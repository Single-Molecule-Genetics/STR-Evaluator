import argparse
import pickle
import gc
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from STR_Functions import *
from TierFunctions import *

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="STR calling from PE bam files"
    )

    # Required inputs
    parser.add_argument(
        "--str_finder_pkl",
        dest="str_finder_pkl_path",
        required=True,
        help="Path to STR-Finder pickle (required).",
    )
    parser.add_argument(
        "--dcsbam",
        dest="dcs_bam_path",
        required=True,
        help="Path to DCS BAM file (required).",
    )
    parser.add_argument(
        "--pebam",
        dest="pe_bam_path",
        required=True,
        help="Path to PE BAM file (required).",
    )
    parser.add_argument(
        "--regions",
        required=True,
        help="Comma-separated list of chromosomes to process (e.g., 'chr1,chr2,chr3'). "
             "Must match names in the BAM files (required).",
    )

    # Optional: cross table + strains (must be provided together)
    parser.add_argument(
        "--crosstable",
        dest="cross_df_path",
        default="",
        help="Path to genotype cross-table Excel (optional; if provided, must be used with --strain1 and --strain2).",
    )
    parser.add_argument(
        "--strain1",
        default="",
        help="Name of the first strain in the cross-table (optional; requires --cross-table and --strain2).",
    )
    parser.add_argument(
        "--strain2",
        default="",
        help="Name of the second strain in the cross-table (optional; requires --cross-table and --strain1).",
    )

    # Other options
    parser.add_argument(
        "--minseqlength",
        dest="min_seq_length",
        type=int,
        default=60,
        help="Minimum sequence length to consider (default: 60).",
    )

    args = parser.parse_args()

    # Normalize/validate regions
    args.regions = [r.strip() for r in args.regions.split(",") if r.strip()]
    if not args.regions:
        parser.error("No valid regions were provided via --regions.")

    # Cross-table/strain consistency checks
    has_cross = bool(args.cross_df_path)
    has_s1 = bool(args.strain1)
    has_s2 = bool(args.strain2)

    if has_cross and (not has_s1 or not has_s2):
        parser.error("When --cross-table is provided, both --strain1 and --strain2 must also be provided.")
    if (has_s1 or has_s2) and not has_cross:
        parser.error("--strain1/--strain2 require --cross-table to be provided.")

    return args


def main():
    args = parse_arguments()
    res_dir = "results"
    dec = 'favor'

    def fast_concat(dfs):
        round_num = 0
        while len(dfs) > 1:
            if round_num % 3 == 0:
                tqdm.write(f"Concatenation round {round_num}: {len(dfs)} DataFrames")
                round_num += 1
            dfs = [
                pd.concat(dfs[i:i + 2], ignore_index=True)
                for i in tqdm(range(0, len(dfs), 2), desc=f"Concat round {round_num}")
            ]
        return dfs[0] if dfs else pd.DataFrame()

    def del_variable(x_v):
        del x_v
        gc.collect()


    res_dir = "results"
    stutter_dir = os.path.join(res_dir, "StutterVsTierPlots")

    os.makedirs(stutter_dir, exist_ok=True)
    for region in args.regions:
        os.makedirs(os.path.join(stutter_dir, region), exist_ok=True)
    if args.cross_df_path:

        cross_sheets = pd.ExcelFile(args.cross_df_path).sheet_names
        unified_regions_cross_df = pd.concat([pd.read_excel(args.cross_df_path, sheet) for sheet in cross_sheets])
        except_strs_df = unified_regions_cross_df[unified_regions_cross_df['type'] == 'STR']
        except_strs = [[f"{row['Chr']}-{row['Position']}-{row['REF']}-{row[args.strain1]}",
                        f"{row['Chr']}-{row['Position']}-{row['REF']}-{row[args.strain2]}"] for _, row in
                       except_strs_df.iterrows()]
        except_strs = sum(except_strs, [])
    else:
        except_strs = []

    with open(args.str_finder_pkl_path, 'rb') as f:
        chromosome_STRdict = pickle.load(f)
    chromosome_STRdict = chromosome_STRdict['chromosome_STRdict']

    dcs_tags = {read.qname for read in pysam.AlignmentFile(args.dcs_bam_path, "rb").fetch(until_eof=True)}

    flat_strs = []

    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        futures = [
            executor.submit(process_region, region, chromosome_STRdict, args.pe_bam_path, args.min_seq_length, dcs_tags)
            for region in args.regions
        ]
        for f in tqdm(as_completed(futures), total=len(futures), desc="Processing flat_strs of regions"):
            flat_strs.extend(f.result())

    del_variable(dcs_tags)

    flat_strs = [lst for lst in flat_strs if lst]
    # flat_strs = list(chain(*flat_strs))iss

    flat_strs_gen = ((i[0], int(i[2]) + 1, i[1], '.'.join(i[4]), (i[1], int(i[2]) + 1)) for i in flat_strs)

    res_df_noTiers = pd.DataFrame(flat_strs_gen, columns=['id', 'start_pos', 'tag', 'tag_nr', 'combined'])

    res_df_noTiers.to_pickle(f"{res_dir}/_{dec}_STRs_dfRaw.pkl", compression="gzip")
    # del_variable(flat_strs)

    groups = res_df_noTiers.groupby(['start_pos', 'tag'])
    dfs_with_ties = []
    dfs_without_ties = []
    for start_pos_tag, group_df in tqdm(groups, desc=f'isolating ties'):
        if len(group_df['id'].unique()) > 1:
            dfs_with_ties.append(group_df)
        else:
            dfs_without_ties.append(group_df)

    del_variable(res_df_noTiers)

    in_args = [(tie_df, except_strs, dec, start_pos_tag) for tie_df in dfs_with_ties]

    with ProcessPoolExecutor() as executor:
        corrected_ties = list(
            tqdm(executor.map(clean_tied_data, in_args), total=len(in_args), desc=f"correcting ties mode: -> {dec}"))

    del_variable(groups)

    all_dfs = dfs_without_ties + corrected_ties

    del_variable(dfs_without_ties)
    del_variable(dfs_with_ties)

    if all_dfs:

        tie_free_df = fast_concat(all_dfs)
    else:
        print(f'NO STRs found')
        # continue
    # print('adding alt column')

    del_variable(all_dfs)

    tie_free_df['alt'] = tie_free_df['id'].str.split('-').str[-1]

    print('groupping tie_free_df')
    grouped = tie_free_df.groupby('tag')
    print('done groupping tie_free_df')
    tag_groups = list(grouped)
    del_variable(tag_groups)
    with ProcessPoolExecutor() as executor:
        results = list(
            tqdm(executor.map(process_tag_group, tag_groups), total=len(tag_groups), desc="processing tag_str"))

    # Convert list of tuples back to dict
    del_variable(tag_groups)
    tag_str = dict(results)

    tag_str_percent = {tag: {repeat: None for repeat in tag_str[tag]} for tag in
                       tag_str}  # this contains the % of an STR presence on a strand (how many ab1 of total ab1 have the STR (for those that actually have a call)
    all_reads = read_generator(args.pe_bam_path, args.min_seq_length, None, dcs_tags)
    reads_dict = len_read_types(all_reads)
    for tag in tqdm(tag_str, desc="processing tag_str_percent"):
        repeat_mInfo = tag_str[tag]

        # Original counts
        ab1 = (reads_dict['ab.1'][tag] if tag in reads_dict['ab.1'] else 0)
        ba2 = (reads_dict['ba.2'][tag] if tag in reads_dict['ba.2'] else 0)
        ab2 = (reads_dict['ab.2'][tag] if tag in reads_dict['ab.2'] else 0)
        ba1 = (reads_dict['ba.1'][tag] if tag in reads_dict['ba.1'] else 0)

        for repeat in repeat_mInfo:
            ab1P = (0 if repeat_mInfo[repeat]['ab1'] == 0 else round(repeat_mInfo[repeat]['ab1'] / ab1, 2))
            ba2P = (0 if repeat_mInfo[repeat]['ba2'] == 0 else round(repeat_mInfo[repeat]['ba2'] / ba2, 2))
            ab2P = (0 if repeat_mInfo[repeat]['ab2'] == 0 else round(repeat_mInfo[repeat]['ab2'] / ab2, 2))
            ba1P = (0 if repeat_mInfo[repeat]['ba1'] == 0 else round(repeat_mInfo[repeat]['ba1'] / ba1, 2))
            tag_str_percent[tag][repeat] = (ab1P, ba2P, ab2P, ba1P)

    flat_str_tiers = []
    for tag in tqdm(tag_str, desc='calculating frequencies and tiers'):
        oab1 = (reads_dict['ab.1'][tag] if tag in reads_dict['ab.1'] else 0)
        oba2 = (reads_dict['ba.2'][tag] if tag in reads_dict['ba.2'] else 0)
        oab2 = (reads_dict['ab.2'][tag] if tag in reads_dict['ab.2'] else 0)
        oba1 = (reads_dict['ba.1'][tag] if tag in reads_dict['ba.1'] else 0)

        for repeat in tag_str[tag]:
            fs = sum(tag_str[tag][repeat].values())
            counts = tuple(tag_str[tag][repeat].values())
            frequencies = tag_str_percent[tag][repeat]
            ab1, ba2, ab2, ba1 = counts
            tag_nr = (1 if ab1 > 0 and ba2 > 0 and ab2 == 0 and ba1 == 0 else 2)

            res = (repeat, int(repeat.split('-')[1]), tag, tier_assigner(counts, frequencies), ab1, ba2, ab2, ba1,
                   frequencies[0], frequencies[1], frequencies[2], frequencies[3], oab1, oba2, oab2, oba1, tag_nr)
            flat_str_tiers.append(res)

            #       0       1   2           3       4       5       6       7      8    9     10      11      12     13    14      15     16      17    18     19        20
    output_df_columns = ['chr', 'id', 'start_pos', 'tag', 'tier', '#ab1', '#ba2', '#ab2', '#ba1', '%ab1', '%ba2',
                         '%ab2', '%ba1', 'oab1', 'oba2', 'oab2', 'oba1', 'fs', 'oFs', 'tag_nr', 'cnr']

    # flat_str_tiers = [[i[0].split('-')[0]] + list(i[:-1]) + [sum(i[4:8]), sum(i[12:16])  , i[16] ,int(i[0].split('-')[0].split('chr')[-1]) ] for i in flat_str_tiers ]
    try:
        flat_str_tiers = [
            [i[0].split('-')[0]] + list(i[:-1]) + [max(i[4:6]) + max(i[6:8]), max(i[12:14]) + max(i[14:16]), i[16],
                                                   int(i[0].split('-')[0].split('chr')[-1])] for i in flat_str_tiers]
    except ValueError:
        flat_str_tiers = [
            [i[0].split('-')[0]] + list(i[:-1]) + [max(i[4:6]) + max(i[6:8]), max(i[12:14]) + max(i[14:16]), i[16],
                                                   int(1)] for i in flat_str_tiers]

    # Create a nested dictionary
    # stutter_rates = (defaultdict(lambda: defaultdict(dict)))

    # Iterate through each row and populate the dictionary
    # for _, row in stutter_df.iterrows():
    #     lib_type = row['lib_type']
    #     repeat_type = row['repeat_type']
    #     repeat_len = row['repeat_len']
    #     kde_value = row['KDE']

    # Add values without overwriting existing ones
    #     stutter_rates[lib_type][repeat_type][repeat_len] = kde_value
    #
    # # Convert back to a regular dictionary (optional)
    # stutter_rates = dict(stutter_rates)

    # MAX_ROWS_PER_SHEET = 1040000
    # print('writing summary')
    # with (pd.ExcelWriter(f"{res_dir}/_{dec}_STRs.xlsx") as writer):
    #     cat = 'STRs'
    #     # Create the output DataFrame

    output_df = pd.DataFrame(flat_str_tiers, columns=output_df_columns)

    output_df.sort_values(['cnr', 'start_pos', 'tier'], ascending=True, inplace=True)
    output_df['chr-pos'] = output_df['id'].apply(lambda x: "-".join(x.split("-")[:2]))
    pos_cvg = output_df['chr-pos'].value_counts().to_dict()
    output_df['pos_coverage'] = output_df['chr-pos'].map(pos_cvg)
    # if min_cvg:
    #     output_df = output_df[output_df['pos_coverage'] >= min_cvg]

    id_counts = output_df['id'].value_counts().to_dict()
    output_df['vaf'] = output_df['id'].map(id_counts) / output_df['chr-pos'].map(pos_cvg)
    output_df['ref'] = [i.split('-')[-2] for i in output_df['id']]
    output_df['alt'] = [i.split('-')[-1] for i in output_df['id']]
    output_df['tier_score'] = [get_tier_score(tier) for tier in output_df['tier']]
    # med_lib_cvg = [med_region_cvg[reg] for reg in output_df['chr']]
    # output_df['med_lib_cvg'] = med_lib_cvg
    # avg_lib_cvg = [avg_region_cvg[reg] for reg in output_df['chr']]
    # output_df['avg_lib_cvg'] = avg_lib_cvg
    # output_df['repeat_type'] = [get_motif_type(i) for i in output_df['alt']]
    output_df['heterology'] = output_df.groupby('start_pos')['vaf'].transform(lambda x: (x > 0.4).any())
    output_df['heterology'] = output_df.groupby('start_pos')['vaf'].transform(lambda x: ((x > 0.4) & (x < 0.7)).any())
    # lib_type = 'TargetedF' if lib[:2]  in ['01','02','03'] else 'RandomF'

    output_df['repeat_size'] = [int(return_digit(x.split('-')[-1]).split('-')[0])
                                for x in output_df['alt']]
    output_df['repeat_size_type'] = ['long' if i > 8 else 'short' for i in output_df['repeat_size']]

    # output_df['rate'] = [
    #     stutter_rates.get(lib_type, {}).get(row['repeat_type'], {}).get(row['repeat_size_type'], 0.02115)
    #     if not row['heterology']
    #     else stutter_rates.get(lib_type, {}).get(row['repeat_type'], {}).get(row['repeat_size_type'], 0.021155) / 2
    #     for _, row in output_df.iterrows()
    # ]
    # output_df['Bayes'] = [int(calculate_posterior(row['tier_score'],row['pos_coverage'],row['med_lib_cvg'],row['vaf'],row['rate'], min_bay )[0]) for _,row in tqdm(output_df.iterrows(), total=len(output_df))]

    # output_df['P(True)'] = [float(calculate_posterior(row['tier_score'],row['pos_coverage'],row['med_lib_cvg'],row['vaf'],row['rate'],min_bay)[1]) for _,row in tqdm(output_df.iterrows(), total=len(output_df))]

    heterology_df = output_df[output_df['heterology'] == True]['id'].unique()
    vaf_dic = {row['id']: row['vaf'] for _, row in output_df.iterrows()}

    tie_free_df['vaf'] = [vaf_dic[id] if id in vaf_dic else 0 for id in tie_free_df['id']]
    tie_free_df['heterology'] = [True if id_ in heterology_df else False for id_ in tie_free_df['id']]
    tie_free_df['alt'] = [i.split('-')[-1] for i in tie_free_df['id']]
    tie_free_df.to_pickle(f"{res_dir}/_{dec}_{'' if 'SSCS' not in args.pe_bam_path else 'sscs'}STRs_dfRaw.pkl",
                          compression="gzip")

    del_variable(tie_free_df)
    # print('-->EXCEL')
    # Split the output_df DataFrame into multiple sheets if necessary
    # num_sheets = math.ceil(len(output_df) / MAX_ROWS_PER_SHEET)
    # for i in range(num_sheets):
    #     start = i * MAX_ROWS_PER_SHEET
    #     end = min(start + MAX_ROWS_PER_SHEET, len(output_df))
    #     sheet_name = f"{cat}({i+1})"
    #     output_df[start:end].to_excel(writer, sheet_name=sheet_name, index=False)

    print('saving pkl')
    output_df.to_pickle(f"{res_dir}/_STRSummaryCash.pkl", compression="gzip")
    print('compressing summary')
    compressed_df = compress_strs2(output_df, 'library', res_dir)

    output_df['type'] = [i.split('-')[-1] for i in output_df['id']]

    fsd_dic = {repeat: {typ: Counter([(tier, fsd) for tier, fsd in
                                      zip(output_df[(output_df['id'] == repeat) & (output_df['type'] == typ)]['tier'],
                                          output_df[(output_df['id'] == repeat) & (output_df['type'] == typ)]['fs'])])
                        for typ in output_df[output_df['id'] == repeat]['type'].unique()} for repeat in
               output_df['id'].unique()}

    del_variable(output_df)
    repeat_figDf = {repeat: [] for repeat in fsd_dic}

    for repeat in fsd_dic:
        for typ in fsd_dic[repeat]:
            temp_df = pd.DataFrame(columns=['Repeat', 'type', 'fsd', 'Tier', 'count'])
            temp_df['fsd'] = [i[1] if i[1] < 20 else 20 for i in fsd_dic[repeat][typ]]
            temp_df['Tier'] = [i[0] for i in fsd_dic[repeat][typ]]
            temp_df['count'] = [fsd_dic[repeat][typ][i] for i in fsd_dic[repeat][typ]]
            temp_df['Repeat'] = [f"{typ}" for i in fsd_dic[repeat][typ]]
            temp_df['type'] = [typ for i in fsd_dic[repeat][typ]]
            repeat_figDf[repeat].append(temp_df)
        repeat_figDf[repeat] = sorted(repeat_figDf[repeat], key=custom_sort_key)

    legend_handles = [Rectangle((0, 0), 1, 1, color=color) for color in base_colors.values()]
    legend_labels = [f"{key}: {tier_description[key]}" for key in base_colors.keys()]
    for repeat in tqdm(repeat_figDf, desc='Plotting stutter'):

        repeat_chr, repeat_pos, _, _ = repeat.split('-')
        typs_nr = len(repeat_figDf[repeat])
        if typs_nr > 1:

            n_cols = (typs_nr if typs_nr <= 3 else 3)
            t_rows = 0
            if typs_nr > 3:
                while typs_nr > 0:
                    typs_nr -= n_cols
                    t_rows += 1
            else:
                t_rows = 1
            n_rows = t_rows
            fig, axs = plt.subplots(n_rows, n_cols, figsize=(30, 70))
            axe = axs.ravel()

            for idx, df in enumerate(repeat_figDf[repeat]):

                df_grouped = df.groupby(['Repeat', 'fsd', 'Tier'], as_index=False).agg({'count': 'sum'})

                df_pivot = df_grouped.pivot(index=['Repeat', 'fsd'], columns='Tier', values='count').fillna(0)
                cols_order = {tier: get_tier_order(tier) for tier in df_pivot.columns}
                sorted_cols = list([i[0] for i in sorted(list(cols_order.items()), key=lambda x: x[1])])
                df_pivot = df_pivot[sorted_cols]

                fsd_total = df_pivot.sum(axis=1)

                df_pct = df_pivot.divide(fsd_total, axis=0)
                colors = {i: get_adjusted_color(int(str(i)[0]), i) for i in df_pct.columns}

                n_bars = len(df_pct.index.tolist())

                # df_pct.plot(kind='barh', stacked=True, color=colors, width=0.9, ax=axe[idx], legend=False, grid=False)
                df_pct.plot(kind='barh', stacked=True,
                            # color=colors,
                            width=0.9, ax=axe[idx], legend=False, grid=False)

                axe[idx].set_xlim(right=1.1)
                axe[idx].set_title(f"{df_pct.index.tolist()[0][0]}  ({int(sum(fsd_total))})", fontsize=26)

                for i, v in enumerate(fsd_total):
                    axe[idx].text(1.05, i, str(int(v)), color='black', fontsize=18, ha='center', va='center')
                axe[idx].tick_params(labelsize=15)

                y_labels = [f"FS {entry[1]}" for entry in df_pct.index.tolist()]
                axe[idx].set_yticklabels(y_labels, fontsize=20)
                axe[idx].set_ylabel("")


        else:
            df = repeat_figDf[repeat][0]
            fig, axs = plt.subplots(figsize=(13, 10))

            df_grouped = df.groupby(['Repeat', 'fsd', 'Tier'], as_index=False).agg({'count': 'sum'})

            df_pivot = df_grouped.pivot(index=['Repeat', 'fsd'], columns='Tier', values='count').fillna(0)

            fsd_total = df_pivot.sum(axis=1)
            cols_order = {tier: get_tier_order(tier) for tier in df_pivot.columns}
            sorted_cols = list([i[0] for i in sorted(list(cols_order.items()), key=lambda x: x[1])])
            df_pivot = df_pivot[sorted_cols]

            df_pct = df_pivot.divide(fsd_total, axis=0)
            colors = {i: get_adjusted_color(int(str(i)[0]), i) for i in df_pct.columns}

            df_pct.plot(kind='barh', stacked=True,
                        color=[get_adjusted_color(int(str(i)[0]), i) for i in df_pct.columns],
                        width=0.8, ax=axs, legend=False)

            for i, v in enumerate(fsd_total):
                axs.text(1.05, i, str(int(v)), color='black', fontsize=18, ha='center', va='center')

            y_labels = [f"FS {entry[1]}" for entry in df_pct.index.tolist()]
            axs.set_yticklabels(y_labels, fontsize=20)
            axs.set_ylabel("")
            axs.set_title(f"{df_pct.index.tolist()[0][0]}  ({int(sum(fsd_total))})", fontsize=26)

            axs.set_xlim(right=1.1)
            axs.tick_params(labelsize=20)
            axs.grid(False)

        fig.legend(handles=legend_handles, labels=legend_labels, title='Tier',
                   bbox_to_anchor=(0.94, 0.5), loc='center', fontsize=10)
        plt.suptitle(repeat, fontsize=30)
        plt.savefig(f"{stutter_dir}/{repeat_chr}/_{dec}_{repeat}_fsTd.pdf")
        plt.close(fig)

    for region in args.regions:
        merge_pdfs(sorted(glob.glob(stutter_dir + '/' + region + '/*')),
                   stutter_dir + '/' + region + '_library' + '_' + region + '_StutterVsFS_Tiers')
        shutil.rmtree((stutter_dir + '/' + region))


if __name__ == "__main__":
    main()
