# %%
import os
import random
from bisect import bisect_left, bisect_right
from collections import Counter
from itertools import chain
from typing import Dict, Optional
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from tqdm import tqdm
import glob
from Shared_Functions import merge_pdfs
import shutil
import math
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import re
from pyfasta import Fasta
from Shared_Functions import return_digit
from TierFunctions import tier_assigner, tier_score, tier_assigner_original

min_mono = 5
min_di = 6
min_tri = 4

mono_regex = r"(?:A{5,}|T{5,}|C{5,}|G{5,}|a{5,}|t{5,}|c{5,}|g{5,})"
di_regex = r"((?:([ACGTacgt])(?!\2)[ACGTacgt]))(?:\1){" + str(min_di - 1) + ",}"
tri_regex = r"([ACGT]{3})\1{" + str(min_tri - 1) + ",}"


# %%
def color_column(col):
    color = STR_Tiers_Colors.get(col.name)
    if color:
        return [f'background-color: {color}' if val > 0 else '' for val in col]
    return ['' for _ in col]


# %%
ordered_tiers = []  # All tier combinations
for A in range(1, 6):  # First digit ranges from 1 to 5
    for B in range(9, -1, -1):  # Second digit ranges from 9 down to 0
        for C in range(B, -1, -1):  # Third digit ranges from current B down to 0
            for D in range(5, 0, -1):  # Fourth digit ranges from 5 down to 1
                for E in range(D, -1, -1):  # Fifth digit ranges from 5 down to 1
                    # Apply the compound condition
                    if ((A <= 3 and (B == 0 or C == 0))
                            or (A > 3 and (B > 0 and C > 0))
                            or (B == C == 0)
                            or (A <= 3 and ((D == 0) or (E == 0)))
                            or (A > 3 and ((D > 0) and (E > 0)))
                    ):
                        continue
                    # Construct the number and add to the list
                    if A > 3:
                        number = float(f"{A}.{B}{C}0{max(E, D)}")

                    else:
                        number = float(f"{A}.{B}{C}{D}{E}")
                    ordered_tiers.append(number)
# %%
STR_Tiers_Colors = {1.1: '#006400', 1.2: '#008000', 1.3: '#4C944C',
                    2.1: '#85EE1D', 2.2: '#99E54E', 2.3: '#A9E86A', 2.4: '#B8EC84',
                    3.1: '#FFEA22', 3.2: '#FFED41', 3.3: '#FFEF5A', 3.4: '#FFF16B', 3.5: '#FFF37C', 3.6: '#FFF7A2',
                    3.7: '#FFF9B9',
                    4.1: '#FF8C00', 4.2: '#FFA024', 4.3: '#FFAA39', 4.4: '#FFB555', 4.5: '#FDC47A',
                    5.1: '#7F00FF', 5.2: '#8D1AFF', 5.3: '#9830FF', 5.4: '#A54AFF', 5.5: '#B05FFF', 5.6: '#BD7AFF',
                    5.7: '#CB97FF',
                    7.0: '#FF0000'}

labels_colors_dic = {
    'NR1': 'green',
    'NR2': 'orange',
    'CO_R1': 'pink',
    'CO_R2': 'pink',
    'NCO': '#00CCCC',
    'L_SNPS0': 'blue',
    'R_SNPS0': 'blue',
    'SNPs L|R 1': '#FF6868',
    'other': 'grey',
    'All': 'grey',
    'NR1_HQ': 'green',
    'NR2_HQ': 'orange',
    'CO_R1_HQ': 'pink',
    'CO_R2_HQ': 'pink',
    'NCO_HQ': '#00CCCC',
    'L_SNPS0_HQ': 'blue',
    'R_SNPS0_HQ': 'blue',
    'SNPs L|R 1_HQ': '#FF6868',
    'other_HQ': 'grey',
    'All_HQ': 'grey',
    'Tier6': 'purple'
}

# %%
from collections import Counter


def compress_strs(str_file_path, name=None, res_dir=None):
    rows_info = []
    compress_df_group = str_file_path.groupby('id')
    cvg_df = str_file_path.groupby(['chr', 'start_pos'])

    cvg_values = {chro: {} for chro in str_file_path['chr'].unique()}
    for (chrom, pos), sub_df in cvg_df:
        cvg_values[chrom][pos] = len(sub_df.index)

    for repeat_id, repeat_df in compress_df_group:

        region, pos, ref, alt = repeat_id.split('-')
        pos = int(pos)

        AC = len(repeat_df['tag'].unique())
        cvg = cvg_values[region][pos]
        AF = AC / cvg

        tier_counts = Counter(repeat_df['tier'])
        for tier in list(STR_Tiers_Colors):
            if tier not in tier_counts:
                tier_counts[tier] = 0
        rows_info.append([repeat_id, ('REF' if ref == alt else 'ALT'), cvg, AC, AF] + [tier_counts[tier] for tier in
                                                                                       STR_Tiers_Colors] + [pos,
                                                                                                            int(
                                                                                                                region.split(
                                                                                                                    'chr')[
                                                                                                                    -1])])
    compressed_df = pd.DataFrame(rows_info,
                                 columns=['STR ID', 'type', 'cvrg', 'AC', 'AF'] + list(STR_Tiers_Colors) + ['start_pos',
                                                                                                            'chr'])
    compressed_df.sort_values(by=['chr', 'start_pos'], ascending=True, inplace=True)

    if name:
        output_file = f"{res_dir}/{name}_STRFreqFile.xlsx"
    else:
        output_file = f"{str_file_path.split('/')[-1][:-15]}STRFreqFile.xlsx"

    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        # Write the DataFrame to the Excel file
        compressed_df.to_excel(writer, sheet_name='Sheet1', index=False)

        workbook = writer.book
        worksheet = writer.sheets['Sheet1']

        for col_name, color in STR_Tiers_Colors.items():
            if col_name in compressed_df.columns:
                col_letter = compressed_df.columns.get_loc(col_name)  # + 1
                worksheet.conditional_format(1, col_letter, compressed_df.shape[0], col_letter, {
                    'type': 'cell',
                    'criteria': '>',
                    'value': 0,
                    'format': workbook.add_format({'bg_color': color}),
                })
    return compressed_df


# Function to calculate the color adjustment
# def get_adjusted_color(group, tier_value):
#     # Get the base color for the main group
#     base_color = base_colors.get(group)
#     if not base_color:
#         raise ValueError("Group not recognized")
#
#     # Convert base color to RGB
#     base_rgb = mcolors.hex2color(base_color)
#
#     # Split tier_value into main group (A), frequency (BC), and family size (D)
#     main_tier, decimals = divmod(tier_value, 1)
#     group = int(main_tier)  # A is the main group
#     bc_d = str(int(decimals * 1000)).zfill(3)  # Get the last three digits (BCD)
#     b, c, d = int(bc_d[0]), int(bc_d[1]), int(bc_d[2])  # Separate B, C, and D
#
#     # Determine BC confidence: treat as high if neither B nor C is 0, otherwise low
#     bc_confidence = 1 if (b != 0 and c != 0) else 0  # 1 for high confidence, 0 for low
#
#     # Adjust based on D (family size): Higher D makes it darker (closer to the main color), lower D makes it lighter
#     d_brightness_factor = 0.5 - (d - 1) * 0.1  # D=1 gives 0.5 (lighter), D=5 gives 0.1 (darker)
#
#     # Apply the color adjustment
#     adjusted_rgb = [
#         min(1, channel + d_brightness_factor * (1 - bc_confidence))  # Increase brightness if low confidence
#         for channel in base_rgb
#     ]
#
#     # Convert back to hex color
#     return mcolors.to_hex(adjusted_rgb)


# %%
# Define the function to get the GC content using a 40-base window (20 before and after)


# %%


# %%
def clean_tied_data(tpl):
    df = tpl[0]
    expected_ids = tpl[1]  #
    decision = tpl[2]
    current_pos = tpl[3]
    id_counts = df['id'].value_counts()
    max_id_count = id_counts.max()

    if (id_counts == max_id_count).sum() > 1:
        tied_ids = id_counts[id_counts == max_id_count].index.tolist()
        random_id = [random.choice(tied_ids)]

        if decision == 'random':
            df_tied_ids = df[df['id'].isin(random_id)]
            return df_tied_ids

        elif decision == 'favor':
            common_ids = list(set(tied_ids).intersection(expected_ids))

            if common_ids:
                if len(common_ids) > 1:
                    random_common_id = random.choice(common_ids)
                    df_tied_ids = df[df['id'] == random_common_id]
                else:
                    df_tied_ids = df[df['id'].isin(common_ids)]

                return df_tied_ids
            else:
                df_tied_ids = df[df['id'].isin(random_id)]
                return df_tied_ids

        elif decision == 'remove':
            # Remove the information for tied IDs
            df_tied_ids = df[df['start_pos'] != current_pos]
            return df_tied_ids
    else:
        max_id = id_counts.idxmax()
        df_max_id = df[df['id'] == max_id]
        return df_max_id


# %%
def summary_file_prep(xfile, chromosome, region_span, file_type):
    xfile = xfile.fillna(0)
    xfile = xfile[xfile['variant ID'] != 0]  # remove nans
    xfile['chr'] = [var_id.split('-')[0] for var_id in xfile['variant ID']]
    if file_type == 'SNPs':
        xfile['start_pos'] = [
            int(var_id.split('-')[1]) for var_id in xfile['variant ID']]
    else:
        xfile['start_pos'] = [
            int(var_id.split('-')[1]) + 1 if len(var_id.split('-')[-2]) > len(var_id.split('-')[-1]) else int(
                var_id.split('-')[1]) for var_id in xfile['variant ID']]

    xfile = xfile[
        (xfile['chr'] == chromosome) & (xfile['start_pos'].isin(
            region_span))]  # -Limit to one chromosome of interest ( this is important for haplotypes)

    grouped = xfile.groupby(['tag', 'variant ID'])
    to_remove = []
    for _, group in grouped:
        alt_exists = 'alt' in group['allele'].values
        ref_exists = 'ref' in group['allele'].values

        # Check if both 'ref' and 'alt' alleles exist in the group
        if alt_exists and ref_exists:
            # Find the index of 'ref' allele and add it to the indices to be removed
            to_remove.append(group[group['allele'] == 'ref'].index[0])
    xfile = xfile.drop(to_remove).reset_index(drop=True)
    new_alt = []
    for var_id, tier, allele in zip(xfile['variant ID'], xfile['tier'], xfile['allele']):
        temp_ref = var_id.split('-')[-2]
        temp_alt = var_id.split('-')[-1]

        if allele == 'ref':
            print('AA')
            if len(temp_ref) > len(temp_alt):  # Deletions
                if tier < 3:
                    new_alt.append(temp_ref[1])
                else:
                    new_alt.append(temp_ref[1].lower())

            else:  # point mutation
                if tier < 3:
                    new_alt.append(temp_ref)
                else:
                    new_alt.append(temp_ref.lower())

        else:  # alt
            if file_type == 'SNPs':
                if tier < 3:
                    new_alt.append(temp_alt)
                else:
                    new_alt.append(temp_alt.lower())

            else:
                if len(temp_ref) > len(temp_alt):  # Deletions
                    if tier < 3:
                        new_alt.append('D')
                    else:
                        new_alt.append('d')

                elif len(temp_ref) < len(temp_alt):  # Insertions
                    if tier < 3:
                        new_alt.append(temp_alt[1:])
                    else:
                        new_alt.append(temp_alt[1:].lower())

                else:  # point mutation
                    if tier < 3:
                        new_alt.append(temp_alt)
                    else:
                        new_alt.append(temp_alt.lower())
    xfile['alt'] = new_alt
    return xfile

# %%
def compare_w_strain(indexes, row_lst, target):  # takes indexes of intrest and compare them from 2 lists
    target_elements = [target[idx] for idx in indexes]
    return sum([1 for t, r in zip(target_elements, row_lst) if t == r])


# %%
def crossing_infoV2(genotype_df, cross, given_genotype=None):
    # print(genotype_df)

    strain1_name = cross[:2]
    strain2_name = cross[2:]
    if given_genotype is None:
        genotype_df = genotype_df[~genotype_df['type'].isin(['ISM', 'AMP'])]
        given_genotype = genotype_df['Position'].tolist()

    strain_1_bases = genotype_df[strain1_name]
    strain_2_bases = genotype_df[strain2_name]

    strain_1_refs = {pos: ref for pos, ref in zip(given_genotype, strain_1_bases)}  # compact
    strain_2_refs = {pos: ref for pos, ref in zip(given_genotype, strain_2_bases)}

    result_dic = {
        strain1_name: (strain_1_bases, strain_1_refs,),
        strain2_name: (strain_2_bases, strain_2_refs),
        # add more when needed
    }

    legend_image = f"/Users/shehabmoukbel/Desktop/11_Bioinformatics Analysis/7_Resources/9_CrossesPngs/{cross}.png"

    return (result_dic[strain1_name][0],
            result_dic[strain1_name][1],
            strain1_name,
            result_dic[strain2_name][0],
            result_dic[strain2_name][1],
            strain2_name,
            legend_image)


def crossing_info(cross, given_genotype=None):
    if given_genotype is None:
        given_genotype = sorted([695296,
                                 695330,
                                 695271,
                                 695368,
                                 695530,
                                 695406,
                                 695481,
                                 695515,
                                 695319,
                                 695389,
                                 695485])
    w1 = ['A', 'C', '6A', 'G', 'T', '5T', 'T', 'T', '5A', 'G', 'A']  # raw
    s1 = ['G', 'T', '5A', 'C', 'G', '5T', 'G', 'C', '13A', 'T', 'C']
    w2 = ['A', 'C', '13A', 'G', 'T', '6T', 'T', 'T', '6A', 'G', 'A']
    s2 = ['G', 'T', '5A', 'C', 'G', '5T', 'G', 'C', '6A', 'T', 'C']

    # add more when needed

    w1_r = {pos: ref for pos, ref in zip(given_genotype, w1)}  # compact
    s1_r = {pos: ref for pos, ref in zip(given_genotype, s1)}
    w2_r = {pos: ref for pos, ref in zip(given_genotype, w2)}
    s2_r = {pos: ref for pos, ref in zip(given_genotype, s2)}

    result_dic = {
        'w1': (w1, w1_r,),
        's1': (s1, s1_r),
        'w2': (w2, w2_r),
        's2': (s2, s2_r)
        # add more when needed
    }
    strain1 = cross[:2]
    strain2 = cross[2:]
    legend_image = f"/Users/shehabmoukbel/Desktop/STR_Calling/STR Presentations and Posters/{cross}.png"

    return (result_dic[strain1][0],
            result_dic[strain1][1],
            strain1,
            result_dic[strain2][0],
            result_dic[strain2][1],
            strain2,
            legend_image)


# %%
def plot_tiers_counts(files):
    # Initialize an empty DataFrame to store combined data
    combined_df = pd.DataFrame(columns=['File', 'Tiers', 'Counts'])

    for file in files:
        # Read the Excel file into a pandas DataFrame
        df = pd.read_excel(file)
        # Extract tiers and counts from the DataFrame
        tiers = [str(i) for i in df['Unnamed: 0']]
        counts = df['counts']

        # Create a new DataFrame for the current file
        file_df = pd.DataFrame({'File': os.path.basename(file).split('_')[-3],
                                'Tiers': tiers,
                                'Counts': counts})
        # Append the data from the current file to the combined DataFrame
        combined_df = combined_df.append(file_df, ignore_index=True)

    # Calculate the total count for all files
    total_count = combined_df['Counts'].sum()

    # Pivot the DataFrame to get the data in a suitable format for grouped bars
    pivot_df = combined_df.pivot(index='Tiers', columns='File', values='Counts').reset_index()

    # Plot the grouped bars and add percentages to the bars
    fig, ax = plt.subplots(figsize=(12, 6))
    pivot_df.plot(x='Tiers', kind='bar', ax=ax)

    ax.set_title(f"Tiers counts")
    ax.set_xlabel('Tiers')
    ax.set_ylabel('Counts')
    ax.set_yscale('log')
    ax.legend(title='File', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Display the plot
    plt.tight_layout()
    plt.savefig(''.join(file.split('xlsx')) + 'pdf')
    # plt.show()


# %%
def plot_stutter(df, res_dir, name, fsize=(10, 6)):
    try:
        groups = df.groupby(['tag', 'start_pos'])
    except AttributeError:
        df = pd.read_pickle(df)
        groups = df.groupby(['tag', 'start_pos'])
    tag_idsPerPos = {tag: {} for tag in df['tag'].unique()}
    for group in groups:
        tag = group[0][0]
        pos = group[0][1]
        ids_size = len(group[1]['id'].unique())
        if ids_size > 1:
            tag_idsPerPos[tag][pos] = ids_size
    tag_idsPerPos = {tag: tag_idsPerPos[tag] for tag in tag_idsPerPos if tag_idsPerPos[tag]}
    multi_calls_poses = Counter(sum([list(tag_idsPerPos[tag]) for tag in tag_idsPerPos], []))

    #####
    groups = df.groupby(['start_pos', 'tag'])
    groups_dict = {tpl: df for tpl, df in groups}
    tag_pos_fs = {tag: {} for tag in tag_idsPerPos}
    for tag in tag_idsPerPos:
        poses = list(tag_idsPerPos[tag])
        for p in poses:
            tag_pos_fs[tag][p] = len(groups_dict.get((p, tag)))

    tie_tags = []
    for tpl, df in groups:
        id_counts = df['id'].value_counts()
        max_id_count = id_counts.max()
        if (id_counts == max_id_count).sum() > 1:
            tie_tags.append(tpl)
    tie_tags = set(tie_tags)

    # plotting
    stutter_counts = Counter(sum([list(tag_pos_fs[tag].values()) for tag in tag_pos_fs], []))
    tie_stutter_counts = Counter([tag_pos_fs[tpl[1]][tpl[0]] for tpl in tie_tags])
    X_stutter = [str(i) for i in sorted(stutter_counts)] + ['']
    Y_stutter = [stutter_counts[x] for x in sorted(stutter_counts)]
    sum_stutter = sum(Y_stutter)
    Y_stutter_pct = [(i / sum_stutter) * 100 for i in Y_stutter] + [0]

    X_stutter_tie = [str(i) for i in sorted(tie_stutter_counts)] + ['']
    Y_stutter_tie = [tie_stutter_counts[x] for x in sorted(tie_stutter_counts)]
    sum_stutter_tie = sum(Y_stutter_tie)
    Y_stutter_tie_pct = [(i / sum_stutter_tie) * 100 for i in Y_stutter_tie] + [0]

    X_multi_calls = sorted(multi_calls_poses)
    Y_multi_calls = [multi_calls_poses[x] for x in X_multi_calls]
    sum_multi_calls = sum(Y_multi_calls)
    Y_multi_calls_pct = [(i / sum_multi_calls) * 100 for i in Y_multi_calls]

    plt.figure(figsize=fsize)
    plt.plot(X_stutter, Y_stutter_pct, color='blue')
    plt.plot(X_stutter_tie, Y_stutter_tie_pct, color='red')
    plt.bar(X_multi_calls, Y_multi_calls_pct)

    plt.title(name + ' STR Stutter/Multicalls FS')

    plt.xticks(X_stutter + X_multi_calls, size=15, rotation=90)

    plt.ylabel('%', size=15)
    plt.xlabel(f"FS - Positions on which multicalls occur, total multicalls {sum_stutter}", size=11)
    custom_labels = ['Stutter % vs. FS', 'Tie within stutter']
    plt.legend(custom_labels, title='Custom Legend Title')

    plt.tight_layout()
    plt.savefig(f"{res_dir}/{name}_MulticallsPos.pdf")


def custom_sort_key(df):
    repeat_value = df['Repeat'].iloc[0]  # Assuming the 'repeat' column is consistent within each DataFrame
    letter = ''.join([i for i in repeat_value if not i.isdigit()])
    number = int(''.join([i for i in repeat_value if i.isdigit()]))
    return letter, number

def process_tag_group(tag_group):
    tag, group = tag_group
    temp_counts = group.groupby(['id', 'tag_nr']).size().to_dict()
    temp_dic = {
        repeat: {
            'ab1': temp_counts.get((repeat, 'ab.1'), 0),
            'ba2': temp_counts.get((repeat, 'ba.2'), 0),
            'ab2': temp_counts.get((repeat, 'ab.2'), 0),
            'ba1': temp_counts.get((repeat, 'ba.1'), 0),
        }
        for repeat in set(group['id'])
    }
    return tag, temp_dic

def compute_tag_percent_chunk(tag_str_chunk, reads_dict, tag_str_percent):
    out = []
    for tag, repeat_mInfo in tag_str_chunk.items():
        ab1 = reads_dict['ab.1'].get(tag, 0)
        ba2 = reads_dict['ba.2'].get(tag, 0)
        ab2 = reads_dict['ab.2'].get(tag, 0)
        ba1 = reads_dict['ba.1'].get(tag, 0)

        result = {}
        for repeat in repeat_mInfo:
            r_info = repeat_mInfo[repeat]
            ab1P = 0 if r_info['ab1'] == 0 or ab1 == 0 else round(r_info['ab1'] / ab1, 2)
            ba2P = 0 if r_info['ba2'] == 0 or ba2 == 0 else round(r_info['ba2'] / ba2, 2)
            ab2P = 0 if r_info['ab2'] == 0 or ab2 == 0 else round(r_info['ab2'] / ab2, 2)
            ba1P = 0 if r_info['ba1'] == 0 or ba1 == 0 else round(r_info['ba1'] / ba1, 2)
            result[repeat] = (ab1P, ba2P, ab2P, ba1P)
        out.append((tag, result))
    return out

    return tag, result
def process_read(tpl):
    pref_suf = 2

    tag, tag_typ, chromosome, pos_lst, seq, nums, chr_dict = tpl

    mono = [repeat for repeat in re.finditer(mono_regex, seq)
            if (repeat.start() >= pref_suf and len(seq) - repeat.end() >= pref_suf and 'N' not in repeat.group())]
    di = [repeat for repeat in re.finditer(di_regex, seq)
          if (repeat.start() >= pref_suf and len(seq) - repeat.end() >= pref_suf and 'N' not in repeat.group())]
    tri_repeats = [repeat for repeat in re.finditer(tri_regex, seq)
                   if
                   (repeat.start() >= pref_suf and len(seq) - repeat.end() >= pref_suf and 'N' not in repeat.group())]
    tri = [match for match in tri_repeats if len(set(match.group())) != 1]
    repeats_before = sum([mono if mono else [], di if di else [], tri if tri else []], [])

    if repeats_before:
        repeats_se = [(repeat.start(), repeat.end()) for repeat in repeats_before]

        repeat_positions = [list(i for i in pos_lst[tpl[0]:tpl[1]] if isinstance(i, int) if i in nums) for tpl in
                            repeats_se]
        non_ref_strs_positions = [list(i for i in pos_lst[tpl[0]:tpl[1]] if isinstance(i, int)) for tpl in
                                  repeats_se]  # same as above but here without if statment. just to see the positions

        def non_ref_strs_func(tag, missing, tag_typ, r):
            return tag, missing, tag_typ, r

        non_ref_strs = [non_ref_strs_func(tag, missing, tag_typ, r) for r, s, missing in
                        zip(repeats_before, repeat_positions, non_ref_strs_positions) if not s]
        repeats = [r for r, s, missing in zip(repeats_before, repeat_positions, non_ref_strs_positions) if s]
        subsets = [s for r, s, missing in zip(repeats_before, repeat_positions, non_ref_strs_positions) if s]

        repeats_refs = [chr_dict[nums[0]].split('-') for nums in subsets if nums]

        annotations = [
            f"{lst[0]}-{int(lst[1]) + 1}-{lst[3]}-{len(repeat.group()) // len(lst[3].split(return_digit(lst[3]).strip('--'))[-1])}{repeat.group()[0:len(lst[3].split(return_digit(lst[3]).strip('--'))[-1])]}"
            for lst, repeat in zip(repeats_refs, repeats)]
        return [(annotation, tag, lst[1], seq, tag_typ) for annotation, lst in zip(annotations, repeats_refs)]


def process_region(region, chromosome_STRdict, PeBamPath, min_seq_length, dcs_tags):
    current_nums = set(chromosome_STRdict[region])
    current_chr_dic = chromosome_STRdict[region]
    flat_strs_local = []
    for read in read_generator(PeBamPath, min_seq_length, region, dcs_tags):
        result = process_read(read_extractor(read, current_nums, current_chr_dic))
        if result:
            flat_strs_local.extend(result)
    return flat_strs_local

def read_extractor(read, current_nums, current_chr_dic):
    # current_nums = set(chromosome_STRdict[read.reference_name])
    # current_chr_dic = chromosome_STRdict[read.reference_name]
    return [read.query_name.split('.')[0], read.query_name.split('.')[1:], read.reference_name,
            tuple(read.get_reference_positions(True)), read.seq, current_nums, current_chr_dic]


def read_generator(bam_file_path: str, min_seq_length: int, region: Optional[str] = None, dcs_tags={}):
    with pysam.AlignmentFile(bam_file_path, 'rb') as bamfile:
        if region:
            fetch_generator = bamfile.fetch(region)
        else:
            fetch_generator = bamfile.fetch(until_eof=True)
        for read in fetch_generator:
            if read.seq and read.reference_name and len(read.seq) > min_seq_length and read.qname.split('.')[0] in dcs_tags:
                yield read


def len_read_types(reads, bam_type=None):
    if bam_type == 'dcs':
        reads_dict: Dict[str, Dict[str, int]] = {'dcs': {}}
    else:
        reads_dict: Dict[str, Dict[str, int]] = {'ab.1': {}, 'ab.2': {}, 'ba.1': {}, 'ba.2': {}}
    for read in reads:
        if bam_type != 'dcs':
            tag, typ = read.qname.split('.', 1)
        else:
            tag, typ = read.qname, 'dcs'
        # If this part is one of the keys in reads_dict, increment its count
        if reads_dict.get(typ) is not None:
            reads_dict[typ][tag] = reads_dict[typ].get(tag, 0) + 1

    return reads_dict
