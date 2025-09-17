import glob
import json
import os
import pickle
import re
from collections import Counter
import PyPDF2
import pandas as pd
import pysam
from pyfasta import Fasta
from Bio.Seq import Seq


def calculate_new_alt(ref, alt, is_deletion=True):
    # Extract the repeat count and sequence from ref
    ref_match = re.match(r'(\d+)([A-Z]+)', ref)
    if not ref_match:
        raise ValueError("Invalid reference format")

    ref_count, ref_sequence = ref_match.groups()
    ref_count = int(ref_count)

    # Calculate the number of times ref_sequence occurs in alt
    alt_sequence_length = len(ref_sequence)

    # Sliding window approach to count occurrences and identify additional bases
    alt_count = 0
    additional_bases = []
    i = 0

    while i <= len(alt) - alt_sequence_length:
        if alt[i:i + alt_sequence_length] == ref_sequence:
            alt_count += 1
            i += alt_sequence_length
        else:
            additional_bases.append(alt[i])
            i += 1

    # Include any remaining bases after the last repeat
    additional_bases.extend(alt[i:])

    # Calculate the difference in repeat counts
    if is_deletion:
        new_alt_count = ref_count - alt_count
    else:
        new_alt_count = ref_count + alt_count

    # Count additional bases
    additional_counts = Counter(additional_bases)
    additional_str = ''.join(f"{count}{base}" for base, count in additional_counts.items())

    # Generate the new alt in the required format
    if new_alt_count == 0:
        new_alt = additional_str or ref_sequence
    elif new_alt_count > 0:
        new_alt = f"{new_alt_count}{ref_sequence}" + additional_str
    else:
        new_alt = f"{abs(new_alt_count)}{ref_sequence}" + additional_str

    return new_alt


# %%
def sum_numbers_in_string(s):
    #This script sums all integers in a string, e.g., '16A1' --> 17
    # Use regular expression to find all numbers in the string
    numbers = re.findall(r'\d+', s)
    # Convert the numbers to integers and sum them
    total_sum = sum(int(num) for num in numbers)
    return total_sum


# %%
def find_clusters(variants, window_size, min_hits):
    # This script finds clustered variants ('chr16-123-G-T') it looks if within a given window size, n number of
    # variants exits [('chr16-123-G-T'),('chr16-124-G-T'),('chr16-125-G-T')('chr16-126-G-T')].
    # Extract chromosome and positions from the variant strings
    variant_data = [(variant, variant.split('-')[0], int(variant.split('-')[1])) for variant in variants]
    # Sort the variants first by chromosome, then by position
    variant_data.sort(key=lambda x: (x[1], x[2]))

    clusters = []
    n = len(variant_data)

    for i in range(n):
        cluster = [variant_data[i][0]]
        for j in range(i + 1, n):
            if variant_data[j][1] == variant_data[i][1] and variant_data[j][2] - variant_data[i][2] <= window_size:
                cluster.append(variant_data[j][0])
            else:
                break
        if len(cluster) >= min_hits:
            clusters.append(cluster)

    return clusters


# %%
def split_coveragefile_by_chromosome(chr_regions: list, cvg_path: str):
    out_path = '/'.join(cvg_path.split('/')[:-1])
    lib = cvg_path.split('/')[-2]
    cvg_file = pd.read_csv(cvg_path, sep='\t', header=None, names=['chr', 'start', 'count'])
    groups = cvg_file.groupby(['chr'])
    for chrom, chr_df in groups:
        if chrom[0] in chr_regions:
            lst = [f"{start}\t{cvg}\n" for start, cvg in zip(chr_df['start'], chr_df['count'])]
            with open(f"{out_path}/{lib}_{chrom[0]}_DCS_coverage.tabular", 'w') as outfile:
                outfile.writelines(lst)


# %%
def bamtofasta(path, name):
    reads = [read for read in
             pysam.AlignmentFile(path).fetch(
                 until_eof=True)]

    reads_dict = {read.qname: [] for read in reads}
    for read in reads:
        reads_dict[read.qname].append(read.seq)

    to = '/'.join(path.split('/')[:-1]) + '/' + name + '_bamto.txt'
    with open(to, "w+") as f:
        for read in reads_dict:
            sequences = reads_dict[read]
            for num, seq in enumerate(sequences):
                f.write(f">{read}.{num + 1}\n")
                f.write(f"{seq}\n")


# %%

def excel2pkl(xlsx_file_path, lib, pckl_path, snp, strr):
    if strr:
        file_pkl_path = f"{pckl_path}/{lib}_STRSummaryCash.pkl"
    else:
        if snp:
            file_pkl_path = f"{pckl_path}/{lib}_SNPSummaryCash.pkl"
        else:
            file_pkl_path = f"{pckl_path}/{lib}_SummaryCash.pkl"

    summary_file_sheets = pd.ExcelFile(xlsx_file_path).sheet_names
    print(summary_file_sheets)
    temp_summary_dfs = []
    for sheet in summary_file_sheets:
        print(sheet)
        print()
        temp_summary = pd.read_excel(xlsx_file_path, sheet_name=sheet)
        temp_summary_dfs.append(temp_summary)
    summary_file = pd.concat(temp_summary_dfs)
    summary_file.to_pickle(file_pkl_path,compression="gzip")

    print('DONE!')


# %%
def generate_html_with_buttons(plot_files, output_file):
    with open(output_file, 'w') as f:
        f.write(
            '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="UTF-8">\n<title>Plot Viewer</title>\n</head>\n<body>\n')
        for i, plot_file in enumerate(plot_files):
            button_id = f'button_{i}'

            # add inline style for the button
            f.write(
                f'<button id="{button_id}" onclick="openPlot(\'{plot_file}\')" style="font-size: 30px; padding: 10px;">{plot_file.split("/")[-1].split("_")[-1][:-5]} </button>\n')

        f.write('<script>\nfunction openPlot(plotFile) {\nwindow.open(plotFile);\n}\n</script>\n')
        f.write('</body>\n</html>')


# %%
def seqNucleotides(galaxy_folder, regions_path, wf_date):
    main_analysis_folder = '/Users/shehabmoukbel/Desktop/11_Bioinformatics Analysis/0_LibrariesAnalysis'
    current_analysis_folder = f"{main_analysis_folder}/{wf_date}_Analysis"
    galaxy_name = galaxy_folder.split('/')[-2]
    libs = [i for i in sorted(glob.glob(galaxy_folder))]
    dcs_bam_paths = sorted([bam for i in libs for bam in glob.glob(i + '/*') if
                            'DCS.bam' in bam and 'bai' not in bam and 'extracted' not in bam and 'sai' not in bam])
    region_df = pd.read_excel(regions_path)[:-2].fillna('')
    regions_dic = {}
    region_df_group = region_df.groupby(['chr', 'type'])
    for (chro, type), temp_df in region_df_group:
        regions_dic[f"{type}"] = (int(list(temp_df['start'])[0]), int(list(temp_df['end'])[0]) + 1)

    libs_n = {lib.split('/')[-1]: 0 for lib in libs}

    for Dcs_BamPath in dcs_bam_paths:
        lib = Dcs_BamPath.split('/')[-2]

        reads_lengths = {region: 0 for region in regions_dic}
        reads_lengths['cassette'] = 0
        for tpl in regions_dic.items():
            chromosome, (start, end) = tpl
            temp_chr = chromosome
            if 'chr16' in chromosome:
                temp_chr = 'chr16'
            for read in pysam.AlignmentFile(Dcs_BamPath).fetch(temp_chr, start - 1, end + 1):
                reads_lengths[chromosome] += read.query_alignment_length

            for read in pysam.AlignmentFile(Dcs_BamPath).fetch('chr16', 695270, 695531):
                reads_lengths['cassette'] += read.query_alignment_length

        libs_n[lib] = reads_lengths
    # all_n = list(libs_n.values()) + [sum(list(libs_n.values()))]

    res = pd.DataFrame(index=list(libs_n) + ['Sequenced nucleotides'],
                       columns=['chr3-SeqBp', 'chr15-SeqBp', 'chr16PS-SeqBp', 'chr16FK-SeqBp', 'cassette-SeqBp',
                                'total'])
    for lib in libs_n:
        res.loc[lib, 'chr3-SeqBp'] = libs_n[lib]['chr3']
        res.loc[lib, 'chr15-SeqBp'] = libs_n[lib]['chr15']
        res.loc[lib, 'chr16FK-SeqBp'] = libs_n[lib]['chr16FK']
        res.loc[lib, 'chr16PS-SeqBp'] = libs_n[lib]['chr16PS']
        res.loc[lib, 'cassette-SeqBp'] = libs_n[lib]['cassette']
        res.loc[lib, 'total'] = sum(list(libs_n[lib].values())[:-1])

    res.loc['Sequenced nucleotides', 'total'] = sum(res['total'][:-1])

    res.to_excel(f'{current_analysis_folder}/{galaxy_name}_NB_SequencedBp.xlsx', sheet_name='entire_bam 5')
    print('Done!')


# %%
def var_region_finder(regions_file_path: str, variants_list: list, version: int):
    regions_file = pd.read_excel(regions_file_path, sheet_name='Probe regions')[:-2]
    region_name_not_int = False
    if version == 3:
        present_combinations = [(int(i), k, s, e) for i, k, s, e in
                                zip(regions_file['chr'], regions_file['region'], regions_file['start'],
                                    regions_file['end'])]
    else:
        try:
            present_combinations = [(int(i), k, s, e) for i, k, s, e in
                                    zip(regions_file['chr'], regions_file['type'], regions_file['start'],
                                        regions_file['end'])]
        except ValueError:
            present_combinations = [(i, k, s, e) for i, k, s, e in
                                    zip(regions_file['chr'], regions_file['type'], regions_file['start'],
                                        regions_file['end'])]
            region_name_not_int = True




    regions_dic = {int(i[0]): {} for i in present_combinations} if not region_name_not_int else {i[0]: {} for i in present_combinations}
    for tpl in present_combinations:
        regions_dic[tpl[0]][tpl[1]] = range(int(tpl[2]), int(tpl[3] + 1))

    regions_res = []
    for var in variants_list:
        temp_lst = var.split('-')[:-2]
        chrom, start = int(temp_lst[0].split('chr')[1]) if not region_name_not_int else temp_lst[0] , int(temp_lst[1])
        found = False
        if chrom not in regions_dic:
            if version == 3:
                regions_res.append(f"OffTarget")
            else:
                regions_res.append(f"O{temp_lst[0]}")
        else:
            for possibility in regions_dic[chrom]:
                if start in regions_dic[chrom][possibility]:
                    if version == 1:
                        regions_res.append(possibility)
                    elif version == 2:
                        regions_res.append(f'T{possibility}')
                    elif version == 3:
                        regions_res.append(f'{possibility.split("-")[0]}')
                    found = True
                    break
            if not found:
                if version == 3:
                    regions_res.append(f"OffTarget")
                else:
                    regions_res.append(f"O{temp_lst[0]}")
    return regions_res


# %%
def get_key(val, dic):  # This function returns the key of a value of a dictionary
    for key, value in dic.items():
        if val in value:
            return key


# %%

def merge_pdfs(files: list, file_name):  # This function merges PDF files
    pdf_files = files

    # Create a PDF merger object
    pdf_merger = PyPDF2.PdfMerger()

    try:
        # Append each PDF to the merger object
        for pdf_file in pdf_files:
            pdf_merger.append(pdf_file)

        # Output file name for the merged PDF
        output_pdf = f"{file_name}.pdf"

        # Write the merged PDF to the output file
        with open(output_pdf, 'wb') as output_file:
            pdf_merger.write(output_file)

        print(f'Merged PDF saved as {output_pdf}')

    except Exception as e:
        print(f'Error: {str(e)}')

    finally:
        # Close the PDF merger object
        pdf_merger.close()


# %%

def rename_files_and_folders(path):
    # Iterate through the contents of the folder
    for item in os.listdir(path):
        item_path = os.path.join(path, item)

        # Check if it's a file or folder
        if os.path.isfile(item_path) or os.path.isdir(item_path):
            # Check if the item's name contains ":"
            if ":" in item:
                new_name = item.replace(":", "-")
                new_path = os.path.join(path, new_name)
                os.rename(item_path, new_path)
                print(f'Renamed: {item} to {new_name}')

            # If it's a folder, recursively process its contents
            if os.path.isdir(item_path):
                rename_files_and_folders(item_path)


# %%
def find_str_pos(pos_lst, ref_str_genome):  # extracts STR positions from a list of positions
    str_only = [pos for pos in pos_lst if any(pos in i for i in ref_str_genome)]
    return sorted(str_only)


# %%
def find_digit(string):  # returns True if any char of a string is a number
    for character in string:
        if character.isdigit():
            return True
    return False


# %%
def return_digit(string):  # this will return the numbers of a string and separate them by -
    digit = ''
    for character in string:
        if character.isdigit():
            digit += character
        else:
            digit += '-'
    if not all(i == '-' for i in digit):
        return digit
    else:
        return False


# %%
def return_chromosome(variants):
    return [i.split('-')[0] for i in variants]


# %%

def type_finder(lst, simple=True):
    if simple:
        types_lst = ['Str' if find_digit(i.split('-')[-2]) else 'Indel' if len(i.split('-')[-1]) > len(
            i.split('-')[-2]) else 'Indel' if len(i.split('-')[-1]) < len(i.split('-')[-2]) else 'pM' for i in lst]
    else:
        types_lst = ['Str' if find_digit(i.split('-')[-2]) else 'Ins' if len(i.split('-')[-1]) > len(
            i.split('-')[-2]) and not find_digit(i.split('-')[-2]) else 'Del' if len(i.split('-')[-1]) < len(
            i.split('-')[-2]) and not find_digit(i.split('-')[-2]) else 'pM' for i in lst]
    return types_lst


# %%
import re


def classify_strs(variants):
    classifications = []

    for variant in variants:
        # Extract the sequences from the variant
        match = re.search(r'(\d+[A-Z]+)-(\d+[A-Z]+)', variant)
        if match:
            seq1, seq2 = match.groups()

            # Extract the repeat counts and sequences
            count1, base_seq1 = re.match(r'(\d+)([A-Z]+)', seq1).groups()
            count2, base_seq2 = re.match(r'(\d+)([A-Z]+)', seq2).groups()

            # Convert repeat counts to integers
            count1 = int(count1)
            count2 = int(count2)

            if base_seq1 == base_seq2:
                if count1 < count2:
                    classifications.append('repeat-gain')
                elif count1 > count2:
                    classifications.append('repeat-loss')
                else:
                    classifications.append('no change')
            else:
                classifications.append('repeat-pM')

    return classifications


# %%

def find_genome_strs(fasta_path: str, str_l: int, present_chromosome=None, identifier=None):
    genome = Fasta(fasta_path)
    save_path = '/'.join(fasta_path.split('/')[:-1])
    regex = r"(?:A){5,}|(?:G){5,}|(?:C){5,}|(?:T){5,}|(?:g){5,}|(?:t){5,}|(?:c){5,}|(?:a){5,}"
    # contains each Chromosome as key and its sequence
    if present_chromosome:
        ref_info = {i[0].split(' ')[0]: {'seq': str(i[1])} for i in list(genome.items()) if
                    i[0].split(' ')[0] in present_chromosome}
    else:
        ref_info = {i[0].split(' ')[0]: {'seq': str(i[1])} for i in list(genome.items())}
    # contains each chromosome as key and the STR coordinates
    ref_str_genome = {i: [] for i in ref_info}

    # here you can define the minimum repeat size must match called STR parameters

    for chrom in ref_str_genome:
        seq = ref_info[chrom]['seq']
        matches = [repeat for repeat in list(re.finditer(regex, seq)) if len(repeat.group()) >= str_l]

        for match in matches:
            str_start = match.start()
            str_end = match.end()
            ref_str_genome[chrom].append(range(str_start, str_end + 1))

    res_name = f"{save_path}/GenomeSTRs_L{str_l}.pkl" if not identifier else f"{save_path}/{identifier}_GenomeSTRs_L{str_l}.pkl"
    with open(res_name, "wb") as f:
        pickle.dump(ref_str_genome, f)
    return f'Done!, Saved in {save_path}'


# %%


def find_tags_by_haplotype(excel_hap: str,
                           results):  # helps to find the tags of a haplotype for quick check for rows in Excel file
    hap = excel_hap.replace('\t', '')
    return {i: results[i] for i in results if hap == results[i]}


def extract_regions(input_bam, regions_lst, lib):  # exacts a subset bamfile of a given region
    def extract_region(input_samfile, output_bam, chromosome, start, end):
        output_samfile = pysam.AlignmentFile(output_bam, "wb", header=input_samfile.header)
        for read in input_samfile.fetch(chromosome, start, end):
            output_samfile.write(read)
        output_samfile.close()

    input_samfile = pysam.AlignmentFile(input_bam, "rb")

    for region in regions_lst:
        region_name = region[0]

        # Extract region and generate output BAM file based on region name
        output_bam_file = f"{lib}_{region[1]}.bam"
        chromosome, start_position, end_position = region_name.split('-')
        extract_region(input_samfile, output_bam_file, chromosome, int(start_position) - 1, int(end_position) + 1)
        pysam.index(output_bam_file)

    input_samfile.close()


# %%
def flatten_sequences(lst):
    result = []
    for sublist in lst:
        flattened_sublist = ''
        i = 0
        while i < len(sublist):
            if sublist[i].isdigit():
                # Multiply the next element with the number
                flattened_sublist += sublist[i + 1] * int(sublist[i])
                i += 2
            else:
                # Just append the element to the flattened_sublist
                flattened_sublist += sublist[i]
                i += 1
        result.append(flattened_sublist)
    return result


def flatten_sequence(seq):
    result = ''
    prev_is_number = False
    c = 0
    for char in seq:
        if char.isdigit():

            if not prev_is_number:
                result += f'-{char}'
            else:
                result += char
            prev_is_number = True
            c += 1
        else:
            if prev_is_number:
                result += f"-{char}"
                prev_is_number = False
            else:
                result += char

    return flatten_sequences([[result] if c == 0 else result[1:].split('-')])[0]


def adjust_seq_length(seq):
    remainder = len(seq) % 3
    if remainder:
        seq += Seq('N' * (3 - remainder))
    return seq


def get_protein_substitution(wild_type, mutant):
    wild_type = adjust_seq_length(wild_type)
    mutant = adjust_seq_length(mutant)
    wild_type_protein = wild_type.translate(to_stop=True)
    mutant_protein = mutant.translate(to_stop=True)

    for i, (wt_aa, mu_aa) in enumerate(zip(wild_type_protein, mutant_protein)):
        if wt_aa != mu_aa:
            return f"p.{wt_aa}{i + 1}{mu_aa}"
