from collections import defaultdict, Counter
from functools import lru_cache
import matplotlib.colors as mcolors
import numpy as np
import math
import pandas as pd

## important note: this list does not contain the mirror values (3,2911 3,9211)
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

tier_cat = defaultdict(list)  # tier classes counts 1:225 ..
for tier in ordered_tiers:
    tier_cat[int(str(tier)[0])].append(tier)
tier_cat = {i: len(tier_cat[i]) for i in tier_cat}
#bkup
# new_tier_scale = [
#     *np.linspace(1, 0.84, num=tier_cat[1]),
#     *np.linspace(0.84, 0.68, num=tier_cat[2]),
#     *np.linspace(0.68, 0.52, num=tier_cat[3]),
#     *np.linspace(0.52, 0.27, num=tier_cat[4]),
#     *np.linspace(0.27, 0.02, num=tier_cat[5])
# ]
# Generate the scaled values for each category
new_tier_scale = [
    *np.linspace(1, 0.84, num=tier_cat[1]),
    *np.linspace(0.84, 0.68, num=tier_cat[2]),
    *np.linspace(0.68, 0.60, num=tier_cat[3]),
    *np.linspace(0.6, 0.30, num=tier_cat[4]),
    *np.linspace(0.30, 0.02, num=tier_cat[5])
]

new_tier_color_scale = [
    *np.linspace(1, 0, num=tier_cat[1]),
    *np.linspace(1, 0, num=tier_cat[2]),
    *np.linspace(1, 0, num=tier_cat[3]),
    *np.linspace(1, 0, num=tier_cat[4]),
    *np.linspace(1, 0, num=tier_cat[5])
]

# tier_score = {tier: score for tier, score in zip(ordered_tiers, new_tier_scale)}  # tiers scores
# Parameters for exponential decay
lambda_decay = 0.00038  # Adjust this to control the decay rate
num_elements = len(ordered_tiers)

# Create a dictionary to store scores for each unique number
tier_score = {}

# Assign scores using exponential decay
for index, number in enumerate(ordered_tiers):
    # Calculate the score
    score = math.exp(-lambda_decay * index)
    # Round the score to avoid floating-point precision issues
    score = round(score, 5)  # Adjust precision as needed
    # Store the score in the dictionary
    tier_score[number] = score

tier_color_score = {tier: score for tier, score in zip(ordered_tiers, new_tier_color_scale)}  # tiers scores
tiers_order = {k: i for i, k in enumerate(ordered_tiers)}
tier_description = {
    1: 'DCS-MV2',
    2: 'DCS-MV1',
    3: 'DCS-MV0',
    4: 'SSCS-MV1',
    5: 'SSCS-MV0',
}
# Define base colors for each main group
base_colors = {
    1: '#006400',  # Dark Green for Group 1
    2: '#85EE1D',  # Light Green for Group 2
    3: '#B8EC84',  # Yellow for Group 3
    4: '#FFEA22',  # Orange for Group 4
    5: '#FF8C00',
}


@lru_cache(maxsize=None)
def get_tier_score(tier_value):
    if len(str(tier_value)) == 5 and str(tier_value)[0] in ['4', '5']:
        tier_value = str(tier_value)
        tier_value = float(f"{tier_value[0]}.{tier_value[2]}{tier_value[3]}0{tier_value[4]}")
    try:
        return tier_score[tier_value]
    except KeyError:
        tier_value = str(tier_value)
        tier_value0 = float(f"{tier_value[0]}.{tier_value[3]}{tier_value[2]}{tier_value[4]}{tier_value[5]}")
        try:
            return tier_score[tier_value0]
        except KeyError:
            tier_value = str(tier_value)
            tier_value1 = float(f"{tier_value[0]}.{tier_value[2]}{tier_value[3]}{tier_value[5]}{tier_value[4]}")
            try:
                return tier_score[tier_value1]
            except KeyError:
                tier_value = str(tier_value)
                tier_value2 = float(f"{tier_value[0]}.{tier_value[3]}{tier_value[2]}{tier_value[5]}{tier_value[4]}")
                try:
                    return tier_score[tier_value2]
                except KeyError:
                    return 'N/A'


def get_tier_color_score(tier_value):
    if len(str(tier_value)) == 5 and str(tier_value)[0] in ['4', '5']:
        tier_value = str(tier_value)
        tier_value = float(f"{tier_value[0]}.{tier_value[2]}{tier_value[3]}0{tier_value[4]}")
    try:
        return tier_color_score[tier_value]
    except KeyError:
        tier_value = str(tier_value)
        tier_value0 = float(f"{tier_value[0]}.{tier_value[3]}{tier_value[2]}{tier_value[4]}{tier_value[5]}")
        try:
            return tier_color_score[tier_value0]
        except KeyError:
            tier_value = str(tier_value)
            tier_value1 = float(f"{tier_value[0]}.{tier_value[2]}{tier_value[3]}{tier_value[5]}{tier_value[4]}")
            try:
                return tier_color_score[tier_value1]
            except KeyError:
                tier_value = str(tier_value)
                tier_value2 = float(f"{tier_value[0]}.{tier_value[3]}{tier_value[2]}{tier_value[5]}{tier_value[4]}")
                try:
                    return tier_color_score[tier_value2]
                except KeyError:
                    return 'N/A'


def get_tier_order(tier_value):
    if len(str(tier_value)) == 5 and str(tier_value)[0] in ['4', '5']:
        tier_value = str(tier_value)
        tier_value = float(f"{tier_value[0]}.{tier_value[2]}{tier_value[3]}0{tier_value[4]}")
    try:
        return tiers_order[tier_value]
    except KeyError:
        tier_value = str(tier_value)
        tier_value0 = float(f"{tier_value[0]}.{tier_value[3]}{tier_value[2]}{tier_value[4]}{tier_value[5]}")
        try:
            return tiers_order[tier_value0]
        except KeyError:
            tier_value = str(tier_value)
            tier_value1 = float(f"{tier_value[0]}.{tier_value[2]}{tier_value[3]}{tier_value[5]}{tier_value[4]}")
            try:
                return tiers_order[tier_value1]
            except KeyError:
                tier_value = str(tier_value)
                tier_value2 = float(f"{tier_value[0]}.{tier_value[3]}{tier_value[2]}{tier_value[5]}{tier_value[4]}")
                try:
                    return tiers_order[tier_value2]
                except KeyError:
                    return 'N/A'


# %%


# Function to calculate the color adjustment
def get_adjusted_color(group, tier_value):
    frequency = get_tier_color_score(tier_value)
    # Get the base color for the main group
    base_color = base_colors.get(group)
    if not base_color:
        raise ValueError("Group not recognized")

    # Convert base color to RGB
    base_rgb = mcolors.hex2color(base_color)

    # Calculate brightness factor based on frequency (0 = brightest, 1 = base color)
    brightness_factor = 1 - frequency  # Lower frequency = higher brightness adjustment

    # Apply the brightness adjustment to each RGB channel
    adjusted_rgb = [
        min(1, channel + brightness_factor * (1 - channel))  # Increase brightness proportionally
        for channel in base_rgb
    ]

    # Convert back to hex color
    return mcolors.to_hex(adjusted_rgb)


def tier_assigner(counts, frequencies):
    def find_tier_pct(num):
        if 0.9 <= num <= 1:
            return 9
        elif 0.8 <= num < 0.9:
            return 8
        elif 0.6 <= num < 0.8:
            return 6
        elif 0.4 <= num < 0.6:
            return 4
        elif 0.2 <= num < 0.4:
            return 2
        elif 0.1 <= num < 0.2:
            return 1
        else:
            return 0

    def fs_reduction(true_fs):
        return int(true_fs) if true_fs < 5 else 5

    sscs1 = True if counts[0] > 0 and counts[2] > 0 else False
    sscs2 = True if counts[1] > 0 and counts[3] > 0 else False
    fab = find_tier_pct(max([frequencies[0], frequencies[2]]))
    fba = find_tier_pct(max([frequencies[1], frequencies[3]]))
    fs1 = max([counts[0], counts[2]])
    fs2 = max([counts[1], counts[3]])
    tier_res = ''
    if fab == fba == 0:
        return 5.1001
    elif sscs1 and sscs2 and fab > 0 and fba > 0:
        tier_res += '1.'  # both strands both mates
    elif (sscs1 and sscs2) and ((fab == fba == 0) or (fab > 0 and fba == 0) or (fab == 0 and fba > 0)):
        tier_res += '4.'
    elif ((sscs1 and not sscs2) or (sscs2 and not sscs1)) and (fab > 0 and fba > 0):
        tier_res += '2.'  # one mate and one strand
    elif (not sscs1 and not sscs2) and (fab > 0 and fba > 0):
        tier_res += '3.'  # one mate
    elif ((sscs1 and not sscs2) or (not sscs1 and sscs2)) and (fab > 0 or fba > 0):
        tier_res += '4.'  # one strand
    elif (not sscs1 and not sscs2) and ((fab > 0 and fba == 0) or (fab == 0 and fba > 0)):
        tier_res += '5.'  # half strand  / half mate
    if tier_res in ['4.', '5.']:
        fab = max(fab, fba)
        fba = 0
        fs2 = max(fs1, fs2)
        fs1 = 0
        latter_info = f"{fab}{fba}{fs_reduction(fs1)}{fs_reduction(fs2)}"
    else:
        latter_info = f"{fab}{fba}{fs_reduction(fs1)}{fs_reduction(fs2)}"
    return float(tier_res + latter_info)


def tier_assigner_original(counts, frequencies):  # original
    ab = max([counts[0], counts[2]])
    ba = max([counts[1], counts[3]])
    fab = max([frequencies[0], frequencies[2]])
    fba = max([frequencies[1], frequencies[3]])

    if ab >= 3 and ba >= 3 and fab == fba == 1:
        tier = 1.1
    elif (ab >= 3 and ba >= 1 and fab == fba == 1) or (ab >= 1 and ba >= 3 and fab == fba == 1):
        tier = 1.2
    elif ab >= 2 and ba >= 2 and fab == fba == 1:
        tier = 1.3

    elif ab >= 3 and ba >= 3 and fab >= 0.6 and fba >= 0.6:
        tier = 2.1
    elif (ab >= 3 and ba >= 1 and fab >= 0.6 and fba >= 0.6) or (ab >= 1 and ba >= 3 and fab >= 0.6 and fba >= 0.6):
        tier = 2.2
    elif (ab >= 3 and ba >= 3) and (fab >= 0.6 or fba >= 0.6):
        tier = 2.3
    elif ((ab >= 3 and ba >= 1) or (ab >= 1 and ba >= 3)) and (fab >= 0.6 or fba >= 0.6):
        tier = 2.4

    elif ab >= 3 and ba >= 3 and fab >= 0.5 and fba >= 0.5:
        tier = 3.1
    elif (ab >= 3 and ba >= 1 and fab >= 0.5 and fba >= 0.5) or (ab >= 1 and ba >= 3 and fab >= 0.5 and fba >= 0.5):
        tier = 3.2
    elif (ab >= 1 and ba >= 2 and fab >= 0.5 and fba >= 0.5) or (ab >= 2 and ba >= 1 and fab >= 0.5 and fba >= 0.5):
        tier = 3.3
    elif (ab == 1 and ba == 1) and (fab == fba == 1):
        tier = 3.4
    elif (ab >= 3 and ba >= 3) and (fab >= 0.5 or fba >= 0.5):
        tier = 3.5
    elif ((ab >= 3 and ba >= 1) or (ab >= 1 and ba >= 3)) and (fab >= 0.5 or fba >= 0.5):
        tier = 3.6
    elif ((ab >= 2 and ba >= 1) or (ab >= 1 and ba >= 2)) and (fab >= 0.5 or fba >= 0.5):
        tier = 3.7

    elif (ab >= 5 and ba == 0 and fab == 1) or (ba >= 5 and ab == 0 and fba == 1):
        tier = 4.1
    elif (ab >= 5 and ba == 0 and fab >= 0.8) or (ba >= 5 and ab == 0 and fba >= 0.8):
        tier = 4.2
    elif (ab >= 5 and ba == 0 and fab >= 0.6) or (ba >= 5 and ab == 0 and fba >= 0.6):
        tier = 4.3

    elif (ab >= 3 and ba == 0 and fab == 1) or (ba >= 3 and ab == 0 and fba == 1):
        tier = 5.1
    elif (ab >= 3 and ba == 0 and fab >= 0.8) or (ba >= 3 and ab == 0 and fba >= 0.8):
        tier = 5.2
    elif (ab >= 3 and ba == 0 and fab >= 0.6) or (ba >= 3 and ab == 0 and fba >= 0.6):
        tier = 5.3
    elif (ab >= 2 and ba == 0 and fab >= 1) or (ba >= 2 and ab == 0 and fba >= 1):
        tier = 5.4
    elif (ab >= 1 and ba == 0 and fab >= 1) or (ba >= 1 and ab == 0 and fba >= 1):
        tier = 5.5
    elif (ab >= 2 and ba == 0 and fab >= 0.5) or (ba >= 2 and ab == 0 and fba >= 0.5):
        tier = 5.6
    elif (ab >= 2 and ba == 0 and fab >= 0.3) or (ba >= 2 and ab == 0 and fba >= 0.3):
        tier = 5.7
    else:
        tier = 7.0

    return tier


def compress_strs2(str_file_path, name=None, res_dir=None):
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
        main_tiers_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
        for t in tier_counts:
            main_tiers_counts[int(t)] += tier_counts[t]

        try:
            rows_info.append(
                [repeat_id, ('REF' if ref == alt else 'ALT'), cvg, AC, AF] + [main_tiers_counts[tier] for tier in
                                                                              main_tiers_counts] + [pos,
                                                                                                    int(
                                                                                                        region.split(
                                                                                                            'chr')[
                                                                                                            -1])])
        except ValueError:
            rows_info.append(
                [repeat_id, ('REF' if ref == alt else 'ALT'), cvg, AC, AF] + [main_tiers_counts[tier] for tier in
                                                                              main_tiers_counts] + [pos,
                                                                                                    int(1)])

    compressed_df = pd.DataFrame(rows_info,
                                 columns=['STR ID', 'type', 'cvrg', 'AC', 'AF'] + list(main_tiers_counts) + [
                                     'start_pos',
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
        tier_color_score = {i: get_adjusted_color(int(i), (float(f"{i}.9955") if i <4 else float(f"{i}.9005"))) for i in list(range(1,6))}
        for col_name, color in tier_color_score.items():
            if col_name in compressed_df.columns:
                col_letter = compressed_df.columns.get_loc(col_name)  # + 1
                worksheet.conditional_format(1, col_letter, compressed_df.shape[0], col_letter, {
                    'type': 'cell',
                    'criteria': '>',
                    'value': 0,
                    'format': workbook.add_format({'bg_color': color}),
                })
    return compressed_df
