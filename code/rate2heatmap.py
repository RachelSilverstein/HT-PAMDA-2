import os
import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import skew
from tqdm.autonotebook import tqdm

from inputs import *
from functions import *


def rate2heatmap(run_name,
                 barcode_csv,
                 pam_length,
                 pam_start,
                 spacers,
                 control_spacers,
                 pam1_nucleotide_rank={1: 'A', 2: 'C', 3: 'G', 4: 'T'},
                 pam2_nucleotide_rank={1: 'A', 2: 'C', 3: 'G', 4: 'T'},
                 pam1_index_rank=None,
                 pam2_index_rank=None,
                 avg_spacer=True,
                 heatmap_fixed_min=False,
                 heatmap_fixed_max=False,
                 log_scale_heatmap=True,
                 input_csv=None):
    """
    generate a heatmap representation of PAM preference
    """

    if input_csv == None:
        csv_input = 'output/%s/PAM_start_%s_length_%s/PAMDA_3_rates.csv' % (run_name, pam_start, pam_length)
    else:
        csv_input = input_csv

    plt.switch_backend('agg')

    variant_ids = pd.read_csv(barcode_csv)  # sample barcode file input
    variants = variant_ids['sample'].tolist()
    variant_name_dict = {}
    for index, row in variant_ids.iterrows():
        variant_name_dict[row['sample']] = row['description']

    if not os.path.exists('figures/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length)):
        os.makedirs('figures/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length))

    df_input = pd.read_csv(csv_input)
    # Remove control only spacers if used
    if control_spacers is not None:
        ind = ~df_input["PAM"].isin(control_spacers.keys())
        df_input = df_input.loc[ind, :]
        df_input.reset_index(inplace=True)

    # Set very low rate constants to heatmap min so they don't pull down the average
    df_input["Rate_Constant_k"] = np.maximum(df_input["Rate_Constant_k"], heatmap_fixed_min)

    if (pam1_index_rank == None or pam2_index_rank == None):
        # if not specified, split PAM in the middle for x and y-axis
        # default ordering is higher priority for "inner" bases of the PAM
        split_pam_index = np.floor_divide(pam_length, 2)
        pam1_index_rank = [x for x in range(0, split_pam_index)][::-1]
        pam2_index_rank = [x for x in range(split_pam_index, pam_length)]
    else:
        split_pam_index = len(pam1_index_rank)

    df_input['PAM_pt1'] = [x[:split_pam_index] for x in df_input['PAM'].tolist()]
    df_input['PAM_pt2'] = [x[split_pam_index:] for x in df_input['PAM'].tolist()]

    spacers = list(spacers.keys())

    pam_length = len(df_input['PAM'].tolist()[0])

    numbers = ['1', '2', '3', '4']
    pam_space = [''.join(p) for p in itertools.product(numbers, repeat=pam_length)]

    # Sort PAMs according to rules
    pam1_ids = []
    pam2_ids = []
    for pam in pam_space:
        pam1_ids.append(int(pam[:split_pam_index]))
        pam2_ids.append(int(pam[split_pam_index:]))

    pam1_ids = sorted([int(x) for x in set(pam1_ids)])
    pam2_ids = sorted([int(x) for x in set(pam2_ids)])

    columns = []
    indices = []
    tmp = ''
    for pam in pam1_ids:
        tmp_dict = {i: j for i, j in zip(pam1_index_rank, list(str(pam)))}
        for i in sorted(tmp_dict):
            tmp += pam1_nucleotide_rank[int(tmp_dict[i])]
        indices.append(tmp)
        tmp = ''

    tmp = ''
    for pam in pam2_ids:
        tmp_dict = {i: j for i, j in zip(pam2_index_rank, list(str(pam)))}
        for i in sorted(tmp_dict):
            tmp += pam2_nucleotide_rank[int(tmp_dict[i])]
        columns.append(tmp)
        tmp = ''

    # loop through variants and make heatmaps
    pbar = tqdm(desc='samples: ', total=df_input['Sample'].nunique(),
                file=sys.stdout)
    new_columns = []
    for variant in df_input['Sample'].unique().tolist():
        for column in columns:
            new_columns.append(str(column) + '-' + str(variant))
        df_output = pd.DataFrame(columns=new_columns, index=indices)

        if avg_spacer:
            for row in df_output.index:
                for column in df_output.columns:
                    pam2 = str(column.split('-')[0].strip('\n'))
                    sample = str(column.split('-')[1].strip('\n'))
                    rate_avg = 0
                    for spacer in sorted(spacers):
                        rate_avg += float(df_input['Rate_Constant_k'][(df_input['PAM_pt1'] == str(row)) &
                                                                      (df_input['PAM_pt2'] == pam2) &
                                                                      (df_input['Sample'] == sample) &
                                                                      (df_input['Spacer'] == spacer)].tolist()[0])
                    rate_avg = rate_avg / len(spacers)
                    if log_scale_heatmap:
                        df_output.loc[row, column] = np.log10(rate_avg)
                    else:
                        df_output.loc[row, column] = rate_avg
            df_output = df_output[df_output.columns].astype(float)
            if (pam_length <= 5) and (pam_length >= 2):
                plot_nice_heatmap(df_output, heatmap_fixed_min, heatmap_fixed_max, log_scale_heatmap, pam_length, split_pam_index, columns, indices, avg_spacer, run_name, pam_start, variant, variant_name_dict)
            else:
                plot_default_heatmap(df_output, avg_spacer, heatmap_fixed_min, heatmap_fixed_max, variant, variant_name_dict, run_name, pam_start, pam_length, pam2_nucleotide_rank, pam2_ids)
            if len(spacers) == 2:
                spacer_correlation(df_input, variant, spacers[0], spacers[1], run_name, pam_start, pam_length, variant, variant_name_dict)

        else:
            for spacer in sorted(spacers):
                for row in df_output.index:
                    for column in df_output.columns:
                        pam2 = str(column.split('-')[0].strip('\n'))
                        sample = str(column.split('-')[1].strip('\n'))
                        rate = float(df_input['Rate_Constant_k'][(df_input['PAM_pt1'] == str(row)) &
                                                                 (df_input['PAM_pt2'] == pam2) &
                                                                 (df_input['Sample'] == sample) &
                                                                 (df_input['Spacer'] == spacer)].tolist()[0])
                        if log_scale_heatmap:
                            df_output.loc[row, column] = np.log10(rate)
                        else:
                            df_output.loc[row, column] = rate
                df_output = df_output[df_output.columns].astype(float)
                if (pam_length <= 5) and (pam_length >= 2):
                    plot_nice_heatmap(df_output, heatmap_fixed_min, heatmap_fixed_max, log_scale_heatmap, pam_length, split_pam_index, columns, indices, avg_spacer, run_name, pam_start, variant, variant_name_dict)
                else:
                    plot_default_heatmap(df_output, avg_spacer, heatmap_fixed_min, heatmap_fixed_max, variant, variant_name_dict, run_name, pam_start, pam_length, pam2_nucleotide_rank, pam2_ids)
                if len(spacers) == 2:
                    spacer_correlation(df_input, variant, spacers[0], spacers[1], run_name, pam_start, pam_length, variant, variant_name_dict)
        pbar.update()
        new_columns = []
    pbar.close()



if __name__ == "__main__":

    os.chdir('../')
    sys.path.append('code')

    # check inputs
    check_inputs(RUN_NAME,
                 BARCODE_CSV,
                 FASTQ_DIR,
                 TIMEPOINT_FASTQ,
                 PAM_ORIENTATION,
                 PAM_LENGTH,
                 PAM_START,
                 CONTROL_RAW_COUNT_CSV,
                 CONTROL_SAMPLE,
                 CONTROL_SAMPLE_TIMEPOINT_FASTQ,
                 TIMEPOINTS,
                 MAX_PAM_LENGTH,
                 SPACERS,
                 CONTROL_SPACERS,
                 P5_SAMPLE_BARCODE_START,
                 P7_SAMPLE_BARCODE_START,
                 USE_TIMEPOINTS,
                 TOP_N_NORMALIZE,
                 INIT_RATE_EST,
                 READ_SUM_MIN,
                 TPS_SUM_MIN,
                 PAM1_NT_RANK,
                 PAM2_NT_RANK,
                 PAM1_INDEX_RANK,
                 PAM2_INDEX_RANK,
                 AVERAGE_SPACER,
                 HEATMAP_FIXED_MIN,
                 HEATMAP_FIXED_MAX,
                 LOG_SCALE_HEATMAP)

    # generate heatmap representations of PAM preference from rates
    print("CONTROL SPACERS in rate2heatmap:", CONTROL_SPACERS)  # delete later
    rate2heatmap(RUN_NAME,
                 BARCODE_CSV,
                 PAM_LENGTH,
                 PAM_START,
                 SPACERS,
                 CONTROL_SPACERS,
                 pam1_nucleotide_rank=PAM1_NT_RANK,
                 pam2_nucleotide_rank=PAM2_NT_RANK,
                 pam1_index_rank=PAM1_INDEX_RANK,
                 pam2_index_rank=PAM2_INDEX_RANK,
                 avg_spacer=AVERAGE_SPACER,
                 heatmap_fixed_min=HEATMAP_FIXED_MIN,
                 heatmap_fixed_max=HEATMAP_FIXED_MAX,
                 log_scale_heatmap=LOG_SCALE_HEATMAP
                 )


