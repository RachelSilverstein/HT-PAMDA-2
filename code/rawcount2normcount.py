from inputs import *
from functions import *

import os
import sys
import glob
import gzip
import itertools
import os
import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import skew
from tqdm.autonotebook import tqdm



def rawcount2normcount(run_name,
                       control_rawcount_csv,
                       control_sample,
                       control_sample_timepoint_fastq,
                       pam_orientation,
                       pam_length,
                       pam_start,
                       spacers={'SPACER1': 'GGGCACGGGCAGCTTGCCGG', 'SPACER2': 'GTCGCCCTCGAACTTCACCT'},
                       control_spacers=None,
                       control_spacer_mapping=None,
                       timepoints=[0, 60, 480, 1920],
                       max_pam_length=8,
                       top_n=5,
                       input_csv=None):
    """
    generate normalized PAM read counts from raw counts
    """
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    nucleotides = ['A', 'T', 'C', 'G']
    total_pam_space = [''.join(p) for p in itertools.product(nucleotides, repeat=pam_length)]

    if input_csv is None:
        df_input = pd.read_csv('output/%s/PAMDA_1_raw_counts.csv.gz' % (run_name))
    else:
        df_input = pd.read_csv(input_csv)

    if control_rawcount_csv is not None:
        df_control = pd.read_csv(control_rawcount_csv)
        df_input = pd.concat([df_input, df_control], sort=False)
        control_sample_timepoint_fastq = 1

    df_input.reset_index(inplace=True, drop=True)

    # reformat the input dataframe so we count the uncleaved control spacers as extra "PAMs"
    if control_spacers is not None:
        df_input["PAM"] = df_input["PAM"].where(df_input["PAM"] != 'fixed', other=df_input["Spacer"])
        for i in range(len(df_input)):
            spacer = df_input.loc[i, "Spacer"]
            for main_spacer in control_spacer_mapping:
                if spacer in control_spacer_mapping[main_spacer]:
                    df_input.loc[i, "Spacer"] = main_spacer


    # group the long PAMs into shorter PAMs
    df_input = group_pams(df_input, pam_start, pam_length, max_pam_length, pam_orientation, control_spacers)


    # add t0 read counts to each sample from control sample
    df_input = add_t0_raw_counts(df_input, control_sample, control_sample_timepoint_fastq)

    # add fractional counts for each sample/spacer combo
    df_input = add_fractional_counts(df_input, timepoints)

    # choose control PAMs for each sample (most enriched ones or designated control spacers)
    # correct for increased counts in later timepoints
    if top_n > 0: # do not perform uptrend corrections if the top n is set to 0
        df_input = correct_uptrends(df_input, control_spacers, timepoints, run_name, pam_start, pam_length, top_n)

    # normalizing read counts
    df_input = norm_to_t0_abundance(df_input, timepoints)

    if not os.path.exists('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length)):
        os.makedirs('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length))

    df_input.to_csv('output/%s/PAM_start_%s_length_%s/PAMDA_2_norm_counts.csv' %
                    (run_name, pam_start, pam_length), index=False)

#-----------------------------------------------------------------------------------------------------------------------------#

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

    # calculate normalized read counts from raw read counts
    rawcount2normcount(RUN_NAME,
                       CONTROL_RAW_COUNT_CSV,
                       CONTROL_SAMPLE,
                       CONTROL_SAMPLE_TIMEPOINT_FASTQ,
                       PAM_ORIENTATION,
                       PAM_LENGTH,
                       PAM_START,
                       SPACERS,
                       CONTROL_SPACERS,
                       CONTROL_SPACER_MAPPING,
                       TIMEPOINTS,
                       MAX_PAM_LENGTH,
                       TOP_N_NORMALIZE)

