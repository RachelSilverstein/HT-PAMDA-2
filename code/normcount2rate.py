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
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.stats import linregress
from scipy.stats import pearsonr
from scipy.stats import skew
from tqdm.autonotebook import tqdm

from inputs import *
from functions import *

def normcount2rate(run_name,
                   pam_length,
                   pam_start,
                   timepoints,
                   init_rate_est=[0.0001, 0.001, 0.01],
                   read_sum_minimum=10,
                   tps_sum_minimum=2,
                   use_timepoints=None,
                   input_csv=None):
    """
    calculate rates (k) of PAM depletion, modeled as exponential decay y(t)=e^-(k*t)
    """

    if input_csv == None:
        df = pd.read_csv('output/%s/PAM_start_%s_length_%s/PAMDA_2_norm_counts.csv' %
                         (run_name, pam_start, pam_length))
    else:
        df = pd.read_csv(input_csv)

    if use_timepoints is None:
        use_timepoints = [n for n in range(len(timepoints))]

    # exponential decay
    def func(x, a, b):
        return a * np.exp(-b * x)

    ks = []
    timepoint_indices = []
    for i in use_timepoints:
        timepoint_indices.append([i, timepoints[i]])

    # Find rates:
    print('calculating rate constants')

    # loop through one row at a time and nonlinear fit rates

    pbar = tqdm(desc='samples: ', total=df['Sample'].nunique(),
                file=sys.stdout)
    previous_row = 'n/a'
    for index, row in df.iterrows():

        current_row = row['Sample']
        if current_row != previous_row:
            pbar.update()

        tps = [x[1] for x in timepoint_indices]
        obs_raw = [0.00000000001] + [row['Raw_Counts_' + str(x[0])]
                                     for x in timepoint_indices if str(x[0]) != '0']
        obs_norm = [row['Norm_Counts_' + str(x[0])] for x in timepoint_indices]

        zero_indices = [i for i, e in enumerate(obs_raw) if e == 0]
        for zero_index in sorted(zero_indices)[::-1]:
            del tps[zero_index]
            del obs_norm[zero_index]
            del obs_raw[zero_index]

        if sum(obs_raw) >= read_sum_minimum and len(obs_norm) >= tps_sum_minimum:

            # map initial conditions to fit error
            min_search = []

            # loop across different initial conditions for k
            for j in init_rate_est:
                p0 = [1.0, j]
                # nonlinear curve fit
                try:
                    popt, pcov = curve_fit(func, tps, obs_norm, p0=p0)
                    # prediction from outputted paramters
                    pred = [func(x, popt[0], popt[1]) for x in tps]
                    # error
                    error = sum([(x[0] - x[1]) ** 2 for x in zip(pred, obs_norm)])
                    # append into dictionary mapping initial conditions to fit error
                    min_search.append([error, list(popt)])
                except:
                    continue
            if len(min_search) != 0:
                ks.append(sorted(min_search)[0][1][1])
            else:
                ks.append('NaN')

        else:
            ks.append('NaN')
        previous_row = current_row
    pbar.close()

    print('appending rate constants')

    min_k = min([100 if ((isinstance(x, float) and x <= 0) or x == 'NaN') else x for x in ks])

    df['Rate_Constant_k'] = [x if ((isinstance(x, float) and x > 0) or x == 'NaN') else min_k for x in ks]

    print('output to CSV')

    if not os.path.exists('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length)):
        os.makedirs('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length))

    df.to_csv('output/%s/PAM_start_%s_length_%s/PAMDA_3_rates.csv' %
              (run_name, pam_start, pam_length), index=False)

#-----------------------------------------------------------------------------------------------------------------------------#

if __name__ == "__main__":

    os.chdir('../')
    sys.path.append('code')

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

    # calculate rate constants from normalized read counts
    normcount2rate(RUN_NAME,
                   PAM_LENGTH,
                   PAM_START,
                   TIMEPOINTS,
                   INIT_RATE_EST,
                   READ_SUM_MIN,
                   TPS_SUM_MIN,
                   USE_TIMEPOINTS)
