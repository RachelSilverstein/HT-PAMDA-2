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


def reverse_complement(seq):
    nt_complement = dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})
    return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])


def find_BCs_and_spacer(spacers, control_spacers, read_sequenceR1, read_sequenceR2,
                        P5_sample_BC_start, P5_sample_BC_len,
                        P7_sample_BC_start, P7_sample_BC_len):
    """
    find the sample barcodes and the spacer, orientation, and location
    more complicated than necessary for standard HT-PAMDA library prep
    intended for compatibility with other library prep methods
    such as cases where P5 and P7 ends are not defined

    return:
        top_read: read containing spacer, oriented 5' to 3'. If the spacer is found in both reads,
            top read is the read where the spacer occurs earlier in the sequence
        bot_read: the other read
        spacer: the spacer that was found in the input 'spacers' list
        spacer_loc: the location of the 5' end of the spacer in the top_read
        P5_sample_BC: the P5 sample barcode
        P7_sample_BC: the P7 sample barcode
    """

    spacer_loc = -1  # location of 5' end of spacer
    spacer_loc_rc = -1
    spacer = None
    top_read = None
    bot_read = None
    P5_sample_BC = None
    P7_sample_BC = None

    if control_spacers is not None:
        spacers.update(control_spacers)

    for sp in spacers:
        spacer_loc = read_sequenceR1.find(spacers[sp])
        spacer_loc_rc = reverse_complement(read_sequenceR2).find(spacers[sp])
        if (spacer_loc != -1 or spacer_loc_rc != -1):
            spacer = sp
            P5_sample_BC = read_sequenceR2[P5_sample_BC_start:
                                           P5_sample_BC_start + P5_sample_BC_len]
            P7_sample_BC = read_sequenceR1[P7_sample_BC_start:
                                           P7_sample_BC_start + P7_sample_BC_len]
            if spacer_loc == -1:
                spacer_loc = spacer_loc_rc
                top_read = reverse_complement(read_sequenceR2)
                bot_read = read_sequenceR1
                break
            elif spacer_loc_rc == -1:
                top_read = read_sequenceR1
                bot_read = reverse_complement(read_sequenceR2)
                break
            elif spacer_loc > spacer_loc_rc:
                spacer_loc = spacer_loc_rc
                top_read = reverse_complement(read_sequenceR2)
                bot_read = read_sequenceR1
                break
            else:
                top_read = read_sequenceR1
                bot_read = reverse_complement(read_sequenceR2)
                break
        else:
            spacer_loc = read_sequenceR2.find(spacers[sp])
            spacer_loc_rc = reverse_complement(read_sequenceR1).find(spacers[sp])
            if (spacer_loc != -1 or spacer_loc_rc != -1):
                spacer = sp
                P5_sample_BC = read_sequenceR1[P5_sample_BC_start:
                                               P5_sample_BC_start + P5_sample_BC_len]
                P7_sample_BC = read_sequenceR2[P7_sample_BC_start:
                                               P7_sample_BC_start + P7_sample_BC_len]
                if spacer_loc == -1:
                    spacer_loc = spacer_loc_rc
                    top_read = reverse_complement(read_sequenceR1)
                    bot_read = read_sequenceR2
                    break
                elif spacer_loc_rc == -1:
                    top_read = read_sequenceR2
                    bot_read = reverse_complement(read_sequenceR1)
                    break
                elif spacer_loc > spacer_loc_rc:
                    spacer_loc = spacer_loc_rc
                    top_read = reverse_complement(read_sequenceR1)
                    bot_read = read_sequenceR2
                    break
                else:
                    top_read = read_sequenceR2
                    bot_read = reverse_complement(read_sequenceR1)
                    break

    return top_read, bot_read, spacer, spacer_loc, P5_sample_BC, P7_sample_BC

def raw_count_summary(run_name, input_csv=None):
    """
    just sum across all PAMs for each sample and timepoint to summarize raw counts
    """
    if input_csv is None:
        df_input = pd.read_csv('output/%s/PAMDA_1_raw_counts.csv.gz' % run_name)
    else:
        df_input = pd.read_csv(input_csv)

    df_output = df_input.groupby(['Sample', 'Spacer']).sum()

    if not os.path.exists('output/%s' % run_name):
        os.makedirs('output/%s/' % run_name)

    df_output.to_csv('output/%s/PAMDA_1_raw_counts_summary.csv' %
                     run_name)


def group_pams(input_df, pam_start, pam_length, max_pam_length, pam_orientation, control_spacers):
    """generate counts for the PAMs defined by indicated start position and PAM length"""
    print('grouping counts by indicated PAM bases')
    # make a copy to avoid mutating the original
    output_df = input_df.copy(deep=True)
    # make a new column for the PAM length we are interested in
    if control_spacers is None:
        control_spacers = []
    if pam_orientation == 'three_prime':
        selected_pams = output_df["PAM"].where(output_df.PAM.isin(control_spacers),
                                                           other=output_df['PAM'].str[pam_start:pam_start + pam_length])
    else:
        selected_pams = output_df["PAM"].where(output_df.PAM.isin(control_spacers),
                                                           other=output_df['PAM'].str[max_pam_length - pam_length - pam_start:max_pam_length - pam_start])
    output_df.loc[:, "PAM"] = selected_pams
    # group by sample and spacer and add up all PAMs that are equivalent
    output_df = output_df.groupby(by=["Sample", "Spacer", "PAM"], as_index=False).sum()
    output_df.reset_index(drop=True)
    return output_df

def add_t0_raw_counts(input_df, control_sample, control_sample_timepoint_fastq, timepoints):
    """Convert the control sample raw counts to t0 raw counts for each sample"""
    # make a copy
    output_df = input_df.copy(deep=True)

    # make sure the order of the spacers and PAMs is the same between the control samples and the other samples
    sample_names = pd.unique(output_df["Sample"])
    control_spacers = output_df.loc[output_df["Sample"] == control_sample, "Spacer"]
    control_pams = output_df.loc[output_df["Sample"] == control_sample, "PAM"]
    for sample_name in sample_names:
        sample_spacers = output_df.loc[output_df["Sample"] == sample_name, "Spacer"]
        sample_pams = output_df.loc[output_df["Sample"] == sample_name, "PAM"]
        assert(len(sample_spacers) == len(control_spacers)), "The following sample contains a different number of spacers than the control sample: " + sample_name
        assert all(np.equal(sample_pams.values, control_pams.values)), "Sample spacers do not match control spacers for this sample: " + sample_name
        assert all(np.equal(sample_spacers.values, control_spacers.values)), "Sample PAMs do not match control PAMS for this sample: " + sample_name

    # separate the control from the rest of the dataframe
    control_df = output_df.loc[output_df["Sample"] == control_sample, :]
    output_df = output_df.loc[output_df["Sample"] != control_sample, :]

    # make the control counts the t0 counts for each sample
    control_counts = control_df["Raw_Counts_" + str(control_sample_timepoint_fastq)].values
    control_counts = [control_counts for _ in range(len(sample_names)-1)]
    control_counts = np.concatenate(control_counts)
    output_df["Raw_Counts_0"] = control_counts
    raw_count_colnames = ["Raw_Counts_" + str(i) for i in range(len(timepoints))]
    colnames = ["Sample", "Spacer", "PAM"].extend(raw_count_colnames)
    output_df = output_df.reindex(columns=colnames)
    return output_df


def add_fractional_counts(input_df, timepoints):
    """For each sample spacer TP combo in the dataframe,
    express raw read counts as a fraction of the total reads for that sample-spacer-tp"""
    # make a copy
    output_df = input_df.copy(deep=True)

    # make a df of the total counts for each sample-spacer-tp combo
    total_counts_df = output_df.drop(["PAM"], axis=1)
    total_counts_df = total_counts_df.groupby(by=["Sample", "Spacer"]).sum()  # this has a multi index of sample and spacer

    # initialize the new columns
    fract_column_names = []
    for i in range(len(timepoints)):
        col = "Fractional_Counts_" + str(i)
        fract_column_names.append(col)
        output_df[col] = np.zeros(len(output_df))

    # get the column names of the raw counts
    raw_column_names = []
    for i in range(len(timepoints)):
        col = "Raw_Counts_" + str(i)
        raw_column_names.append(col)

    # divide each row of our data by the corresponding total counts to make fractional counts
    for index, row in output_df.iterrows():
        sample = row["Sample"]
        spacer = row["Spacer"]
        divisors = total_counts_df.loc[(sample, spacer)]
        new = list(np.divide(row.loc[raw_column_names], divisors))
        output_df.loc[index, fract_column_names] = new

    return output_df


def make_control_df(input_df, timepoints, top_n):
    """Choose the top n rows of the dataframe with largest uptrends to use as controls"""

    input_df = input_df.copy(deep=True)
    print("choosing control samples")

    # initialize new columns
    input_df["slopes"] = np.zeros(len(input_df))

    # populate new columns
    x = range(len(timepoints))
    for index, row in input_df.iterrows():
        y = [row['Fractional_Counts_' + str(i)] for i in range(len(timepoints))]
        slope = linregress(x, y)[0]
        input_df.loc[index, 'slopes'] = slope

    # take only the rows with the top n slopes for each sample and each spacer
    sample_spacer_dataframes = input_df.groupby(by=["Sample", "Spacer"], as_index=False)

    control_dfs_list = []
    for index, df in sample_spacer_dataframes:
        df = df.sort_values(by="slopes", ascending=False)
        df.reset_index(inplace=True, drop=True)
        df = df.loc[:top_n, :]
        control_dfs_list.append(df)

    control_df = pd.concat(control_dfs_list)
    control_df.reset_index(inplace=True, drop=True)

    return control_df


def correct_uptrends(input_df, control_spacers, timepoints, run_name, pam_start, pam_length, top_n):
    """Correct later timepoints when control samples trend upwards"""
    # make a copy
    output_df = input_df.copy(deep=True)

    # df of just the control spacers so that we can calculate the increases for these
    if control_spacers is not None:
        sel = np.isin(output_df["PAM"], list(control_spacers.keys()))
        control_df = output_df.loc[sel, :]
    else:
        control_df = make_control_df(input_df, timepoints, top_n)

    if not os.path.exists('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length)):
        os.makedirs('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length))
    control_df.to_csv('output/%s/PAM_start_%s_length_%s/control_sample_slopes.csv' %
                               (run_name, pam_start, pam_length), index=True)

    print("normalizing read counts")

    # initialize new columns
    correction_colnames = []
    for i in range(len(timepoints)):
        new = "Correction_Factor_" + str(i)
        control_df[new] = np.zeros(len(control_df))
        correction_colnames.append(new)

    fract_column_names = []
    for i in range(len(timepoints)):
        col = "Fractional_Counts_" + str(i)
        fract_column_names.append(col)

    # calculate median increase factors for the controls since TP0
    for i in range(len(timepoints)):
        control_df["Correction_Factor_" + str(i)] = control_df["Fractional_Counts_" + str(i)] / control_df["Fractional_Counts_0"] # increase factor since t0
    control_df = control_df.groupby(by=["Sample", "Spacer"]).median(numeric_only=True)
    correction_factors_df = control_df[correction_colnames]

    if not os.path.exists('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length)):
        os.makedirs('output/%s/PAM_start_%s_length_%s' % (run_name, pam_start, pam_length))
    correction_factors_df.to_csv('output/%s/PAM_start_%s_length_%s/count_correction_factors.csv' %
              (run_name, pam_start, pam_length), index=True)

    # divide by the correction factors
    for index, row in output_df.iterrows():
        sample = row.loc["Sample"]
        spacer = row.loc["Spacer"]
        correction_factors = list(correction_factors_df.loc[(sample, spacer), :])
        output_df.loc[index, fract_column_names] = np.divide(output_df.loc[index, fract_column_names], correction_factors)
    return output_df




def norm_to_t0_abundance(input_df, timepoints):
    """Express the fractional counts as a fraction on abundance of each PAM at t0"""
    # make a copy
    output_df = input_df.copy(deep=True)

    fract_column_names = []
    for i in range(len(timepoints)):
        col = "Fractional_Counts_" + str(i)
        fract_column_names.append(col)

    norm_column_names = []
    for i in range(len(timepoints)):
        col = "Norm_Counts_" + str(i)
        norm_column_names.append(col)
        output_df[col] = np.zeros(len(output_df))

    for i in range(len(timepoints)):
        output_df["Norm_Counts_" + str(i)] = np.divide(output_df["Fractional_Counts_" + str(i)],
                                                       output_df["Fractional_Counts_0"])

    return output_df




def check_inputs(RUN_NAME,
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
                 LOG_SCALE_HEATMAP):
    """
    perform some checks for input parameters
    """
    if not os.path.exists(BARCODE_CSV):
        raise Exception('BARCODE_CSV "%s" not found' % BARCODE_CSV)
    if not os.path.exists(FASTQ_DIR):
        raise Exception('fastq directory "%s" not found' % FASTQ_DIR)
    if CONTROL_RAW_COUNT_CSV != None:
        if not os.path.exists(CONTROL_RAW_COUNT_CSV):
            raise Exception('CONTROL_RAW_COUNT_CSV "%s" not found' % CONTROL_RAW_COUNT_CSV)

    fastqs = glob.glob(FASTQ_DIR + '/**/*R1*.fastq.gz', recursive=True)
    if len(fastqs) == 0:
        raise Exception('no fastq files found')
    for fastq in fastqs:
        fastqR1 = fastq
        fastqR2 = fastq.replace('R1', 'R2')
        fastq_name = fastqR1.split('/')[-1]
        fastq_name = fastq_name.split('_L00')[0]
        try:
            TIMEPOINT_FASTQ[fastq_name]
        except:
            warnings.warn('%s not found in TIMEPOINT_FASTQ. This fastq will be ignored.' % fastq_name)
    if not isinstance(MAX_PAM_LENGTH, int):
        raise Exception('MAX_PAM_LENGTH should be an integer value, you entered: %s' % MAX_PAM_LENGTH)
    if not isinstance(P5_SAMPLE_BARCODE_START, int):
        raise Exception('P5_SAMPLE_BARCODE_START should be an integer value, you entered: %s' % P5_SAMPLE_BARCODE_START)
    if not isinstance(P7_SAMPLE_BARCODE_START, int):
        raise Exception('P7_SAMPLE_BARCODE_START should be an integer value, you entered: %s' % P7_SAMPLE_BARCODE_START)
    if not isinstance(PAM_LENGTH, int):
        raise Exception('PAM_LENGTH should be an integer value, you entered: %s' % PAM_LENGTH)
    if not isinstance(PAM_START, int):
        raise Exception('PAM_START should be an integer value, you entered: %s' % PAM_START)
    if (CONTROL_RAW_COUNT_CSV == None) and (CONTROL_SAMPLE_TIMEPOINT_FASTQ == None):
        raise Exception('Either CONTROL_RAW_COUNT_CSV or CONTROL_SAMPLE_TIMEPOINT_FASTQ must be specified')
    if not (isinstance(CONTROL_SAMPLE_TIMEPOINT_FASTQ, int) or (CONTROL_SAMPLE_TIMEPOINT_FASTQ == None)):
        raise Exception(
            'CONTROL_SAMPLE_TIMEPOINT_FASTQ should be "None" or an integer value, you entered: %s' % CONTROL_SAMPLE_TIMEPOINT_FASTQ)
    if not isinstance(TOP_N_NORMALIZE, int):
        raise Exception('TOP_N_NORMALIZE should be an integer value, you entered: %s' % TOP_N_NORMALIZE)
    if not isinstance(READ_SUM_MIN, int):
        raise Exception('READ_SUM_MIN should be an integer value, you entered: %s' % READ_SUM_MIN)
    if not isinstance(TPS_SUM_MIN, int):
        raise Exception('TPS_SUM_MIN should be an integer value, you entered: %s' % TPS_SUM_MIN)
    if not isinstance(TIMEPOINTS, list):
        raise Exception('TIMEPOINTS should be a list, you entered: %s' % TIMEPOINTS)
    if not isinstance(INIT_RATE_EST, list):
        raise Exception('INIT_RATE_EST should be a list, you entered: %s' % INIT_RATE_EST)
    if not isinstance(AVERAGE_SPACER, bool):
        raise Exception("AVERAGE_SPACER should be 'True' or 'False', you entered: %s" % AVERAGE_SPACER)
    if not isinstance(LOG_SCALE_HEATMAP, bool):
        raise Exception("LOG_SCALE_HEATMAP should be 'True' or 'False', you entered: %s" % LOG_SCALE_HEATMAP)
    if (HEATMAP_FIXED_MIN != False and not (
            isinstance(HEATMAP_FIXED_MIN, int) or isinstance(HEATMAP_FIXED_MIN, float))):
        raise Exception(
            "HEATMAP_FIXED_MIN should be 'False' or a float or an integer, you entered: %s" % HEATMAP_FIXED_MIN)
    if (HEATMAP_FIXED_MAX != False and not (
            isinstance(HEATMAP_FIXED_MAX, int) or isinstance(HEATMAP_FIXED_MAX, float))):
        raise Exception(
            "HEATMAP_FIXED_MAX should be 'False' or a float or an integer, you entered: %s" % HEATMAP_FIXED_MAX)
    if PAM_ORIENTATION not in ['three_prime', 'five_prime']:
        raise Exception("please enter 'three_prime' or 'five_prime' for PAM_ORIENTATION")
    if PAM_LENGTH > 8:
        warnings.warn('PAM lengths longer than 8 are not recommended')
    if PAM_LENGTH < 1:
        raise Exception('please choose a PAM length >0')
    if PAM_START + PAM_LENGTH > MAX_PAM_LENGTH:
        raise Exception('PAM_START (%s) + PAM_LENGTH (%s) is greater than MAX_PAM_LENGTH (%s)'
                        % (PAM_START, PAM_LENGTH, MAX_PAM_LENGTH))
    if (PAM1_INDEX_RANK != None and PAM2_INDEX_RANK != None):
        if not isinstance(PAM1_INDEX_RANK, list):
            raise Exception('PAM1_INDEX_RANK should be a list, you entered: %s' % PAM1_INDEX_RANK)
        if not isinstance(PAM2_INDEX_RANK, list):
            raise Exception('PAM2_INDEX_RANK should be a list, you entered: %s' % PAM2_INDEX_RANK)
        if len(PAM1_INDEX_RANK) + len(PAM2_INDEX_RANK) != PAM_LENGTH:
            raise Exception(
                'The number of ranked PAM positions in PAM1_INDEX_RANK and PAM2_INDEX_RANK is not equal to PAM_LENGTH')
    elif ((PAM1_INDEX_RANK == None and PAM2_INDEX_RANK != None) or
          (PAM1_INDEX_RANK != None and PAM2_INDEX_RANK == None)):
        raise Exception('Please specify both PAM1_INDEX_RANK and PAM2_INDEX_RANK or leave both as "None".')


def translate_pam(pam_as_numbers, tranlsate_dict):
    pam_as_bases = ''
    for n in str(pam_as_numbers):
        pam_as_bases += tranlsate_dict[int(n)]
    return pam_as_bases

def plot_default_heatmap(df_output, avg_spacer, heatmap_fixed_min, heatmap_fixed_max, variant, variant_name_dict, run_name, pam_start, pam_length, pam2_nucleotide_rank, pam2_ids, spacer=None):
    """
    plot heatmap without fancy formatting
    """
    heatmap_min = None
    heatmap_max = None
    if heatmap_fixed_min:
        heatmap_min = heatmap_fixed_min
    if heatmap_fixed_max:
        heatmap_max = heatmap_fixed_max

    # generate heatmaps and save
    sns.set(font_scale=1)
    fig, ax = plt.subplots()
    plt.title(variant + '(' + variant_name_dict[variant] + ')', y=1)
    ax = sns.heatmap(df_output,
                     vmin=heatmap_min,
                     vmax=heatmap_max,
                     square=True,
                     cmap='Blues',
                     cbar=True,
                     cbar_kws={"shrink": 0.5},
                     linewidth=0.2,
                     linecolor="White",
                     xticklabels=[translate_pam(pam, pam2_nucleotide_rank) for pam in pam2_ids])
    fig.tight_layout(pad=4)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    save_heatmap(avg_spacer, run_name, pam_start, pam_length, variant, variant_name_dict, df_output, spacer)
    plt.close()


def plot_nice_heatmap(df_output, heatmap_fixed_min, heatmap_fixed_max, log_scale_heatmap, pam_length, split_pam_index,
                      columns, indices, avg_spacer, run_name, pam_start, variant, variant_name_dict, spacer=None):
    """
    plot pretty heatmaps with color-formatted axes labels (only for 2-5 nt PAMs)
    """
    heatmap_min = None
    heatmap_max = None
    if heatmap_fixed_min:
        heatmap_min = heatmap_fixed_min
    if heatmap_fixed_max:
        heatmap_max = heatmap_fixed_max

    if log_scale_heatmap:
        cbar_label = 'Log10(rate)'
    else:
        cbar_label = 'rate'

    # heatmap plotting
    axes = [4 * (pam_length - split_pam_index), 4 * (split_pam_index)]
    fig, ax1 = plt.subplots(1, figsize=(axes[0], axes[1]))
    sns.heatmap(df_output,
                ax=ax1,
                vmin=heatmap_min,
                vmax=heatmap_max,
                square=True,
                cmap='Blues',
                cbar=True,
                cbar_kws={'shrink': axes[1] / axes[0] / 2,
                          'label': cbar_label,
                          'aspect': 8},
                linewidth=0.2,
                linecolor="White",
                xticklabels=False,
                yticklabels=False)
    heatmap_size = fig.get_size_inches()

    # format the axis labels

    # colors of the bases of the axis labels
    colors = {'A': '#feb2b1', 'C': '#14c7fe', 'T': '#14f485', 'G': '#f8ffa3'}
    # scaling dict manually scales the colored axis labels
    # dict structure: {pam_len:{split_index:[x_width,x_height,y_width,y_height]}}
    scaling_dict = {2: {1: [1, 3.5, 4, 1.07]},
                    3: {1: [1, 1.7, 1, 2.15], 2: [1, 2, 4, 3.53]},
                    4: {1: [1, 0.7, 0.3, 5.5], 2: [1, 1.7, 1, 4.3], 3: [1, 1, 4, 14.1]},
                    5: {1: [1, 0.25, 0.08, 17], 2: [1, 0.6, 0.3, 11.5],
                        3: [1, 1, 1, 14.15], 4: [1, 0.3, 3, 56.5]}}
    x_text = [[columns[n][m] for n in range(len(columns))] for m in range(len(columns[0]))]
    x_text = x_text[::-1]
    x_color = [[colors[x_text[n][m]] for m in range(len(x_text[0]))] for n in range(len(x_text))]
    xtable = ax1.table(cellText=x_text, cellColours=x_color,
                       cellLoc='center', loc='top')
    for key, cell in xtable.get_celld().items():
        cell.set_linewidth(0)
    y_text = [[indices[n][m] for m in range(len(indices[0]))] for n in range(len(indices))]
    y_color = [[colors[y_text[n][m]] for m in range(len(y_text[0]))] for n in range(len(y_text))]
    y_widths = [0.06 for n in enumerate(y_text)]
    ytable = ax1.table(cellText=y_text, cellColours=y_color, colWidths=y_widths,
                       cellLoc='center', loc='left')
    for key, cell in ytable.get_celld().items():
        cell.set_linewidth(0)
    xtable.set_fontsize(8)
    ytable.set_fontsize(8)
    xtable.scale(scaling_dict[pam_length][split_pam_index][0],
                 scaling_dict[pam_length][split_pam_index][1])
    ytable.scale(scaling_dict[pam_length][split_pam_index][2],
                 heatmap_size[1] / scaling_dict[pam_length][split_pam_index][3])
    plt.tight_layout(pad=6)

    # save the heatmap
    save_heatmap(avg_spacer, run_name, pam_start, pam_length, variant, variant_name_dict, df_output, spacer)
    plt.close()


def spacer_correlation(df_input, sample, spacer1, spacer2, run_name, pam_start, pam_length, variant, variant_name_dict):
    """
    if len(spacers)==2, the correlation between the two spacers may be plotted
    """
    x = df_input['Rate_Constant_k'][(df_input['Sample'] == sample) & (df_input['Spacer'] == spacer1)].fillna(0)
    y = df_input['Rate_Constant_k'][(df_input['Sample'] == sample) & (df_input['Spacer'] == spacer2)].fillna(0)
    corr, pval = pearsonr(x, y)
    plt.figure(figsize=(6, 6))
    plt.scatter(x, y, color='Black')
    plt.title('Correlation %s versus %s \n Pearson r: %s' % (spacer1, spacer2, round(corr, 4)))
    plt.xlabel('Rate constant %s' % spacer1)
    plt.ylabel('Rate constant %s' % spacer2)
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.tight_layout()
    plt.savefig('figures/%s/PAM_start_%s_length_%s/PAMDA_spacer_correlation_%s_%s.pdf' %
                (run_name, pam_start, pam_length, variant, variant_name_dict[variant]))
    plt.close()

def save_heatmap(avg_spacer, run_name, pam_start, pam_length, variant, variant_name_dict, df_output, spacer):
    """
    save heatmaps
    """
    if avg_spacer:
        plt.savefig('figures/%s/PAM_start_%s_length_%s/PAMDA_HEATMAP_%s_%s.pdf' %
                    (run_name, pam_start, pam_length, variant, variant_name_dict[variant]))
        df_output.to_csv('figures/%s/PAM_start_%s_length_%s/PAMDA_HEATMAP_%s_%s.csv' %
                    (run_name, pam_start, pam_length, variant, variant_name_dict[variant]))
    else:
        plt.savefig('figures/%s/PAM_start_%s_length_%s/PAMDA_HEATMAP_%s_%s_%s.pdf' %
                    (run_name, pam_start, pam_length, variant, variant_name_dict[variant], spacer))
        df_output.to_csv('figures/%s/PAM_start_%s_length_%s/PAMDA_HEATMAP_%s_%s_%s.csv' %
                    (run_name, pam_start, pam_length, variant, variant_name_dict[variant], spacer))
