from functions import *
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import skew
from fastq2count import fastq2count


def library_QC_check_inputs(RUN_NAME,
                            CONTROL_BARCODE_CSV,
                            CONTROL_FASTQ_DIR,
                            CONTROL_FASTQ,
                            PAM_ORIENTATION,
                            PAM_LENGTH,
                            PAM_START,
                            MAX_PAM_LENGTH,
                            SPACERS,
                            CONTROL_SPACERS,
                            P5_SAMPLE_BARCODE_START,
                            P7_SAMPLE_BARCODE_START):
    """
    perform some checks for input parameters
    """
    if not os.path.exists(CONTROL_BARCODE_CSV):
        raise Exception('CONTROL_BARCODE_CSV "%s" not found' % CONTROL_BARCODE_CSV)
    if not os.path.exists(CONTROL_FASTQ_DIR):
        raise Exception('CONTROL_FASTQ_DIR "%s" not found' % CONTROL_FASTQ_DIR)
    fastqs = glob.glob(CONTROL_FASTQ_DIR + '/**/*R1*.fastq.gz', recursive=True)
    if len(fastqs) == 0:
        raise Exception('no fastq files found')
    fastq_names = []
    for fastq in fastqs:
        fastqR1 = fastq
        fastqR2 = fastq.replace('R1', 'R2')
        fastq_name = fastqR1.split('/')[-1]
        fastq_name = fastq_name.split('_L00')[0]
        fastq_names.append(fastq_name)
        if fastq_name != CONTROL_FASTQ:
            warnings.warn('%s is not the CONTROL_FASTQ. This fastq will be ignored.' % fastq_name)
    if CONTROL_FASTQ not in fastq_names:
        raise Exception('CONTROL_FASTQ %s not found' % CONTROL_FASTQ)
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
    if PAM_ORIENTATION not in ['three_prime', 'five_prime']:
        raise Exception("please enter 'three_prime' or 'five_prime' for PAM_ORIENTATION")
    if PAM_LENGTH > 8:
        warnings.warn('PAM lengths longer than 8 are not recommended')
    if PAM_LENGTH < 1:
        raise Exception('please choose a PAM length >0')
    if PAM_START + PAM_LENGTH > MAX_PAM_LENGTH:
        raise Exception('PAM_START (%s) + PAM_LENGTH (%s) is greater than MAX_PAM_LENGTH (%s)'
                        % (PAM_START, PAM_LENGTH, MAX_PAM_LENGTH))


def control_fastq2count(run_name,
                        barcode_csv,
                        fastq_dir,
                        control_fastq,
                        pam_orientation,
                        pam_length,
                        pam_start,
                        max_pam_len=8,
                        spacers={'SPACER1': 'GGGCACGGGCAGCTTGCCGG', 'SPACER2': 'GTCGCCCTCGAACTTCACCT'},
                        control_spacers=None,
                        P5_sample_BC_start=2,
                        P7_sample_BC_start=2):
    """
    generate raw read counts from a fastq file containing control reads
    just runs fastq2count with a single timepoint_fastq and a single timepoint
    """

    timepoint_fastq = {control_fastq: 0}
    timepoints = [0, 1]

    fastq2count(run_name,
                barcode_csv,
                fastq_dir,
                timepoint_fastq,
                pam_orientation,
                timepoints=timepoints,
                max_pam_len=max_pam_len,
                spacers=spacers,
                control_spacers=control_spacers,
                P5_sample_BC_start=P5_sample_BC_start,
                P7_sample_BC_start=P7_sample_BC_start)


def rawcount2PAMcount(run_name,
                      pam_orientation,
                      pam_length,
                      pam_start,
                      max_pam_length=8):
    """
    generate counts for the PAMs defined by indicated start position and PAM length
    generate QC plots
    """

    df_input = pd.read_csv('output/%s/PAMDA_1_raw_counts.csv.gz' % (run_name))

    # generate counts for the PAMs defined by indicated start position and PAM length
    print('grouping counts by indicated PAM bases')

    column_sort = {'Sample': 'first', 'Spacer': 'first'}
    count_columns = df_input.columns.values[3:]
    new_columns = ['Sample', 'Spacer', 'PAM']
    for count_column in count_columns:
        column_sort[count_column] = 'sum'
        new_columns.append(count_column)
    if pam_orientation == 'three_prime':
        df_input['selected_PAM'] = df_input['PAM'].where(df_input['PAM'] == 'fixed',
                                                         other=df_input['PAM'].str[pam_start:pam_start + pam_length])
    else:
        df_input['selected_PAM'] = df_input['PAM'].where(df_input['PAM'] == 'fixed',
                                                         other=df_input['PAM'].str[max_pam_length - pam_length - pam_start:max_pam_length - pam_start])
    df_list = []
    for sample in df_input['Sample'].unique().tolist():
        for spacer in df_input['Spacer'].unique().tolist():
            temp_df = df_input[(df_input['Sample'] == sample) & (df_input['Spacer'] == spacer)] \
                .groupby(['selected_PAM'], as_index=False).agg(column_sort)
            df_list.append(temp_df)
    df = pd.concat(df_list)
    df = df.rename(columns={'selected_PAM': 'PAM'})
    df = df.loc[:, new_columns]
    df = df.reset_index(drop=True)

    for sample in df['Sample'].unique().tolist():
        QC_metrics(run_name, df[df['Sample'] == sample])

    if not os.path.exists('output/%s' % run_name):
        os.makedirs('output/%s' % run_name)

    df.to_csv('output/%s/PAMDA_1_raw_counts_PAM_start_%s_length_%s.csv' %
              (run_name, pam_start, pam_length), index=False)


def QC_metrics(run_name, df):
    """
    Simple QC plots for PAM library composition
    """

    plt.switch_backend('agg')

    if not os.path.exists('figures/%s' % run_name):
        os.makedirs('figures/%s' % run_name)

    for spacer in df['Spacer'].unique().tolist():
        df_spacer = df[df['Spacer'] == spacer]
        counts = df_spacer['Raw_Counts_1'].sort_values(ascending=False).tolist()
        PAM_count = len(df_spacer)
        top10 = int(0.1 * PAM_count)
        bot10 = int(0.9 * PAM_count)
        val90 = counts[top10]
        val10 = counts[bot10]
        if val10 == 0:
            ratio9010 = 'nan (10th percentile is zero)'
        else:
            ratio9010 = round(val90 / val10, 4)
        skewness = round(skew(counts), 4)

        if counts[-1] == 0:
            ratiomaxmin = 'nan (min is zero)'
        else:
            ratiomaxmin = round(counts[0] / counts[-1], 4)
        normalized_reads = df_spacer['Raw_Counts_1'].tolist() / df_spacer['Raw_Counts_1'].sum()

        plt.figure(figsize=(12, 4))
        plt.subplot(1, 2, 1)
        plt.hist(df_spacer['Raw_Counts_1'], color='Black')
        plt.title('PAM read count histogram for %s \n max:min ratio: %s \n 90:10 ratio: %s \n skewness: %s' %
                  (spacer, ratiomaxmin, ratio9010, skewness))
        plt.xlabel('Read counts')
        plt.ylabel('PAM count')
        plt.subplot(1, 2, 2)
        plt.plot([x + 1 for x in range(PAM_count)], sorted(normalized_reads), color='Black')
        plt.ylim(bottom=0)
        plt.title('PAM representation for %s \n max:min ratio: %s \n 90:10 ratio: %s \n skewness: %s' %
                  (spacer, ratiomaxmin, ratio9010, skewness))
        plt.xlabel('PAMs')
        plt.ylabel('Proportion of library')
        plt.tight_layout()
        plt.savefig('figures/%s/%s_PAMDA_library_QC_%s.pdf' % (run_name, run_name, spacer))
        plt.close()

        print('PAM library %s max:min ratio: %s' % (spacer, ratiomaxmin))
        print('PAM library %s 90:10 ratio: %s' % (spacer, ratio9010))
        print('PAM library %s skewness: %s' % (spacer, skewness))