import pandas as pd
from inputs import *
from functions import *

def fastq2count(run_name,
                barcode_csv,
                fastq_dir,
                timepoint_fastq,
                pam_orientation,
                timepoints,
                max_pam_len,
                spacers,
                control_spacers=None,
                control_spacer_mapping=None,
                P5_sample_BC_start=2,
                P7_sample_BC_start=2):
    """
    generate raw PAM read counts from fastq files
    """
    # check inputs
    try:
        variant_ids = pd.read_csv(barcode_csv)
    except:
        raise Exception('BARCODE_CSV "%s" not found' % barcode_csv)

    fastqs = glob.glob(fastq_dir + '/**/*R1*.fastq.gz', recursive=True)
    if len(fastqs) == 0:
        raise Exception('no fastq files found')

    if pam_orientation not in ['three_prime', 'five_prime']:
        raise Exception("please enter 'three_prime' or 'five_prime' for PAM_ORIENTATION")

    P5_sample_BCs = variant_ids['P5_sample_barcode'].tolist()
    P5_sample_BC_len = len(P5_sample_BCs[0])
    P7_sample_BCs = variant_ids['P7_sample_barcode'].tolist()
    P7_sample_BC_len = len(P5_sample_BCs[0])

    nt_complement = dict({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '_': '_', '-': '-'})

    nucleotides = ['A', 'T', 'C', 'G']
    total_pam_space = [''.join(p) for p in itertools.product(nucleotides, repeat=max_pam_len)]

    variant_dict = {}
    for index, row in variant_ids.iterrows():
        variant_dict[str(row['P5_sample_barcode']) + '_' + str(row['P7_sample_barcode'])] = row['sample']

    store_all_data = {}
    store_all_control_data = {}
    norm_counts_scale = {}


    for sample in variant_ids['sample']:
        store_all_data[sample] = {spacer: {x: [0] * (len(timepoints) - 1) for x in total_pam_space}
                                  for spacer in spacers}
        if control_spacers is not None:
            store_all_control_data[sample] = {control_spacer: {'fixed': [0] * (len(timepoints) - 1)}
                                              for control_spacer in control_spacers}

    pbar1 = tqdm(desc='fastq files: ', total=len(fastqs))
    pbar2 = tqdm(desc='reads: ')

    for fastq in fastqs:

        fastqR1 = fastq
        fastqR2 = fastq.replace('R1', 'R2')
        fastq_name = fastqR1.split('/')[-1]
        fastq_name = fastq_name.split('_L00')[0]

        try:
            timepoint = timepoint_fastq[fastq_name]
        except:
            continue

        if fastqR1.endswith('.gz'):
            infileR1 = gzip.open(fastqR1, 'rt')
            infileR2 = gzip.open(fastqR2, 'rt')
        else:
            infileR1 = open(fastqR1, 'r')
            infileR2 = open(fastqR2, 'r')

        wrong_barcode = 0
        wrong_spacer = 0
        total_reads = 0
        counted_reads = 0
        level1 = 0
        level2 = 0
        level1rc = 0
        level2rc = 0

        while infileR1.readline() and infileR2.readline():
            read_sequenceR1 = infileR1.readline().strip()
            infileR1.readline()
            infileR1.readline()
            read_sequenceR2 = infileR2.readline().strip()
            infileR2.readline()
            infileR2.readline()

            total_reads += 1

            top_read, bot_read, spacer, spacer_loc, P5_sample_BC, P7_sample_BC = \
                find_BCs_and_spacer(spacers, control_spacers, read_sequenceR1, read_sequenceR2,
                                    P5_sample_BC_start, P5_sample_BC_len,
                                    P7_sample_BC_start, P7_sample_BC_len)

            if spacer_loc == -1:
                wrong_spacer += 1
                continue

            if P5_sample_BC in P5_sample_BCs and P7_sample_BC in P7_sample_BCs:
                barcode_pair = P5_sample_BC + '_' + P7_sample_BC
                if barcode_pair in variant_dict.keys():
                    if pam_orientation == 'three_prime':
                        spacer3p = spacer_loc + len(spacers[spacer])
                        PAM = top_read[spacer3p: spacer3p + max_pam_len]
                        try:
                            store_all_data[variant_dict[barcode_pair]][spacer][PAM][timepoint] += 1
                            counted_reads += 1
                        except:
                            try: # try to see if it goes in the control dictionary if there is not a key for it in the all_data dictionary
                                store_all_control_data[variant_dict[barcode_pair]][spacer]["fixed"][timepoint] += 1
                                counted_reads += 1
                            except:
                                pass
                    elif pam_orientation == 'five_prime':
                        PAM = top_read[spacer_loc - max_pam_len: spacer_loc]
                        try:
                            store_all_data[variant_dict[barcode_pair]][spacer][PAM][timepoint] += 1
                            counted_reads += 1
                        except:
                            try:
                                store_all_control_data[variant_dict[barcode_pair]][spacer]["fixed"][timepoint] += 1
                                counted_reads += 1
                            except:
                                pass
            else:
                wrong_barcode += 1

            pbar2.update()

        pbar2.reset()
        pbar1.update()

        write_out = str(round(float(counted_reads) / float(total_reads) * 100, 2)) \
                    + '% of reads mapped from ' + str(fastq_name) + ' (' + str(counted_reads) + ' reads)'
        tqdm.write(write_out, file=sys.stdout)

    pbar1.close()
    pbar2.close()

    # output raw count results as a csv

    print('writing compressed CSV output')

    if not os.path.exists('output/%s' % run_name):
        os.makedirs('output/%s' % run_name)

    with gzip.open('output/%s/PAMDA_1_raw_counts.csv.gz' % (run_name), mode='wb') as f_out:
        f_out.write((','.join(map(str, ['Sample', 'Spacer', 'PAM'] +
                                  ['Raw_Counts_' + str(x)
                                   for x in range(1, len(timepoints))])) + '\n').encode('utf-8'))
        for fastq in store_all_data:
            for spacer in store_all_data[fastq]:
                for pam in store_all_data[fastq][spacer]:
                    total_info = [fastq, spacer, pam] + store_all_data[fastq][spacer][pam]
                    f_out.write((','.join(map(str, total_info)) + '\n').encode('utf-8'))
            if control_spacers is not None:
                for spacer in store_all_control_data[fastq]:
                    for pam in store_all_control_data[fastq][spacer]:
                        total_info = [fastq, spacer, pam] + store_all_control_data[fastq][spacer][pam]
                        f_out.write((','.join(map(str, total_info)) + '\n').encode('utf-8'))

    # also output summary csv file
    print('summarizing raw read counts')
    raw_count_summary(run_name)

#-------------------------------------

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

    # run the complete pamda pipeline:
    fastq2count(RUN_NAME,
                BARCODE_CSV,
                FASTQ_DIR,
                TIMEPOINT_FASTQ,
                PAM_ORIENTATION,
                TIMEPOINTS,
                MAX_PAM_LENGTH,
                SPACERS,
                CONTROL_SPACERS,
                CONTROL_SPACER_MAPPING,
                P5_SAMPLE_BARCODE_START,
                P7_SAMPLE_BARCODE_START)
