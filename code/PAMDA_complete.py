
from inputs import *
from functions import *
from fastq2count import fastq2count
from rawcount2normcount import rawcount2normcount
from normcount2rate import normcount2rate
from rate2heatmap import rate2heatmap



def PAMDA_complete(RUN_NAME,
                   BARCODE_CSV,
                   FASTQ_DIR,
                   TIMEPOINT_FASTQ,
                   PAM_ORIENTATION,
                   PAM_LENGTH,
                   PAM_START,
                   CONTROL_RAW_COUNT_CSV,
                   CONTROL_SAMPLE,
                   CONTROL_SAMPLE_TIMEPOINT_FASTQ=None,
                   TIMEPOINTS=[0, 60, 480, 1920],
                   MAX_PAM_LENGTH=8,
                   SPACERS={'SPACER1': 'GGGCACGGGCAGCTTGCCGG',
                            'SPACER2': 'GTCGCCCTCGAACTTCACCT'},
                   CONTROL_SPACERS=None,
                   P5_SAMPLE_BARCODE_START=2,
                   P7_SAMPLE_BARCODE_START=2,
                   USE_TIMEPOINTS=None,
                   TOP_N_NORMALIZE=5,
                   INIT_RATE_EST=[0.0001, 0.001, 0.01],
                   READ_SUM_MIN=4,
                   TPS_SUM_MIN=1,
                   PAM1_NT_RANK={1: 'A', 2: 'C', 3: 'G', 4: 'T'},
                   PAM2_NT_RANK={1: 'A', 2: 'C', 3: 'G', 4: 'T'},
                   PAM1_INDEX_RANK=None,
                   PAM2_INDEX_RANK=None,
                   AVERAGE_SPACER=True,
                   HEATMAP_FIXED_MIN=False,
                   HEATMAP_FIXED_MAX=False,
                   LOG_SCALE_HEATMAP=True):
    """
    Runs the complete PAMDA analysis from fastq files to heat map visualization.
    """

    # perform some checks
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

    # run the complete pipeline

    print('BEGIN: generate counts from fastqs')
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
    print('FINISHED: generate counts from fastqs \n')

    print('BEGIN: convert raw counts to normalized counts')
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
    print('FINISHED: convert raw counts to normalized counts \n')

    print('BEGIN: calculate rates from normalized counts')
    normcount2rate(RUN_NAME,
                   PAM_LENGTH,
                   PAM_START,
                   TIMEPOINTS,
                   INIT_RATE_EST,
                   READ_SUM_MIN,
                   TPS_SUM_MIN,
                   USE_TIMEPOINTS)
    print('FINISHED: calculate rates from normalized counts \n')

    print('BEGIN: plot heat maps')
    rate2heatmap(RUN_NAME,
                 BARCODE_CSV,
                 PAM_LENGTH,
                 PAM_START,
                 SPACERS,
                 PAM1_NT_RANK,
                 PAM2_NT_RANK,
                 PAM1_INDEX_RANK,
                 PAM2_INDEX_RANK,
                 AVERAGE_SPACER,
                 HEATMAP_FIXED_MIN,
                 HEATMAP_FIXED_MAX,
                 LOG_SCALE_HEATMAP,
                 control_spacers=CONTROL_SPACERS)
    print('FINISHED: plot heat maps')

    print(" _______  _______  __   __  ______   _______ ")
    print("|       ||   _   ||  |_|  ||      | |   _   |")
    print("|    _  ||  |_|  ||       ||  _    ||  |_|  |")
    print("|   |_| ||       ||       || | |   ||       |")
    print("|    ___||       ||       || |_|   ||       |")
    print("|   |    |   _   || ||_|| ||       ||   _   |")
    print("|___|    |__| |__||_|   |_||______| |__| |__|")
    print('complete')

#-----------------------------------------------------------------------------------------------------------------------------#

if __name__ == "__main__":
    print("HELLO from PAMDA complete")
    print(os.getcwd())
    os.chdir('../')
    sys.path.append('code')

    PAMDA_complete(RUN_NAME,
                   BARCODE_CSV,
                   FASTQ_DIR,
                   TIMEPOINT_FASTQ,
                   PAM_ORIENTATION,
                   PAM_LENGTH,
                   PAM_START,
                   CONTROL_RAW_COUNT_CSV,
                   CONTROL_SAMPLE,
                   CONTROL_SAMPLE_TIMEPOINT_FASTQ = CONTROL_SAMPLE_TIMEPOINT_FASTQ,
                   TIMEPOINTS = TIMEPOINTS,
                   MAX_PAM_LENGTH = MAX_PAM_LENGTH,
                   SPACERS = SPACERS,
                   CONTROL_SPACERS = CONTROL_SPACERS,
                   P5_SAMPLE_BARCODE_START = P5_SAMPLE_BARCODE_START,
                   P7_SAMPLE_BARCODE_START = P7_SAMPLE_BARCODE_START,
                   USE_TIMEPOINTS = USE_TIMEPOINTS,
                   TOP_N_NORMALIZE = TOP_N_NORMALIZE,
                   INIT_RATE_EST = INIT_RATE_EST,
                   READ_SUM_MIN = READ_SUM_MIN,
                   TPS_SUM_MIN = TPS_SUM_MIN,
                   PAM1_NT_RANK = PAM1_NT_RANK,
                   PAM2_NT_RANK = PAM2_NT_RANK,
                   PAM1_INDEX_RANK = PAM1_INDEX_RANK,
                   PAM2_INDEX_RANK = PAM2_INDEX_RANK,
                   AVERAGE_SPACER = AVERAGE_SPACER,
                   HEATMAP_FIXED_MIN = HEATMAP_FIXED_MIN,
                   HEATMAP_FIXED_MAX = HEATMAP_FIXED_MAX,
                   LOG_SCALE_HEATMAP = LOG_SCALE_HEATMAP)