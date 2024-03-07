import os
import sys
from library_QC_inputs import *
from library_QC_functions import *



def library_QC(RUN_NAME,
               CONTROL_BARCODE_CSV,
               CONTROL_FASTQ_DIR,
               CONTROL_FASTQ,
               PAM_ORIENTATION,
               PAM_LENGTH,
               PAM_START,
               MAX_PAM_LENGTH=8,
               SPACERS={'SPACER1': 'GGGCACGGGCAGCTTGCCGG', 'SPACER2': 'GTCGCCCTCGAACTTCACCT'},
               CONTROL_SPACERS=None,
               P5_SAMPLE_BARCODE_START=2,
               P7_SAMPLE_BARCODE_START=2):
    print('Begin library QC')

    library_QC_check_inputs(RUN_NAME,
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
                            P7_SAMPLE_BARCODE_START)

    control_fastq2count(RUN_NAME,
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
                        P7_SAMPLE_BARCODE_START)

    rawcount2PAMcount(RUN_NAME,
                      PAM_ORIENTATION,
                      PAM_LENGTH,
                      PAM_START,
                      MAX_PAM_LENGTH)

    print('Library QC complete')


if __name__ == "__main__":
    os.chdir('../')
    sys.path.append('code')

    library_QC(RUN_NAME,
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
               P7_SAMPLE_BARCODE_START)