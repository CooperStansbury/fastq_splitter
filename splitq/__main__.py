import sys
import os
from datetime import datetime
import argparse
import splitq

if __name__ == '__main__':
    
    ############################################################################################
    # INPUT ARGUMENT DEFINITIONS
    ############################################################################################
    
    desc = """A Python3 commandline tool process Pore-C-SnakeMake outputs for PaohVis"""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-d", 
                        help="The path to a directory of .fastq files.")
    parser.add_argument("-od", 
                        nargs='?',
                        type=str, 
                        default="results/",
                        help="The path to a directory for all output files.")
    parser.add_argument("-overwrite", 
                    nargs='?',
                    type=bool, 
                    default=True,
                    help="If True, overwrite the directory at `-od`")
  
    ############################################################################################
    # INPUT ARGUMENT PARSING
    ############################################################################################
    args = parser.parse_args()
    
    # argument parsing - may need to handle more robustly
    FASTQ_DIR = args.d
    OUTPUT_DIR = args.od
    OVERWRITE = args.overwrite
    TODAY = datetime.today().strftime('%Y-%m-%d')
    
    # make output dir
    os.makedirs(OUTPUT_DIR, mode=0o777, exist_ok=OVERWRITE)
    
    splitter = splitq.SplitQ(FASTQ_DIR, OUTPUT_DIR)
    
#     print(TODAY)

    ############################################################################################
    # DATA PROCESSING
    ############################################################################################
    
    
    ############################################################################################
    # OUTPUTS
    ############################################################################################
  