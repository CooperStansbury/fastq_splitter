import sys
import os
from datetime import datetime
import argparse


def split_fastq(infile):
    """A function to digest the reads by a cut site prior to mapping 
    against the reference.

    This is adapted from:

    https://github.com/dekkerlab/MC-3C_scripts

    Parameters:
    -----------------------------
        : infile (str): path to the input fastq file
        : outfile (str): path to the output fastq file
    """
    output_handle = open(outfile, "w")
    handle = open(infile, "r")
    records = list(SeqIO.parse(handle, "fastq"))

    handle.close()

    for record in records:
        restriction_site_len = len(Restriction.NlaIII.site)
        cut_sites = Restriction.NlaIII.search(record.seq)

        print(cut_sites)
        break

        if len(cut_sites) == 0:
            output_handle.write(record.format("fastq"))
        else :
            subread_counter = 0

            for i in range(0, len(cut_sites)):
                if i == 0:
                    subread = record[0 : cut_sites[i]-1]
                else:
                    subread = record[cut_sites[i-1]-1 : cut_sites[i]-1]

                if len(subread.seq) <= 4:
                    continue

                subread.id += "/0" + str(subread_counter)
                subread.name += "/0" + str(subread_counter)
                subread.description += "/0" + str(subread_counter)
                subread_counter += 1
                output_handle.write(subread.format("fastq"))

            subread = record[cut_sites[i]-1 :]
            subread.id += "/0" + str(subread_counter)
            subread.name += "/0" + str(subread_counter)
            subread.description += "/0" + str(subread_counter)

            if len(subread.seq) <= 4:
                continue

            output_handle.write(subread.format("fastq"))
    output_handle.close()


if __name__ == '__main__':
    
    ############################################################################################
    # INPUT ARGUMENT DEFINITIONS
    ############################################################################################
    
    desc = """A Python3 commandline tool process Pore-C-SnakeMake outputs for PaohVis"""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-i", 
                        help="Input merged fastq file.")
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
    FASTQ = args.i
    OUTPUT_DIR = args.od
    OVERWRITE = args.overwrite
    
    # make output dir, handle overwriting flag
    os.makedirs(OUTPUT_DIR, mode=0o777, exist_ok=OVERWRITE)
    
    ############################################################################################
    # DATA PROCESSING
    ############################################################################################
    
    """TODO: finish this function, perform QA, and write out new fastq file to the newly created
    directory"""
    
    split_fastq(FASTQ)
    
    ############################################################################################
    # OUTPUTS
    ############################################################################################
    
    """TODO: store outputs"""