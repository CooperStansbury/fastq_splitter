import sys
import os
from datetime import datetime
import argparse
import pandas as pd

from Bio import Restriction
from Bio.Seq import Seq
from Bio import SeqIO

def partition(alist, indices):
    """A function to split a list based on item indices 
    
    Parameters:
    -----------------------------
        : alist (list): a list to be split
        : indices (list): list of indices on which to divide the input list
        
    Returns:
    -----------------------------
        : splits (list): a list of subreads based on cut sites
    """
    return [alist[i:j] for i, j in zip([0]+indices, indices+[None])]


def split_fastq(infile, N=300):
    """A function to digest the reads by a cut site prior to mapping 
    against the reference.

    This is adapted from:

    https://github.com/dekkerlab/MC-3C_scripts

    Parameters:
    -----------------------------
        : infile (str): path to the input fastq file
        : N (int): number of base pairs surrounding the cut-site to define as pairs
        
        
    Returns:
    -----------------------------
        : pairs_table (pd.DataFrame): a dataframe with subread pairs and metadata
    """
    handle = open(infile, "r")
    records = list(SeqIO.parse(handle, "fastq"))
    handle.close()
    
    new_rows = []

    # iterate through each reaa of the fastq file
    for record in records:
        
        # digest read on NlaIII sites
        restriction_site_len = len(Restriction.NlaIII.site)
        cut_sites = Restriction.NlaIII.search(record.seq)
        cut_sites = [x-1 for x in cut_sites] # cut site indices are the bp AFTER the CATG enzyme
       
        if len(cut_sites) > 0:
            
            # split the read on NlaIII cites
            splits = partition(record, cut_sites)
        
            # loop through pairs of subreads (i, i+1)
            subread_count = 0
            for idx in range(len(splits) - 1):
                subread_count += 1
                
                left_read = splits[idx]
                right_read = splits[idx + 1]
                
                # trim thye CATG sequence from the reads
                if left_read.seq.endswith("CATG"):
                    left_read = left_read[:-4]
                    
                if right_read.seq.endswith("CATG"):
                    right_read = right_read[:-4]
                
                # trim subreads based on parameter N
                if len(left_read) > N:
                    left_read = left_read[-N:]
                    
                if len(right_read) > N:
                    right_read = right_read[:N]
                
                # get read metadata from left read - assume consistent between subreads
                metadata = dict(item.split("=") for item in left_read.description.split(" ") if "=" in item)
                
                # construct new data record
                new_record = {
                    'read_id' : left_read.id,
                    'subread_id' : subread_count,
                    'left_read' : "".join(left_read.seq),
                    'right_read' : "".join(right_read.seq),
                    'subread_max_length' : N,
                }
                
                # add metadata
                for k, v in metadata.items():
                    new_record[k] = v
                    
                new_rows.append(new_record)
        
    # build pairs table 
    pairs_table = pd.DataFrame(new_rows)
    return pairs_table
        


if __name__ == '__main__':
    
    ############################################################################################
    # INPUT ARGUMENT DEFINITIONS
    ############################################################################################
    
    desc = """A Python3 commandline tool process Pore-C-SnakeMake outputs for PaohVis"""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-i", 
                        help="Input merged fastq file.")
    parser.add_argument("-n", 
                        nargs='?',
                        default=300,
                        help="Number of base pairs surrounding a cut site.")
    parser.add_argument("-o", 
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
    N = int(args.n)
    OUTPUT_DIR = args.o
    OVERWRITE = args.overwrite
    
    # make output dir, handle overwriting flag
    os.makedirs(OUTPUT_DIR, mode=0o777, exist_ok=OVERWRITE)
    
    ############################################################################################
    # DATA PROCESSING
    ############################################################################################
    
    # generate pairs
    pairs_table = split_fastq(FASTQ, N)
        
    # get input file basename and generate output filename 
    base = os.path.basename(FASTQ)
    basename = os.path.splitext(base)[0]
    new_filename = f"{OUTPUT_DIR}{basename}_PAIRS.csv"
    
    ############################################################################################
    # OUTPUTS
    ############################################################################################
    
    pairs_table.to_csv(new_filename, index=False)
    print(f"Done saving: `{new_filename}`")