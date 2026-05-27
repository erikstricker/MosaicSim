#!/usr/bin/env python3

import sys
import os
import time
import pysam
import subprocess
from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def stream_merge_replace(input_bam, modified_bam, output_bam):
    """
    Synchronized streaming replacement using coordinate-sorted pointers.
    Uses native pysam output for maximum stability and speed.
    """
    print("Running synchronized merge-replace...")
    
    # Open input and modified BAMs
    orig_bam = pysam.AlignmentFile(input_bam, 'rb')
    mod_bam = pysam.AlignmentFile(modified_bam, 'rb')
    
    # Open output BAM using the original as a template (preserves headers/index)
    out_bam = pysam.AlignmentFile(output_bam, 'wb', template=orig_bam)
    
    # Setup iterator for modified reads
    mod_iter = iter(mod_bam)
    curr_mod = next(mod_iter, None)
    
    replaced_count = 0
    
    # Stream through the original BAM
    for orig_aln in orig_bam:
        # Check if the current modified read matches the original read's physical properties
        if curr_mod and (curr_mod.reference_id == orig_aln.reference_id and 
                         curr_mod.pos == orig_aln.pos and 
                         curr_mod.query_name == orig_aln.query_name):
            
            # Inject modified sequences into the original read structure
            orig_aln.cigarstring = curr_mod.cigarstring
            orig_aln.query_sequence = curr_mod.query_sequence
            orig_aln.query_qualities = curr_mod.query_qualities
            orig_aln.set_tags(curr_mod.get_tags())
            
            out_bam.write(orig_aln)
            replaced_count += 1
            curr_mod = next(mod_iter, None)
        else:
            out_bam.write(orig_aln)
            
    out_bam.close()
    orig_bam.close()
    mod_bam.close()
    
    print(f"\n--- SUCCESS: Replaced {replaced_count} reads using pointer-sync! ---\n")

def main(args):
    if args.time_run:
        start_time = time.time()
        print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Ensure output directory exists
    out_dir = os.path.dirname(args.out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        
    stream_merge_replace(args.input_bam, args.modified_bam, args.out_file)
    
    print("Indexing final BAM...")
    # The output is sorted because we kept the original read's position/reference_id
    subprocess.run(["samtools", "index", args.out_file], check=True)

    if args.time_run:
        end_time = time.time()
        print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Total run time: {end_time - start_time:.2f} seconds")
    
    print("Done.")

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-b", "--input_bam", required=True, help="Unmodified input bam file")
    parser.add_argument("-m", "--modified_bam", required=True, help="BAM file with modified reads only")
    parser.add_argument("-o", "--out_file", required=True, help="Output file name")
    parser.add_argument("--time_run", action="store_true", default=False, help="Export start time, end time, and total run time")
    
    args = parser.parse_args()
    main(args)
