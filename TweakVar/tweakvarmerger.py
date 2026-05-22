#!/usr/bin/env python3

import sys
import os
import time
import shutil
import pysam
import subprocess
from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def check_dependencies():
    """Ensure samtools is installed and available in the system PATH."""
    if shutil.which("samtools") is None:
        sys.exit("Error: 'samtools' is not installed or not in your PATH. It is required for fast streaming.")

def load_modified_reads(modified_bam_path):
    """Loads modified reads directly from the BAM into a dictionary as raw SAM strings."""
    print("Loading modified reads into memory...")
    modified_reads = {}
    
    with pysam.AlignmentFile(modified_bam_path, 'rb') as mbam:
        for aln in mbam:
            # We use (QNAME, string(FLAG)) as the unique key to prevent paired-end collisions
            key = (aln.query_name, str(aln.flag))
            # pysam's to_string() gives the SAM format, we just append the newline
            modified_reads[key] = aln.to_string() + "\n"
            
    return modified_reads

def stream_and_replace(input_bam, output_bam, modified_reads):
    """Coordinates the samtools background processes and streams the substitution."""
    print("Streaming and substituting reads via subprocess pipeline...")
    
    # 1. Spawn background samtools to decompress the input BAM to stdout
    p_in = subprocess.Popen(
        ["samtools", "view", "-h", input_bam],
        stdout=subprocess.PIPE,
        text=True # Ensures we read strings, not bytes
    )
    
    # 2. Spawn background samtools to compress stdin to the output BAM
    p_out = subprocess.Popen(
        ["samtools", "view", "-b", "-o", output_bam],
        stdin=subprocess.PIPE,
        text=True
    )
    
    # 3. Read from p_in, substitute if needed, and write to p_out
    for line in p_in.stdout:
        if line.startswith('@'):
            p_out.stdin.write(line)
            continue
            
        # Split only the first two columns (QNAME and FLAG)
        cols = line.split('\t', 2)
        key = (cols[0], cols[1])
        
        # Exact mate-to-mate substitution
        if key in modified_reads:
            p_out.stdin.write(modified_reads[key])
        else:
            p_out.stdin.write(line)

    # 4. Clean up and wait for processes to finish cleanly
    p_in.stdout.close()
    p_out.stdin.close()
    p_in.wait()
    p_out.wait()
    
    if p_in.returncode != 0 or p_out.returncode != 0:
        sys.exit("Error: A samtools subprocess failed during the stream.")

def main(args):
    check_dependencies()
    
    if args.time_run:
        start_time = time.time()
        print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Ensure output directory exists
    out_dir = os.path.dirname(args.out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Step 1: Load the modified dictionary
    modified_dict = load_modified_reads(args.modified_bam)
    
    # Step 2: Run the fast stream replacement
    stream_and_replace(args.input_bam, args.out_file, modified_dict)
    
    # Step 3: Index the newly created, inherently sorted BAM
    print("Indexing final BAM...")
    subprocess.run(["samtools", "index", args.out_file], check=True)

    if args.time_run:
        end_time = time.time()
        print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Total run time: {end_time - start_time:.2f} seconds")
    
    print("Done.")

if __name__ == '__main__':
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Fast BAM read replacer using subprocess streaming."
    )
    parser.add_argument("-b", "--input_bam", required=True, help="Unmodified input bam file")
    parser.add_argument("-m", "--modified_bam", required=True, help="BAM file with modified reads only")
    parser.add_argument("-o", "--out_file", help="Output file name", default=os.path.join(os.getcwd(), "tweakvar_modified.bam"))
    parser.add_argument("--time_run", action="store_true", default=False, help="Export start time, end time, and total run time")
    
    args = parser.parse_args()
    main(args)
