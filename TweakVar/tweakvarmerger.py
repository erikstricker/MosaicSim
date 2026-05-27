#!/usr/bin/env python3

import sys
import os
import time
import shutil
import subprocess
from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def check_dependencies():
    """Ensure samtools is installed and available in the system PATH."""
    if shutil.which("samtools") is None:
        sys.exit("Error: 'samtools' is not installed or not in your PATH. It is required for fast streaming.")

def clean_qname(qname):
    """Strips any paired-end suffixes that simulators occasionally append."""
    if qname.endswith('/1') or qname.endswith('/2'):
        return qname[:-2]
    return qname

def get_next_mod(stdout_stream, total_count_wrapper):
    """Fetches the next available modified read from the stream in O(1) time."""
    line = stdout_stream.readline()
    if not line:
        return None
        
    total_count_wrapper[0] += 1  # Track total reads present in the modified file
    cols = line.rstrip('\n').split('\t')
    
    qname = cols[0]
    if qname.endswith('/1') or qname.endswith('/2'):
        qname = qname[:-2]
        
    return {
        'qname': qname,
        'chrom': cols[2],  # Column 3 is RNAME (Chrom)
        'pos': cols[3],    # Column 4 is POS
        'cigar': cols[5],  # Column 6 is CIGAR
        'seq': cols[9],    # Column 10 is SEQ
        'qual': cols[10]   # Column 11 is QUAL (Required for Indel length matching)
    }

def stream_and_replace(input_bam, modified_bam, output_bam):
    """Coordinates parallel sorted streams using an O(1) synchronized tracking pointer."""
    print("Streaming and substituting reads via synchronized parallel pipeline...")
    
    # 1. Decompress the original input BAM
    p_in = subprocess.Popen(
        ["samtools", "view", "-h", input_bam],
        stdout=subprocess.PIPE,
        text=True
    )
    
    # 2. Stream the pre-sorted modified BAM on a parallel track (no headers)
    p_mod = subprocess.Popen(
        ["samtools", "view", modified_bam],
        stdout=subprocess.PIPE,
        text=True
    )
    
    # 3. Compress directly to final destination
    p_out = subprocess.Popen(
        ["samtools", "view", "-b", "-o", output_bam],
        stdin=subprocess.PIPE,
        text=True
    )
    
    total_mod_count = [0]
    replaced_count = 0
    
    # Initialize pointer to the first modified read
    curr_mod = get_next_mod(p_mod.stdout, total_mod_count)
    
    # Main loop over the original BAM stream
    for line in p_in.stdout:
        if line[0] == '@':
            p_out.stdin.write(line)
            continue
            
        # FAST-PATH: Split only the first 4 columns to get coordinates
        # This completely avoids parsing the massive SEQ and QUAL strings for unaltered reads
        cols = line.split('\t', 4)
        qname_raw = cols[0]
        chrom = cols[2]
        pos = cols[3]
        
        if qname_raw.endswith('/1') or qname_raw.endswith('/2'):
            qname = qname_raw[:-2]
        else:
            qname = qname_raw
            
        # INTRA-CHROMOSOME CATCH-UP: If the modified stream pointer lags behind the original position
        # on the same chromosome, advance it. This eliminates any chance of a pointer lock.
        while curr_mod and chrom == curr_mod['chrom'] and int(pos) > int(curr_mod['pos']):
            curr_mod = get_next_mod(p_mod.stdout, total_mod_count)
            
        # If positions and names match perfectly, execute substitution
        if curr_mod and chrom == curr_mod['chrom'] and pos == curr_mod['pos'] and qname == curr_mod['qname']:
            # Slow Path: Only triggered for target modifications
            orig_cols = line.rstrip('\n').split('\t')
            
            # FIXED: Overwrite CIGAR (5), SEQ (9), and QUAL (10) to support Indels cleanly
            orig_cols[5] = curr_mod['cigar']
            orig_cols[9] = curr_mod['seq']
            orig_cols[10] = curr_mod['qual']
            
            p_out.stdin.write('\t'.join(orig_cols) + '\n')
            replaced_count += 1
            
            # Advance the modified pointer
            curr_mod = get_next_mod(p_mod.stdout, total_mod_count)
        else:
            # FAST-PATH: Keep original read entirely intact
            p_out.stdin.write(line)

    # Exhaust remaining lines in modified stream if any are left to ensure total file count is exact
    while curr_mod:
        curr_mod = get_next_mod(p_mod.stdout, total_mod_count)

    # Clean up streams
    p_in.stdout.close()
    p_mod.stdout.close()
    p_out.stdin.close()
    
    p_in.wait()
    p_mod.wait()
    p_out.wait()
    
    print(f"Total reads available in modified BAM: {total_mod_count[0]}")
    print(f"Total reads successfully replaced: {replaced_count}")
    
    if p_in.returncode != 0 or p_out.returncode != 0:
        sys.exit("Error: A samtools subprocess failed during the stream.")

def main(args):
    check_dependencies()
    
    if args.time_run:
        start_time = time.time()
        print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    out_dir = os.path.dirname(args.out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    stream_and_replace(args.input_bam, args.modified_bam, args.out_file)
    
    print("Indexing final BAM...")
    subprocess.run(["samtools", "index", args.out_file], check=True)

    if args.time_run:
        end_time = time.time()
        print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        total_seconds = end_time - start_time
        hours = int(total_seconds // 3600)
        minutes = int((total_seconds % 3600) // 60)
        seconds = total_seconds % 60
        
        print(f"Total run time: {hours:02d}:{minutes:02d}:{seconds:05.2f}")
    
    print("Done.")

if __name__ == '__main__':
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Fast BAM read replacer using sorted synchronized parallel streams."
    )
    parser.add_argument("-b", "--input_bam", required=True, help="Unmodified input bam file")
    parser.add_argument("-m", "--modified_bam", required=True, help="BAM file with modified reads only")
    parser.add_argument("-o", "--out_file", help="Output file name", default=os.path.join(os.getcwd(), "tweakvar_modified.bam"))
    parser.add_argument("--time_run", action="store_true", default=False, help="Export start time, end time, and total run time")
    
    args = parser.parse_args()
    main(args)
