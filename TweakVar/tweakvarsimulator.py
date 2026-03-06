import sys
import pysam
from pysam import depth
import argparse
import os
import numpy as np
from numpy.random import choice, uniform, seed as npseed
from numpy import random as nran
from math import ceil
import random

# Initialize parser
parser = argparse.ArgumentParser(description="TweakVarSimulator: A tool for simulating variants")

# Required arguments
parser.add_argument("-i", "--input", required=True, help="Path to input BAM file or txt file containing BAMs listed from lowest to highest coverage. If you are using a text file, ensure that the bam files are listed from top to bottom in increasing file size.")
parser.add_argument("-T", "--reference", required=True, help="Path to reference genome")
parser.add_argument("-o", "--output", required=True, help="Output file prefix or txt file containing prefixes")

# Optional arguments
parser.add_argument("-s", "--seed", type=int, required=False, help="Seed for random number generation")
parser.add_argument("-minAFsv", "--minimum_allele_frequency_sv", type=float, required=False)
parser.add_argument("-maxAFsv", "--maximum_allele_frequency_sv", type=float, required=False)
parser.add_argument("-minAFsnv", "--minimum_allele_frequency_snv", type=float, required=False)
parser.add_argument("-maxAFsnv", "--maximum_allele_frequency_snv", type=float, required=False)
parser.add_argument("-numsv", "--number_of_svs", type=int, required=False)
parser.add_argument("-numsnv", "--number_of_snvs", type=int, required=False)
parser.add_argument("-maxsnvl", "--maximum_snv_length", type=int, required=False)
parser.add_argument("-minsvl", "--minimum_sv_length", type=int, required=False)
parser.add_argument("-maxsvl", "--maximum_sv_length", type=int, required=False)
parser.add_argument("-sub", "--substitution_rate", type=float, required=False)
parser.add_argument("-insdelsnv", "--insdel_snv_rate", type=float, required=False)
parser.add_argument("-insdel", "--insdel_sv_rate", type=float, required=False)
parser.add_argument("--ref_chroms", nargs='+', help="Optional list of reference chromosomes")
parser.add_argument("-snv", "--SNV_truth_file", required=False, help="Optional SNV truth VCF file")
parser.add_argument("-sv", "--SV_truth_file", required=False, help="Optional SV truth VCF file")
parser.add_argument("-bps", "--bp_shift", type=int,  required=False, help="Optionally allows you to shift all bp positions up or downstream by the indicated number of bps")
parser.add_argument("-imc", "--ignore_minimum_cov", type=bool, required=False, help="If True, loci will not be evaluated for minimum coverage")

args = parser.parse_args()

### bam_paths length determination
# determine if input is a single bam file or txt file
input_files = args.input

bam_paths = []
file_sizes = []

if input_files.endswith(".bam"):
    bam_paths.append(input_files)
    file_sizes.append(os.path.getsize(input_files))

elif input_files.endswith(".txt"):
    with open(input_files, 'r') as f:
        for line in f:
            path = line.strip()
            if path:  # Skip empty lines
                bam_paths.append(path)
                # Check if file exists before getting size to avoid crashing
                if os.path.exists(path):
                    file_sizes.append(os.path.getsize(path))
                else:
                    print(f"Warning: File {path} not found.")

    # Check for increasing order
    # Using 'all' with a generator is a very "Pythonic" way to check sorting
    is_sorted = all(file_sizes[i] <= file_sizes[i+1] for i in range(len(file_sizes)-1))
    
    if not is_sorted:
        print("Files are not ordered by increasing file size. "
              "Note: this may affect time for variant simulation.")

else:
    print("No .bam or .txt file detected. Check your extensions.")

# ensure all references align in the bam files to the reference file
ref_path = args.reference

# assign lengths to each sequence name in fasta file
fasta_contigs = {}
with pysam.FastaFile(ref_path) as fasta:
    # fasta.references and fasta.lengths are lists of names and lengths
    for name, length in zip(fasta.references, fasta.lengths):
        fasta_contigs[name] = length  # Use [] for assignment

# extract reference sequence information from bam files
for bam_path in bam_paths:
    bam_contigs = {}
    print(f"Checking alignment of {bam_path} with reference...")
    
    # Pass the variable bam_path, not the string "bam_path"
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # bam.header.get("SQ") returns a list of dictionaries
        for sq in bam.header.get("SQ", []):
            bam_contigs[sq["SN"]] = sq["LN"]

    # check if both dictionaries match (keys and values)
    if bam_contigs != fasta_contigs:
        print(f"Error: {bam_path} header does not match the reference sequence.", file=sys.stderr)
        sys.exit(1)

print("Proceeding...")

output_prefix = args.output
SVvcf_paths = []
SNVvcf_paths = []

if input_files.endswith(".txt"):
    # Open the input list, not the output prefix
    with open(input_files, 'r') as f:
        # Load lines and remove whitespace/newlines
        input_lines = [line.strip() for line in f if line.strip()]
    
    # Check if the number of rows matches bam_paths
    if len(input_lines) != len(bam_paths):
        raise ValueError(f"Mismatch: {len(input_lines)} inputs found, but {len(bam_paths)} BAM paths provided.")

    for item in input_lines:
        SVvcf_paths.append(f"{item}_SV.vcf")
        SNVvcf_paths.append(f"{item}_SNV.vcf")
        
else:
    # Single file mode
    SVvcf_paths.append(f"{output_prefix}_SV.vcf")
    SNVvcf_paths.append(f"{output_prefix}_SNV.vcf")

snv_truth_file = args.SNV_truth_file
sv_truth_file = args.SV_truth_file

bp_shift = args.bp_shift
ignore_minimum_cov = args.ignore_minimum_cov

if ignore_minimum_cov is None:
    ignore_minimum_cov = False

if bp_shift is None:
    bp_shift = 0
else:
    bp_shift = int(bp_shift)

seed = args.seed if args.seed is not None else 0
npseed(seed)
random.seed(seed)

minAFsv = args.minimum_allele_frequency_sv if args.minimum_allele_frequency_sv is not None else 0.01
maxAFsv = args.maximum_allele_frequency_sv if args.maximum_allele_frequency_sv is not None else 0.05
numsv = args.number_of_svs if args.number_of_svs is not None else 50
minsvl = args.minimum_sv_length if args.minimum_sv_length is not None else 50
maxsvl = args.maximum_sv_length if args.maximum_sv_length is not None else 10000
insdel = args.insdel_sv_rate if args.insdel_sv_rate is not None else 0.7

minAFsnv = args.minimum_allele_frequency_snv if args.minimum_allele_frequency_snv is not None else 0.01
maxAFsnv = args.maximum_allele_frequency_snv if args.maximum_allele_frequency_snv is not None else 0.05
numsnv = args.number_of_snvs if args.number_of_snvs is not None else 200
maxsnvl = args.maximum_snv_length if args.maximum_snv_length is not None else 100
sub = args.substitution_rate if args.substitution_rate is not None else 1
insdelsnv = args.insdel_snv_rate if args.insdel_snv_rate is not None else 0.5


def count_decimals(x):
    s = str(x)
    if '.' in s:
        return len(s.split('.')[-1].rstrip('0'))
    else:
        return 0

def parse_truth_vcf(vcf_path):
    records = []
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, _, ref, alt, _, _, info, *_ = parts
            af_val = 0.05
            for field in info.split(';'):
                if field.startswith('AF='):
                    try:
                        af_val = float(field.split('=')[1])
                    except:
                        pass
            records.append((chrom, pos, ref, alt, af_val))
    return records

def get_chrom_lengths(bam_path, ref_chroms=None):
    with pysam.AlignmentFile(bam_path) as bam:
        chromol = dict(zip(bam.references, bam.lengths))

    if ref_chroms is not None:
        # Filter only if ref_chroms is specified
        chromol = {c: l for c, l in chromol.items() if c in ref_chroms}

    return chromol, tuple(chromol.keys())

def genlocSNV(num, bam_paths, mincov=20):
    # Assuming npseed and nran are aliases for numpy.random
    npseed(seed)
    ref_chroms = args.ref_chroms
    
    # Load reference fasta file
    ref_fasta = pysam.FastaFile(ref_path)
    
    # We can use the first BAM to get the chromosome lengths
    chromol, chrom = get_chrom_lengths(bam_paths[0], ref_chroms=ref_chroms)
    
    locations = []
    
    for i in range(num):
        print(f"Creating SNV {i+1}/{num}...")
        
        while True:
            ranchrom = nran.choice(chrom)
            loc_int = nran.randint(0, chromol[ranchrom]) + bp_shift
            loc = str(loc_int)
            
            # Fetch the reference base
            ref_base = ref_fasta.fetch(ranchrom, loc_int - 1, loc_int).upper()
            if ref_base not in ['A', 'C', 'G', 'T']:
                continue  # It's an N or a gap, try a new location
                
            # --- NEW: Get coverage for ALL BAM files at this locus ---
            covers = []
            for bam_path in bam_paths:
                try:
                    cover_val = int(depth(bam_path, '-r', f"{ranchrom}:{loc}-{loc}").rstrip("\n").split("\t")[-1])
                except (ValueError, IndexError):
                    cover_val = 0
                covers.append(cover_val)
            
            # Test for minimum coverage requirement
            if ignore_minimum_cov:
                break
            else:
                # Require ALL BAM files to meet the minimum coverage.
                # Change 'all' to 'any' if you only need at least one BAM to pass.
                if all(c >= mincov for c in covers):
                    break
                    
        # Append chrom, loc, the ref_base, and the list of coverages
        locations.append((ranchrom, loc, ref_base, covers))
        
    return locations

def genlocSV(num, bam_paths, mincov=20):
    npseed(seed)
    ref_chroms = args.ref_chroms
    # Use the first BAM to get chrom lengths
    chromol, chrom = get_chrom_lengths(bam_paths[0], ref_chroms=ref_chroms)
    
    locations = []
    for i in range(num):
        print(f"Creating SV {i+1}/{num}...")
        while True:
            ranchrom = nran.choice(chrom)
            loc = str(nran.randint(0, chromol[ranchrom]) + bp_shift)
            
            # Check for overlapping/nearby SVs 
            # (Changed loop variable from 'i' to 'existing' to avoid overwriting the outer loop variable)
            too_close = False
            for existing in locations:
                if existing[0] == ranchrom and abs(int(existing[1])) - int(loc) < 50000:
                    too_close = True
                    break
            
            if too_close:
                continue
                
            # --- Get coverage for ALL BAM files ---
            covers = []
            for bam_path in bam_paths:
                try:
                    cover_val = int(depth(bam_path, '-r', f"{ranchrom}:{loc}-{loc}").rstrip("\n").split("\t")[-1])
                except (ValueError, IndexError):
                    cover_val = 0
                covers.append(cover_val)
                
            ## Test for minimum coverage requirement
            if ignore_minimum_cov:
                break
            else:
                # Require ALL BAM files to meet the minimum coverage
                if all(c >= mincov for c in covers):
                    break
                    
        # Append chrom, loc, and the list of coverages
        locations.append((ranchrom, loc, covers))
        
    return tuple(locations)

def genseq(minl,maxl):
    npseed(seed)
    nucl=tuple(["A","T","C","G"])
    res=''
    for i in range(choice(range(minl,maxl+1))):
        res+=choice(nucl)
    return res

def gensnps(maxsnvl,sub,insdelsnv=insdelsnv, snplist=[]):
    npseed(seed)
    result=[]
    nucl=["A","T","C","G"]
    for i in range(len(snplist)):
        draw = choice(tuple(['snp','indel']), 1, p= [sub,1-sub])    
        if draw=='snp':
            nucl=["A","T","C","G"]
            nucl.remove(snplist[i][3][0].upper())
            result.append(tuple(list(snplist[i][0:3])+[snplist[i][3][0].upper()]+[choice(nucl,1)[0]]))
        else:
            draw = choice(tuple(['in','del'],), 1, p= [insdelsnv,1-insdelsnv])    
            if draw=='in':
                result.append(tuple(list(snplist[i][0:3])+[snplist[i][3][0].upper()]+[snplist[i][3][0].upper()+genseq(1,maxsnvl-1).upper()]))
            else:
                rmlen=choice(range(maxsnvl+1)[1:])
                result.append(tuple(list(snplist[i][0:3])+[snplist[i][3].upper()[:rmlen]]+[snplist[i][3].upper()[0]]))

    return tuple(result)
def getrefsnp(reffile,snplist=1):
    with open(reffile,'r') as f:
        data=[]
        cta=[]
        for i in snplist:
            cta.append(i[0])
        cta=set(cta)
        start=1
        for line in f:
            if line.startswith(">"):
                if start==1:
                    start=0
                else:
                    if analyse==1:
                        data.append(tuple([chromosome,seq]))
                if "chr"==snplist[0][0][0:3]:
                    chromosome=line.lstrip(">").split(" ")[0].rstrip("\n")
                else:
                    chromosome=line.lstrip(">").split(" ")[0].lstrip("chr").rstrip("\n")
                seq=''
                analyse=0
                if chromosome in cta:
                    analyse=1
            else:
                if analyse==0:
                    continue
                seq+=line.rstrip("\n") 
        if analyse==1:
            data.append(tuple([chromosome,seq]))
        data=tuple(data)

    f.close()
    snplistmod=[]
    for i in snplist:
        for j in data:
            if i[0]==j[0]:
                snplistmod.append(tuple(list(i)+[j[1][int(i[1])-1:int(i[1])+maxsnvl]]))
                break
    return tuple(snplistmod)
    


def genseq(minl,maxl):
    npseed(seed)
    nucl=tuple(["A","T","C","G"])
    res=''
    for i in range(choice(range(minl,maxl+1))):
        res+=choice(nucl)
    return res

def main():
    # Only proceed if we have a truth file
    if snv_truth_file is not None:
        raw_snvloc = parse_truth_vcf(snv_truth_file)

        # Loop over every BAM file
        for idx, bam_path in enumerate(bam_paths):
            print(f"\n--- Processing BAM {idx+1}/{len(bam_paths)}: {bam_path} ---")
            
            # 1. Determine the correct output VCF path
            if len(SNVvcf_paths) > 1:
                # Use the directly corresponding path if a list was provided
                current_SNVvcf = SNVvcf_paths[idx]
            else:
                # If only one path exists, enumerate it based on the BAM index
                base_path = SNVvcf_paths[0]
                if base_path.endswith('.vcf'):
                    current_SNVvcf = f"{base_path[:-4]}_{idx+1}.vcf"
                else:
                    current_SNVvcf = f"{base_path}_{idx+1}.vcf"

            # 2. Get chromosome lengths for the current BAM
            chromol, chrom = get_chrom_lengths(bam_path)

            # 3. Initialize fresh VCF headers for this specific BAM iteration
            vcfsnv = [
                '##fileformat=VCFv4.2',
                '##FILTER=<ID=PASS,Description="All filters passed">',
                '##FILTER=<ID=FAIL,Description="Failed minimum coverage">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
                '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">',
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Target Allele Frequency">',
                '##INFO=<ID=EAF,Number=A,Type=Float,Description="Exact Allele Frequency (variant reads / total depth)">'
            ]
            vcfsnv.extend([f'##contig=<ID={c},length={l}>' for c, l in chromol.items()])
            vcfsnv.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))

            vcfsnv_failed = vcfsnv.copy()
            header_len = len(vcfsnv) 

            # 4. Process each SNV location
            for i in range(len(raw_snvloc)):
                # Optional: If you have thousands of SNVs, consider printing this only every 100 or 1000 iterations to avoid console spam.
                print(f"Retrieving SNV {i+1}/{len(raw_snvloc)} from truth vcf...")
                
                chrom, pos, ref, alt, af_val = raw_snvloc[i]
                
                # update the position
                pos = str(int(pos) + bp_shift)
                
                # Failsafe coverage calculation using the current bam_path
                try:
                    cover = int(depth(bam_path, '-r', f"{chrom}:{pos}-{pos}").rstrip("\n").split("\t")[-1])
                except (ValueError, IndexError): # Added IndexError just in case depth() returns empty
                    cover = 0 
                    
                pos = int(pos)
                
                # Test for minimum coverage requirement
                if not ignore_minimum_cov and cover < mincov:
                    vcfsnv_failed.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1500\tFAIL\tAF={af_val};EAF=0.0\tGT:AD:DV\t0/0:{cover}:0")
                else:
                    from math import ceil # Ensure ceil is imported
                    readnum = ceil(af_val * cover)
                    exact_af = round(readnum / cover, 4) if cover > 0 else 0.0
                    
                    vcfsnv.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1500\tPASS\tAF={af_val};EAF={exact_af}\tGT:AD:DV\t0/0:{cover-readnum}:{readnum}")

            # 5. Write the passing VCF for this BAM
            with open(current_SNVvcf, "w") as f:
                f.write('\n'.join(vcfsnv) + '\n')
            print(f"New SNV truth vcf written to {current_SNVvcf}.")

            # 6. Write the failed VCF for this BAM (if any failed)
            if len(vcfsnv_failed) > header_len:
                if current_SNVvcf.endswith('.vcf'):
                    failed_SNVvcf = f"{current_SNVvcf[:-4]}_failed.vcf"
                else:
                    failed_SNVvcf = f"{current_SNVvcf}_failed.vcf"
                    
                with open(failed_SNVvcf, "w") as f:
                    f.write('\n'.join(vcfsnv_failed) + '\n')
                    
                num_failed = len(vcfsnv_failed) - header_len
                print(f"Found {num_failed} SNVs below minimum coverage. Written to {failed_SNVvcf}.")          
else:
    random.seed(seed)
    
    # 1. Generate the random locations (which now includes a list of covers for all BAMs)
    snvloc = genlocSNV(numsnv, bam_paths, ceil(1/maxAFsnv))

    print("Writing the SNV output files...")

    if numsnv > 0:
        snps = getrefsnp(ref_path, snvloc)
        
        # 2. Materialize the variants FIRST so they are consistent across all BAM VCFs
        simulated_variants = []
        for i in gensnps(maxsnvl=maxsnvl, sub=sub, snplist=snps):
            # i[0]=chrom, i[1]=pos, i[2]=list of covers, i[3]=ref, i[4]=alt
            AFnum = round(random.uniform(minAFsnv, maxAFsnv), count_decimals(minAFsnv))
            simulated_variants.append({
                'chrom': i[0],
                'pos': i[1],
                'covers': i[2], # This is the list of coverages from genlocSNV
                'ref': i[3],
                'alt': i[4],
                'target_af': AFnum
            })

        # 3. Loop over each BAM file to create its unique VCF
        for idx, bam_path in enumerate(bam_paths):
            
            # Determine the correct output VCF path
            if len(SNVvcf_paths) > 1:
                current_SNVvcf = SNVvcf_paths[idx]
            else:
                base_path = SNVvcf_paths[0]
                if base_path.endswith('.vcf'):
                    current_SNVvcf = f"{base_path[:-4]}_{idx+1}.vcf"
                else:
                    current_SNVvcf = f"{base_path}_{idx+1}.vcf"

            # Initialize headers for this specific VCF
            vcfsnv = [
                '##fileformat=VCFv4.2',
                '##FILTER=<ID=PASS,Description="All filters passed">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
                '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">',
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Target Allele Frequency">',
                '##INFO=<ID=EAF,Number=A,Type=Float,Description="Exact Allele Frequency (variant reads / total depth)">'
            ]

            # Dynamically add contig lines (assuming chromol was assigned earlier in main)
            vcfsnv.extend([f'##contig=<ID={c},length={l}>' for c, l in chromol.items()])
            vcfsnv.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))

            # Process the pre-generated variants using the specific coverage for this BAM
            for var in simulated_variants:
                # Grab the coverage matching the current BAM index
                cover = var['covers'][idx] 
                readnum = ceil(var['target_af'] * cover)
                
                # Calculate exact AF (prevent division by zero)
                exact_af = round(readnum / cover, 4) if cover > 0 else 0.0
                
                vcfsnv.append(f"{var['chrom']}\t{var['pos']}\t.\t{var['ref']}\t{var['alt']}\t1500\tPASS\tAF={var['target_af']};EAF={exact_af}\tGT:AD:DV\t0/0:{cover-readnum}:{readnum}")
            
            # Write the file for this BAM
            with open(current_SNVvcf, "w") as f:
                f.write('\n'.join(vcfsnv) + '\n')

            print(f"New SNV truth vcf for BAM {idx+1} written to {current_SNVvcf}.")

    # SV processing
    if sv_truth_file is not None:
        print(f"Retrieving SVs from truth vcf...")
        svloc = parse_truth_vcf(sv_truth_file)
        
        # Loop over every BAM file for SV truth parsing
        for idx, bam_path in enumerate(bam_paths):
            print(f"\n--- Processing BAM {idx+1}/{len(bam_paths)} for SVs: {bam_path} ---")
            
            # 1. Determine the correct output VCF path
            if len(SVvcf_paths) > 1:
                current_SVvcf = SVvcf_paths[idx]
            else:
                base_path = SVvcf_paths[0]
                if base_path.endswith('.vcf'):
                    current_SVvcf = f"{base_path[:-4]}_{idx+1}.vcf"
                else:
                    current_SVvcf = f"{base_path}_{idx+1}.vcf"

            chromol, chrom = get_chrom_lengths(bam_path)
            
            vcfsv = [
                '##fileformat=VCFv4.2',
                '##ALT=<ID=INS,Description="Insertion">',
                '##ALT=<ID=DEL,Description="Deletion">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">', 
                '##FILTER=<ID=PASS,Description="All filters passed">',
                '##FILTER=<ID=FAIL,Description="Failed minimum coverage">',
                '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">',
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">',
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">',
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">',
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'
            ]
            vcfsv.extend([f'##contig=<ID={c},length={l}>' for c, l in chromol.items()])
            vcfsv.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))

            vcfsv_failed = vcfsv.copy()
            header_len = len(vcfsv)

            for rec_idx, record in enumerate(svloc):
                chrom, pos, ref, alt, af = record
                pos = str(int(pos) + bp_shift)

                ## Failsafe coverage calculation tied to specific BAM
                try:
                    cover = int(depth(bam_path, '-r', f"{chrom}:{pos}-{pos}").rstrip("\n").split("\t")[-1])
                except (ValueError, IndexError):
                    cover = 0
                
                svtype = 'INS' if len(alt) > len(ref) else 'DEL'
                svlen = len(alt) - len(ref) if svtype == 'INS' else -(len(ref) - len(alt))
                end = int(pos) + abs(svlen)
                ID = f"SV{rec_idx+1}"

                ## Test for minimum coverage requirement
                if not ignore_minimum_cov and cover < mincov:
                    vcfsv_failed.append(f"{chrom}\t{pos}\t{ID}\t{ref}\t{alt}\t60\tFAIL\tPRECISE;SVTYPE={svtype};SVLEN={svlen};END={end};AF={af}\tGT:GQ:DP\t0/0:60:{cover}")
                else:
                    vcfsv.append(f"{chrom}\t{pos}\t{ID}\t{ref}\t{alt}\t60\tPASS\tPRECISE;SVTYPE={svtype};SVLEN={svlen};END={end};AF={af}\tGT:GQ:DP\t0/0:60:{cover}")

            # --- File Writing Logic ---
            with open(current_SVvcf, "w") as f:
                f.write('\n'.join(vcfsv) + '\n')

            print(f"New SV truth vcf written to {current_SVvcf}.")

            if len(vcfsv_failed) > header_len:
                if current_SVvcf.endswith('.vcf'):
                    failed_SVvcf = f"{current_SVvcf[:-4]}_failed.vcf"
                else:
                    failed_SVvcf = f"{current_SVvcf}_failed.vcf"
                    
                with open(failed_SVvcf, "w") as f:
                    f.write('\n'.join(vcfsv_failed) + '\n')
                    
                num_failed = len(vcfsv_failed) - header_len
                print(f"Found {num_failed} SVs below minimum coverage. Written to {failed_SVvcf}.")

    else:
        # Pass the list of BAMs to genlocSV
        svloc = genlocSV(numsv, bam_paths, ceil(1/minAFsv))
        random.seed(seed)

        print("Writing the SV output files...")
        if numsv > 0:
            # 1. Materialize SV properties FIRST so they map perfectly across all VCFs
            simulated_svs = []
            insertnum = 1
            delnum = 1
            
            for i in svloc:
                # i[0] = chrom, i[1] = pos, i[2] = list of covers
                draw = choice(tuple(['in','del']), 1, p=[insdel, 1-insdel]) 
                af = round(random.uniform(minAFsv, maxAFsv), count_decimals(minAFsv))
                
                if draw == 'in':
                    seq = genseq(minsvl, maxsvl)
                    simulated_svs.append({
                        'chrom': i[0], 'pos': i[1], 'covers': i[2],
                        'id': f"HackIns{insertnum}", 'ref': 'N', 'alt': seq,
                        'svtype': 'INS', 'svlen': len(seq), 'end': int(i[1]) + 1, 'af': af
                    })
                    insertnum += 1
                else:
                    dellen = choice(range(minsvl, maxsvl))
                    simulated_svs.append({
                        'chrom': i[0], 'pos': i[1], 'covers': i[2],
                        'id': f"HackDel{delnum}", 'ref': 'N', 'alt': '<DEL>',
                        'svtype': 'DEL', 'svlen': -dellen, 'end': int(i[1]) + dellen, 'af': af
                    })
                    delnum += 1

            # 2. Generate a VCF for each BAM file
            for idx, bam_path in enumerate(bam_paths):
                if len(SVvcf_paths) > 1:
                    current_SVvcf = SVvcf_paths[idx]
                else:
                    base_path = SVvcf_paths[0]
                    if base_path.endswith('.vcf'):
                        current_SVvcf = f"{base_path[:-4]}_{idx+1}.vcf"
                    else:
                        current_SVvcf = f"{base_path}_{idx+1}.vcf"

                chromol, chrom = get_chrom_lengths(bam_path)
                
                vcfsv = [
                    '##fileformat=VCFv4.2',
                    '##ALT=<ID=INS,Description="Insertion">',
                    '##ALT=<ID=DEL,Description="Deletion">',
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
                    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">', 
                    '##FILTER=<ID=PASS,Description="All filters passed">',
                    '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">',
                    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">',
                    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">',
                    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">',
                    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'
                ]

                vcfsv.extend([f'##contig=<ID={c},length={l}>' for c, l in chromol.items()])
                vcfsv.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))

                # Map materialized variants to VCF format string
                for var in simulated_svs:
                    # Fetch specific coverage for this BAM
                    cover = var['covers'][idx]
                    vcfsv.append(
                        f"{var['chrom']}\t{var['pos']}\t{var['id']}\t{var['ref']}\t{var['alt']}\t60\tPASS\t"
                        f"PRECISE;SVTYPE={var['svtype']};SVLEN={var['svlen']};END={var['end']};AF={var['af']}\t"
                        f"GT:GQ:DP\t0/0:60:{cover}"
                    )

                # Write out final files
                # Note: changed output_prefix to avoid errors if missing
                os.makedirs(os.path.dirname(current_SVvcf) or '.', exist_ok=True)

                with open(current_SVvcf, "w") as f:
                    f.write('\n'.join(vcfsv) + '\n')
                
                print(f"New SV truth vcf for BAM {idx+1} written to {current_SVvcf}.")
        

if __name__ == "__main__":
    main()
