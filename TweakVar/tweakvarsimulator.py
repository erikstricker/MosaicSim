import os
import sys
import time
import argparse
import numpy
import pysam
from datetime import datetime
from numpy import random
from numpy.random import choice, randint, seed
from pysam import depth, AlignmentFile, FastaFile
from math import ceil
import random

MAX_ATTEMPTS = 100000
SV_MIN_GAP = 50000


def str2bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in ("true", "t", "yes", "y", "1"):
        return True
    if value.lower() in ("false", "f", "no", "n", "0"):
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got '{value}'.")

# Initialize parser
parser = argparse.ArgumentParser(description="TweakVarSimulator: A tool for simulating variants")

# Required arguments
parser.add_argument("-i", "--input", required=True, help="Required. Path to BAM file or txt file containing BAM (BAMs in a txt file should be listed from lowest to highest coverage)")
parser.add_argument("-T", "--reference", required=True, help="Required. Path to reference genome")
parser.add_argument("-o", "--output", required=True, help="Required. Output file prefix, or a .txt file with one output prefix per BAM (same order as the BAM list). If a multi-BAM .txt input is used with a plain prefix, SNV vcf files will be generated in the same directory as each of the BAM files.")

# Optional arguments
parser.add_argument("-s", "--seed", type=int, required=False, help="Optional. Random seed for reproducibility. If not set, the results will vary between runs. The same seed with different input but identical parameters will lead to the same list of loci.")
parser.add_argument("-minAFsv", "--minimum_allele_frequency_sv", type=float, required=False, help="Optional. Minimum allele frequency for simulated structural variants (SVs).")
parser.add_argument("-maxAFsv", "--maximum_allele_frequency_sv", type=float, required=False, help="Optional. Maximum allele frequency for simulated structural variants (SVs).")
parser.add_argument("-minAFsnv", "--minimum_allele_frequency_snv", type=float, required=False, help="Optional. Minimum allele frequency for simulated single nucleotide variants (SNVs). May also be used to create a minimum value for the allele frequency when generating variants with a SNV truth file.")
parser.add_argument("-maxAFsnv", "--maximum_allele_frequency_snv", type=float, required=False, help="Optional. Maxiumum allele frequency for simulated single nucleotide variants (SNVs.) May also be used to limit the allele frequency when generating variants with a SNV truth file.")
parser.add_argument("-numsv", "--number_of_svs", type=int, required=False, help="Optional. Number of structual variants (SVs) to simulate.")
parser.add_argument("-numsv_align", "--number_of_svs_to_align", type=int, required=False, help="Optional. Number of structural variants (SVs) to match the same positions and AF as the truth dataset.")
parser.add_argument("-numsnv", "--number_of_snvs", type=int, required=False, help="Optional. Number of single nucleotide variants (SNVs) to simulate.")
parser.add_argument("-numsnv_align", "--number_of_snvs_to_align", type=int, required=False, help="Optional. Number of single nucleotide variants (SNVs) to match the same positions and AF as the truth dataset.")
parser.add_argument("-maxsnvl", "--maximum_snv_length", type=int, required=False)
parser.add_argument("-minsvl", "--minimum_sv_length", type=int, required=False, help="Optional. Minimum length of simulated structural variants (SVs) (in base pairs).")
parser.add_argument("-maxsvl", "--maximum_sv_length", type=int, required=False, help="Optional. Maximum length of simulated structural variants (SVs) (in base pairs).")
parser.add_argument("-sub", "--substitution_rate", type=float, required=False, help="Optional. Probability of generating a SNP versus an indel (in SNVs). A value of 1.0 means only SNPs will be generated.")
parser.add_argument("-insdelsnv", "--insdel_snv_rate", type=float, required=False, help="Optional. Probability of generating an insertion vs. a deletion in SNVs. A value of 0.5 means equal chances of insertion and deletion.")
parser.add_argument("-insdel", "--insdel_sv_rate", type=float, required=False, help="Optional. Probability of generating an insertion vs. a deletion in SVs. A value of 0.5 means equal chances of insertion and deletion.")
parser.add_argument("--ref_chroms", nargs='+', help="Optional. List of chromosomes to include in the simulation (e.g., chr1 chr2 ... chrY). Only these chromosomes will be used for variant generation. may also be used for limiting the chromosomes to match to a provided SNV truth dataset.")
parser.add_argument("-snv", "--SNV_truth_file", required=False, help="Optional. SNV truth file of which variant locations, allele frequencies, and alternate alleles are used as the foundation to create a new SNV VCF file with updated coverage and read count information.")
parser.add_argument("-sv", "--SV_truth_file", required=False, help="Optional. SV truth file of which variant locations, allele frequencies, and alternate alleles are used as the foundation to create a new SV VCF file with updated coverage and read count information.")
parser.add_argument("-bps", "--bp_shift", type=int,  required=False, help="Optional. Allows the shift of all bp positions up- or downstream by the indicated number of bps.")
parser.add_argument("-imc", "--ignore_minimum_cov", type=str2bool, nargs="?", const=True, default=False, help="Optional. Ignores minimum coverage requirement for the selection of variant loci. Useful when selecting identical regions in down-sampled files. Use as a flag (-imc) or with a value (-imc True / -imc False).")

args = parser.parse_args()

## Initial sanity checks
snv_truth_file = args.SNV_truth_file
sv_truth_file = args.SV_truth_file
input_files = args.input
output_prefix = args.output

# Declaration of bam_path elements
bam_paths = []

# Ensure the BAM file path exists
if input_files.endswith(".bam"):
    path = args.input
    if os.path.exists(path):
        bam_paths.append(path)
    else:
        print(f"File not found: {path}")
        sys.exit()
# Ensure files are sorted by BAM size given a .txt file
elif input_files.endswith(".txt"):
    # Add each line in the txt file to the 'bam_paths' list and grab each BAM files' size
    txtfile = args.input
    if os.path.exists(txtfile):
        with open(input_files, 'r') as bam_files:
            file_sizes = []
            for line in bam_files:
                path = line.strip()
                if path:
                    if path.endswith(".bam"):
                        path = os.path.expandvars(path)
                        if os.path.exists(path):
                            bam_paths.append(path)
                            file_sizes.append(os.path.getsize(path))
                        else:
                            print(f"File not found: {path}")
                            sys.exit()
                    else:
                        print(f"File not a bam file: {path}")
                        sys.exit()
        # Notify the user that the output prefix will do nothing when a truth dataset is provided.
        if output_prefix.endswith(".txt"):
            print(f"Notice: Multi BAM input detected. VCF files will be named using the prefixes listed in {output_prefix}.")
        elif snv_truth_file:
            print(f"Notice: Multi BAM input detected. Provided output prefix will be ignored. VCF files will be generated in the same directory as each BAM file with the same file name.")
        else:
            print(f"Notice: Multi BAM input detected. VCF files will be generated in the same directory as each BAM file with the same file name.")
    else:
        print(f"File not found: {txtfile}")
        sys.exit()

    # Ensure BAM files are in increasing file size order
    wrong = 0
    for i in range(len(file_sizes)-1):
        if file_sizes[i] > file_sizes[i+1]:
            wrong = 1
    if wrong == 1:
        print("Files are not ordered by increasing file size in the provided text file. Note that this may affect time for variant simulation.")
else:
    print("No .bam or .txt file detected. Check extensions.")
    sys.exit()

output_prefixes = None
if output_prefix.endswith(".txt"):
    if not os.path.exists(output_prefix):
        print(f"File not found: {output_prefix}")
        sys.exit()
    with open(output_prefix, 'r') as prefix_file:
        output_prefixes = [os.path.expandvars(line.strip()) for line in prefix_file if line.strip()]
    if len(output_prefixes) != len(bam_paths):
        print(f"Mismatch: {len(output_prefixes)} output prefixes found in {output_prefix}, but {len(bam_paths)} BAM paths provided.")
        sys.exit()


## Remaining variable declarations
ref_path = args.reference
ref_chroms = args.ref_chroms
seed = args.seed if args.seed is not None else 0
sub = args.substitution_rate if args.substitution_rate is not None else 1
bp_shift = 0 if args.bp_shift is None else int(args.bp_shift)
ignore_minimum_cov = args.ignore_minimum_cov if args.ignore_minimum_cov is not None else False
sv_rng = numpy.random.default_rng(seed + 1)

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
insdelsnv = args.insdel_snv_rate if args.insdel_snv_rate is not None else 0.5


## Methods
# Grabs the lengths of a single BAM file or FASTA file
def get_chrom_lengths(path, ref_chroms=None):
    # Get chromosomes and lengths for a single BAM file
    if path.endswith(".bam"):
        with AlignmentFile(path) as bam:
            bam_chroms = dict(zip(bam.references, bam.lengths))

        # Filter if ref_chroms is specified and verify alignment
        if ref_chroms is not None:
            if not set(ref_chroms).issubset(bam_chroms.keys()):
                print(f"Warning: One or more values from \'--ref_chroms\' does not match with the BAM file: {path}\nEnsure you are using the correct reference name.")
                print("\nBAM references: ", bam_chroms.keys())
                print("\nRef_chroms: ", ref_chroms)
                return 1
            bam_chroms = {chrom: length for chrom, length in bam_chroms.items() if chrom in ref_chroms}
        return bam_chroms

    # Get chroms & lengths for a FASTA file
    if path.endswith(".fasta") or path.endswith(".fa"):
        with FastaFile(path) as fasta:
            fasta_chroms = dict(zip(fasta.references, fasta.lengths))

        # Filter if ref_chroms is specified and verify alignment
        if ref_chroms is not None:
            if not set(ref_chroms).issubset(fasta_chroms.keys()):
                print("Warning: One or more values from \'--ref_chroms\' does not match with the FASTA file. Ensure you are using the correct reference name.")
                print("\nFASTA references: ", fasta_chroms.keys())
                print("\nRef_chroms: ", ref_chroms)
                return 1
            fasta_chroms = {chrom: length for chrom, length in fasta_chroms.items() if chrom in ref_chroms}
        return fasta_chroms


    # Get chroms & lengths for a single BAM file
    if path.endswith(".bam"):
        with AlignmentFile(path) as bam:
            bam_chroms = dict(zip(bam.references, bam.lengths))

        # Filter if ref_chroms is specified and verify alignment
        if ref_chroms is not None:
            if not set(ref_chroms).issubset(bam_chroms.keys()):
                print(f"Warning: One or more values from \'--ref_chroms\' does not match with the BAM file: {path}\nEnsure you are using the correct reference name.")
                print("\nBAM references: ", bam_chroms.keys())
                print("\nRef_chroms: ", ref_chroms)
                return 1
            bam_chroms = {chrom: length for chrom, length in bam_chroms.items() if chrom in ref_chroms}
        return bam_chroms

    # Get chroms & lengths for a FASTA file
    if path.endswith(".fasta") or path.endswith(".fa"):
        with FastaFile(path) as fasta:
            fasta_chroms = dict(zip(fasta.references, fasta.lengths))

        # Filter if ref_chroms is specified and verify alignment
        if ref_chroms is not None:
            if not set(ref_chroms).issubset(fasta_chroms.keys()):
                print("Warning: One or more values from \'--ref_chroms\' does not match with the FASTA file. Ensure you are using the correct reference name.")
                print("\nFASTA references: ", fasta_chroms.keys())
                print("\nRef_chroms: ", ref_chroms)
                return 1

            fasta_chroms = {chrom: length for chrom, length in fasta_chroms.items() if chrom in ref_chroms}
        return fasta_chroms

def get_chrom_lengths_for_sv(bam_path, ref_chroms=None):
    with pysam.AlignmentFile(bam_path) as bam:
        chromol = dict(zip(bam.references, bam.lengths))

    if ref_chroms is not None:
        # Filter only if ref_chroms is specified
        chromol = {c: l for c, l in chromol.items() if c in ref_chroms}

    return chromol, tuple(chromol.keys())

# Counts the number of decimal places in the minAFsnv or minAFsv variables
def count_decimals(num):
    s = str(num)
    if '.' in s:
        # Return only significant digits (including leading zeros after the decimal)
        return len(s.split('.')[-1].rstrip('0'))
    else:
        return 0

# Create a list of lists ('records') of all variants with their chromosome, position, reference, alt, and AF value in the truth VCF file
def parse_truth_vcf(vcf_path):
    records = []
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, _, ref, alt, _, _, info, *_ = parts
            # Find AF value, set a default value if not acceptable
            af_val = 0.05
            for field in info.split(';'):
                if field.startswith('AF='):
                    try:
                        af_val = float(field.split('=')[1])
                    except:
                        pass

            records.append((chrom, pos, ref, alt, af_val))
    return records

def parse_truth_sv_vcf(vcf_path):
    records = []
    with open(vcf_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            chrom, pos, _, ref, alt, _, _, info, *_ = parts
            info_dict = {}
            for field in info.split(';'):
                if '=' in field:
                    key, value = field.split('=', 1)
                    info_dict[key] = value
                elif field:
                    info_dict[field] = True
            af_val = 0.05
            try:
                af_val = float(info_dict.get('AF', 0.05))
            except (TypeError, ValueError):
                pass
            records.append((chrom, pos, ref, alt, af_val, info_dict))
    return records

def describe_truth_sv(pos, ref, alt, info):
    symbolic = alt.startswith('<')

    svtype = info.get('SVTYPE')
    if not isinstance(svtype, str):
        if symbolic:
            svtype = alt.strip('<>').split(':')[0]
        else:
            svtype = 'INS' if len(alt) > len(ref) else 'DEL'

    if isinstance(info.get('SVLEN'), str):
        svlen = int(info['SVLEN'].split(',')[0])
    elif not symbolic:
        svlen = len(alt) - len(ref)
    elif isinstance(info.get('END'), str):
        svlen = int(info['END']) - pos
    else:
        svlen = 0

    svlen = -abs(svlen) if svtype == 'DEL' else abs(svlen)
    end = pos + abs(svlen) if svtype == 'DEL' else pos
    return svtype, svlen, end

# Generate a location to generate a SNV based on the size of each chromosome
def genlocSNV(bam_paths, mincov=20):
    # Start the generation time
    start_time = time.time()
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Se the random seed here for reproducibility
    numpy.random.seed(seed)
    ref_fasta = pysam.FastaFile(ref_path)

    # Calculate total genome length based on the BAM file or reference chromosomes provided
    bam_chroms = get_chrom_lengths(bam_paths[0], ref_chroms=ref_chroms)
    genome_length = sum(bam_chroms.values())

    # Calculate weights for distributed variant simulation
    p = []
    for chrom in bam_chroms.keys():
        p.append(bam_chroms[chrom] / genome_length)

    # Create a coverage dictionary to store all coverages in all BAM files
    # key: BAM file
    # value: list of coverages at each variant locus
    callable_depth_dict = {}
    coverage_dict = {}
    for path in bam_paths:
        callable_depth_dict[path] = []
        coverage_dict[path] = []
    ## Generate variants
    num_bams = len(bam_paths)
    chrom_and_loc = []
    taken = set()

    for snv in range(numsnv):
        print(f"Creating SNV {snv+1}/{numsnv}...")
        attempts = 0

        # Exits once an acceptable SNP is found
        while True:
            attempts += 1
            if attempts > MAX_ATTEMPTS:
                print(f"Error: Could not find a valid SNV locus after {MAX_ATTEMPTS} attempts. Check coverage (minimum {mincov}), --ref_chroms, or use -imc.")
                sys.exit()

            # Select a random chromosome and locus based on the chromosome ratio to total genome length
            #rand_chrom = numpy.random.choice(bam_chroms.keys(), p=p)
            rand_chrom = str(numpy.random.choice(list(bam_chroms.keys()), 1, p=p)[0])
            rand_locus_int = numpy.random.randint(0,bam_chroms[rand_chrom])
            if rand_locus_int == 0:
                rand_locus_int = 1
            rand_locus = str(rand_locus_int+bp_shift)
            if int(rand_locus) < 1 or int(rand_locus) > bam_chroms[rand_chrom]:
                continue
            if (rand_chrom, rand_locus) in taken:
                continue

            # Try a new location if the reference nt is 'N' or a gap
            ref_base = ref_fasta.fetch(rand_chrom, int(rand_locus)-1, int(rand_locus)).upper()
            if ref_base not in ['A', 'C', 'G', 'T']:
                continue

            # Check for minimum coverage requirement in each BAM file at 1 specific variant locus
            failed = 0
            num_above_mincov = 0

            # Create a list of numbers that represents only the number of reads that a variant caller is likely to see
            callable_depths = []
            coverages = []
            for path in bam_paths:
                indiv_depth = 0
                try:
                    # Grab the depth at a single location
                    indiv_depth = int(depth(path,'-r',rand_chrom+":"+rand_locus+"-"+rand_locus).rstrip("\n").split("\t")[-1])
                except (ValueError, IndexError, pysam.utils.SamtoolsError):
                    indiv_depth = 0

                # Keep track of the number of BAMs that meet the minimum coverage requirement for each locus
                if ignore_minimum_cov or indiv_depth >= mincov:
                    start_pos = int(rand_locus)
                    start_pos -= 1
                    end_pos = start_pos+1
                    bamfile = AlignmentFile(path, "rb")
                    num_reads_left = indiv_depth
                    num_good_reads = 0

                    # Check that there are an appropriate number of good reads to modify to match the VAF
                    for read in bamfile.fetch(rand_chrom, start_pos, end_pos):
                        if read.query_sequence is not None and read.flag in {0, 16, 99, 147, 83, 163}:
                            # Check that there is not a gap in the read at a given position. Excludes all positions that are insertions, deletions, or gaps in the read
                            # NOTEE: this does not take into consideration any SNVs that are more than 1 bp long
                            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                                if ref_pos == start_pos and query_pos is None:
                                    continue
                                if ref_pos == start_pos and query_pos is not None:
                                    num_good_reads += 1
                        num_reads_left -= 1

                    bamfile.close()
                    # Check that it's possible to modify the reads according to the specified min and max VAF values
                    max_num_reads_to_modify = round(num_good_reads * maxAFsnv)
                    if num_good_reads == 0 or max_num_reads_to_modify < 1:
                        failed = 1
                        break

                    # Ensure that the EAF won't be smaller than 1%
                    temp_EAF = max_num_reads_to_modify / num_good_reads
                    if temp_EAF < 0.01:
                        failed = 1
                        break

                    callable_depths.append(num_good_reads)
                    coverages.append(indiv_depth)
                    num_above_mincov += 1
                # If the current BAM file does not meet the minimum coverage requirment, try another locus
                else:
                    failed = 1
                    break

            # Try another locus if none of the AFs work at that position
            if failed == 1 or num_above_mincov != num_bams:
                continue

            # If all BAM files meet the coverage requirement, move to the next SNP
            if num_above_mincov == len(bam_paths):
                chrom_and_loc.append(tuple([rand_chrom,rand_locus]))
                taken.add((rand_chrom, rand_locus))
                for i, path in enumerate(bam_paths):
                    callable_depth_dict[path].append(callable_depths[i])
                    # RETRIEVE AND STORE COVERAGE HERE:
                    coverage_dict[path].append(coverages[i])
                break

    end_time = time.time()
    print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total run time: {end_time - start_time:.2f} seconds")

    return tuple(chrom_and_loc), callable_depth_dict, coverage_dict

# Generate a location to generate a SV
def genlocSV(num, bam_paths, mincov=20):
    ref_chroms = args.ref_chroms
    # Use the first BAM to get chrom lengths
    chromol, chrom = get_chrom_lengths_for_sv(bam_paths[0], ref_chroms=ref_chroms)
    chrom_names = list(chrom)
    weights = numpy.array([chromol[c] for c in chrom_names], dtype=float)
    weights /= weights.sum()
    ref_fasta = pysam.FastaFile(ref_path)

    locations = []
    for i in range(num):
        print(f"Creating SV {i+1}/{num}...")
        attempts = 0
        while True:
            attempts += 1
            if attempts > MAX_ATTEMPTS:
                print(f"Error: Could not find a valid SV locus after {MAX_ATTEMPTS} attempts. Try fewer SVs, a smaller -maxsvl, a larger -minAFsv (minimum coverage {mincov}), or -imc.")
                sys.exit()

            ranchrom = chrom_names[sv_rng.choice(len(chrom_names), p=weights)]
            loc_int = int(sv_rng.integers(1, chromol[ranchrom] + 1)) + bp_shift

            if loc_int < 1 or loc_int + maxsvl > chromol[ranchrom]:
                continue

            # Check for overlapping/nearby SVs
            too_close = False
            for existing in locations:
                if existing[0] == ranchrom and abs(int(existing[1]) - loc_int) < SV_MIN_GAP:
                    too_close = True
                    break

            if too_close:
                continue

            # Try a new location if the reference nt is 'N' or a gap
            ref_base = ref_fasta.fetch(ranchrom, loc_int - 1, loc_int).upper()
            if ref_base not in ['A', 'C', 'G', 'T']:
                continue

            # --- Get coverage for ALL BAM files ---
            loc = str(loc_int)
            covers = []
            for bam_path in bam_paths:
                try:
                    cover_val = int(depth(bam_path, '-r', f"{ranchrom}:{loc}-{loc}").rstrip("\n").split("\t")[-1])
                except (ValueError, IndexError, pysam.utils.SamtoolsError):
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
        locations.append((ranchrom, loc, covers, ref_base))

    ref_fasta.close()
    return tuple(locations)

# Generates SNVs (the 'ALT' column in the output VCF file)
# snplist format: ((chr1,586528,76,T),...)
def gensnps(ref_nts = []):
    numpy.random.seed(seed)
    result = []
    nucl = ["A","T","C","G"]
    for i in range(len(ref_nts)):
        draw = choice(tuple(['snp','indel']), 1, p=[sub,1-sub])
        if draw == 'snp':
            nucl = ["A","T","C","G"]
            # Remove the nt already in the reference genome
            nucl.remove(ref_nts[i][3][0].upper())
            alt_nt = str(choice(nucl,1)[0])
            result.append(tuple(list(ref_nts[i][0:3])+[ref_nts[i][3][0].upper()]+[alt_nt]))
        else:
            draw = choice(tuple(['ins','del'],), 1, p= [insdelsnv,1-insdelsnv])
            if draw == 'ins':
                result.append(tuple(list(ref_nts[i][0:3])+[ref_nts[i][3][0].upper()]+[ref_nts[i][3][0].upper()+genseq(1,maxsnvl-1).upper()]))
            else:
                rmlen = choice(range(maxsnvl+1)[1:])
                result.append(tuple(list(ref_nts[i][0:3])+[ref_nts[i][3].upper()[:rmlen]]+[ref_nts[i][3].upper()[0]]))

    return tuple(result)

# Grabs the reference sequence with the reference chromosomes
def addrefnt(snplist=1):
    with open(ref_path,'r') as fasta:
        data = []
        chroms = []
        for i in snplist:
            chroms.append(i[0])
        chroms = set(chroms)
        start = 1
        for line in fasta:
            if line.startswith(">"):
                if start==1:
                    start=0
                else:
                    if check==1:
                        data.append(tuple([chromosome,seq]))

                if "chr"==snplist[0][0][0:3]:
                    chromosome=line.lstrip(">").split(" ")[0].rstrip("\n")
                else:
                    chromosome=line.lstrip(">").split(" ")[0].lstrip("chr").rstrip("\n")
                seq=''
                check=0
                if chromosome in chroms:
                    check=1
            else:
                if check == 0:
                    continue
                seq += line.rstrip("\n")
        if check == 1:
            data.append(tuple([chromosome,seq]))
        data=tuple(data)
    fasta.close()

    # Modify the input snplist variable
    snplistmod = []
    for i in snplist:
        for j in data:
            if i[0] == j[0]:
                snplistmod.append(tuple(list(i)+[j[1][int(i[1])-1:int(i[1])+maxsnvl]]))
                break

    # Return an updated snplist with the reference nt (to be a part of the output VCF file)
    return tuple(snplistmod)

# Generates a sequence of random nt's
def genseq(minl,maxl):
    numpy.random.seed(seed)
    rand_length = choice(range(minl,maxl+1))

    # Generate a list of random nucleotides
    nucl = tuple(["A","T","C","G"])
    result = ""
    for i in range(rand_length):
        result += choice(nucl)
    return result

#  Writes passing and failing VCF files
def writevcf(path, vcf_header_snv, variants, failed_variants, num_failed):
    if output_prefixes is not None:
        prefix = output_prefixes[bam_paths.index(path)]
        vcf = f"{prefix}_SNV.vcf"
        new_vcf = os.path.expandvars(vcf)
        if os.path.exists(new_vcf):
            counter = 2
            while os.path.exists(os.path.expandvars(f"{prefix}_SNV_{counter}.vcf")):
                counter += 1
            new_vcf = os.path.expandvars(f"{prefix}_SNV_{counter}.vcf")

    # Write the passing VCF for a single BAM input file and ensure a file with the same name does not get overwritten
    elif len(bam_paths) == 1:
        vcf = f"{output_prefix}_SNV.vcf"
        new_vcf = os.path.expandvars(vcf)
        if os.path.exists(new_vcf):
            counter = 2
            while os.path.exists(os.path.expandvars(f"{output_prefix}_SNV_{counter}.vcf")):
                counter += 1
            new_vcf = os.path.expandvars(f"{output_prefix}_SNV_{counter}.vcf")

    # Write the passing VCF for a multiple BAM input files and ensure a file with the same name does not get overwritten
    # The path to the VCF files is the same path to each individual BAM file
    else:
        vcf = f"{path[:-4]}_SNV.vcf"
        new_vcf = os.path.expandvars(vcf)
        # BUGFIX: this previously checked os.path.exists(path) -- "path" is the
        # INPUT BAM file, which always exists (we just read from it), so this
        # branch was *always* taken and every multi-BAM run unconditionally
        # skipped "<name>_SNV.vcf" in favor of "<name>_SNV_2.vcf" (or higher),
        # regardless of whether a prior output actually existed. The correct
        # check is whether the intended *output* VCF (new_vcf) already exists.
        if os.path.exists(new_vcf):
            counter = 2
            while os.path.exists(os.path.expandvars(f"{path[:-4]}_SNV_{counter}.vcf")):
                counter += 1
            new_vcf = os.path.expandvars(f"{path[:-4]}_SNV_{counter}.vcf")

    os.makedirs(os.path.dirname(new_vcf) or ".", exist_ok=True)
    with open(new_vcf, "w") as f:
        f.write('\n'.join(vcf_header_snv) + "\n")
        for v in variants:
            f.write(v + '\n')
        print(f"New truth VCF file written to {new_vcf}")

    # Write the failing VCF if at least one variant failed to meet specified parameters given a truth dataset and ensure a file with the same name does not get overwritten
    if num_failed > 0:
        failed_prefix = output_prefixes[bam_paths.index(path)] if output_prefixes is not None else output_prefix
        vcf_failed = f"{failed_prefix}_SNV_FAILED.vcf"
        new_vcf_failed = os.path.expandvars(vcf_failed)
        if os.path.exists(new_vcf_failed):
            counter = 2
            while os.path.exists(os.path.expandvars(f"{failed_prefix}_SNV_FAILED_{counter}.vcf")):
                counter += 1
            new_vcf_failed = os.path.expandvars(f"{failed_prefix}_SNV_FAILED_{counter}.vcf")

        os.makedirs(os.path.dirname(new_vcf_failed) or ".", exist_ok=True)
        with open(new_vcf_failed, "w") as f:
            f.write('\n'.join(vcf_header_snv) + '\n')
            for v in failed_variants:
                f.write(v + '\n')
        print(f"Failed VCF written to {new_vcf_failed}")

    if snv_truth_file is not None:
        if num_failed != 0:
            print(f"{num_failed} variant(s) could not be aligned to do an insufficient quantity of variants that meet the provided input requirement(s) in the truth VCF file.")


## Main function
def main():
    # Set global seed values
    numpy.random.seed(seed)
    random.seed(seed)

    # Verify the BAM(s) match the reference FASTA and gather lengths of the BAMs
    fasta_chroms = get_chrom_lengths(ref_path, ref_chroms)
    if isinstance(fasta_chroms, int): sys.exit()

    bam_chroms = {}
    error = 0
    # Ensure that each BAM file is aligned with the reference file
    for path in bam_paths:
        with AlignmentFile(path, "rb") as bam:
            if not bam.has_index():
                print(f"Error: {path} is not indexed. Run: samtools index {path}")
                sys.exit()
        bam_chroms = get_chrom_lengths(path, ref_chroms)
        if isinstance(bam_chroms, int): sys.exit()
        if set(bam_chroms.keys()).issubset(fasta_chroms.keys()):
            for chrom in bam_chroms.keys():
                if bam_chroms[chrom] != fasta_chroms[chrom]:
                    print(f"Warning: {path} is not aligned with the reference file.")
                    error = 1
        else:
            print(f"Warning: {path} contains chromosomes that are not in the reference file.")
            error = 1

    if error == 1:
        sys.exit()

    # Calculate the genome length off the last BAM file if all match to the reference
    genome_length = 0
    for length in bam_chroms.values():
        genome_length += length

    # VCF header information
    vcf_header_snv = [
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=FAIL,Description="Failed minimum coverage">',
        '##FORMAT=<ID=CO,Number=1,Type=Integer,Description="Coverage; Number of reads at a specific locus">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of reads that pass quality filter; callable depth">',
        '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Target Allele Frequency">',
        '##INFO=<ID=EAF,Number=A,Type=Float,Description="Exact Allele Frequency (variant reads / callable depth)">'
    ]

    # If a truth file is provided, simulate aligned variants using the BAM file and generate 2 VCF files of variants that successfully aligned and failed to align
    # SNV processing
    if snv_truth_file:
        snvs_in_truth = parse_truth_vcf(snv_truth_file)

        # Check that there are no more variants being requested for alignment than are in the truth dataset
        numsnv_align = args.number_of_snvs_to_align if args.number_of_snvs_to_align is not None else len(snvs_in_truth)
        if numsnv_align > len(snvs_in_truth):
            print(f"Warning: The number of variants you requested to align ({numsnv_align}) is greater than the number of variants in the provided truth VCF file ({len(snvs_in_truth)}).")
            sys.exit()

        if len(bam_paths) == 1:
            # Extend the header
            vcf_header_snv.extend([f'##contig=<ID={chrom},length={length}>' for chrom, length in bam_chroms.items()])
            vcf_header_snv.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))

            # Align BAM with VCF (maybe change loop description)
            num_failed = 0
            variants = []
            failed_variants = []
            for i in range(numsnv_align):
                print(f"\tRetrieving SNV {i+1}/{numsnv_align} from truth vcf...")
                chrom, pos, ref, alt, af_val = snvs_in_truth[i]

                # Update the position if a bp shift is specified
                pos = str(int(pos)+bp_shift)

                try:
                    cover = int(depth(bam_paths[0], '-r', chrom + ":" + pos + "-" + pos).rstrip("\n").split("\t")[-1])
                except (ValueError, IndexError, pysam.utils.SamtoolsError):
                    cover = 0

                # Change data type for position back to integer
                pos = int(pos)

                # Ensure the user wants to use this variant
                if ((args.minimum_allele_frequency_snv is None) or af_val >= minAFsnv) and ((args.maximum_allele_frequency_snv is None) or af_val <= maxAFsnv) and ((ref_chroms is None) or chrom in ref_chroms):
                    # Find the callable depth at that locus
                    bamfile = AlignmentFile(bam_paths[0], "rb")
                    start_pos = int(pos)
                    start_pos -= 1
                    end_pos = start_pos+1
                    num_good_reads = 0

                    for read in bamfile.fetch(chrom, start_pos, end_pos):
                        if read.query_sequence is not None and read.flag in {0, 16, 99, 147, 83, 163}:
                            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                                if ref_pos == start_pos and query_pos is None:
                                    continue
                                if ref_pos == start_pos and query_pos is not None:
                                    num_good_reads += 1
                    bamfile.close()

                else:
                    failed_variants.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1500\tFAIL\tAF={af_val};EAF=0.0\tCO:AD:DV\t{cover}:.:0")
                    continue

                # Check that it's possible to modify the reads according to the specified min and max VAF values
                num_read_to_modify = round(af_val * num_good_reads)
                if num_read_to_modify > 0:
                    exact_af = round(num_read_to_modify / num_good_reads, 4)
                    variants.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1500\tPASS\tAF={af_val};EAF={exact_af}\tCO:AD:DV\t{cover}:{num_good_reads}:{num_read_to_modify}")
                else:
                    num_failed += 1
                    exact_af = round(num_read_to_modify / num_good_reads, 4) if num_good_reads > 0 else 0.0
                    failed_variants.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1500\tFAIL\tAF={af_val};EAF={exact_af}\tCO:AD:DV\t{cover}:{num_good_reads}:0")

            # Write the final VCF files for passing and failing variants
            writevcf(bam_paths[0], vcf_header_snv, variants, failed_variants, num_failed)

        else:
            vcf_header_snvs = []
            for path in range(len(bam_paths)):
                # Extend the header
                vcf_header_snv_copy = vcf_header_snv.copy()
                vcf_header_snv_copy.extend([f'##contig=<ID={chrom},length={length}>' for chrom, length in bam_chroms.items()])
                vcf_header_snv_copy.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))
                vcf_header_snvs.append(vcf_header_snv_copy)

            # Align BAMs with VCF
            for b, path in enumerate(bam_paths):
                num_failed = 0
                variants_processed = 0
                variants = []
                failed_variants = []
                for i in range(numsnv_align):
                    print(f"\tRetrieving SNV {i+1}/{numsnv_align} from truth vcf...")
                    chrom, pos, ref, alt, af_val = snvs_in_truth[i]

                    # Update the position if a bp shift is specified
                    pos = str(int(pos)+bp_shift)

                    # Calculate the coverage at a single position
                    cover = 0
                    try:
                        cover = int(pysam.depth(bam_paths[b], '-r', chrom + ":" + pos + "-" + pos).rstrip("\n").split("\t")[-1])
                    except (ValueError, IndexError, pysam.utils.SamtoolsError):
                        cover = 0

                    # Change data type for position back to integer
                    pos = int(pos)

                     # Ensure the user wants to use this variant
                    if ((args.minimum_allele_frequency_snv is None) or af_val >= minAFsnv) and ((args.maximum_allele_frequency_snv is None) or af_val <= maxAFsnv) and ((ref_chroms is None) or chrom in ref_chroms):
                        # Find the callable depth at that locus
                        bamfile = AlignmentFile(bam_paths[b], "rb")
                        start_pos = int(pos)
                        start_pos -= 1
                        end_pos = start_pos+1
                        num_good_reads = 0

                        for read in bamfile.fetch(chrom, start_pos, end_pos):
                            if read.query_sequence is not None and read.flag in {0, 16, 99, 147, 83, 163}:
                                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                                    if ref_pos == start_pos and query_pos is None:
                                        continue
                                    if ref_pos == start_pos and query_pos is not None:
                                        num_good_reads += 1
                        bamfile.close()

                    else:
                        failed_variants.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1500\tFAIL\tAF={af_val};EAF=0.0\tCO:AD:DV\t{cover}:.:0")
                        continue

                    # Check that it's possible to modify the reads according to the specified min and max VAF values
                    num_read_to_modify = round(af_val * num_good_reads)
                    if num_read_to_modify > 0:
                        exact_af = round(num_read_to_modify / num_good_reads, 4)
                        variants.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1500\tPASS\tAF={af_val};EAF={exact_af}\tCO:AD:DV\t{cover}:{num_good_reads}:{num_read_to_modify}")
                    else:
                        num_failed += 1
                        exact_af = round(num_read_to_modify / num_good_reads, 4) if num_good_reads > 0 else 0.0
                        failed_variants.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1500\tFAIL\tAF={af_val};EAF={exact_af}\tCO:AD:DV\t{cover}:{num_good_reads}:0")

                # Write the final VCF files for passing and failing variants
                writevcf(bam_paths[b], vcf_header_snvs[b], variants, failed_variants, num_failed)
    else:
        if numsnv > 0:
            vcf_header_snvs = []
            tot_bam_len = str(len(bam_paths))

            # SNV location is generated based on coverage of all BAM files
            """
            Format:
            chrom_and_loc: ((chr1, 326542683), (chr1, 75495683),...)
            depth_dict: {path1: [54,98,65,...],
                            path2: [76,45,97,...],...}
            """
            chrom_and_loc, callable_depth_dict, covers_dict = genlocSNV(bam_paths, ceil(1/maxAFsnv))

            # Rewrite snvloc how it was originally written
            # snvloc = ((chr1, 68745375, 98), (chr1, 6883452, 87),...)
            snvloc = []
            for j, depth_ in enumerate(callable_depth_dict[path]):
                snvloc.append([chrom_and_loc[j][0], chrom_and_loc[j][1], depth_])

            ref_nts = addrefnt(snvloc)
            elements = gensnps(ref_nts)

            # Set the allele frequencies that are compatable with all BAMs
            AFs = []
            random.seed(seed)

            # itereate through each variant that was generated ('v' acts as an index)
            # iterate the-number-of-variants-you-need-to-generate times
            mindepth = round(1/maxAFsnv)
            for variant in range(numsnv):
                af_attempts = 0
                while True:
                    af_attempts += 1
                    if af_attempts > MAX_ATTEMPTS:
                        print(f"Error: Could not find an allele frequency between {minAFsnv} and {maxAFsnv} that gives at least 1 variant read for SNV {variant+1} in every BAM file.")
                        sys.exit()
                    failed = 0
                    differences = []
                    EAFs = []
                    AFnum = round(random.uniform(minAFsnv, maxAFsnv), count_decimals(minAFsnv))

                    for key in callable_depth_dict.keys():
                        readnum = round(AFnum * callable_depth_dict[key][variant])
                        if readnum < 1:
                            failed = 1
                            break

                        if callable_depth_dict[key][variant] > 0:
                            EAF = round(readnum / callable_depth_dict[key][variant], 4)
                        else:
                            EAF = 0.0
                        EAFs.append(EAF)

                    if failed == 1:
                        continue
                    if failed == 0:
                        break

                AFs.append([AFnum, tuple(EAFs)])

            for i, path in enumerate(bam_paths):
                print(f"Writing SNV output file {i+1}/{tot_bam_len}:")
                # Extend the header
                vcf_header_snv_copy = vcf_header_snv.copy()
                vcf_header_snv_copy.extend([f'##contig=<ID={chrom},length={length}>' for chrom, length in bam_chroms.items()])
                vcf_header_snv_copy.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))
                vcf_header_snvs.append(vcf_header_snv_copy)
                variants = []
                for d, (k, p, c) in enumerate(zip(elements, AFs, callable_depth_dict[path])):
                    reads_to_mod = round(p[0] * callable_depth_dict[path][d])
                    variants.append(f"{k[0]}\t{k[1]}\t.\t{k[3]}\t{k[4]}\t1500\tPASS\tAF={p[0]};EAF={p[1][i]}\tCO:AD:DV\t{covers_dict[path][d]}:{callable_depth_dict[path][d]}:{reads_to_mod}")

                failed_variants = 0
                variants_processed = 0
                num_failed = 0
                writevcf(path, vcf_header_snvs[i], variants, failed_variants, num_failed)

    # SV processing (in beta)
    if sv_truth_file is not None:
        print(f"Retrieving SVs from truth vcf...")
        svloc = parse_truth_sv_vcf(sv_truth_file)

        # Check that there are no more variants being requested for alignment than are in the truth dataset
        numsv_align = args.number_of_svs_to_align if args.number_of_svs_to_align is not None else len(svloc)
        if numsv_align > len(svloc):
            print(f"Warning: The number of SVs you requested to align ({numsv_align}) is greater than the number of SVs in the provided truth VCF file ({len(svloc)}).")
            sys.exit()
        svloc = svloc[:numsv_align]

        # Loop over every BAM file for SV truth parsing
        for idx, bam_path in enumerate(bam_paths):
            print(f"\n--- Processing BAM {idx+1}/{len(bam_paths)} for SVs: {bam_path} ---")

            # 1. Determine the correct output VCF path
            if output_prefixes is not None:
                current_SVvcf = f"{output_prefixes[idx]}_SV.vcf"
            elif len(bam_paths) > 1:
                if output_prefix.endswith('.vcf'):
                    current_SVvcf = f"{output_prefix[:-4]}_SV_{idx+1}.vcf"
                else:
                    current_SVvcf = f"{output_prefix}_SV_{idx+1}.vcf"
            else:
                current_SVvcf = f"{output_prefix}_SV.vcf"

            chromol, chrom = get_chrom_lengths_for_sv(bam_path)

            vcfsv = [
                '##fileformat=VCFv4.2',
                '##ALT=<ID=INS,Description="Insertion">',
                '##ALT=<ID=DEL,Description="Deletion">',
                '##FILTER=<ID=PASS,Description="All filters passed">',
                '##FILTER=<ID=FAIL,Description="Failed minimum coverage">',
                '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">',
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">',
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">',
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                '##FORMAT=<ID=CO,Number=1,Type=Integer,Description="Coverage; Number of reads at a specific locus">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of reads that pass quality filter; callable depth">',
                '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">',
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Target Allele Frequency">',
                '##INFO=<ID=EAF,Number=A,Type=Float,Description="Exact Allele Frequency (variant reads / callable depth)">'
            ]
            vcfsv.extend([f'##contig=<ID={c},length={l}>' for c, l in chromol.items()])
            vcfsv.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))

            vcfsv_failed = vcfsv.copy()
            header_len = len(vcfsv)

            for rec_idx, record in enumerate(svloc):
                chrom, pos, ref, alt, af, info = record
                pos = str(int(pos) + bp_shift)

                ## Failsafe coverage calculation tied to specific BAM
                try:
                    cover = int(depth(bam_path, '-r', f"{chrom}:{pos}-{pos}").rstrip("\n").split("\t")[-1])
                except (ValueError, IndexError, pysam.utils.SamtoolsError):
                    cover = 0

                # Find the callable depth at that locus
                bamfile = AlignmentFile(bam_path, "rb")
                start_pos = int(pos)
                start_pos -= 1
                end_pos = start_pos+1
                num_good_reads = 0

                for read in bamfile.fetch(chrom, start_pos, end_pos):
                    if read.query_sequence is not None and read.flag in {0, 16, 99, 147, 83, 163}:
                        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                            if ref_pos == start_pos and query_pos is None:
                                continue
                            if ref_pos == start_pos and query_pos is not None:
                                num_good_reads += 1
                bamfile.close()

                svtype, svlen, end = describe_truth_sv(int(pos), ref, alt, info)
                ID = f"SV{rec_idx+1}"

                ## Test for minimum coverage requirement
                if not ignore_minimum_cov and cover < ceil(1/minAFsv):
                    vcfsv_failed.append(f"{chrom}\t{pos}\t{ID}\t{ref}\t{alt}\t60\tFAIL\tPRECISE;SVTYPE={svtype};SVLEN={svlen};END={end};AF={af}\tDP\t{cover}")
                else:
                    vcfsv.append(f"{chrom}\t{pos}\t{ID}\t{ref}\t{alt}\t60\tPASS\tPRECISE;SVTYPE={svtype};SVLEN={svlen};END={end};AF={af}\tDP\t{cover}")

            # --- File Writing Logic ---
            os.makedirs(os.path.dirname(current_SVvcf) or '.', exist_ok=True)
            with open(current_SVvcf, "w") as f:
                f.write('\n'.join(vcfsv) + '\n')

            print(f"New SV truth VCF file written to {current_SVvcf}.")

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
                ref_base = i[3]
                af = round(float(sv_rng.uniform(minAFsv, maxAFsv)), count_decimals(minAFsv))

                if sv_rng.random() < insdel:
                    seq_len = int(sv_rng.integers(minsvl, maxsvl + 1))
                    seq = ''.join(sv_rng.choice(["A","T","C","G"], size=seq_len))
                    simulated_svs.append({
                        'chrom': i[0], 'pos': i[1], 'covers': i[2],
                        'id': f"HackIns{insertnum}", 'ref': ref_base, 'alt': ref_base + seq,
                        'svtype': 'INS', 'svlen': len(seq), 'end': int(i[1]), 'af': af
                    })
                    insertnum += 1
                else:
                    dellen = int(sv_rng.integers(minsvl, maxsvl + 1))
                    simulated_svs.append({
                        'chrom': i[0], 'pos': i[1], 'covers': i[2],
                        'id': f"HackDel{delnum}", 'ref': ref_base, 'alt': '<DEL>',
                        'svtype': 'DEL', 'svlen': -dellen, 'end': int(i[1]) + dellen, 'af': af
                    })
                    delnum += 1

            contig_order = {c: n for n, c in enumerate(get_chrom_lengths_for_sv(bam_paths[0])[0])}
            simulated_svs.sort(key=lambda var: (contig_order[var['chrom']], int(var['pos'])))
            insertnum = 1
            delnum = 1
            for var in simulated_svs:
                if var['svtype'] == 'INS':
                    var['id'] = f"HackIns{insertnum}"
                    insertnum += 1
                else:
                    var['id'] = f"HackDel{delnum}"
                    delnum += 1

            # 2. Generate a VCF for each BAM file
            for idx, bam_path in enumerate(bam_paths):
                # 1. Determine the correct output VCF path
                if output_prefixes is not None:
                    current_SVvcf = f"{output_prefixes[idx]}_SV.vcf"
                elif len(bam_paths) > 1:
                    if output_prefix.endswith('.vcf'):
                        current_SVvcf = f"{output_prefix[:-4]}_SV_{idx+1}.vcf"
                    else:
                        current_SVvcf = f"{output_prefix}_SV_{idx+1}.vcf"
                else:
                    current_SVvcf = f"{output_prefix}_SV.vcf"

                chromol, chrom = get_chrom_lengths_for_sv(bam_path)

                vcfsv = [
                    '##fileformat=VCFv4.2',
                    '##ALT=<ID=INS,Description="Insertion">',
                    '##ALT=<ID=DEL,Description="Deletion">',
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
                        f"DP\t{cover}"
                    )

                # Write out final files
                # Note: changed output_prefix to avoid errors if missing
                os.makedirs(os.path.dirname(current_SVvcf) or '.', exist_ok=True)

                with open(current_SVvcf, "w") as f:
                    f.write('\n'.join(vcfsv) + '\n')

                print(f"New SV truth vcf for BAM {idx+1} written to {current_SVvcf}.")

if __name__ == "__main__":
    main()
