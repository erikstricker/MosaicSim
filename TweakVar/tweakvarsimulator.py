import sys
import os
import argparse
from math import ceil

import numpy as np
import pysam

MAX_ATTEMPTS = 100000   # cap on location draws per variant before giving up
SV_MIN_GAP = 50000      # minimum spacing between simulated SVs on the same contig

# ============================== Arguments ==============================
parser = argparse.ArgumentParser(description="TweakVarSimulator: A tool for simulating variants")

# Required arguments
parser.add_argument("-i", "--input", required=True,
                    help="Path to an input BAM file, or a .txt file listing BAM paths (one per line, "
                         "ordered from lowest to highest coverage / file size)")
parser.add_argument("-T", "--reference", required=True, help="Path to reference genome (FASTA)")
parser.add_argument("-o", "--output", required=True,
                    help="Output prefix, or a .txt file with one output prefix per BAM (same order as the BAM list)")

# Optional arguments
parser.add_argument("-s", "--seed", type=int, help="Seed for random number generation")
parser.add_argument("-minAFsv", "--minimum_allele_frequency_sv", type=float)
parser.add_argument("-maxAFsv", "--maximum_allele_frequency_sv", type=float)
parser.add_argument("-minAFsnv", "--minimum_allele_frequency_snv", type=float)
parser.add_argument("-maxAFsnv", "--maximum_allele_frequency_snv", type=float)
parser.add_argument("-numsv", "--number_of_svs", type=int)
parser.add_argument("-numsnv", "--number_of_snvs", type=int)
parser.add_argument("-maxsnvl", "--maximum_snv_length", type=int)
parser.add_argument("-minsvl", "--minimum_sv_length", type=int)
parser.add_argument("-maxsvl", "--maximum_sv_length", type=int)
parser.add_argument("-sub", "--substitution_rate", type=float)
parser.add_argument("-insdelsnv", "--insdel_snv_rate", type=float)
parser.add_argument("-insdel", "--insdel_sv_rate", type=float)
parser.add_argument("--ref_chroms", nargs='+', help="Optional list of reference chromosomes to simulate on")
parser.add_argument("-snv", "--SNV_truth_file", help="Optional SNV truth VCF file")
parser.add_argument("-sv", "--SV_truth_file", help="Optional SV truth VCF file")
parser.add_argument("-bps", "--bp_shift", type=int, default=0,
                    help="Shift all bp positions up or downstream by the indicated number of bps")
parser.add_argument("-imc", "--ignore_minimum_cov", action="store_true",
                    help="If set, loci will not be evaluated for minimum coverage")

args = parser.parse_args()


def die(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(1)


def read_list_file(path):
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


# ============================== Input BAMs ==============================
if args.input.endswith(".bam"):
    bam_paths = [args.input]
elif args.input.endswith(".txt"):
    bam_paths = read_list_file(args.input)
else:
    die("Input must be a .bam file or a .txt file listing BAM paths.")

if not bam_paths:
    die(f"No BAM paths found in {args.input}.")

missing = [p for p in bam_paths if not os.path.exists(p)]
if missing:
    die("BAM file(s) not found: " + ", ".join(missing))

file_sizes = [os.path.getsize(p) for p in bam_paths]
if any(file_sizes[i] > file_sizes[i + 1] for i in range(len(file_sizes) - 1)):
    print("Warning: BAMs are not ordered by increasing file size. "
          "This may affect time for variant simulation.")

# ============================== Reference check ==============================
ref_path = args.reference
with pysam.FastaFile(ref_path) as _fa:
    fasta_contigs = dict(zip(_fa.references, _fa.lengths))
contig_order = {c: i for i, c in enumerate(fasta_contigs)}

for bam_path in bam_paths:
    print(f"Checking alignment of {bam_path} with reference...")
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        bam_contigs = dict(zip(bam.references, bam.lengths))
        indexed = bam.has_index()
    if bam_contigs != fasta_contigs:
        die(f"{bam_path} header does not match the reference sequence.")
    if not indexed:
        die(f"{bam_path} is not indexed. Run: samtools index {bam_path}")

print("Proceeding...")

if args.ref_chroms:
    unknown = [c for c in args.ref_chroms if c not in fasta_contigs]
    if unknown:
        die("Contig(s) in --ref_chroms not found in reference: " + ", ".join(unknown))
    sample_contigs = {c: l for c, l in fasta_contigs.items() if c in args.ref_chroms}
else:
    sample_contigs = dict(fasta_contigs)

# ============================== Output paths ==============================
if args.output.endswith(".txt"):
    prefixes = read_list_file(args.output)
    if len(prefixes) != len(bam_paths):
        die(f"Mismatch: {len(prefixes)} output prefixes found, but {len(bam_paths)} BAM paths provided.")
elif len(bam_paths) == 1:
    prefixes = [args.output]
else:
    prefixes = [f"{args.output}_{i + 1}" for i in range(len(bam_paths))]

SNVvcf_paths = [f"{p}_SNV.vcf" for p in prefixes]
SVvcf_paths = [f"{p}_SV.vcf" for p in prefixes]
for p in SNVvcf_paths + SVvcf_paths:
    os.makedirs(os.path.dirname(p) or ".", exist_ok=True)

# ============================== Parameters ==============================
snv_truth_file = args.SNV_truth_file
sv_truth_file = args.SV_truth_file
bp_shift = args.bp_shift
ignore_minimum_cov = args.ignore_minimum_cov
seed = args.seed if args.seed is not None else 0


def opt(value, default):
    return value if value is not None else default


minAFsv = opt(args.minimum_allele_frequency_sv, 0.01)
maxAFsv = opt(args.maximum_allele_frequency_sv, 0.05)
numsv = opt(args.number_of_svs, 50)
minsvl = opt(args.minimum_sv_length, 50)
maxsvl = opt(args.maximum_sv_length, 10000)
insdel = opt(args.insdel_sv_rate, 0.7)

minAFsnv = opt(args.minimum_allele_frequency_snv, 0.01)
maxAFsnv = opt(args.maximum_allele_frequency_snv, 0.05)
numsnv = opt(args.number_of_snvs, 200)
maxsnvl = opt(args.maximum_snv_length, 100)
sub = opt(args.substitution_rate, 1.0)
insdelsnv = opt(args.insdel_snv_rate, 0.5)

if not 0 < minAFsnv <= maxAFsnv <= 1:
    die("SNV allele frequencies must satisfy 0 < minAFsnv <= maxAFsnv <= 1.")
if not 0 < minAFsv <= maxAFsv <= 1:
    die("SV allele frequencies must satisfy 0 < minAFsv <= maxAFsv <= 1.")
if not 1 <= minsvl <= maxsvl:
    die("SV lengths must satisfy 1 <= minsvl <= maxsvl.")
if maxsnvl < 1:
    die("maxsnvl must be at least 1.")
for _name, _val in (("substitution_rate", sub), ("insdel_snv_rate", insdelsnv), ("insdel_sv_rate", insdel)):
    if not 0 <= _val <= 1:
        die(f"{_name} must be between 0 and 1.")

snv_mincov = ceil(1 / maxAFsnv)
sv_mincov = ceil(1 / minAFsv)

# Independent RNG streams so SNV and SV draws don't land on identical loci
snv_rng = np.random.default_rng(seed)
sv_rng = np.random.default_rng(seed + 1)

# ============================== Helpers ==============================
VCF_COLUMNS = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])

SNV_HEADER = [
    '##fileformat=VCFv4.2',
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FILTER=<ID=FAIL,Description="Failed minimum coverage">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
    '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Target Allele Frequency">',
    '##INFO=<ID=EAF,Number=A,Type=Float,Description="Exact Allele Frequency (variant reads / total depth)">',
]

SV_HEADER = [
    '##fileformat=VCFv4.2',
    '##ALT=<ID=INS,Description="Insertion">',
    '##ALT=<ID=DEL,Description="Deletion">',
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FILTER=<ID=FAIL,Description="Failed minimum coverage">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
    '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
]


def write_vcf(path, header, records):
    contigs = [f'##contig=<ID={c},length={l}>' for c, l in fasta_contigs.items()]
    with open(path, "w") as f:
        f.write('\n'.join(header + contigs + [VCF_COLUMNS] + records) + '\n')


def failed_path(path):
    return f"{path[:-4]}_failed.vcf" if path.endswith(".vcf") else f"{path}_failed.vcf"


def count_decimals(x):
    s = str(x)
    return len(s.split('.')[-1].rstrip('0')) if '.' in s else 0


def draw_af(lo, hi, rng):
    decimals = max(count_decimals(lo), count_decimals(hi))
    return round(float(rng.uniform(lo, hi)), decimals)


def random_seq(minl, maxl, rng):
    n = int(rng.integers(minl, maxl + 1))
    return ''.join(rng.choice(list("ACGT"), size=n))


def get_depth(bam_path, chrom, pos):
    try:
        out = pysam.depth(bam_path, "-r", f"{chrom}:{pos}-{pos}")
        return int(out.rstrip("\n").split("\t")[-1])
    except (ValueError, IndexError, pysam.utils.SamtoolsError):
        return 0  # samtools depth prints nothing for zero-coverage sites


def get_depths(chrom, pos):
    return [get_depth(b, chrom, pos) for b in bam_paths]


def passes_coverage(covers, mincov):
    return ignore_minimum_cov or all(c >= mincov for c in covers)


_sample_names = list(sample_contigs)
_sample_weights = np.array([sample_contigs[c] for c in _sample_names], dtype=float)
_sample_weights /= _sample_weights.sum()


def draw_position(rng):
    """Pick a contig weighted by length, then a 1-based position on it (with bp_shift applied)."""
    chrom = _sample_names[rng.choice(len(_sample_names), p=_sample_weights)]
    pos = int(rng.integers(1, sample_contigs[chrom] + 1)) + bp_shift
    return chrom, pos


def parse_truth_vcf(vcf_path):
    records = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, _, _, info = parts[:8]
            info_d = {}
            for field in info.split(';'):
                if '=' in field:
                    k, v = field.split('=', 1)
                    info_d[k] = v
                elif field:
                    info_d[field] = True
            try:
                af = float(info_d.get('AF', 0.05))
            except (TypeError, ValueError):
                af = 0.05
            records.append({'chrom': chrom, 'pos': int(pos), 'id': vid,
                            'ref': ref, 'alt': alt, 'af': af, 'info': info_d})
    return records


# ============================== SNV logic ==============================
def snv_record(chrom, pos, ref, alt, af, cover, passed=True):
    readnum = min(ceil(af * cover), cover) if passed else 0
    eaf = round(readnum / cover, 4) if cover > 0 else 0.0
    filt = "PASS" if passed else "FAIL"
    return (f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t1500\t{filt}\tAF={af};EAF={eaf}"
            f"\tGT:AD:DV\t0/0:{cover - readnum}:{readnum}")


def gen_snv_locations(num, mincov, rng, fasta):
    locations, taken = [], set()
    for i in range(num):
        print(f"Creating SNV {i + 1}/{num}...")
        for _ in range(MAX_ATTEMPTS):
            chrom, pos = draw_position(rng)
            length = fasta_contigs[chrom]
            if pos < 1 or pos > length or (chrom, pos) in taken:
                continue
            # Window = ref base plus up to maxsnvl downstream bases (used for deletions)
            window = fasta.fetch(chrom, pos - 1, min(pos + maxsnvl, length)).upper()
            if window[0] not in "ACGT":
                continue  # N or gap
            covers = get_depths(chrom, pos)
            if passes_coverage(covers, mincov):
                break
        else:
            die(f"Could not find a valid SNV locus after {MAX_ATTEMPTS} attempts. "
                f"Check coverage (min {mincov}) or use -imc.")
        taken.add((chrom, pos))
        locations.append((chrom, pos, window, covers))
    return locations


def make_snv(window, rng):
    ref_base = window[0]
    if rng.random() < sub:
        alt = rng.choice([n for n in "ACGT" if n != ref_base])
        return ref_base, str(alt)
    if rng.random() < insdelsnv or len(window) < 2:
        return ref_base, ref_base + random_seq(1, max(1, maxsnvl - 1), rng)
    dellen = int(rng.integers(1, len(window)))  # always deletes >= 1 base
    return window[:dellen + 1], ref_base


def process_snvs(fasta):
    if snv_truth_file is not None:
        records = parse_truth_vcf(snv_truth_file)
        for idx, bam_path in enumerate(bam_paths):
            print(f"\n--- Processing SNVs for BAM {idx + 1}/{len(bam_paths)}: {bam_path} ---")
            passed, failed = [], []
            for n, rec in enumerate(records, 1):
                if n % 100 == 0 or n == len(records):
                    print(f"Retrieved {n}/{len(records)} SNVs from truth vcf...")
                pos = rec['pos'] + bp_shift
                cover = get_depth(bam_path, rec['chrom'], pos)
                ok = ignore_minimum_cov or cover >= snv_mincov
                (passed if ok else failed).append(
                    snv_record(rec['chrom'], pos, rec['ref'], rec['alt'], rec['af'], cover, passed=ok))

            write_vcf(SNVvcf_paths[idx], SNV_HEADER, passed)
            print(f"New SNV truth vcf written to {SNVvcf_paths[idx]}.")
            if failed:
                write_vcf(failed_path(SNVvcf_paths[idx]), SNV_HEADER, failed)
                print(f"Found {len(failed)} SNVs below minimum coverage. "
                      f"Written to {failed_path(SNVvcf_paths[idx])}.")
        return

    if numsnv <= 0:
        return

    locations = gen_snv_locations(numsnv, snv_mincov, snv_rng, fasta)
    print("Writing the SNV output files...")

    # Materialize variants once so they're identical across all BAM VCFs
    variants = []
    for chrom, pos, window, covers in locations:
        ref, alt = make_snv(window, snv_rng)
        af = draw_af(minAFsnv, maxAFsnv, snv_rng)
        variants.append((chrom, pos, ref, alt, af, covers))
    variants.sort(key=lambda v: (contig_order[v[0]], v[1]))

    for idx in range(len(bam_paths)):
        records = [snv_record(c, p, r, a, af, covers[idx]) for c, p, r, a, af, covers in variants]
        write_vcf(SNVvcf_paths[idx], SNV_HEADER, records)
        print(f"New SNV truth vcf for BAM {idx + 1} written to {SNVvcf_paths[idx]}.")


# ============================== SV logic ==============================
def sv_record(chrom, pos, vid, ref, alt, svtype, svlen, end, af, cover, passed=True):
    filt = "PASS" if passed else "FAIL"
    return (f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t60\t{filt}\t"
            f"PRECISE;SVTYPE={svtype};SVLEN={svlen};END={end};AF={af}\tGT:GQ:DP\t0/0:60:{cover}")


def describe_truth_sv(rec):
    """Work out SVTYPE/SVLEN/END, handling symbolic ALTs like <DEL> via INFO fields."""
    info, ref, alt = rec['info'], rec['ref'], rec['alt']
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
        svlen = int(info['END']) - rec['pos']
    else:
        svlen = 0

    if svtype == 'DEL':
        svlen = -abs(svlen)
    else:
        svlen = abs(svlen)

    pos = rec['pos'] + bp_shift
    end = pos + abs(svlen) if svtype == 'DEL' else pos
    return pos, svtype, svlen, end


def gen_sv_locations(num, mincov, rng):
    locations = []
    for i in range(num):
        print(f"Creating SV {i + 1}/{num}...")
        for _ in range(MAX_ATTEMPTS):
            chrom, pos = draw_position(rng)
            if pos < 1 or pos + maxsvl > fasta_contigs[chrom]:
                continue  # keep the whole SV on the contig
            if any(c == chrom and abs(p - pos) < SV_MIN_GAP for c, p, _ in locations):
                continue
            covers = get_depths(chrom, pos)
            if passes_coverage(covers, mincov):
                break
        else:
            die(f"Could not find a valid SV locus after {MAX_ATTEMPTS} attempts. "
                f"Try fewer SVs, a smaller -maxsvl, or -imc.")
        locations.append((chrom, pos, covers))
    return locations


def process_svs(fasta):
    if sv_truth_file is not None:
        print("Retrieving SVs from truth vcf...")
        records = parse_truth_vcf(sv_truth_file)
        for idx, bam_path in enumerate(bam_paths):
            print(f"\n--- Processing SVs for BAM {idx + 1}/{len(bam_paths)}: {bam_path} ---")
            passed, failed = [], []
            for n, rec in enumerate(records, 1):
                pos, svtype, svlen, end = describe_truth_sv(rec)
                cover = get_depth(bam_path, rec['chrom'], pos)
                vid = rec['id'] if rec['id'] != '.' else f"SV{n}"
                ok = ignore_minimum_cov or cover >= sv_mincov
                (passed if ok else failed).append(
                    sv_record(rec['chrom'], pos, vid, rec['ref'], rec['alt'],
                              svtype, svlen, end, rec['af'], cover, passed=ok))

            write_vcf(SVvcf_paths[idx], SV_HEADER, passed)
            print(f"New SV truth vcf written to {SVvcf_paths[idx]}.")
            if failed:
                write_vcf(failed_path(SVvcf_paths[idx]), SV_HEADER, failed)
                print(f"Found {len(failed)} SVs below minimum coverage. "
                      f"Written to {failed_path(SVvcf_paths[idx])}.")
        return

    if numsv <= 0:
        return

    locations = gen_sv_locations(numsv, sv_mincov, sv_rng)
    print("Writing the SV output files...")

    variants = []
    ins_n = del_n = 1
    for chrom, pos, covers in locations:
        ref_base = fasta.fetch(chrom, pos - 1, pos).upper()
        af = draw_af(minAFsv, maxAFsv, sv_rng)
        if sv_rng.random() < insdel:
            seq = random_seq(minsvl, maxsvl, sv_rng)
            variants.append((chrom, pos, f"HackIns{ins_n}", ref_base, ref_base + seq,
                             "INS", len(seq), pos, af, covers))
            ins_n += 1
        else:
            dellen = int(sv_rng.integers(minsvl, maxsvl + 1))
            variants.append((chrom, pos, f"HackDel{del_n}", ref_base, "<DEL>",
                             "DEL", -dellen, pos + dellen, af, covers))
            del_n += 1
    variants.sort(key=lambda v: (contig_order[v[0]], v[1]))

    for idx in range(len(bam_paths)):
        records = [sv_record(c, p, vid, r, a, t, l, e, af, covers[idx])
                   for c, p, vid, r, a, t, l, e, af, covers in variants]
        write_vcf(SVvcf_paths[idx], SV_HEADER, records)
        print(f"New SV truth vcf for BAM {idx + 1} written to {SVvcf_paths[idx]}.")


# ============================== Main ==============================
def main():
    with pysam.FastaFile(ref_path) as fasta:
        process_snvs(fasta)
        process_svs(fasta)


if __name__ == "__main__":
    main()
