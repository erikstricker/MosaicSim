import sys
import pysam
from pysam import depth
import argparse
import os
import numpy as np
from numpy.random import choice, uniform, seed as npseed
from numpy import random as nran
from math import ceil, log

# Initialize parser
parser = argparse.ArgumentParser(description="TykeVarSimulator: A tool for simulating variants")

# Define required arguments
parser.add_argument("-i", "--input", required=True, help="Path to input BAM file")
parser.add_argument("-T", "--reference", required=True, help="Path to reference genome")
parser.add_argument("-o", "--output", required=True, help="Output file prefix")

# Optional arguments
parser.add_argument("-s", "--seed", type=int, required=False, help="Seed for random number generation")
parser.add_argument("-minAFsv", "--minimum-allele-frequency-sv", type=float, required=False, help="Minimum allele frequency for structural variants")
parser.add_argument("-maxAFsv", "--maximum-allele-frequency-sv", type=float, required=False, help="Maximum allele frequency for structural variants")
parser.add_argument("-minAFsnv", "--minimum-allele-frequency-snv", type=float, required=False, help="Minimum allele frequency for single nucleotide variants")
parser.add_argument("-maxAFsnv", "--maximum-allele-frequency-snv", type=float, required=False, help="Maximum allele frequency for single nucleotide variants")

# Newly added arguments
parser.add_argument("-numsv", "--number-of-svs", type=int, required=False, help="Number of structural variants to simulate")
parser.add_argument("-numsnv", "--number-of-snvs", type=int, required=False, help="Number of single nucleotide variants to simulate")
parser.add_argument("-maxsnvl", "--maximum-snv-length", type=int, required=False, help="Maximum length of single nucleotide variants indels")
parser.add_argument("-minsvl", "--minimum-sv-length", type=int, required=False, help="Minimum length of structural variants")
parser.add_argument("-maxsvl", "--maximum-sv-length", type=int, required=False, help="Maximum length of structural variants")
parser.add_argument("-sub", "--substitution-rate", type=float, required=False, help="Probability of producing a SNP vs indel in SNV")
parser.add_argument("-insdelsnv", "--insdel-snv-rate", type=float, required=False, help="Probability of producing an insertion vs deletion in SNV")
parser.add_argument("-insdel", "--insdel-sv-rate", type=float, required=False, help="Probability of producing an insertion vs deletion in SV")

# Parse arguments
args = parser.parse_args()

# Check for missing arguments
if len(sys.argv) < 2:
    print("Usage: python vcfgen.py -i <path_to_bam> -T <path_to_ref> -o <output_path_prefix> [optional arguments]")
    print("\nExample: python vcfgen.py -i chr22.bam -T hs37d5.fa -o chr22")
    print("The above generates a chr22_SV.vcf and chr22_SNV.vcf file")
    exit(0)


# Access variables
bam_path = args.input
ref_path = args.reference
output_prefix = args.output

## Test variables
if args.seed is None:
    seed=0
else:
    seed=args.seed

if args.minimum_allele_frequency_sv is None:
    minAFsv=0.01 # minimum allele frequency SV 
else:
    minAFsv=args.minimum_allele_frequency_sv  

if args.maximum_allele_frequency_sv is None:
    maxAFsv=0.05  # maximum allele frequency SV
else:
    maxAFsv=args.maximum_allele_frequency_sv 

if args.minimum_allele_frequency_snv is None:
    minAFsnv=0.01 # minimum allele frequency SNV 
else:
    minAFsnv=args.minimum_allele_frequency_snv  

if args.maximum_allele_frequency_snv is None:
    maxAFsnv=0.05  # maximum allele frequency SNV
else:
    maxAFsnv=args.maximum_allele_frequency_snv 

if args.number_of_svs is None:
    numsv=50  # number of SVs to simulate
else:
    numsv=args.number_of_svs 

if args.number_of_snvs is None:
    numsnv=200  # number of SNVs to simulate
else:
    numsnv=args.number_of_snvs  

if args.maximum_snv_length is None:
    maxsnvl=100  # maximum SNV length
else:
    maxsnvl=args.maximum_snv_length 

if args.minimum_sv_length is None:
    minsvl=50  # minimum SV length
else:
    minsvl=args.minimum_sv_length 

if args.maximum_sv_length is None:
    maxsvl=10000  # maximum SV length
else:
    maxsvl=args.maximum_sv_length 

if args.substitution_rate is None:
    sub=1  # probability of producing a SNP vs indel in SNV
else:
    sub=args.substitution_rate 

if args.insdel_snv_rate is None:
    insdelsnv=0.5  # probability of producing INS vs DEL in SNV
else:
    insdelsnv=args.insdel_snv_rate 

if args.insdel_sv_rate is None:
    insdel=0.7  # probability of producing INS vs DEL in SV
else:
    insdel=args.insdel_sv_rate 


SVvcf=f"{output_prefix}_SV.vcf" # name of output vcf file for SV
SNVvcf=f"{output_prefix}_SNV.vcf" # name of output vcf file for SNV

CHROMS = [str(i) for i in range(1,23)] + ["X","Y"]
CHROMS += [f'chr{c}' for c in CHROMS]

def get_chrom_lengths(bam_path, ref_chroms=CHROMS):
    with pysam.AlignmentFile(bam_path) as bam:
        chromol = dict((c,l) for c,l in zip(bam.references, bam.lengths) if c in ref_chroms)
    return chromol, tuple(chromol.keys())

def genlocSNV(num,bam_path,mincov=20):
    npseed(seed)
    chromol, chrom = get_chrom_lengths(bam_path)
    locations=[]
    for i in range(num):
        while True:
            ranchrom=nran.choice(chrom)
            loc=str(nran.randint(0,chromol[ranchrom]))
            try:
                cover=int(depth(bam_path,'-r',ranchrom+":"+loc+"-"+loc).rstrip("\n").split("\t")[-1])
            except ValueError:
                continue
            if cover>=mincov:
                break
        locations.append(tuple([ranchrom,loc,cover]))
    return tuple(locations)

def genlocSV(num,bam_path,mincov=20):
    npseed(seed)
    chromol, chrom = get_chrom_lengths(bam_path)
    locations=[]
    for i in range(num):
        while True:
            ranchrom=nran.choice(chrom)
            loc=str(nran.randint(0,chromol[ranchrom]))
            for i in locations:
                if i[0]==ranchrom and abs(int(i[1]))-int(loc)<50000:
                    continue
            try:
                cover=int(depth(bam_path,'-r',ranchrom+":"+loc+"-"+loc).rstrip("\n").split("\t")[-1])
            except ValueError:
                continue
            if cover>=mincov:
                break
        locations.append(tuple([ranchrom,loc,cover]))
    return tuple(locations)

def genseq(minl,maxl):
    npseed(seed)
    nucl=tuple(["A","T","C","G"])
    res=''
    for i in range(choice(range(minl,maxl+1))):
        res+=choice(nucl)
    return res

def gensnps(maxsnvl=maxsnvl,sub=sub,insdelsnv=insdelsnv, snplist=[]):
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
    


def main():

    svloc=genlocSV(numsv,bam_path,ceil(1/minAFsv))
    snvloc=genlocSNV(numsnv,bam_path,ceil(1/minAFsnv))
    
    chromol, chrom = get_chrom_lengths(bam_path)

    if numsv > 0:
        #print(snvloc)
        vcfsv=[
            '##fileformat=VCFv4.2',
            '##ALT=<ID=INS,Description="Insertion">',
            '##ALT=<ID=DEL,Description="Deletion">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">',
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'
        ]

        # Dynamically add contig lines
        vcfsv.extend([f'##contig=<ID={c},length={l}>' for c, l in chromol.items()])

        # Add VCF column headers
        vcfsv.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))

        insertnum=1
        delnum=1
        for i in svloc:
            draw = choice(tuple(['in','del']), 1, p=[insdel,1-insdel])    
            if draw=='in':
                seq=genseq(minsvl,maxsvl)
                vcfsv.append(str(i[0])+"\t"+str(i[1])+"\tHackIns"+str(insertnum)+"\tN\t"+seq+"\t60\tPASS\tPRECISE;SVTYPE=INS;SVLEN="+str(len(seq))+";END="+str(int(i[1])+1)+";AF="+str(round(uniform(minAFsv,maxAFsv),2))+"\tGT:GQ\t0/0:60")
                #PRECISE;SVTYPE=INS;SVLEN=333;END=748218 AF \t GT:GQ:DR:DV \t	0/0:28:28:5
                insertnum+=1
            else:
                dellen=choice(range(minsvl,maxsvl))
                vcfsv.append(str(i[0])+"\t"+str(i[1])+"\tHackDel"+str(delnum)+"\tN\t<DEL>\t60\tPASS\tPRECISE;SVTYPE=DEL;SVLEN=-"+str(dellen)+";END="+str(int(i[1])+dellen)+";AF="+str(round(uniform(minAFsv,maxAFsv),2))+"\tGT:GQ\t0/0:60")
                delnum+=1

        # Ensure parent directories exist
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

        with open(SVvcf,"w") as f:
            for i in tuple(vcfsv)[:-1]:
                f.write(i+'\n')
            f.write(tuple(vcfsv)[-1])
        f.close()

    if numsnv > 0:
        vcfsnv=[
            '##fileformat=VCFv4.2',
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
            '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'
        ]

        # Dynamically add contig lines
        vcfsnv.extend([f'##contig=<ID={c},length={l}>' for c, l in chromol.items()])

        # Add VCF column headers
        vcfsnv.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']))

        snps=getrefsnp(ref_path,snvloc)
        
        
        for i in gensnps(maxsnvl=maxsnvl,sub=sub, snplist=snps):
            AFnum=(round(uniform(minAFsnv,maxAFsnv),2))
            readnum=ceil(AFnum*i[2])
            vcfsnv.append(str(i[0])+'\t'+str(i[1])+'\t.\t'+str(i[3])+'\t'+str(i[4])+"\t1500\tPASS\tAF="+str(AFnum)+"\tGT:AD:DV\t0/0:"+str(i[2]-readnum)+":"+str(readnum))
        with open(SNVvcf,"w") as f:
            for i in tuple(vcfsnv)[:-1]:
                f.write(i+'\n')
            f.write(tuple(vcfsnv)[-1])
        f.close()
if __name__=="__main__":
    main()
