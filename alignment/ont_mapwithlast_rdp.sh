#!/usr/bin/env sh

if [ $# -lt 5 ] ; then
    echo "Usage: ont_mapwithlast.py readsfastq readsfmt readgroup /path/to/reference.fasta outdir a m"
    echo "       Map the long reads to a reference using the last program."
    echo "       where"
    echo "           readsfastq is the input file of reads"
    echo "           readsfmt must be fastq"
    echo "           readgroup contains just the ID tag (e.g., htgburnin1)"
    echo "           a is used as the lastal a parameter"
    echo "           m is used as the lastal m parameter"
    echo
    exit 1
fi

readspath=$1
readsfmt=$2
readgroup=$3
ref_fasta=$4
outdir=$5
lastal_a=$6
lastal_m=$7


# Extract useful things from command-line args.
bindir=`dirname $0`
reads_basename=`basename ${readspath} | sed "s,.fastq,,g"`.${lastal_a}.${lastal_m}
ref_dirname=`dirname ${ref_fasta} | sed "s,.fasta,,g"`
ref_basename=`basename ${ref_fasta} | sed "s,.fasta,,g"`

# Create the last db index files (if necessary).
if [ ! -f ${ref_dirname}/${ref_basename}.lastindex.bck ] ; then
    echo "Erro: last db files do not exist ${ref_dirname}/${ref_basename}.lastindex.*"
    exit 2
#    echo "lastdb -Q 0 ${outdir}/${ref_basename}.lastindex ${ref_fasta}"
#          lastdb -Q 0 ${outdir}/${ref_basename}.lastindex ${ref_fasta}
#else
#    echo "Info: ${outdir}/${ref_basename}.lastindex.* files already exist"
fi

# Align the reads with last (if necessary).
if [ ! -f ${outdir}/${reads_basename}.last.txt ] ; then
    if [ $readsfmt = 'fasta' ] ; then
    echo "lastal -m ${lastal_m} -Q 0 -r 1 -q 1 -a ${lastal_a} -b 1 -s 2 ${ref_dirname}/${ref_basename}.lastindex ${readspath} > ${outdir}/${reads_basename}.last.txt"
          lastal -m ${lastal_m} -Q 0 -r 1 -q 1 -a ${lastal_a} -b 1 -s 2 ${ref_dirname}/${ref_basename}.lastindex ${readspath} > ${outdir}/${reads_basename}.last.txt
    elif [ $readsfmt = 'fastq' ] ; then
    echo "lastal -m ${lastal_m} -Q 1 -r 1 -q 1 -a ${lastal_a} -b 1 -s 2 ${ref_dirname}/${ref_basename}.lastindex ${readspath} > ${outdir}/${reads_basename}.last.txt"
          lastal -m ${lastal_m} -Q 1 -r 1 -q 1 -a ${lastal_a} -b 1 -s 2 ${ref_dirname}/${ref_basename}.lastindex ${readspath} > ${outdir}/${reads_basename}.last.txt
    else
        print 'Unrecognised readsfmt *$readsfmt*'
    fi
else
    echo "Info: last output already exists " + ${outdir}/${reads_basename}.last.txt
fi

# Convert output to sorted bam format.
if [ ! -f ${outdir}/${reads_basename}.last.maf.sam ] ; then
    echo "maf-convert.py sam ${outdir}/${reads_basename}.last.txt > ${outdir}/${reads_basename}.last.maf.sam"
          maf-convert.py sam ${outdir}/${reads_basename}.last.txt > ${outdir}/${reads_basename}.last.maf.sam
else
    echo "Info: mapped reads in sam format already exist " + ${outdir}/${reads_basename}.last.maf.sam
fi

# if [ ! -f ${outdir}/${reads_basename}.last.sam ] ; then
#     echo "fixmafsam.py --fastqreads ${readspath} --lastmaf ${outdir}/${reads_basename}.last.txt --lastsam ${outdir}/${reads_basename}.last.maf.sam --readgroup htgburnin1 --fastaref ${ref_fasta} --outsampath ${outdir}/${reads_basename}.last.sam > ${outdir}/${reads_basename}.last.sam.out 2>${outdir}/${reads_basename}.last.sam.err"
#     ${bindir}/fixmafsam_rdp.py \
#         --fastqreads ${readspath} \
#         --lastmaf ${outdir}/${reads_basename}.last.txt \
#         --lastsam ${outdir}/${reads_basename}.last.maf.sam \
#         --readgroup htgburnin1 \
#         --fastaref ${ref_fasta} \
#         --outsampath ${outdir}/${reads_basename}.last.sam \
#         > ${outdir}/${reads_basename}.last.sam.out \
#         2>${outdir}/${reads_basename}.last.sam.err
# else
#     echo "Info: fixed sam file already exists " + ${outdir}/${reads_basename}.last.sam
# fi
#
if [ ! -f ${outdir}/${reads_basename}.last.sorted.bam ] ; then
    echo "samtools view -T ${ref_fasta} -bS ${outdir}/${reads_basename}.last.maf.sam | samtools sort - ${outdir}/${reads_basename}.last.sorted"
          samtools view -T ${ref_fasta} -bS ${outdir}/${reads_basename}.last.maf.sam | samtools sort - ${outdir}/${reads_basename}.last.sorted
else
    echo "Info: sorted bam file already exists " + ${outdir}/${reads_basename}.last.sorted.bam
fi
if [ ! -f ${outdir}/${reads_basename}.last.sorted.bam.bai ] ; then
    echo "samtools index ${outdir}/${reads_basename}.last.sorted.bam"
          samtools index ${outdir}/${reads_basename}.last.sorted.bam

else
    echo "Info: sorted bam index file already exists " + ${outdir}/${reads_basename}.last.sorted.bam.bai
fi

# Print success message
echo "Finished successfully : $*"

