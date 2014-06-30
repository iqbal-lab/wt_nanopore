#! /bin/bash 

## Takes the output of lastal and converts it into a samfile. 
## Must have maf-convert.py in PATH (it comes with the last install)

## USAGE : ./last2Sam.sh /path/to/lastal/output/file.txt


if [ "$1" == "" ] || ["$1" == "-h"] || ["$1" == "-help" ] ; then
    echo "Usage: `basename $0` /path/to/lastal/output/file.txt"
    exit 0
fi

# if [ ! -x  "maf-convert.py"  ]; then
# 	echo "Can't find maf-convert.py. Add it to your path (it should have come with the last install)."
# 	exit 0
# fi


f=$1
b=$(basename $f)
directory=$(dirname ${f})
filename="${b%.*}"
echo "Converting $f to samfile"
~/tools/banyan_last/scripts/maf-convert.py sam $f > "$directory"/"$filename".sam
echo "Finished converting $f to samfile. Output in $directory/$filename.sam"
echo "Converting $directory/$filename.sam to sorted bamfile"
~/tools/samtools-0.1.17/samtools view -T /Net/wombat/dipro/mmm/data/Nanopore/processed_data/Run_1/ref/BX571856.1.fasta -bS "$directory"/"$filename".sam | ~/tools/samtools-0.1.17/samtools sort - "$directory"/"$filename".sam.sorted
~/tools/samtools-0.1.17/samtools index "$directory"/"$filename".sam.sorted.bam
mv "$directory"/"$filename".sam.sorted.bam "$directory"/"$filename".sorted.bam
mv "$directory"/"$filename".sam.sorted.bam.bai "$directory"/"$filename".sorted.bam.bai
echo "Finished converting $directory/$filename.sam to sorted bamfile. Output in $directory/$filename.sorted.bam " 
