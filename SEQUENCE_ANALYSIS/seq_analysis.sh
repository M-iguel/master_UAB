#!/bin/bash

FILE=$1
STRFASTA='fasta'
N_LINES=$(bc <<< $(cat $FILE | wc -l)'+'1)

if [[ $FILE == *$STRFASTA ]]
then
	echo 'Running program'
	PDOWN=2
	PUP=2
	echo '' > trimers.txt
	# Get every sequence one by one from the FASTA file
	while [ $PDOWN -lt $N_LINES ]
	do		
		head -n $PDOWN $FILE | tail -n 2 > fasta.txt
		cat fasta.txt | tail -n 1 > temp.txt
		# Split sequences by 3 nt and adding them to trimers.txt file #
		cat temp.txt | tr [:lower:] [:upper:] | fold -w 3 > trimers1.txt ; echo 'I' >> trimers1.txt ; cat trimers1.txt | tr 'I' '\n' >> trimers.txt
		cat temp.txt | fold -w 1 | tail -n +2 | tr -d '\n' | tr [:lower:] [:upper:] | fold -w 3 > trimers2.txt ; echo 'I' >> trimers2.txt ; cat trimers2.txt | tr 'I' '\n' >> trimers.txt
		cat temp.txt | fold -w 1 | tail -n +3 | tr -d '\n' | tr [:lower:] [:upper:] | fold -w 3 > trimers3.txt ; echo 'I' >> trimers3.txt ; cat trimers3.txt | tr 'I' '\n' >> trimers.txt

		# Save every fasta seq in a new file #
		NUMRECORD=$(bc <<< $PDOWN'/'2)
		FILEO=$NUMRECORD"_"$FILE
		cat fasta.txt > $FILEO

		# Trim sequences by right side and save them into new files#
		HEADER=$(cat fasta.txt | head -n 1)
		SEQ=$(cat fasta.txt | tail -n 1)
		TRIMO=$NUMRECORD"_trim_"$FILE
		echo $HEADER > $TRIMO
		echo "${SEQ:0:-20}" >> $TRIMO

		PDOWN=$(bc <<< $PDOWN'+'2)
	done
	# Getting rid of monomers and trimers, sorting them, counting them and sorting the file numerically #
	cat trimers.txt | awk '{ if (length($0) == 3) print }' | sort | uniq -c | sort -rnk1 > $FILE.stats

	# Removing temporal files
	rm fasta.txt temp.txt trimers*
fi



# Same proces than before but including few changes to work with FASTQ files
STRFASTQ='fastq'
if [[ $FILE == *$STRFASTQ ]]
then
	echo 'Running program'
	PDOWN=4
	echo '' > trimers.txt
	# Get every sequence one by one from the FASTQ file
	while [ $PDOWN -lt $N_LINES ]
	do		
		head -n $PDOWN $FILE | tail -n 4 > fastq.txt
		cat fastq.txt | head -n 2 | tail -n 1 > temp.txt
		# Split sequences by 3 nt and adding them to trimers.txt file #
		cat temp.txt | tr [:lower:] [:upper:] | fold -w 3 > trimers1.txt ; echo 'I' >> trimers1.txt ; cat trimers1.txt | tr 'I' '\n' >> trimers.txt
		cat temp.txt | fold -w 1 | tail -n +2 | tr -d '\n' | tr [:lower:] [:upper:] | fold -w 3 > trimers2.txt ; echo 'I' >> trimers2.txt ; cat trimers2.txt | tr 'I' '\n' >> trimers.txt
		cat temp.txt | fold -w 1 | tail -n +3 | tr -d '\n' | tr [:lower:] [:upper:] | fold -w 3 > trimers3.txt ; echo 'I' >> trimers3.txt ; cat trimers3.txt | tr 'I' '\n' >> trimers.txt

		# Save every fasta seq in a new file #
		NUMRECORD=$(bc <<< $PDOWN'/'4)
		FILEO=$NUMRECORD"_"$FILE
		cat fastq.txt > $FILEO

		# Trim sequences by right side and save them into new files#
		HEADER=$(cat fastq.txt | head -n 1)
		SEQ=$(cat fastq.txt | head -n 2 | tail -n 1)
		PLUS=$(cat fastq.txt | head -n 3 | tail -n 1)
		QUALI=$(cat fastq.txt | tail -n 1)
		TRIMO=$NUMRECORD"_trim_"$FILE
		echo $HEADER > $TRIMO
		echo "${SEQ:0:-20}" >> $TRIMO
		echo $PLUS >> $TRIMO
		echo "${QUALI:0:-20}" >> $TRIMO

		PDOWN=$(bc <<< $PDOWN'+'4)
	done
	# Getting rid of monomers and trimers, sorting them, counting them and sorting the file numerically #
	cat trimers.txt | awk '{ if (length($0) == 3) print }' | sort | uniq -c | sort -rnk1 > $FILE.stats

	# Removing temporal files
	rm fasta.txt temp.txt trimers*
fi
echo 'Individualized sequences files and trimmed sequences files have been created'

# MAPING AND ALIGNING READS #

# Creation of index file
bwa index s_cerevisiae.fna
samtools faidx s_cerevisiae.fna

# iteration over FASTA files to align them against the reference genome
if [[ $FILE == *$STRFASTA ]]
then
	for i in $(ls *_trim_$FILE)
	do
	# Alignment of the sequence
	bwa mem s_cerevisiae.fna $i > $i.sam
	done
fi

# iteration over FASTQ files to align them against the reference genome
if [[ $FILE == *$STRFASTQ ]]
then
	for i in $(ls *_trim_$FILE)
	do
	# Convertion of FASTQ file to FASTA file to ease alignment
	sed -n '1~4s/^@/>/p;2~4p' $i > converted.fasta
	cat converted.fasta
	# Alignment of the sequence and saving the outcome into a .sam file
	bwa mem s_cerevisiae.fna converted.fasta > $i.sam
	done
	rm converted.fasta
fi

# Create or empty full_aln.sam file to save all alingment results into the same file
echo '' > full_aln.sam
# iteration over individualized .sam files
for i in $(ls *$FILE.sam)
do
cat $i | tail -n 1 >> full_aln.sam
done

# Sort full_aln.sam file
cat full_aln.sam | sort -k3,3 -k4n,4 | sed -r '/^\s*$/d' > full_aln_sorted.sam

# return number of sequences well alingned
ALN=$(cut full_aln_sorted.sam -f4 | grep -v '*' | wc -l)
echo '##################################'
echo "# $ALN sequences have been aligned #"
echo '##################################'
