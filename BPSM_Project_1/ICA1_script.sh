#!/bin/bash
destination=$(pwd)  #The user can type in their desired destination to store the raw data and processing files, otherwise the script do that in the current directory

# 0. We will start by copying the raw data files to users local machine
echo "Copying the raw data files to your current directory..."
cp -R /localdisk/data/BPSM/ICA1 $destination 

echo "!Download complete!"
echo "Running FastQC Analysis..."


# 1. Run the fastqc analysis
mkdir $destination/ICA1/processing #Creating a directory to store the fastqc analysis files

no_files=$(ls -l $destination/ICA1/fastq/*.gz | wc -l) #We need to get the number of files we are working with, to know how many fastqc processes to run in parallel
find $destination/ICA1/fastq/*.gz | parallel -j $no_files "fastqc --extract -t 20" #This command will run FASTQC analysis on all the files in 'fastq' directory and create 2 files per each file (html. and .zip)

# !!!Only use the next line if you don't have Parallel/GNU installed. If that applies to you, then replace the line before with this one:
# fastqc --extract -t 2 $destination/ICA1/fastq/*gz 

#find $destination/ICA1/fastq/*.html -delete # We do not need the html files, so we can delete them straight away
#find $destination/ICA1/fastq/*.zip -delete

mv $destination/ICA1/fastq/*fastqc $destination/ICA1/processing # Move the folders into a separate directory


#This block will create a text file which will containt all the PASS/FAIL results from fastqc analysis and delete the source directory
cd $destination/ICA1/processing
touch fastqc_summaries.txt
for dir in */;
do
	cat ./$dir/summary.txt >> ./fastqc_summaries.txt
	rm -rf ./$dir &
done

echo "########################################################"
echo "!FastQC analysis complete!"








# 2. Print the analysis (PASS/FAIL/WARN) so the user can check if that satisfies him enough to continue

#There are 10lines for each sequence, we are mostly interested in "Per Base Sequence Quality" and "Overrepresented sequences". The first will show confidence in the bases retrieved and the latter will let us know if there is any possible contamination.

while true; do
    read -p "Do you wish to check the 'per base sequence quality? (y/n): " yn 
    case $yn in
        [Yy]* ) 
		#Code to count pass/fail/warn for "per base seq quality'
			count_p=0
			count_f=0
			count_w=0
			i=0
			line_interest=2
			while read line
			do
				[ -z "$line" ] && continue
				i=$(( $i + 1 ))
				if test $i == $line_interest 
				then
					if test ${line:0:1} == "P"
					then
					count_p=$(( $count_p+1 ))
					elif test ${line:0:1} == "F"
					then
					count_f=$(( $count_f+1 ))
					elif test ${line:0:1} == "W"
					then
					count_w=$(( $count_w+1 ))
					fi
					line_interest=$(( $line_interest + 10 ))	# 10 lines per sequence, therefore needs to jump +10 lines to keep looking at the per base quality results.
				fi
			done < fastqc_summaries.txt
      echo ""
			echo "!For Per base sequence quality there are:"
			echo "$count_p Passes"
			echo "$count_f Fails"
			echo "$count_w Warnings"
      echo ""
		break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done

while true; do
    read -p "Do you wish to check for the 'Overrepresented sequences? (y/n): " yn
    case $yn in
        [Yy]* ) 
		#Code to count pass/fail/warn for "Overrepresented sequences'
			count_p=0
			count_f=0
			count_w=0
			i=0
			line_interest=9
			while read line
			do
				[ -z "$line" ] && continue
				i=$(( $i + 1 ))
				if test $i == $line_interest
				then
					if test ${line:0:1} == "P"
					then
					count_p=$(( $count_p+1 ))
					elif test ${line:0:1} == "F"
					then
					count_f=$(( $count_f+1 ))
					elif test ${line:0:1} == "W"
					then
					count_w=$(( $count_w+1 ))
					fi
					line_interest=$(( $line_interest + 10 ))	
				fi
			done < fastqc_summaries.txt
      echo ""
			echo "!For the overrepresented sequences there are:!"
			echo "$count_p Passes"
			echo "$count_f Fails"
			echo "$count_w Warnings"
      echo ""
		break;;
        [Nn]* ) break;;
        * ) echo "Please answer yes or no.";;
    esac
done
cd ..
read -p "Press enter to start the alignment process" # I need to give time for the reader to analyse the results, so the script does not continue running straight away.










# 3. Next step is to index the Tcongo reference files so the aligning software can read it
echo "Indexing the reference file..."
bowtie2-build ./Tcongo_genome/* ./fastq/bowtie
echo "!Indexing of the bowtie2 reference file complete!"

# Now let's align the files with bowtie2
cd ./fastq

echo "Aligning the sequences..."
for file in *1.fq.gz # I only need to run it once per file name (Tco-XXXX), not for both file 1 and 2. Therefore I need to skip the loop for every other file ending with .fq.gz
  do
      bowtie2 --very-fast-local -p32 -x ./bowtie -1 ${file:0:9}'1.fq.gz' -2 ${file:0:9}'2.fq.gz' | samtools view -Sb -o ${file:0:8}.bam & #Pipeline to avoid storing the .sam files on local disk, instead let's keep it in RAM to save time 
done
wait 

ct=0
no_files=$(ls -l *1.fq.gz | wc -l)
# I decided to run the indexing on a different loop, so I can run bowtie2 in parallel, which will speed up the process significantly. (couldn't use a single loop as indexing would start on a file that does not exist yet)
for file in *1.fq.gz # I only need to run it once per file name (Tco-XXXX), not for both file 1 and 2. Therefore I need to skip the loop for every other file ending with .fq.gz
  do
      ct=$(( $ct+1 ))
      echo "Sorting and indexing the .bam files... ("$ct"/"$no_files")"
      samtools sort -@$no_files ${file:0:8}.bam > ${file:0:8}_sorted.bam 
      samtools index ${file:0:8}_sorted.bam &
done
wait
mv *sorted.bam* ${destination}/ICA1 #Move all the output files to the same directory where the .bed file is. Preparing for the step 4.







# 4. Generating counts data
cd $destination/ICA1
echo ""
echo "Generating counts data..."
find *bam | parallel 'bedtools intersect -a *.bed -b {} -c >> {}"_bed_res.txt"' #Running bedtools and -c option adds a tab separated column with the number of reads that align to the regions of the genome

# ! Only use the next 'for' loop if you don't have Parallel/GNU installed. If that applies to you, then replace the line before with this loop:
#for f in *.bam
#do
#  bedtools intersect -a *.bed -b $f -c >> ${f:0:8}"_bed_res.txt" & 
#done








# 5.
# This part is all about dividing the Tco.fqfiles into different sample groups (12 total). I've done it in a more complicated way, so it can still work if new clones are added in the future
tail -n +2 "$destination/ICA1/fastq/Tco.fqfiles" >> $destination/ICA1/fastq/Tco.txt

while read line
do 
  awk -v clone=2 'BEGIN{FS="\t";} 
  {
     if($4 == "0")
     {
       if ($5 == "Induced"){print $0 > $clone"_0_Indsum.txt";}
       else print $0 > $clone"_0_Unisum.txt";
     }
     
     else if($4 == "24")
      {
       if ($5 == "Induced"){print $0 > $clone"_24_Indsum.txt";}
       else print $0 > $clone"_48_Unisum.txt";
      } 
      
     else if($4 == "48")
      {
       if ($5 == "Induced"){print $0 > $clone"_48_Indsum.txt";}
       else print $0 > $clone"_48_Unisum.txt";
      }  
  }' 
    
done<$destination/ICA1/fastq/Tco.txt

mv *.txt $destination/ICA1/processing
cd $destination/ICA1/processing

for u in *bed_res.txt
do
rep=$(cat $destination/ICA1/fastq/Tco.txt | grep ${u:0:8} | awk 'BEGIN{FS="\t";}{print $3}')
cat $u | awk -v r=$rep 'BEGIN{FS="\t"; OFS = "\t";}{print $1,$2,$2,$4,$5,$6=$6/r}' >> ${u:0:8}"_sorted.bam_bed_res_rep.txt" #This is a crucial line, because it divides counts of every gene by the number of replicates carried out on the sample
done

for item in *sum.txt
do
  while read line
  do
    echo $line | awk 'BEGIN { ORS = " "; OFS = "\t"}; { print substr($6,0,8)"_sorted.bam_bed_res_rep.txt" }' >> $item"filenames.tsv" #This gives me the files with file names for each sample with their replicates.
    

  done<$item
done

  
for plik in *names.tsv
do  
paste $(cat $plik) >> $plik"_merged.tsv"
done






#Now let's create a tab separated word file
for tsv in *merged.tsv
do
cat $tsv | awk 'BEGIN{ FS="\t"; OFS = "\t";}{print $1=$6+$12+$18+$24+$30+$36+$42+$48,$4,$5,$6,$12,$18,$24,$30,$36,$42,$48,$54}' >> ${tsv:0:12}"prefi.tsv" # I need that middle step in order to no divide the final (count/replicate) value ($1) by the no. of member in the group which I will work out in the next while loop
done



for prefi in *prefi.tsv
do
  while read line
  do
    sed 's/\t\+/\t/g;s/^\t//' >>$prefi"_full.tsv"
  done<$prefi
done


for full in *full.tsv
do
  
  fields=$(head -1 $full | awk --field-separator="\t" '{ print NF }') #Number of fields per each file
  members=$(expr $fields - 4)
  
  cat $full | awk -v m=$members -v fu=$full 'BEGIN{ FS="\t"; OFS = "\t";}{print $1=$1/m,$2,$3 > (fu"mean_expr.txt")}' #This line is outputting the final tab separated txt file
  
done

mkdir $destination/ICA1/Step6
mv *mean_expr.txt $destination/ICA1/Step6
cd $destination/ICA1/Step6

echo "! Mean Averages of the counts per gene calculated! Correspodning files created!" 
echo "You can find them in: "$destination"/ICA1/Step6" 
echo ""
echo "Unfortunately, the script was not finished before the deadline :("
echo ""

#Step 6

#Unfortunately, I did not anticipate that part 5 will take me so long. For that reason I did not have time to code the last part of the script




















