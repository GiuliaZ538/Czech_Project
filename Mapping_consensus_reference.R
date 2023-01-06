#####Mapping damaged reads to consensus reference of mitogenomes########

##Tools
bowtie2-build, bowtie2-align (v. 2.3.2)
mafft (v7.427)
raxml-ng (v. 0.8.1 BETA)
angsd (v. 0.928)
Geneious Prime (v. 2021.2.2 )
Jalview (v. 2.11.0)

1. ###Create fastq file with taxa of interest

#Print sample list based on name files
ll *lca | awk '{print $}' > sample.list
ll *lca | cut -f2 -d: | cut -f2 -d " " > sample.list

##Taxalist
nano name_genus

##Create txt file with all readIDs
while read -r line
do
arr=($line)
#if [ "${arr[1]}" = "$(basename $folder)" ]
lib=${arr[0]}
echo $lib
cat taxa.list | parallel -j20 "grep {} $lib | cut -f1,2,3,4,5,6,7 -d: > $lib.{}.readID.txt"
done < sample.list

##Remove file were not found name_genus
wc -l *.readID.txt| awk '$1 == 0' | awk '{print $2}' > rm.list
cat rm.list | parallel -j20 "rm {}"

#Create file tot sequences found per sample
wc -l *.readID.txt| paste > tot_genus_sequences.txt

#Create fastq from readIDs
for infile in *readID.txt
do
bname=$(basename $infile)
bname1=$(echo $bname | sed 's/.readID.txt*/.fq/')   
bname2=$(echo $bname | cut -f2 -d: | cut -f2 -d" " | cut -f1 -d.)
echo $bname1
echo $bname2
seqtk subseq fq_link/*$bname2*rmdup.fq $bname > $bname1 &
  #seqtk seq -a $bname1 > $bname1.fa ##NOTHING IN IT conversion into fasta didnt work
  done

##Crete symbolic links of fastq files to new folder
for file in PATHfqfiles_extracted/*.readID.fq
do
bname=$(basename $file)
echo $bname
ln -s $file PATH_new_folder
done

2. ##Prepare fasta reference files
#Check duplicates with Jalview software and remove them if any.

###Be sure to have removed tab and , " , ' , : from sequence!
sed "s/:/_/g" in.fa > out.fa
sed "s/ /_/g" in.fa > out.fa
sed "s/;/_/g" in.fa > out.fa 
sed "s/__/_/g" in.fa > out.fa
sed "s/(/./g" in.fa > out.fa 
sed "s/)/./g" in.fa > out.fa

#Alignment references
mafft --thread 40 --leavegappyregion in_to_align.fa > out_aligned.fa

#if needed a consensus of the references
Make Consensus Sequence with Geneious software: use the most common base, and ignore gaps

#Checking tree of of the reference
raxml-ng --search1 --msa in_aligned.fa --model GTR+G --prefix Out_name --threads 2 --seed 2
#If some warnings, run the Reduced.phy files to get a clean Tree.

#Build indexes for database
bowtie2-build -f in_aligned.fa out_db
bowtie2-inspect -s in.fa #you can inspect file the Database 

3. ##Mapping readID.fq files to the reference database 
for file in *.fq
do
db=OVISclean2c_db
fasta=OVISclean2c_aligned.fa
bowtie2 -x $db -U $file --threads 60 --no-unal | samtools view -bS -> $file.bam
angsd -doFasta 2 -doCounts 1  -i $file.bam -out $file.consensus.fa
done
----------------------
3.1 # for making individual consensus split the fastq from each genome into one read fastq
for file in *.fq
do
split -l 4 $file $file.prefix.fq &
  done

# then map all individual reads against the reference
for file in *prefix*
  do
db=in_db
fasta=in_aligned.fa
bowtie2 -x $db -U $file --threads 60 --no-unal | samtools view -bS -> $file_split.bam
angsd -doFasta 2 -doCounts 1  -i $file_split.bam -out $file.consensus.fa
done
-------------------------
4. ##Bamcov outputs
for file in *.bam
do
samtools sort -O BAM -o sort_$file $file
done

##To get bamcov histograms
for file in sort*.fq.bam
do
echo $file
bamcov -m -w0 $file | paste > bamcov_$file.txt 
done

##To get Bamcov tables
for file in sort_*.bam
do
bname=$(basename $file)
echo $bname
bamcov -Q 30 $file | paste > bamcov_table30filter_$file.txt
done

5. ####GET the consensus reads 
for file in *.bam
do
angsd -doFasta 2 -doCounts 1 -i $file -out $file.consensus.fa
done

#GUNZIP MANY READS
for file in *.consensus.fa.fa.gz
do
bname=$(basename $file)
echo $bname
gunzip $file 
done

##DELETE SIZE 0
for file in /willerslev/users-shared/science-snm-willerslev-mxd270/MetaDMG_CzechALL/LCA/VM/Bovinae/*.consensus.fa.fa
do
bname=$(basename $file)
echo $bname
find $file -size 0 -delete
done

##Rename headers consensus reads
for infile in *.fa.fa
do
bname=$(basename $infile)
bname2=$(ll $bname | cut -f9-10 -d-)
echo $bname
echo $bname2
/willerslev/software/bbmap/rename.sh in=$bname out=$bname.fasta prefix=$bname2
done

6. ###Alignment and Tree building/CHECK BEAST
cat *consensus.fasta Ref_not_aligned.fa > Allref_mycons.fa

mafft --thread 40 --leavegappyregion Ref_Cons_reads.fa > Ref_Cons_reads_aligned.fa

raxml-ng --search1 --msa Ref_Cons_reads_aligned.fa --model GTR+G --prefix OVis_CR_Cybt_ALL --threads 2 --seed 2 

###REFERENCES Tools

