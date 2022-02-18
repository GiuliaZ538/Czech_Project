#####Mapping damaged reads to consensus reference of mitogenomes########

1. #EXTRACT READS FROM FQ FILES 

#Print sample list based on name files
ll *lca | awk '{print $}' > sample.list
ll *lca | cut -f2 -d: | cut -f2 -d " " > sample.list

##Taxalist
nano name_genus

###CREATE FILES based on SAMPLE.list and Taxa.list
while read -r line
do
arr=($line)
#if [ "${arr[1]}" = "$(basename $folder)" ]
lib=${arr[0]}
echo $lib
cat taxa.list | parallel -j20 "grep {} $lib | cut -f1,2,3,4,5,6,7 -d: > $lib.{}.readID.txt"
done < sample.list

##RMOVE 0 lines in file were not found name_genus
wc -l *.readID.txt| awk '$1 == 0' | awk '{print $2}' > rm.list
cat rm.list | parallel -j20 "rm {}"

#CREATE file tot sequences found per sample
wc -l *.readID.txt| paste > tot_genus_sequences.txt

#EXTRACTION 
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

###IF IN BLANK FOUND SEQUENCE: REMOVE IT from all!!
for infile in /willerslev/users-shared/science-snm-willerslev-mxd270/MetaDMG_CzechALL/VM/LCA_95/*.Canis.fq
do
bname=$(basename $infile)
echo $bname
fastq-grep -v 'GATAGAGTGTTCCTAACCACTTCTGGTTTGGTGCTACCCACCAGCTGGAATTAATTTTT' $bname > grep_$bname
done

OR

fastq-grep -i read_ID $bname > grep2_$bname

##Crete symbolic links of FQ files to new folder
for file in PATHfqfiles_extracted/*.readID.fq
do
bname=$(basename $file)
echo $bname
ln -s $file PATH_new_folder
done

2. ##PREPARE FASTA REFERENCE files
###Check duplicates with Jalview and remove them if any.

###Be sure to have removed tab and , " , ' , : from sequence!
sed "s/:/_/g" in.fa > out.fa
sed "s/ /_/g" in.fa > out.fa
sed "s/;/_/g" in.fa > out.fa 
sed "s/__/_/g" in.fa > out.fa
sed "s/(/./g" in.fa > out.fa 
sed "s/)/./g" in.fa > out.fa

##ALignment 
mafft --thread 40 --leavegappyregion in_to_align.fa > out_aligned.fa

##CHecking Tree of Reference
raxml-ng --search1 --msa in_aligned.fa --model GTR+G --prefix Out_name --threads 2 --seed 2
#If some WARNINGS, run the Reduced.phy files to get clean Tree.

##Build indexes for DATABASE
bowtie2-build -f in_aligned.fa out_db

#you can inspect file the Database with 
#willerslev/software/bowtie2/bowtie2-inspect -s in.fa

3. ##MAPPING readID.fq files with the Reference database just created
for file in *.fq
do
db=OVISclean2c_db
fasta=OVISclean2c_aligned.fa
bowtie2 -x $db -U $file --threads 60 --no-unal | samtools view -bS -> $file.bam
angsd -doFasta 2 -doCounts 1  -i $file.bam -out $file.consensus.fa
done
----------------------
  3.1 ############## for making individual consensus split the fastq from each genome into one read fastq's THIS IS JUST SPLITTING READS TO MAP EACH READ in the TRee.
for file in *.fq
do
split -l 4 $file $file.prefix.fq &
  done

######## then map all individual reads against the reference
for file in *prefix*
  do
db=in_db
fasta=in_aligned.fa
bowtie2 -x $db -U $file --threads 60 --no-unal | samtools view -bS -> $file_split.bam
#angsd -doFasta 2 -doCounts 1  -i $file_split.bam -out $file.consensus.fa
done
-------------------------
  
  4. ####BAMCOV outputs
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

#then (to have both 95 and 98 simil. in alignment) - change name consensus reads 98 sim.
for infile in *_98.fa
do
bname=$(basename $infile)
echo $bname
/willerslev/software/bbmap/rename.sh in=$bname out=$bname.fasta prefix=$bname
done

6. ###Alignment and Tree building
cat *consensus.fasta Ref_not_aligned.fa > Allref_mycons.fa

mafft --thread 40 --leavegappyregion Ref_Cons_reads.fa > Ref_Cons_reads_aligned.fa

raxml-ng --search1 --msa Ref_Cons_reads_aligned.fa --model GTR+G --prefix OVis_CR_Cybt_ALL --threads 2 --seed 2 

####TO GET A FILE WITH REDUCED TREE: WITHOUT identical sequences and/or undetermined columns.  
raxml-ng --search1 --msa reduced.phy --model GTR+G --prefix Reduced --threads 2 --seed 2  
