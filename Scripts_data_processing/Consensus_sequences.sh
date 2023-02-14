#####Read extraction by readID and consensus fasta files

##Tools
bowtie2-build, bowtie2-align (v. 2.3.2)
mafft (v7.427)
raxml-ng (v. 0.8.1 BETA)
angsd (v. 0.928)
Geneious Prime (v. 2021.2.2 )

1. ##Create fastq file with taxa of interest

#Print sample list based on name files
ll *lca | awk '{print $}' > sample.list
ll *lca | cut -f2 -d: | cut -f2 -d " " > sample.list

#Taxalist
nano name_genus or name_family

#Create txt file with all readIDs
while read -r line
do
arr=($line)
#if [ "${arr[1]}" = "$(basename $folder)" ]
lib=${arr[0]}
echo $lib
cat taxa.list | parallel -j20 "grep {} $lib | cut -f1,2,3,4,5,6,7 -d: > $lib.{}.readID.txt"
done < sample.list

#Remove file were not found name_genus
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

#Create symbolic links of fastq files to new folder
for file in PATHfqfiles_extracted/*.readID.fq
do
bname=$(basename $file)
echo $bname
ln -s $file PATH_new_folder
done

2. ##Prepare fasta reference files
#Be sure to have removed tab and , " , ' , : from sequence!
sed "s/:/_/g" in.fa > out.fa
sed "s/ /_/g" in.fa > out.fa
sed "s/;/_/g" in.fa > out.fa 
sed "s/__/_/g" in.fa > out.fa
sed "s/(/./g" in.fa > out.fa 
sed "s/)/./g" in.fa > out.fa

#Alignment references
mafft --thread 40 --leavegappyregion in_to_align.fa > out_aligned.fa

#if needed a consensus of the references(Geneious software: the most common base rule, and ignore gaps)

#Checking tree of of the reference
raxml-ng --search1 --msa in_aligned.fa --model GTR+G --prefix Out_name --threads 2 --seed 2
#If some warnings, run the Reduced.phy files to get a clean Tree.

#Build indexes for database
bowtie2-build -f in_aligned.fa out_db
bowtie2-inspect -s in.fa #you can inspect file the Database 

3. ##Mapping readID.fq files to the reference database 
for file in *.fq
do
db=Reference_db
fasta=Reference_aligned.fa
bowtie2 -x $db -U $file --threads 60 --no-unal | samtools view -bS -> $file.bam
angsd -doFasta 2 -doCounts 1  -i $file.bam -out $file.consensus.fa
done

4. ##Bamcov outputs
for file in *.bam
do
samtools sort -O BAM -o sort_$file $file
done

#To get bamcov histograms
for file in sort*.fq.bam
do
echo $file
bamcov -m -w0 $file | paste > bamcov_$file.txt 
done

#To get Bamcov tables
for file in sort_*.bam
do
bname=$(basename $file)
echo $bname
bamcov -Q 30 $file | paste > bamcov_table30filter_$file.txt
done

5. ##GeT the consensus reads 
for file in *.bam
do
angsd -doFasta 2 -doCounts 1 -i $file -out $file.consensus.fa
done

#Gunzip
for file in *.consensus.fa.fa.gz
do
bname=$(basename $file)
echo $bname
gunzip $file 
done

#Delete Size 0
for file in *.consensus.fa.fa
do
bname=$(basename $file)
echo $bname
find $file -size 0 -delete
done

#Rename headers consensus reads
for infile in *.fa.fa
do
bname=$(basename $infile)
bname2=$(ll $bname | cut -f9-10 -d-)
echo $bname
echo $bname2
/willerslev/software/bbmap/rename.sh in=$bname out=$bname.fasta prefix=$bname2
done

