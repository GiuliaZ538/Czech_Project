######pathPhynder commands used 
### bwa mapping of samples to cons from references ###
bwa aln -l 1024 -n 0.001 -t 10 /path_to_folder/Consensus_Mafft_All_references.fa /path_to_folder/library.fq
| bwa samse /path_to_folder/Consensus_Mafft_All_references.fa  - /path_to_folder/library.fq
| samtools view -F 4 -q 25 -@ 10 -uS - | samtools sort -@ 10 -o library.sort.bam

### PathPhynder command1: ###
phynder -B -o /path_to_folder/branches.snp /path_to_folder/Mafft_All_reference.nwk /path_to_folder/Mafft_All_references_fixed_consensus.vcf

### PathPhynder command2: ###
pathPhynder -s prepare -i /path_to_folder/Mafft_All_reference.nwk -p taxa_pathphynder_tree -f /path_to_folder/branches.snp -r /path_to_folder/Consensus_Mafft_All_references.fa

### PathPhynder command3: ###
pathPhynder -s all -t 100 -m transversions -i /path_to_folder/Mafft_All_references.nwk -p /path_to_folder/taxa_pathphynder_tree -l /path_to_folder/bamlist.txt -r /path_to_folder/Consensus_Mafft_All_references.fa

### PathPhynder command3: ###transitions and transvertions
pathPhynder -s all -t 100 -i /path_to_folder/Mafft_All_references.nwk -p /path_to_folder/taxa_pathphynder_tree -l /path_to_folder/bamlist.txt -r /path_to_folder/Consensus_Mafft_All_references.fa
