########python script Sourcetracker2##########

#making sourcetracker input files
import pandas as pd

####parameter
infile1="subdata.csv"
infile2="VM.csv"
outfile1s="species.tsv"
outfile1g="genus.tsv"
outfile2s="species_z2_D5.tsv"
outfile2g="genus_z2_D5.tsv"
outfile3s="species_lessD5.tsv"
outfile3g="genus_lessD5.tsv"
outfile4s="species_D5.tsv"
outfile4g="genus_D5.tsv"


#output column names
col_names=["tax_name","sample","N_reads"]

######

a1=pd.read_csv(infile1)
a2=pd.read_csv(infile2)

#species level
a1_s=a1[a1["tax_rank"]=="species"]
a2_s=a2[a2["tax_rank"]=="species"]

#genus level
a1_g=a1[a1["tax_rank"]=="genus"]
a2_g=a2[a2["tax_rank"]=="genus"]

##nofilter
b1_s=a1_s[col_names]
b2_s=a2_s[col_names]
b_s_all=pd.concat([b1_s,b2_s])
c_s=pd.pivot_table(b_s_all,index="tax_name",columns="sample",values="N_reads",fill_value=0)
c_s.to_csv(outfile1s,index=True,sep="\t")

b1_g=a1_g[col_names]
b2_g=a2_g[col_names]
b_g_all=pd.concat([b1_g,b2_g])
c_g=pd.pivot_table(b_g_all,index="tax_name",columns="sample",values="N_reads",fill_value=0)
c_g.to_csv(outfile1g,index=True,sep="\t")

##Baysian Z>2 & Baysian Bayesian_D_max>5%
a2_s_ZD=a2_s[(a2_s["Bayesian_D_max"]>0.05) & (a2_s["Bayesian_z"]>2)]
b1_s_ZD=a1_s[col_names]
b2_s_ZD=a2_s_ZD[col_names]
b_s_ZD_all=pd.concat([b1_s_ZD,b2_s_ZD])
c_s_ZD=pd.pivot_table(b_s_ZD_all,index="tax_name",columns="sample",values="N_reads",fill_value=0)
c_s_ZD.to_csv(outfile2s,index=True,sep="\t")

a2_g_ZD=a2_g[(a2_g["Bayesian_D_max"]>0.05) & (a2_g["Bayesian_z"]>2)]
b1_g_ZD=a1_g[col_names]
b2_g_ZD=a2_g_ZD[col_names]
b_g_ZD_all=pd.concat([b1_g_ZD,b2_g_ZD])
c_g_ZD=pd.pivot_table(b_g_ZD_all,index="tax_name",columns="sample",values="N_reads",fill_value=0)
c_g_ZD.to_csv(outfile2g,index=True,sep="\t")

#Baysian Bayesian_D_max=<5%
a2_s_LD=a2_s[(a2_s["Bayesian_D_max"]<=0.05)]
b1_s_LD=a1_s[col_names]
b2_s_LD=a2_s_LD[col_names]
b_s_LD_all=pd.concat([b1_s_LD,b2_s_LD])
c_s_LD=pd.pivot_table(b_s_LD_all,index="tax_name",columns="sample",values="N_reads",fill_value=0)
c_s_LD.to_csv(outfile3s,index=True,sep="\t")

a2_g_LD=a2_g[(a2_g["Bayesian_D_max"]<=0.05)]
b1_g_LD=a1_g[col_names]
b2_g_LD=a2_g_LD[col_names]
b_g_LD_all=pd.concat([b1_g_LD,b2_g_LD])
c_g_LD=pd.pivot_table(b_g_LD_all,index="tax_name",columns="sample",values="N_reads",fill_value=0)
c_g_LD.to_csv(outfile3g,index=True,sep="\t")

#Baysian Bayesian_D_max>5%
a2_s_D5=a2_s[(a2_s["Bayesian_D_max"]>0.05)]
b1_s_D5=a1_s[col_names]
b2_s_D5=a2_s_D5[col_names]
b_s_D5_all=pd.concat([b1_s_D5,b2_s_D5])
c_s_D5=pd.pivot_table(b_s_D5_all,index="tax_name",columns="sample",values="N_reads",fill_value=0)
c_s_D5.to_csv(outfile4s,index=True,sep="\t")

a2_g_D5=a2_g[(a2_g["Bayesian_D_max"]>0.05)]
b1_g_D5=a1_g[col_names]
b2_g_D5=a2_g_D5[col_names]
b_g_D5_all=pd.concat([b1_g_D5,b2_g_D5])
c_g_D5=pd.pivot_table(b_g_D5_all,index="tax_name",columns="sample",values="N_reads",fill_value=0)
c_g_D5.to_csv(outfile4g,index=True,sep="\t")


###count the number of read and taxon


####python script end##################
#sourcetracker2 analysis with downsampling(--sink_rarefaction_depth 100)

for table in *.tsv;
do
name=$(basename "${table}" .tsv)
#########making biom table
conda activate biom
biom convert -i ${table} -o table.from_txt_hdf5.biom --table-type="OTU table" --to-hdf5
echo created_${name}_biom_table
#########sourcetracker
conda activate st2
sourcetracker2 gibbs -i table.from_txt_hdf5.biom -m map_rm.txt -o result_${name} --sink_rarefaction_depth 100
echo created_${name}_sourcetracker_result
cp result_${name}/mixing_proportions.txt ${name}_mixing_proportions.txt
done
