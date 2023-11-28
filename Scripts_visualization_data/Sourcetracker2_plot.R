##Sourcetracker2 plots

#!/usr/bin/env python
# coding: utf-8

from plotnine import *
from plotnine.data import mtcars
import pandas as pd
import os
from plotnine import scales


def stat_out(st2_infile, biom_infile,outfile): #output mix_proportion with taxon and read count
    st2=pd.read_csv(st2_infile,sep="\t",index_col=0)
    sample_order=["VM-28","VM-26","VM-24","VM-22","VM-19","VM-17","VM-15","VM-14","VM-11","VM-3","VM-2"]
    sample_order2=sample_order[::-1]
    st2=st2.reindex(index=sample_order2) #reorder by sample order
    biom=pd.read_csv(biom_infile,sep="\t")
    taxon_count=biom.astype(bool).sum(axis=0)
    read_count=biom.sum(axis=0)
    st2["taxon_count"]=taxon_count
    st2["read_count"]=read_count
    st2=st2.round(3) #To the third decimal place
    st2.to_csv(outfile,index=True, index_label="header")
    print("writing output to " + outfile)

def stat_out_combined(st2_infile, biom_infile,outfile): #output mix_proportion with taxon and read count after collapsing each category
    st2=pd.read_csv(st2_infile,sep="\t",index_col=0)
    sample_order=["VM-28","VM-26","VM-24","VM-22","VM-19","VM-17","VM-15","VM-14","VM-11","VM-3","VM-2"]
    sample_order2=sample_order[::-1]
    st2=st2.reindex(index=sample_order2) #reorder by sample order
    biom=pd.read_csv(biom_infile,sep="\t")
    taxon_count=biom.astype(bool).sum(axis=0)
    read_count=biom.sum(axis=0)
    st2.loc["category"]=pd.Series({"Bos_taurus_dung":"Bovine_rumen","Bovine_rumen":"Bovine_rumen","Cow_Rumen":"Bovine_rumen","Homo_sapiens_fecal":"Human_fecal","Ovis_aries_fecal":"Sheep_fecal","agricultural_soil":"Other_soil","biocrust":"Other_soil","desert_fairy_circles":"Other_soil","desert_varnish":"Other_soil","freshwater_sediment_lake":"Other_soil","gut_metagenome":"Sheep_fecal","melting_permafrost_tundra":"Other_soil","permafrost":"Other_soil","pig_fecal":"Pig_fecal","river_lime_clay":"River_lime_clay","salt_marsh_metagenome":"Other_soil","soil_pH.4.5_agricultural":"Other_soil","soil_pH.5.0_agricultural":"Other_soil","temperate_forest_soil_mineral_layer":"Temperate_forest","temperate_forest_soil_organic_layer":"Temperate_forest","temperate_prairie_grassland":"Temperate_grassland","wetland_acidic_soil":"Wetland_acidic_soil","Unknown":"Unknown"})
    st2=st2.T.groupby("category").sum().T.reset_index().set_index("SampleID")
    ##reorder category
    category_order=["Bovine_rumen","Sheep_fecal","Pig_fecal","Human_fecal","Temperate_forest","Wetland_acidic_soil","River_lime_clay","Temperate_grassland","Other_soil","Unknown"]
    st2=st2.reindex(columns=category_order) #reorder by sample order
    st2["taxon_count"]=taxon_count
    st2["read_count"]=read_count
    st2=st2.round(3) #To the third decimal place
    st2.to_csv(outfile,index=True, index_label="header")
    print("writing output to " + outfile)


files = os.listdir()
mix_files = [file for file in files if file.endswith("_mixing_proportions.txt")]
name_files = [x.split('_mixing_proportions.', 1)[0] for x in mix_files]

#apply stat_out for each file
for name_file in name_files:
    st2_name=name_file + "_mixing_proportions.txt"
    biom_name=name_file + ".tsv"
    outfile_name=name_file + "_result.csv"
    stat_out(st2_name, biom_name,outfile_name)

#apply stat_out combined for each file
for name_file in name_files:
    st2_name=name_file + "_mixing_proportions.txt"
    biom_name=name_file + ".tsv"
    outfile_name=name_file + "_result_combined.csv"
    stat_out_combined(st2_name, biom_name,outfile_name)


for name_file in name_files:
    st2_name=name_file + "_mixing_proportions.txt"
    biom_name=name_file + ".tsv"
    outfile_pdf=name_file + "_result.pdf"
    st2=pd.read_csv(st2_name,sep="\t",index_col=0)
    biom=pd.read_csv(
        biom_name,sep="\t")

    taxon_count=biom.astype(bool).sum(axis=0)
    read_count=biom.sum(axis=0)
    st2.loc["category"]=pd.Series({"Bos_taurus_dung":"Bovine_rumen","Bovine_rumen":"Bovine_rumen","Cow_Rumen":"Bovine_rumen","Homo_sapiens_fecal":"Human_fecal","Ovis_aries_fecal":"Sheep_fecal","agricultural_soil":"Other_soil","biocrust":"Other_soil","desert_fairy_circles":"Other_soil","desert_varnish":"Other_soil","freshwater_sediment_lake":"Other_soil","gut_metagenome":"Sheep_fecal","melting_permafrost_tundra":"Other_soil","permafrost":"Other_soil","pig_fecal":"Pig_fecal","river_lime_clay":"River_lime_clay","salt_marsh_metagenome":"Other_soil","soil_pH.4.5_agricultural":"Other_soil","soil_pH.5.0_agricultural":"Other_soil","temperate_forest_soil_mineral_layer":"Temperate_forest","temperate_forest_soil_organic_layer":"Temperate_forest","temperate_prairie_grassland":"Temperate_grassland","wetland_acidic_soil":"Wetland_acidic_soil","Unknown":"Unknown"})
    st2_c=st2.T.groupby("category").sum().T.reset_index().set_index("SampleID")

    st2_c["taxon_count"]=taxon_count.astype(str)
    st2_c["read_count"]=read_count.astype(str)
    st2_c["taxon_read"]=st2_c["taxon_count"]+ ":" + st2_c["read_count"]
    st2_c2=st2_c.reset_index()

    ##reorder
    sample_order=["VM-28","VM-26","VM-24","VM-22","VM-19","VM-17","VM-15","VM-14","VM-11","VM-3","VM-2"]
    sample_order2=sample_order[::-1]
    st2_c2["SampleID"] = pd.Categorical(st2_c2["SampleID"], categories=sample_order2)

    ##reorder for facet_grid
    st2_c2=st2_c2.sort_values(by = 'SampleID')
    st2_c2["taxon_read"]=pd.Categorical(st2_c2["taxon_read"], categories=st2_c2["taxon_read"].tolist())

    #reshape for plot
    st2_td=st2_c2.melt(id_vars=['SampleID',"taxon_count","read_count","taxon_read"])
    
    #removing below 1% category for plot
    st2_td=st2_td[st2_td["value"]>0.01]

    ##reorder category
    category_order=["Bovine_rumen","Sheep_fecal","Pig_fecal","Human_fecal","Temperate_forest","Wetland_acidic_soil","River_lime_clay","Temperate_grassland","Other_soil","Unknown"]
    st2_td["category"]=st2_td["category"] = pd.Categorical(st2_td["category"], categories=category_order)
    #color
    category_color=["#8B2323", "#EE7600",  "#db7093", "#6495ED","#458B00","#32cd32","#458B74", "#9acd32", "#b8860b","#EBEBEB"]

    result=(ggplot(st2_td,aes(x="SampleID", y="value", fill = "category"))
     + theme_bw()
     + theme(panel_background = element_rect(fill="white"),panel_border = element_rect(colour = "white", fill="white"),panel_grid_major=element_blank(),panel_grid_minor=element_blank()) 
     + geom_bar(stat="identity", position = "fill", width=0.7)
     + scale_fill_manual(values=category_color)
     + facet_grid("taxon_read~", scales="free", space="free")
     + theme(strip_background=element_rect(colour = "white",fill="white"),strip_text = element_text(angle=0,hjust=0, vjust=0.5, size=10)) 
     + theme(axis_text_y = element_text(angle=0, hjust=0.5, vjust=0.5, size=10,color="black"))  
     + theme(axis_ticks = element_blank())  
     + theme(axis_title_x= element_blank(),axis_title_y= element_blank())
     + theme(legend_position = "bottom")
     + ggtitle(name_file)
     + coord_flip()
    )
    print(result)
    result.save(outfile_pdf)
