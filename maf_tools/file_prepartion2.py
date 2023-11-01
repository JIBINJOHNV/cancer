
import pandas as pd
import os,glob
import subprocess

samples_1=['T_01','T_02','T_03','T_04','T_05','T_06','T_07','T_08','T_09','T_10','T_11','T_13','Met_13','T_14','T_15','T_16',
           'T_18','T_19','Met_19','T_20','Met_20_a','Met_20_b', 'T_23', 'Met_24', 'Met_25']

vcf_files1=[x+".vcf" for x in samples_1 ]

sample_names_1=['Sample_01_Tumour_01', 'Sample_02_Tumour_02', 'Sample_03_Tumour_03', 'Sample_04_Tumour_04', 'Sample_05_Tumour_05', 'Sample_06_Tumour_06', 
                'Sample_07_Tumour_07', 'Sample_08_Tumour_08', 'Sample_09_Tumour_09', 'Sample_10_Tumour_10', 'Sample_11_Tumour_11', 'Sample_13_Tumour_13', 
                'Sample_13_Met_13', 'Sample_14_Tumour_14', 'Sample_15_Tumour_15', 'Sample_16_Tumour_16', 'Sample_18_Tumour_18', 'Sample_19_Tumour_19', 
                'Sample_19_Met_19','Sample_20_Tumour_20','Sample_20_Met_20_a','Sample_20_Met_20_b','Sample_23_Tumour_23','Sample_24_Met_24','Sample_25_Met_25']

samples_1_dict=dict(zip(vcf_files1,sample_names_1))


for vcf,name in samples_1_dict.items():
        os.system(f'''perl /mnt/disks/sdb/VCF-Ap/vcf2maf/vcf2maf.pl --input-vcf {vcf} --output-maf {name}.maf --ref-fasta /mnt/disks/sdb/Resourses/GATK_Rsource/Homo_sapiens_assembly38.fasta --tumor-id {name} --vep-path /home/jjohn41/anaconda3/envs/mamba/bin/ --vep-data /mnt/disks/sdb/Resourses/ --ncbi-build GRCh38 ''')


samples_2=['T_12','T_17','T_21','T_22']
vcf_files2=[x+".vcf" for x in samples_2 ]
sample_names_2=[['Sample_12_Normal_12','Sample_12_Tumour_12'],['Sample_17_Normal_17', 'Sample_17_Tumour_17'],['Sample_21_Normal_21','Sample_21_Tumour_21'],['Sample_22_Normal_22','Sample_22_Tumour_22']]

samples_2_dict=dict(zip(vcf_files2,sample_names_2))

for vcf,name in samples_2_dict.items():
        os.system(f"sed -i 's/\\r/\\n/' {vcf}")
        os.system(f"sed -i '/^$/d' {vcf}")
        vep=vcf.replace(".vcf","")
        os.system(f"rm {vep}.vep.vcf")
        os.system(f'''perl /mnt/disks/sdb/VCF-Ap/vcf2maf/vcf2maf.pl --input-vcf {vcf} --output-maf {name[1]}.maf \
                  --ref-fasta /mnt/disks/sdb/Resourses/GATK_Rsource/Homo_sapiens_assembly38.fasta \
                  --tumor-id {name[1]} --normal-id {name[0]} --vep-path /home/jjohn41/anaconda3/envs/mamba/bin/ --vep-data /mnt/disks/sdb/Resourses/ --ncbi-build GRCh38 ''')


### Control samples
for vcf,name in samples_2_dict.items():
        os.system(f"sed -i 's/\\r/\\n/' {vcf}")
        os.system(f"sed -i '/^$/d' {vcf}")
        vep=vcf.replace(".vcf","")
        os.system(f"rm {vep}.vep.vcf")
        os.system(f'''perl /mnt/disks/sdb/VCF-Ap/vcf2maf/vcf2maf.pl --input-vcf {vcf} --output-maf {name[0]}.maf \
                  --ref-fasta /mnt/disks/sdb/Resourses/GATK_Rsource/Homo_sapiens_assembly38.fasta \
                  --normal-id {name[0]} --vep-path /home/jjohn41/anaconda3/envs/mamba/bin/ --vep-data /mnt/disks/sdb/Resourses/ --ncbi-build GRCh38 ''')



### Mergee maf files
df=pd.read_csv("Sample_13_Tumour_13.maf",sep="\t",quotechar='"',quoting=3,skiprows=1 )
df.to_csv("Sample_13_Tumour_13.maf",sep="\t")



## Creating Tumor MAF file by merging all tumor sample wise MAF files
maf_files=glob.glob("*.maf")
maf_files=sorted(maf_files)
master_df=pd.DataFrame()

normal_files=[x for x in maf_files if "_Normal_" in x]
tumor_files=[x for x in maf_files if "_Normal_" not in x]

for file in tumor_files:
    df=pd.read_csv(file, sep="\t", quotechar='"',quoting=3,skiprows=1)
    df["MaxMAF"]=df[[x for x in df.columns  if "_AF" in x]].max(axis=1)
    df['Population_MAxMAF'] = df[[x for x in df.columns  if "_AF" in x]].idxmax(axis=1)
    df2=df.drop([x for x in df.columns  if "_AF" in x],axis=1)
    output = subprocess.check_output(f"wc -l {file}", shell=True).decode().strip().split()[0]
    print(file,df2.shape[0],output)
    master_df=pd.concat([master_df,df2])

## Calculating VAF and removing panel_of_normals, germline and common_variant
master_df['VAF']=master_df['t_alt_count'] / master_df['t_depth']
master_df2=master_df[~master_df['FILTER'].str.contains("panel_of_normals|germline|common_variant")]
master_df2=master_df2[master_df2["t_depth"]>=10]
#master_df3=master_df2[["Hugo_Symbol","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification",
#                       "Tumor_Sample_Barcode","IMPACT","CLIN_SIG",'VAF','FILTER']]



## Creating normal MAF file by merging all normal sample wise MAF files

normal_master_df=pd.DataFrame()

for file in normal_files:
    df=pd.read_csv(file, sep="\t", quotechar='"',quoting=3,skiprows=1)
    df["MaxMAF"]=df[[x for x in df.columns  if "_AF" in x]].max(axis=1)
    df['Population_MAxMAF'] = df[[x for x in df.columns  if "_AF" in x]].idxmax(axis=1)
    df2=df.drop([x for x in df.columns  if "_AF" in x],axis=1)
    output = subprocess.check_output(f"wc -l {file}", shell=True).decode().strip().split()[0]
    print(file,df2.shape[0],output)
    normal_master_df=pd.concat([normal_master_df,df2])

normal_master_df['VAF']=normal_master_df['n_alt_count'] / normal_master_df['n_depth']
normal_master_df=normal_master_df[normal_master_df["n_depth"]>=10]
normal_master_df=normal_master_df[normal_master_df["VAF"]>0.01]
normal_master_df2=normal_master_df[['Chromosome', 'Start_Position', 'End_Position','Reference_Allele', 'Tumor_Seq_Allele1','Tumor_Seq_Allele2']]
normal_master_df2["Type"]="Normal"

## merging Tumor and Normal MAF files
merged=pd.merge(master_df2,normal_master_df2,on=['Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2'],how="left")
## Remove the variants if it is reeporteed in the normal samples
merged=merged[merged["Type"].isna()].drop("Type",axis=1)
#Remove the variants if it has AF >0.01 in public data bases
merged=merged[ (merged["MaxMAF"]<0.01) |(merged["MaxMAF"].isna()) ]

##Removee the variants if VAF <0.01 or VAF>0.4
merged_2=merged[(merged["VAF"]>0.01) & (merged["VAF"]<0.4) ]

target_strings = ["KRAS", "GNAS", "MUC16", "BRAF", "TTN", "TGFBR2", "SMAD4", "TP53"]
selected=merged_2[merged_2["Hugo_Symbol"].isin(target_strings)]


clinical_df=pd.read_csv("/mnt/disks/sdb/back_VCF-Ap/maf/back_TumorOnly_clinical_details.tsv",sep="\t")

merged_2=pd.merge(merged_2,clinical_df,on="Tumor_Sample_Barcode")
merged_2.to_csv("TumorOnly_noGermline_NoCommonvariant_VAFGT0.01_LT0.4.tsv", sep = "\t",index=None)

selected=pd.merge(selected,clinical_df,on="Tumor_Sample_Barcode")
selected.to_csv("TumorOnly_noGermline_NoCommonvariant_VAFGT0.01_LT0.4_SelectedGenes.tsv", sep = "\t",index=None)



for group in clinical_df["Group"].unique():
      merged_2[merged_2["Group"]==group].to_csv(f"{group}_noGermline_NoCommonvariant_VAFGT0.01_LT0.4.tsv", sep = "\t",index=None)
      selected[selected["Group"]==group].to_csv(f"{group}_noGermline_NoCommonvariant_VAFGT0.01_LT0.4_SelectedGenes.tsv", sep = "\t",index=None)

##Split the clinical file
for group in clinical_df["Group"].unique():
      df=clinical_df[clinical_df["Group"]==group]
      df.to_csv(f"{group}_clinical_details.tsv",sep="\t",index=None)
      


