#!/usr/bin/env python

import sys
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def generate_gc_content_data(input_file, class_code):
	df = []
	
	for seq_record in SeqIO.parse(input_file, "fasta"):
		name = str(seq_record.id)
		seq = str(seq_record.seq)
		seq=seq.upper()
		gc=float(seq.count("G") + seq.count("C")) / float(len(seq))
		result = pd.Series(data=[name, class_code, gc])
		df.append(result)
	df = pd.concat(df,axis=1).T
	df.columns = ["Name", "Type", "GC"]
	return df
	


###################### Plot GC

lncrna_fa, nc_fa, inter_fa, coding_fa, lncrna_list, nc_list, coding_list, stringtie_bed, transcript_tpm, outdir = sys.argv[1:] 

lncrna_gc = generate_gc_content_data(lncrna_fa, "lncRNA")
nc_gc = generate_gc_content_data(nc_fa, "non-coding")
int_gc = generate_gc_content_data(inter_fa,"intergenic")
cod_gc = generate_gc_content_data(coding_fa, "coding")


combined_gc = pd.concat([lncrna_gc, nc_gc, int_gc, cod_gc])
combined_gc = combined_gc[combined_gc["GC"] != 0]
combined_gc["logGC"] = np.log2(combined_gc["GC"].astype(float))

fig, ax = plt.subplots(figsize= (12, 10))
sns.boxplot(y='logGC', x='Type', data=combined_gc, ax=ax)

plt.tight_layout()
fig.savefig(f"{outdir}/GCDistribution.png", dpi=300)

###################### Plot Length


lncrnas = pd.read_csv(lncrna_list,header=None)[0].tolist()
non_coding_genes = pd.read_csv(nc_list,header=None)[0].tolist()
coding_genes = pd.read_csv(coding_list,header=None)[0].tolist()


bed_file = pd.read_csv(stringtie_bed,header=None,sep="\t")

bed_file = bed_file.iloc[:,:6]
bed_file.columns = ["Chr","Start","End","TransID","Score","Strand"]

lncrna_len = bed_file.loc[bed_file["TransID"].isin(lncrnas),:].copy()
lncrna_len.loc[:,"Type"] = "lncRNA"

cod_len = bed_file.loc[bed_file["TransID"].str.contains('|'.join(coding_genes)),:].copy()
cod_len.loc[:,"Type"] = "coding"

ncd_len = bed_file.loc[bed_file["TransID"].str.contains('|'.join(non_coding_genes)),:].copy()
ncd_len.loc[:,"Type"] = "non-coding"

combined_len = pd.concat([cod_len,lncrna_len, ncd_len])
combined_len["Length"]=abs(combined_len["End"]-combined_len["Start"])
combined_len["LogLength"]=np.log2(combined_len["Length"])

fig, ax = plt.subplots(figsize= (12, 10))
sns.boxplot(y='LogLength', x='Type', data=combined_len, ax=ax)

plt.tight_layout()
fig.savefig(f"{outdir}/LengthDistribution.png", dpi=300)

###################### Plot Expression

TPM = pd.read_csv(transcript_tpm,sep="\t")

TPM = TPM.melt("Name")


lncrna_tpm = TPM.loc[TPM["Name"].isin(lncrnas),:].copy()
lncrna_tpm.loc[:,"Type"] = "lncRNA"

cod_tpm = TPM.loc[TPM["Name"].str.contains('|'.join(coding_genes)),:].copy()
#cod_tpm = cod_tpm.loc[~(TPM["Transcript_ID"].str.contains("MSTRG")),:].copy()
cod_tpm.loc[:,"Type"] = "coding"

ncd_tpm = TPM.loc[TPM["Name"].str.contains('|'.join(non_coding_genes)),:].copy()
ncd_tpm.loc[:,"Type"] = "non-coding"


combined_tpm = pd.concat([cod_tpm, lncrna_tpm]) #, ncd_tpm])
combined_tpm["LogExpression"]=np.log2(combined_tpm["value"] +0.01)



fig, ax = plt.subplots(figsize= (12, 5))
sns.boxplot(y='LogExpression', x='Type', data=combined_tpm, ax=ax)

plt.tight_layout()
fig.savefig(f"{outdir}/ExpressionDistribution.png", dpi=300)


###################### Plot All

fig, axes = plt.subplots(ncols=3, figsize=(12,6),sharex=False, gridspec_kw={'width_ratios': [4,3,3]})
all_combined = pd.concat([combined_gc.loc[:,["Type", "logGC"]], combined_len.loc[:, ["Type", "LogLength"]], combined_tpm.loc[:,["Type", "LogExpression"]]])
all_combined = all_combined.melt("Type")

for i, label in enumerate(["logGC", "LogLength", "LogExpression"]):
	if i ==0:
		order = ["coding", "lncRNA", "non-coding","intergenic"]
	elif i == 2:
		order = ["coding", "lncRNA"]
	else:
		order = ["coding", "lncRNA", "non-coding"]
	sns.boxplot(x='Type', y='value', data=all_combined.loc[all_combined["variable"]==label],ax=axes[i], order=order)
	axes[i].set_ylabel(["log2(GC)","log2(Length)", "log2(TPM)"][i],fontsize=15)
	axes[i].set_xlabel("")
	axes[i].tick_params(labelsize=15)
	axes[i].tick_params("x", rotation=45)
	axes[i].spines['right'].set_visible(False)
	axes[i].spines['top'].set_visible(False)


axes[1].set_xlabel("Type",fontsize=15)
plt.tight_layout()
plt.subplots_adjust(wspace=0.2)

fig.savefig(f"{outdir}/CatDistribution.png", dpi=300)


