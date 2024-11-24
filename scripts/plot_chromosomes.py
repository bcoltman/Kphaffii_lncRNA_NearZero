#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

rename_dict = {"cbs7435_chr1":"Chr. 1", "cbs7435_chr2":"Chr. 2", "cbs7435_chr3":"Chr. 3","cbs7435_chr4":"Chr. 4","cbs7435_mt":"Chr. Mit"}

def import_gtf(gtf_file, biotype):
    # Read the GTF file
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
    gtf.columns = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gtf.chr = gtf.chr.apply(lambda x: rename_dict[x])
    #Extracting relevant info from attributes
    gtf['transcript_id'] = gtf['attribute'].str.extract('transcript_id "([^"]+)"')
    gtf['gene_id'] = gtf['attribute'].str.extract('gene_id "([^"]+)"')
    gtf['biotype'] = biotype
    # Filter for only transcripts 
    gtf = gtf[gtf['feature'] == 'transcript']
    return gtf

def import_gaps(gtf_file, biotype):
    # Read the GTF file
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
    gtf.columns = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gtf.chr = gtf.chr.apply(lambda x: rename_dict[x])
    #Extracting relevant info from attributes
    #gtf['transcript_id'] = gtf['attribute'].str.extract('transcript_id "([^"]+)"')
    #gtf['gene_id'] = gtf['attribute'].str.extract('gene_id "([^"]+)"')
    gtf['biotype'] = biotype
    gtf.set_index("chr",inplace=True)
    # Filter for only transcripts 
    #gtf = gtf[gtf['feature'] == 'transcript']
    return gtf



def plot_gtf_transcripts(gtf_df, chrom_size_file, gaps_df, name, outdir):
    # Import chromosome size
    chrom_sizes = pd.read_csv(chrom_size_file, sep='\t', header=None)
    chrom_sizes.columns = ['chr', 'size']
    chrom_sizes.chr = chrom_sizes.chr.apply(lambda x: rename_dict[x])
    chrom_sizes = chrom_sizes.set_index("chr")['size'].to_dict()
    # Plotting
    palette = {'lncRNA':'red', 'coding':'lightblue'}
    fig, ax = plt.subplots(figsize= (12, 10))
    sns.swarmplot(y='chr', x='start', hue='biotype', data=gtf_df, dodge=False, ax=ax, palette=palette)
    # Draw a vertical line for each chromosome end
    for i, ylab in enumerate(ax.get_yticklabels()):
        x, y = ax.collections[i].get_offsets().T
        ymin = y.min()
        ymax = y.max()
        #print(ylab.get_text())
        #if ylab.get_text() in gaps_df.index:
            #chr_gaps = gaps_df.loc[[ylab.get_text()],:]
            #for j, gap in chr_gaps.iterrows():
                #ax.add_patch(plt.Rectangle((gap["start"], ymin), gap["end"] - gap["start"],ymax-ymin, facecolor="green",lw=5, edgecolor="green"))
#	ax.vlines(x=gaps_df.loc[ylab.get_text(), "start"], ymin=ymin, ymax=ymax, colors="black", linestyles="--")
        ax.vlines(x=chrom_sizes[ylab.get_text()], ymin=ymin, ymax=ymax, colors="black", linestyles="--")
        #, row in chrom_sizes.iterrows():
        #if row['chr'] in gtf_df['chr'].values:
          #  ax.axvline(x=row['size'], color='black', linestyle='--')
    #ax.set_ylabel('Chromosome', fontsize=15)
    #ax.set_xlabel('Start Coordinate', fontsize=15)
    ax.set_ylabel('', fontsize=15)
    ax.set_xlabel('', fontsize=15)
    #ax.set_title('Swarm plot of transcripts by chromosome', fontsize=15)
    #ax.ticklabel_format(axis="x", style="scientific")
    ax.tick_params(labelsize=15, labelbottom=False, bottom=False) #, 
    ax.legend(fontsize=20, markerscale=3, frameon=False, bbox_to_anchor=(0.97, 0.97), loc="upper right")
    for i in ['left', 'top', 'right', 'bottom']:
        ax.spines[i].set_visible(False)
    #ax.set_yticklabels(ax.get_xticklabels(),rotation=45)
    plt.tight_layout()
   
    fig.savefig(f"{outdir}/{name}Chromosomes_lncRNA.png", dpi=300)

def plot_gtf_transcripts_vertical(gtf_df, chrom_size_file, gaps_df, name, outdir):
    # Import chromosome size
    chrom_sizes = pd.read_csv(chrom_size_file, sep='\t', header=None)
    chrom_sizes.columns = ['chr', 'size']
    chrom_sizes.chr = chrom_sizes.chr.apply(lambda x: rename_dict[x])
    chrom_sizes = chrom_sizes.set_index("chr")['size'].to_dict()
    # Plotting
    palette = {'lncRNA':'red', 'coding':'lightblue'}
    fig, ax = plt.subplots(figsize= (12, 10))
    sns.swarmplot(x='chr', y='start', hue='biotype', data=gtf_df, dodge=False, ax=ax, palette=palette)
    # Draw a vertical line for each chromosome end
    for i, ylab in enumerate(ax.get_xticklabels()):
        x, y = ax.collections[i].get_offsets().T
        ymin = x.min()
        ymax = x.max()
        #print(ylab.get_text())
        if ylab.get_text() in gaps_df.index:
            chr_gaps = gaps_df.loc[[ylab.get_text()],:]
            #for j, gap in chr_gaps.iterrows():
                #ax.add_patch(plt.Rectangle((ymin, gap["start"]), ymax-ymin, gap["end"] - gap["start"], facecolor="green",lw=5, edgecolor="green"))

        ax.hlines(y=chrom_sizes[ylab.get_text()], xmin=ymin, xmax=ymax, colors="black", linestyles="--")
       
    ax.set_ylabel('', fontsize=15)
    ax.set_xlabel('', fontsize=15)
   
    ax.tick_params(labelsize=15, labelleft=False, left=False) #, 
    ax.legend(fontsize=20, markerscale=3, frameon=False, bbox_to_anchor=(0.97, 0.97), loc="upper right")
    for i in ['left', 'top', 'right', 'bottom']:
        ax.spines[i].set_visible(False)
    #ax.set_yticklabels(ax.get_xticklabels(),rotation=45)
    plt.tight_layout()
   
    fig.savefig(f"{outdir}/{name}Chromosomes_lncRNA.png", dpi=300)


######

coding_gtf, lncrna_gtf, gaps_gtf, chrom_sizes, outdir = sys.argv[1:]

cod_gtf = import_gtf(coding_gtf, "coding")
lncRNA_gtf = import_gtf(lncrna_gtf, "lncRNA")
gaps_gtf = import_gaps(gaps_gtf, "Gaps")	
combined = pd.concat([cod_gtf, lncRNA_gtf])

#no_mt = combined[combined.chr != "cbs7435_mt"]
no_mt = combined[combined.chr != "Chr. Mit"]
plot_gtf_transcripts(no_mt, chrom_sizes, gaps_gtf, "NoMT", outdir)
plot_gtf_transcripts(combined, chrom_sizes, gaps_gtf, "All", outdir)
plot_gtf_transcripts_vertical(no_mt, chrom_sizes, gaps_gtf, "Vertical_NoMT", outdir)
plot_gtf_transcripts_vertical(combined, chrom_sizes, gaps_gtf, "Vertical_All", outdir)    


