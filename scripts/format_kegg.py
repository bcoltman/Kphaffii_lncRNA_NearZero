#!/usr/bin/env python
# Format kegg file for gprofiler

import pandas as pd
import json

todf = []

cbs7435_gs155 = pd.read_csv("data/cbs7435_gs115.tsv", sep="\t", header=None, index_col=[1])
kegg = json.load(open('data/kp_kegg.json'))

for main in kegg['children']:
    for broad in main['children']:
        for sub in broad['children']:
            for gn in sub.get('children', [None]):
                if gn:
                    nm = gn['name'].split('[\t;]')[0].split(' ')
                    locus = gn['name'].split(' ')[0]
                    nm = gn['name'].split(' ',maxsplit=1)[1].split('\t')
                    description = nm[0]
                    genename = nm[1].split(';')[0]
                    extra = nm[1].split(';')[1]

                    todf.append([locus, description, genename, extra, sub['name'], broad['name'], main['name']])
                else:
                    # some terms have no children
                    pass
                    #print(sub)



todf = pd.DataFrame(todf)
todf.columns=['GS115Locus', 'Description', 'GeneName', ' Extra', 'Cat1', 'Cat2', 'Cat3']


todf["CBS7435 Locus"] = todf["GS115Locus"].replace(cbs7435_gs155[0].to_dict())

todf.to_csv('reference_genome/GO_analysis/Kegg_Ontology.csv', index=False)

result_df = todf.groupby('Cat1')['CBS7435 Locus'].apply(list).reset_index()
result_df = pd.concat([result_df["Cat1"].str.split(" ", n=1,expand=True), pd.DataFrame(result_df["CBS7435 Locus"].tolist())],axis=1)

# Replace None with empty string
result_df.fillna('', inplace=True)

# Custom function to format rows
def format_row(row):
    # Remove trailing empty strings
    while row and row[-1] == '':
        row.pop()
    return '\t'.join(str(item) for item in row)

# Apply custom function and save to TSV
with open("reference_genome/GO_analysis/Kegg_Ontology.gmt", 'w') as file:
    # Write rows
    for index, row in result_df.iterrows():
        file.write(format_row(row.tolist()) + '\n')



