#!/usr/bin/python3

#############################
# This script takes a genes.gff3 file (GFF3 file where only genes are present) and creates two set of BED files:
# 1. co_genes.bed file with the genes codirecred with replication and rev_genes file with the genes directed opposite to replication.
# 2. Four BED files with intergenic intervals (that are shorter than 500 bp) plus 100 nt from each side: conv.bed conteins intervals between the two convergent genes, div.bed conteins intervals between the two divergent genes, plus.bed conteins intervals between the two genes in + strand and minus.bed conteins intervals between the two genes in - strand.
#############################


import pandas as pd

# Read the genes.gff3 file.
gff = pd.read_csv('../ref_tmp/genes.gff3', sep='\t', header=None)
gff.columns = ['chr', 'source', 'type', 'start', 'end', 'empty_1', 'strand', 'empty_2', 'description']
gff = gff[['chr', 'start', 'end', 'strand']]
genes = pd.DataFrame({'chr': gff['chr'], 'start': gff['start'], 'end': gff['end'], 'name': gff['start'], 'score': '0', 'strand': gff['strand']})

# Save the genes codirected with the replication in co_genes.bed.
co_genes = genes.query('(start > 1700000 and end < 3815000 and strand == "-") or (start > 3815000 and strand == "+") or (end < 1200000 and strand == "+")')
co_genes.to_csv('../genes_coverage/co_genes.bed', sep="\t", header=False, index=False)

# Save the genes directed opposite to replication in rev_genes.bed.
rev_genes = genes.query('(start > 1700000 and end < 3815000 and strand == "+") or (start > 3815000 and strand == "-") or (end < 1200000 and strand == "-")')
rev_genes.to_csv('../genes_coverage/rev_genes.bed', sep="\t", header=False, index=False)

# Create empty DataFrames.
column_names = ["chr", "start", "end"]
conv = pd.DataFrame(columns = column_names)
div = pd.DataFrame(columns = column_names)
plus = pd.DataFrame(columns = column_names)
minus = pd.DataFrame(columns = column_names)

# Select intergenic intervals between convergent genes. Write them to the conv DataFrame.
for index in range(len(genes)-1):
    if (genes.loc[index, 'strand'] == "+") and (genes.loc[index+1, 'strand'] == "-"):
        new_row = {'chr': genes.loc[index, 'chr'], 'start': genes.loc[index, 'end'], 'end': genes.loc[index+1, 'start']}
        conv = conv.append(new_row, ignore_index=True)

conv = conv.query('(end - start) < 500 and (end - start) > -200')
conv['start'] = conv['start'] - 100
conv['end'] = conv['end'] + 100
conv['name'] = conv['start']

# Select intergenic intervals between divvergent genes. Write them to the div DataFrame.
for index in range(len(genes)-1):
    if (genes.loc[index, 'strand'] == "-") and (genes.loc[index+1, 'strand'] == "+"):
        new_row = {'chr': genes.loc[index, 'chr'], 'start': genes.loc[index, 'end'], 'end': genes.loc[index+1, 'start']}
        div = div.append(new_row, ignore_index=True)

div = div.query('(end - start) < 500 and (end - start) > -200')
div['start'] = div['start'] - 100
div['end'] = div['end'] + 100
div['name'] = div['start']

# Select intergenic intervals between the genes in the + strand. Write them to the plus DataFrame.
for index in range(len(genes)-1):
    if (genes.loc[index, 'strand'] == "+") and (genes.loc[index+1, 'strand'] == "+"):
        new_row = {'chr': genes.loc[index, 'chr'], 'start': genes.loc[index, 'end'], 'end': genes.loc[index+1, 'start']}
        plus = plus.append(new_row, ignore_index=True)
        
plus = plus.query('(end - start) < 500 and (end - start) > -200')
plus['start'] = plus['start'] - 100
plus['end'] = plus['end'] + 100
plus['name'] = plus['start']

# Select intergenic intervals between the genes in the - strand. Write them to the minus DataFrame.
for index in range(len(genes)-1):
    if (genes.loc[index, 'strand'] == "-") and (genes.loc[index+1, 'strand'] == "-"):
        new_row = {'chr': genes.loc[index, 'chr'], 'start': genes.loc[index, 'end'], 'end': genes.loc[index+1, 'start']}
        minus = minus.append(new_row, ignore_index=True)
        
minus = minus.query('(end - start) < 500 and (end - start) > -200')
minus['start'] = minus['start'] - 100
minus['end'] = minus['end'] + 100
minus['name'] = minus['start']

# Save DataFrames to BED files.
conv.to_csv('../intergenic_coverage/conv.bed', sep="\t", header=False, index=False)
div.to_csv('../intergenic_coverage/div.bed', sep="\t", header=False, index=False)
plus.to_csv('../intergenic_coverage/plus.bed', sep="\t", header=False, index=False)
minus.to_csv('../intergenic_coverage/minus.bed', sep="\t", header=False, index=False)
