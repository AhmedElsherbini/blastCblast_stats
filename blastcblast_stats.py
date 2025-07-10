# -*- coding: utf-8 -*-
"""
Original author : Ahmed Elsherbini 
Date 20-11-2024
Spyder Editor
"""
###########################################################

import pandas as pd
from Bio import Entrez
import ete3
from ete3 import NCBITaxa, TreeStyle, PieChartFace, faces
import argparse
import warnings
import math
###########################################################
warnings.filterwarnings("ignore")

my_parser = argparse.ArgumentParser(description='Hello!')
my_parser.add_argument('-i', '--input', action='store', metavar='input', type=str, help='the path to your file')
my_parser.add_argument('-og', '--outgroup', action='store', metavar='outgroup', type=str, help='the name to your outgroup species')

args = my_parser.parse_args()
f_name = args.input
og = args.outgroup

f_name = "2029_30_mi.csv"
og = "deinococcus_radiodurans"
try:
    Entrez.email = 'drahmedelsherbini@gmail.com'
    df = pd.read_csv(f_name, header=None,error_bad_lines=False)
    #Take care of error_bad_line and numpy update
    
    if df[0].iloc[0] == 'Description':
        df.rename(columns=df.iloc[0], inplace = True)
        df['Scientific Name'] = df['Scientific Name'].apply(lambda x: ' '.join(str(x).split()[:2]) if isinstance(x, str) else str(x))    
        df = df['Scientific Name'].value_counts().reset_index()
        df.columns = ['Species', 'count']
        substring = 'Scientific Name'
        filterx = df['Species'].str.contains(substring)
        df = df[~filterx]
        substring = 'sp.'
        filterx = df['Species'].str.contains(substring)
        substring = 'nan'
        filterx = df['Species'].str.contains(substring)
        df = df[~filterx]
        
        assmebly = []
        print("I think this is a NCBI blast input")
    
    else:
        df = df.iloc[:, 0].to_frame()
        df = df[0].str.split().str[:2].str.join(' ').to_frame()
        df = df[0].value_counts().to_frame()
        df = df.reset_index()
        df.columns = ['Species', 'count']
        substring = 'sp.'
        filterx = df['Species'].str.contains(substring)
        substring = 'nan'
        filterx = df['Species'].str.contains(substring)
        df = df[~filterx]
        assmebly = []
        print("I think this is a Cblaster input")
    
    
    print("And this is blastCblast_stats (Ahmed Elsherbini)")
    
    for row in df['Species']:
        species_name = str(row)
        handle = Entrez.esearch(db="assembly", term=species_name, retmode="xml")
        record = Entrez.read(handle)
        count = int(record['Count'])
        print(f"Number of {species_name} occurrences in NCBI assembly database: {count}")
        assmebly.append(count)
    
    df['assembly'] = assmebly
    df['%_in_assembly_db'] = df['count'] / df['assembly'] * 100
    file_name = 'database_percentage_%s.csv' % (f_name[:-4])
    df.to_csv(file_name, index=False)
    print("Done for the database file")
    
    # Taxonomy retrieval and tree generation
    ncbi = NCBITaxa()
    species_taxids = {}
    species_names = pd.Series(df["Species"])
    species_names = species_names[~species_names.str.contains("sp.", regex=False)]
    
    # Add outlier if provided
    if og:
        og = og.replace("_", " ")
        species_names.loc[len(species_names)] = og
    
    for species_name in species_names:
        taxid = ncbi.get_name_translator([species_name])
        species_taxids[species_name] = taxid[species_name][0] if taxid else None
    
    taxa_ids = [taxid for taxid in species_taxids.values() if taxid]
    tree = ncbi.get_topology(taxa_ids)
    #qc,this to help me, shall print numbers
    #print(tree)
    
    def annotate_tree_with_scientific_names(tree):
        for node in tree.traverse():
            if node.is_leaf():
                taxid = int(node.name)
                scientific_name = ncbi.get_taxid_translator([taxid]).get(taxid)
                node.name = scientific_name if scientific_name else "Unknown"
    
    annotate_tree_with_scientific_names(tree)
    print(tree)
 
    output_file = "%s_tree.nwk" % (f_name[:-4])
    tree.write(outfile=output_file)
    
    # Prepare data for pie charts
    df1 = df[['Species', '%_in_assembly_db']]
    pie_data = df1.set_index('Species')['%_in_assembly_db'].to_dict()
    pie_data = {k: v for k, v in pie_data.items() if "sp." not in k}
    pie_data = {k: v for k, v in pie_data.items() if "Species" not in k}
    pie_data = {k: v for k, v in pie_data.items() if not math.isinf(v)}
    
    def layout(node):
        if node.is_leaf() and node.name in pie_data:
            pie_values = [pie_data[node.name], 100 - pie_data[node.name]]
            colors = ["green", "lightgray"]
            pie_face = PieChartFace(pie_values, colors=colors, width=50, height=50)
            faces.add_face_to_node(pie_face, node, column=1)
    
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = True
    #ts.title.add_face(faces.TextFace("Phylogenetic Tree with Pie Charts", fsize=12), column=0)

    tree.render("%s_tree_with_pies.pdf" % (f_name[:-4]), tree_style=ts)
    print("Tree with pie charts rendered successfully!")

except Exception as e:
   print("An error occurred:", e)
