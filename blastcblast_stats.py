# -*- coding: utf-8 -*-
"""
Original author: Ahmed Elsherbini
Date: 23-10-2025
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

# ----------------------- ARGUMENT PARSING ----------------------- #
my_parser = argparse.ArgumentParser(description='NCBI Species Tree Builder with Duplicates Filtering')
my_parser.add_argument('-i', '--input', metavar='input', type=str, required=False, help='Path to input CSV file')
my_parser.add_argument('-og', '--outgroup', metavar='outgroup', type=str, required=False, help='Outgroup species (optional)')
my_parser.add_argument('-m', '--mode', metavar='mode', type=str, choices=['strict', 'relaxed'], default='relaxed', help='Mode to handle duplicates: strict or relaxed')
args = my_parser.parse_args()

# ----------------------- DEV MODE (Optional for Debugging) ----------------------- #
# Uncomment the following lines to test without command-line arguments
# args.input = "2029_30_mi.csv"
# args.outgroup = "deinococcus_radiodurans"
# args.mode = "strict"

# ----------------------- MAIN VARIABLES ----------------------- #
f_name = args.input
og = args.outgroup
mode = args.mode

#

#f_name = "novel.csv"
#og = "deinococcus_radiodurans"
#mode = "strict"



    

#
try:
    Entrez.email = 'drahmedelsherbini@gmail.com'

    # ----------------------- READ CSV FILE ----------------------- #
    if mode == 'strict':
        try:
            df = pd.read_csv(f_name, header=None, error_bad_lines=False)
            df.drop_duplicates(subset=df.columns[0], keep='first', inplace=True)
            print("We are in the strict mode with old pandas")
        except TypeError:
            df = pd.read_csv(f_name, header=None, on_bad_lines='skip')
            df.drop_duplicates(subset=df.columns[0], keep='first', inplace=True)
            print("We are in the strict mode with the new pandas ")


#no drop of duplicates here 
    elif mode == 'relaxed':
        try:
            df = pd.read_csv(f_name, header=None, error_bad_lines=False)
            print("We are in the relaxed mode with old pandas")


        except TypeError:
            df = pd.read_csv(f_name, header=None, on_bad_lines='skip')
            print("We are in the relaxed mode with new pandas")



    # ----------------------- PROCESS DATA ----------------------- #
    #as a marker for blastN 
    if df[0].iloc[0] == 'Description':
        df.rename(columns=df.iloc[0], inplace=True)
        df['Scientific Name'] = df['Scientific Name'].apply(lambda x: ' '.join(str(x).split()[:2]) if isinstance(x, str) else str(x))
        df = df['Scientific Name'].value_counts().reset_index()
        df.columns = ['Species', 'count']
        df = df[~df['Species'].str.contains('Scientific Name', regex=False)]
        df = df[~df['Species'].str.contains('sp.', regex=False)]
        df = df[~df['Species'].str.contains('nan', regex=False)]
        assmebly = []
        print("Detected NCBI (BlastN) input format.")
    
    #as a marker for blastp results

    elif df[0].iloc[0] == 'Cluster Composition':
       df.rename(columns=df.iloc[0], inplace=True)
       df["Species"] = df["Representative sequence"].str.extract(r"\[([^\]]+)\]")
       df = df['Species'].value_counts().reset_index()
       df.columns = ['Species', 'count']
       df = df[~df['Species'].str.contains('Representative sequence', regex=False)]
       df = df[~df['Species'].str.contains('sp.', regex=False)]
       df = df[~df['Species'].str.contains('nan', regex=False)]
       df = df[~df['Species'].str.contains('bacterium', regex=False)]
       assmebly = []
       print("Detected NCBI (Blastp) input format.")


    else:
        df = df.iloc[:, 0].to_frame()
        df = df[0].str.split().str[:2].str.join(' ').to_frame()
        df = df[0].value_counts().reset_index()
        df.columns = ['Species', 'count']
        df = df[~df['Species'].str.contains('sp.', regex=False)]
        df = df[~df['Species'].str.contains('nan', regex=False)]
        assmebly = []
        print("Detected (Cblaster) input format.")

    print("Running blastCblast_stats 2025 (Ahmed Elsherbini)")

    # ----------------------- NCBI Assembly Query ----------------------- #
    for row in df['Species']:
        species_name = str(row)
        handle = Entrez.esearch(db="assembly", term=species_name, retmode="xml")
        record = Entrez.read(handle)
        count = int(record['Count'])
        print(f"Number of {species_name} occurrences in NCBI assembly database: {count}")
        assmebly.append(count)

    df['assembly'] = assmebly
    df['%_in_assembly_db'] = df['count'] / df['assembly'] * 100
    file_name = f'database_percentage_{f_name[:-4]}.csv'
    df.to_csv(file_name, index=False)
    print("Saved database statistics to:", file_name)

    # ----------------------- TAXONOMY TREE ----------------------- #
    ncbi = NCBITaxa()
    species_taxids = {}
    species_names = pd.Series(df["Species"])
    species_names = species_names[~species_names.str.contains("sp.", regex=False)]

    if og:
        og = og.replace("_", " ")
        species_names.loc[len(species_names)] = og

    for species_name in species_names:
        taxid = ncbi.get_name_translator([species_name])
        species_taxids[species_name] = taxid[species_name][0] if taxid else None

    taxa_ids = [taxid for taxid in species_taxids.values() if taxid]
    tree = ncbi.get_topology(taxa_ids)

    def annotate_tree_with_scientific_names(tree):
        for node in tree.traverse():
            if node.is_leaf():
                taxid = int(node.name)
                scientific_name = ncbi.get_taxid_translator([taxid]).get(taxid)
                node.name = scientific_name if scientific_name else "Unknown"

    annotate_tree_with_scientific_names(tree)
    print(tree)

    output_file = f"{f_name[:-4]}_tree.nwk"
    tree.write(outfile=output_file)

    # ----------------------- PIE CHART TREE ----------------------- #
    df1 = df[['Species', '%_in_assembly_db']]
    pie_data = df1.set_index('Species')['%_in_assembly_db'].to_dict()
    pie_data = {k: v for k, v in pie_data.items() if "sp." not in k and "Species" not in k and not math.isinf(v)}

    def layout(node):
        if node.is_leaf() and node.name in pie_data:
            pie_values = [pie_data[node.name], 100 - pie_data[node.name]]
            colors = ["green", "lightgray"]
            pie_face = PieChartFace(pie_values, colors=colors, width=50, height=50)
            faces.add_face_to_node(pie_face, node, column=1)

    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = True

    tree.render(f"{f_name[:-4]}_tree_with_pies.pdf", tree_style=ts)
    print("Tree with pie charts rendered successfully!")

except Exception as e:
    print("An error occurred:", e)
