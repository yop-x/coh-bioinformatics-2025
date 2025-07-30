import os
import tarfile
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import re
import pyranges as pr 


# the updated function should be able 
# to take EITHER dna_coordinate or 
# gene name as input 'query' to generate 
# the corresponding co-expression network 

def get_coexp_network_updated(query, iMARGI_files, freq):
    df = pd.read_csv(iMARGI_files, sep='\t', comment='#', header=None)
    gene_names = pd.read_csv("./genes_df.csv")
    genes = gene_names[['Chromosome','Start','End','gene_name']]

    # if feeding the function with all 10 columns 
    if df.shape[1] == 10:
        
        # assign column names 
        df.columns = [
            "RNA_chr", "RNA_start", "RNA_end",
            "DNA_chr", "DNA_start", "DNA_end",
            "name", "score", "RNA_strand", "DNA_strand"
        ]
        
        
        # left join RNA names with hg38 reference 
        rna_annot = genes.rename(columns={
        'Chromosome':'RNA_chr',
        'Start'     :'RNA_start',
        'End'       :'RNA_end',
        'gene_name' :'RNA_gene_name'})
    
        df = df.merge(
            rna_annot,
            how = 'left',
            on = ['RNA_chr', 'RNA_start', 'RNA_end']
        )
        
        # doing the same for DNA names 
        dna_annot = genes.rename(columns={
            'Chromosome':'DNA_chr',
            'Start'     :'DNA_start',
            'End'       :'DNA_end',
            'gene_name' :'DNA_gene_name'})
        
        df = df.merge(
            dna_annot,
            how = 'left',
            on = ['DNA_chr', 'DNA_start', 'DNA_end']
        )
        
        
        # when query is DNA_coordinate 
        if query[:3]== 'chr':
            
            query_type = 'coordinate'
            query_chr = query.split(':')[0]
            coord = query.split(':')[1].replace(',','').replace('–', '-')
            query_start, query_end = map(int, coord.split('-'))
            
            # getting the df with all the DNA interactions
            # that match the query RNA coordinate 
            query_df = df[(df['RNA_chr'] == query_chr) &
            (df['RNA_start'] <= query_end) &
            (df['RNA_end'] >= query_start)
            ]
            
            # add 'DNA_coord'
            query_df['DNA_coord'] = (
                query_df['DNA_chr'].astype(str)
                + ':'
                + query_df['DNA_start'].astype(str)
                + '-'
                + query_df['DNA_end'].astype(str)
            )
            
            # due to lack of time/effort
            # for now only checking the 
            # occurances of different 'DNA_start' to showcase vc 
        
            vc = query_df[['RNA_chr', 'DNA_coord']].value_counts().reset_index(name='count')
            vc = vc.sort_values(by='count', ascending=False).head(100)
            
            G = nx.Graph()
            G.add_node(query)

            for _, row in vc.iterrows():
                dna_label = row['DNA_coord']
                G.add_node(dna_label, count=row['count'])
                G.add_edge(query, dna_label, weight=row['count'])

            plt.figure(figsize=(14, 12))
            pos = nx.spring_layout(G, k=0.6, seed=42)

            node_sizes = []
            node_colors = []

            for node in G.nodes():
                if node == query:
                    node_sizes.append(800)
                    node_colors.append('red')
                else:
                    count = G.nodes[node].get("count", 1)
                    node_sizes.append(100 + count * 10)
                    # Highlight if count > frequency
                    if count > freq:
                        node_colors.append('orange')  # Highlight color
                    else:
                        node_colors.append('skyblue')  # Default
                
        # when query is RNA_gene_name    
        else:
            if query not in df['RNA_gene_name'].values:
                print(f'Sorry, gene name {query} not found in dataset.')
                
            else:
                query_df = df[df['RNA_gene_name'] == query][['RNA_gene_name', 'DNA_gene_name']]
                
                
                vc = query_df.value_counts().reset_index(name='count')
                vc = vc.sort_values(by='count', ascending=False).head(100)
                
                G = nx.Graph()
                G.add_node(query)

                for _, row in vc.iterrows():
                    dna_label = row['DNA_gene_name']
                    G.add_node(dna_label, count=row['count'])
                    G.add_edge(query, dna_label, weight=row['count'])

                plt.figure(figsize=(14, 12))
                pos = nx.spring_layout(G, k=0.6, seed=42)

                node_sizes = []
                node_colors = []

                for node in G.nodes():
                    if node == query:
                        node_sizes.append(800)
                        node_colors.append('red')
                    else:
                        count = G.nodes[node].get("count", 1)
                        node_sizes.append(100 + count * 10)
                        # Highlight if count > frequency
                        if count > freq:
                            node_colors.append('orange')  # Highlight color
                        else:
                            node_colors.append('skyblue')  
            

        # Draw the graph
        plt.figure(figsize=(14, 12))
        pos = nx.spring_layout(G, k=0.6, seed=42)

        nx.draw(
            G, pos,
            with_labels=True,
            node_size=node_sizes,
            node_color=node_colors,
            edge_color='gray',
            font_size=8
        )

        plt.title(f"{query} Interaction Network (Highlighting count > {freq})", fontsize=16)
        plt.axis('off')
        
        return plt.gcf() 
        
    # if feeding the dataset with gene info already 
    else:
        
        # iMARGI datasets with more than 10 cols 
        # usually contains the additional info 
        # this is to locate the targeted columns
        # if the datasets contain more/less than 10 columns    
        
        
        gene_info_cols = []
        # save the column names that contain 
        # addtional DNA/RNA info into a list  
        
        for col in df.columns:
            sample_values = df[col].dropna().astype(str).head(20)
            
            match_ratio = sum(
                bool(re.match(r'^(?:[^|]*\|){2,}[^|]*$', val)) for val in sample_values
            ) / len(sample_values)
            
            # setting the matchin pattern to be
            # at least 70% matching 
            if match_ratio > 0.7:
                gene_info_cols.append(col)
            
            
        # rename the columns 
        df.rename(
            columns=dict(zip(gene_info_cols, ["RNA_gene_info", "DNA_gene_info"])),
            inplace=True
            )
        
        
        # the first 6 columns usually 
        # follow this pattern 
        df.columns = [
            "RNA_chr", "RNA_start", "RNA_end",
            "DNA_chr", "DNA_start", "DNA_end",
            *df.columns[6:].astype(str)
            ]
    
        # create a new df with the preferred columns 
        df = df[["RNA_chr", "RNA_start", "RNA_end",
            "DNA_chr", "DNA_start", "DNA_end",
            "RNA_gene_info", "DNA_gene_info"]]
        
        # filter rows where both gene info fields have exactly 2 pipe characters (i.e., 3 parts)
        
        df_rna_filtered = df[df['RNA_gene_info'].str.count(r'\|')==2]
        
        # now do the same for DNA 
        df_all_filtered = df_rna_filtered[df_rna_filtered['DNA_gene_info'].str.count(r'\|') == 2]
        
        df = df_all_filtered 
        
        # Split RNA_gene_info into 3 new columns
        df[['RNA_gene_id', 'RNA_gene_name', 'RNA_gene_type']] = df['RNA_gene_info'].str.split('|', expand=True)

        # Split DNA_gene_info into 3 new columns
        df[['DNA_gene_id', 'DNA_gene_name', 'DNA_gene_type']] = df['DNA_gene_info'].str.split('|', expand=True)


        # if the query matches dna_coordinates 
        if query[:3]== 'chr':
            query_type = 'coordinate'
            query_chr = query.split(':')[0]
            coord = query.split(':')[1].replace(',','').replace('–', '-')
            query_start, query_end = map(int, coord.split('-'))
        
            # getting the df with all the DNA interactions
            # that match the query RNA coordinate 
            query_df = df[(df['RNA_chr'] == query_chr) &
            (df['RNA_start'] <= query_end) &
            (df['RNA_end'] >= query_start)
            ]
            
            vc = query_df['DNA_gene_name'].value_counts().reset_index(name='count')
            vc = vc.sort_values(by='count', ascending=False).head(100)
        
            G = nx.Graph()
            G.add_node(query)

            for _, row in vc.iterrows():
                dna_label = f"{row['DNA_gene_name']}"
                G.add_node(dna_label, count=row['count'])

                G.add_edge(query, dna_label, weight=row['count'])

            plt.figure(figsize=(14, 12))
            pos = nx.spring_layout(G, k=0.6, seed=42)

            node_sizes = []
            node_colors = []

            for node in G.nodes():
                if node == query:
                    node_sizes.append(800)
                    node_colors.append('red')
                else:
                    count = G.nodes[node].get("count", 1)
                    node_sizes.append(100 + count * 10)
                    # Highlight if count > frequency
                    if count > freq:
                        node_colors.append('orange')  # Highlight color
                    else:
                        node_colors.append('skyblue')  # Default


        # if the query is RNA_gene_name   
        else:
            if query not in df['RNA_gene_name'].values:
                print(f'Sorry, gene name {query} not found in dataset.')
                
            else:
                
                query_df = df[df['RNA_gene_name'] == query][['DNA_gene_name']]
                
                
                vc = query_df.value_counts().reset_index(name='count')
                vc = vc.sort_values(by='count', ascending=False).head(100)
                
                G = nx.Graph()
                G.add_node(query)

                for _, row in vc.iterrows():
                    dna_label = row['DNA_gene_name']
                    G.add_node(dna_label, count=row['count'])
                    G.add_edge(query, dna_label, weight=row['count'])

                plt.figure(figsize=(14, 12))
                pos = nx.spring_layout(G, k=0.6, seed=42)

                node_sizes = []
                node_colors = []

                for node in G.nodes():
                    if node == query:
                        node_sizes.append(800)
                        node_colors.append('red')
                    else:
                        count = G.nodes[node].get("count", 1)
                        node_sizes.append(100 + count * 10)
                        # Highlight if count > frequency
                        if count > freq:
                            node_colors.append('orange')  # Highlight color
                        else:
                            node_colors.append('skyblue')  
                

            
        
        
        # Draw the graph
        plt.figure(figsize=(14, 12))
        pos = nx.spring_layout(G, k=0.6, seed=42)

        nx.draw(
            G, pos,
            with_labels=True,
            node_size=node_sizes,
            node_color=node_colors,
            edge_color='gray',
            font_size=8
        )

        plt.title(f"{query} Interaction Network (Highlighting count > {freq})", fontsize=16)
        plt.axis('off')
        
        return plt.gcf()
        
        
        