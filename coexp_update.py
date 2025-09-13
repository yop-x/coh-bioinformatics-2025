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


def build_norm_one_file(df):
    # aggregate once
    edges = (df.groupby(["RNA_gene_name","DNA_gene_name"])
               .size().reset_index(name="count"))
    N = edges["count"].sum()
    edges["EPM"] = 1e6 * edges["count"] / N                 # library-size scaling
    edges["p_pair_given_RNA"] = edges["count"] / edges.groupby("RNA_gene_name")["count"].transform("sum")
    return edges, N




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
        
        df = df.rename(columns={
            'RNA_chr':'Chromosome',
            'RNA_start':'Start',
            'RNA_end':'End'
        })
        
        df = df[['Chromosome','Start','End', 'DNA_chr', 'DNA_start', 'DNA_end']]
        
        
        # left join RNA names with hg38 reference 
        
        
        pr_df = pr.PyRanges(df)
        pr_genes = pr.PyRanges(genes)
        
        joined = pr_df.join(pr_genes, how='left')
        jdf = joined.df
        
        contained = jdf[(jdf["Start"] >= jdf["Start_b"]) & (jdf["End"] <= jdf["End_b"])]
        contained["container_len"] = contained["End_b"] - contained["Start_b"]
        
        best = (contained
        .sort_values(["Chromosome","Start","End","container_len"])
        .drop_duplicates(subset=["Chromosome","Start","End"], keep = 'first'))
        
        out = df.merge(best[["Chromosome","Start","End","gene_name"]],
               on=["Chromosome","Start","End"], how="left")
          
          
          
          
        
          
          
        # doing the same for DNA names 
        out = out.rename(columns={
            'Chromosome' : 'RNA_chr',
            'Start' : 'RNA_start',
            'End' : 'RNA_end',
            'gene_name' : 'RNA_gene_name'
        })
        
        
        out = out.rename(columns={
            'DNA_chr' : 'Chromosome',
            'DNA_start' : 'Start',
            'DNA_end' : 'End',
        })
        
        out_dna = out[['Chromosome', 'Start', 'End']]
        
        pr_dna = pr.PyRanges(out_dna)
        
        joined = pr_dna.join(pr_genes, how='left')
        
        jdf = joined.df
        contained = contained = jdf[(jdf['Start']>= jdf['Start_b']) & (jdf['End']<= jdf['End_b'])]
        contained['container_len'] = contained['End_b'] - contained['Start_b']
        
        best = (contained
        .sort_values(['Chromosome', 'Start', 'End', 'container_len'])
        .drop_duplicates(subset=['Chromosome', 'Start', 'End'], keep='first'))
        
        
        
        out2 = out2 = out.merge(best[['Chromosome', 'Start', 'End', 'gene_name']],
               on=['Chromosome', 'Start', 'End'], how='left') 
        
        out2 = out2.rename(columns={
            'gene_name' : 'DNA_gene_name'
        })
        
        df = out2
        
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
            
            RNA_gene_name = query_df['RNA_gene_name'].value_counts().index[0]
            
        # when query is RNA_gene_name    
        else:
            if query not in df['RNA_gene_name'].values:
                print(f'Sorry, gene name {query} not found in dataset.')
                
            else:
                RNA_gene_name = query 
        
        query_df = build_norm_one_file(df)[0]
        query_df = query_df[query_df['RNA_gene_name'] == RNA_gene_name]
        
        
        
        vc = (query_df[["DNA_gene_name", "EPM"]]
                .groupby("DNA_gene_name", as_index=False)["EPM"].sum()
                .rename(columns={"EPM": "count"})
                .sort_values(by="count", ascending=False)
                .head(200))
    
        G = nx.Graph()
        G.add_node(query)
        
        
        for _, row in vc.iterrows():
            dna_label = f"{row['DNA_gene_name']}"
            G.add_node(dna_label, count=row['count'])
            G.add_edge(query, dna_label, weight=row['count'])

        plt.figure(figsize=(14, 12))
        pos = nx.spring_layout(G, k=0.6, seed=42)

        node_size_query = 500 
        node_size_other = 120 
        
        node_sizes = [node_size_query if n == query else node_size_other for n in G.nodes()]
        node_colors = [
            'red' if n == query
            else ('orange' if G.nodes[n].get('count', 1)> freq else 'skyblue')
            for n in G.nodes()
        ]

        

        nx.draw(
            G, pos,
            with_labels=True,
            node_size=node_sizes,
            node_color=node_colors,
            edge_color='gray',
            font_size=8
        )

        plt.title(f"{query} Interaction Network (Highlighting EPM > {freq})", fontsize=16)
        plt.axis('off')
        plt.tight_layout()
        plt.show()

        
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
            # at least 50% matching 
            if match_ratio > 0.5:
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
            
            RNA_gene_name = query_df['RNA_gene_name'].value_counts().index[0]
            
            
        # if the query is RNA_gene_name   
        else:
            if query not in df['RNA_gene_name'].values:
                print(f'Sorry, gene name {query} not found in dataset.')
                
            else:
                RNA_gene_name = query
                
        query_df = build_norm_one_file(df)[0]
        query_df = query_df[query_df['RNA_gene_name'] == RNA_gene_name]
        
                
                
        vc = (query_df[["DNA_gene_name", "EPM"]]
                .groupby("DNA_gene_name", as_index=False)["EPM"].sum()
                .rename(columns={"EPM": "count"})
                .sort_values(by="count", ascending=False)
                .head(200))
    
        G = nx.Graph()
        G.add_node(query)
        
        
        for _, row in vc.iterrows():
            dna_label = f"{row['DNA_gene_name']}"
            G.add_node(dna_label, count=row['count'])
            G.add_edge(query, dna_label, weight=row['count'])

        plt.figure(figsize=(14, 12))
        pos = nx.spring_layout(G, k=0.6, seed=42)

        node_size_query = 500 
        node_size_other = 120 
        
        node_sizes = [node_size_query if n == query else node_size_other for n in G.nodes()]
        node_colors = [
            'red' if n == query
            else ('orange' if G.nodes[n].get('count', 1)> freq else 'skyblue')
            for n in G.nodes()
        ]

        

        nx.draw(
            G, pos,
            with_labels=True,
            node_size=node_sizes,
            node_color=node_colors,
            edge_color='gray',
            font_size=8
        )

        plt.title(f"{query} Interaction Network (Highlighting EPM > {freq})", fontsize=16)
        plt.axis('off')
        plt.tight_layout()
        plt.show()