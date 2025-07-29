import os
import tarfile
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def get_coexp_network(dna_coordinate, iMARGI_files, freq):
  df = pd.read_csv(iMARGI_files, sep='\t', comment='#', header=None)

  # truncate for 10 columns
  if df.shape[1] > 10:
    df = df.iloc[:, :10]


  df.columns = [
    "RNA_chr", "RNA_start", "RNA_end",
    "DNA_chr", "DNA_start", "DNA_end",
    "name", "score", "RNA_strand", "DNA_strand"
  ]


  linc_chr = dna_coordinate.split(':')[0]
  coord = dna_coordinate.split(':')[1].replace(',', '').replace('–', '-')
  linc_start, linc_end = map(int, coord.split('-'))


  linc_df = df[
    (df['DNA_chr'] == linc_chr) &
    (df['DNA_start'] <= linc_end) &
    (df['DNA_end'] >= linc_start)
  ]

  vc = linc_df.value_counts().reset_index(name='count')
  vc = vc.sort_values(by='count', ascending=False).head(100)

  G = nx.Graph()
  G.add_node(dna_coordinate)

  for _, row in vc.iterrows():
    dna_label = f"{row['DNA_chr']}:{row['DNA_start']}-{row['DNA_end']}"
    G.add_node(dna_label, count=row['count'])

    G.add_edge(dna_coordinate, dna_label, weight=row['count'])

  plt.figure(figsize=(14, 12))
  pos = nx.spring_layout(G, k=0.6, seed=42)

  node_sizes = []
  node_colors = []

  for node in G.nodes():
      if node == dna_coordinate:
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

  plt.title(f"{dna_coordinate} Interaction Network (Highlighting count > {freq})", fontsize=16)
  plt.axis('off')
  plt.show()
