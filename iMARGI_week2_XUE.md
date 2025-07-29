# iMARGI Week 2 

1. ## Function overview `get_coexp_network_update`

   	Last week, I was able to draft a preliminary function that takes very limited data (6 columns) from the provided iMARGI datasets (usually range from 10-15 columns) to generate the co-expression network of the given **RNA genomic coordinates** and targeted **iMARGI file** and highlight the high-frequency connections with a customizable `freq` parameter.   
   	This week, I continued working on the function and implemented several algorithms that improve the function's performance. The function is now able to take either **RNA genomic coordinate** or **RNA gene name** as the input `query`. The function will identify which type of query is entered and work accordingly.   
   For **RNA genomic coordinate**, if the iMARGI function only contains the simple RNA/DNA start and end data, the function will then create new columns specifying the gene names from the supplementary file `gene_df.csv` by merging the `gene_df.csv` with existing DNAs and RNAs in our datasets, then match the gene name to the RNA genomic coordinate interval to fetch the data needed for the network graph generation. If the input query is the **RNA gene name**, the function still merges the dataset with the gene name but omits the matching part by directly calling the gene name to fetch the data.   
   The datasets with 15 columns generally work better in this case, thanks to the additional columns already including the gene names we need. Therefore, the function does not go through the step to merge the dataset with `gene_df.csv`. Instead, the function creates three new columns, each for DNA and RNA data, which specify their gene IDs, gene names, and gene types. As a result, the generated network graph provides a much clearer, more intuitive depiction of the RNA–DNA co-expression patterns captured by iMARGI, making downstream analysis significantly easier.  
 


2. Transcript-level bound VS. Gene-level bound 

   Many transcripts appear different despite being the same gene, making value counts in the dataframes challenging. The counts are therefore very much fragmented due to the mismatch of the transcript and its corresponding gene names, causing potential harm to the performance of the function in practice. 

3. ## Missing column/feature names 

   	The column names for the features of the iMARGI dataset are missing from all of the provided datasets. Had to do research and assign them manually. (already implemented in the function to automate the process, might be tricky when inputting datasets with slightly different structure/columns.   
   

   `df.columns = [`  
              `"RNA_chr", "RNA_start", "RNA_end",`  
              `"DNA_chr", "DNA_start", "DNA_end",`  
              `"name", "score", "RNA_strand", "DNA_strand"`  
          `]`  
     
   

   For datasets with 15 columns, they contain columns that congregate all the metadata of the DNA/RNA. An algorithm was implemented to locate the column position in the dataset first by regular expression, then assign the names `RNA_gene_info` and `DNA_gene_info` for further manipulation. 

   

4. ## Extracting useful columns 

For datasets with 10 columns, only `["RNA_chr", "RNA_start", "RNA_end", "DNA_chr", "DNA_start", "DNA_end"]` were taken as useful. 

For datasets with 15 columns, assuming they have the DNA/RNA gene info provided like given iMARGI datasets, three new columns were created for easier RNA/DNA name look-up/matching. (`DNA_gene_id`, `DNA_gene_name`, `DNA_gene_type,RNA_gene_id`, `RNA_gene_name`, `RNA_gene_type)` 

5. ## GENCODE 

   I’m using the `GENCODE` GTF file to map genomic intervals (e.g., `chr2:12343–12345`) back to their gene names, saved in a separate file in the directory as `genes_df.csv`. Right now, my function only returns a hit when the interval exactly matches a transcript’s start and end, so any fragment that falls just outside a specific transcript ends up unannotated. Because gene-level boundaries often extend beyond individual transcripts, this strict exact-match logic can miss perfectly valid gene overlaps.

   Using the `genes_df.csv`, I was still able to annotate DNA and RNA in the dataset despite a great deal of information loss. I implemented that in the function, so that the datasets without gene name information could append this additional information when doing the gene name lookup. 

   

6. ## Challenges 

   1. Gene transcript inconsistency 

      One of the biggest challenges that still needs to be solved is the inconsistency of the genomic region intervals from those of the 10-column datasets. It cannot be easily done by simply parsing them with the `hg38` human gene reference dataset, which I had prepared in the CSV file `genes_df.csv`. One approach is to implement a **fuzzy lookup** to assign all the genes that fall into a certain interval a consistent name.

      

      

   2. Clustering of the genes found in one row.   
        
      While exploring the datasets, especially with the 15-column iMARGI files, I found that the columns `DNA_gene_info` and `RNA_gene_info` sometimes contain more than one set of information about the gene, suggesting that the coexpression might have happened among multiple genes. (even though they are presented in one row) I omitted the clustering coexpression in the function for the simplicity of the analysis, yet it remains an interesting topic for future exploration. 