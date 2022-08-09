# Orrantia, et al. 2022 Repository

This is the repository for differential expression analyzes and functional analyzes of "In vivo expansion of a functionally different CD9 decidual-like NK cell subset following autologous hematopoietic stem cell transplantation"
Orrantia et al. 

You can find the raw RNA-Seq data at GEO under the ID: GSE199608.

The data after the quantification by RSEM can also be found in the Raw_Data/ folder.

To calculate the DEGs you just have to execute the python script in DEG/ maintaining the directory structure. It needs you to unpack the GTF of the human transcriptome (from ENSEMBL) present in Clean_GTF/.

For functional analyses, after calculating the DEGs, it is necessary to execute the R script present in FunctionalAnalysis/.

Please, if you use this data or the information in the paper, cite it.
