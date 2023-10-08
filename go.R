# Install requier libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("gProfileR")

library(gprofiler2)

#Set working directory to file with gene file lists
## Change this PATH !!!!!!!
path = "......................"
setwd(path)

# Read files
genes_background <- read.csv(file="background.txt", header=TRUE, sep=",")
genes_up <- read.csv(file="up.txt", header=TRUE, sep=",")
genes_down <- read.csv(file="down.txt", header=TRUE, sep=",")

# Store genes as a vector
genes_background_v = as.vector(genes_background[,1])
genes_up_v = as.vector(genes_up[,1])
genes_down_y = as.vector(genes_down[,1])

# Run gene ontology
# up
goterms_up <- gost(query = genes_up_v, 
                   custom_bg = genes_background_v, 
                   organism = "hsapiens", 
                   significant = TRUE, 
                   exclude_iea = TRUE, 
                   user_threshold = 0.05, 
                   correction_method = "g_SCS")

# down
goterms_down <- gost(query = genes_down_v, 
                   custom_bg = genes_background_v, 
                   organism = "hsapiens", 
                   significant = TRUE, 
                   exclude_iea = TRUE, 
                   user_threshold = 0.05, 
                   correction_method = "g_SCS")

# Export results
write.csv(goterms_up, file ='gene_ontology_up.csv')
write.csv(goterms_down, file = 'gene_ontology_down.csv')

# Visualization
gostplot(goterms_up, capped = TRUE, interactive = TRUE)
gostplot(goterms_down, capped = TRUE, interactive = TRUE)
