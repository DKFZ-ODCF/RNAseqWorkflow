args <- commandArgs(trailingOnly = TRUE)

full_matrix <- read.table(args[1], sep="\t", header=T, row.name=1, comment.char="")

# gene_names_matrix <- full_matrix[,1]
# gene_names_matrix <- as.data.frame(gene_names_matrix)

# rownames(gene_names_matrix) <- rownames(full_matrix)

expression_matrix <- full_matrix[,2:(length(full_matrix[1,])-1)]

source(args[2])
qc.dir <- args[3]
makeQC_plot(expression_matrix, args[4], qc.dir)