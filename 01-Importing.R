# If not intalled Biocounductor, install it with the following command
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install recount
BiocManager::install("recount")

# Load recount3 package
library(recount3)

# List available projects
human_projects <- available_projects()

# Filter SRP150320 project
proj_info <- subset(
  human_projects,
  project == "SRP150320" & project_type == "data_sources"
)

# Create a RangedSummarizedExperiment object for the SRP150320 project
rse_gene_SRP150320 <- create_rse(proj_info)

# Visualize and explore data
rse_gene_SRP150320

metadata(rse_gene_SRP150320)

## Number of genes by number of samples
dim(rse_gene_SRP150320)

rowData(rse_gene_SRP150320)


# Transform the recount3 data from counts per nucleotide into counts per lecture
assay(rse_gene_SRP150320, "counts") <- compute_read_counts(rse_gene_SRP150320)


rse_gene_SRP150320 <- expand_sra_attributes(rse_gene_SRP150320)

# Maintain only the columns that start with sra_attribute
colData(rse_gene_SRP150320)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP150320)))
]

# Procesar los datos, convertir las columnas de sra_attribute en factores
rse_gene_SRP150320$sra_attribute.source_name <- factor(rse_gene_SRP150320$sra_attribute.source_name)
rse_gene_SRP150320$sra_attribute.sirna <- factor(rse_gene_SRP150320$sra_attribute.sirna)
rse_gene_SRP150320$sra_attribute.replicate <- factor(rse_gene_SRP150320$sra_attribute.replicate)
rse_gene_SRP150320$sra_attribute.cell_line <- factor(rse_gene_SRP150320$sra_attribute.cell_line)


# Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP150320)[
  ,
  grepl("^sra_attribute.[source_name|sirna|replicate|cell_line]", colnames(colData(rse_gene_SRP150320)))
]))
# Debería salir como 'factor', pero aparece caracter con el summary


#
rse_gene_SRP150320$assigned_gene_prop <- rse_gene_SRP150320$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP150320$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP150320$assigned_gene_prop)

with(colData(rse_gene_SRP150320), tapply(assigned_gene_prop, sra_attribute.sirna, summary))

rse_gene_SRP150320_unfiltered <- rse_gene_SRP150320

hist(rse_gene_SRP150320$assigned_gene_prop)

table(rse_gene_SRP150320$assigned_gene_prop < 0.3)

# Calcular los nivelos medios de expresión con las cuentas
gene_means <- rowMeans(assay(rse_gene_SRP150320, "counts"))

summary(gene_means)

rse_gene_SRP150320 <- rse_gene_SRP150320[gene_means > 0.1, ]
summary(rse_gene_SRP150320)
dim(rse_gene_SRP150320)


round(nrow(rse_gene_SRP150320) / nrow(rse_gene_SRP150320_unfiltered) * 100, 2)

# Normalización de los Datos
library("edgeR")
dge <- DGEList(
  counts = assay(rse_gene_SRP150320, "counts"),
  genes = rowData(rse_gene_SRP150320)
)
dge <- calcNormFactors(dge)


library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP150320)), aes(y = assigned_gene_prop, x = sra_attribute.sirna)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("siRNA Group")
)

mod <- model.matrix(~ sra_attribute.sirna + assigned_gene_prop,
                    data = colData(rse_gene_SRP150320)
)
colnames(mod)


# Análisis estadístico con limma
library("limma")
vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 1,
  number = nrow(rse_gene_SRP150320),
  sort.by = "none"
)
dim(de_results)
head(de_results)

## Genes diferencialmente expresados entre pre y post natal con FDR < 5%
table(de_results$adj.P.Val < 0.05)

## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 2)

volcanoplot(eb_results, coef = 1, highlight = 5, names = de_results$gene_name)
de_results[de_results$gene_name %in% c("ZSCAN2", "VASH2", "KIAA0922"), ]


# 50 genes mayormente expresados

## Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene_SRP150320)[, c("sra_attribute.sirna","sra_attribute.replicate")])
colnames(df) <- c("siRNA", "Replicate")

## Hagamos un heatmap
library("pheatmap")
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)


## Para colores
library("RColorBrewer")

## Conviertiendo los grupos de edad a colores
col.group <- df$siRNA
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

## MDS por grupos de edad
plotMDS(vGene$E, labels = df$siRNA, col = col.group)

## Conviertiendo los valores de Sex a colores
col.rep <- df$Replicate
levels(col.rep) <- brewer.pal(nlevels(col.rep), "Dark2")
col.rep <- as.character(col.rep)

## MDS por sexo
plotMDS(vGene$E, labels = df$Replicate, col = col.rep)

