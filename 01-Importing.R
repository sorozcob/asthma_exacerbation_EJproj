library(recount3)

## SRP150456

recount3::create_rse_manual(
  project = "SRP150456",
  project_home = "data_sources/sra",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)
