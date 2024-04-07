library(tidyverse)
library(KEGGREST)
library(stringi)

# loading my functions
source("my-functions.R")

all_cov <- read.table("../combined-outputs/All-combined-KO-function-coverages.tsv", sep = "\t", header = TRUE, quote = "")

# making KEGG summary table with all info for the KOs we have
all_our_KOs <- all_cov %>% pull(KO_ID) %>% as.vector()

Master_KEGG_tab <- make_KO_summary_tab(all_our_KOs)

write.table(Master_KEGG_tab, "kegg-tables/Master-kegg-annotations.tsv", sep="\t", quote=F, row.names=F)

# making N-focused one

  # searching KEGG pathways for Nitrogen
N_pathway_info <- keggFind(database = "pathway", query = "nitrogen")

  # getting "Nitrogen metabolism" pathway ID
N_pathway_ID <- sub(names(N_pathway_info), pattern = "path:", replacement = "")

  # getting info for all KOs associated with "Nitrogen metabolism" pathway
N_pathway_KOs_info <- keggLink("ko", N_pathway_ID)

  # getting KO IDs alone
N_pathway_KOs <- as.character(sub(N_pathway_KOs_info, pattern="ko:", replacement=""))

length(N_pathway_KOs) # 65

Nitrogen_metabolism_KEGG_tab <- make_KO_summary_tab(N_pathway_KOs)

write.table(Nitrogen_metabolism_KEGG_tab, "kegg-tables/Kegg-N-annotations.tsv", sep="\t", quote=F, row.names=F)


# making P-focused one

  # searching KEGG pathways for Nitrogen
P_pathway_info <- keggFind(database = "pathway", query = "phosphorus")

  # getting "Nitrogen metabolism" pathway ID
N_pathway_ID <- sub(names(N_pathway_info), pattern = "path:", replacement = "")

  # getting info for all KOs associated with "Nitrogen metabolism" pathway
N_pathway_KOs_info <- keggLink("ko", N_pathway_ID)

  # getting KO IDs alone
N_pathway_KOs <- as.character(sub(N_pathway_KOs_info, pattern="ko:", replacement=""))

length(N_pathway_KOs) # 65

Nitrogen_metabolism_KEGG_tab <- make_KO_summary_tab(N_pathway_KOs)

write.table(Nitrogen_metabolism_KEGG_tab, "kegg-tables/Kegg-N-annotations.tsv", sep="\t", quote=F, row.names=F)
