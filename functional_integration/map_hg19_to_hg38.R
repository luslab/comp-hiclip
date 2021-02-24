library(rtracklayer)
library(data.table)
library(tidyverse)
library(stringr)
library(rio)


# ==========
# Load data
# ==========

# Load data (mapped to hg19) and get gene ids from (Mukherjee et al. 2017)
work.dir <- "/camp/home/iosubi/home/users/iosubi/projects/comp_hiclip/gene_id_mapping/"
genomes.dir <- "/camp/lab/luscomben/home/users/chakraa2/projects/flora/ref/human/icount/"

data_list <- import_list(paste0(work.dir,"41594_2017_BFnsmb3325_MOESM3_ESM.xlsx"), which = c("3b","3c","3d"))
data.df <- data_list %>% reduce(inner_join, by = c("Gene", "Simple"))

# hg19 annotation from (Mukherjee et al. 2017)
hg19.gr <- import.gff2(paste0(work.dir,"gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz"))
hg19.df <- as.data.frame(hg19.gr)
hg19.df <- hg19.df %>% filter(type == "gene")

# hg38, gencode v33 used in STAU1 hiCLIP annotation 
regions.gr <- import.gff2(paste0(genomes.dir,"regions.gtf.gz"))
regions.df <- as.data.frame(regions.gr)


# ==========
# Mapping: based on matching ensembl IDs and gene names
# ==========

hg19.df <- semi_join(hg19.df, data.df, by = c("gene_id" = "Gene")) # filter hg19 for genes present in the Mukherjee et al. 2017 data

hg19.map <- hg19.df %>%
  select(gene_id, gene_name) %>%
  mutate(gene_id_trimmed = str_remove(gene_id, "\\.[^.]*$"),
         gene_name_trimmed = str_remove(gene_name, "\\.[^.]*$"))

hg38.map <- regions.df %>% 
  select(gene_id, gene_name) %>%
  mutate(gene_id_trimmed = str_remove(gene_id, "\\.[^.]*$"),
         gene_name_trimmed = str_remove(gene_name, "\\.[^.]*$"))

# matched by gene_id
hg19to38map_id.df <- inner_join(hg19.map, hg38.map, by = "gene_id_trimmed", suffix = c(".hg19", ".hg38")) %>%
  distinct()


# matched by gene_name (for some genes the gene name is more stable between assemblies than the ensembl id)

# antijoin to match by name only the unmapped hg19 ids
hg19.map_name <- anti_join(hg19.map, hg19to38map_id.df, by = c("gene_id" = "gene_id.hg19"))

hg19to38map_name.df <- inner_join(hg19.map_name, hg38.map, by = "gene_name_trimmed", suffix = c(".hg19", ".hg38")) %>%
  distinct()

hg19to38map_id.df <- hg19to38map_id.df %>%
  select(-gene_id_trimmed, -gene_name_trimmed.hg19, -gene_name_trimmed.hg38)

hg19to38map_name.df <- hg19to38map_name.df %>%
  select(-gene_id_trimmed.hg19, -gene_name_trimmed, -gene_id_trimmed.hg38)

hg19to38map.df <- rbind(hg19to38map_id.df, hg19to38map_name.df)

# ==========
# Unmatched hg19 IDs
# ==========

# Extract ids & names of the umapped
hg19.unmapped <- anti_join(hg19.map, hg19to38map.df, by = c("gene_id" = "gene_id.hg19"))

hg19.unmapped <- hg19.unmapped %>%
  rename_at(vars(gene_name, gene_id),function(x) paste0(x,".hg19")) %>%
  select(gene_id.hg19, gene_name.hg19) %>%
  mutate(gene_id.hg38 = NA, gene_name.hg38 = NA)

# Append the ids & names of the unmapped
hg19to38map.df <- rbind(hg19to38map.df, hg19.unmapped)

fwrite(hg19to38map.df, sep = "\t", paste0(genomes.dir,"hg19_to_hg38_map.tsv"))
