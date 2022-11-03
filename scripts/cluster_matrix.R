library(tidyverse)
library(pheatmap)
library(RColorBrewer)

cluster_df <- read_tsv("output/clusters_0.70.tsv", col_names=c("cluster", "target"))
hsp60_df <- read_csv("output/amplicon_sequences.csv")

cycl_strains <- select(filter(hsp60_df, genotype=="cycl_bloomer"), strain)

day_annotation <- tribble( ~strain, ~day,
    "10N.286.45.B7","286",
    "10N.286.46.B3","286",
    "10N.286.46.B6","286",
    "10N.286.46.C7","286",
    "10N.286.48.C5","286",
    "10N.286.48.E9","286",
   "10N.286.49.E10","286",
           "1F-111",NA,
           "1F-175",NA,
           "1F-273",NA,
            "1F-53",NA,
            "1F-97",NA,
        "237312C11","237",
        "237312D05","237",
        "237312G08","237",
        "239311A09","239",
        "239311B03","239",
        "239311C11","239",
        "239311D07","239",
        "239311D05","239",
        "239311E03","239",
        "239311F08","239",
        "239311G07","239",
        "239311H11","239",
        "239312A01","239",
        "239312H01","239",
         "24731027","247",
         "24731032","247",
         "24731124","247",
         "24731162","247",
         "24731185","247",
           "FF-160",NA,
         "24731082","247"
)
color_structure <- structure(c("#D55E00","#984ea3","#56B4E9","#009E73"),
  names=c(237,239,247,286))
  annotation_color_data <- list(day=color_structure)

cluster_matrix <- cluster_df %>%
    separate(target, into=c("strain"), extra="drop", sep="_") %>%
    mutate(presence=1) %>%
    group_by(strain, cluster) %>%
    filter(row_number()==1) %>%
    ungroup() %>%
    semi_join(cycl_strains) %>%
    group_by(cluster) %>%
    {filter(., n() >= 1)} %>%
    mutate(cluster_group = group_indices()) %>%
    mutate(cluster_label = paste0("cluster_", str_pad(as.character(cluster_group), width=5, pad="0"))) %>%
    ungroup() %>%
    select(-cluster, -cluster_group) %>%
    group_by(cluster_label) %>%
    filter(sum(presence) > 2, sum(presence) < nrow(cycl_strains) - 1) %>%
    ungroup() %>%
    spread(cluster_label, presence, fill=0) %>%
    as.data.frame %>%
    column_to_rownames("strain")

annotation_row_data <- tibble(strain=rownames(cluster_matrix)) %>% left_join(select(day_annotation, strain, day)) %>%
    as.data.frame() %>% column_to_rownames("strain")

pheatmap(cluster_matrix,
  annotation_row=annotation_row_data,
  annotation_colors=annotation_color_data,
  show_colnames=FALSE,
  treeheight_col=0,
  color=brewer.pal("Greys", n=3),
  border_color="#ffffff",
  main="V. cyclitrophicus bloomer hsp60 flexible genome", filename=snakemake@output[["genes"]], cellheight = 7, cellwidth=0.25, fontsize=8)
