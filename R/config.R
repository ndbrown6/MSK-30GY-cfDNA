#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
suppressPackageStartupMessages(library("gdata"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("fuzzyjoin"))
suppressPackageStartupMessages(library("pander"))
suppressPackageStartupMessages(library("ggsignif"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggforce"))
suppressPackageStartupMessages(library("ggord"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("ggridges"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("superheat"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("grid"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("drc"))
suppressPackageStartupMessages(library("preseqR"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("maftools"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("pROC"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("rpart"))
suppressPackageStartupMessages(library("klaR"))


registerDoMC(8)

FRAGMENT_LENGTH_THRESHOLD = 37
N_TILES = 250
USE_ALL = TRUE

url_manifest <- "../data/manifest.txt"
url_primers <- "../data/HPV_primers.bed"
url_bed_file <- "../data/hpv_31_33_35H_16_18.bed"
url_preanalytical_conidtions <- "../data/pre-analytical_conditions.txt"
url_mrd_summary <- "../data/TM-2021-MSK-002_MRD_Results_Summary_20220825.txt"
url_mutation_summary <- "../data/mutation_summary.maf"
url_hpv_type <- "../data/TM-2021-MSK-002_HPV_Typing_WES_WGS.txt"
url_no_node_dissection <- "../data/patients_0_nodal_dissection.txt"
url_clinical <- "../data/clinical.txt"

url_aln_metrics <- "../data/aln_metrics.txt"
url_aln_metrics_ft <- "../data/aln_metrics_ft.txt"

url_hs_metrics <- "../data/hs_metrics.txt"
url_hs_metrics_ft <- "../data/hs_metrics_ft.txt"

url_hs_target_metrics <- "../data/hs_metrics_target.txt"
url_hs_target_metrics_ft <- "../data/hs_metrics_target_ft.txt"

url_idx_metrics <- "../data/idx_metrics.txt"
url_idx_metrics_ft <- "../data/idx_metrics_ft.txt"

url_insert_metrics <- "../data/insert_metrics.txt"
url_insert_metrics_ft <- "../data/insert_metrics_ft.txt"
url_insert_metrics_by_gene <- "../data/insert_metrics_by_gene.txt"

url_insert_summary <- "../data/insert_summary.txt"
url_insert_summary_ft <- "../data/insert_summary_ft.txt"

url_gatk_summary <- "../data/gatk_summary.txt"

target_contigs <- c("HPV-16" = "NC001526.4",
		    "HPV-18" = "NC001357.1",
		    "HPV-31" = "J04353.1",
		    "HPV-33" = "M12732.1",
		    "HPV-35" = "X74477.1")

'scientific_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^", scales::scientific_format()(x)), , fixed = TRUE))
}

'scientific_1e' <- function(x) {
	parse(text=gsub("+", "", gsub("1e", " 10^", scales::scientific_format()(x)), , fixed = TRUE))
}

'str_split' <- function(x, split, n) {
	return(invisible(unlist(strsplit(x, split = split))[n]))
}