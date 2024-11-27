#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert()

mrd_smry = readr::read_tsv(file = url_mrd_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert()

preanalytical_conditions = readr::read_tsv(file = url_preanalytical_conidtions, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			   readr::type_convert()

hpv_smry = readr::read_tsv(file = url_hpv_type, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::mutate(patient_id_mskcc = gsub(pattern = "-T", replacement = "", x= patient_name, fixed = TRUE)) %>%
	   dplyr::mutate(hpv_type_wes = gsub("HPV", "HPV-", hpv_type_wes, fixed = TRUE)) %>%
	   dplyr::mutate(hpv_type_wgs = gsub("HPV", "HPV-", hpv_type_wgs, fixed = TRUE)) %>%
	   dplyr::mutate(hpv_type_wes_wgs = gsub("HPV", "HPV-", hpv_type_wes_wgs, fixed = TRUE)) %>%
	   dplyr::mutate(hpv_type_wes_wgs = case_when(
		   is.na(hpv_type_wes_wgs) ~ "Unknown",
		   TRUE ~ hpv_type_wes_wgs
	   ))

manifest = manifest %>%
	   dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	   dplyr::mutate(patient_id_mskcc = case_when(
		   is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
		   TRUE ~ patient_id_mskcc
	   )) %>%
	   dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	   dplyr::left_join(hpv_smry, by = "patient_id_mskcc")

aln_metrics = readr::read_tsv(file = url_aln_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::filter(CATEGORY == "PAIR")

aln_metrics_ft = readr::read_tsv(file = url_aln_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      	 readr::type_convert() %>%
	      	 dplyr::filter(CATEGORY == "PAIR")

idx_metrics = readr::read_tsv(file = url_idx_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::filter(CHROMOSOME %in% target_contigs) %>%
	      reshape2::dcast(SAMPLE_NAME ~ CHROMOSOME, value.var = "ALIGNED_READS") %>%
	      dplyr::left_join(manifest, by = "SAMPLE_NAME") %>%
	      dplyr::rename(`HPV-31` = `J04353.1`,
			    `HPV-33` = `M12732.1`,
			    `HPV-18` = `NC001357.1`,
			    `HPV-16` = `NC001526.4`,
			    `HPV-35` = `X74477.1`)

idx_metrics_ft = readr::read_tsv(file = url_idx_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	         readr::type_convert() %>%
		 dplyr::filter(FRAGMENT_LENGTH == 0) %>%
	         dplyr::filter(CHROMOSOME %in% target_contigs) %>%
	         reshape2::dcast(SAMPLE_NAME ~ CHROMOSOME, value.var = "ALIGNED_READS") %>%
		 dplyr::left_join(manifest, by = "SAMPLE_NAME") %>%
		 dplyr::rename(`HPV-31` = `J04353.1`,
			       `HPV-33` = `M12732.1`,
			       `HPV-18` = `NC001357.1`,
			       `HPV-16` = `NC001526.4`,
			       `HPV-35` = `X74477.1`)

plot_ = aln_metrics %>%
	dplyr::left_join(manifest, by = "SAMPLE_NAME") %>%
	dplyr::filter(!is.na(concentration_ng_uL)) %>%
	dplyr::filter(concentration_ng_uL>.01) %>%
	dplyr::filter(TOTAL_READS < (2*1E4)) %>%
	ggplot(aes(x = concentration_ng_uL, y = TOTAL_READS/2)) +
	geom_point(stat = "identity", shape = 21, fill = NA, alpha = .55, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x, se = .90, color = "goldenrod3", alpha = .55, size = 1.5) +
	xlab(expression("cfDNA concentration  (ng"~.~mu~L^-1~")")) +
	ylab("Number of Read Pairs Aligned") +
	labs(title = "HPV Assay\n") +
	scale_x_log10(limits = c(NA, 1),
		      labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      plot.title = element_text(hjust = 0.5))

pdf("../res/Number_Read_Pairs_by_cfDNA_Concentration_ng_muL_HPV_Assay.pdf", width = 3.25, height = 3.5)
print(plot_)
dev.off()

plot_ = mrd_smry %>%
	dplyr::mutate(SAMPLE_NAME = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
	dplyr::left_join(manifest, by = "SAMPLE_NAME") %>%
	dplyr::filter(!is.na(concentration_ng_uL)) %>%
	dplyr::filter(concentration_ng_uL>.01) %>%
	ggplot(aes(x = concentration_ng_uL, y = Total_Sample_Complexity)) +
	geom_point(stat = "identity", shape = 21, fill = NA, alpha = .55, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x, se = .90, color = "goldenrod3", alpha = .55, size = 1.5) +
	xlab(expression("cfDNA concentration  (ng"~.~mu~L^-1~")")) +
	ylab("Number of Read Pairs Aligned") +
	labs(title = "PCM Assay\n") +
	scale_x_log10(limits = c(NA, 1),
		      labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      plot.title = element_text(hjust = 0.5))

pdf("../res/Number_Read_Pairs_by_cfDNA_Concentration_ng_muL_PCM_Assay.pdf", width = 3.25, height = 3.5)
print(plot_)
dev.off()
