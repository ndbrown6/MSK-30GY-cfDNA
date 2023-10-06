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
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      plot.title = element_text(hjust = 0.5))

pdf("../res/Number_Read_Pairs_by_cfDNA_Concentration_ng_muL__HPV.pdf", width = 5, height = 5)
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
	labs(title = "MRD Assay\n") +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      plot.title = element_text(hjust = 0.5))

pdf("../res/Number_Read_Pairs_by_cfDNA_Concentration_ng_muL__MRD.pdf", width = 5, height = 5)
print(plot_)
dev.off()

aln_metrics %>%
dplyr::left_join(manifest, by = "SAMPLE_NAME") %>%
dplyr::left_join(mrd_smry %>%
		 dplyr::mutate(SAMPLE_NAME = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)),
		 by = "SAMPLE_NAME") %>%
dplyr::mutate(hpv_fraction = (TOTAL_READS/2)/(Total_Sample_Complexity+(TOTAL_READS/2))) %>%
dplyr::filter(SAMPLE_NAME == "21-020-04087-GER2107031") %>%
.[["hpv_fraction"]] * 100

plot_ = idx_metrics %>%
	dplyr::mutate(is_hpv_yes_no = case_when(
		hpv_type_wes_wgs == "HPV-16" ~ "HPV-16 +ve",
		TRUE ~ "HPV-16 -ve"
	)) %>%
	ggplot(aes(x = `HPV-16`)) +
	geom_histogram(stat = "bin", bins = 25, fill = "#2c7fb8", color = NA, alpha = .35) +
	xlab("Read Pairs Aligned to HPV-16 Contig") +
	ylab("Number of Samples") +
	scale_x_log10(expand = c(0, 0),
		      labels = scientific_10) +
	scale_y_continuous() +
	theme_minimal() +
	facet_wrap(~is_hpv_yes_no) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      panel.spacing = unit(2, "lines"))
	

pdf(file = "../res/Number_Read_Pairs_Aligned_by_Contigs_HPV-16.pdf", width = 7, height = 3)
print(plot_)
dev.off()

plot_ = idx_metrics %>%
	dplyr::mutate(is_hpv_yes_no = case_when(
		hpv_type_wes_wgs == "HPV-18" ~ "HPV-18 +ve",
		TRUE ~ "HPV-18 -ve"
	)) %>%
	ggplot(aes(x = `HPV-18`)) +
	geom_histogram(stat = "bin", bins = 25, fill = "#7fcdbb", color = NA, alpha = .55) +
	xlab("Read Pairs Aligned to HPV-18 Contig") +
	ylab("Number of Samples") +
	scale_x_log10(expand = c(0, 0),
		      labels = scientific_10) +
	scale_y_continuous() +
	theme_minimal() +
	facet_wrap(~is_hpv_yes_no) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      panel.spacing = unit(2, "lines"))

pdf(file = "../res/Number_Read_Pairs_Aligned_by_Contigs_HPV-18.pdf", width = 7, height = 3)
print(plot_)
dev.off()

plot_ = idx_metrics %>%
	dplyr::mutate(is_hpv_yes_no = case_when(
		hpv_type_wes_wgs == "HPV-33" ~ "HPV-33 +ve",
		TRUE ~ "HPV-33 -ve"
	)) %>%
	ggplot(aes(x = `HPV-33`)) +
	geom_histogram(stat = "bin", bins = 25, fill = "#dd1c77", color = NA, alpha = .35) +
	xlab("Read Pairs Aligned to HPV-33 Contig") +
	ylab("Number of Samples") +
	scale_x_log10(expand = c(0, 0),
		      labels = scientific_10) +
	scale_y_continuous() +
	theme_minimal() +
	facet_wrap(~is_hpv_yes_no) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      panel.spacing = unit(2, "lines"))

pdf(file = "../res/Number_Read_Pairs_Aligned_by_Contigs_HPV-33.pdf", width = 7, height = 3)
print(plot_)
dev.off()
	
plot_ = idx_metrics %>%
	dplyr::mutate(is_hpv_yes_no = case_when(
		hpv_type_wes_wgs == "HPV-35" ~ "HPV-35 +ve",
		TRUE ~ "HPV-35 -ve"
	)) %>%
	ggplot(aes(x = `HPV-35`)) +
	geom_histogram(stat = "bin", bins = 25, fill = "#de2d26", color = NA, alpha = .35) +
	xlab("Read Pairs Aligned to HPV-35 Contig") +
	ylab("Number of Samples") +
	scale_x_log10(expand = c(0, 0),
		      labels = scientific_10) +
	scale_y_continuous() +
	theme_minimal() +
	facet_wrap(~is_hpv_yes_no) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      panel.spacing = unit(2, "lines"))

pdf(file = "../res/Number_Read_Pairs_Aligned_by_Contigs_HPV-35.pdf", width = 7, height = 3)
print(plot_)
dev.off()

ii = 1
data_ = list()
for (i in 1:(length(target_contigs)-1)) {
	for (j in (i+1):length(target_contigs)) {
		data_[[ii]] = idx_metrics %>%
			      dplyr::select(all_of(c("SAMPLE_NAME", "hpv_type_wes_wgs", names(target_contigs)[i], names(target_contigs)[j]))) %>%
			      dplyr::rename(x = names(target_contigs)[i],
					    y = names(target_contigs)[j]) %>%
			      dplyr::mutate(xlab = names(target_contigs)[i],
					    ylab = names(target_contigs)[j])
		ii = ii + 1
	}
}

plot_ = do.call(bind_rows, data_) %>%
	dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
	dplyr::left_join(manifest %>%
			 dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
			 dplyr::group_by(hpv_type_wes_wgs) %>%
			 dplyr::summarize(n = n()), by = "hpv_type_wes_wgs") %>%
	dplyr::arrange(desc(n)) %>%
	dplyr::mutate(hpv_type_wes_wgs = factor(hpv_type_wes_wgs, levels = rev(c("Unknown", "HPV-18", "HPV-58", "HPV-35", "HPV-33", "HPV-16")), ordered = TRUE)) %>%
	ggplot(aes(x +1, y +1, color = hpv_type_wes_wgs)) +
	geom_abline(slope = 1, intercept = 0, linetype = 1, size = .75, alpha = .55, color = "goldenrod3") +
	geom_vline(xintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_hline(yintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_point(stat = "identity", shape = 21, fill = NA, alpha = .55, size = 2) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("Number of Read Pairs Aligned") +
	ylab("Number of Read Pairs Aligned") +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	facet_grid(ylab ~ xlab) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(color = guide_legend(title = "Tumor\nHPV Type", override.aes = list(alpha = 1)))

pdf(file = "../res/Number_Read_Pairs_Aligned_by_Contigs.pdf", width = 9, height = 7)
print(plot_)
dev.off()


plot_ = do.call(bind_rows, data_) %>%
	dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
	dplyr::left_join(manifest %>%
			 dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
			 dplyr::group_by(hpv_type_wes_wgs) %>%
			 dplyr::summarize(n = n()), by = "hpv_type_wes_wgs") %>%
	dplyr::arrange(desc(n)) %>%
	dplyr::mutate(hpv_type_wes_wgs = factor(hpv_type_wes_wgs, levels = rev(c("Unknown", "HPV-18", "HPV-58", "HPV-35", "HPV-33", "HPV-16")), ordered = TRUE)) %>%
	dplyr::filter(hpv_type_wes_wgs != "HPV-18" & hpv_type_wes_wgs != "HPV-31" & hpv_type_wes_wgs != "Unknown") %>%
	dplyr::filter(xlab == "HPV-18" & ylab == "HPV-31") %>%
	ggplot(aes(x +1, y +1, color = hpv_type_wes_wgs)) +
	geom_abline(slope = 1, intercept = 0, linetype = 1, size = 1.5, alpha = .55, color = "#7fcdbb") +
	geom_smooth(method = "glm", formula = y ~ x, linetype = 1, size = 1.5, alpha = .55, color = "goldenrod3", se = FALSE, fullrange = TRUE) +
	geom_vline(xintercept = 1, color = "black", linetype = 3, size = .75, alpha = .75) +
	geom_hline(yintercept = 1, color = "black", linetype = 3, size = .75, alpha = .75) +
	geom_point(stat = "identity", shape = 21, fill = NA, alpha = .45, size = 2.5) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("Read Pairs Aligned to HPV-18 Contig") +
	ylab("Read Pairs Aligned to HPV-31 Contig") +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(color = guide_legend(title = "Tumor\nHPV Type", override.aes = list(alpha = 1)))

pdf(file = "../res/Number_Read_Pairs_Aligned_HPV-18_HPV_31.pdf", width = 6, height = 5)
print(plot_)
dev.off()

ii = 1
data_ = list()
for (i in 1:(length(target_contigs)-1)) {
	for (j in (i+1):length(target_contigs)) {
		data_[[ii]] = idx_metrics_ft %>%
			      dplyr::select(all_of(c("SAMPLE_NAME", "hpv_type_wes_wgs", names(target_contigs)[i], names(target_contigs)[j]))) %>%
			      dplyr::rename(x = names(target_contigs)[i],
					    y = names(target_contigs)[j]) %>%
			      dplyr::mutate(xlab = names(target_contigs)[i],
					    ylab = names(target_contigs)[j])
		ii = ii + 1
	}
}

plot_ = do.call(bind_rows, data_) %>%
	dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
	dplyr::left_join(manifest %>%
			 dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
			 dplyr::group_by(hpv_type_wes_wgs) %>%
			 dplyr::summarize(n = n()), by = "hpv_type_wes_wgs") %>%
	dplyr::arrange(desc(n)) %>%
	dplyr::mutate(hpv_type_wes_wgs = factor(hpv_type_wes_wgs, levels = rev(c("Unknown", "HPV-18", "HPV-58", "HPV-35", "HPV-33", "HPV-16")), ordered = TRUE)) %>%
	ggplot(aes(x +1, y +1, color = hpv_type_wes_wgs)) +
	geom_abline(slope = 1, intercept = 0, linetype = 1, size = .75, alpha = .55, color = "goldenrod3") +
	geom_vline(xintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_hline(yintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_point(stat = "identity", shape = 21, fill = NA, alpha = .55, size = 2) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("Number of Read Pairs Aligned") +
	ylab("Number of Read Pairs Aligned") +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	facet_grid(ylab ~ xlab) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(color = guide_legend(title = "Tumor\nHPV Type", override.aes = list(alpha = 1)))

pdf(file = "../res/Number_Read_Pairs_Aligned_by_Contigs_Primer_Filtered.pdf", width = 9, height = 7)
print(plot_)
dev.off()

plot_ = idx_metrics %>%
	reshape2::melt(id.vars = c("SAMPLE_NAME", "hpv_type_wes_wgs"), measure.vars = names(target_contigs)) %>%
	dplyr::rename(contig = variable, aligned_reads = value) %>%
	dplyr::full_join(idx_metrics_ft %>%
			 reshape2::melt(id.vars = "SAMPLE_NAME", measure.vars = names(target_contigs)) %>%
			 dplyr::rename(contig = variable, aligned_reads_ft = value), by = c("SAMPLE_NAME", "contig")) %>%
	dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
	dplyr::left_join(manifest %>%
			 dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
			 dplyr::group_by(hpv_type_wes_wgs) %>%
			 dplyr::summarize(n = n()), by = "hpv_type_wes_wgs") %>%
	dplyr::arrange(desc(n)) %>%
	dplyr::mutate(hpv_type_wes_wgs = factor(hpv_type_wes_wgs, levels = rev(c("Unknown", "HPV-18", "HPV-58", "HPV-35", "HPV-33", "HPV-16")), ordered = TRUE)) %>%
	ggplot(aes(x = aligned_reads+1, y = aligned_reads_ft+1, color = hpv_type_wes_wgs)) +
	geom_abline(slope = 1, intercept = 0, linetype = 1, size = .75, alpha = .55, color = "goldenrod3") +
	geom_vline(xintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_hline(yintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_point(stat = "identity", shape = 21, fill = NA, alpha = .55, size = 2) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("Total Number of Read Pairs") +
	ylab("Number of Read Pairs Filtered") +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	facet_grid(hpv_type_wes_wgs~contig) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(color = FALSE)
	
pdf(file = "../res/Number_Read_Pairs_Aligned_by_Contigs__Compare_Primer_Filtered.pdf", width = 8, height = 8)
print(plot_)
dev.off()
