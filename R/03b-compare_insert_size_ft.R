#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(!is.na(bam_file_name_hpv)) %>%
	   dplyr::mutate(sample_uuid = paste0(sample_id_mskcc, "-", sample_id_invitae))

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

nodal_dissection_smry = readr::read_tsv(file = url_no_node_dissection, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::rename(patient_id_mskcc = sample)

manifest = manifest %>%
	   dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	   dplyr::mutate(patient_id_mskcc = case_when(
		   is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
		   TRUE ~ patient_id_mskcc
	   )) %>%
	   dplyr::mutate(sample_name = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	   dplyr::left_join(hpv_smry, by = "patient_id_mskcc")

insert_size_metrics = readr::read_tsv(file = url_insert_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert()

insert_size_metrics_ft = readr::read_tsv(file = url_insert_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      	 readr::type_convert()

insert_size_smry = readr::read_tsv(file = url_insert_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		   readr::type_convert()

insert_size_smry_ft = readr::read_tsv(file = url_insert_summary_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert()

insert_size_smry = insert_size_smry %>%
		   dplyr::group_by(SAMPLE_NAME) %>%
		   dplyr::summarize(TOTAL_READS = sum(READ_COUNT)) %>%
		   dplyr::ungroup() %>%
		   dplyr::right_join(insert_size_smry, by = "SAMPLE_NAME") %>%
		   dplyr::mutate(`%_READS` = READ_COUNT/TOTAL_READS)

insert_size_smry_ft = insert_size_smry_ft %>%
		      dplyr::group_by(SAMPLE_NAME, FRAGMENT_LENGTH) %>%
		      dplyr::summarize(TOTAL_READS = sum(READ_COUNT)) %>%
		      dplyr::ungroup() %>%
		      dplyr::right_join(insert_size_smry_ft, by = c("SAMPLE_NAME", "FRAGMENT_LENGTH")) %>%
		      dplyr::mutate(`%_READS` = READ_COUNT/TOTAL_READS)

primer_set = readr::read_tsv(file = url_primers, col_names = FALSE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert() %>%
	     dplyr::mutate(insert_size = X3 - X2)

mutation_smry = readr::read_tsv(file = url_mutation_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		readr::type_convert()

clinical = readr::read_tsv(file = url_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert()

#########################################################
# (+) Samples with mean AF > 5% from MRD assay
# (+) Samples with HPV subtype in assay
# (+) Samples with +ve MRD assay
#########################################################
smry_t_pos = mutation_smry %>%
	     dplyr::filter(FILTER == "PASS") %>%
	     dplyr::group_by(Tumor_Sample_Barcode) %>%
	     dplyr::summarize(mean_af = mean(t_maf)) %>%
	     dplyr::ungroup() %>%
	     dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
	     ## Mean AF > 5%
	     dplyr::filter(mean_af > (5/100)) %>%
	     dplyr::left_join(manifest, by = "sample_name") %>%
	     ## Known HPV
	     dplyr::filter(hpv_type_wes_wgs %in% names(target_contigs)) %>%
	     ## Positive MRD assay
	     dplyr::left_join(mrd_smry %>%
			      dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			      dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
	     dplyr::filter(`MRD-Landmark_Result` == "PRESENT") %>%
	     dplyr::mutate(Is_ctDNA = "+ve")

#########################################################
# (-) Samples with mean AF < 0.01% from MRD assay
# (-) Samples with max AF < 0.1% from MRD assay
# (-) Patients with no nodal dissection ≥ 2 years
# (-) No duplicate patients
# (+) Samples with -ve MRD assay
#########################################################
smry_t_neg = mutation_smry %>%
	     dplyr::filter(FILTER == "PASS") %>%
	     dplyr::group_by(Tumor_Sample_Barcode) %>%
	     dplyr::summarize(mean_af = mean(t_maf),
			      max_af = max(t_maf)) %>%
	     dplyr::ungroup() %>%
	     dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
	     ## Mean AF < 0.01%
	     dplyr::filter(mean_af < (0.01/100)) %>%
	     ## Max AF < 0.1%
	     dplyr::filter(max_af < (0.1/100)) %>%
	     dplyr::left_join(manifest, by = "sample_name") %>%
	     dplyr::left_join(nodal_dissection_smry, by = "patient_id_mskcc") %>%
	     ## No nodal dissection ≥ 2 years
	     dplyr::filter(nd_event == 0 & timepoint_days_since_start_of_RT > 730) %>%
	     ## No duplicate patients
	     dplyr::left_join(mutation_smry %>%
	     		      dplyr::filter(FILTER == "PASS") %>%
			      dplyr::group_by(Tumor_Sample_Barcode) %>%
			      dplyr::summarize(mean_af = mean(t_maf),
					       max_af = max(t_maf)) %>%
			      dplyr::ungroup() %>%
			      dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
			      ## Mean AF < 0.01%
			      dplyr::filter(mean_af < (0.01/100)) %>%
			      ## Max AF < 0.1%
			      dplyr::filter(max_af < (0.1/100)) %>%
			      dplyr::left_join(manifest, by = "sample_name") %>%
			      dplyr::left_join(nodal_dissection_smry, by = "patient_id_mskcc") %>%
			      ## No nodal dissection ≥ 2 years
			      dplyr::filter(nd_event == 0 & timepoint_days_since_start_of_RT > 730) %>%
			      dplyr::group_by(patient_id_mskcc) %>%
			      dplyr::summarize(timepoint_days_since_start_of_RT = max(timepoint_days_since_start_of_RT)) %>%
			      dplyr::mutate(is_selected = TRUE), by = c("patient_id_mskcc", "timepoint_days_since_start_of_RT")) %>%
	     dplyr::filter(!is.na(is_selected)) %>%
	     dplyr::select(-is_selected) %>%
	     ## Negative MRD assay
	     dplyr::left_join(mrd_smry %>%
			      dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			      dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
	     dplyr::filter(`MRD-Landmark_Result` == "ABSENT") %>%
	     dplyr::mutate(Is_ctDNA = "-ve")

smry_ft = dplyr::bind_rows(smry_t_pos %>%
			   dplyr::select(all_of(intersect(colnames(smry_t_pos), colnames(smry_t_neg)))),
			   smry_t_neg %>%
			   dplyr::select(all_of(intersect(colnames(smry_t_pos), colnames(smry_t_neg)))))

plot_ = insert_size_smry %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::mutate(Is_ctDNA = case_when(
		Is_ctDNA == "+ve" ~ "ctDNA +ve",
		Is_ctDNA == "-ve" ~ "ctDNA -ve"
	)) %>%
	ggplot(aes(x = INSERT_SIZE, y = READ_COUNT, group = sample_name)) +
	geom_step(stat = "identity", size = .25, alpha = .55, color = "#377eb8") +
	scale_x_sqrt(limits = c(10, 500),
		     breaks = c(10, 25, 50, 100, 200, 300, 400, 500),
		     labels = c(10, 25, 50, 100, 200, 300, 400, 500)) +
	scale_y_continuous(limits = c(0, NA),
		     	   labels = scientific_10) +
	xlab("Insert Size (bp)") +
	ylab("Number of Aligned Read Pairs") +
	facet_wrap(~Is_ctDNA, ncol = 1, scales = "free") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      strip.background = element_rect(colour="white", fill="white")) +
	geom_rug(data = primer_set, mapping = aes(x = insert_size), size = .2, , sides = "b", inherit.aes = FALSE)

pdf(file = "../res/Insert_Size_Distribution_no_Filtering.pdf", width = 3.25, height = 3.25)
print(plot_)
dev.off()

insert_size_metrics = readr::read_tsv(file = url_insert_metrics_by_gene, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert() %>%
		      dplyr::rename(sample_name = SAMPLE_NAME) %>%
		      dplyr::filter(TARGET_REGION != "NC001526.4:3372-4794") %>%
		      dplyr::filter(grepl("NC001526.4", TARGET_REGION, fixed = TRUE)) %>%
		      dplyr::group_by(sample_name) %>%
		      dplyr::summarize(MEAN_INSERT_SIZE = mean(MEAN_INSERT_SIZE)) %>%
		      dplyr::left_join(smry_ft %>%
				       dplyr::select(sample_name, Is_ctDNA),
				       by = "sample_name") %>%
		      dplyr::mutate(Is_ctDNA = case_when(
			      is.na(Is_ctDNA) ~ "Unknown",
			      TRUE ~ Is_ctDNA
		      )) %>%
		      dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve", "Unknown"), ordered = TRUE)) %>%
		      dplyr::left_join(manifest, by = "sample_name") %>%
		      dplyr::filter(
			      (Is_ctDNA %in% c("+ve", "-ve") & hpv_type_wes_wgs == "HPV-16") |
			      (Is_ctDNA == "Unknown" & hpv_type_wes_wgs %in% c("HPV-18", "HPV-33", "HPV-35", "HPV-58"))
		      )

plot_ = insert_size_metrics %>%
	ggplot(aes(x = Is_ctDNA, y = MEAN_INSERT_SIZE, fill = Is_ctDNA)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white") +
	geom_jitter(stat = "identity", width = .15, height = 0, shape = 21, alpha = .85, size = 3) +
	scale_fill_manual(values = c("+ve" = "#377eb8",
				     "-ve" = "#377eb8",
				     "Unknown" = "#e41a1c")) +
	xlab("") +
	ylab("HPV-16 Mean Insert Size (bp)") +
	scale_x_discrete() +
	scale_y_continuous(limits = c(20, 130),
		      	   breaks = c(25, 50, 75, 100, 125),
			   labels = c(25, 50, 75, 100, 125)) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = 130,
		    tip_length = 0.01) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20), size = 14),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = FALSE)

pdf(file = "../res/Insert_Size_Distribution_Fraction_Primer_Filtered_with_Unknown.pdf", width = 3.25, height = 3.25)
print(plot_)
dev.off()


insert_size_metrics = readr::read_tsv(file = url_insert_metrics_by_gene, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert() %>%
		      dplyr::filter(TARGET_REGION != "NC001526.4:3372-4794") %>%
		      dplyr::rename(sample_name = SAMPLE_NAME) %>%
		      dplyr::mutate(contig = unlist(lapply(TARGET_REGION, function(x) { strsplit(x, ":", fixed = TRUE)[[1]][1] }))) %>%
		      dplyr::group_by(sample_name, contig) %>%
		      dplyr::summarize(mean_insert_size = mean(MEAN_INSERT_SIZE)) %>%
		      dplyr::left_join(dplyr::tibble(contig = target_contigs,
				       		     chromosome = names(target_contigs)),
				       by = "contig") %>%
		      dplyr::left_join(manifest, by = "sample_name") %>%
		      dplyr::filter(chromosome == hpv_type_wes_wgs) %>%
		      dplyr::filter(timepoint_days_since_start_of_RT<0) %>%
		      dplyr::left_join(mutation_smry %>%
				       dplyr::filter(FILTER == "PASS") %>%
				       dplyr::group_by(Tumor_Sample_Barcode) %>%
				       dplyr::summarize(mean_af = mean(t_maf)) %>%
				       dplyr::ungroup() %>%
				       dplyr::rename(sample_name = Tumor_Sample_Barcode),
				       by = "sample_name")

idx_metrics_ft = readr::read_tsv(file = url_idx_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      	 readr::type_convert() %>%
		 dplyr::select(sample_name = SAMPLE_NAME, contig = CHROMOSOME, aligned_reads = ALIGNED_READS, fragment_length = FRAGMENT_LENGTH) %>%
		 dplyr::left_join(dplyr::tibble(contig = target_contigs,
				       		chromosome = names(target_contigs)),
				  by = "contig") %>%
		 dplyr::filter(fragment_length == 37) %>%
		 dplyr::filter(!is.na(chromosome)) %>%
		 dplyr::left_join(manifest, by = "sample_name") %>%
		 dplyr::filter(chromosome == hpv_type_wes_wgs) %>%
		 dplyr::filter(timepoint_days_since_start_of_RT<0) %>%
		 dplyr::left_join(mutation_smry %>%
				       dplyr::filter(FILTER == "PASS") %>%
				       dplyr::group_by(Tumor_Sample_Barcode) %>%
				       dplyr::summarize(mean_af = mean(t_maf)) %>%
				       dplyr::ungroup() %>%
				       dplyr::rename(sample_name = Tumor_Sample_Barcode),
				       by = "sample_name")

plot_ = insert_size_metrics %>%
	ggplot(aes(x = mean_insert_size, y = (mean_af*100)+(1E-3), color = chromosome)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, color = "goldenrod3", size = 1.75, se = FALSE) +
	geom_point(stat = "identity", shape = 21, fill = "white", alpha = .75, size = 2.5) +
	stat_cor(data = insert_size_metrics,
		 mapping = aes(x = mean_insert_size, y = (mean_af*100)+(1E-3)),
		 method = "spearman", size = 4, inherit.aes = FALSE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_continuous(limits = c(25, 125)) +
	scale_y_log10(limits = c(1e-3, 110),
		      labels = scientific_10) +
	xlab("Mean HPV Insert Size (bp)") +
	ylab("Mean PCM AF (%)") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = guide_legend(title = " "))

pdf(file = "../res/Mean_Insert_Size_Mean_AF.pdf", width = 4.5, height = 3.25)
print(plot_)
dev.off()

plot_ = idx_metrics_ft %>%
	ggplot(aes(x = aligned_reads, y = (mean_af*100)+(1E-3), color = chromosome)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, color = "goldenrod3", size = 1.75, se = FALSE) +
	geom_point(stat = "identity", shape = 21, fill = "white", alpha = .75, size = 2.5) +
	stat_cor(data = idx_metrics_ft,
		 mapping = aes(x = aligned_reads, y = (mean_af*100)+(1E-3)),
		 method = "spearman", size = 4, inherit.aes = FALSE) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(limits = c(1e-3, 110),
		      labels = scientific_10) +
	xlab("Aligned HPV Read Pairs") +
	ylab("Mean PCM AF (%)") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = guide_legend(title = " "))

pdf(file = "../res/Total_Reads_Mean_AF.pdf", width = 4.5, height = 3.25)
print(plot_)
dev.off()

fit_ = lm(formula = aligned_reads ~ mean_af, data = idx_metrics_ft %>%
	  					    dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	  					    dplyr::mutate(mean_af = log10((mean_af*100)+1e-3),
								  aligned_reads = log10(aligned_reads)) %>%
	  					    dplyr::select(aligned_reads, mean_af) %>%
	  					    as.data.frame())

resi_ = dplyr::tibble(patient_id_mskcc = idx_metrics_ft %>%
		      			 dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
		      			 .[["patient_id_mskcc"]],
		      mean_af = fit_$model$mean_af,
		      aligned_reads = fit_$model$aligned_reads,
		      aligned_reads_residuals = fit_$residuals,
		      hpv_type_wes_wgs = idx_metrics_ft %>%
		      			 dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
		      			 .[["hpv_type_wes_wgs"]])

#########################################################
# Primary Tumor HPV copy number - Panel
#########################################################
plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	ggplot(aes(x = hpv_panel_copynumber, y = aligned_reads)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, color = "goldenrod3", size = 1.5) +
	geom_point(stat = "identity", shape = 21, color = "#1b9e77", fill = "white", alpha = .85, size = 2.5) +
	scale_x_log10() +
	scale_y_continuous(limits = c(2, 7),
			   breaks = c(2, 3, 4, 5, 6, 7),
			   labels = scientific_10(10^c(2, 3, 4, 5, 6, 7))) +
	xlab("Primary Tumor HPV Copies") +
	ylab("cfDNA Aligned HPV Read Pairs") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/HPV_Copy_Number_Total_Read_Pairs.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	ggplot(aes(x = hpv_panel_copynumber, y = (10^mean_af))) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, color = "goldenrod3", size = 1.5) +
	geom_point(stat = "identity", shape = 21, color = "#1b9e77", fill = "white", alpha = .85, size = 2.5) +
	scale_x_log10() +
	scale_y_log10(limits = c(1e-3, 110),
		      breaks = c(.001, .01, .1, 1, 10, 100),
		      labels = c("ND", ".01", ".1", "1", "10", "100")) +
	xlab("Primary Tumor HPV Copies") +
	ylab("Mean PCM AF (%)\n") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/HPV_Copy_Number_Mean_AF.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	ggplot(aes(x = hpv_panel_copynumber, y = aligned_reads_residuals)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, color = "goldenrod3", size = 1.5) +
	geom_point(stat = "identity", shape = 21, color = "#1b9e77", fill = "white", alpha = .85, size = 2.5) +
	scale_x_log10() +
	scale_y_continuous(limits = c(-2.0, 2.5),
			   breaks = c(-2, -1, 0, 1, 2),
			   labels = scientific_10(10^c(-2, -1, 0, 1, 2))) +
	xlab("Primary Tumor HPV copies") +
	ylab("cfDNA Adjusted HPV Read Pairs") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/HPV_Copy_Number_Total_Reads_Residuals.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

#########################################################
# Primary Tumor HPV copy number - WGS
#########################################################
plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::mutate(hpv_wgs_copynumber = ifelse(hpv_wgs_copynumber<.1, .1, hpv_wgs_copynumber)) %>%
	ggplot(aes(x = hpv_wgs_copynumber, y = aligned_reads)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, color = "goldenrod3", size = 1.5) +
	geom_point(stat = "identity", shape = 21, color = "#1b9e77", fill = "white", alpha = .85, size = 2.5) +
	scale_x_log10(breaks = c(.1, 1, 10, 100),
		      labels = c(".1", "1", "10", "100")) +
	scale_y_continuous(limits = c(2, 7),
			   breaks = c(2, 3, 4, 5, 6, 7),
			   labels = scientific_10(10^c(2, 3, 4, 5, 6, 7))) +
	xlab("Primary Tumor HPV copies") +
	ylab("cfDNA Aligned HPV Read Pairs") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/HPV_Copy_Number_Total_Reads_Whole_Genome_Sequencing.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::mutate(hpv_wgs_copynumber = ifelse(hpv_wgs_copynumber<.1, .1, hpv_wgs_copynumber)) %>%
	ggplot(aes(x = hpv_wgs_copynumber, y = (10^mean_af))) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, color = "goldenrod3", size = 1.5) +
	geom_point(stat = "identity", shape = 21, color = "#1b9e77", fill = "white", alpha = .85, size = 2.5) +
	scale_x_log10(breaks = c(.1, 1, 10, 100),
		      labels = c(".1", "1", "10", "100")) +
	scale_y_log10(limits = c(1e-3, 110),
		      labels = scientific_10) +
	xlab("Primary Tumor HPV copies") +
	ylab("Mean PCM AF (%)\n") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/HPV_Copy_Number_Mean_AF_Whole_Genome_Sequencing.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::mutate(hpv_wgs_copynumber = ifelse(hpv_wgs_copynumber<.1, .1, hpv_wgs_copynumber)) %>%
	ggplot(aes(x = hpv_wgs_copynumber, y = aligned_reads_residuals)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, color = "goldenrod3", size = 1.5) +
	geom_point(stat = "identity", shape = 21, color = "#1b9e77", fill = "white", alpha = .85, size = 2.5) +
	scale_x_log10(breaks = c(.1, 1, 10, 100),
		      labels = c(".1", "1", "10", "100")) +
	scale_y_continuous(limits = c(-2.0, 2.5),
			   breaks = c(-2, -1, 0, 1, 2),
			   labels = scientific_10(10^c(-2, -1, 0, 1, 2))) +
	xlab("Primary Tumor HPV copies") +
	ylab("cfDNA Adjusted HPV Read Pairs") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/HPV_Copy_Number_Total_Reads_Residuals_Whole_Genome_Sequencing.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

#########################################################
# Age
#########################################################
plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::mutate(age_cat = case_when(
		Age < quantile(Age, 1/3, na.rm = TRUE) ~ "1",
		Age >= quantile(Age, 1/3, na.rm = TRUE) & Age < quantile(Age, 2/3, na.rm = TRUE) ~ "2",
		Age > quantile(Age, 2/3, na.rm = TRUE) ~ "3",
		TRUE ~ ""
	)) %>%
	dplyr::filter(age_cat != "") %>%
	ggplot(aes(x = age_cat, y = aligned_reads)) +
	geom_violin(stat = "ydensity", draw_quantiles = c(.25, .5, .75), trim = FALSE, scale = "area") +
	scale_x_discrete(breaks = c(1, 2, 3),
			 labels = c(expression("1"^st~"Tertile"), expression("2"^nd~"Tertile"), expression("3"^rd~"Tertile"))) +
	scale_y_continuous(limits = c(1, 11),
			   breaks = c(2, 3, 4, 5, 6, 7),
			   labels = scientific_10(10^c(2, 3, 4, 5, 6, 7))) +
	xlab("Age") +
	ylab("cfDNA Aligned HPV Read Pairs") +
	geom_signif(stat = "signif",
		    comparisons = list(c(1, 2),
				       c(1, 3),
				       c(2, 3)),
		    test = "wilcox.test",
		    test.args = list(alternative = "greater"),
		    y_position = c(8, 9, 10)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Age_Total_Read_Pairs.pdf", width = 3.35, height = 3.25)
print(plot_)
dev.off()

#########################################################
# Primary Tumor Volume
#########################################################
plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	ggplot(aes(x = PlanVol, y = aligned_reads)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, color = "goldenrod3", size = 1.5) +
	geom_point(stat = "identity", shape = 21, color = "#1b9e77", fill = "white", alpha = .85, size = 2.5) +
	scale_x_log10() +
	scale_y_continuous(limits = c(2, 7),
			   breaks = c(2, 3, 4, 5, 6, 7),
			   labels = scientific_10(10^c(2, 3, 4, 5, 6, 7))) +
	xlab(expression("Primary Tumor Volume ("~cm^3~")")) +
	ylab("cfDNA Aligned HPV Read Pairs") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Tumor_Volume_Total_Reads.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

#########################################################
# Tumor ssGSEA Cell Cycle
#########################################################
plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	ggplot(aes(x = ssGSEA_pathway_CELL_CYCLE, y = aligned_reads)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, color = "goldenrod3", size = 1.5) +
	geom_point(stat = "identity", shape = 21, color = "#1b9e77", fill = "white", alpha = .85, size = 2.5) +
	scale_x_continuous(limits = c(.35, .60)) +
	scale_y_continuous(limits = c(2, 7),
			   breaks = c(2, 3, 4, 5, 6, 7),
			   labels = scientific_10(10^c(2, 3, 4, 5, 6, 7))) +
	xlab("Primary Tumor Cell Cycle Score") +
	ylab("cfDNA Aligned HPV Read Pairs") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Tumor_ssGSEA_Total_Reads.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

data_ = idx_metrics_ft %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::mutate(mean_af = mean_af*100,
		      aligned_reads = log10(aligned_reads)) %>%
        dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
        dplyr::select(aligned_reads,
		      mean_af,
		      cfdna_concentration = concentration_ng_uL,
		      hpv_copynumber = hpv_panel_copynumber,
		      tumor_purity = purity.wes,
		      tumor_volume = PlanVol,
		      tumor_size = primary_tumor_size_cm,
		      age = Age,
		      t_stage = TStage,
		      n_stage = NStage,
		      smoking_status = Smoking2,
		      hypoxia = Simplified_hypoxia_group,
		      ssgsea_cell_cycle = ssGSEA_pathway_CELL_CYCLE,
		      ssgsea_cell_death = ssGSEA_pathway_CELL_DEATH) %>%
	readr::type_convert() %>%
	dplyr::mutate(mean_af = scale(mean_af),
		      cfdna_concentration = scale(cfdna_concentration),
		      hpv_copynumber = scale(hpv_copynumber),
		      tumor_volume = scale(tumor_volume),
		      tumor_size = scale(tumor_size),
		      age = scale(age),
		      ssgsea_cell_cycle = scale(ssgsea_cell_cycle),
		      ssgsea_cell_death = scale(ssgsea_cell_death)) %>%
        tidyr::drop_na() %>%
        as.data.frame()

fit_ = lm(formula = aligned_reads ~ ., data = data_)
res_ = summary(fit_)
aic_ = MASS::stepAIC(object = lm(formula = aligned_reads ~ ., data = data_), direction = "backward")

plot_ = res_ %>%
	.[["coefficients"]] %>%
	as.data.frame() %>%
	tibble::rownames_to_column("variable") %>%
	dplyr::as_tibble() %>%
	dplyr::filter(variable != "(Intercept)") %>%
	dplyr::mutate(is_significant = ifelse(`Pr(>|t|)`<.1, "Yes", "No")) %>%
	dplyr::mutate(variable_cat = case_when(
		variable == "mean_af" ~ "cfDNA",
		variable == "cfdna_concentration" ~ "cfDNA",
		variable == "hpv_copynumber" ~ "Tumor",
		variable == "tumor_purity" ~ "Tumor",
		variable == "tumor_volume" ~ "Tumor",
		variable == "tumor_size" ~ "Tumor",
		variable == "age" ~ "Patient",
		variable == "t_stageT2" ~ "Patient",
		variable == "t_stageT2" ~ "Patient",
		variable == "n_stageN2a" ~ "Patient",
		variable == "n_stageN2b" ~ "Patient",
		variable == "n_stageN2c" ~ "Patient",
		variable == "smoking_statusYes" ~ "Patient",
		variable == "hypoxianever_hypoxic" ~ "Tumor",
		variable == "hypoxiapersistent" ~ "Tumor",
		variable == "ssgsea_cell_cycle" ~ "Tumor",
		variable == "ssgsea_cell_death" ~ "Tumor"
	)) %>%
	dplyr::mutate(variable_cat = factor(variable_cat, levels = rev(c("cfDNA", "Tumor", "Patient")), ordered = TRUE)) %>%
	dplyr::arrange(variable_cat, Estimate) %>%
	dplyr::mutate(variable = factor(variable, levels = variable, ordered = TRUE)) %>%
	ggplot(aes(x = variable, ymin = 0, ymax = Estimate, y = Estimate, fill = is_significant, size = -log10(`Pr(>|t|)`), shape = variable_cat)) +
	geom_linerange(stat = "identity", size = .5) +
	geom_hline(yintercept = 0, size = 1) +
	geom_point(stat = "identity") +
	xlab("") +
	ylab("") +
	scale_fill_manual(values = c("#bdbdbd", "#e41a1c")) +
	scale_shape_manual(values = c("cfDNA" = 21, "Tumor" = 22, "Patient" = 23)) +
	scale_x_discrete() +
	scale_y_continuous() +
	coord_flip() +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = FALSE,
	       shape = guide_legend(title = "Variable Category", order = 1),
	       size = guide_legend(title = expression(-Log[10]~"p-value")))

pdf(file = "../res/Linear_Regression_Coefficients.pdf", width = 6, height = 3.5)
print(plot_)
dev.off()

variables_filtered = c("All", "-hypoxia", "-tumor_size", "-tumor_purity", "-smoking_status",
		       "-cfdna_concentration", "-n_stage", "-t_stage", "-ssgsea_cell_death", "-ssgsea_cell_cycle",
		       "-age", "-hpv_copynumber", "-mean_af")
extracted_aic = vector(mode = "numeric", length = length(variables_filtered))
extracted_aic[1] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -hypoxia)
extracted_aic[2] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -tumor_size)
extracted_aic[3] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -tumor_purity)
extracted_aic[4] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -smoking_status)
extracted_aic[5] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -cfdna_concentration)
extracted_aic[6] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -n_stage)
extracted_aic[7] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -t_stage)
extracted_aic[8] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -ssgsea_cell_death)
extracted_aic[9] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -ssgsea_cell_cycle)
extracted_aic[10] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -age)
extracted_aic[11] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -hpv_copynumber)
extracted_aic[12] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -mean_af)
extracted_aic[13] = extractAIC(fit_)[2]

plot_ = dplyr::tibble(variables = variables_filtered,
		      aic = extracted_aic) %>%
	dplyr::mutate(index = 1:n()) %>%
	ggplot(aes(x = index, y = aic)) +
	geom_line(stat = "identity", size = .95) +
	geom_point(stat = "identity", color = "black", fill = "white", shape = 21, size = 3) +
	scale_x_continuous(breaks = 1:length(variables_filtered),
			   labels = paste0(variables_filtered, " ", 13:1)) +
	scale_y_continuous() +
	xlab("") +
	ylab("AIC") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust=1),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Linear_Regression_Stepwise_AIC.pdf", width = 3.35, height = 4)
print(plot_)
dev.off()


idx_metrics_ft = readr::read_tsv(file = url_idx_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      	 readr::type_convert() %>%
		 dplyr::select(sample_name = SAMPLE_NAME, contig = CHROMOSOME, aligned_reads = ALIGNED_READS, fragment_length = FRAGMENT_LENGTH) %>%
		 dplyr::left_join(dplyr::tibble(contig = target_contigs,
				       		chromosome = names(target_contigs)),
				  by = "contig") %>%
		 dplyr::filter(fragment_length == 37) %>%
		 dplyr::filter(!is.na(chromosome)) %>%
		 dplyr::left_join(manifest, by = "sample_name") %>%
		 dplyr::filter(chromosome == hpv_type_wes_wgs) %>%
		 dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
			 timepoint_days_since_start_of_RT >= 0 & timepoint_days_since_start_of_RT < 7 ~ "wk0",
			 timepoint_days_since_start_of_RT >= 7 & timepoint_days_since_start_of_RT < 14 ~ "wk1",
			 timepoint_days_since_start_of_RT >= 14 & timepoint_days_since_start_of_RT < 21 ~ "wk2",
			 timepoint_days_since_start_of_RT >= 21 & timepoint_days_since_start_of_RT < 28 ~ "wk3",
			 timepoint_days_since_start_of_RT >= 28 & timepoint_days_since_start_of_RT < 35 ~ "wk4",
			 timepoint_days_since_start_of_RT >= 35 & timepoint_days_since_start_of_RT < 42 ~ "wk5",
			 timepoint_days_since_start_of_RT >= 42 & timepoint_days_since_start_of_RT < 49 ~ "wk6",
			 timepoint_days_since_start_of_RT >= 49 ~ "wk7+",
			 TRUE ~ "Pre-treatment"
		 )) %>%
		 dplyr::left_join(mutation_smry %>%
				  dplyr::filter(FILTER == "PASS") %>%
				  dplyr::group_by(Tumor_Sample_Barcode) %>%
				  dplyr::summarize(mean_af = mean(t_maf)) %>%
				  dplyr::ungroup() %>%
				  dplyr::rename(sample_name = Tumor_Sample_Barcode),
				  by = "sample_name")

fit_ = lm(formula = aligned_reads ~ mean_af, data = idx_metrics_ft %>%
	  					    dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	  					    dplyr::mutate(mean_af = log10((mean_af*100 + 1e-3)),
								  aligned_reads = log10(aligned_reads + 1e-3)) %>%
	  					    dplyr::select(aligned_reads, mean_af) %>%
	  					    as.data.frame())

resi_ = dplyr::tibble(patient_id_mskcc = idx_metrics_ft %>%
		      			 dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
		      			 .[["patient_id_mskcc"]],
		      mean_af = fit_$model$mean_af,
		      aligned_reads = fit_$model$aligned_reads,
		      aligned_reads_residuals = fit_$residuals,
		      hpv_type_wes_wgs = idx_metrics_ft %>%
		      			 dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
		      			 .[["hpv_type_wes_wgs"]],
		      timepoint_weeks_since_start_of_RT = idx_metrics_ft %>%
		      					  dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
		      					  .[["timepoint_weeks_since_start_of_RT"]])

plot_ = resi_ %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("wk1", "wk2", "wk3", "wk5")) %>%
	ggplot(aes(x = hpv_panel_copynumber, y = aligned_reads_residuals)) +
	geom_smooth(stat = "smooth", method = "glm", formula = y ~ x, se = FALSE, color = "goldenrod3", size = 1.5) +
	geom_point(stat = "identity", shape = 21, color = "#1b9e77", fill = "white", alpha = .85, size = 2.5) +
	scale_x_log10(limits = c(.4, 250),
		      breaks = c(1, 10, 100),
		      labels = c("1", "10", "100")) +
	scale_y_continuous(limits = c(-2.0, 2.5),
			   breaks = c(-2, -1, 0, 1, 2),
			   labels = scientific_10(10^c(-2, -1, 0, 1, 2))) +
	xlab("Primary Tumor HPV Copies") +
	ylab("cfDNA Adjusted HPV Read Pairs") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10),
	      panel.spacing = unit(1, "lines"),
	      strip.background = element_rect(colour="white", fill="white")) +
	facet_wrap(~timepoint_weeks_since_start_of_RT, scales = "free")

pdf(file = "../res/HPV_Copy_Number_Total_Reads_Residuals_bywk.pdf", width = 3.35*1.5, height = 3*1.5)
print(plot_)
dev.off()
