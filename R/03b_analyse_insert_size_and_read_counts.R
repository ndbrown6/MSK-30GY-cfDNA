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
	   readr::type_convert() %>%
	   dplyr::mutate(sample_uuid = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer))

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

idx_metrics_ft = readr::read_tsv(file = url_idx_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      	 readr::type_convert()

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

plot_ = insert_size_smry_ft %>%
	dplyr::filter(FRAGMENT_LENGTH==0) %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft %>%
			 dplyr::select(sample_name, Is_ctDNA), by = "sample_name") %>%
	dplyr::filter(is.na(Is_ctDNA)) %>%
	dplyr::select(-Is_ctDNA) %>%
	dplyr::left_join(readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		    	 readr::type_convert(),
			 by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::rename(Is_ctDNA__HPV = Is_ctDNA) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::mutate(Is_ctDNA__MRD = case_when(
		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
		TRUE ~ "-ve"
	)) %>%
	dplyr::filter(Is_ctDNA__HPV == "+ve" & Is_ctDNA__MRD == "-ve") %>%
	ggplot(aes(x = INSERT_SIZE, y = READ_COUNT, group = sample_name)) +
	geom_rect(aes(xmin = min(primer_set$insert_size), xmax = max(primer_set$insert_size), ymin = 0, ymax = 495),
		  fill = "grey90", color = NA, alpha = .75) +
	geom_smooth(stat = "smooth", method = "loess", formula = y ~x, color = "#1b9e77", se = FALSE, size = .25, n = 250, span = .05) +
	scale_x_sqrt(limits = c(10, 500),
		     breaks = c(10, 25, 50, 100, 200, 300, 400, 500),
		     labels = c(10, 25, 50, 100, 200, 300, 400, 500)) +
	scale_y_continuous(limits = c(0, 500),
		     	   labels = scientific_10) +
	xlab("Insert Size (bp)") +
	ylab("Number of Aligned Read Pairs") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      strip.background = element_rect(colour="white", fill="white"))

pdf(file = "../res/Insert_Size_Distribution_Primer_Filtered_apparent_False_Positive.pdf", width = 3.25, height = 3.25/1.25)
print(plot_)
dev.off()

plot_ = insert_size_smry_ft %>%
	dplyr::filter(FRAGMENT_LENGTH==0) %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft %>%
			 dplyr::select(sample_name, Is_ctDNA), by = "sample_name") %>%
	dplyr::filter(is.na(Is_ctDNA)) %>%
	dplyr::select(-Is_ctDNA) %>%
	dplyr::left_join(readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		    	 readr::type_convert(),
			 by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::rename(Is_ctDNA__HPV = Is_ctDNA) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::mutate(Is_ctDNA__MRD = case_when(
		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
		TRUE ~ "-ve"
	)) %>%
	dplyr::filter(Is_ctDNA__HPV == "-ve" & Is_ctDNA__MRD == "+ve") %>%
	ggplot(aes(x = INSERT_SIZE, y = READ_COUNT, group = sample_name)) +
	geom_rect(aes(xmin = min(primer_set$insert_size), xmax = max(primer_set$insert_size), ymin = 0, ymax = 495),
		  fill = "grey90", color = NA, alpha = .75) +
	geom_smooth(stat = "smooth", method = "loess", formula = y ~x, color = "#1b9e77", se = FALSE, size = .25, n = 250, span = .15) +
	scale_x_sqrt(limits = c(10, 500),
		     breaks = c(10, 25, 50, 100, 200, 300, 400, 500),
		     labels = c(10, 25, 50, 100, 200, 300, 400, 500)) +
	scale_y_continuous(limits = c(0, 500),
		     	   labels = scientific_10) +
	xlab("Insert Size (bp)") +
	ylab("Number of Aligned Read Pairs") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      strip.background = element_rect(colour="white", fill="white"))

pdf(file = "../res/Insert_Size_Distribution_Primer_Filtered_apparent_False_Negative.pdf", width = 3.25, height = 3.25/1.25)
print(plot_)
dev.off()

plot_ = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	readr::type_convert() %>%
	dplyr::rename(Is_ctDNA__HPV = Is_ctDNA) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::mutate(Is_ctDNA__MRD = case_when(
		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
		TRUE ~ "-ve"
	)) %>%
	dplyr::filter(!(sample_name %in% smry_ft$sample_name)) %>%
	dplyr::mutate(class = case_when(
		Is_ctDNA__HPV == "-ve" & Is_ctDNA__MRD == "+ve" ~ "FN",
		Is_ctDNA__HPV == "-ve" & Is_ctDNA__MRD == "-ve" ~ "TN",
		Is_ctDNA__HPV == "+ve" & Is_ctDNA__MRD == "-ve" ~ "FP",
		Is_ctDNA__HPV == "+ve" & Is_ctDNA__MRD == "+ve" ~ "TP"
	)) %>%
	dplyr::filter(class %in% c("FP", "FN")) %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(!duplicated(paste0(patient_id_mskcc, ":", class))) %>%
	ggplot(aes(x = class, y = hpv_panel_copy_number)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	geom_jitter(stat = "identity", width = .15, height = 0, shape = 21, alpha = .85, size = 3, color = "#1b9e77", fill = "white") +
	scale_x_discrete(breaks = c("FN", "FP"),
			 labels = c("\nMRD +ve/\nHPV -ve", "\nMRD -ve/\nHPV +ve")) +
	scale_y_continuous(limits = c(0, 300),
			   breaks = seq(0, 300, by = 50),
			   labels = seq(0, 300, by = 50)) +
	xlab("") +
	ylab("HPV Copy Number") +
	geom_signif(stat = "signif",
		    comparisons = list(c("FP", "FN")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = 275,
		    tip_length = 0.01) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)))

pdf(file = "../res/MRD_HPV_HPV_Copy_number_Misclassifications.pdf", width = 2.75, height = 3.25)
print(plot_)
dev.off()

plot_ = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	readr::type_convert() %>%
	dplyr::rename(Is_ctDNA__HPV = Is_ctDNA) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::mutate(Is_ctDNA__MRD = case_when(
		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
		TRUE ~ "-ve"
	)) %>%
	dplyr::filter(!(sample_name %in% smry_ft$sample_name)) %>%
	dplyr::mutate(class = case_when(
		Is_ctDNA__HPV == "-ve" & Is_ctDNA__MRD == "+ve" ~ "FN",
		Is_ctDNA__HPV == "-ve" & Is_ctDNA__MRD == "-ve" ~ "TN",
		Is_ctDNA__HPV == "+ve" & Is_ctDNA__MRD == "-ve" ~ "FP",
		Is_ctDNA__HPV == "+ve" & Is_ctDNA__MRD == "+ve" ~ "TP"
	)) %>%
	dplyr::filter(class %in% c("TP", "TN")) %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::filter(!duplicated(paste0(patient_id_mskcc, ":", class))) %>%
	ggplot(aes(x = class, y = hpv_panel_copy_number)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	geom_jitter(stat = "identity", width = .15, height = 0, shape = 21, alpha = .85, size = 3, color = "#1b9e77", fill = "white") +
	scale_x_discrete(breaks = c("TN", "TP"),
			 labels = c("\nMRD -ve/\nHPV -ve", "\nMRD +ve/\nHPV +ve")) +
	scale_y_continuous(limits = c(0, 300),
			   breaks = seq(0, 300, by = 50),
			   labels = seq(0, 300, by = 50)) +
	xlab("") +
	ylab("HPV Copy Number") +
	geom_signif(stat = "signif",
		    comparisons = list(c("TP", "TN")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = 275,
		    tip_length = 0.01) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)))

pdf(file = "../res/MRD_HPV_HPV_Copy_number_Correct_Classifications.pdf", width = 2.75, height = 3.25)
print(plot_)
dev.off()

plot_ = readr::read_tsv(file = url_insert_metrics_by_gene, col_names = TRUE, col_types = cols(.default = col_character())) %>%
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
	) %>%
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

pdf(file = "../res/Mean_Insert_Size_Primer_Filtered_with_Unknown.pdf", width = 3.25, height = 3.25)
print(plot_)
dev.off()

#########################################################
# Linear Regression
# Aligned Reads ~ AF + HPV copy number
#########################################################
idx_metrics_ft = idx_metrics_ft %>%
		 dplyr::select(sample_name = SAMPLE_NAME, contig = CHROMOSOME, aligned_reads = ALIGNED_READS, fragment_length = FRAGMENT_LENGTH) %>%
		 dplyr::left_join(dplyr::tibble(contig = target_contigs,
				       		chromosome = names(target_contigs)),
				  by = "contig") %>%
		 dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD) %>%
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
				  by = "sample_name") %>%
		dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
		dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
		dplyr::filter(!is.na(hpv_panel_copy_number)) %>%
		dplyr::filter(!is.na(aligned_reads)) %>%
		dplyr::filter(!is.na(mean_af))


fit_ = stats::lm(formula = aligned_reads ~ hpv_panel_copy_number,
		 data = idx_metrics_ft %>%
		        dplyr::mutate(aligned_reads = log10(aligned_reads)) %>%
		        dplyr::mutate(hpv_panel_copy_number = log10(hpv_panel_copy_number)) %>%
		 	as.data.frame())

res_ = summary(object = lm(formula = aligned_reads ~ hpv_panel_copy_number + mean_af,
			   data = idx_metrics_ft %>%
			   	  dplyr::mutate(aligned_reads = log10(aligned_reads)) %>%
			   	  dplyr::mutate(hpv_panel_copy_number = log10(hpv_panel_copy_number))))

aic_ = MASS::stepAIC(object = lm(formula = aligned_reads ~ .,
				 data = idx_metrics_ft %>%
				 	dplyr::mutate(aligned_reads = log10(aligned_reads+1E-3)) %>%
				 	dplyr::mutate(hpv_panel_copy_number = log10(hpv_panel_copy_number)) %>%
				 	dplyr::select(aligned_reads, hpv_panel_copy_number, mean_af)), direction = "backward")

plot_ = idx_metrics_ft %>%
	ggplot(aes(x = (mean_af*100)+(1E-3), y = aligned_reads)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, color = "goldenrod3", size = 1.25, se = FALSE) +
	geom_point(stat = "identity", shape = 21, fill = "white", color = "#1b9e77", alpha = .75, size = 2.5) +
	stat_cor(data = idx_metrics_ft,
		 mapping = aes(x = (mean_af*100)+(1E-3), y = aligned_reads),
		 method = "spearman", size = 4, inherit.aes = FALSE) +
	scale_x_log10(limits = c(1e-3, 110),
		      breaks = c(.001, .01, .1, 1, 10, 100),
		      labels = c("0", ".01", ".1", "1", "10", "100")) +
	scale_y_log10(limits = 10^c(2, 7),
		      breaks = 10^c(2, 3, 4, 5, 6, 7),
		      labels = scientific_10) +
	xlab("Mean AF (%)") +
	ylab("cfDNA Aligned HPV Read Pairs") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = guide_legend(title = " "))

pdf(file = "../res/Mean_AF_Number_Read_Pairs.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

plot_ = idx_metrics_ft %>%
	dplyr::mutate(aligned_reads_residuals = 10^fit_$residuals) %>%
	ggplot(aes(x = (mean_af*100)+(1E-3), y = aligned_reads_residuals)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, color = "goldenrod3", size = 1.25, se = FALSE) +
	geom_point(stat = "identity", shape = 22, fill = "white", color = "#1b9e77", alpha = .75, size = 2.75) +
	stat_cor(data = idx_metrics_ft %>%
		 	dplyr::mutate(aligned_reads_residuals = 10^fit_$residuals),
		 mapping = aes(x = (mean_af*100)+(1E-3), y = aligned_reads_residuals),
		 method = "spearman", size = 4, inherit.aes = FALSE) +
	scale_x_log10(limits = c(1e-3, 110),
		      breaks = c(.001, .01, .1, 1, 10, 100),
		      labels = c("0", ".01", ".1", "1", "10", "100")) +
	scale_y_log10(limits = c(1E-2, 200),
		      labels = scientific_10) +
	xlab("Mean AF (%)") +
	ylab("cfDNA Adjusted HPV Read Pairs") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = guide_legend(title = " "))

pdf(file = "../res/Mean_AF_Number_Read_Pairs_Adjusted.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

plot_ = res_ %>%
	.[["coefficients"]] %>%
	as.data.frame() %>%
	tibble::rownames_to_column("variable") %>%
	dplyr::as_tibble() %>%
	dplyr::mutate(variable = case_when(
		variable == "mean_af" ~ "Mean AF (%)",
		variable == "hpv_panel_copy_number" ~ "HPV copy number",
		variable == "(Intercept)" ~ "Y intercept"
	)) %>%
	dplyr::mutate(is_significant = ifelse(`Pr(>|t|)`<.1, "Yes", "No")) %>%
	dplyr::arrange(Estimate) %>%
	dplyr::mutate(variable = factor(variable, levels = variable, ordered = TRUE)) %>%
	ggplot(aes(x = variable, ymin = 0, ymax = Estimate, y = Estimate, fill = is_significant, size = -log10(`Pr(>|t|)`))) +
	geom_linerange(stat = "identity", size = .5) +
	geom_hline(yintercept = 0, size = 1) +
	geom_point(stat = "identity", shape = 21) +
	xlab("") +
	ylab("") +
	scale_fill_manual(values = c("No" = "#bdbdbd", "Yes" = "#e41a1c")) +
	scale_x_discrete() +
	scale_y_continuous() +
	coord_flip() +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = FALSE,
	       size = guide_legend(title = expression(-Log[10]~"p-value")))

pdf(file = "../res/Linear_Regression_Coefficients_Mini.pdf", width = 6, height = 1.75)
print(plot_)
dev.off()

plot_ = dplyr::tibble(variables = c("None", "- HPV copy\nnumber", "Mean\nAF (%)"),
		      aic = c(-15.6000, -6.9827, -0.6432)) %>%
	dplyr::mutate(index = 1:n()) %>%
	ggplot(aes(x = index, y = aic)) +
	geom_line(stat = "identity", size = .95) +
	geom_point(stat = "identity", color = "black", fill = "white", shape = 21, size = 3) +
	scale_x_continuous(breaks = 1:3,
			   labels = c("None", "HPV copy\nnumber", "Mean\nAF (%)")) +
	scale_y_continuous(limits = c(-20, NA)) +
	xlab("") +
	ylab("AIC") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust=1),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Linear_Regression_Stepwise_AIC_Mini.pdf", width = 2, height = 2)
print(plot_)
dev.off()

#########################################################
# Liner Regression
# Aligned Reads ~ .
#########################################################
data_ = idx_metrics_ft %>%
	dplyr::left_join(readr::read_tsv(file = "../res/GSEA.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert(),
			 by = "patient_id_mskcc") %>%
	dplyr::mutate(mean_af = (mean_af*100)+1E-3) %>%
        dplyr::select(aligned_reads,
		      mean_af,
		      cfdna_concentration = concentration_ng_uL,
		      hpv_copynumber = hpv_panel_copy_number,
		      tumor_purity = purity_wes,
		      tumor_volume = plan_volume,
		      tumor_size = primary_tumor_size_cm,
		      age = age,
		      t_stage = t_stage,
		      n_stage = n_stage,
		      smoking_status = smoking_category_yes_never,
		      hypoxia = simplified_hypoxia_group,
		      ssGSEA_apoptosis_pos = GOBP_REGULATION_OF_MITOCHONDRIAL_MEMBRANE_PERMEABILITY_INVOLVED_IN_APOPTOTIC_PROCESS,
		      ssGSEA_apoptosis_neg = GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS,
		      ssGSEA_necrosis = REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS,
		      ssGSEA_mitosis = GOBP_MITOTIC_G1_S_TRANSITION_CHECKPOINT_SIGNALING) %>%
	readr::type_convert() %>%
	dplyr::mutate(mean_af = log10(mean_af),
		      aligned_reads = log10(aligned_reads),
		      cfdna_concentration = scale(cfdna_concentration),
		      hpv_copynumber = scale(hpv_copynumber),
		      tumor_volume = scale(tumor_volume),
		      tumor_size = scale(tumor_size),
		      age = scale(age),
		      ssGSEA_apoptosis_pos = scale(ssGSEA_apoptosis_pos),
		      ssGSEA_apoptosis_neg = scale(ssGSEA_apoptosis_neg),
		      ssGSEA_necrosis = scale(ssGSEA_necrosis),
		      ssGSEA_mitosis = scale(ssGSEA_mitosis)) %>%
        tidyr::drop_na() %>%
        as.data.frame()
	
smry_reads = summary(lm(formula = aligned_reads ~ ., data = data_ %>% dplyr::select(-mean_af)))
smry_af = summary(lm(formula = mean_af ~ ., data = data_ %>% dplyr::select(-aligned_reads, -hpv_copynumber)))

plot_ = smry_reads$coefficients %>%
	as.data.frame() %>%
	tibble::rownames_to_column("variable") %>%
	dplyr::select(variable, coefficient_reads = `Estimate`, p_value_reads = `Pr(>|t|)`) %>%
	dplyr::left_join(smry_af$coefficients %>%
			 as.data.frame() %>%
			 tibble::rownames_to_column("variable") %>%
			 dplyr::select(variable, coefficient_af = `Estimate`, p_value_af = `Pr(>|t|)`),
			 by = "variable") %>%
	dplyr::filter(variable != "(Intercept)") %>%
	dplyr::mutate(variable_cat = case_when(
		variable == "cfdna_concentration" ~ "Genomics",
		variable == "hpv_copynumber" ~ "Genomics",
		variable == "tumor_purity" ~ "Genomics",
		variable == "tumor_volume" ~ "Imaging",
		variable == "tumor_size" ~ "Clinical",
		variable == "age" ~ "Clinical",
		variable == "t_stageT2" ~ "Clinical",
		variable == "t_stageT2" ~ "Clinical",
		variable == "n_stageN2a" ~ "Clinical",
		variable == "n_stageN2b" ~ "Clinical",
		variable == "n_stageN2c" ~ "Clinical",
		variable == "smoking_statusyes" ~ "Clinical",
		variable == "hypoxianever_hypoxic" ~ "Imaging",
		variable == "hypoxiapersistent" ~ "Clinical",
		variable == "ssGSEA_apoptosis_pos" ~ "Transcriptomics",
		variable == "ssGSEA_apoptosis_neg" ~ "Transcriptomics",
		variable == "ssGSEA_necrosis" ~ "Transcriptomics",
		variable == "ssGSEA_mitosis" ~ "Transcriptomics"
	)) %>%
	dplyr::mutate(is_significant = p_value_reads<.1 | p_value_af<.1) %>%
	dplyr::mutate(variable_cat = factor(variable_cat, levels = rev(c("Genomics", "Transcriptomics", "Imaging", "Clinical")), ordered = TRUE)) %>%
	ggplot(aes(x = coefficient_reads, y = coefficient_af, shape = variable_cat, fill = is_significant)) +
	geom_abline(intercept = 0, slope = 1, color = "goldenrod3", size = 1, alpha = .85, linetype = 2) +
	geom_point(stat = "identity", color = "black", size = 3, alpha = .75) +
	scale_fill_manual(values = c("#ffffff", "#e41a1c")) +
	scale_shape_manual(values = c("Genomics" = 21, "Transcriptomics" = 22, "Imaging" = 23, "Clinical" = 24)) +
	scale_x_continuous(limits = c(-1, 2)) +
	scale_y_continuous(limits = c(-1, 2)) +
	xlab("Coefficients variables\nregressed number of reads") +
	ylab("Coefficients variables\nregressed Mean AF") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 10)),
	      axis.title.y = element_text(margin = margin(r = 10)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = FALSE,
	       shape = guide_legend(title = "Variable\nCategory", order = 1, override.aes = list(size = 3)))

pdf(file = "../res/Compare_Coefficients_Mean_AF_Number_Read_Pairs.pdf", width = 4.35, height = 3)
print(plot_)
dev.off()

data_ = idx_metrics_ft %>%
	dplyr::left_join(readr::read_tsv(file = "../res/GSEA.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert(),
			 by = "patient_id_mskcc") %>%
	dplyr::mutate(mean_af = mean_af*100,
		      aligned_reads = log10(aligned_reads)) %>%
        dplyr::select(aligned_reads,
		      mean_af,
		      cfdna_concentration = concentration_ng_uL,
		      hpv_copynumber = hpv_panel_copy_number,
		      tumor_purity = purity_wes,
		      tumor_volume = plan_volume,
		      tumor_size = primary_tumor_size_cm,
		      age = age,
		      t_stage = t_stage,
		      n_stage = n_stage,
		      smoking_status = smoking_category_yes_never,
		      hypoxia = simplified_hypoxia_group,
		      ssGSEA_apoptosis_pos = GOBP_REGULATION_OF_MITOCHONDRIAL_MEMBRANE_PERMEABILITY_INVOLVED_IN_APOPTOTIC_PROCESS,
		      ssGSEA_apoptosis_neg = GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS,
		      ssGSEA_necrosis = REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS,
		      ssGSEA_mitosis = GOBP_MITOTIC_G1_S_TRANSITION_CHECKPOINT_SIGNALING) %>%
	readr::type_convert() %>%
	dplyr::mutate(mean_af = scale(mean_af),
		      cfdna_concentration = scale(cfdna_concentration),
		      hpv_copynumber = scale(hpv_copynumber),
		      tumor_volume = scale(tumor_volume),
		      tumor_size = scale(tumor_size),
		      age = scale(age),
		      ssGSEA_apoptosis_pos = scale(ssGSEA_apoptosis_pos),
		      ssGSEA_apoptosis_neg = scale(ssGSEA_apoptosis_neg),
		      ssGSEA_necrosis = scale(ssGSEA_necrosis),
		      ssGSEA_mitosis = scale(ssGSEA_mitosis)) %>%
        tidyr::drop_na() %>%
        as.data.frame()

fit_ = lm(formula = aligned_reads ~ ., data = data_ %>% dplyr::select(-mean_af))
res_ = summary(fit_)
aic_ = MASS::stepAIC(object = lm(formula = aligned_reads ~ ., data = data_ %>% dplyr::select(-mean_af)), direction = "backward")

plot_ = res_ %>%
	.[["coefficients"]] %>%
	as.data.frame() %>%
	tibble::rownames_to_column("variable") %>%
	dplyr::as_tibble() %>%
	dplyr::filter(variable != "(Intercept)") %>%
	dplyr::mutate(is_significant = ifelse(`Pr(>|t|)`<.1, "Yes", "No")) %>%
	dplyr::mutate(variable_cat = case_when(
		variable == "cfdna_concentration" ~ "Genomics",
		variable == "hpv_copynumber" ~ "Genomics",
		variable == "tumor_purity" ~ "Genomics",
		variable == "tumor_volume" ~ "Imaging",
		variable == "tumor_size" ~ "Clinical",
		variable == "age" ~ "Clinical",
		variable == "t_stageT2" ~ "Clinical",
		variable == "t_stageT2" ~ "Clinical",
		variable == "n_stageN2a" ~ "Clinical",
		variable == "n_stageN2b" ~ "Clinical",
		variable == "n_stageN2c" ~ "Clinical",
		variable == "smoking_statusyes" ~ "Clinical",
		variable == "hypoxianever_hypoxic" ~ "Imaging",
		variable == "hypoxiapersistent" ~ "Imaging",
		variable == "ssGSEA_apoptosis_pos" ~ "Transcriptomics",
		variable == "ssGSEA_apoptosis_neg" ~ "Transcriptomics",
		variable == "ssGSEA_necrosis" ~ "Transcriptomics",
		variable == "ssGSEA_mitosis" ~ "Transcriptomics"
	)) %>%
	dplyr::mutate(variable_cat = factor(variable_cat, levels = rev(c("Genomics", "Transcriptomics", "Imaging", "Clinical")), ordered = TRUE)) %>%
	dplyr::arrange(variable_cat, Estimate) %>%
	dplyr::mutate(variable = factor(variable, levels = variable, ordered = TRUE)) %>%
	ggplot(aes(x = variable, ymin = 0, ymax = Estimate, y = Estimate, fill = is_significant, size = -log10(`Pr(>|t|)`), shape = variable_cat)) +
	geom_linerange(stat = "identity", size = .5) +
	geom_hline(yintercept = 0, size = 1) +
	geom_point(stat = "identity") +
	xlab("") +
	ylab("") +
	scale_fill_manual(values = c("#bdbdbd", "#e41a1c")) +
	scale_shape_manual(values = c("Genomics" = 21, "Transcriptomics" = 22, "Imaging" = 23, "Clinical" = 24)) +
	scale_x_discrete() +
	scale_y_continuous() +
	coord_flip() +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = FALSE,
	       shape = guide_legend(title = "Variable Category", order = 1, override.aes = list(size = 3)),
	       size = guide_legend(title = expression(-Log[10]~"p-value")))

pdf(file = "../res/Linear_Regression_Coefficients_Maxi.pdf", width = 6, height = 4)
print(plot_)
dev.off()

variables_filtered = c("All", "n_stage", "hypoxia",
		       "smoking_status", "ssGSEA_mitosis", "t_stage",
		       "cfdna_concentration", "ssGSEA_apoptosis_pos", "ssGSEA_apoptosis_neg",
		       "tumor_purity", "tumor_size", "hpv_copynumber", "ssGSEA_necrosis",
		       "age")
extracted_aic = vector(mode = "numeric", length = length(variables_filtered))
names(extracted_aic) = variables_filtered

extracted_aic["All"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -n_stage)
extracted_aic["n_stage"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -hypoxia)
extracted_aic["hypoxia"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -smoking_status)
extracted_aic["smoking_status"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -ssGSEA_mitosis)
extracted_aic["ssGSEA_mitosis"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -t_stage)
extracted_aic["t_stage"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -cfdna_concentration)
extracted_aic["cfdna_concentration"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -ssGSEA_apoptosis_pos)
extracted_aic["ssGSEA_apoptosis_pos"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -ssGSEA_apoptosis_neg)
extracted_aic["ssGSEA_apoptosis_neg"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -tumor_purity)
extracted_aic["tumor_purity"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -tumor_size)
extracted_aic["tumor_size"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -hpv_copynumber)
extracted_aic["hpv_copynumber"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -ssGSEA_necrosis)
extracted_aic["ssGSEA_necrosis"] = extractAIC(fit_)[2]

fit_ = update(fit_, . ~ . -age)
extracted_aic["age"] = extractAIC(fit_)[2]

plot_ = dplyr::tibble(variables = variables_filtered,
		      aic = extracted_aic) %>%
	dplyr::mutate(index = 1:n()) %>%
	ggplot(aes(x = index, y = aic)) +
	geom_line(stat = "identity", size = .95) +
	geom_point(stat = "identity", color = "black", fill = "white", shape = 21, size = 3) +
	scale_x_continuous(breaks = 1:length(variables_filtered),
			   labels = paste0("- ", variables_filtered, " ", length(variables_filtered):1)) +
	scale_y_continuous() +
	xlab("") +
	ylab("AIC") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust=1),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Linear_Regression_Stepwise_AIC_Maxi.pdf", width = 3.35, height = 4)
print(plot_)
dev.off()

data_ = idx_metrics_ft %>%
	dplyr::left_join(readr::read_tsv(file = "../res/GSEA.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert(),
			 by = "patient_id_mskcc") %>%
	dplyr::mutate(mean_af = mean_af*100,
		      aligned_reads = log10(aligned_reads)) %>%
        dplyr::select(aligned_reads,
		      mean_af,
		      cfdna_concentration = concentration_ng_uL,
		      hpv_copynumber = hpv_panel_copy_number,
		      tumor_purity = purity_wes,
		      tumor_volume = plan_volume,
		      tumor_size = primary_tumor_size_cm,
		      age = age,
		      t_stage = t_stage,
		      n_stage = n_stage,
		      smoking_status = smoking_category_yes_never,
		      hypoxia = simplified_hypoxia_group,
		      ssGSEA_apoptosis_pos = GOBP_REGULATION_OF_MITOCHONDRIAL_MEMBRANE_PERMEABILITY_INVOLVED_IN_APOPTOTIC_PROCESS,
		      ssGSEA_apoptosis_neg = GOBP_INFLAMMATORY_CELL_APOPTOTIC_PROCESS,
		      ssGSEA_necrosis = REACTOME_RIPK1_MEDIATED_REGULATED_NECROSIS,
		      ssGSEA_mitosis = GOBP_MITOTIC_G1_S_TRANSITION_CHECKPOINT_SIGNALING) %>%
	readr::type_convert() %>%
	tidyr::drop_na() %>%
        as.data.frame()

#########################################################
# Tumor volume
#########################################################
plot_ = data_ %>%
	ggplot(aes(x = tumor_volume, y = aligned_reads)) +
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

pdf(file = "../res/Total_Reads_Tumor_Volume.pdf", width = 3.35, height = 3)
print(plot_)
dev.off()

#########################################################
# Age
#########################################################
plot_ = data_ %>%
	dplyr::mutate(age_cat = case_when(
		age < quantile(age, 1/3, na.rm = TRUE) ~ paste0("<", round(quantile(age, 1/3, na.rm = TRUE))),
		age >= quantile(age, 1/3, na.rm = TRUE) & age < quantile(age, 2/3, na.rm = TRUE) ~ paste0(round(quantile(age, 1/3, na.rm = TRUE)), "-", round(quantile(age, 2/3, na.rm = TRUE))),
		age > quantile(age, 2/3, na.rm = TRUE) ~ paste0(">", round(quantile(age, 2/3, na.rm = TRUE))),
		TRUE ~ ""
	)) %>%
	dplyr::filter(age_cat != "") %>%
	dplyr::mutate(age_cat = factor(age_cat, levels = c("<56", "56-60", ">60"), ordered = TRUE)) %>%
	ggplot(aes(x = age_cat, y = aligned_reads)) +
	geom_violin(stat = "ydensity", draw_quantiles = c(.25, .5, .75), trim = FALSE, scale = "area") +
	scale_x_discrete() +
	scale_y_continuous(limits = c(1, 11),
			   breaks = c(2, 3, 4, 5, 6, 7),
			   labels = scientific_10(10^c(2, 3, 4, 5, 6, 7))) +
	xlab("Age Tertiles (yrs)") +
	ylab("cfDNA Aligned HPV Read Pairs") +
	geom_signif(stat = "signif",
		    comparisons = list(c(1, 2),
				       c(1, 3),
				       c(2, 3)),
		    test = "wilcox.test",
		    test.args = list(alternative = "greater"),
		    y_position = c(9, 10, 11)-.5) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Total_Reads_Age.pdf", width = 3.35, height = 3.25)
print(plot_)
dev.off()
