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

idx_metrics = readr::read_tsv(file = url_idx_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::select(sample_name = SAMPLE_NAME, chromosome = CHROMOSOME, aligned_reads = ALIGNED_READS) %>%
	      dplyr::filter(chromosome %in% target_contigs)

idx_metrics_ft = readr::read_tsv(file = url_idx_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      	 readr::type_convert() %>%
		 dplyr::select(sample_name = SAMPLE_NAME, chromosome = CHROMOSOME, aligned_reads = ALIGNED_READS, fragment_length = FRAGMENT_LENGTH) %>%
	      	 dplyr::filter(chromosome %in% target_contigs)

insert_size_metrics = readr::read_tsv(file = url_insert_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert()

insert_size_metrics_ft = readr::read_tsv(file = url_insert_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      	 readr::type_convert()

mutation_smry = readr::read_tsv(file = url_mutation_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
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

plot_ = idx_metrics %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::group_by(sample_name, Is_ctDNA) %>%
	dplyr::summarize(aligned_reads = sum(aligned_reads)) %>%
	dplyr::ungroup() %>%
	dplyr::left_join(insert_size_metrics %>%
			 dplyr::rename(sample_name = SAMPLE_NAME), by = "sample_name") %>%
	ggplot(aes(x = Is_ctDNA, y = aligned_reads+1, size = MEAN_INSERT_SIZE)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous() +
	scale_x_discrete() +
	scale_y_log10(labels = scientific_10) +
	xlab("ctDNA") +
	ylab("Number of Read Pairs") +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = log10(1E7),
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Mean Insert\nSize (bp)"))

pdf(file = "../res/Number_of_Read_Pairs_by_Insert_Size_no_Filtering.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = idx_metrics %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = case_when(
		is.na(Is_ctDNA) ~ "?",
		TRUE ~ Is_ctDNA
	)) %>%
	dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve", "?"), ordered = TRUE)) %>%
	dplyr::group_by(sample_name, Is_ctDNA) %>%
	dplyr::summarize(aligned_reads = sum(aligned_reads)) %>%
	dplyr::ungroup() %>%
	dplyr::left_join(insert_size_metrics %>%
			 dplyr::rename(sample_name = SAMPLE_NAME), by = "sample_name") %>%
	ggplot(aes(x = Is_ctDNA, y = aligned_reads+1, size = MEAN_INSERT_SIZE)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous() +
	scale_x_discrete() +
	scale_y_log10(labels = scientific_10) +
	xlab("ctDNA") +
	ylab("Number of Read Pairs") +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = log10(1E7),
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Mean Insert\nSize (bp)"))

pdf(file = "../res/Number_of_Read_Pairs_by_Insert_Size_no_Filtering_with_Unknown.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = idx_metrics_ft %>%
	dplyr::filter(fragment_length == 0) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::group_by(sample_name, Is_ctDNA) %>%
	dplyr::summarize(aligned_reads = sum(aligned_reads)) %>%
	dplyr::ungroup() %>%
	dplyr::left_join(insert_size_metrics_ft %>%
			 dplyr::filter(FRAGMENT_LENGTH == 0) %>%
			 dplyr::rename(sample_name = SAMPLE_NAME), by = "sample_name") %>%
	ggplot(aes(x = Is_ctDNA, y = aligned_reads+1, size = MEAN_INSERT_SIZE)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous() +
	scale_x_discrete() +
	scale_y_log10(labels = scientific_10) +
	xlab("ctDNA") +
	ylab("Number of Read Pairs") +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = log10(1E7),
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Mean Insert\nSize (bp)"))

pdf(file = "../res/Number_of_Read_Pairs_by_Insert_Size_Primer_Filtered.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = idx_metrics_ft %>%
	dplyr::filter(fragment_length == 0) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = case_when(
		is.na(Is_ctDNA) ~ "?",
		TRUE ~ Is_ctDNA
	)) %>%
	dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve", "?"), ordered = TRUE)) %>%
	dplyr::group_by(sample_name, Is_ctDNA) %>%
	dplyr::summarize(aligned_reads = sum(aligned_reads)) %>%
	dplyr::ungroup() %>%
	dplyr::left_join(insert_size_metrics_ft %>%
			 dplyr::filter(FRAGMENT_LENGTH == 0) %>%
			 dplyr::rename(sample_name = SAMPLE_NAME), by = "sample_name") %>%
	ggplot(aes(x = Is_ctDNA, y = aligned_reads+1, size = MEAN_INSERT_SIZE)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous() +
	scale_x_discrete() +
	scale_y_log10(labels = scientific_10) +
	xlab("ctDNA") +
	ylab("Number of Read Pairs") +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = log10(1E7),
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Mean Insert\nSize (bp)"))

pdf(file = "../res/Number_of_Read_Pairs_by_Insert_Size_Primer_Filtered_with_Unknown.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = idx_metrics_ft %>%
	dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::group_by(sample_name, Is_ctDNA) %>%
	dplyr::summarize(aligned_reads = sum(aligned_reads)) %>%
	dplyr::ungroup() %>%
	dplyr::left_join(insert_size_metrics_ft %>%
			 dplyr::filter(FRAGMENT_LENGTH == FRAGMENT_LENGTH_THRESHOLD) %>%
			 dplyr::rename(sample_name = SAMPLE_NAME), by = "sample_name") %>%
	ggplot(aes(x = Is_ctDNA, y = aligned_reads+1, size = MEAN_INSERT_SIZE)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous() +
	scale_x_discrete() +
	scale_y_log10(labels = scientific_10) +
	xlab("ctDNA") +
	ylab("Number of Read Pairs") +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = log10(1E7),
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Mean Insert\nSize (bp)"))

pdf(file = "../res/Number_of_Read_Pairs_by_Insert_Size_Primer_Filtered_Fragment_Filtered.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = idx_metrics_ft %>%
	dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = case_when(
		is.na(Is_ctDNA) ~ "?",
		TRUE ~ Is_ctDNA
	)) %>%
	dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve", "?"), ordered = TRUE)) %>%
	dplyr::group_by(sample_name, Is_ctDNA) %>%
	dplyr::summarize(aligned_reads = sum(aligned_reads)) %>%
	dplyr::ungroup() %>%
	dplyr::left_join(insert_size_metrics_ft %>%
			 dplyr::filter(FRAGMENT_LENGTH == FRAGMENT_LENGTH_THRESHOLD) %>%
			 dplyr::rename(sample_name = SAMPLE_NAME), by = "sample_name") %>%
	ggplot(aes(x = Is_ctDNA, y = aligned_reads+1, size = MEAN_INSERT_SIZE)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous() +
	scale_x_discrete() +
	scale_y_log10(labels = scientific_10) +
	xlab("ctDNA") +
	ylab("Number of Read Pairs") +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = log10(1E7),
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Mean Insert\nSize (bp)"))

pdf(file = "../res/Number_of_Read_Pairs_by_Insert_Size_Primer_Filtered_Fragment_Filtered_with_Unknown.pdf", width = 5, height = 5)
print(plot_)
dev.off()


fit = idx_metrics_ft %>%
      dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD) %>%
      dplyr::left_join(mutation_smry %>%
		       dplyr::filter(FILTER == "PASS") %>%
		       dplyr::group_by(Tumor_Sample_Barcode) %>%
		       dplyr::summarize(mean_af = mean(t_maf),
					median_af = median(t_maf),
					max_af = max(t_maf)) %>%
		       dplyr::ungroup() %>%
		       dplyr::rename(sample_name = Tumor_Sample_Barcode), by = "sample_name") %>%
       dplyr::left_join(manifest, by = "sample_name") %>%
       dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
       dplyr::mutate(chromosome = case_when(
		chromosome == "J04353.1" ~ "HPV-31",
		chromosome == "M12732.1" ~ "HPV-33",
	   	chromosome == "NC001357.1" ~ "HPV-18",
		chromosome == "NC001526.4" ~ "HPV-16",
		chromosome == "X74477.1" ~ "HPV-35"
	)) %>%
	dplyr::left_join(smry_ft %>%
			 dplyr::select(sample_name, Is_ctDNA), by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = ifelse(is.na(Is_ctDNA), "?", Is_ctDNA)) %>%
	dplyr::filter(Is_ctDNA != "?") %>%
	dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve"), ordered = TRUE)) %>%
	dplyr::filter(chromosome == "HPV-16") %>%
	dplyr::mutate(Is_ctDNA = case_when(
		  hpv_type_wes_wgs == "HPV-16" & Is_ctDNA == "+ve" ~ "+ve",
		  TRUE ~ "-ve"
	 )) %>%
	rpart::rpart(formula = Is_ctDNA ~ aligned_reads, data = ., method = "class")

smry = dplyr::tibble(chromosome = rep("HPV-16", length(unique(smry_ft$hpv_type_wes_wgs))),
		     hpv_type_wes_wgs = unique(smry_ft$hpv_type_wes_wgs),
		     xintercept = fit$splits[1,"index"])

plot_ = idx_metrics_ft %>%
	dplyr::left_join(mutation_smry %>%
	 		 dplyr::filter(FILTER == "PASS") %>%
			 dplyr::group_by(Tumor_Sample_Barcode) %>%
			 dplyr::summarize(mean_af = mean(t_maf),
					  median_af = median(t_maf),
				   	  max_af = max(t_maf)) %>%
			 dplyr::ungroup() %>%
			 dplyr::rename(sample_name = Tumor_Sample_Barcode), by = "sample_name") %>%
	dplyr::left_join(manifest, by = "sample_name") %>%
	dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
	dplyr::mutate(chromosome = case_when(
		chromosome == "J04353.1" ~ "HPV-31",
		chromosome == "M12732.1" ~ "HPV-33",
	   	chromosome == "NC001357.1" ~ "HPV-18",
		chromosome == "NC001526.4" ~ "HPV-16",
		chromosome == "X74477.1" ~ "HPV-35"
	)) %>%
	dplyr::left_join(smry_ft %>%
			 dplyr::select(sample_name, Is_ctDNA), by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = ifelse(is.na(Is_ctDNA), "?", Is_ctDNA)) %>%
	dplyr::filter(Is_ctDNA != "?") %>%
	dplyr::mutate(mean_af = ifelse(mean_af == 0, .000001, mean_af)) %>%
	dplyr::mutate(median_af = ifelse(median_af == 0, .000001, median_af)) %>%
	dplyr::mutate(max_af = ifelse(max_af == 0, .000001, max_af)) %>%
	dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD) %>%
	dplyr::arrange(desc(Is_ctDNA)) %>%
	ggplot(aes(x = aligned_reads+1, y = 100*mean_af, color = Is_ctDNA, shape = hpv_type_wes_wgs)) +
	geom_vline(data = smry, aes(xintercept = xintercept), alpha = 1, linetype = 3) +
	geom_point(stat = "identity", fill = NA, alpha = 1, size = 1.5) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_shape_manual(values = c(1, 2, 3, 4, 5, 6)) +
	xlab("Number of Read Pairs Aligned") +
	ylab("Mean Allele Fraction (%)") +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(limits = c(1e-4, 100),
		      breaks = c(1e-4, 1e-3, 1e-1, 1, 10, 100),
		      labels = c("ND", ".001", ".01", "1", "10", "100")) +
	facet_grid(hpv_type_wes_wgs ~ chromosome) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(shape = guide_legend(title = "Tumor\nHPV Type"),
	       color = guide_legend(title = "ctDNA", override.aes = list(alpha = 1)))

pdf(file = "../res/Number_Read_Pairs_Aligned_by_Mean_AF.pdf", width = 12/1.15, height = 10/1.15)
print(plot_)
dev.off()
