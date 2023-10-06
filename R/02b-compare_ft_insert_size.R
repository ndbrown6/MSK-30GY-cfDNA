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

plot_ = insert_size_metrics %>%
	dplyr::select(sample_name = SAMPLE_NAME, mean_insert_size = MEAN_INSERT_SIZE, read_pairs = READ_PAIRS) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::arrange(desc(read_pairs)) %>%
	ggplot(aes(x = factor(Is_ctDNA), y = mean_insert_size, size = read_pairs)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous(breaks = c(1E3, 1E4, 1E5, 2E5, 3E5),
			      label = scientific_10) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 175)) +
	xlab("ctDNA") +
	ylab("Mean Insert Size (bp)") +
	geom_rug(data = primer_set, mapping = aes(y = insert_size), sides = "r", size= .2, inherit.aes = FALSE) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Number of\nRead Pairs"))
	
pdf(file = "../res/Mean_Insert_Size_no_Filtering.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_metrics %>%
	dplyr::select(sample_name = SAMPLE_NAME, mean_insert_size = MEAN_INSERT_SIZE, read_pairs = READ_PAIRS) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = case_when(
		is.na(Is_ctDNA) ~ "?",
		TRUE ~ Is_ctDNA
	)) %>%
	dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve", "?"), ordered = TRUE)) %>%
	dplyr::arrange(desc(read_pairs)) %>%
	ggplot(aes(x = factor(Is_ctDNA), y = mean_insert_size, size = read_pairs)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous(breaks = c(1E3, 1E4, 1E5, 2E5, 3E5),
			      label = scientific_10) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 175)) +
	xlab("ctDNA") +
	ylab("Mean Insert Size (bp)") +
	geom_rug(data = primer_set, mapping = aes(y = insert_size), sides = "r", size= .2, inherit.aes = FALSE) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Number of\nRead Pairs"))
	
pdf(file = "../res/Mean_Insert_Size_no_Filtering_with_Unknown.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_metrics_ft %>%
	dplyr::filter(FRAGMENT_LENGTH == 0) %>%
	dplyr::select(sample_name = SAMPLE_NAME, mean_insert_size = MEAN_INSERT_SIZE, read_pairs = READ_PAIRS) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::arrange(desc(read_pairs)) %>%
	ggplot(aes(x = factor(Is_ctDNA), y = mean_insert_size, size = read_pairs)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous(breaks = c(1E3, 1E4, 1E5, 2E5, 3E5),
			      label = scientific_10) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 175)) +
	xlab("ctDNA") +
	ylab("Mean Insert Size (bp)") +
	geom_rug(data = primer_set, mapping = aes(y = insert_size), sides = "r", size= .2, inherit.aes = FALSE) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Number of\nRead Pairs"))
	
pdf(file = "../res/Mean_Insert_Size_Primer_Filtered.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_metrics_ft %>%
	dplyr::filter(FRAGMENT_LENGTH == 0) %>%
	dplyr::select(sample_name = SAMPLE_NAME, mean_insert_size = MEAN_INSERT_SIZE, read_pairs = READ_PAIRS) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = case_when(
		is.na(Is_ctDNA) ~ "?",
		TRUE ~ Is_ctDNA
	)) %>%
	dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve", "?"), ordered = TRUE)) %>%
	dplyr::arrange(desc(read_pairs)) %>%
	ggplot(aes(x = factor(Is_ctDNA), y = mean_insert_size, size = read_pairs)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous(breaks = c(1E3, 1E4, 1E5, 2E5, 3E5),
			      label = scientific_10) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 175)) +
	xlab("ctDNA") +
	ylab("Mean Insert Size (bp)") +
	geom_rug(data = primer_set, mapping = aes(y = insert_size), sides = "r", size= .2, inherit.aes = FALSE) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Number of\nRead Pairs"))
	
pdf(file = "../res/Mean_Insert_Size_Primer_Filtered_with_Unknown.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_smry %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::mutate(Is_ctDNA = case_when(
		Is_ctDNA == "+ve" ~ "ctDNA +ve",
		Is_ctDNA == "-ve" ~ "ctDNA -ve"
	)) %>%
	ggplot(aes(x = INSERT_SIZE, y = READ_COUNT, group = sample_name)) +
	geom_step(stat = "identity", size = .5, alpha = .25, color = "#cb181d") +
	scale_x_sqrt(limits = c(5, NA),
		     breaks = c(10, 15, 25, 50, 100, 200, 300, 400, 500),
		     labels = c(10, 15, 25, 50, 100, 200, 300, 400, 500)) +
	scale_y_continuous(limits = c(0, NA),
		     	   labels = scientific_10) +
	xlab("Insert size (bp)") +
	ylab("Number of Read Pairs") +
	facet_wrap(~Is_ctDNA, ncol = 1, scales = "free_y") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_rug(data = primer_set, mapping = aes(x = insert_size), size = .2, , sides = "b", inherit.aes = FALSE)

pdf(file = "../res/Insert_Size_Distribution_no_Filtering.pdf", width = 6, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_smry_ft %>%
	dplyr::filter(FRAGMENT_LENGTH == 0) %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::mutate(Is_ctDNA = case_when(
		Is_ctDNA == "+ve" ~ "ctDNA +ve",
		Is_ctDNA == "-ve" ~ "ctDNA -ve"
	)) %>%
	ggplot(aes(x = INSERT_SIZE, y = READ_COUNT, group = sample_name)) +
	geom_step(stat = "identity", size = .5, alpha = .25, color = "#cb181d") +
	scale_x_sqrt(limits = c(5, NA),
		     breaks = c(10, 15, 25, 50, 100, 200, 300, 400, 500),
		     labels = c(10, 15, 25, 50, 100, 200, 300, 400, 500)) +
	scale_y_continuous(limits = c(0, NA),
		     	   labels = scientific_10) +
	xlab("Insert size (bp)") +
	ylab("Number of Read Pairs") +
	facet_wrap(~Is_ctDNA, ncol = 1, scales = "free_y") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_rug(data = primer_set, mapping = aes(x = insert_size), size = .2, , sides = "b", inherit.aes = FALSE)

pdf(file = "../res/Insert_Size_Distribution_Primer_Filtered.pdf", width = 6, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_smry %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::mutate(Is_ctDNA = case_when(
		Is_ctDNA == "+ve" ~ "ctDNA +ve",
		Is_ctDNA == "-ve" ~ "ctDNA -ve"
	)) %>%
	ggplot(aes(x = INSERT_SIZE, y = 100*`%_READS`, group = sample_name)) +
	geom_step(stat = "identity", size = .5, alpha = .25, color = "#cb181d") +
	scale_x_sqrt(limits = c(5, NA),
		     breaks = c(10, 15, 25, 50, 100, 200, 300, 400, 500),
		     labels = c(10, 15, 25, 50, 100, 200, 300, 400, 500)) +
	scale_y_continuous(limits = c(0, NA)) +
	xlab("Insert size (bp)") +
	ylab("Fraction of Read Pairs (%)") +
	facet_wrap(~Is_ctDNA, ncol = 1, scales = "free_y") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_rug(data = primer_set, mapping = aes(x = insert_size), size = .2, , sides = "b", inherit.aes = FALSE)

pdf(file = "../res/Insert_Size_Distribution_Fraction_no_Filtering.pdf", width = 6, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_smry %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::filter(INSERT_SIZE == 25) %>%
	ggplot(aes(x = Is_ctDNA, y = 100*`%_READS`)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75, size = 3) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 100)) +
	xlab("ctDNA") +
	ylab("Fraction of Total Reads (%)") +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 100,
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)))

pdf(file = "../res/Insert_Size_25bp_Fraction_no_Filtering.pdf", width = 4, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_smry_ft %>%
	dplyr::filter(FRAGMENT_LENGTH == 0) %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::mutate(Is_ctDNA = case_when(
		Is_ctDNA == "+ve" ~ "ctDNA +ve",
		Is_ctDNA == "-ve" ~ "ctDNA -ve"
	)) %>%
	ggplot(aes(x = INSERT_SIZE, y = 100*`%_READS`, group = sample_name)) +
	geom_step(stat = "identity", size = .5, alpha = .25, color = "#cb181d") +
	scale_x_sqrt(limits = c(5, NA),
		     breaks = c(10, 15, 25, 50, 100, 200, 300, 400, 500),
		     labels = c(10, 15, 25, 50, 100, 200, 300, 400, 500)) +
	scale_y_continuous(limits = c(0, NA)) +
	xlab("Insert size (bp)") +
	ylab("Fraction of Read Pairs (%)") +
	facet_wrap(~Is_ctDNA, ncol = 1, scales = "free_y") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_rug(data = primer_set, mapping = aes(x = insert_size), size = .2, , sides = "b", inherit.aes = FALSE)

pdf(file = "../res/Insert_Size_Distribution_Fraction_Primer_Filtered.pdf", width = 6, height = 5)
print(plot_)
dev.off()

smry_ = insert_size_smry %>%
	reshape2::dcast(SAMPLE_NAME ~ INSERT_SIZE, value.var = "%_READS", fill = 0) %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft %>% dplyr::select(sample_name, Is_ctDNA), by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::select(-sample_name) %>%
	tidyr::gather(key = variable, value = value, -Is_ctDNA) %>% 
	dplyr::group_by(Is_ctDNA, variable) %>% 
  	dplyr::summarise(value = list(value)) %>% 
	tidyr::spread(Is_ctDNA, value) %>%
	dplyr::group_by(variable) %>% 
  	dplyr::mutate(p_value = t.test(unlist(`+ve`), unlist(`-ve`))$p.value,
		      t_value = t.test(unlist(`+ve`), unlist(`-ve`))$statistic) %>%
	dplyr::ungroup() %>%
	readr::type_convert() %>%
	dplyr::mutate(Filtering = "None") %>%
	dplyr::bind_rows(
		insert_size_smry_ft %>%
		dplyr::filter(FRAGMENT_LENGTH == 0) %>%
		reshape2::dcast(SAMPLE_NAME ~ INSERT_SIZE, value.var = "%_READS", fill = 0) %>%
		dplyr::rename(sample_name = SAMPLE_NAME) %>%
		dplyr::left_join(smry_ft %>% dplyr::select(sample_name, Is_ctDNA), by = "sample_name") %>%
		dplyr::filter(!is.na(Is_ctDNA)) %>%
		dplyr::select(-sample_name) %>%
		tidyr::gather(key = variable, value = value, -Is_ctDNA) %>% 
		dplyr::group_by(Is_ctDNA, variable) %>% 
		dplyr::summarise(value = list(value)) %>% 
		tidyr::spread(Is_ctDNA, value) %>%
		dplyr::group_by(variable) %>% 
		dplyr::mutate(p_value = t.test(unlist(`+ve`), unlist(`-ve`))$p.value,
			      t_value = t.test(unlist(`+ve`), unlist(`-ve`))$statistic) %>%
		dplyr::ungroup() %>%
		readr::type_convert() %>%
		dplyr::mutate(Filtering = "Primer\nFiltered"))

plot_ = smry_ %>%
	ggplot(aes(x = variable, y = t_value, color = Filtering)) +
	geom_hline(yintercept = 0, linetype = 3, color = "black", size = 1) +
	geom_step(stat = "identity", alpha = .75, size = 1) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_continuous() +
	scale_y_continuous(labels = scientific_10) +
	xlab("Insert size (bp)") +
	ylab(expression(tau~"statistic")) +
	theme(axis.text = element_text(size = 9),
	      axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20)),
	      plot.title = element_text(hjust = 0.5),
	      axis.ticks.x = element_blank(),
	      axis.ticks.y = element_blank(),
	      panel.background = element_rect(fill = "white", colour = "white", size = 2, linetype = "solid"),
  	      panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "#ebebebff"), 
  	      panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "#ebebebff"),
	      legend.key = element_rect(colour = NA, fill = NA)) +
	geom_rug(data = primer_set, mapping = aes(x = insert_size), size = .2, , sides = "b", inherit.aes = FALSE) +
	facet_zoom(xlim = c(19, 50))

pdf(file = "../res/Insert_Size_Distribution_by_ctDNA.pdf", width = 9, height = 9)
print(plot_)
dev.off()

plot_ = insert_size_metrics_ft %>%
	dplyr::filter(FRAGMENT_LENGTH == FRAGMENT_LENGTH_THRESHOLD) %>%
	dplyr::select(sample_name = SAMPLE_NAME, mean_insert_size = MEAN_INSERT_SIZE, read_pairs = READ_PAIRS) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::arrange(desc(read_pairs)) %>%
	ggplot(aes(x = factor(Is_ctDNA), y = mean_insert_size, size = read_pairs)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous(breaks = c(1E3, 1E4, 1E5, 2E5, 3E5),
			      label = scientific_10) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 175)) +
	xlab("ctDNA") +
	ylab("Mean Insert Size (bp)") +
	geom_rug(data = primer_set, mapping = aes(y = insert_size), sides = "r", size= .2, inherit.aes = FALSE) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Number of\nRead Pairs"))
	
pdf(file = "../res/Mean_Insert_Size_Primer_Filtered_Fragment_Filtered.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_metrics_ft %>%
	dplyr::filter(FRAGMENT_LENGTH == FRAGMENT_LENGTH_THRESHOLD) %>%
	dplyr::select(sample_name = SAMPLE_NAME, mean_insert_size = MEAN_INSERT_SIZE, read_pairs = READ_PAIRS) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = case_when(
		is.na(Is_ctDNA) ~ "?",
		TRUE ~ Is_ctDNA
	)) %>%
	dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve", "?"), ordered = TRUE)) %>%
	dplyr::arrange(desc(read_pairs)) %>%
	ggplot(aes(x = factor(Is_ctDNA), y = mean_insert_size, size = read_pairs)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(stat = "identity", fill = "salmon", width = .1, height = 0, shape = 21, alpha = .75) +
	scale_size_continuous(breaks = c(1E3, 1E4, 1E5, 2E5, 3E5),
			      label = scientific_10) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 175)) +
	xlab("ctDNA") +
	ylab("Mean Insert Size (bp)") +
	geom_rug(data = primer_set, mapping = aes(y = insert_size), sides = "r", size= .2, inherit.aes = FALSE) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    y_position = 175,
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
 	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(size = guide_legend(title = "Number of\nRead Pairs"))
	
pdf(file = "../res/Mean_Insert_Size_Primer_Filtered_Fragment_Filtered_with_Unknown.pdf", width = 5, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_smry_ft %>%
	dplyr::filter(FRAGMENT_LENGTH == FRAGMENT_LENGTH_THRESHOLD) %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::mutate(Is_ctDNA = case_when(
		Is_ctDNA == "+ve" ~ "ctDNA +ve",
		Is_ctDNA == "-ve" ~ "ctDNA -ve"
	)) %>%
	ggplot(aes(x = INSERT_SIZE, y = READ_COUNT, group = sample_name)) +
	geom_step(stat = "identity", size = .5, alpha = .25, color = "#cb181d") +
	scale_x_sqrt(limits = c(5, NA),
		     breaks = c(10, 15, 25, 50, 100, 200, 300, 400, 500),
		     labels = c(10, 15, 25, 50, 100, 200, 300, 400, 500)) +
	scale_y_continuous(limits = c(0, NA),
		     	   labels = scientific_10) +
	xlab("Insert size (bp)") +
	ylab("Number of Read Pairs") +
	facet_wrap(~Is_ctDNA, ncol = 1, scales = "free_y") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_rug(data = primer_set, mapping = aes(x = insert_size), size = .2, , sides = "b", inherit.aes = FALSE)

pdf(file = "../res/Insert_Size_Distribution_Primer_Filtered_Fragment_Filtered.pdf", width = 6, height = 5)
print(plot_)
dev.off()

plot_ = insert_size_smry_ft %>%
	dplyr::filter(FRAGMENT_LENGTH == FRAGMENT_LENGTH_THRESHOLD) %>%
	dplyr::rename(sample_name = SAMPLE_NAME) %>%
	dplyr::left_join(smry_ft, by = "sample_name") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::mutate(Is_ctDNA = case_when(
		Is_ctDNA == "+ve" ~ "ctDNA +ve",
		Is_ctDNA == "-ve" ~ "ctDNA -ve"
	)) %>%
	ggplot(aes(x = INSERT_SIZE, y = 100*`%_READS`, group = sample_name)) +
	geom_step(stat = "identity", size = .5, alpha = .25, color = "#cb181d") +
	scale_x_sqrt(limits = c(5, NA),
		     breaks = c(10, 15, 25, 50, 100, 200, 300, 400, 500),
		     labels = c(10, 15, 25, 50, 100, 200, 300, 400, 500)) +
	scale_y_continuous(limits = c(0, NA)) +
	xlab("Insert size (bp)") +
	ylab("Fraction of Read Pairs (%)") +
	facet_wrap(~Is_ctDNA, ncol = 1, scales = "free_y") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	geom_rug(data = primer_set, mapping = aes(x = insert_size), size = .2, , sides = "b", inherit.aes = FALSE)

pdf(file = "../res/Insert_Size_Distribution_Fraction_Primer_Filtered_Fragment_Filtered.pdf", width = 6, height = 5)
print(plot_)
dev.off()
