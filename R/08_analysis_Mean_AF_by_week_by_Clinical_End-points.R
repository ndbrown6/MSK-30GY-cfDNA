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

posterior_probability = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::mutate(timepoint_weeks_since_start_of_RT = floor(timepoint_days_since_start_of_RT/7)) %>%
			dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
				timepoint_weeks_since_start_of_RT < 0 ~ "Pre-treatment",
				TRUE ~ paste0("wk", timepoint_weeks_since_start_of_RT)
			))

idx_metrics_ft = readr::read_tsv(file = url_idx_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      	 readr::type_convert() %>%
		 dplyr::select(sample_name = SAMPLE_NAME, contig = CHROMOSOME, aligned_reads = ALIGNED_READS, fragment_length = FRAGMENT_LENGTH) %>%
		 dplyr::left_join(dplyr::tibble(contig = target_contigs,
				       		chromosome = names(target_contigs)),
				  by = "contig") %>%
		 dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD) %>%
		 dplyr::filter(!is.na(chromosome)) %>%
		 dplyr::left_join(manifest, by = "sample_name") %>%
		 dplyr::filter(chromosome == hpv_type_wes_wgs) %>%
		 dplyr::left_join(mutation_smry %>%
				       dplyr::filter(FILTER == "PASS") %>%
				       dplyr::group_by(Tumor_Sample_Barcode) %>%
				       dplyr::summarize(mean_af = mean(t_maf)) %>%
				       dplyr::ungroup() %>%
				       dplyr::rename(sample_name = Tumor_Sample_Barcode),
				       by = "sample_name") %>%
		dplyr::mutate(timepoint_weeks_since_start_of_RT = floor(timepoint_days_since_start_of_RT/7)) %>%
		dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
			timepoint_weeks_since_start_of_RT < 0 ~ "Pre-treatment",
			TRUE ~ paste0("wk", timepoint_weeks_since_start_of_RT)
		))

#==================================================
# Relative Mean AF Hclust
#==================================================
data_ = idx_metrics_ft %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(mean_af = mean(mean_af+1E-5, na.rm = TRUE)) %>%
	reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	tidyr::drop_na() %>%
	dplyr::mutate(rel_wk0 = 0,
		      rel_wk1 = log2(wk1/`Pre-treatment`),
		      rel_wk2 = log2(wk2/`Pre-treatment`),
		      rel_wk3 = log2(wk3/`Pre-treatment`)) %>%
	dplyr::rename(raw_wk0 = `Pre-treatment`,
		      raw_wk1 = wk1,
		      raw_wk2 = wk2,
		      raw_wk3 = wk3)

clust_ = data_ %>%
	 dplyr::select(rel_wk1, rel_wk2, rel_wk3) %>%
	 dist(method = "euclidean") %>%
	 hclust(method = "ward.D") %>%
	 cutree(k = 4) %>%
	 dplyr::as_tibble() %>%
	 dplyr::rename(canonical_groups = value) %>%
	 dplyr::bind_cols(data_) %>%
	 dplyr::select(patient_id_mskcc, canonical_groups)

#==================================================
# Logistic Regeression Relative/ Raw Mean AF
#==================================================
data_ = idx_metrics_ft %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(mean_af = mean(mean_af+1e-5, na.rm = TRUE)) %>%
	reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	tidyr::drop_na() %>%
	dplyr::mutate(rel_wk1 = wk1/`Pre-treatment`,
		      rel_wk2 = wk2/`Pre-treatment`,
		      rel_wk3 = wk3/`Pre-treatment`) %>%
	dplyr::rename(raw_wk0 = `Pre-treatment`,
		      raw_wk1 = wk1,
		      raw_wk2 = wk2,
		      raw_wk3 = wk3) %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::left_join(clust_, by = "patient_id_mskcc") %>%
	dplyr::select(patient_id_mskcc,
		      composite_end_point,
		      canonical_groups,
		      raw_wk0, raw_wk1, raw_wk2, raw_wk3,
		      rel_wk1, rel_wk2, rel_wk3,
		      hpv_copynumber = hpv_panel_copy_number,
		      tumor_volume = plan_volume) %>%
	dplyr::mutate(composite_end_point = case_when(
		composite_end_point ~ 1,
		!composite_end_point ~ 0
	))

tmp_ <- data_ %>%
	dplyr::select(composite_end_point, contains("raw"), contains("rel")) %>%
	dplyr::mutate(raw_wk0 = scale(raw_wk0),
		      raw_wk1 = scale(raw_wk1),
		      raw_wk2 = scale(raw_wk2),
		      raw_wk3 = scale(raw_wk3),
		      rel_wk1 = scale(rel_wk1),
		      rel_wk2 = scale(rel_wk2),
		      rel_wk3 = scale(rel_wk3))

smry_ = stats::glm(composite_end_point ~ ., data = tmp_) %>%
	summary()

plot_ = smry_ %>%
	.[["coefficients"]] %>%
	as.data.frame() %>%
	tibble::rownames_to_column("variable") %>%
	dplyr::as_tibble() %>%
	dplyr::filter(variable != "(Intercept)") %>%
	dplyr::mutate(is_significant = ifelse(`Pr(>|t|)`<.1, "Yes", "No")) %>%
	dplyr::mutate(variable_cat = case_when(
		grepl("rel", variable) ~ "Relative",
		grepl("raw", variable) ~ "Absolute",
	)) %>%
	dplyr::arrange(desc(variable_cat), desc(variable)) %>%
	dplyr::mutate(variable_cat = factor(variable_cat, levels = rev(c("Relative", "Absolute")), ordered = TRUE)) %>%
	dplyr::mutate(variable = factor(variable, levels = unique(variable), ordered = TRUE)) %>%
	ggplot(aes(x = variable, ymin = 0, ymax = Estimate, y = Estimate, fill = is_significant, size = -log10(`Pr(>|t|)`), shape = variable_cat)) +
	geom_linerange(stat = "identity", size = .5) +
	geom_hline(yintercept = 0, size = 1) +
	geom_point(stat = "identity") +
	xlab("") +
	ylab("") +
	scale_fill_manual(values = c("#bdbdbd", "#e41a1c")) +
	scale_shape_manual(values = c("Relative" = 21, "Absolute" = 22)) +
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

pdf(file = "../res/Logistic_Regression_Coefficients_Composite_End-Points.pdf", width = 5, height = 4)
print(plot_)
dev.off()

tmp_ <- data_ %>%
	dplyr::select(composite_end_point, contains("raw"), canonical_groups) %>%
	dplyr::mutate(raw_wk0 = log2(raw_wk0),
		      raw_wk1 = log2(raw_wk1),
		      raw_wk2 = log2(raw_wk2),
		      raw_wk3 = log2(raw_wk3)) %>%
        dplyr::mutate(canonical_groups = case_when(
		canonical_groups == "4" ~ "FastSlow",
		canonical_groups == "3" ~ "FastSlow",
		canonical_groups == "2" ~ "Noclearance",
		canonical_groups == "1" ~ "Intermediate"
	)) %>%
	dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("FastSlow", "Intermediate", "Noclearance"), ordered = FALSE))

smry_ = stats::glm(composite_end_point ~ .:canonical_groups - canonical_groups, data = tmp_) %>%
	summary()

plot_ = smry_ %>%
	.[["coefficients"]] %>%
	as.data.frame() %>%
	tibble::rownames_to_column("variable") %>%
	dplyr::as_tibble() %>%
	dplyr::filter(variable != "(Intercept)") %>%
	dplyr::mutate(is_significant = ifelse(`Pr(>|t|)`<.1, "Yes", "No")) %>%
	dplyr::mutate(variable_cat = case_when(
		grepl("FastSlow", variable) ~ "Fast/Slow",
		grepl("Intermediate", variable) ~ "Intermediate",
		grepl("Noclearance", variable) ~ "No clearance"
	)) %>%
	dplyr::arrange(desc(variable_cat), desc(variable)) %>%
	dplyr::mutate(variable_cat = factor(variable_cat, levels = c("Fast/Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
	dplyr::mutate(variable = factor(variable, levels = unique(variable), ordered = TRUE)) %>%
	ggplot(aes(x = variable, ymin = 0, ymax = Estimate, y = Estimate, fill = is_significant, size = -log10(`Pr(>|t|)`), shape = variable_cat)) +
	geom_linerange(stat = "identity", size = .5) +
	geom_hline(yintercept = 0, size = 1) +
	geom_point(stat = "identity") +
	xlab("") +
	ylab("") +
	scale_fill_manual(values = c("#bdbdbd", "#e41a1c")) +
	scale_shape_manual(values = c("Fast/Slow" = 21, "Intermediate" = 23, "No clearance" = 24)) +
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

pdf(file = "../res/Logistic_Regression_Coefficients_Composite_End-Points_by_Relative_Classification.pdf", width = 8, height = 4)
print(plot_)
dev.off()

params = rpart.control(minsplit = 6, minbucket = 4, maxdepth = 6, cp = 0.001)
tree = rpart(composite_end_point ~ .,
	     data = data_ %>%
	     	    dplyr::select(-patient_id_mskcc) %>%
	     	    dplyr::mutate(raw_wk0 = ifelse(raw_wk0<=1e-3, 0, raw_wk0)) %>%
	     	    dplyr::mutate(raw_wk1 = ifelse(raw_wk1<=1e-3, 0, raw_wk1)) %>%
	     	    dplyr::mutate(raw_wk2 = ifelse(raw_wk2<=1e-3, 0, raw_wk2)) %>%
	     	    dplyr::mutate(raw_wk3 = ifelse(raw_wk3<=1e-3, 0, raw_wk3)),
	     method = "class", control = params)
rpart.plot(tree)