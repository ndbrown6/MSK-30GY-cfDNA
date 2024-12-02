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
	         ))

#==================================================
# Number of Read Pairs by Week
#==================================================
plot_ = idx_metrics_ft %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(aligned_reads = mean(log10(aligned_reads))) %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3", "wk5")) %>%
	dplyr::mutate(timepoint_weeks_since_start_of_RT = factor(x = timepoint_weeks_since_start_of_RT,
								 levels = c("wk5", "wk3", "wk2", "wk1", "Pre-treatment"),
								 ordered = TRUE)) %>%
	ggplot(aes(x = aligned_reads, y = timepoint_weeks_since_start_of_RT,
		   fill = timepoint_weeks_since_start_of_RT)) +
	geom_density_ridges2(stat = "density_ridges", bandwidth = .35, scale = 1.5, color = "black", alpha = 1) +
	scale_fill_viridis(discrete = TRUE) +
	scale_x_continuous(limits = c(1, 7),
			   breaks = c(2, 3, 4, 5, 6, 7),
			   labels = scientific_1e(10^c(2, 3, 4, 5, 6, 7))) +
	scale_y_discrete(breaks = c("wk5", "wk3", "wk2", "wk1", "Pre-treatment"),
			 labels = c("Week 5", "Week 3", "Week 2", "Week 1", "Pre-treatment")) +
	xlab("cfDNA Aligned HPV Read Pairs") +
	ylab("") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10)) +
	guides(fill = FALSE) +
	coord_cartesian(clip = "off")

pdf(file = "../res/Number_Read_Pairs_by_Wk.pdf", width = 3*1.15, height = 2.5*1.25)
print(plot_)
dev.off()

#==================================================
# Mean AF by Week
#==================================================
plot_ = idx_metrics_ft %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(mean_af = mean(mean_af)) %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3", "wk5")) %>%
	dplyr::mutate(timepoint_weeks_since_start_of_RT = factor(x = timepoint_weeks_since_start_of_RT,
								 levels = c("wk5", "wk3", "wk2", "wk1", "Pre-treatment"),
								 ordered = TRUE)) %>%
	ggplot(aes(x = (mean_af*100)+(1E-3), y = timepoint_weeks_since_start_of_RT,
		   fill = timepoint_weeks_since_start_of_RT)) +
	geom_density_ridges2(stat = "density_ridges", bandwidth = .35, scale = 1.5, color = "black", alpha = 1) +
	scale_fill_viridis(discrete = TRUE) +
	scale_x_log10(limits = c(NA, 110),
		      breaks = c(.001, .01, .1, 1, 10, 100),
		      labels = c(".001", ".01", ".1", "1", "10", "100")) +
	scale_y_discrete(breaks = c("wk5", "wk3", "wk2", "wk1", "Pre-treatment"),
			 labels = c("Week 5", "Week 3", "Week 2", "Week 1", "Pre-treatment")) +
	xlab("Mean PCM AF (%)") +
	ylab("") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10)) +
	guides(fill = FALSE) +
	coord_cartesian(clip = "off")

pdf(file = "../res/Mean_AF_by_Wk.pdf", width = 3*1.15, height = 2.5*1.25)
print(plot_)
dev.off()

#==================================================
# MRI Volume by Week
#==================================================
plot_ = clinical %>%
	dplyr::select(patient_id_mskcc,
		      `Pre-treatment` = MRI_rawdata_wk0,
		      `wk1` = MRI_rawdata_wk1,
		      `wk2` = MRI_rawdata_wk2,
		      `wk3` = MRI_rawdata_wk3,
		      `wk5` = MRI_rawdata_wk4) %>%
	reshape2::melt(variable.name = "timepoint_weeks_since_start_of_RT", value.name = "MRI_volume") %>%
	dplyr::left_join(manifest %>%
			 dplyr::group_by(patient_id_mskcc) %>%
			 dplyr::summarize(hpv_type_wes_wgs = unique(hpv_type_wes_wgs)),
			 by = "patient_id_mskcc") %>%
	tidyr::drop_na() %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3", "wk5")) %>%
	dplyr::mutate(timepoint_weeks_since_start_of_RT = factor(x = timepoint_weeks_since_start_of_RT,
								 levels = c("wk5", "wk3", "wk2", "wk1", "Pre-treatment"),
								 ordered = TRUE)) %>%
	ggplot(aes(x = MRI_volume, y = timepoint_weeks_since_start_of_RT,
		   fill = timepoint_weeks_since_start_of_RT)) +
	geom_density_ridges2(stat = "density_ridges", bandwidth = .095, scale = 1.5, color = "black", alpha = 1) +
	scale_fill_viridis(discrete = TRUE) +
	scale_x_log10(labels = scientific_10) +
	scale_y_discrete(breaks = c("wk5", "wk3", "wk2", "wk1", "Pre-treatment"),
			 labels = c("Week 5", "Week 3", "Week 2", "Week 1", "Pre-treatment")) +
	xlab(expression("MRI Volume ("~mm^3~")")) +
	ylab("") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10)) +
	guides(fill = FALSE) +
	coord_cartesian(clip = "off")

pdf(file = "../res/MRI_Volume_by_Wk.pdf", width = 3*1.15, height = 2.5*1.25)
print(plot_)
dev.off()

#==================================================
# Posterior Probability by Week
#==================================================
plot_ = posterior_probability %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(`Pr(x=1)` = mean(`Pr(x=1)`, na.rm=TRUE)) %>%
	tidyr::drop_na() %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3", "wk5")) %>%
	dplyr::mutate(timepoint_weeks_since_start_of_RT = factor(x = timepoint_weeks_since_start_of_RT,
								 levels = c("wk5", "wk3", "wk2", "wk1", "Pre-treatment"),
								 ordered = TRUE)) %>%
	ggplot(aes(x = `Pr(x=1)`, y = timepoint_weeks_since_start_of_RT,
		   fill = timepoint_weeks_since_start_of_RT)) +
	geom_density_ridges2(stat = "density_ridges", bandwidth = .095, scale = 1.5, color = "black", alpha = 1) +
	scale_fill_viridis(discrete = TRUE) +
	scale_x_log10(labels = scientific_10) +
	scale_y_discrete(breaks = c("wk5", "wk3", "wk2", "wk1", "Pre-treatment"),
			 labels = c("Week 5", "Week 3", "Week 2", "Week 1", "Pre-treatment")) +
	xlab("Posterior Probability") +
	ylab("") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10)) +
	guides(fill = FALSE) +
	coord_cartesian(clip = "off")

pdf(file = "../res/Posterior_Probability_by_Wk.pdf", width = 3*1.15, height = 2.5*1.25)
print(plot_)
dev.off()

#==================================================
# tSNE Mean AF
#==================================================
smry_ = idx_metrics_ft %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(mean_af = mean(mean_af+1E-5, na.rm = TRUE)) %>%
	reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	dplyr::select(-wk0, -wk4, -wk5, -wk6, -`wk7+`) %>%
	tidyr::drop_na() %>%
	dplyr::mutate(wk1 = log2(wk1/`Pre-treatment`),
		      wk2 = log2(wk2/`Pre-treatment`),
		      wk3 = log2(wk3/`Pre-treatment`),
		      `Pre-treatment` = 0)

fit_ = smry_ %>%
       dplyr::select(where(is.numeric)) %>%
       Rtsne::Rtsne(perplexity = 10, theta = 0, pca = TRUE, pca_center = FALSE, pca_scale = FALSE, normalize = FALSE, exaggeration_factor = 5)

plot_ = smry_ %>%
	dplyr::bind_cols(fit_$Y %>%
			 dplyr::as_tibble()) %>%
	dplyr::mutate(drop = case_when(
		(V1 <= (-7) & V2 >= (3)) ~ "Fast",
		(V1 > (-7) & V1 <= (1) & V2 < (5)) & wk1 < 0 ~ "Slow",
		(V1 > (-1) & V2 > (5)) ~ "Initial increase",
		(V1 > (10) & wk3 > 0) ~ "No clearance",
		TRUE ~ "No discernible pattern"
	)) %>%
	dplyr::mutate(drop = factor(drop, levels = c("Fast", "Slow", "Initial increase", "No clearance", "No discernible pattern"))) %>%
	ggplot(aes(x = V1, y = V2, color = drop, shape = drop)) +
	geom_smooth(mapping = aes(x = V1, y = V2),
		    stat = "smooth", formula = y ~ x, method = "loess",
		    se = FALSE, span = 1, linetype = 2, color = "#636363", size = .25, inherit.aes = FALSE) +
	geom_point(stat = "identity", size = 3) +
	scale_shape_manual(values = c("Fast" = 21, "Slow" = 22, "Initial increase" = 23, "No clearance" = 24, "No discernible pattern" = 3)) +
	scale_color_brewer(type = "qual", palette = 2) +
	scale_x_continuous() +
	scale_y_continuous() +
	xlab("tSNE 1") +
	ylab("tSNE 2") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 7)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = guide_legend(title = "ctDNA clearance"),
	       shape = guide_legend(title = "ctDNA clearance"))

pdf(file = "../res/tSNE_Mean_AF.pdf", width = 3*1.95, height = 2.5*1.5)
print(plot_)
dev.off()

plot_ = smry_ %>%
	dplyr::bind_cols(fit_$Y %>%
			 dplyr::as_tibble()) %>%
	dplyr::mutate(drop = case_when(
		(V1 <= (-7) & V2 >= (3)) ~ "Fast",
		(V1 > (-7) & V1 <= (1) & V2 < (5)) & wk1 < 0 ~ "Slow",
		(V1 > (-1) & V2 > (5)) ~ "Initial increase",
		(V1 > (10) & wk3 > 0) ~ "No clearance",
		TRUE ~ "No discernible pattern"
	)) %>%
	dplyr::mutate(drop = factor(drop, levels = c("Fast", "Slow", "Initial increase", "No clearance", "No discernible pattern"))) %>%
	dplyr::filter(drop != "No discernible pattern") %>%
	dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3, drop) %>%
	reshape2::melt() %>%
	ggplot(aes(x = variable, y = value, group = patient_id_mskcc, color = drop, shape = drop)) +
	geom_hline(yintercept = 0, color = "goldenrod3", size = .5, linetype = 3) +
	geom_path(stat = "identity", alpha = 1, size = .35) +
	geom_point(stat = "identity", fill = "white", size = 2.5, alpha = 1) +
	scale_shape_manual(values = c("Fast" = 21, "Slow" = 22, "Initial increase" = 23, "No clearance" = 24, "No discernible pattern" = 3)) +
	scale_color_brewer(type = "qual", palette = 2) +
	scale_x_discrete(breaks = c("Pre-treatment", "wk1", "wk2", "wk3"),
			 labels = c("Pre-treat\n-ment", "wk1", "wk2", "wk3")) +
	scale_y_continuous(limits = c(-15, 10)) +
	xlab("") +
	ylab(expression(Log[2]~"Ratio Mean AF")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 7)),
	      axis.text.x = element_text(size = 8),
	      axis.text.y = element_text(size = 8),
	      strip.background = element_blank()) +
	guides(color = FALSE, shape = FALSE) +
	facet_wrap(~drop, nrow = 1, scales = "free_y")

pdf(file = "../res/tSNE_Mean_AF_Canonical.pdf", width = 3*4, height = 2.5*1.0)
print(plot_)
dev.off()

centroids = smry_ %>%
	    dplyr::bind_cols(fit_$Y %>%
			     dplyr::as_tibble()) %>%
	    dplyr::mutate(drop = case_when(
		    (V1 <= (-7) & V2 >= (3)) ~ "Fast",
		    (V1 > (-7) & V1 <= (1) & V2 < (5)) & wk1 < 0 ~ "Slow",
		    (V1 > (-1) & V2 > (5)) ~ "Initial increase",
		    (V1 > (10) & wk3 > 0) ~ "No clearance",
		    TRUE ~ "No discernible pattern"
	    )) %>%
	    dplyr::filter(drop != "No discernible pattern") %>%
	    dplyr::group_by(drop) %>%
	    dplyr::summarize(mean_v1 = mean(V1),
			     mean_v2 = mean(V2))

predictions = smry_ %>%
	      dplyr::bind_cols(fit_$Y %>%
			       dplyr::as_tibble()) %>%
	      dplyr::mutate(drop = case_when(
		    (V1 <= (-7) & V2 >= (3)) ~ "Fast",
		    (V1 > (-7) & V1 <= (1) & V2 < (5)) & wk1 < 0 ~ "Slow",
		    (V1 > (-1) & V2 > (5)) ~ "Initial increase",
		    (V1 > (10) & wk3 > 0) ~ "No clearance",
		    TRUE ~ "No discernible pattern"
	      )) %>%
	      dplyr::filter(drop == "No discernible pattern") %>%
	      dplyr::left_join(centroids %>%
			       dplyr::rename(centroid = drop) %>%
			       dplyr::mutate(drop = "No discernible pattern")) %>%
	      dplyr::mutate(d = sqrt((V1 - mean_v1)^2 + (V2 - mean_v2)^2)) %>%
	      dplyr::group_by(patient_id_mskcc) %>%
	      dplyr::summarize(nearest_centroid = centroid[which.min(d)])

plot_ = smry_ %>%
	dplyr::bind_cols(fit_$Y %>%
			 dplyr::as_tibble()) %>%
	dplyr::left_join(predictions, by = "patient_id_mskcc") %>%
	dplyr::filter(!is.na(nearest_centroid)) %>%
	dplyr::rename(drop = nearest_centroid) %>%
	dplyr::mutate(drop = factor(drop, levels = c("Fast", "Slow", "Initial increase", "No clearance", "No discernible pattern"))) %>%
	dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3, drop) %>%
	reshape2::melt() %>%
	ggplot(aes(x = variable, y = value, group = patient_id_mskcc, color = drop, shape = drop)) +
	geom_hline(yintercept = 0, color = "goldenrod3", size = .5, linetype = 3) +
	geom_path(stat = "identity", alpha = 1, size = .35) +
	geom_point(stat = "identity", fill = "white", size = 2.5, alpha = 1) +
	scale_shape_manual(values = c("Fast" = 21, "Slow" = 22, "Initial increase" = 23, "No clearance" = 24, "No discernible pattern" = 3)) +
	scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a")) +
	scale_x_discrete(breaks = c("Pre-treatment", "wk1", "wk2", "wk3"),
			 labels = c("Pre-treat\n-ment", "wk1", "wk2", "wk3")) +
	scale_y_continuous(limits = c(-15, 10)) +
	xlab("") +
	ylab(expression(Log[2]~"Ratio Mean AF")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 7)),
	      axis.text.x = element_text(size = 8),
	      axis.text.y = element_text(size = 8),
	      strip.background = element_blank()) +
	guides(color = FALSE, shape = FALSE) +
	facet_wrap(~drop, nrow = 1, scales = "free_y")

pdf(file = "../res/tSNE_Mean_AF_Unknown.pdf", width = 3*3, height = 2.5*1.0)
print(plot_)
dev.off()

#==================================================
# Log2 Ratio cfDNA Read Pairs / Mean AF
#==================================================
smry_af = idx_metrics_ft %>%
	  dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	  dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	  dplyr::summarize(mean_af = mean(mean_af+1E-5, na.rm = TRUE)) %>%
	  reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			  fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	  dplyr::select(-wk0, -wk4, -wk5, -wk6, -`wk7+`) %>%
	  dplyr::mutate(wk1 = log2(wk1/`Pre-treatment`),
		        wk2 = log2(wk2/`Pre-treatment`),
		        wk3 = log2(wk3/`Pre-treatment`),
		        `Pre-treatment` = 0)

smry_reads = idx_metrics_ft %>%
	     dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	     dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	     dplyr::summarize(aligned_reads = mean(aligned_reads, na.rm = TRUE)) %>%
	     reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			     fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "aligned_reads") %>%
	     dplyr::select(-wk0, -wk4, -wk5, -wk6, -`wk7+`) %>%
	     dplyr::mutate(wk1 = log2(wk1/`Pre-treatment`),
			   wk2 = log2(wk2/`Pre-treatment`),
			   wk3 = log2(wk3/`Pre-treatment`),
			   `Pre-treatment` = 0)

plot_ = smry_af %>%
	reshape2::melt(variable.name = "week", value.name = "log2_af") %>%
	dplyr::left_join(smry_reads %>%
			 reshape2::melt(variable.name = "week", value.name = "log2_reads"),
			 by = c("patient_id_mskcc", "week")) %>%
	dplyr::filter(week != "Pre-treatment") %>%
	dplyr::mutate(week = case_when(
		week == "wk1" ~ "Week 1",
		week == "wk2" ~ "Week 2",
		week == "wk3" ~ "Week 3"
	)) %>%
	ggplot(aes(x = log2_af, y = log2_reads, color = week, shape = week)) +
	geom_abline(intercept = 0, slope = 1, color = "#969696", size = .5, linetype = 1) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x +0, se = FALSE, color = "goldenrod3", size = .5, linetype = 1) +
	geom_point(stat = "identity", size = 2, fill = "white", alpha = .75) +
	scale_shape_manual(values = c("Week 1" = 21, "Week 2" = 22, "Week 3" = 23)) +
	scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a")) +
	geom_hline(data = smry_reads %>%
		   	  reshape2::melt(variable.name = "week", value.name = "log2_reads") %>%
		   	  dplyr::filter(week != "Pre-treatment") %>%
		   	  dplyr::group_by(week) %>%
		   	  dplyr::summarize(mean_log2_reads = mean(log2_reads, na.rm=TRUE)) %>%
		   	  dplyr::mutate(week = case_when(
						week == "wk1" ~ "Week 1",
				  		week == "wk2" ~ "Week 2",
				  		week == "wk3" ~ "Week 3"
			  )),
		   aes(yintercept = mean_log2_reads), color = "#74c476", size = .5, linetype = 3) +
	geom_vline(data = smry_af %>%
		   	  reshape2::melt(variable.name = "week", value.name = "log2_af") %>%
		   	  dplyr::filter(week != "Pre-treatment") %>%
		   	  dplyr::group_by(week) %>%
		   	  dplyr::summarize(mean_log2_af = mean(log2_af, na.rm=TRUE)) %>%
		   	  dplyr::mutate(week = case_when(
				  		week == "wk1" ~ "Week 1",
				  		week == "wk2" ~ "Week 2",
				  		week == "wk3" ~ "Week 3"
			  )),
		   aes(xintercept = mean_log2_af), color = "#74c476", size = .5, linetype = 3) +
	scale_x_continuous(limits = c(-15, 10)) +
	scale_y_continuous(limits = c(-15, 10)) +
	xlab("Relative Mean AF") +
	ylab("Relative cfDNA HPV Read Pairs") +
	stat_cor(method = "spearman") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(shape = FALSE, color = FALSE)

pdf(file = "../res/Relative_Change_Mean_AF_cfDNA_Read_Pairs.pdf", width = 3*3, height = 2.5*1.0)
print(plot_)
dev.off()

plot_ = smry_af %>%
	reshape2::melt(variable.name = "week", value.name = "log2_af") %>%
	dplyr::left_join(smry_reads %>%
			 reshape2::melt(variable.name = "week", value.name = "log2_reads"),
			 by = c("patient_id_mskcc", "week")) %>%
	dplyr::filter(week != "Pre-treatment") %>%
	reshape2::melt(variable.name = "category", value.name = "log2_ratio") %>%
	dplyr::mutate(week = case_when(
		week == "wk1" ~ "Week 1",
		week == "wk2" ~ "Week 2",
		week == "wk3" ~ "Week 3"
	)) %>%
	tidyr::drop_na() %>%
	ggplot(aes(x = category, y = 2^log2_ratio)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	scale_x_discrete(breaks = c("log2_af", "log2_reads"),
			 labels = c("\nMean\nAF", "\ncfDNA HPV\nReads")) +
	scale_y_log10(breaks = c(.001, .01, .1, 1, 10, 100, 1000),
		      labels = scientific_10) +
	xlab(" ") +
	ylab("Relative Change") +
	geom_signif(stat = "signif",
		    comparisons = list(c("log2_af", "log2_reads")),
		    test = "wilcox.test",
		    test.args = list(alternative = "less")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	facet_wrap(~week)

