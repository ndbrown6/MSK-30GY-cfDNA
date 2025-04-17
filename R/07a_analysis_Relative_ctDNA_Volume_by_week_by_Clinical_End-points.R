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
# Relative cfDNA Read Pairs / Mean AF / MRI Volume
#==================================================
smry_af = idx_metrics_ft %>%
	  dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	  dplyr::summarize(mean_af = mean(mean_af+1E-5, na.rm = TRUE)) %>%
	  reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			  fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	  dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	  dplyr::mutate(wk1 = log2(wk1/`Pre-treatment`),
		        wk2 = log2(wk2/`Pre-treatment`),
		        wk3 = log2(wk3/`Pre-treatment`),
		        `Pre-treatment` = 0)

smry_reads = idx_metrics_ft %>%
	     dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	     dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	     dplyr::summarize(aligned_reads = mean(aligned_reads+1, na.rm = TRUE)) %>%
	     reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			     fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "aligned_reads") %>%
	     dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	     dplyr::mutate(wk1 = log2(wk1/`Pre-treatment`),
			   wk2 = log2(wk2/`Pre-treatment`),
			   wk3 = log2(wk3/`Pre-treatment`),
			   `Pre-treatment` = 0)

smry_mri = clinical %>%
	   dplyr::select(patient_id_mskcc,
		        `Pre-treatment` = MRI_rawdata_wk0,
		        `wk1` = MRI_rawdata_wk2,
		        `wk2` = MRI_rawdata_wk3,
		        `wk3` = MRI_rawdata_wk4) %>%
	   dplyr::mutate(wk1 = log2(wk1/`Pre-treatment`),
		         wk2 = log2(wk2/`Pre-treatment`),
		         wk3 = log2(wk3/`Pre-treatment`),
			 `Pre-treatment` = 0)
	   
plot_ = smry_af %>%
	reshape2::melt(variable.name = "week", value.name = "log2_af") %>%
	dplyr::full_join(smry_reads %>%
			 reshape2::melt(variable.name = "week", value.name = "log2_reads"),
			 by = c("patient_id_mskcc", "week")) %>%
	dplyr::full_join(smry_mri %>%
			 reshape2::melt(variable.name = "week", value.name = "log2_volume"),
			 by = c("patient_id_mskcc", "week")) %>%
	dplyr::filter(week != "Pre-treatment") %>%
	dplyr::mutate(week = case_when(
		week == "wk1" ~ "Week 1",
		week == "wk2" ~ "Week 2",
		week == "wk3" ~ "Week 3"
	)) %>%
	dplyr::select(log2_af, log2_reads, week) %>%
	tidyr::drop_na() %>%
	ggplot(aes(x = log2_af, y = log2_reads, color = week, shape = week)) +
	geom_abline(intercept = 0, slope = 1, color = "#969696", size = .5, linetype = 1) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x +0, se = FALSE, color = "goldenrod3", size = .5, linetype = 1) +
	geom_point(stat = "identity", size = 2, fill = "white", alpha = .75) +
	scale_shape_manual(values = c("Week 1" = 21, "Week 2" = 22, "Week 3" = 23)) +
	scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a")) +
	geom_hline(yintercept = 0, color = "#74c476", size = .5, linetype = 3) +
	geom_vline(xintercept = 0, color = "#74c476", size = .5, linetype = 3) +
	scale_x_continuous(limits = c(-15, 10)) +
	scale_y_continuous(limits = c(-15, 10)) +
	xlab("Relative Mean AF") +
	ylab("Relative cfDNA HPV Read Pairs") +
	stat_cor(method = "spearman", label.x = -15, label.y = 9) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(shape = FALSE, color = FALSE)

pdf(file = "../res/Relative_Change_Mean_AF_cfDNA_Read_Pairs.pdf", width = 3*2.75, height = 2.95*1.0)
print(plot_)
dev.off()

plot_ = smry_af %>%
	reshape2::melt(variable.name = "week", value.name = "log2_af") %>%
	dplyr::full_join(smry_reads %>%
			 reshape2::melt(variable.name = "week", value.name = "log2_reads"),
			 by = c("patient_id_mskcc", "week")) %>%
	dplyr::full_join(smry_mri %>%
			 reshape2::melt(variable.name = "week", value.name = "log2_volume"),
			 by = c("patient_id_mskcc", "week")) %>%
	dplyr::filter(week != "Pre-treatment") %>%
	dplyr::mutate(week = case_when(
		week == "wk1" ~ "Week 1",
		week == "wk2" ~ "Week 2",
		week == "wk3" ~ "Week 3"
	)) %>%
	dplyr::select(patient_id_mskcc, log2_af, log2_volume, week) %>%
	tidyr::drop_na() %>%
	ggplot(aes(x = log2_af, y = log2_volume, color = week, shape = week)) +
	geom_abline(intercept = 0, slope = 1, color = "#969696", size = .5, linetype = 1) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x +0, se = FALSE, color = "goldenrod3", size = .5, linetype = 1) +
	geom_point(stat = "identity", size = 2, fill = "white", alpha = .75) +
	scale_shape_manual(values = c("Week 1" = 21, "Week 2" = 22, "Week 3" = 23)) +
	scale_color_manual(values = c("#1b9e77", "#d95f02", "#e7298a")) +
	geom_hline(yintercept = 0, color = "#74c476", size = .5, linetype = 3) +
	geom_vline(xintercept = 0, color = "#74c476", size = .5, linetype = 3) +
	scale_x_continuous(limits = c(-15, 10)) +
	scale_y_continuous(limits = c(-3, 1)) +
	xlab("Relative Mean AF") +
	ylab("Relative MRI Volume") +
	stat_cor(method = "spearman", label.x = -15, label.y = .9) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(shape = FALSE, color = FALSE)

pdf(file = "../res/Relative_Change_Mean_AF_MRI_Volume.pdf", width = 3*2.75, height = 2.95*1.0)
print(plot_)
dev.off()

#==================================================
# Relative/ Raw Mean AF tSNE, Hclust
#==================================================
smry_ = idx_metrics_ft %>%
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

set.seed(3)

tsne_ = smry_ %>%
        dplyr::select(rel_wk1, rel_wk2, rel_wk3) %>%
        Rtsne::Rtsne(perplexity = 10, theta = 0, pca = TRUE, pca_center = FALSE, pca_scale = FALSE, normalize = FALSE, exaggeration_factor = 5)

clust_ = smry_ %>%
	 dplyr::select(rel_wk1, rel_wk2, rel_wk3) %>%
	 dist(method = "euclidean") %>%
	 hclust(method = "ward.D") %>%
	 cutree(k = 4)

plot_ = smry_ %>%
	dplyr::bind_cols(tsne_$Y %>%
			 dplyr::as_tibble() %>%
			 dplyr::select(tSNE1 = V1, tSNE2 = V2)) %>%
	dplyr::mutate(canonical_groups = as.character(clust_)) %>%
	dplyr::mutate(canonical_groups = case_when(
		canonical_groups == "4" ~ "Fast",
		canonical_groups == "3" ~ "Slow",
		canonical_groups == "2" ~ "No clearance",
		canonical_groups == "1" ~ "Intermediate"
	)) %>%
	dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
	ggplot(aes(x = tSNE1, y = tSNE2, color = canonical_groups, shape = canonical_groups)) +
	geom_smooth(mapping = aes(x = tSNE1, y = tSNE2),
		    stat = "smooth", formula = y ~ x, method = "loess",
		    se = FALSE, span = 1, linetype = 2, color = "#636363", size = .25, inherit.aes = FALSE) +
	geom_point(stat = "identity", size = 3, fill = "white") +
	scale_shape_manual(values = c("Fast" = 21, "Slow" = 22, "Intermediate" = 23, "No clearance" = 24)) +
	scale_color_manual(values = c("Fast" = "#1b9e77", "Slow" = "#d95f02", "Intermediate" = "#7570b3", "No clearance" = "#e7298a")) +
	scale_x_continuous() +
	scale_y_continuous() +
	xlab("Component 1") +
	ylab("Component 2") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 7)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = guide_legend(title = "ctDNA clearance"),
	       shape = guide_legend(title = "ctDNA clearance"))

pdf(file = "../res/tSNE_Relative_Mean_AF.pdf", width = 3*1.95, height = 2.5*1.5)
print(plot_)
dev.off()

plot_ = smry_ %>%
	dplyr::select(rel_wk1, rel_wk2, rel_wk3) %>%
	dist(method = "euclidean") %>%
	hclust(method = "ward.D") %>%
	ggdendrogram(rotate = FALSE, size = 2) +
	geom_hline(yintercept = 35, color = "goldenrod3", linetype = 2) +
	geom_hline(yintercept = 25, color = "goldenrod3", linetype = 2) +
	geom_hline(yintercept = 15, color = "goldenrod3", linetype = 2)

pdf(file = "../res/Hclust_Relative_Mean_AF.pdf", width = 3*1.95, height = 2.5*1.5)
print(plot_)
dev.off()

plot_ = smry_ %>%
	dplyr::select(patient_id_mskcc, rel_wk0, rel_wk1, rel_wk2, rel_wk3) %>%
	dplyr::mutate(canonical_groups = as.character(clust_)) %>%
	dplyr::mutate(canonical_groups = case_when(
		canonical_groups == "4" ~ "Fast",
		canonical_groups == "3" ~ "Slow",
		canonical_groups == "2" ~ "No clearance",
		canonical_groups == "1" ~ "Intermediate"
	)) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups")) %>%
	dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
	ggplot(aes(x = variable, y = value, group = patient_id_mskcc)) +
	geom_hline(yintercept = 0, color = "goldenrod3", size = .5, linetype = 3) +
	geom_path(stat = "identity", size = .55, color = "grey75", alpha = .5) +
	geom_path(data = smry_ %>%
		  	 dplyr::select(patient_id_mskcc, rel_wk0, rel_wk1, rel_wk2, rel_wk3) %>%
		  	 dplyr::mutate(canonical_groups = as.character(clust_)) %>%
		  	 dplyr::mutate(canonical_groups = case_when(
				 		canonical_groups == "4" ~ "Fast",
				 		canonical_groups == "3" ~ "Slow",
				 		canonical_groups == "2" ~ "No clearance",
				 		canonical_groups == "1" ~ "Intermediate"
			 )) %>%
		  	 reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups")) %>%
		  	 dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
		  	 dplyr::group_by(canonical_groups, variable) %>%
		  	 dplyr::summarize(value = mean(value)) %>%
		  	 dplyr::ungroup(),
		  mapping = aes(x = variable, y = value, color = canonical_groups, group = canonical_groups),
		  size = 1.05, alpha = 1, inherit.aes = FALSE) +
	geom_point(data = smry_ %>%
		   	  dplyr::select(patient_id_mskcc, rel_wk0, rel_wk1, rel_wk2, rel_wk3) %>%
		  	  dplyr::mutate(canonical_groups = as.character(clust_)) %>%
		  	  dplyr::mutate(canonical_groups = case_when(
				 		canonical_groups == "4" ~ "Fast",
				 		canonical_groups == "3" ~ "Slow",
				 		canonical_groups == "2" ~ "No clearance",
				 		canonical_groups == "1" ~ "Intermediate"
			  )) %>%
		  	  reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups")) %>%
		  	  dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
		  	  dplyr::group_by(canonical_groups, variable) %>%
		  	  dplyr::summarize(value = mean(value)) %>%
		  	  dplyr::ungroup(),
		   mapping = aes(x = variable, y = value, color = canonical_groups, shape = canonical_groups, group = canonical_groups),
		   fill = "white", size = 3, alpha = 1) +
	scale_shape_manual(values = c("Fast" = 21, "Slow" = 22, "Intermediate" = 23, "No clearance" = 24)) +
	scale_color_manual(values = c("Fast" = "#1b9e77", "Slow" = "#d95f02", "Intermediate" = "#7570b3", "No clearance" = "#e7298a")) +
	scale_x_discrete(breaks = c("rel_wk0", "rel_wk1", "rel_wk2", "rel_wk3"),
			 labels = c("Pre-treat\n-ment", "wk1", "wk2", "wk3")) +
	scale_y_continuous(limits = c(-15, 10)) +
	xlab("") +
	ylab("Relative Mean AF") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 9, margin = margin(t = 5)),
	      axis.text.y = element_text(size = 9),
	      strip.background = element_blank()) +
	guides(color = FALSE, shape = FALSE) +
	facet_wrap(~canonical_groups, nrow = 1, ncol = 4, scales = "free")

pdf(file = "../res/Relative_Mean_AF_by_Wk_Hclust.pdf", width = 3*4, height = 2.5*1.0)
print(plot_)
dev.off()

f_ = smry_ %>%
     dplyr::mutate(canonical_groups = as.character(clust_)) %>%
     dplyr::mutate(canonical_groups = case_when(
			canonical_groups == "4" ~ "Fast",
			canonical_groups == "3" ~ "Slow",
			canonical_groups == "2" ~ "No clearance",
			canonical_groups == "1" ~ "Intermediate"
     )) %>%
     dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
     dplyr::select(hypoxia_resolution, neck_dissection_yes_no, composite_end_point, canonical_groups)

f_ = f_ %>%
     dplyr::group_by(hypoxia_resolution, canonical_groups) %>%
     dplyr::summarize(n = n()) %>%
     dplyr::ungroup() %>%
     reshape2::dcast(formula = canonical_groups ~ hypoxia_resolution, fill = 0, value.var = "n") %>%
     dplyr::mutate(f = persistent+resolved) %>%
     dplyr::mutate(persistent = 100*persistent/f,
		   resolved = 100*resolved/f) %>%
     dplyr::select(-f) %>%
     dplyr::left_join(f_ %>%
    		      dplyr::group_by(neck_dissection_yes_no, canonical_groups) %>%
		      dplyr::summarize(n = n()) %>%
		      dplyr::ungroup() %>%
		      reshape2::dcast(formula = canonical_groups ~ neck_dissection_yes_no, fill = 0, value.var = "n") %>%
		      dplyr::mutate(f = yes+no) %>%
		      dplyr::mutate(yes = 100*yes/f,
				    no = 100*no/f) %>%
			dplyr::select(-f), by = "canonical_groups") %>%
     dplyr::left_join(f_ %>%
		      dplyr::group_by(composite_end_point, canonical_groups) %>%
		      dplyr::summarize(n = n()) %>%
		      dplyr::ungroup() %>%
		      reshape2::dcast(formula = canonical_groups ~ composite_end_point, fill = 0, value.var = "n") %>%
		      dplyr::mutate(f = `FALSE`+`TRUE`) %>%
		      dplyr::mutate(`FALSE` = 100*`FALSE`/f,
				    `TRUE` = 100*`TRUE`/f) %>%
		      dplyr::select(-f), by = "canonical_groups") %>%
     dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
     dplyr::arrange(canonical_groups) %>%
     tibble::column_to_rownames("canonical_groups")

pdf(file = "../res/Relative_Classification_by_Clinical_End_Points.pdf", width = 7, height = 7)
corrplot(corr = as.matrix(f_),
	 method = "circle",
	 type = "full",
	 is.corr = FALSE,
	 outline = FALSE,
	 addgrid.col = NA,
	 col = COL2('RdBu', 200))
dev.off()

f_ = smry_ %>%
     dplyr::mutate(canonical_groups = as.character(clust_)) %>%
     dplyr::mutate(canonical_groups = case_when(
			canonical_groups == "4" ~ "Fast",
			canonical_groups == "3" ~ "Slow",
			canonical_groups == "2" ~ "No clearance",
			canonical_groups == "1" ~ "Intermediate"
     )) %>%
     dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
     dplyr::select(hypoxia_resolution, neck_dissection_yes_no, composite_end_point, canonical_groups)

f_ %>%
dplyr::group_by(hypoxia_resolution, canonical_groups) %>%
dplyr::summarize(n = n()) %>%
dplyr::ungroup() %>%
reshape2::dcast(formula = canonical_groups ~ hypoxia_resolution, fill = 0, value.var = "n") %>%
dplyr::left_join(f_ %>%
    		 dplyr::group_by(neck_dissection_yes_no, canonical_groups) %>%
		 dplyr::summarize(n = n()) %>%
		 dplyr::ungroup() %>%
		 reshape2::dcast(formula = canonical_groups ~ neck_dissection_yes_no, fill = 0, value.var = "n"),
		 by = "canonical_groups") %>%
dplyr::left_join(f_ %>%
		 dplyr::group_by(composite_end_point, canonical_groups) %>%
		 dplyr::summarize(n = n()) %>%
		 dplyr::ungroup() %>%
		 reshape2::dcast(formula = canonical_groups ~ composite_end_point, fill = 0, value.var = "n"),
		 by = "canonical_groups") %>%
dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
dplyr::arrange(canonical_groups) %>%
tibble::column_to_rownames("canonical_groups") %>%
pander::pander(caption = "Relative Classification by Clinical End Points")

plot_ = smry_ %>%
	dplyr::select(patient_id_mskcc, raw_wk0, raw_wk1, raw_wk2, raw_wk3) %>%
	dplyr::mutate(canonical_groups = as.character(clust_)) %>%
	dplyr::mutate(canonical_groups = case_when(
		canonical_groups == "4" ~ "Fast",
		canonical_groups == "3" ~ "Slow",
		canonical_groups == "2" ~ "No clearance",
		canonical_groups == "1" ~ "Intermediate"
	)) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups")) %>%
	dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
	ggplot(aes(x = variable, y = value, group = patient_id_mskcc)) +
	geom_hline(yintercept = 0, color = "goldenrod3", size = .5, linetype = 3) +
	geom_path(stat = "identity", size = .55, color = "grey75", alpha = .5) +
	geom_path(data = smry_ %>%
		  	 dplyr::select(patient_id_mskcc, raw_wk0, raw_wk1, raw_wk2, raw_wk3) %>%
		  	 dplyr::mutate(canonical_groups = as.character(clust_)) %>%
		  	 dplyr::mutate(canonical_groups = case_when(
				 		canonical_groups == "4" ~ "Fast",
				 		canonical_groups == "3" ~ "Slow",
				 		canonical_groups == "2" ~ "No clearance",
				 		canonical_groups == "1" ~ "Intermediate"
			 )) %>%
		  	 reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups")) %>%
		  	 dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
		  	 dplyr::group_by(canonical_groups, variable) %>%
		  	 dplyr::summarize(value = mean(value)) %>%
		  	 dplyr::ungroup(),
		  mapping = aes(x = variable, y = value, color = canonical_groups, group = canonical_groups),
		  size = 1.05, alpha = 1, inherit.aes = FALSE) +
	geom_point(data = smry_ %>%
		   	  dplyr::select(patient_id_mskcc, raw_wk0, raw_wk1, raw_wk2, raw_wk3) %>%
		  	  dplyr::mutate(canonical_groups = as.character(clust_)) %>%
		  	  dplyr::mutate(canonical_groups = case_when(
				 		canonical_groups == "4" ~ "Fast",
				 		canonical_groups == "3" ~ "Slow",
				 		canonical_groups == "2" ~ "No clearance",
				 		canonical_groups == "1" ~ "Intermediate"
			  )) %>%
		  	  reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups")) %>%
		  	  dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
		  	  dplyr::group_by(canonical_groups, variable) %>%
		  	  dplyr::summarize(value = mean(value)) %>%
		  	  dplyr::ungroup(),
		   mapping = aes(x = variable, y = value, color = canonical_groups, shape = canonical_groups, group = canonical_groups),
		   fill = "white", size = 3, alpha = 1) +
	scale_shape_manual(values = c("Fast" = 21, "Slow" = 22, "Intermediate" = 23, "No clearance" = 24)) +
	scale_color_manual(values = c("Fast" = "#1b9e77", "Slow" = "#d95f02", "Intermediate" = "#7570b3", "No clearance" = "#e7298a")) +
	scale_x_discrete(breaks = c("raw_wk0", "raw_wk1", "raw_wk2", "raw_wk3"),
			 labels = c("Pre-treat\n-ment", "wk1", "wk2", "wk3")) +
	scale_y_log10(limits = c(5e-6, 50),
		      labels = scientific_10) +
	xlab("") +
	ylab("Mean AF (%)") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 9, margin = margin(t = 5)),
	      axis.text.y = element_text(size = 9),
	      strip.background = element_blank()) +
	guides(color = FALSE, shape = FALSE) +
	facet_wrap(~canonical_groups, nrow = 1, ncol = 4, scales = "free")

pdf(file = "../res/Raw_Mean_AF_by_Wk_Hclust.pdf", width = 3*4, height = 2.5*1.0)
print(plot_)
dev.off()

plot_ = smry_ %>%
	dplyr::mutate(canonical_groups = as.character(clust_)) %>%
	dplyr::mutate(canonical_groups = case_when(
			canonical_groups == "4" ~ "Fast/Slow",
			canonical_groups == "3" ~ "Fast/Slow",
			canonical_groups == "2" ~ "No clearance",
			canonical_groups == "1" ~ "Intermediate"
	)) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups"), variable.name = "week", value.name = "AF") %>%
	tidyr::drop_na() %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::select(patient_id_mskcc, canonical_groups, week, AF, hypoxia_resolution, neck_dissection_yes_no, composite_end_point) %>%
	dplyr::mutate(hypoxia_resolution = case_when(
			hypoxia_resolution == "persistent" ~ "Persistent",
			hypoxia_resolution == "resolved" ~ "Resolved",
			TRUE ~ "NA"
	)) %>%
	dplyr::mutate(neck_dissection_yes_no = case_when(
			neck_dissection_yes_no == "yes" ~ "Yes",
			neck_dissection_yes_no == "no" ~ "No",
			TRUE ~ "NA"
	)) %>%
	dplyr::mutate(composite_end_point = case_when(
			composite_end_point == "TRUE" ~ "1",
			composite_end_point == "FALSE" ~ "0",
			TRUE ~ "NA"
	)) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups", "week", "AF"), variable.name = "clinical_variable", value.name = "end_point") %>%
	readr::type_convert() %>%
	tidyr::drop_na() %>%
	dplyr::filter(grepl("raw", week)) %>%
	dplyr::filter(week != "rel_wk0") %>%
	dplyr::mutate(clinical_variable = case_when(
		clinical_variable == "hypoxia_resolution" ~ "Hypoxia resolution",
		clinical_variable == "neck_dissection_yes_no" ~ "Neck dissection",
		clinical_variable == "composite_end_point" ~ "Composite end-point"
	)) %>%
	dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast/Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
	dplyr::mutate(week = case_when(
		week == "raw_wk0" ~ "Pre-treatment",
		week == "raw_wk1" ~ "Week 1",
		week == "raw_wk2" ~ "Week 2",
		week == "raw_wk3" ~ "Week 3",
	)) %>%
	dplyr::mutate(week = factor(week, levels = c("Pre-treatment", "Week 1", "Week 2", "Week 3"), ordered = TRUE)) %>%
	dplyr::mutate(end_point = factor(end_point)) %>%
	dplyr::filter(clinical_variable == "Composite end-point") %>%
	ggplot(aes(x = end_point, y = 100*AF, color = end_point)) +
	geom_boxplot(stat = "boxplot", width = .7, outlier.shape = NA, fill = "white") +
	geom_jitter(stat = "identity", shape = 21, fill = "white", width = .15, height = 0, size = 2.5) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_x_discrete() +
	scale_y_log10(limits = c(0.001, 500),
		      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
		      labels = c(".001", ".01", ".1", "1", "10", "100")) +
	xlab("") +
	ylab("Mean AF (%)") +
	geom_signif(stat = "signif",
		    comparisons = list(c("0", "1")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided"),
		    size = .5, textsize = 3, color = "grey25", vjust = -1, 
		    y_position = log10(100)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 9, margin = margin(t = 5)),
	      axis.text.y = element_text(size = 9),
	      strip.background = element_rect(color = "white", fill = "white")) +
	guides(color = FALSE, shape = FALSE) +
	facet_wrap(~week, nrow = 1, scales = "free_y")
	
pdf(file = "../res/Raw_Mean_AF_by_Wk_by_Clinical_End_Points.pdf", width = 3*1.75, height = 3)
print(plot_)
dev.off()

plot_ = smry_ %>%
	dplyr::mutate(canonical_groups = as.character(clust_)) %>%
	dplyr::mutate(canonical_groups = case_when(
			canonical_groups == "4" ~ "Fast/Slow",
			canonical_groups == "3" ~ "Fast/Slow",
			canonical_groups == "2" ~ "No clearance",
			canonical_groups == "1" ~ "Intermediate"
	)) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups"), variable.name = "week", value.name = "AF") %>%
	tidyr::drop_na() %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::select(patient_id_mskcc, canonical_groups, week, AF, hypoxia_resolution, neck_dissection_yes_no, composite_end_point) %>%
	dplyr::mutate(hypoxia_resolution = case_when(
			hypoxia_resolution == "persistent" ~ "Persistent",
			hypoxia_resolution == "resolved" ~ "Resolved",
			TRUE ~ "NA"
	)) %>%
	dplyr::mutate(neck_dissection_yes_no = case_when(
			neck_dissection_yes_no == "yes" ~ "Yes",
			neck_dissection_yes_no == "no" ~ "No",
			TRUE ~ "NA"
	)) %>%
	dplyr::mutate(composite_end_point = case_when(
			composite_end_point == "TRUE" ~ "1",
			composite_end_point == "FALSE" ~ "0",
			TRUE ~ "NA"
	)) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups", "week", "AF"), variable.name = "clinical_variable", value.name = "end_point") %>%
	readr::type_convert() %>%
	tidyr::drop_na() %>%
	dplyr::filter(grepl("raw", week)) %>%
	dplyr::filter(week != "rel_wk0") %>%
	dplyr::mutate(clinical_variable = case_when(
		clinical_variable == "hypoxia_resolution" ~ "Hypoxia resolution",
		clinical_variable == "neck_dissection_yes_no" ~ "Neck dissection",
		clinical_variable == "composite_end_point" ~ "Composite end-point"
	)) %>%
	dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast/Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
	dplyr::mutate(week = case_when(
		week == "raw_wk0" ~ "Pre-treatment",
		week == "raw_wk1" ~ "Week 1",
		week == "raw_wk2" ~ "Week 2",
		week == "raw_wk3" ~ "Week 3",
	)) %>%
	dplyr::mutate(week = factor(week, levels = c("Pre-treatment", "Week 1", "Week 2", "Week 3"), ordered = TRUE)) %>%
	dplyr::mutate(end_point = factor(end_point)) %>%
	dplyr::filter(clinical_variable == "Composite end-point") %>%
	ggplot(aes(x = end_point:canonical_groups, y = 100*AF, color = canonical_groups)) +
	geom_boxplot(stat = "boxplot", width = .85, outlier.shape = NA, fill = "white") +
	geom_jitter(stat = "identity", shape = 21, fill = "white", width = .15, height = 0, size = 2.15) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_x_discrete() +
	scale_y_log10(limits = c(0.001, 500),
		      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
		      labels = c(".001", ".01", ".1", "1", "10", "100")) +
	xlab("") +
	ylab("Mean AF (%)") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 9, margin = margin(t = 5), angle = 90, vjust = 0.5, hjust = 1),
	      axis.text.y = element_text(size = 9),
	      strip.background = element_rect(color = "white", fill = "white")) +
	guides(color = FALSE, shape = FALSE) +
	facet_wrap(~week, nrow = 1, scales = "free_y")
	
pdf(file = "../res/Raw_Mean_AF_by_Wk_by_Clinical_End_Points_by_ctDNA_Clearance.pdf", width = 3*2.2, height = 3*1.3)
print(plot_)
dev.off()

#==================================================
# Absolute MRI volume wk1, wk2, wk3, wk5
#==================================================
plot_ = smry_ %>%
	dplyr::bind_cols(tsne_$Y %>%
			 dplyr::as_tibble() %>%
			 dplyr::select(tSNE1 = V1, tSNE2 = V2)) %>%
	dplyr::mutate(canonical_groups = as.character(clust_)) %>%
	dplyr::mutate(canonical_groups = case_when(
		canonical_groups == "4" ~ "Fast",
		canonical_groups == "3" ~ "Slow",
		canonical_groups == "2" ~ "No clearance",
		canonical_groups == "1" ~ "Intermediate"
	)) %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::select(patient_id_mskcc, canonical_groups, MRI_rawdata_wk0, MRI_rawdata_wk2, MRI_rawdata_wk3, MRI_rawdata_wk4) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "canonical_groups"), variable.name = "week", value.name = "MRI volume") %>%
	dplyr::mutate(canonical_groups = factor(canonical_groups, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE)) %>%
	dplyr::mutate(week = gsub(pattern = "MRI_rawdata_wk", replacement = "", x = week)) %>%
	ggplot(aes(x = week, y = `MRI volume`, color = canonical_groups, shape = canonical_groups)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "black") +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", size = 3, alpha = .95) +
	scale_shape_manual(values = c("Fast" = 21, "Slow" = 22, "Intermediate" = 23, "No clearance" = 24)) +
	scale_color_manual(values = c("Fast" = "#1b9e77", "Slow" = "#d95f02", "Intermediate" = "#7570b3", "No clearance" = "#e7298a")) +
	scale_x_discrete(breaks = c("0", "2", "3", "4"),
			 labels = c("0", "1", "2", "3")) +
	scale_y_continuous(limits = c(0, 5E4),
			   labels = scientific_10) +
	xlab("") +
	ylab("MRI Volume") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 9.5),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	guides(color = FALSE, shape = FALSE) +
	facet_wrap(~canonical_groups, nrow = 1)

pdf(file = "../res/ctDNA_Clearance_MRI_Volume_All_weeks.pdf", width = 2.85*2.5, height = 2.5*1.5)
print(plot_)
dev.off()
