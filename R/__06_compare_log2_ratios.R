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
		 
aligned_reads = idx_metrics_ft %>%
		dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
		dplyr::summarize(aligned_reads = mean(aligned_reads, na.rm=TRUE),
				 hpv_type_wes_wgs = unique(hpv_type_wes_wgs)) %>%
	        dplyr::ungroup() %>%
		reshape2::dcast(patient_id_mskcc + hpv_type_wes_wgs ~ timepoint_weeks_since_start_of_RT, value.var = "aligned_reads", fill = NA)

mean_af = idx_metrics_ft %>%
	  dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	  dplyr::summarize(mean_af = mean(mean_af, na.rm=TRUE),
	  		   hpv_type_wes_wgs = unique(hpv_type_wes_wgs)) %>%
	  dplyr::ungroup() %>%
	  reshape2::dcast(patient_id_mskcc + hpv_type_wes_wgs ~ timepoint_weeks_since_start_of_RT, value.var = "mean_af", fill = NA)

mri_volume = clinical %>%
	     dplyr::select(patient_id_mskcc,
			   `Pre-treatment` = MRI_rawdata_wk0,
			   wk1 = MRI_rawdata_wk1,
			   wk2 = MRI_rawdata_wk2,
			   wk3 = MRI_rawdata_wk3,
			   wk5 = MRI_rawdata_wk4)

plot_ = aligned_reads %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::select(patient_id_mskcc, wk1, wk2, wk3, wk5) %>%
	reshape2::melt(variable.name = "Week", value.name = "aligned_reads") %>%
	dplyr::full_join(mean_af %>%
			 dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
			 dplyr::select(patient_id_mskcc, wk1, wk2, wk3, wk5) %>%
			 reshape2::melt(variable.name = "Week", value.name = "mean_af"),
			 by = c("patient_id_mskcc", "Week")) %>%
	tidyr::drop_na() %>%
	ggplot(aes(x = (mean_af*100)+(1E-3), y = aligned_reads)) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, color = "goldenrod3", size = 1.75, se = FALSE) +
	geom_point(stat = "identity", shape = 21, fill = "white", color = "salmon", alpha = .75, size = 2.5) +
	stat_cor(method = "spearman", size = 3) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_log10(limits = c(1e-3, 110),
		      labels = scientific_10) +
	scale_y_log10(limits = c(1e1, 1e6),
		      labels = scientific_10) +
	xlab("Mean PCM AF (%)") +
	ylab("cfDNA Aligned HPV Read Pairs") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      panel.spacing = unit(1, "lines"),
	      strip.background = element_rect(colour="white", fill="white")) +
	facet_wrap(~Week, nrow = 2, ncol = 2, scales = "free")

pdf(file = "../res/Total_Reads_Mean_AF_bywk.pdf", width = 5.25, height = 5.25)
print(plot_)
dev.off()
	
plot_ = aligned_reads %>%
	dplyr::mutate(delta_aligned_reads = log2(`wk0`/`Pre-treatment`)) %>%
	dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
	dplyr::full_join(mean_af %>%
			 dplyr::mutate(delta_mean_af = log2(`wk0`/`Pre-treatment`)) %>%
			 dplyr::select(patient_id_mskcc, delta_mean_af),
			 by = "patient_id_mskcc") %>%
	dplyr::mutate(delta = "wk0") %>%
	dplyr::bind_rows(
		aligned_reads %>%
		dplyr::mutate(delta_aligned_reads = log2(`wk1`/`Pre-treatment`)) %>%
		dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
		dplyr::full_join(mean_af %>%
				 dplyr::mutate(delta_mean_af = log2(`wk1`/`Pre-treatment`)) %>%
				 dplyr::select(patient_id_mskcc, delta_mean_af),
				 by = "patient_id_mskcc") %>%
		dplyr::mutate(delta = "wk1")) %>%
	dplyr::bind_rows(
		aligned_reads %>%
		dplyr::mutate(delta_aligned_reads = log2(`wk2`/`Pre-treatment`)) %>%
		dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
		dplyr::full_join(mean_af %>%
				 dplyr::mutate(delta_mean_af = log2(`wk2`/`Pre-treatment`)) %>%
				 dplyr::select(patient_id_mskcc, delta_mean_af),
				 by = "patient_id_mskcc") %>%
		dplyr::mutate(delta = "wk2")) %>%
	dplyr::bind_rows(
		aligned_reads %>%
		dplyr::mutate(delta_aligned_reads = log2(`wk3`/`Pre-treatment`)) %>%
		dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
		dplyr::full_join(mean_af %>%
				 dplyr::mutate(delta_mean_af = log2(`wk3`/`Pre-treatment`)) %>%
				 dplyr::select(patient_id_mskcc, delta_mean_af),
				 by = "patient_id_mskcc") %>%
		dplyr::mutate(delta = "wk3")) %>%
	dplyr::bind_rows(
		aligned_reads %>%
		dplyr::mutate(delta_aligned_reads = log2(`wk5`/`Pre-treatment`)) %>%
		dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
		dplyr::full_join(mean_af %>%
				 dplyr::mutate(delta_mean_af = log2(`wk5`/`Pre-treatment`)) %>%
				 dplyr::select(patient_id_mskcc, delta_mean_af),
				 by = "patient_id_mskcc") %>%
		dplyr::mutate(delta = "wk5")) %>%
	dplyr::bind_rows(
		aligned_reads %>%
		dplyr::mutate(delta_aligned_reads = log2(`wk6`/`Pre-treatment`)) %>%
		dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
		dplyr::full_join(mean_af %>%
				 dplyr::mutate(delta_mean_af = log2(`wk6`/`Pre-treatment`)) %>%
				 dplyr::select(patient_id_mskcc, delta_mean_af),
				 by = "patient_id_mskcc") %>%
		dplyr::mutate(delta = "wk6")) %>%
	dplyr::bind_rows(
		aligned_reads %>%
		dplyr::mutate(delta_aligned_reads = log2(`wk7+`/`Pre-treatment`)) %>%
		dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
		dplyr::full_join(mean_af %>%
				 dplyr::mutate(delta_mean_af = log2(`wk7+`/`Pre-treatment`)) %>%
				 dplyr::select(patient_id_mskcc, delta_mean_af),
				 by = "patient_id_mskcc") %>%
		dplyr::mutate(delta = "wk7+")) %>%
	tidyr::drop_na() %>%
	dplyr::filter(delta %in% c("wk1", "wk2", "wk3", "wk5")) %>%
	ggplot(aes(x = delta_mean_af, y = delta_aligned_reads)) +
	geom_abline(intercept = 0, slope = 1, color = "lightgrey", alpha = .75, size = 1) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, color = "goldenrod3", size = 1.75, se = FALSE) +
	geom_point(stat = "identity", shape = 21, fill = "white", color = "salmon", alpha = .75, size = 2.5) +
	stat_cor(method = "spearman", size = 3) +
	scale_x_continuous(limits = c(-15, 10)) +
	scale_y_continuous(limits = c(-15, 10)) +
	xlab(expression(Log[2]~"Ratio Mean PCM AF")) +
	ylab(expression(Log[2]~"Ratio cfDNA Aligned HPV Reads Pairs")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      panel.spacing = unit(1, "lines"),
	      strip.background = element_rect(colour="white", fill="white")) +
	facet_wrap(~delta, nrow = 2, ncol = 2, scales = "free")

pdf(file = "../res/Log2_Ratio_Total_Reads_Mean_AF_bywk.pdf", width = 5.25, height = 5.25)
print(plot_)
dev.off()

plot_ = aligned_reads %>%
	dplyr::mutate(delta_aligned_reads = log2(`wk1`/`Pre-treatment`)) %>%
	dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
	dplyr::full_join(mri_volume %>%
			 dplyr::mutate(delta_mri = log2(`wk1`/`Pre-treatment`)) %>% 
			 dplyr::select(patient_id_mskcc, delta_mri)) %>%
	dplyr::mutate(delta = "wk1") %>%
	dplyr::bind_rows(aligned_reads %>%
			 dplyr::mutate(delta_aligned_reads = log2(`wk2`/`Pre-treatment`)) %>%
			 dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
			 dplyr::full_join(mri_volume %>%
					  dplyr::mutate(delta_mri = log2(`wk2`/`Pre-treatment`)) %>% 
					  dplyr::select(patient_id_mskcc, delta_mri)) %>%
			 dplyr::mutate(delta = "wk2")) %>%
	dplyr::bind_rows(aligned_reads %>%
			 dplyr::mutate(delta_aligned_reads = log2(`wk3`/`Pre-treatment`)) %>%
			 dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
			 dplyr::full_join(mri_volume %>%
					  dplyr::mutate(delta_mri = log2(`wk3`/`Pre-treatment`)) %>% 
					  dplyr::select(patient_id_mskcc, delta_mri)) %>%
			 dplyr::mutate(delta = "wk3")) %>%
	dplyr::bind_rows(aligned_reads %>%
			 dplyr::mutate(delta_aligned_reads = log2(`wk5`/`Pre-treatment`)) %>%
			 dplyr::select(patient_id_mskcc, delta_aligned_reads) %>%
			 dplyr::full_join(mri_volume %>%
					  dplyr::mutate(delta_mri = log2(`wk5`/`Pre-treatment`)) %>% 
					  dplyr::select(patient_id_mskcc, delta_mri)) %>%
			 dplyr::mutate(delta = "wk5")) %>%
	tidyr::drop_na() %>%
	dplyr::filter(delta %in% c("wk1", "wk2", "wk3", "wk5")) %>%
	ggplot(aes(x = delta_mri, y = delta_aligned_reads)) +
	geom_abline(intercept = 0, slope = 1, color = "lightgrey", alpha = .75, size = 1) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, color = "goldenrod3", size = 1.75, se = FALSE) +
	geom_point(stat = "identity", shape = 21, fill = "white", color = "salmon", alpha = .75, size = 2.5) +
	stat_cor(method = "spearman", size = 3) +
	scale_x_continuous(limits = c(-2, 1)) +
	scale_y_continuous(limits = c(-15, 10)) +
	xlab(expression(Log[2]~"Ratio MRI Volume")) +
	ylab(expression(Log[2]~"Ratio cfDNA Aligned HPV Reads Pairs")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      panel.spacing = unit(1, "lines"),
	      strip.background = element_rect(colour="white", fill="white")) +
	facet_wrap(~delta, nrow = 2, ncol = 2, scales = "free")

pdf(file = "../res/Log2_Ratio_Total_Reads_MRI_Volume_bywk.pdf", width = 5.25, height = 5.25)
print(plot_)
dev.off()

plot_ = mean_af %>%
	dplyr::mutate(delta_mean_af = log2(`wk1`/`Pre-treatment`)) %>%
	dplyr::select(patient_id_mskcc, delta_mean_af) %>%
	dplyr::full_join(mri_volume %>%
			 dplyr::mutate(delta_mri = log2(wk1/`Pre-treatment`)) %>% 
			 dplyr::select(patient_id_mskcc, delta_mri)) %>%
	dplyr::mutate(delta = "wk1") %>%
	dplyr::bind_rows(mean_af %>%
			 dplyr::mutate(delta_mean_af = log2(`wk2`/`Pre-treatment`)) %>%
			 dplyr::select(patient_id_mskcc, delta_mean_af) %>%
			 dplyr::full_join(mri_volume %>%
					  dplyr::mutate(delta_mri = log2(`wk2`/`Pre-treatment`)) %>% 
					  dplyr::select(patient_id_mskcc, delta_mri)) %>%
			 dplyr::mutate(delta = "wk2")) %>%
	dplyr::bind_rows(mean_af %>%
			 dplyr::mutate(delta_mean_af = log2(`wk3`/`Pre-treatment`)) %>%
			 dplyr::select(patient_id_mskcc, delta_mean_af) %>%
			 dplyr::full_join(mri_volume %>%
					  dplyr::mutate(delta_mri = log2(`wk3`/`Pre-treatment`)) %>% 
					  dplyr::select(patient_id_mskcc, delta_mri)) %>%
			 dplyr::mutate(delta = "wk3")) %>%
	dplyr::bind_rows(mean_af %>%
			 dplyr::mutate(delta_mean_af = log2(`wk5`/`Pre-treatment`)) %>%
			 dplyr::select(patient_id_mskcc, delta_mean_af) %>%
			 dplyr::full_join(mri_volume %>%
					  dplyr::mutate(delta_mri = log2(`wk5`/`Pre-treatment`)) %>% 
					  dplyr::select(patient_id_mskcc, delta_mri)) %>%
			 dplyr::mutate(delta = "wk5")) %>%
	tidyr::drop_na() %>%
	dplyr::filter(delta %in% c("wk1", "wk2", "wk3", "wk5")) %>%
	ggplot(aes(x = delta_mri, y = delta_mean_af)) +
	geom_abline(intercept = 0, slope = 1, color = "lightgrey", alpha = .75, size = 1) +
	geom_smooth(stat = "smooth", method = "lm", formula = y ~ x, color = "goldenrod3", size = 1.75, se = FALSE) +
	geom_point(stat = "identity", shape = 21, fill = "white", color = "salmon", alpha = .75, size = 2.5) +
	stat_cor(method = "spearman", size = 3) +
	scale_x_continuous(limits = c(-2, 1)) +
	scale_y_continuous(limits = c(-15, 10)) +
	xlab(expression(Log[2]~"Ratio MRI Volume")) +
	ylab(expression(Log[2]~"Ratio Mean PCM AF")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      panel.spacing = unit(1, "lines"),
	      strip.background = element_rect(colour="white", fill="white")) +
	facet_wrap(~delta, nrow = 2, ncol = 2, scales = "free")

pdf(file = "../res/Log2_Ratio_Mean_AF_MRI_Volume_bywk.pdf", width = 5.25, height = 5.25)
print(plot_)
dev.off()

plot_ = mean_af %>%
	dplyr::mutate(wk0 = log2(wk0/`Pre-treatment`),
		      wk1 = log2(wk1/`Pre-treatment`),
		      wk2 = log2(wk2/`Pre-treatment`),
		      wk3 = log2(wk3/`Pre-treatment`),
		      wk4 = log2(wk4/`Pre-treatment`),
		      wk5 = log2(wk5/`Pre-treatment`),
		      wk6 = log2(wk6/`Pre-treatment`),
		      `wk7+` = log2(`wk7+`/`Pre-treatment`)) %>%
	dplyr::select(patient_id_mskcc, wk0, wk1, wk2, wk3, -wk4, wk5, wk6, `wk7+`) %>%
	readr::type_convert() %>%
	reshape2::melt(id.vars = "patient_id_mskcc") %>%
	tidyr::drop_na() %>%
	dplyr::mutate(value = case_when(
		is.infinite(value) ~ -10,
		TRUE ~ value
	)) %>%
	dplyr::filter(variable!= "wk0") %>%
	dplyr::filter(variable!= "wk7+") %>%
	ggplot(aes(x = variable, y = value)) +
	geom_boxplot(stat = "boxplot", position = "dodge2", outlier.shape = NA, color = "grey5", fill = "white") +
	geom_beeswarm(shape = 21, dodge.width = 0.75, fill = "white", alpha = .75, size = 2.5, cex = 1) +
	geom_hline(yintercept = 0, color = "goldenrod3", size = 1, linetype = 3, alpha = .85) +
	scale_x_discrete(breaks = c("wk1", "wk2", "wk3", "wk4", "wk5", "wk6"),
			 labels = c("wk1", "wk2", "wk3", "wk4", "wk5", "wk6")) +
	scale_y_continuous(limits = c(-15, 10)) +
	xlab("Weeks after CRT") +
	ylab(expression(Log[2]~"Ratio Mean PCM AF")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Log2_Ratio_Mean_AF_bywk.pdf", width = 5.25, height = 5.25)
print(plot_)
dev.off()

data_ = aligned_reads %>%
	dplyr::mutate(wk0 = log2(wk0/`Pre-treatment`),
		      wk1 = log2(wk1/`Pre-treatment`),
		      wk2 = log2(wk2/`Pre-treatment`),
		      wk3 = log2(wk3/`Pre-treatment`),
		      wk4 = log2(wk4/`Pre-treatment`),
		      wk5 = log2(wk5/`Pre-treatment`),
		      wk6 = log2(wk6/`Pre-treatment`),
		      wk7 = log2(`wk7+`/`Pre-treatment`)) %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	readr::type_convert() %>%
        dplyr::select(wk0, wk1, wk2, wk3, -wk4, wk5, wk6, `wk7`,
		      crt_randomization, primary_tumor_size_cm, neck_dissection_yes_no,
		      age, sex, t_stage, n_stage, smoking_category_yes_never, simplified_hypoxia_group,
		      plan_volume, hypoxia_resolution) %>%
       readr::type_convert()

fit_ = list()

fit_[[1]] = data_ %>%
	    dplyr::select(-wk1, -wk2, -wk3, -wk5, -wk6, -wk7) %>%
	    tidyr::drop_na() %>%
	    dplyr::filter(!is.infinite(wk0)) %>%
	    lm(formula = wk0 ~ ., data = .) %>%
	    summary()

fit_[[2]] = data_ %>%
	    dplyr::select(-wk0, -wk2, -wk3, -wk5, -wk6, -wk7) %>%
	    tidyr::drop_na() %>%
	    dplyr::filter(!is.infinite(wk1)) %>%
	    lm(formula = wk1 ~ ., data = .) %>%
	    summary()

fit_[[3]] = data_ %>%
	    dplyr::select(-wk0, -wk1, -wk3, -wk5, -wk6, -wk7) %>%
	    tidyr::drop_na() %>%
	    dplyr::filter(!is.infinite(wk2)) %>%
	    lm(formula = wk2 ~ ., data = .) %>%
	    summary()

fit_[[4]] = data_ %>%
	    dplyr::select(-wk0, -wk1, -wk2, -wk5, -wk6, -wk7) %>%
	    tidyr::drop_na() %>%
	    dplyr::filter(!is.infinite(wk3)) %>%
	    lm(formula = wk3 ~ ., data = .) %>%
	    summary()

fit_[[5]] = data_ %>%
	    dplyr::select(-wk0, -wk1, -wk2, -wk3, -wk6, -wk7) %>%
	    tidyr::drop_na() %>%
	    dplyr::filter(!is.infinite(wk5)) %>%
	    lm(formula = wk5 ~ ., data = .) %>%
	    summary()

fit_[[6]] = data_ %>%
	    dplyr::select(-wk0, -wk1, -wk2, -wk3, -wk5, -wk7) %>%
	    tidyr::drop_na() %>%
	    dplyr::filter(!is.infinite(wk6)) %>%
	    lm(formula = wk6 ~ ., data = .) %>%
	    summary()

fit_[[7]] = data_ %>%
	    dplyr::select(-wk0, -wk1, -wk2, -wk3, -wk5, -wk6) %>%
	    tidyr::drop_na() %>%
	    dplyr::filter(!is.infinite(wk7)) %>%
	    lm(formula = wk7 ~ ., data = .) %>%
	    summary()

p_values = do.call(cbind, lapply(fit_, function(x) { (x$coefficients[,"Pr(>|t|)"])})) %>%
	   as.data.frame() %>%
	   tibble::rownames_to_column("variable") %>%
	   dplyr::as_tibble() %>%
	   dplyr::filter(variable != "(Intercept)") %>%
	   dplyr::rename(wk0 = V1, wk1 = V2, wk2 = V3, wk3 = V4, wk5 = V5, wk6 = V6, wk7 = V7) %>%
	   reshape2::melt(variable.name = "week", value.name = "p_value") %>%
	   dplyr::mutate(fdr = p.adjust(p_value, "fdr"))

pdf(file = "../res/Linear_Regression_Coefficients_Log2_Ratio_bywk.pdf", width = 4.15, height = 3.5)
draw(Heatmap(matrix = p_values %>%
	     	      reshape2::dcast(variable ~ week, value.var = "fdr") %>%
	     	      dplyr::select(-wk0, -wk5, -wk7) %>%
	     	      tibble::column_to_rownames("variable") %>%
	     	      as.matrix(),
	     col = viridis(n = 10),
	     name = "P-value",
	     rect_gp = gpar(col = "white", lwd = 1),
	     border = NA,
	     
	     cluster_rows = TRUE,
	     clustering_distance_rows = "euclidean",
	     clustering_method_rows = "complete",
	     show_row_names = TRUE,
	     row_names_side = "right",
	     row_names_gp = gpar(fontsize = 8),
	     
	     cluster_columns = FALSE,
	     show_column_names = TRUE,
	     column_title_side = "bottom",
	     column_title_gp = gpar(fontsize = 8),

	     use_raster = FALSE,
	     show_heatmap_legend = TRUE))
dev.off()

plot_ = data_ %>%
	dplyr::select(neck_dissection_yes_no, wk6) %>%
	tidyr::drop_na() %>%
	dplyr::filter(!is.infinite(wk6)) %>%
	ggplot(aes(x = neck_dissection_yes_no, y = wk6)) +
	geom_boxplot(stat = "boxplot", position = "dodge2", outlier.shape = NA, color = "grey5", fill = "white") +
	geom_beeswarm(shape = 21, dodge.width = 0.75, fill = "white", alpha = .75, size = 2.5, cex = 1) +
	geom_hline(yintercept = 0, color = "goldenrod3", size = 1, linetype = 3, alpha = .85) +
	scale_x_discrete() +
	scale_y_continuous() +
	xlab("Neck dissection") +
	ylab(expression(Log[2]~"Ratio Mean PCM AF")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

