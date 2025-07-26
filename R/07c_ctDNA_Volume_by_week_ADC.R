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


smry_pcm = idx_metrics_ft %>%
	   dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	   dplyr::summarize(mean_af = mean(mean_af+1E-5, na.rm = TRUE)) %>%
	   reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			   fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	   dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	   readr::type_convert()

smry_hpv = idx_metrics_ft %>%
	   dplyr::filter(chromosome == "HPV-16") %>%
	   dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	   dplyr::summarize(aligned_reads = mean(aligned_reads+1, na.rm = TRUE)) %>%
	   reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			   fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "aligned_reads") %>%
	   dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	   readr::type_convert()

smry_mri = clinical %>%
	   dplyr::select(patient_id_mskcc,
		        `Pre-treatment` = MRI_rawdata_wk0,
		        `wk1` = MRI_rawdata_wk2,
		        `wk2` = MRI_rawdata_wk3,
		        `wk3` = MRI_rawdata_wk4) %>%
	   readr::type_convert()

smry_adc = clinical %>%
	   dplyr::select(patient_id_mskcc,
		        `Pre-treatment` = ADC_Mean_wk0,
		        `wk1` = ADC_Mean_wk2,
		        `wk2` = ADC_Mean_wk3,
		        `wk3` = ADC_Mean_wk4) %>%
	   readr::type_convert()

plot_ = smry_pcm %>%
	reshape2::melt(variable.name = "week", value.name = "AF") %>%
	dplyr::full_join(smry_hpv %>%
			 reshape2::melt(variable.name = "week", value.name = "HPV"),
			 by = c("patient_id_mskcc", "week")) %>%
	tidyr::drop_na() %>%
	dplyr::filter(week != "Pre-treatment") %>%
	dplyr::mutate(week = gsub(pattern = "wk", replacement = "Week ", x = week, fixed = TRUE)) %>%
	ggplot(aes(x = AF, y = HPV, shape = week, color = week)) +
	geom_point(stat = "identity", fill = "white", alpha = 1, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x,
		    se = FALSE, fullrange = TRUE, color = "goldenrod3", size = .85) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22, 23, 24)) +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(limits = c(1E0, 1E7),
		      breaks = c(1E0, 1E2, 1E4, 1E6),
		      labels = scientific_10) +
	xlab("ctDNA Fraction (%)") +
	ylab("cfDNA HPV Aligned Read Pairs") +
	stat_cor(method = "spearman", color = "black", label.x = log10(1E-5), label.y = log10(5E6)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(color = FALSE, shape = FALSE)

pdf(file = "../res/Correlation_AF_HPV_by_week.pdf", width = 3*2.75, height = 2.95*1.0)
print(plot_)
dev.off()

plot_ = smry_pcm %>%
	reshape2::melt(variable.name = "week", value.name = "AF") %>%
	dplyr::full_join(smry_mri %>%
			 reshape2::melt(variable.name = "week", value.name = "MRI"),
			 by = c("patient_id_mskcc", "week")) %>%
	tidyr::drop_na() %>%
	dplyr::filter(week != "Pre-treatment") %>%
	dplyr::mutate(week = gsub(pattern = "wk", replacement = "Week ", x = week, fixed = TRUE)) %>%
	ggplot(aes(x = AF, y = MRI/1000, shape = week, color = week)) +
	geom_point(stat = "identity", fill = "white", alpha = 1, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x,
		    se = FALSE, fullrange = TRUE, color = "goldenrod3", size = .85) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22, 23, 24)) +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(limits = c(1E3, 7E4)/1000,
		      labels = scientific_10) +
	xlab("ctDNA Fraction (%)") +
	ylab(expression("MRI Volume "(cm^3))) +
	stat_cor(method = "spearman", color = "black", label.x = log10(1E-5), label.y = log10(6E4/1000)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(color = FALSE, shape = FALSE)

pdf(file = "../res/Correlation_AF_MRI_Vol_by_week.pdf", width = 3*2.75, height = 2.95*1.0)
print(plot_)
dev.off()

plot_ = smry_pcm %>%
	reshape2::melt(variable.name = "week", value.name = "AF") %>%
	dplyr::full_join(smry_adc %>%
			 reshape2::melt(variable.name = "week", value.name = "ADC"),
			 by = c("patient_id_mskcc", "week")) %>%
	tidyr::drop_na() %>%
	dplyr::filter(week != "Pre-treatment") %>%
	dplyr::mutate(week = gsub(pattern = "wk", replacement = "Week ", x = week, fixed = TRUE)) %>%
	ggplot(aes(x = AF, y = ADC, shape = week, color = week)) +
	geom_point(stat = "identity", fill = "white", alpha = 1, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x,
		    se = FALSE, fullrange = TRUE, color = "goldenrod3", size = .85) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22, 23, 24)) +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(limits = c(0.5, 2.5),
		      labels = scientific_10) +
	xlab("ctDNA Fraction (%)") +
	ylab(expression("ADC ("%.%10^-3~mm^2~sec^-1~")")) +
	stat_cor(method = "spearman", color = "black", label.x = log10(1E-5), label.y = log10(2.5)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(color = FALSE, shape = FALSE)

pdf(file = "../res/Correlation_AF_ADC_Vol_by_week.pdf", width = 3*2.75, height = 2.95*1.0)
print(plot_)
dev.off()

plot_ = smry_hpv %>%
	reshape2::melt(variable.name = "week", value.name = "HPV") %>%
	dplyr::full_join(smry_mri %>%
			 reshape2::melt(variable.name = "week", value.name = "MRI"),
			 by = c("patient_id_mskcc", "week")) %>%
	tidyr::drop_na() %>%
	dplyr::filter(week != "Pre-treatment") %>%
	dplyr::mutate(week = gsub(pattern = "wk", replacement = "Week ", x = week, fixed = TRUE)) %>%
	ggplot(aes(x = HPV, y = MRI/1000, shape = week, color = week)) +
	geom_point(stat = "identity", fill = "white", alpha = 1, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x,
		    se = FALSE, fullrange = TRUE, color = "goldenrod3", size = .85) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22, 23, 24)) +
	scale_x_log10(limits = c(1E1, 1E7),
		      labels = scientific_10) +
	scale_y_log10(limits = c(1E3, 7E4)/1000,
		      labels = scientific_10) +
	xlab("cfDNA HPV Aligned Read Pairs") +
	ylab(expression("MRI Volume "(cm^3))) +
	stat_cor(method = "spearman", color = "black", label.x = log10(1E1), label.y = log10(6E4/1000)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(color = FALSE, shape = FALSE)

pdf(file = "../res/Correlation_HPV_MRI_Vol_by_week.pdf", width = 3*2.75, height = 2.95*1.0)
print(plot_)
dev.off()

plot_ = smry_adc %>%
	reshape2::melt(variable.name = "week", value.name = "ADC") %>%
	dplyr::full_join(smry_mri %>%
			 reshape2::melt(variable.name = "week", value.name = "MRI"),
			 by = c("patient_id_mskcc", "week")) %>%
	tidyr::drop_na() %>%
	dplyr::filter(week != "Pre-treatment") %>%
	dplyr::mutate(week = gsub(pattern = "wk", replacement = "Week ", x = week, fixed = TRUE)) %>%
	ggplot(aes(x = ADC, y = MRI/1000, shape = week, color = week)) +
	geom_point(stat = "identity", fill = "white", alpha = 1, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x,
		    se = FALSE, fullrange = TRUE, color = "goldenrod3", size = .85) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22, 23, 24)) +
	scale_x_log10(limits = c(0.5, 2.5)) +
	scale_y_log10(limits = c(1E3, 7E4)/1000,
		      labels = scientific_10) +
	xlab(expression("ADC ("%.%10^-3~mm^2~sec^-1~")")) +
	ylab(expression("MRI Volume "(cm^3))) +
	stat_cor(method = "spearman", color = "black", label.x = log10(0.55), label.y = log10(6E4/1000)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 10),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(color = FALSE, shape = FALSE)

pdf(file = "../res/Correlation_ADC_MRI_Vol_by_week.pdf", width = 3*2.75, height = 2.95*1.0)
print(plot_)
dev.off()

fit_ctdna = smry_pcm %>%
	    reshape2::melt(variable.name = "week", value.name = "ctDNA") %>%
	    dplyr::full_join(smry_mri %>%
			     reshape2::melt(variable.name = "week", value.name = "Volume"),
			     by = c("patient_id_mskcc", "week")) %>%
	    dplyr::full_join(smry_adc %>%
			     reshape2::melt(variable.name = "week", value.name = "ADC"),
			     by = c("patient_id_mskcc", "week")) %>%
	    tidyr::drop_na() %>%
	    dplyr::mutate(Volume = scale(Volume, center = TRUE, scale = TRUE)) %>%
	    stats::glm(formula = ctDNA ~ Volume + ADC:week, data = .)

fit_hpv = smry_hpv %>%
	  reshape2::melt(variable.name = "week", value.name = "HPV") %>%
	  dplyr::full_join(smry_mri %>%
			   reshape2::melt(variable.name = "week", value.name = "Volume"),
			   by = c("patient_id_mskcc", "week")) %>%
	  dplyr::full_join(smry_adc %>%
			   reshape2::melt(variable.name = "week", value.name = "ADC"),
			   by = c("patient_id_mskcc", "week")) %>%
	  dplyr::mutate(HPV = scale(HPV, center = TRUE, scale = TRUE)) %>%
	  dplyr::mutate(Volume = scale(Volume, center = TRUE, scale = TRUE)) %>%
	  stats::glm(formula = HPV ~ Volume + ADC:week, data = .)

plot_ = summary(fit_ctdna) %>%
	.[["coefficients"]] %>%
	as.data.frame() %>%
	tibble::rownames_to_column("variable") %>%
	dplyr::as_tibble() %>%
	dplyr::filter(variable != "(Intercept)") %>%
	dplyr::mutate(assay = "PCM") %>%
	dplyr::mutate(Estimate = Estimate*20) %>%
	dplyr::bind_rows(summary(fit_hpv) %>%
			 .[["coefficients"]] %>%
			 as.data.frame() %>%
			 tibble::rownames_to_column("variable") %>%
			 dplyr::as_tibble() %>%
			 dplyr::filter(variable != "(Intercept)") %>%
			 dplyr::mutate(assay = "HPV")) %>%
	dplyr::mutate(is_significant = ifelse(`Pr(>|t|)`<.1, "Yes", "No")) %>%
	dplyr::mutate(variable = gsub(pattern = "week", replacement = "", x= variable, fixed = TRUE)) %>%
	dplyr::mutate(variable = gsub(pattern = "wk", replacement = "Week ", x= variable, fixed = TRUE)) %>%
	dplyr::mutate(variable = factor(variable, levels = rev(unique(variable)), ordered = TRUE)) %>%
	ggplot(aes(x = variable, ymin = 0, ymax = Estimate, y = Estimate, fill = is_significant, size = -log10(`Pr(>|t|)`))) +
	geom_linerange(stat = "identity", size = .5) +
	geom_hline(yintercept = 0, size = 1) +
	geom_point(stat = "identity", shape = 21) +
	xlab("") +
	ylab("Standardized Coefficients") +
	scale_fill_manual(values = c("#bdbdbd", "#e41a1c")) +
	scale_size_continuous(breaks = c(1,3,5)) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(-.75, .35)) +
	coord_flip() +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 8),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank(),
	      panel.spacing = unit(2, "lines")) +
	guides(fill = guide_legend(title = expression(p<0.1), order = 1),
	       size = guide_legend(title = expression(-Log[10]~"p-value"), override.aes = list(shape = 21, fill = "black"))) +
	facet_wrap(~assay, ncol = 2, scales = "free_x")

pdf(file = "../res/Linear_Regression_Coefficients_ADC.pdf", width = 6.5, height = 3)
print(plot_)
dev.off()
