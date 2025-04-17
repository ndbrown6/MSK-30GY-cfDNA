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
	dplyr::full_join(smry_mri %>%
			 reshape2::melt(variable.name = "week", value.name = "MRI"),
			 by = c("patient_id_mskcc", "week")) %>%
	tidyr::drop_na() %>%
	ggplot(aes(x = AF, y = MRI, shape = week, color = week)) +
	geom_point(stat = "identity", fill = "white", alpha = .75, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x,
		    se = FALSE, fullrange = TRUE, color = "goldenrod3", alpha = .15, size = 1.5) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22, 23, 24)) +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(limits = c(1E3, 4E4),
		      labels = scientific_10) +
	xlab("ctDNA Fraction (%)") +
	ylab(expression("MRI Volume "(mm^3))) +
	stat_cor(method = "spearman", color = "black", label.x = log10(1E-5), label.y = log10(3E4)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = ggplot2::margin(t = 20)),
	      axis.title.y = element_text(margin = ggplot2::margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(color = FALSE, shape = FALSE)

pdf(file = "../res/Correlation_AF_MRI_Vol_by_week.pdf", width = 13, height = 3.5)
print(plot_)
dev.off()

plot_ = smry_hpv %>%
	reshape2::melt(variable.name = "week", value.name = "HPV") %>%
	dplyr::full_join(smry_mri %>%
			 reshape2::melt(variable.name = "week", value.name = "MRI"),
			 by = c("patient_id_mskcc", "week")) %>%
	tidyr::drop_na() %>%
	ggplot(aes(x = HPV, y = MRI, shape = week, color = week)) +
	geom_point(stat = "identity", fill = "white", alpha = .75, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x,
		    se = FALSE, fullrange = TRUE, color = "goldenrod3", alpha = .15, size = 1.5) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22, 23, 24)) +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(limits = c(1E3, 4E4),
		      labels = scientific_10) +
	xlab("cfDNA HPV Aligned Read Pairs") +
	ylab(expression("MRI Volume "(mm^3))) +
	stat_cor(method = "spearman", color = "black", label.x = log10(1E0), label.y = log10(3E4)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = ggplot2::margin(t = 20)),
	      axis.title.y = element_text(margin = ggplot2::margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(color = FALSE, shape = FALSE)

pdf(file = "../res/Correlation_HPV_MRI_Vol_by_week.pdf", width = 13, height = 3.5)
print(plot_)
dev.off()

fit_ctDNA = smry_pcm %>%
	    reshape2::melt(variable.name = "week", value.name = "ctDNA") %>%
	    dplyr::full_join(smry_mri %>%
			     reshape2::melt(variable.name = "week", value.name = "Volume"),
			     by = c("patient_id_mskcc", "week")) %>%
	    dplyr::full_join(smry_adc %>%
			     reshape2::melt(variable.name = "week", value.name = "ADC"),
			     by = c("patient_id_mskcc", "week")) %>%
	    tidyr::drop_na() %>%
	    stats::glm(formula = ctDNA ~ Volume + ADC:week, data = .)

fit_HPV = smry_hpv %>%
	  reshape2::melt(variable.name = "week", value.name = "HPV") %>%
	  dplyr::full_join(smry_mri %>%
			   reshape2::melt(variable.name = "week", value.name = "Volume"),
			   by = c("patient_id_mskcc", "week")) %>%
	  dplyr::full_join(smry_adc %>%
			   reshape2::melt(variable.name = "week", value.name = "ADC"),
			   by = c("patient_id_mskcc", "week")) %>%
	  tidyr::drop_na() %>%
	  stats::glm(formula = HPV ~ Volume + ADC:week, data = .)

data_ = smry_pcm %>%
	reshape2::melt(variable.name = "week", value.name = "ctDNA") %>%
	dplyr::full_join(smry_mri %>%
			 reshape2::melt(variable.name = "week", value.name = "Volume"),
			 by = c("patient_id_mskcc", "week")) %>%
	tidyr::drop_na()

residuals = list()
weeks = c("Pre-treatment", "wk1", "wk2", "wk3")
for (i in 1:4) {
	residuals[[i]] = dplyr::tibble(residuals = data_ %>%
				       		   dplyr::filter(week == weeks[i]) %>%
				       		   stats::lm(formula = ctDNA ~ Volume, data = .) %>%
				       		   .[["residuals"]]) %>%
			 dplyr::mutate(week = weeks[i]) %>%
			 dplyr::mutate(patient_id_mskcc = data_ %>%
				       		   	  dplyr::filter(week == weeks[i]) %>%
				       			  .[["patient_id_mskcc"]])
}

plot_ = do.call(rbind, residuals) %>%
	dplyr::full_join(smry_adc %>%
			 reshape2::melt(variable.name = "week", value.name = "ADC"),
			 by = c("patient_id_mskcc", "week")) %>%
	tidyr::drop_na() %>%
	ggplot(aes(x = residuals, y = ADC, color = week, shape = week)) +
	geom_point(stat = "identity", fill = "white", alpha = .75, size = 2) +
	geom_smooth(stat = "smooth", method = "rlm", formula = y ~ x,
		    se = FALSE, fullrange = TRUE, color = "goldenrod3", alpha = .15, size = 1.5) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22, 23, 24)) +
	scale_x_continuous(labels = scientific_10) +
	scale_y_continuous(#limits = c(.6, 1.5),
			   labels = scientific_10) +
	xlab(expression("Residual MRI Volume "(mm^3))) +
	ylab("ADC") +
	stat_cor(method = "spearman", color = "black", label.x = -2, label.y = 1.5) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = ggplot2::margin(t = 20)),
	      axis.title.y = element_text(margin = ggplot2::margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	facet_wrap(~week, nrow = 1, scales = "free_y") +
	guides(color = FALSE, shape = FALSE)

#pdf(file = "../res/Correlation_MRI_Residuals_Vol_ADC_by_week.pdf", width = 11, height = 3.5)
#print(plot_)
#dev.off()
