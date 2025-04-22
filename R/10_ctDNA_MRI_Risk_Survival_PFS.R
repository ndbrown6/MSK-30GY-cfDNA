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

mrd_smry = readr::read_tsv(file = url_mrd_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::rename(sample_id_invitae = `Invitae Biospecimen ID`) %>%
	   dplyr::left_join(manifest, by = "sample_id_invitae") %>%
	   dplyr::mutate(timepoint_weeks_since_start_of_RT = floor(timepoint_days_since_start_of_RT/7)) %>%
	   dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
			timepoint_weeks_since_start_of_RT < 0 ~ "Pre-treatment",
			TRUE ~ paste0("wk", timepoint_weeks_since_start_of_RT)
	   ))

# ctDNA fraction PCM
smry_pcm = idx_metrics_ft %>%
	   dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	   dplyr::summarize(mean_af = mean(mean_af, na.rm = TRUE)) %>%
	   reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			   fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	   dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	   readr::type_convert()

# MRD +/- PCM
smry_mrd = mrd_smry %>%
	   dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	   dplyr::summarize(landmark_mrd = any(`MRD-Landmark_Result`=="PRESENT")) %>%
	   reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			   fun.aggregate = function(x) { any(x) }, fill = NaN, value.var = "landmark_mrd") %>%
	   dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	   readr::type_convert()

# Number Reads HPV
smry_hpv = idx_metrics_ft %>%
	   dplyr::filter(chromosome == "HPV-16") %>%
	   dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	   dplyr::summarize(aligned_reads = mean(aligned_reads+1, na.rm = TRUE)) %>%
	   reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			   fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "aligned_reads") %>%
	   dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	   readr::type_convert()

# Posterior Probability +/- HPV
smry_ppr = posterior_probability %>%
	   dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	   dplyr::summarize(posterior_probability = any(`Is_ctDNA`=="+ve")) %>%
	   reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			   fun.aggregate = function(x) { any(x) }, fill = NaN, value.var = "posterior_probability") %>%
	   dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	   readr::type_convert()

# Volume MRI
smry_mri = clinical %>%
	   dplyr::select(patient_id_mskcc,
		        `Pre-treatment` = MRI_rawdata_wk0,
		        `wk1` = MRI_rawdata_wk2,
		        `wk2` = MRI_rawdata_wk3,
		        `wk3` = MRI_rawdata_wk4) %>%
	   readr::type_convert()

# ADC MRI
smry_adc = clinical %>%
	   dplyr::select(patient_id_mskcc,
		        `Pre-treatment` = ADC_Mean_wk0,
		        `wk1` = ADC_Mean_wk2,
		        `wk2` = ADC_Mean_wk3,
		        `wk3` = ADC_Mean_wk4) %>%
	   readr::type_convert()


# ctDNA fraction PCM 30Gy arm
data = clinical %>%
       dplyr::filter(crt_randomization == "30Gy") %>%
       dplyr::select(patient_id_mskcc, pfs_time, pfs_event) %>%
       dplyr::left_join(smry_pcm,
		        by = "patient_id_mskcc") %>%
       tidyr::drop_na() %>%
       dplyr::mutate(`Pre-treatment` = case_when(
	       `Pre-treatment` > median(`Pre-treatment`) ~ "high",
	       `Pre-treatment` <= median(`Pre-treatment`) ~ "low"
       )) %>%
       dplyr::mutate(wk1 = case_when(
	       wk1 > median(wk1) ~ "high",
	       wk1 <= median(wk1) ~ "low"
       )) %>%
       dplyr::mutate(wk2 = case_when(
	       wk2 > median(wk2) ~ "high",
	       wk2 <= median(wk2) ~ "low"
       )) %>%
       dplyr::mutate(wk3 = case_when(
	       wk3 > median(wk3) ~ "high",
	       wk3 <= median(wk3) ~ "low"
       )) %>%
       reshape2::melt(id.vars = c("patient_id_mskcc", "pfs_time", "pfs_event"))

pdf(file = "../res/Survival_ctDNA_30Gy.pdf", width = 4, height = 5)
for (i in c("Pre-treatment", "wk1", "wk2", "wk3")) {
	fit = survfit(Surv(pfs_time, pfs_event) ~ value, data = data %>% dplyr::filter(variable == i))
	p = ggsurvplot(fit = fit,
		       data = data %>% dplyr::filter(variable == i),
		       palette = c("#fc8d62", "#8da0cb"),
		       risk.table = TRUE,
		       pval = TRUE,
		       conf.int = FALSE,
		       xlim = c(0, 48),
		       xlab = "Time (months)",
		       ylab = "Survival rate (%)",    
		       break.time.by = 12,
		       ggtheme = theme_classic(),
		       risk.table.y.text.col = TRUE,
		       risk.table.y.text = FALSE,
		       tables.theme = theme_void())
	p$plot = p$plot +
		 ggtitle(gsub("wk", "Week ", i, fixed = TRUE)) +
		 scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1),
				    labels = c(0, 20, 40, 60, 80, 100)) +
		 theme(axis.title.x = element_text(margin = margin(t = 20)),
		       axis.title.y = element_text(margin = margin(r = 20)),
		       axis.text.x = element_text(size = 12),
		       axis.text.y = element_text(size = 12),
		       plot.title = element_text(hjust = 0.5),
		       legend.position = "bottom") +
		guides(color = guide_legend(title = " "))

	if (i=="Pre-treatment") {
		print(p, newpage = FALSE)
	} else {
		print(p, newpage = TRUE)
	}
}
dev.off()

# Landmark MRD ctDNA 30Gy arm
data = clinical %>%
       dplyr::filter(crt_randomization == "30Gy") %>%
       dplyr::select(patient_id_mskcc, pfs_time, pfs_event) %>%
       dplyr::left_join(smry_mrd,
		        by = "patient_id_mskcc") %>%
       tidyr::drop_na() %>%
       reshape2::melt(id.vars = c("patient_id_mskcc", "pfs_time", "pfs_event"))

pdf(file = "../res/Survival_Landmark_MRD_30Gy.pdf", width = 4, height = 5)
for (i in c("Pre-treatment", "wk1", "wk2", "wk3")) {
	fit = survfit(Surv(pfs_time, pfs_event) ~ value, data = data %>% dplyr::filter(variable == i))
	p = ggsurvplot(fit = fit,
		       data = data %>% dplyr::filter(variable == i),
		       palette = c("#8da0cb", "#fc8d62"),
		       risk.table = TRUE,
		       pval = TRUE,
		       conf.int = FALSE,
		       xlim = c(0, 48),
		       xlab = "Time (months)",
		       ylab = "Survival rate (%)",    
		       break.time.by = 12,
		       ggtheme = theme_classic(),
		       risk.table.y.text.col = TRUE,
		       risk.table.y.text = FALSE,
		       tables.theme = theme_void())
	p$plot = p$plot +
		 ggtitle(gsub("wk", "Week ", i, fixed = TRUE)) +
		 scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1),
				    labels = c(0, 20, 40, 60, 80, 100)) +
		 theme(axis.title.x = element_text(margin = margin(t = 20)),
		       axis.title.y = element_text(margin = margin(r = 20)),
		       axis.text.x = element_text(size = 12),
		       axis.text.y = element_text(size = 12),
		       plot.title = element_text(hjust = 0.5),
		       legend.position = "bottom") +
		guides(color = guide_legend(title = " "))

	if (i=="Pre-treatment") {
		print(p, newpage = FALSE)
	} else {
		print(p, newpage = TRUE)
	}
}
dev.off()

# Number HPV Read pairs 30Gy arm
data = clinical %>%
       dplyr::filter(crt_randomization == "30Gy") %>%
       dplyr::select(patient_id_mskcc, pfs_time, pfs_event) %>%
       dplyr::left_join(smry_hpv,
		        by = "patient_id_mskcc") %>%
       tidyr::drop_na() %>%
       dplyr::mutate(`Pre-treatment` = case_when(
	       `Pre-treatment` > median(`Pre-treatment`) ~ "high",
	       `Pre-treatment` <= median(`Pre-treatment`) ~ "low"
       )) %>%
       dplyr::mutate(wk1 = case_when(
	       wk1 > median(wk1) ~ "high",
	       wk1 <= median(wk1) ~ "low"
       )) %>%
       dplyr::mutate(wk2 = case_when(
	       wk2 > median(wk2) ~ "high",
	       wk2 <= median(wk2) ~ "low"
       )) %>%
       dplyr::mutate(wk3 = case_when(
	       wk3 > median(wk3) ~ "high",
	       wk3 <= median(wk3) ~ "low"
       )) %>%
       reshape2::melt(id.vars = c("patient_id_mskcc", "pfs_time", "pfs_event"))

pdf(file = "../res/Survival_HPV_30Gy.pdf", width = 4, height = 5)
for (i in c("Pre-treatment", "wk1", "wk2", "wk3")) {
	fit = survfit(Surv(pfs_time, pfs_event) ~ value, data = data %>% dplyr::filter(variable == i))
	p = ggsurvplot(fit = fit,
		       data = data %>% dplyr::filter(variable == i),
		       palette = c("#fc8d62", "#8da0cb"),
		       risk.table = TRUE,
		       pval = TRUE,
		       conf.int = FALSE,
		       xlim = c(0, 48),
		       xlab = "Time (months)",
		       ylab = "Survival rate (%)",    
		       break.time.by = 12,
		       ggtheme = theme_classic(),
		       risk.table.y.text.col = TRUE,
		       risk.table.y.text = FALSE,
		       tables.theme = theme_void())
	p$plot = p$plot +
		 ggtitle(gsub("wk", "Week ", i, fixed = TRUE)) +
		 scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1),
				    labels = c(0, 20, 40, 60, 80, 100)) +
		 theme(axis.title.x = element_text(margin = margin(t = 20)),
		       axis.title.y = element_text(margin = margin(r = 20)),
		       axis.text.x = element_text(size = 12),
		       axis.text.y = element_text(size = 12),
		       plot.title = element_text(hjust = 0.5),
		       legend.position = "bottom") +
		guides(color = guide_legend(title = " "))

	if (i=="Pre-treatment") {
		print(p, newpage = FALSE)
	} else {
		print(p, newpage = TRUE)
	}
}
dev.off()

# Posterior Probability HPV Read pairs 30Gy arm
data = clinical %>%
       dplyr::filter(crt_randomization == "30Gy") %>%
       dplyr::select(patient_id_mskcc, pfs_time, pfs_event) %>%
       dplyr::left_join(smry_ppr,
		        by = "patient_id_mskcc") %>%
       tidyr::drop_na() %>%
       reshape2::melt(id.vars = c("patient_id_mskcc", "pfs_time", "pfs_event"))

pdf(file = "../res/Survival_Posterior_Probability_HPV_30Gy.pdf", width = 4, height = 5)
for (i in c("Pre-treatment", "wk1", "wk2", "wk3")) {
	fit = survfit(Surv(pfs_time, pfs_event) ~ value, data = data %>% dplyr::filter(variable == i))
	p = ggsurvplot(fit = fit,
		       data = data %>% dplyr::filter(variable == i),
		       palette = c("#8da0cb", "#fc8d62"),
		       risk.table = TRUE,
		       pval = TRUE,
		       conf.int = FALSE,
		       xlim = c(0, 48),
		       xlab = "Time (months)",
		       ylab = "Survival rate (%)",    
		       break.time.by = 12,
		       ggtheme = theme_classic(),
		       risk.table.y.text.col = TRUE,
		       risk.table.y.text = FALSE,
		       tables.theme = theme_void())
	p$plot = p$plot +
		 ggtitle(gsub("wk", "Week ", i, fixed = TRUE)) +
		 scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1),
				    labels = c(0, 20, 40, 60, 80, 100)) +
		 theme(axis.title.x = element_text(margin = margin(t = 20)),
		       axis.title.y = element_text(margin = margin(r = 20)),
		       axis.text.x = element_text(size = 12),
		       axis.text.y = element_text(size = 12),
		       plot.title = element_text(hjust = 0.5),
		       legend.position = "bottom") +
		guides(color = guide_legend(title = " "))

	if (i=="Pre-treatment") {
		print(p, newpage = FALSE)
	} else {
		print(p, newpage = TRUE)
	}
}
dev.off()

# MRI Volume 30Gy arm
data = clinical %>%
       dplyr::filter(crt_randomization == "30Gy") %>%
       dplyr::select(patient_id_mskcc, pfs_time, pfs_event) %>%
       dplyr::left_join(smry_mri,
		        by = "patient_id_mskcc") %>%
       tidyr::drop_na() %>%
       dplyr::mutate(`Pre-treatment` = case_when(
	       `Pre-treatment` > median(`Pre-treatment`) ~ "high",
	       `Pre-treatment` <= median(`Pre-treatment`) ~ "low"
       )) %>%
       dplyr::mutate(wk1 = case_when(
	       wk1 > median(wk1) ~ "high",
	       wk1 <= median(wk1) ~ "low"
       )) %>%
       dplyr::mutate(wk2 = case_when(
	       wk2 > median(wk2) ~ "high",
	       wk2 <= median(wk2) ~ "low"
       )) %>%
       dplyr::mutate(wk3 = case_when(
	       wk3 > median(wk3) ~ "high",
	       wk3 <= median(wk3) ~ "low"
       )) %>%
       reshape2::melt(id.vars = c("patient_id_mskcc", "pfs_time", "pfs_event"))

pdf(file = "../res/Survival_Volume_30Gy.pdf", width = 4, height = 5)
for (i in c("Pre-treatment", "wk1", "wk2", "wk3")) {
	fit = survfit(Surv(pfs_time, pfs_event) ~ value, data = data %>% dplyr::filter(variable == i))
	p = ggsurvplot(fit = fit,
		       data = data %>% dplyr::filter(variable == i),
		       palette = c("#fc8d62", "#8da0cb"),
		       risk.table = TRUE,
		       pval = TRUE,
		       conf.int = FALSE,
		       xlim = c(0, 48),
		       xlab = "Time (months)",
		       ylab = "Survival rate (%)",    
		       break.time.by = 12,
		       ggtheme = theme_classic(),
		       risk.table.y.text.col = TRUE,
		       risk.table.y.text = FALSE,
		       tables.theme = theme_void())
	p$plot = p$plot +
		 ggtitle(gsub("wk", "Week ", i, fixed = TRUE)) +
		 scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1),
				    labels = c(0, 20, 40, 60, 80, 100)) +
		 theme(axis.title.x = element_text(margin = margin(t = 20)),
		       axis.title.y = element_text(margin = margin(r = 20)),
		       axis.text.x = element_text(size = 12),
		       axis.text.y = element_text(size = 12),
		       plot.title = element_text(hjust = 0.5),
		       legend.position = "bottom") +
		guides(color = guide_legend(title = " "))

	if (i=="Pre-treatment") {
		print(p, newpage = FALSE)
	} else {
		print(p, newpage = TRUE)
	}
}
dev.off()

# ctDNA fraction PCM & MRI Volume 30Gy arm
data = clinical %>%
       dplyr::filter(crt_randomization == "30Gy") %>%
       dplyr::select(patient_id_mskcc, pfs_time, pfs_event, composite_end_point) %>%
       dplyr::left_join(smry_pcm,
			by = "patient_id_mskcc") %>%
       dplyr::left_join(smry_mri,
			by = c("patient_id_mskcc")) %>%
       dplyr::mutate(`Pre-treatment.x` = case_when(
	       `Pre-treatment.x` > median(`Pre-treatment.x`, na.rm=TRUE) ~ "high",
	       `Pre-treatment.x` <= median(`Pre-treatment.x`, na.rm=TRUE) ~ "low"
       )) %>%
       dplyr::mutate(wk1.x = case_when(
	       wk1.x > median(wk1.x, na.rm=TRUE) ~ "high",
	       wk1.x <= median(wk1.x, na.rm=TRUE) ~ "low"
       )) %>%
       dplyr::mutate(wk2.x = case_when(
	       wk2.x > median(wk2.x, na.rm=TRUE) ~ "high",
	       wk2.x <= median(wk2.x, na.rm=TRUE) ~ "low"
       )) %>%
       dplyr::mutate(wk3.x = case_when(
	       wk3.x > median(wk3.x, na.rm=TRUE) ~ "high",
	       wk3.x <= median(wk3.x, na.rm=TRUE) ~ "low"
       )) %>%
       dplyr::mutate(`Pre-treatment.y` = case_when(
	       `Pre-treatment.y` > median(`Pre-treatment.y`, na.rm=TRUE) ~ "high",
	       `Pre-treatment.y` <= median(`Pre-treatment.y`, na.rm=TRUE) ~ "low"
       )) %>%
       dplyr::mutate(wk1.y = case_when(
	       wk1.y > median(wk1.y, na.rm=TRUE) ~ "high",
	       wk1.y <= median(wk1.y, na.rm=TRUE) ~ "low"
       )) %>%
       dplyr::mutate(wk2.y = case_when(
	       wk2.y > median(wk2.y, na.rm=TRUE) ~ "high",
	       wk2.y <= median(wk2.y, na.rm=TRUE) ~ "low"
       )) %>%
       dplyr::mutate(wk3.y = case_when(
	       wk3.y > median(wk3.y, na.rm=TRUE) ~ "high",
	       wk3.y <= median(wk3.y, na.rm=TRUE) ~ "low"
       )) %>%
       reshape2::melt(id.vars = c("patient_id_mskcc", "pfs_time", "pfs_event", "composite_end_point")) %>%
       dplyr::mutate(assay = case_when(
	       grepl(".x", variable, fixed = TRUE) ~ "ctDNA",
	       grepl(".y", variable, fixed = TRUE) ~ "Vol"
       )) %>%
       dplyr::mutate(variable = gsub(pattern = ".x", "", variable, fixed = TRUE)) %>%
       dplyr::mutate(variable = gsub(pattern = ".y", "", variable, fixed = TRUE)) %>%
       tidyr::drop_na()

ii = 1
pdf(file = "../res/Survival_ctDNA_Volume_30Gy.pdf", width = 4, height = 5.5)
for (i in c("Pre-treatment", "wk1", "wk2", "wk3")) {
	for (j in c("Pre-treatment", "wk1", "wk2", "wk3")) {
		tmp = data %>%
		      dplyr::filter((variable == i & assay == "ctDNA") |
				    (variable == j & assay == "Vol")) %>%
		      reshape2::dcast(patient_id_mskcc + pfs_time + pfs_event ~ assay, value.var = "value") %>%
		      tidyr::drop_na() %>%
		      dplyr::mutate(value = case_when(
			      Vol == "high" & ctDNA == "high" ~ "high",
			      TRUE ~ "low"
		      ))
		fit = survfit(Surv(pfs_time, pfs_event) ~ value, data = tmp)
		p = ggsurvplot(fit = fit,
			       data = tmp,
			       palette = c("#fc8d62", "#8da0cb"),
			       risk.table = TRUE,
			       pval = TRUE,
			       conf.int = FALSE,
			       xlim = c(0, 48),
			       xlab = "Time (months)",
			       ylab = "Survival rate (%)",    
			       break.time.by = 12,
			       ggtheme = theme_classic(),
			       risk.table.y.text.col = TRUE,
			       risk.table.y.text = FALSE,
			       tables.theme = theme_void())
		p$plot = p$plot +
			 ggtitle(paste0(gsub("wk", "Week ", i, fixed = T), " ctDNA\n", gsub("wk", "Week ", j, fixed = T), " Vol")) +
			 scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1),
					    labels = c(0, 20, 40, 60, 80, 100)) +
			 theme(axis.title.x = element_text(margin = margin(t = 20)),
			       axis.title.y = element_text(margin = margin(r = 20)),
			       axis.text.x = element_text(size = 12),
			       axis.text.y = element_text(size = 12),
			       plot.title = element_text(hjust = 0.5, size = 10),
			       legend.position = "bottom") +
			 guides(color = guide_legend(title = " "))

		if (ii==1) {
			print(p, newpage = FALSE)
		} else {
			print(p, newpage = TRUE)
		}
		ii = ii + 1
	}
}
dev.off()

ii = 1
pdf(file = "../res/Survival_ctDNA_Volume_Risk_group_30Gy.pdf", width = 4, height = 5.5)
for (i in c("Pre-treatment", "wk1", "wk2", "wk3")) {
	for (j in c("Pre-treatment", "wk1", "wk2", "wk3")) {
		tmp = data %>%
		      dplyr::filter((variable == i & assay == "ctDNA") |
				    (variable == j & assay == "Vol")) %>%
		      reshape2::dcast(patient_id_mskcc + pfs_time + pfs_event + composite_end_point ~ assay, value.var = "value") %>%
		      tidyr::drop_na() %>%
		      dplyr::mutate(value = case_when(
			      Vol == "high" & ctDNA == "high" & composite_end_point ~ "high",
			      TRUE ~ "low"
		      ))
		fit = survfit(Surv(pfs_time, pfs_event) ~ value, data = tmp)
		p = ggsurvplot(fit = fit,
			       data = tmp,
			       palette = c("#fc8d62", "#8da0cb"),
			       risk.table = TRUE,
			       pval = TRUE,
			       conf.int = FALSE,
			       xlim = c(0, 48),
			       xlab = "Time (months)",
			       ylab = "Survival rate (%)",    
			       break.time.by = 12,
			       ggtheme = theme_classic(),
			       risk.table.y.text.col = TRUE,
			       risk.table.y.text = FALSE,
			       tables.theme = theme_void())
		p$plot = p$plot +
			 ggtitle(paste0(gsub("wk", "Week ", i, fixed = T), " ctDNA\n", gsub("wk", "Week ", j, fixed = T), " Vol")) +
			 scale_y_continuous(breaks = c(0, .2, .4, .6, .8, 1),
					    labels = c(0, 20, 40, 60, 80, 100)) +
			 theme(axis.title.x = element_text(margin = margin(t = 20)),
			       axis.title.y = element_text(margin = margin(r = 20)),
			       axis.text.x = element_text(size = 12),
			       axis.text.y = element_text(size = 12),
			       plot.title = element_text(hjust = 0.5, size = 10),
			       legend.position = "bottom") +
			 guides(color = guide_legend(title = " "))

		if (ii==1) {
			print(p, newpage = FALSE)
		} else {
			print(p, newpage = TRUE)
		}
		ii = ii + 1
	}
}
dev.off()
