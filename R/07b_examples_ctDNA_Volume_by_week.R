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
# Mean AF by Week
#==================================================
p_mrd = idx_metrics_ft %>%
	dplyr::left_join(mrd_smry %>%
			 dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)),
			 by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = case_when(
		`MRD-Monitoring_Result` == "PRESENT" ~ "+ve",
		`MRD-Monitoring_Result` == "ABSENT" ~ "-ve"
	)) %>%
	dplyr::filter(patient_id_mskcc == "CTMS-161") %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3")) %>%
	dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
		timepoint_weeks_since_start_of_RT == "Pre-treatment" ~ "0",
		timepoint_weeks_since_start_of_RT == "wk1" ~ "1",
		timepoint_weeks_since_start_of_RT == "wk2" ~ "2",
		timepoint_weeks_since_start_of_RT == "wk3" ~ "3"
	)) %>%
	readr::type_convert() %>%
	dplyr::select(patient_id_mskcc, timepoint_weeks_since_start_of_RT, measure = mean_af, Is_ctDNA)
p_mrd = p_mrd %>%
	dplyr::mutate(measure = measure / (p_mrd %>%
					   dplyr::filter(timepoint_weeks_since_start_of_RT==0) %>%
					   .[["measure"]])) %>%
	dplyr::mutate(assay = "PCM")

#==================================================
# Number of Read Pairs by Week
#==================================================
p_hpv = idx_metrics_ft %>%
	dplyr::left_join(posterior_probability %>%
			 dplyr::select(sample_name, Is_ctDNA), by = "sample_name") %>%
	dplyr::filter(patient_id_mskcc == "CTMS-161") %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3")) %>%
	dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
		timepoint_weeks_since_start_of_RT == "Pre-treatment" ~ "0",
		timepoint_weeks_since_start_of_RT == "wk1" ~ "1",
		timepoint_weeks_since_start_of_RT == "wk2" ~ "2",
		timepoint_weeks_since_start_of_RT == "wk3" ~ "3"
	)) %>%
	readr::type_convert() %>%
	dplyr::select(patient_id_mskcc, timepoint_weeks_since_start_of_RT, measure = aligned_reads, Is_ctDNA)
p_hpv = p_hpv %>%
	dplyr::mutate(measure = measure / (p_hpv %>%
					   dplyr::filter(timepoint_weeks_since_start_of_RT==0) %>%
					   .[["measure"]])) %>%
	dplyr::mutate(assay = "HPV")

#==================================================
# MRI Volume by Week
#==================================================
p_mri = clinical %>%
	dplyr::select(patient_id_mskcc,
		      `Pre-treatment` = MRI_rawdata_wk0,
		      `wk1` = MRI_rawdata_wk1,
		      `wk2` = MRI_rawdata_wk2,
		      `wk3` = MRI_rawdata_wk3) %>%
	reshape2::melt(variable.name = "timepoint_weeks_since_start_of_RT", value.name = "MRI_volume") %>%
	dplyr::filter(patient_id_mskcc == "CTMS-161") %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3")) %>%
	dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
		timepoint_weeks_since_start_of_RT == "Pre-treatment" ~ "0",
		timepoint_weeks_since_start_of_RT == "wk1" ~ "1",
		timepoint_weeks_since_start_of_RT == "wk2" ~ "2",
		timepoint_weeks_since_start_of_RT == "wk3" ~ "3"
	)) %>%
	readr::type_convert() %>%
	dplyr::select(patient_id_mskcc, timepoint_weeks_since_start_of_RT, measure = MRI_volume)
p_mri = p_mri %>%
	dplyr::mutate(measure = measure / (p_mri %>%
					   dplyr::filter(timepoint_weeks_since_start_of_RT==0) %>%
					   .[["measure"]])) %>%
	dplyr::mutate(Is_ctDNA = "Not Applicable", assay = "MRI")

plot_ = p_mrd %>%
	dplyr::bind_rows(p_hpv) %>%
	dplyr::bind_rows(p_mri) %>%
	dplyr::mutate(measure = measure * 100) %>%
	dplyr::mutate(assay = factor(assay, levels = c("PCM", "HPV", "MRI"), ordered = TRUE)) %>%
	ggplot(aes(x = timepoint_weeks_since_start_of_RT, y = measure, shape = Is_ctDNA, color = assay)) +
	geom_line(mapping = aes(x = timepoint_weeks_since_start_of_RT, y = measure, color = assay),
		  stat = "identity", size = 1.25, alpha = .75, inherit.aes = FALSE) +
	geom_point(stat = "identity", size = 4, fill = "white") +
	scale_shape_manual(values = c(21, 25, 22)) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_x_continuous(breaks = 0:3,
			   labels = c("Pre-treat\n-ment", "Week 1", "Week 2", "Week 3")) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 25, 50, 75, 100),
			   labels = c(100, 75, 50, 25, 0)) +
	xlab("") +
	ylab("Relative Change from Baseline (%)") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(margin = margin(t = 7), size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = guide_legend(title = "Assay Modality", order = 1),
	       shape = guide_legend(title = "ctDNA Status", order = 2))

pdf(file = "../res/CTMS-30_MRI_HPV_MRD_Combined.pdf", width = 6, height = 3)
print(plot_)
dev.off()

#==================================================
# % Reduction by week
#==================================================
smry_af = idx_metrics_ft %>%
	  dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	  dplyr::summarize(mean_af = mean(mean_af+1E-5, na.rm = TRUE)) %>%
	  reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			  fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	  dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	  dplyr::mutate(wk1 = wk1/`Pre-treatment`,
		        wk2 = wk2/`Pre-treatment`,
		        wk3 = wk3/`Pre-treatment`,
		        `Pre-treatment` = 0) %>%
	  readr::type_convert()

smry_reads = idx_metrics_ft %>%
	     dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	     dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	     dplyr::summarize(aligned_reads = mean(aligned_reads+1, na.rm = TRUE)) %>%
	     reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			     fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "aligned_reads") %>%
	     dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	     dplyr::mutate(wk1 = wk1/`Pre-treatment`,
			   wk2 = wk2/`Pre-treatment`,
			   wk3 = wk3/`Pre-treatment`,
			   `Pre-treatment` = 0) %>%
	     readr::type_convert()

smry_mri = clinical %>%
	   dplyr::select(patient_id_mskcc,
		        `Pre-treatment` = MRI_rawdata_wk0,
		        `wk1` = MRI_rawdata_wk2,
		        `wk2` = MRI_rawdata_wk3,
		        `wk3` = MRI_rawdata_wk4) %>%
	   dplyr::mutate(wk1 = wk1/`Pre-treatment`,
			 wk2 = wk2/`Pre-treatment`,
			 wk3 = wk3/`Pre-treatment`,
			 `Pre-treatment` = 0) %>%
	   readr::type_convert()

plot_ = smry_af %>%
	reshape2::melt(value.name = "PCM", variable.name = "week") %>%
	dplyr::full_join(smry_reads %>%
			 reshape2::melt(value.name = "HPV", variable.name = "week"), by = c("patient_id_mskcc", "week")) %>%
	dplyr::full_join(smry_mri %>%
			 reshape2::melt(value.name = "MRI", variable.name = "week"), by = c("patient_id_mskcc", "week")) %>%
	readr::type_convert() %>%
	dplyr::filter(week != "Pre-treatment") %>%
	dplyr::group_by(week) %>%
	dplyr::summarize(`30_%_pcm` = 100*sum(PCM<=.7, na.rm=TRUE)/sum(!is.na(PCM)),
			 `30_%_hpv` = 100*sum(HPV<=.7, na.rm=TRUE)/sum(!is.na(HPV)),
			 `30_%_mri` = 100*sum(MRI<=.7, na.rm=TRUE)/sum(!is.na(MRI)),
			 `90_%_pcm` = 100*sum(PCM<=.1, na.rm=TRUE)/sum(!is.na(PCM)),
			 `90_%_hpv` = 100*sum(HPV<=.1, na.rm=TRUE)/sum(!is.na(HPV)),
			 `90_%_mri` = 100*sum(MRI<=.1, na.rm=TRUE)/sum(!is.na(MRI))) %>%
	reshape2::melt(value.name = "%", variable.name = "assay") %>%
	dplyr::mutate(`%_reduction` = case_when(
		grepl("30_%", assay, fixed = TRUE) ~ "30%",
		grepl("90_%", assay, fixed = TRUE) ~ "90%"
	)) %>%
	dplyr::mutate(assay = toupper(gsub("30_%_|90_%_", "", assay, perl = TRUE))) %>%
	dplyr::mutate(week = gsub(pattern = "wk", replacement = "", x = week, fixed = TRUE)) %>%
	readr::type_convert() %>%
	dplyr::mutate(assay = factor(assay, levels = c("PCM", "HPV", "MRI"), ordered = TRUE)) %>%
	dplyr::mutate(`%_reduction` = factor(`%_reduction`, levels = c("30%", "90%"), ordered = TRUE)) %>%
	ggplot(aes(x = week, y = `%`, color = assay, shape = assay, linetype = `%_reduction`)) +
	geom_line(stat = "identity", size = .75, alpha = .95) +
	geom_point(stat = "identity", fill = "white", size = 3) +
	scale_shape_manual(values = c(21, 25, 22)) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_linetype_manual(values = c(1, 3)) +
	scale_x_continuous(expand = c(.1,.1),
			   breaks = 1:3,
			   labels = c("1", "2", "3")) +
	scale_y_continuous(limits = c(0, 100),
			   breaks = c(0, 25, 50, 75, 100),
			   labels = c(0, 25, 50, 75, 100)) +
	xlab("Weeks") +
	ylab("Fraction of Patients (%)") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(margin = margin(t = 7), size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	guides(color = guide_legend(title = "Assay Modality"),
	       shape = guide_legend(title = "Assay Modality"),
	       linetype = guide_legend(title = "Baseline\nReduction"))

pdf(file = "../res/Fraction_Patients_30_90_Percentage_Reduction.pdf", width = 4, height = 3)
print(plot_)
dev.off()

#==================================================
# Empirical CDF
#==================================================
ecdf_mri = smry_mri %>%
	   dplyr::select(patient_id_mskcc, wk2) %>%
	   tidyr::drop_na() %>%
	   dplyr::mutate(relative_change = log10(wk2)) %>%
	   dplyr::arrange(relative_change) %>%
	   .[["relative_change"]] %>%
	   ecdf()

ecdf_pcm = smry_af %>%
	   dplyr::select(patient_id_mskcc, wk2) %>%
	   tidyr::drop_na() %>%
	   dplyr::mutate(relative_change = log10(wk2)) %>%
	   dplyr::arrange(relative_change) %>%
	   .[["relative_change"]] %>%
	   ecdf()

plot_ = dplyr::tibble(
		relative_change = knots(ecdf_mri)
	) %>%
	dplyr::mutate(`%` = ecdf_mri(relative_change),
		      assay = "MRI") %>%
	dplyr::bind_rows(
		dplyr::tibble(
			relative_change = knots(ecdf_pcm)
		) %>%
		dplyr::mutate(`%` = ecdf_pcm(relative_change),
			      assay = "PCM")
	) %>%
	ggplot(aes(x = `%`, y = relative_change, color = assay)) +
	geom_segment(aes(x = `%`, xend = `%`, y = 0, yend = relative_change)) +
	geom_hline(yintercept = 0, color = "#bdbdbd", alpha = 1, size = .75) +
	geom_vline(aes(xintercept = xintercept),
		   data = dplyr::tibble(xintercept = c(.03, .03+.68),
					assay = c("MRI", "MRI")),
		   linetype = 3) +
	geom_vline(aes(xintercept = xintercept),
		   data = dplyr::tibble(xintercept = c(.33, .33+.19, .33+.19+.19),
					assay = c("PCM", "PCM", "PCM")),
		   linetype = 3) +
	geom_point(stat = "identity", size = 1, shape = 21, fill = "white") +
	scale_color_manual(values = c("PCM" = "#1b9e77", "MRI" = "#7570b3")) +
	scale_x_reverse() +
	scale_y_continuous(limits = c(NA, 1)) +
	xlab(" ") +
	ylab(expression("Relative Change "(Log[10]))) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      axis.line.x = element_blank(),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank(),
	      strip.text = element_text(color = "white")) +
	facet_wrap(~assay, nrow = 2, ncol = 1, scales = "free_y") +
	guides(color = guide_legend(title = "Assay Modality"))
	
	
pdf(file = "../res/Empirical_Cummulative__Relative_Change_MRI_AF.pdf", width = 6, height = 4)
print(plot_)
dev.off()

smry_af = idx_metrics_ft %>%
	  dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	  dplyr::summarize(mean_af = mean(mean_af+1E-5, na.rm = TRUE)) %>%
	  reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			  fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	  dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	  readr::type_convert()

smry_reads = idx_metrics_ft %>%
	     dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
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

#==================================================
# Absolute MRI volume wk1, wk2, wk3, wk5
#==================================================
plot_ = smry_af %>%
	dplyr::mutate(wk3 = wk3/`Pre-treatment`) %>%
	dplyr::mutate(q4 = case_when(
		wk3 >= quantile(wk3, .75, na.rm=TRUE) ~ "4th",
		wk3 >= quantile(wk3, .5, na.rm=TRUE) & wk3 < quantile(wk3, .75, na.rm=TRUE) ~ "3rd",
		wk3 >= quantile(wk3, .25, na.rm=TRUE) & wk3 < quantile(wk3, .5, na.rm=TRUE) ~ "2nd",
		wk3 < quantile(wk3, .25, na.rm=TRUE) ~ "1st",
		TRUE ~"NA"
	)) %>%
	dplyr::select(patient_id_mskcc, q4) %>%
	dplyr::full_join(smry_mri, by = "patient_id_mskcc") %>%	
	reshape2::melt(id.vars = c("patient_id_mskcc", "q4"), variable.name = "week", value.name = "MRI") %>%
	readr::type_convert() %>%
	tidyr::drop_na() %>%
	ggplot(aes(x = week, y = MRI, shape = q4, color = q4)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "black") +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", size = 3, alpha = .95) +
	scale_shape_manual(values = c("1st" = 21, "2nd" = 22, "3rd" = 23, "4th" = 24)) +
	scale_color_manual(values = c("1st" = "#1b9e77", "2nd" = "#d95f02", "3rd" = "#7570b3", "4th" = "#e7298a")) +
	scale_x_discrete(breaks = c("Pre-treatment", "wk1", "wk2", "wk3"),
			 labels = c("0", "1", "2", "3")) +
	scale_y_continuous(limits = c(0, 5E4),
			   labels = scientific_10) +

	xlab("Weeks") +
	ylab(expression("MRI Volume "(cm^3))) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(margin = margin(t = 7), size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	guides(color = FALSE, shape = FALSE) +
	facet_wrap(~q4, nrow = 1)

pdf(file = "../res/ctDNA_Clearance_MRI_Volume_All_weeks.pdf", width = 2.85*1.75, height = 2.5*1.5)
print(plot_)
dev.off()
