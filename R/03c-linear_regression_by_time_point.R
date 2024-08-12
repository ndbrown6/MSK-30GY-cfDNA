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
		 dplyr::filter(fragment_length == 37) %>%
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

data_ = idx_metrics_ft %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT=="Pre-treatment") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc) %>%
	dplyr::summarize(mean_af = mean(mean_af*100),
			 aligned_reads = mean(log10(aligned_reads)),
			 concentration_ng_uL = mean(concentration_ng_uL)) %>%
        dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
        dplyr::select(aligned_reads,
		      mean_af,
		      cfdna_concentration = concentration_ng_uL,
		      hpv_copynumber = hpv_panel_copynumber,
		      tumor_purity = purity.wes,
		      tumor_volume = PlanVol,
		      tumor_size = primary_tumor_size_cm,
		      age = Age,
		      t_stage = TStage,
		      n_stage = NStage,
		      smoking_status = Smoking2,
		      hypoxia = Simplified_hypoxia_group,
		      ssgsea_cell_cycle = ssGSEA_pathway_CELL_CYCLE,
		      ssgsea_cell_death = ssGSEA_pathway_CELL_DEATH) %>%
	readr::type_convert() %>%
	dplyr::mutate(mean_af = scale(mean_af),
		      cfdna_concentration = scale(cfdna_concentration),
		      hpv_copynumber = scale(hpv_copynumber),
		      tumor_volume = scale(tumor_volume),
		      tumor_size = scale(tumor_size),
		      age = scale(age),
		      ssgsea_cell_cycle = scale(ssgsea_cell_cycle),
		      ssgsea_cell_death = scale(ssgsea_cell_death)) %>%
        tidyr::drop_na() %>%
        as.data.frame()

fit_ = lm(formula = aligned_reads ~ ., data = data_)
wk0 = summary(fit_)

data_ = idx_metrics_ft %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT=="wk1") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc) %>%
	dplyr::summarize(mean_af = mean(mean_af*100),
			 aligned_reads = mean(log10(aligned_reads)),
			 concentration_ng_uL = mean(concentration_ng_uL)) %>%
        dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
        dplyr::select(aligned_reads,
		      mean_af,
		      cfdna_concentration = concentration_ng_uL,
		      hpv_copynumber = hpv_panel_copynumber,
		      tumor_purity = purity.wes,
		      tumor_volume = PlanVol,
		      tumor_size = primary_tumor_size_cm,
		      age = Age,
		      t_stage = TStage,
		      n_stage = NStage,
		      smoking_status = Smoking2,
		      hypoxia = Simplified_hypoxia_group,
		      ssgsea_cell_cycle = ssGSEA_pathway_CELL_CYCLE,
		      ssgsea_cell_death = ssGSEA_pathway_CELL_DEATH) %>%
	readr::type_convert() %>%
	dplyr::mutate(mean_af = scale(mean_af),
		      cfdna_concentration = scale(cfdna_concentration),
		      hpv_copynumber = scale(hpv_copynumber),
		      tumor_volume = scale(tumor_volume),
		      tumor_size = scale(tumor_size),
		      age = scale(age),
		      ssgsea_cell_cycle = scale(ssgsea_cell_cycle),
		      ssgsea_cell_death = scale(ssgsea_cell_death)) %>%
        tidyr::drop_na() %>%
        as.data.frame()

fit_ = lm(formula = aligned_reads ~ ., data = data_)
wk1 = summary(fit_)

data_ = idx_metrics_ft %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT=="wk2") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc) %>%
	dplyr::summarize(mean_af = mean(mean_af*100),
			 aligned_reads = mean(log10(aligned_reads)),
			 concentration_ng_uL = mean(concentration_ng_uL)) %>%
        dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
        dplyr::select(aligned_reads,
		      mean_af,
		      cfdna_concentration = concentration_ng_uL,
		      hpv_copynumber = hpv_panel_copynumber,
		      tumor_purity = purity.wes,
		      tumor_volume = PlanVol,
		      tumor_size = primary_tumor_size_cm,
		      age = Age,
		      t_stage = TStage,
		      n_stage = NStage,
		      smoking_status = Smoking2,
		      hypoxia = Simplified_hypoxia_group,
		      ssgsea_cell_cycle = ssGSEA_pathway_CELL_CYCLE,
		      ssgsea_cell_death = ssGSEA_pathway_CELL_DEATH) %>%
	readr::type_convert() %>%
	dplyr::mutate(mean_af = scale(mean_af),
		      cfdna_concentration = scale(cfdna_concentration),
		      hpv_copynumber = scale(hpv_copynumber),
		      tumor_volume = scale(tumor_volume),
		      tumor_size = scale(tumor_size),
		      age = scale(age),
		      ssgsea_cell_cycle = scale(ssgsea_cell_cycle),
		      ssgsea_cell_death = scale(ssgsea_cell_death)) %>%
        tidyr::drop_na() %>%
        as.data.frame()

fit_ = lm(formula = aligned_reads ~ ., data = data_)
wk2 = summary(fit_)

data_ = idx_metrics_ft %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT=="wk3") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc) %>%
	dplyr::summarize(mean_af = mean(mean_af*100),
			 aligned_reads = mean(log10(aligned_reads)),
			 concentration_ng_uL = mean(concentration_ng_uL)) %>%
        dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
        dplyr::select(aligned_reads,
		      mean_af,
		      cfdna_concentration = concentration_ng_uL,
		      hpv_copynumber = hpv_panel_copynumber,
		      tumor_purity = purity.wes,
		      tumor_volume = PlanVol,
		      tumor_size = primary_tumor_size_cm,
		      age = Age,
		      t_stage = TStage,
		      n_stage = NStage,
		      smoking_status = Smoking2,
		      hypoxia = Simplified_hypoxia_group,
		      ssgsea_cell_cycle = ssGSEA_pathway_CELL_CYCLE,
		      ssgsea_cell_death = ssGSEA_pathway_CELL_DEATH) %>%
	readr::type_convert() %>%
	dplyr::mutate(mean_af = scale(mean_af),
		      cfdna_concentration = scale(cfdna_concentration),
		      hpv_copynumber = scale(hpv_copynumber),
		      tumor_volume = scale(tumor_volume),
		      tumor_size = scale(tumor_size),
		      age = scale(age),
		      ssgsea_cell_cycle = scale(ssgsea_cell_cycle),
		      ssgsea_cell_death = scale(ssgsea_cell_death)) %>%
        tidyr::drop_na() %>%
        as.data.frame()

fit_ = lm(formula = aligned_reads ~ ., data = data_)
wk3 = summary(fit_)

data_ = idx_metrics_ft %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT=="wk6") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc) %>%
	dplyr::summarize(mean_af = mean(mean_af*100),
			 aligned_reads = mean(log10(aligned_reads)),
			 concentration_ng_uL = mean(concentration_ng_uL)) %>%
        dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
        dplyr::select(aligned_reads,
		      mean_af,
		      cfdna_concentration = concentration_ng_uL,
		      hpv_copynumber = hpv_panel_copynumber,
		      tumor_purity = purity.wes,
		      tumor_volume = PlanVol,
		      tumor_size = primary_tumor_size_cm,
		      age = Age,
		      t_stage = TStage,
		      n_stage = NStage,
		      smoking_status = Smoking2,
		      hypoxia = Simplified_hypoxia_group,
		      ssgsea_cell_cycle = ssGSEA_pathway_CELL_CYCLE,
		      ssgsea_cell_death = ssGSEA_pathway_CELL_DEATH) %>%
	readr::type_convert() %>%
	dplyr::mutate(mean_af = scale(mean_af),
		      cfdna_concentration = scale(cfdna_concentration),
		      hpv_copynumber = scale(hpv_copynumber),
		      tumor_volume = scale(tumor_volume),
		      tumor_size = scale(tumor_size),
		      age = scale(age),
		      ssgsea_cell_cycle = scale(ssgsea_cell_cycle),
		      ssgsea_cell_death = scale(ssgsea_cell_death)) %>%
        tidyr::drop_na() %>%
        as.data.frame()

fit_ = lm(formula = aligned_reads ~ ., data = data_)
wk5 = summary(fit_)

p_values = do.call(cbind, lapply(list(wk0, wk1, wk2, wk3, wk5), function(x) { x$coefficients[,"Pr(>|t|)"] })) %>%
	   as.data.frame() %>%
	   tibble::rownames_to_column("variable") %>%
	   dplyr::as_tibble() %>%
	   dplyr::filter(variable != "(Intercept)") %>%
	   dplyr::rename("Pre-treatment" = V1,
			 "wk1" = V2,
			 "wk2" = V3,
			 "wk3" = V4,
			 "wk5" = V5)
	   
pdf(file = "../res/Linear_Regression_Coefficients_bywk.pdf", width = 3.15, height = 3.5)
draw(Heatmap(matrix = p_values %>%
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

plot_ = idx_metrics_ft %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(aligned_reads = mean(log10(aligned_reads))) %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3", "wk5")) %>%
	ggplot(aes(x = timepoint_weeks_since_start_of_RT, y = aligned_reads)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(1, 7),
			   breaks = c(2, 3, 4, 5, 6, 7),
			   labels = scientific_10(10^c(2, 3, 4, 5, 6, 7))) +
	xlab("") +
	ylab("cfDNA Aligned HPV Read Pairs") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Number_Read_Pairs_bywk.pdf", width = 3.75, height = 3.25)
print(plot_)
dev.off()

plot_ = idx_metrics_ft %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(mean_af = mean(mean_af)) %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3", "wk5")) %>%
	ggplot(aes(x = timepoint_weeks_since_start_of_RT, y = (mean_af*100)+(1E-3))) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	scale_x_discrete() +
	scale_y_log10(limits = c(NA, 100),
		      breaks = c(.001, .01, .1, 1, 10, 100),
		      labels = c(".001", ".01", ".1", "1", "10", "100")) +
	xlab("") +
	ylab("Mean PCM AF (%)") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Mean_AF_bywk.pdf", width = 3.75, height = 3.25)
print(plot_)
dev.off()