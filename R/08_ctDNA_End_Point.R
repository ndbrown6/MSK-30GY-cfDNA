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


data_ = idx_metrics_ft %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(mean_af = 100*(mean(mean_af, na.rm = TRUE)+1E-5)) %>%
	dplyr::ungroup() %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3")) %>%
	reshape2::dcast(patient_id_mskcc ~ timepoint_weeks_since_start_of_RT, value.var = "mean_af") %>%
	dplyr::mutate(rwk1 = (wk1/`Pre-treatment`),
		      rwk2 = (wk2/`Pre-treatment`),
		      rwk3 = (wk3/`Pre-treatment`)) %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
        dplyr::select(`Pre-treatment`, wk1, wk2, wk3, rwk1, rwk2, rwk3,
		      tumor_purity = purity_wes,
		      tumor_volume = plan_volume,
		      tumor_size = primary_tumor_size_cm,
		      age = age,
		      t_stage = t_stage,
		      n_stage = n_stage,
		      smoking_status = smoking_category_yes_never,
		      hypoxia = simplified_hypoxia_group,
		      neck_dissection = neck_dissection_yes_no) %>%
	readr::type_convert() %>%
	dplyr::mutate(`Pre-treatment` = log2(`Pre-treatment`),
		      wk1 = log2(wk1),
		      wk2 = log2(wk2),
		      wk3 = log2(wk3),
		      tumor_purity = scale(tumor_purity),
		      tumor_volume = scale(tumor_volume),
		      tumor_size = scale(tumor_size),
		      age = scale(age)) %>%
        as.data.frame()

fit_ = lapply(c("rwk1", "rwk2", "rwk3"), function(x, y) {
		y %>%
		dplyr::select(x, neck_dissection, tumor_purity, tumor_volume, tumor_size, age, t_stage, n_stage, smoking_status) %>%
		dplyr::mutate(neck_dissection = case_when(
			neck_dissection == "yes" ~ "1",
			neck_dissection == "no" ~ "0"
		)) %>%
		readr::type_convert() %>%
		stats::glm(neck_dissection ~ ., family = binomial(link = 'logit'), data = .)}, y = data_)

smry_ = lapply(fit_, summary)



plot_ = fit_wk0 %>%
	dplyr::bind_rows(fit_wk1) %>%
	dplyr::bind_rows(fit_wk2) %>%
	dplyr::bind_rows(fit_wk3) %>%
	dplyr::mutate(is_significant = ifelse(`Pr(>|z|)`<.1, "Yes", "No")) %>%
	ggplot(aes(x = as.factor(variable):as.factor(week), ymin = 0, ymax = Estimate, y = Estimate, fill = is_significant, size = -log10(`Pr(>|z|)`), shape = week)) +
	geom_linerange(stat = "identity", size = .5) +
	geom_hline(yintercept = 0, size = 1) +
	geom_point(stat = "identity") +
	xlab("") +
	ylab("") +
	scale_fill_manual(values = c("#bdbdbd", "#e41a1c")) +
	scale_shape_manual(values = c("Pre-treatment" = 21, "Week 1" = 22, "Week 2" = 23, "Week 3" = 24)) +
	scale_x_discrete() +
	scale_y_continuous() +
	coord_flip() +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = FALSE,
	       shape = guide_legend(title = "Week", order = 1, override.aes = list(size = 3)),
	       size = guide_legend(title = expression(-Log[10]~"p-value")))
	
