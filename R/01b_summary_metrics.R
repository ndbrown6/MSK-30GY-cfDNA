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

mrd_smry = readr::read_tsv(file = url_mrd_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::mutate(sample_uuid = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer))

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

mutation_smry = readr::read_tsv(file = url_mutation_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		readr::type_convert()

nodal_dissection_smry = readr::read_tsv(file = url_no_node_dissection, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::rename(patient_id_mskcc = sample)

#==================================================
## Number of unique patients with MRD assay
#==================================================
smry_ = manifest %>%
	dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	dplyr::mutate(patient_id_mskcc = case_when(
			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
			TRUE ~ patient_id_mskcc
	)) %>%
	dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
	dplyr::mutate(hpv_type_wes_wgs = case_when(
		   is.na(hpv_type_wes_wgs) ~ "Unknown",
		   TRUE ~ hpv_type_wes_wgs
	)) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::group_by(patient_id_mskcc) %>%
	dplyr::summarize(n = n()) %>%
	dplyr::ungroup() %>%
	nrow() %>%
	pander::pander()

#==================================================
## Number of samples per patient
#==================================================
smry_ = manifest %>%
	dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	dplyr::mutate(patient_id_mskcc = case_when(
			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
			TRUE ~ patient_id_mskcc
	)) %>%
	dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
	dplyr::mutate(hpv_type_wes_wgs = case_when(
		   is.na(hpv_type_wes_wgs) ~ "Unknown",
		   TRUE ~ hpv_type_wes_wgs
	)) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::group_by(patient_id_mskcc) %>%
	dplyr::summarize(n = n()) %>%
	dplyr::ungroup() %>%
	dplyr::summarize(mean_n_samples = mean(n),
			 min_n_samples = min(n),
			 max_n_samples = max(n)) %>%
	pander::pander()

plot_ = manifest %>%
	dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	dplyr::mutate(patient_id_mskcc = case_when(
			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
			TRUE ~ patient_id_mskcc
	)) %>%
	dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
	dplyr::mutate(hpv_type_wes_wgs = case_when(
		   is.na(hpv_type_wes_wgs) ~ "Unknown",
		   TRUE ~ hpv_type_wes_wgs
	)) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::group_by(patient_id_mskcc) %>%
	dplyr::summarize(n = n()) %>%
	dplyr::ungroup() %>%
	ggplot(aes(x = n)) +
	geom_bar(stat = "count", fill = "#e41a1c", color = NA, alpha = .55) +
	xlab("Number of Samples") +
	ylab("Number of Patients") +
	scale_x_continuous(breaks = 1:12,
			   labels = 1:12) +
	scale_y_continuous(limits = c(0, 30)) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)))

pdf(file = "../res/Number_Patients_by_Samples.pdf", width = 5, height = 5)
print(plot_)
dev.off()

#==================================================
## Number of patient per HPV subtype
#==================================================
smry_ = manifest %>%
	dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	dplyr::mutate(patient_id_mskcc = case_when(
			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
			TRUE ~ patient_id_mskcc
	)) %>%
	dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
	dplyr::mutate(hpv_type_wes_wgs = case_when(
		   is.na(hpv_type_wes_wgs) ~ "Unknown",
		   TRUE ~ hpv_type_wes_wgs
	)) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::group_by(patient_id_mskcc, hpv_type_wes_wgs) %>%
	dplyr::summarize() %>%
	dplyr::ungroup() %>%
	dplyr::group_by(hpv_type_wes_wgs) %>%
	dplyr::summarize(n_patients = n()) %>%
	dplyr::mutate(`%_patients` = 100*n_patients/sum(n_patients)) %>%
	pander::pander()

plot_ = manifest %>%
	dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	dplyr::mutate(patient_id_mskcc = case_when(
			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
			TRUE ~ patient_id_mskcc
	)) %>%
	dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
	dplyr::mutate(hpv_type_wes_wgs = case_when(
		   is.na(hpv_type_wes_wgs) ~ "Unknown",
		   TRUE ~ hpv_type_wes_wgs
	)) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::filter(!duplicated(patient_id_mskcc)) %>%
	dplyr::mutate(hpv_type_wes_wgs = factor(hpv_type_wes_wgs, levels = sort(unique(hpv_type_wes_wgs)), ordered = TRUE)) %>%
	ggplot(aes(x = hpv_type_wes_wgs)) +
	geom_bar(stat = "count", fill = "#e41a1c", color = NA, alpha = .55) +
	xlab("Tumor HPV Type") +
	ylab("Number of Patients") +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 90)) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)))

pdf(file = "../res/Number_Patients_by_HPV_Type.pdf", width = 5, height = 5)
print(plot_)
dev.off()

#==================================================
## Time since start of RT
#==================================================
plot_ = manifest %>%
	dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	dplyr::mutate(patient_id_mskcc = case_when(
			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
			TRUE ~ patient_id_mskcc
	)) %>%
	dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
	dplyr::mutate(hpv_type_wes_wgs = case_when(
		   is.na(hpv_type_wes_wgs) ~ "Unknown",
		   TRUE ~ hpv_type_wes_wgs
	)) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::filter(!is.na(timepoint_days_since_start_of_RT)) %>%
	dplyr::mutate(hpv_type_wes_wgs = case_when(
		hpv_type_wes_wgs != "HPV-16" & hpv_type_wes_wgs != "Unknown" ~ "Other HPV",
		TRUE ~ hpv_type_wes_wgs
	)) %>%
	dplyr::mutate(hpv_type_wes_wgs = factor(hpv_type_wes_wgs, levels = sort(unique(hpv_type_wes_wgs)), ordered = TRUE)) %>%
	dplyr::arrange(hpv_type_wes_wgs, timepoint_days_since_start_of_RT) %>%
	dplyr::mutate(index = 1:nrow(.)) %>%
	ggplot(aes(x = index, y = timepoint_days_since_start_of_RT)) +
	geom_point(stat = "identity", shape = 21, color = "black", fill = "white", size = 2, alpha = .55) +
	geom_hline(yintercept = 0, color = "goldenrod3", size = .5, alpha = 1, linetype = 2) +
	geom_hline(yintercept = 730, color = "goldenrod3", size = .5, alpha = 1, linetype = 2) +
	xlab("Plasma Sample") +
	ylab("Time (days)") +
	scale_x_continuous() +
	scale_y_continuous(limits = c(-90, 900)) +
	theme_minimal() +
	theme(text = element_text(size = 14),
	      axis.title.x = element_text(margin = margin(t = 20)),
	      axis.text.x = element_text(size = 0),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	facet_grid(0~hpv_type_wes_wgs, scales = "free_x", space = "free_x")

pdf(file = "../res/Sample_by_Time_Point.pdf", width = 9, height = 4)
print(plot_)
dev.off()

#==================================================
## Number of variants per sample
#==================================================
smry_ = manifest %>%
	dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	dplyr::mutate(patient_id_mskcc = case_when(
			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
			TRUE ~ patient_id_mskcc
	)) %>%
	dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
	dplyr::mutate(hpv_type_wes_wgs = case_when(
		   is.na(hpv_type_wes_wgs) ~ "Unknown",
		   TRUE ~ hpv_type_wes_wgs
	)) %>%
	dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::summarize(min_vars_input = min(`MRD-Landmark_Input_Variant_Count`),
			 mean_vars_input = mean(`MRD-Landmark_Input_Variant_Count`),
			 max_vars_input = max(`MRD-Landmark_Input_Variant_Count`),
			 min_vars_pf = min(`MRD-Landmark_Input_Variants_Passing_Filter`),
			 mean_vars_pf = mean(`MRD-Landmark_Input_Variants_Passing_Filter`),
			 max_vars_pf = max(`MRD-Landmark_Input_Variants_Passing_Filter`)) %>%
	pander::pander()

#==================================================
## Pre-treatment ctDNA negative patients
#==================================================
smry_baseline__mrd = manifest %>%
	       	     dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
		     dplyr::mutate(patient_id_mskcc = case_when(
			     			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
			     			TRUE ~ patient_id_mskcc
		     )) %>%
		     dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
		     dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
		     dplyr::mutate(hpv_type_wes_wgs = case_when(
			     			is.na(hpv_type_wes_wgs) ~ "Unknown",
			     			TRUE ~ hpv_type_wes_wgs
		     )) %>%
		     dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
		     dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
		     dplyr::filter(timepoint_days_since_start_of_RT<=0) %>%
		     dplyr::group_by(patient_name) %>%
		     dplyr::summarize(Is_ctDNA = any(`MRD-Landmark_Result` == "PRESENT")) %>%
		     .[["Is_ctDNA"]] %>%
		     sum()

smry_baseline__hpv = readr::read_tsv(file = "../res/posterior_probability_all.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	             readr::type_convert() %>%
		     dplyr::filter(timepoint_days_since_start_of_RT<=0) %>%
		     dplyr::group_by(patient_name) %>%
		     dplyr::summarize(Is_ctDNA = any(Is_ctDNA == "+ve")) %>%
		     .[["Is_ctDNA"]] %>%
		     sum()

smry_baseline__mrd_hpv = manifest %>%
	       	     	 dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
			 dplyr::mutate(patient_id_mskcc = case_when(
				 			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
				 			TRUE ~ patient_id_mskcc
			 )) %>%
			 dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
			 dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
			 dplyr::mutate(hpv_type_wes_wgs = case_when(
				 		is.na(hpv_type_wes_wgs) ~ "Unknown",
			     			TRUE ~ hpv_type_wes_wgs
			 )) %>%
			 dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
			 dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
			 dplyr::filter(timepoint_days_since_start_of_RT<=0) %>%
			 dplyr::group_by(patient_name) %>%
			 dplyr::summarize(Is_ctDNA_MRD = any(`MRD-Landmark_Result` == "PRESENT")) %>%
			 dplyr::left_join(readr::read_tsv(file = "../res/posterior_probability_all.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
	             			  readr::type_convert() %>%
					  dplyr::filter(timepoint_days_since_start_of_RT<=0) %>%
					  dplyr::group_by(patient_name) %>%
					  dplyr::summarize(Is_ctDNA_HPV = any(Is_ctDNA == "+ve")),
					  by = "patient_name") %>%
			 dplyr::mutate(Is_ctDNA = Is_ctDNA_MRD | Is_ctDNA_HPV) %>%
		     	 .[["Is_ctDNA"]] %>%
			 sum()

#==================================================
# ctDNA fraction by time point
#==================================================
smry_ctdna = manifest %>%
	     dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	     dplyr::mutate(patient_id_mskcc = case_when(
		     			is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
		     			TRUE ~ patient_id_mskcc
	     )) %>%
	     dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	     dplyr::left_join(hpv_smry, by = "patient_id_mskcc") %>%
	     dplyr::mutate(hpv_type_wes_wgs = case_when(
		   			is.na(hpv_type_wes_wgs) ~ "Unknown",
		    			TRUE ~ hpv_type_wes_wgs
	     )) %>%
	     dplyr::left_join(mrd_smry, by = "sample_uuid") %>%
	     dplyr::filter(timepoint_days_since_start_of_RT<=0 | timepoint_days_since_start_of_RT>=730) %>%
	     dplyr::mutate(time_point = case_when(
		     timepoint_days_since_start_of_RT<=0 ~ "<=0",
		     timepoint_days_since_start_of_RT>=730 ~ ">=730"
	     ))

plot_ = smry_ctdna %>%
	ggplot(aes(x = time_point, y = `MRD-Monitoring_MRD_Quantification`)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	geom_jitter(stat = "identity", fill = "salmon", width = .15, height = 0, shape = 21, alpha = .75, size = 2) +
	xlab("Time (days)") +
	ylab("ctDNA fraction (?)") +
	scale_x_discrete() +
	scale_y_sqrt(breaks = c(0.00625, .025, .05, .1, .2, .3),
		     labels = scientific_10) +
	geom_signif(stat = "signif",
		    comparisons = list(c("<=0", ">=730")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = sqrt(.4),
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)))

pdf(file = "../res/ctDNA_by_Time_Point.pdf", width = 4, height = 5)
print(plot_)
dev.off()


#########################################################
# (+) Samples with mean AF > 5% from MRD assay
# (+) Samples with HPV subtype in assay
# (+) Samples with +ve MRD assay
# (-) Samples with mean AF < 0.01% from MRD assay
# (-) Samples with max AF < 0.1% from MRD assay
# (-) Patients with no nodal dissection ≥ 2 years
# (-) No duplicate patients
# (+) Samples with -ve MRD assay
#########################################################
manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(!is.na(bam_file_name_hpv)) %>%
	   dplyr::mutate(sample_uuid = paste0(sample_id_mskcc, "-", sample_id_invitae))

manifest = manifest %>%
	   dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	   dplyr::mutate(patient_id_mskcc = case_when(
		   is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
		   TRUE ~ patient_id_mskcc
	   )) %>%
	   dplyr::mutate(sample_name = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	   dplyr::left_join(hpv_smry, by = "patient_id_mskcc")

smry_t_pos = mutation_smry %>%
	     dplyr::filter(FILTER == "PASS") %>%
	     dplyr::group_by(Tumor_Sample_Barcode) %>%
	     dplyr::summarize(mean_af = mean(t_maf),
			      max_af = max(t_maf)) %>%
	     dplyr::ungroup() %>%
	     dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
	     ## Mean AF > 5%
	     dplyr::filter(mean_af > (5/100)) %>%
	     dplyr::left_join(manifest, by = "sample_name") %>%
	     ## Known HPV
	     dplyr::filter(hpv_type_wes_wgs %in% names(target_contigs)) %>%
	     ## Positive MRD assay
	     dplyr::left_join(mrd_smry %>%
			      dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			      dplyr::select(sample_name, `MRD-Landmark_Result`, `MRD-Monitoring_MRD_Quantification`), by = "sample_name") %>%
	     dplyr::filter(`MRD-Landmark_Result` == "PRESENT") %>%
	     dplyr::mutate(Is_ctDNA = "+ve")

smry_t_neg = mutation_smry %>%
	     dplyr::filter(FILTER == "PASS") %>%
	     dplyr::group_by(Tumor_Sample_Barcode) %>%
	     dplyr::summarize(mean_af = mean(t_maf),
			      max_af = max(t_maf)) %>%
	     dplyr::ungroup() %>%
	     dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
	     ## Mean AF < 0.01%
	     dplyr::filter(mean_af < (0.01/100)) %>%
	     ## Max AF < 0.1%
	     dplyr::filter(max_af < (0.1/100)) %>%
	     dplyr::left_join(manifest, by = "sample_name") %>%
	     dplyr::left_join(nodal_dissection_smry, by = "patient_id_mskcc") %>%
	     ## No nodal dissection ≥ 2 years
	     dplyr::filter(nd_event == 0 & timepoint_days_since_start_of_RT > 730) %>%
	     ## No duplicate patients
	     dplyr::left_join(mutation_smry %>%
	     		      dplyr::filter(FILTER == "PASS") %>%
			      dplyr::group_by(Tumor_Sample_Barcode) %>%
			      dplyr::summarize(mean_af = mean(t_maf),
					       max_af = max(t_maf)) %>%
			      dplyr::ungroup() %>%
			      dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
			      ## Mean AF < 0.01%
			      dplyr::filter(mean_af < (0.01/100)) %>%
			      ## Max AF < 0.1%
			      dplyr::filter(max_af < (0.1/100)) %>%
			      dplyr::left_join(manifest, by = "sample_name") %>%
			      dplyr::left_join(nodal_dissection_smry, by = "patient_id_mskcc") %>%
			      ## No nodal dissection ≥ 2 years
			      dplyr::filter(nd_event == 0 & timepoint_days_since_start_of_RT > 730) %>%
			      dplyr::group_by(patient_id_mskcc) %>%
			      dplyr::summarize(timepoint_days_since_start_of_RT = max(timepoint_days_since_start_of_RT)) %>%
			      dplyr::mutate(is_selected = TRUE), by = c("patient_id_mskcc", "timepoint_days_since_start_of_RT")) %>%
	     dplyr::filter(!is.na(is_selected)) %>%
	     dplyr::select(-is_selected) %>%
	     ## Negative MRD assay
	     dplyr::left_join(mrd_smry %>%
			      dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			      dplyr::select(sample_name, `MRD-Landmark_Result`, `MRD-Monitoring_MRD_Quantification`), by = "sample_name") %>%
	     dplyr::filter(`MRD-Landmark_Result` == "ABSENT") %>%
	     dplyr::mutate(Is_ctDNA = "-ve")

smry_ft = dplyr::bind_rows(smry_t_pos %>%
			   dplyr::select(all_of(intersect(colnames(smry_t_pos), colnames(smry_t_neg)))),
			   smry_t_neg %>%
			   dplyr::select(all_of(intersect(colnames(smry_t_pos), colnames(smry_t_neg)))))

plot_ = smry_ft %>%
	ggplot(aes(x = Is_ctDNA, y = 100*mean_af)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	geom_jitter(stat = "identity", fill = "salmon", width = .15, height = 0, shape = 21, alpha = .75, size = 2) +
	xlab("ctDNA") +
	ylab("Mean AF (%)") +
	scale_x_discrete() +
	scale_y_sqrt(breaks = c(0.00625, .025, .05, .1, .2, .3, .4)*100,
		     labels = scientific_10) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = sqrt(.4*100),
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)))

pdf(file = "../res/Mean_AF_by_Time_Point.pdf", width = 4, height = 5)
print(plot_)
dev.off()

plot_ = smry_ft %>%
	ggplot(aes(x = Is_ctDNA, y = 100*max_af)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	geom_jitter(stat = "identity", fill = "salmon", width = .15, height = 0, shape = 21, alpha = .75, size = 2) +
	xlab("ctDNA") +
	ylab("Max AF (%)") +
	scale_x_discrete() +
	scale_y_sqrt(breaks = c(0.00625, .025, .05, .1, .2, .4, .6)*100,
		     labels = scientific_10) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = sqrt(.6*100),
		    tip_length = 0.01) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)))

pdf(file = "../res/Max_AF_by_Time_Point.pdf", width = 4, height = 5)
print(plot_)
dev.off()

plot_ = smry_ft %>%
	ggplot(aes(x = mean_af, y = `MRD-Monitoring_MRD_Quantification`, shape = Is_ctDNA, fill = Is_ctDNA)) +
	geom_abline(intercept = 0, slope = 1, color = "lightgrey", alpha = .65, size = 1.25) +
	geom_smooth(data = smry_ft,
		    mapping = aes(x = mean_af, y = `MRD-Monitoring_MRD_Quantification`),
		    stat = "smooth", method = "glm", formula = y ~ x, color = "goldenrod3", se = FALSE, size = 1.25, alpha = .75, inherit.aes = FALSE) +
	geom_point(stat = "identity", size = 3, alpha = 1) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22)) +
	xlab("Mean AF (%)") +
	ylab("ctDNA farction (?)") +
	scale_x_continuous() +
	scale_y_continuous() +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(shape = guide_legend(title = "ctDNA"),
	       fill =  guide_legend(title = "ctDNA"))

pdf(file = "../res/Mean_AF_by_ctDNA_Fraction.pdf", width = 5.5, height = 5)
print(plot_)
dev.off()

plot_ = smry_ft %>%
	ggplot(aes(x = max_af, y = `MRD-Monitoring_MRD_Quantification`, shape = Is_ctDNA, fill = Is_ctDNA)) +
	geom_abline(intercept = 0, slope = 1, color = "lightgrey", alpha = .65, size = 1.25) +
	geom_smooth(data = smry_ft,
		    mapping = aes(x = max_af, y = `MRD-Monitoring_MRD_Quantification`),
		    stat = "smooth", method = "glm", formula = y ~ x, color = "goldenrod3", se = FALSE, size = 1.25, alpha = .75, inherit.aes = FALSE) +
	geom_point(stat = "identity", size = 3, alpha = 1) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_shape_manual(values = c(21, 22)) +
	xlab("Max AF (%)") +
	ylab("ctDNA farction (?)") +
	scale_x_continuous() +
	scale_y_continuous() +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(shape = guide_legend(title = "ctDNA"),
	       fill =  guide_legend(title = "ctDNA"))

pdf(file = "../res/Max_AF_by_ctDNA_Fraction.pdf", width = 5.5, height = 5)
print(plot_)
dev.off()
