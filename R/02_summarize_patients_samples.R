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

clinical = readr::read_tsv(file = url_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::rename(Tumor_Sample_Barcode = tumor_id_wes) %>%
	   tidyr::drop_na(Tumor_Sample_Barcode) %>%
	   dplyr::select(-additional_crt_post_neck_dissection) %>%
	   dplyr::mutate(is_selected = 1)

tumor_vars = readr::read_tsv(file = url_tumor_variants, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert()
	     
#==================================================
# Number of unique patients with MRD assay
#==================================================
manifest %>%
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
dplyr::mutate(patient_number = 1:n()) %>%
dplyr::select(patient_number, patient_id = patient_id_mskcc, n) %>%
pander::pander(caption = "Number of unique patients with MRD assay")

#==================================================
# Number of samples per patient
#==================================================
manifest %>%
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
pander::pander(caption = "Number of samples per patient")

#==================================================
# Number of patient per HPV subtype
#==================================================
manifest %>%
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
pander::pander(caption = "Number of patient per HPV subtype")

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
	ggplot(aes(x = hpv_type_wes_wgs, fill = hpv_type_wes_wgs)) +
	geom_bar(stat = "count", color = "black", alpha = 1) +
	scale_fill_manual(values = c("HPV-16" = "#377eb8",
				     "HPV-18" = "#e41a1c",
				     "HPV-33" = "#e41a1c",
				     "HPV-35" = "#e41a1c",
				     "HPV-58" = "#e41a1c",
				     "Unknown" = "#737373")) +
	xlab("") +
	ylab("Number of Patients") +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 90),
			   breaks = c(0, 15, 30, 45, 60, 75, 90),
			   labels = c(0, "", 30, "", 60, "", 90)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 11),
	      axis.title.y = element_text(margin = margin(r = 20), size = 11),
	      axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust = 1),
	      axis.text.y = element_text(size = 11)) +
	guides(fill = FALSE)

pdf(file = "../res/Number_Patients_by_HPV_Type.pdf", width = 2.5, height = 3.25)
print(plot_)
dev.off()

#==================================================
# Median number of variants per patient
#==================================================
tumor_vars %>%
dplyr::left_join(clinical, by = "Tumor_Sample_Barcode") %>%
tidyr::drop_na(is_selected) %>%
dplyr::group_by(Tumor_Sample_Barcode) %>%
dplyr::summarize(n = n()) %>%
.[["n"]] %>%
median() %>%
pander::pander(caption = "Median number of varaints per patient")

#==================================================
# Number of variants per sample
#==================================================
manifest %>%
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
		 median_vars_input = median(`MRD-Landmark_Input_Variant_Count`),
		 mean_vars_input = mean(`MRD-Landmark_Input_Variant_Count`),
		 max_vars_input = max(`MRD-Landmark_Input_Variant_Count`),
		 min_vars_pf = min(`MRD-Landmark_Input_Variants_Passing_Filter`),
		 median_vars_pf = median(`MRD-Landmark_Input_Variants_Passing_Filter`),
		 mean_vars_pf = mean(`MRD-Landmark_Input_Variants_Passing_Filter`),
		 max_vars_pf = max(`MRD-Landmark_Input_Variants_Passing_Filter`)) %>%
reshape2::melt() %>%
dplyr::rename(Statistic = variable, `#` = value) %>%
pander::pander(caption = "Number of variants per sample")

#==================================================
# ctDNA HPV-MRD patient cross-tab
#==================================================
smry__mrd = manifest %>%
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
	    dplyr::mutate(timepoint_weeks_since_start_of_RT = floor(timepoint_days_since_start_of_RT/7)) %>%
	    dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
		    		timepoint_weeks_since_start_of_RT < 0 ~ "Pre-treatment",
		    		TRUE ~ paste0("wk", timepoint_weeks_since_start_of_RT)
	    )) %>%
	    dplyr::group_by(timepoint_weeks_since_start_of_RT, patient_name) %>%
	    dplyr::summarize(MRD = any(`MRD-Landmark_Result` == "PRESENT")) %>%
	    dplyr::ungroup()

if (file.exists("../res/Posterior_Probability_ALL.txt")) {
	smry__hpv = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		    readr::type_convert() %>%
		    dplyr::mutate(timepoint_weeks_since_start_of_RT = floor(timepoint_days_since_start_of_RT/7)) %>%
	    	    dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
			    	timepoint_weeks_since_start_of_RT < 0 ~ "Pre-treatment",
			    	TRUE ~ paste0("wk", timepoint_weeks_since_start_of_RT)
		    )) %>%
		    dplyr::group_by(timepoint_weeks_since_start_of_RT, patient_name) %>%
		    dplyr::summarize(HPV = any(Is_ctDNA == "+ve")) %>%
		    dplyr::ungroup()
}

smry_ = smry__mrd %>%
	dplyr::full_join(smry__hpv, by = c("timepoint_weeks_since_start_of_RT", "patient_name")) %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3", "wk5")) %>%
	dplyr::group_by(timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(TP = sum(MRD & HPV),
			 FP = sum(!MRD & HPV),
			 FN = sum(MRD & !HPV),
			 TN = sum(!MRD & !HPV)) %>%
	reshape2::melt() %>%
	reshape2::dcast(variable ~ timepoint_weeks_since_start_of_RT, value.var = "value") %>%
	dplyr::select(-variable) %>%
	as.list() %>%
	lapply(function(x) { y = matrix(x, 2, 2, byrow = TRUE);
			     rownames(y) = colnames(y) = c("+ve", "-ve");
			     return(y)})

pdf(file = "../res/ctDNA_Patient_Table_by_Time_Point_All_Absolute.pdf", width = 15, height = 10)
corrplot(corr = do.call(cbind, smry_),
	 method = "circle",
	 type = "full",
	 col = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))(25),
	 is.corr = FALSE,
	 diag = TRUE,
	 addgrid.col = NA,
	 addCoef.col = "white",
	 order = "original",
	 number.cex = 1,
	 number.font = 1,
	 tl.cex = 1,
	 tl.col = "black",
	 cl.pos = "r",
	 cl.ratio = .1,
	 cl.length = 7)
dev.off()

pdf(file = "../res/ctDNA_Patient_Table_by_Time_Point_All_Percent.pdf", width = 15, height = 10)
corrplot(corr = do.call(cbind, lapply(smry_, function(x) {round(100*x/sum(x))})),
	 method = "circle",
	 type = "full",
	 col = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))(25),
	 is.corr = FALSE,
	 diag = TRUE,
	 addgrid.col = NA,
	 addCoef.col = "white",
	 order = "original",
	 number.cex = 1,
	 number.font = 1,
	 tl.cex = 1,
	 tl.col = "black",
	 cl.pos = "r",
	 cl.ratio = .1,
	 cl.length = 7)
dev.off()

#==================================================
# ctDNA HPV-MRD patient marginal frequency
#==================================================
smry__mrd = manifest %>%
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
	    dplyr::mutate(timepoint_weeks_since_start_of_RT = floor(timepoint_days_since_start_of_RT/7)) %>%
	    dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
		    		timepoint_weeks_since_start_of_RT < 0 ~ "Pre-treatment",
		    		TRUE ~ paste0("wk", timepoint_weeks_since_start_of_RT)
	    )) %>%
	    dplyr::group_by(timepoint_weeks_since_start_of_RT, patient_name) %>%
	    dplyr::summarize(Is_ctDNA = any(`MRD-Landmark_Result` == "PRESENT")) %>%
	    dplyr::ungroup() %>%
	    dplyr::group_by(timepoint_weeks_since_start_of_RT) %>%
	    dplyr::summarize(N = n(),
		`Is_ctDNA_%` = sum(Is_ctDNA)) %>%
	    dplyr::mutate(`Is_ctDNA_%` = 100*`Is_ctDNA_%`/N)

if (file.exists("../res/Posterior_Probability_ALL.txt")) {
	smry__hpv = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		    readr::type_convert() %>%
		    dplyr::mutate(timepoint_weeks_since_start_of_RT = floor(timepoint_days_since_start_of_RT/7)) %>%
	    	    dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
			    	timepoint_weeks_since_start_of_RT < 0 ~ "Pre-treatment",
			    	TRUE ~ paste0("wk", timepoint_weeks_since_start_of_RT)
		    )) %>%
		    dplyr::group_by(timepoint_weeks_since_start_of_RT, patient_name) %>%
		    dplyr::summarize(Is_ctDNA = any(Is_ctDNA == "+ve")) %>%
		    dplyr::ungroup() %>%
		    dplyr::group_by(timepoint_weeks_since_start_of_RT) %>%
		    dplyr::summarize(N = n(),
				`Is_ctDNA_%` = sum(Is_ctDNA)) %>%
		    dplyr::mutate(`Is_ctDNA_%` = 100*`Is_ctDNA_%`/N)

	smry__mrd_hpv = manifest %>%
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
			dplyr::left_join(readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
					 readr::type_convert() %>%
					 dplyr::select(sample_uuid = sample_name, Is_ctDNA), by = "sample_uuid") %>%
			dplyr::mutate(timepoint_weeks_since_start_of_RT = floor(timepoint_days_since_start_of_RT/7)) %>%
	    	        dplyr::mutate(timepoint_weeks_since_start_of_RT = case_when(
					timepoint_weeks_since_start_of_RT < 0 ~ "Pre-treatment",
					TRUE ~ paste0("wk", timepoint_weeks_since_start_of_RT)
			)) %>%
			dplyr::group_by(timepoint_weeks_since_start_of_RT, patient_name) %>%
			dplyr::summarize(Is_ctDNA_MRD = any(`MRD-Landmark_Result` == "PRESENT"),
					 Is_ctDNA_HPV = any(Is_ctDNA == "+ve")) %>%
			dplyr::ungroup() %>%
			dplyr::group_by(timepoint_weeks_since_start_of_RT) %>%
			dplyr::summarize(N = n(),
					 `Is_ctDNA_%` = sum(Is_ctDNA_MRD | Is_ctDNA_HPV)) %>%
			dplyr::mutate(`Is_ctDNA_%` = 100*`Is_ctDNA_%`/N)
}

plot_ = smry__mrd %>%
	dplyr::rename(PCM = `Is_ctDNA_%`) %>%
	dplyr::left_join(smry__hpv %>%
			 dplyr::rename(HPV = `Is_ctDNA_%`),
			 by = c("timepoint_weeks_since_start_of_RT", "N")) %>%
	dplyr::left_join(smry__mrd_hpv %>%
			 dplyr::rename(Combined = `Is_ctDNA_%`),
			 by = c("timepoint_weeks_since_start_of_RT", "N")) %>%
	reshape2::melt() %>%
	dplyr::filter(variable != "N") %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT %in% c("Pre-treatment", "wk1", "wk2", "wk3", "wk5")) %>%
	dplyr::mutate(timepoint_weeks_since_start_of_RT = gsub(pattern = "wk", replacement = "Week ", x = timepoint_weeks_since_start_of_RT, fixed = TRUE)) %>%
	ggplot(aes(x = timepoint_weeks_since_start_of_RT, y = value, fill = variable)) +
	geom_bar(stat = "identity", position = "dodge", color = "white", width = .75) +
	scale_fill_brewer(type = "qual", palette = 7) +
	xlab("") +
	ylab("Fraction of Patients (%)") +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 100),
			   breaks = seq(from = 0, to = 100, by = 10),
			   labels = c(0, "", 20, "", 40, "", 60, "", 80, "", 100)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = guide_legend(title = "Assay"))

pdf(file = "../res/ctDNA_Patient_Summary_by_Time_Point_All.pdf", width = 5.5, height = 3.25)
print(plot_)
dev.off()

plot_ = smry__mrd %>%
	dplyr::rename(PCM = `Is_ctDNA_%`) %>%
	dplyr::left_join(smry__hpv %>%
			 dplyr::rename(HPV = `Is_ctDNA_%`),
			 by = c("timepoint_weeks_since_start_of_RT", "N")) %>%
	dplyr::left_join(smry__mrd_hpv %>%
			 dplyr::rename(Combined = `Is_ctDNA_%`),
			 by = c("timepoint_weeks_since_start_of_RT", "N")) %>%
	reshape2::melt() %>%
	dplyr::filter(variable != "N") %>%
	dplyr::filter(timepoint_weeks_since_start_of_RT == "Pre-treatment") %>%
	ggplot(aes(x = variable, y = value, fill = variable)) +
	geom_bar(stat = "identity", position = "dodge", color = "white", width = 1) +
	scale_fill_brewer(type = "qual", palette = 7) +
	xlab("") +
	ylab("Fraction of Patients (%)") +
	scale_x_discrete() +
	scale_y_continuous(limits = c(0, 100),
			   breaks = seq(from = 0, to = 100, by = 10),
			   labels = c(0, "", 20, "", 40, "", 60, "", 80, "", 100)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = FALSE)

pdf(file = "../res/ctDNA_Patient_Summary_by_Time_Point_Pre-treatment.pdf", width = 1.65, height = 3.25)
print(plot_)
dev.off()
		 
#==================================================
# (+) Samples with mean AF > 5% from MRD assay
# (+) Samples with HPV subtype in assay
# (+) Samples with +ve MRD assay
# (-) Samples with mean AF < 0.01% from MRD assay
# (-) Samples with max AF < 0.1% from MRD assay
# (-) Patients with no nodal dissection ≥ 2 years
# (-) No duplicate patients
# (+) Samples with -ve MRD assay
#==================================================
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
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::bind_rows(
		mutation_smry %>%
		dplyr::filter(FILTER == "PASS") %>%
		dplyr::group_by(Tumor_Sample_Barcode) %>%
		dplyr::summarize(mean_af = mean(t_maf),
				 max_af = max(t_maf)) %>%
		dplyr::ungroup() %>%
		dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
		dplyr::left_join(manifest, by = "sample_name") %>%
		dplyr::filter(!(hpv_type_wes_wgs %in% c("HPV-16", "Unknown"))) %>%
		dplyr::left_join(mrd_smry %>%
				 dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
				 dplyr::select(sample_name, `MRD-Landmark_Result`, `MRD-Monitoring_MRD_Quantification`), by = "sample_name") %>%
		dplyr::filter(`MRD-Landmark_Result` == "PRESENT") %>%
		dplyr::mutate(Is_ctDNA = "Unknown") %>%
		dplyr::select(all_of(colnames(smry_ft)))
	) %>%
	ggplot(aes(x = Is_ctDNA, y = 100*`MRD-Monitoring_MRD_Quantification`, fill = Is_ctDNA)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white") +
	geom_jitter(stat = "identity", width = .15, height = 0, shape = 21, alpha = .85, size = 3) +
	scale_fill_manual(values = c("+ve" = "#377eb8",
				     "-ve" = "#377eb8",
				     "Unknown" = "#e41a1c")) +
	xlab("") +
	ylab("ctDNA Fraction (%)") +
	scale_x_discrete(breaks = c("+ve", "-ve", "Unknown"),
			 labels = c("+ve", "-ve", "Outgroup")) +
	scale_y_continuous(limits = c(0, 32),
			   breaks = c(0, 5, 10, 15, 20, 25, 30),
			   labels = c(0, "", 10, "", 20, "", 30)) +
	geom_signif(stat = "signif",
		    comparisons = list(c("+ve", "-ve")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = 31, vjust = -.25,
		    tip_length = 0.01) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 10),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = FALSE)

pdf(file = "../res/ctDNA_Fraction_by_Time_Point.pdf", width = 2.5, height = 3.25)
print(plot_)
dev.off()
