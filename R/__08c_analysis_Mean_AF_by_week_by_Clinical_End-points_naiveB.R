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
	dplyr::summarize(mean_af = mean(mean_af+1e-5, na.rm = TRUE)) %>%
	reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	dplyr::mutate(rel_AF_wk1 = wk1/`Pre-treatment`,
		      rel_AF_wk2 = wk2/`Pre-treatment`,
		      rel_AF_wk3 = wk3/`Pre-treatment`) %>%
	dplyr::rename(abs_AF_wk0 = `Pre-treatment`,
		      abs_AF_wk1 = wk1,
		      abs_AF_wk2 = wk2,
		      abs_AF_wk3 = wk3) %>%
	dplyr::left_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::select(patient_id_mskcc,
		      composite_end_point,
		      abs_AF_wk0, abs_AF_wk1, abs_AF_wk2, abs_AF_wk3,
		      rel_AF_wk1, rel_AF_wk2, rel_AF_wk3,
		      abs_MRI_wk0 = MRI_rawdata_wk0,
		      abs_MRI_wk1 = MRI_rawdata_wk1,
		      abs_MRI_wk2 = MRI_rawdata_wk2,
		      abs_MRI_wk3 = MRI_rawdata_wk3,
		      abs_MRI_wk4 = MRI_rawdata_wk4,
		      abs_ADC_wk0 = ADC_Mean_wk0,
		      abs_ADC_wk1 = ADC_Mean_wk1,
		      abs_ADC_wk2 = ADC_Mean_wk2,
		      abs_ADC_wk3 = ADC_Mean_wk3,
		      abs_ADC_wk4 = ADC_Mean_wk4,
		      hpv_copynumber = hpv_panel_copy_number,
		      tumor_volume = plan_volume) %>%
	dplyr::mutate(composite_end_point = case_when(
		composite_end_point ~ 1,
		!composite_end_point ~ 0
	)) %>%
	tidyr::drop_na(composite_end_point)


'Class_NoCv' <- function(data) {
	fit = naiveBayes(composite_end_point ~ ., data = data)
	prd = predict(object = fit, newdata = data)
	conf = table(data$composite_end_point, prd)
	mce = 100*(1-sum(diag(conf))/sum(conf))
	return(invisible(mce))
}

'Class_10Cv' <- function(data) {
	mce = list()
	for (j  in 1:100) {
		folds = split(data, cut(sample(1:nrow(data)),10))
		mce[[j]] = rep(NA, length(folds))
		for (i in 1:length(folds)) {
			test = ldply(folds[i], data.frame, .id = NULL)
			train = ldply(folds[-i], data.frame, .id = NULL)
			tmp.fit = naiveBayes(composite_end_point ~ ., data = train)
			tmp.prd = predict(tmp.fit, newdata = test)
			conf = table(test$composite_end_point, tmp.prd)
			mce[[j]][i] <- 1-sum(diag(conf))/sum(conf)
		}
	}
	mce = 100*unlist(lapply(mce, mean))
	return(invisible(mce))
}

'ROC' <- function(data) {
	fit = naiveBayes(composite_end_point ~ ., data = data)
	prd = predict(object = fit, newdata = data, type = "raw")[,"1"]
	res = pROC::roc(predictor = prd, response = data$composite_end_point, ret = "all_coords")
	res = dplyr::tibble(sensitivity = res$sensitivities,
			    specificity = res$specificities)
	return(invisible(res))
}

set.seed(12345)
Training_Error = list()
Test_Error = list()
ROC_Curve = list()

#==================================================
# Relative AF
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, rel_AF_wk1, rel_AF_wk2, rel_AF_wk3) %>%
       dplyr::mutate(rel_AF_wk1 = 100*rel_AF_wk1,
		     rel_AF_wk2 = 100*rel_AF_wk2,
		     rel_AF_wk3 = 100*rel_AF_wk3) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point)) %>%
       tidyr::drop_na()

Training_Error[[1]] = Class_NoCv(data = data)
Test_Error[[1]] = Class_10Cv(data = data)
ROC_Curve[[1]] = ROC(data = data)

#==================================================
# Absolute AF
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, abs_AF_wk0, abs_AF_wk1, abs_AF_wk2, abs_AF_wk3) %>%
       dplyr::mutate(abs_AF_wk0 = 100*abs_AF_wk0,
       		     abs_AF_wk1 = 100*abs_AF_wk1,
		     abs_AF_wk2 = 100*abs_AF_wk2,
		     abs_AF_wk3 = 100*abs_AF_wk3) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point)) %>%
       tidyr::drop_na()

Training_Error[[2]] = Class_NoCv(data = data)
Test_Error[[2]] = Class_10Cv(data = data)
ROC_Curve[[2]] = ROC(data = data)

#==================================================
# MRI
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, abs_MRI_wk0, abs_MRI_wk1, abs_MRI_wk2) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point)) %>%
       tidyr::drop_na()

Training_Error[[3]] = Class_NoCv(data = data)
Test_Error[[3]] = Class_10Cv(data = data)
ROC_Curve[[3]] = ROC(data = data)

#==================================================
# ADC
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, abs_ADC_wk0, abs_ADC_wk1, abs_ADC_wk2) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point)) %>%
       tidyr::drop_na()

Training_Error[[4]] = Class_NoCv(data = data)
Test_Error[[4]] = Class_10Cv(data = data)
ROC_Curve[[4]] = ROC(data = data)

#==================================================
# Relative AF + MRI
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, rel_AF_wk1, rel_AF_wk2, rel_AF_wk3, abs_MRI_wk0, abs_MRI_wk1, abs_MRI_wk2) %>%
       tidyr::drop_na() %>%
       dplyr::mutate(rel_AF_wk1 = 100*scale(rel_AF_wk1, center = TRUE, scale = TRUE),
		     rel_AF_wk2 = 100*scale(rel_AF_wk2, center = TRUE, scale = TRUE),
		     rel_AF_wk3 = 100*scale(rel_AF_wk3, center = TRUE, scale = TRUE),
		     abs_MRI_wk0 = 100*scale(abs_MRI_wk0, center = TRUE, scale = TRUE),
       		     abs_MRI_wk1 = 100*scale(abs_MRI_wk1, center = TRUE, scale = TRUE),
       		     abs_MRI_wk2 = 100*scale(abs_MRI_wk2, center = TRUE, scale = TRUE)) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point))

Training_Error[[5]] = Class_NoCv(data = data)
Test_Error[[5]] = Class_10Cv(data = data)
ROC_Curve[[5]] = ROC(data = data)

#==================================================
# Absolute AF + MRI
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, abs_AF_wk0, abs_AF_wk1, abs_AF_wk2, abs_MRI_wk0, abs_MRI_wk1, abs_MRI_wk2) %>%
       tidyr::drop_na() %>%
       dplyr::mutate(abs_AF_wk0 = 100*scale(abs_AF_wk0, center = TRUE, scale = TRUE),
		     abs_AF_wk1 = 100*scale(abs_AF_wk1, center = TRUE, scale = TRUE),
		     abs_AF_wk2 = 100*scale(abs_AF_wk2, center = TRUE, scale = TRUE),
		     abs_MRI_wk0 = 100*scale(abs_MRI_wk0, center = TRUE, scale = TRUE),
		     abs_MRI_wk1 = 100*scale(abs_MRI_wk1, center = TRUE, scale = TRUE),
		     abs_MRI_wk2 = 100*scale(abs_MRI_wk2, center = TRUE, scale = TRUE)) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point))

Training_Error[[6]] = Class_NoCv(data = data)
Test_Error[[6]] = Class_10Cv(data = data)
ROC_Curve[[6]] = ROC(data = data)

names(Training_Error) = c("Relative Mean AF", "Absolute Mean AF", "MRI volume", "ADC",
			  "Relative Mean\nAF + MRI", "Absolute Mean\nAF + MRI")

names(Test_Error) = c("Relative Mean AF", "Absolute Mean AF", "MRI volume", "ADC",
		      "Relative Mean\nAF + MRI", "Absolute Mean\nAF + MRI")

names(ROC_Curve) = c("Relative Mean AF", "Absolute Mean AF", "MRI volume", "ADC",
		      "Relative Mean\nAF + MRI", "Absolute Mean\nAF + MRI")

errors_nb = list(Training_Error = Training_Error,
		 Test_Error = Test_Error,
		 ROC_Curve = ROC_Curve)

save(errors_nb, file = "../res/Classifier_Nb.RData")