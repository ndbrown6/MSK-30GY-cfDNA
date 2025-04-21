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
	dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(mean_af = 100*mean(mean_af+1e-5, na.rm = TRUE)) %>%
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
	dplyr::full_join(idx_metrics_ft %>%
			 dplyr::filter(hpv_type_wes_wgs == "HPV-16") %>%
			 dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
			 dplyr::summarize(aligned_reads = mean(aligned_reads)) %>%
			 reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
					 fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "aligned_reads") %>%
			 dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
			 dplyr::mutate(rel_HPV_wk1 = wk1/`Pre-treatment`,
		      		       rel_HPV_wk2 = wk2/`Pre-treatment`,
				       rel_HPV_wk3 = wk3/`Pre-treatment`) %>%
			 dplyr::rename(abs_HPV_wk0 = `Pre-treatment`,
				       abs_HPV_wk1 = wk1,
				       abs_HPV_wk2 = wk2,
				       abs_HPV_wk3 = wk3),
			 by = "patient_id_mskcc") %>%
	dplyr::full_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::select(patient_id_mskcc,
		      composite_end_point,
		      abs_AF_wk0, abs_AF_wk1, abs_AF_wk2, abs_AF_wk3,
		      rel_AF_wk1, rel_AF_wk2, rel_AF_wk3,
		      abs_HPV_wk0, abs_HPV_wk1, abs_HPV_wk2, abs_HPV_wk3,
		      rel_HPV_wk1, rel_HPV_wk2, rel_HPV_wk3,
		      abs_MRI_wk0 = MRI_rawdata_wk0,
		      abs_MRI_wk1 = MRI_rawdata_wk2,
		      abs_MRI_wk2 = MRI_rawdata_wk3,
		      abs_MRI_wk3 = MRI_rawdata_wk4) %>%
	dplyr::mutate(rel_MRI_wk1 = abs_MRI_wk1/abs_MRI_wk0,
		      rel_MRI_wk2 = abs_MRI_wk2/abs_MRI_wk0,
		      rel_MRI_wk3 = abs_MRI_wk1/abs_MRI_wk0) %>%
	dplyr::mutate(composite_end_point = case_when(
		composite_end_point ~ 1,
		!composite_end_point ~ 0
	)) %>%
	tidyr::drop_na(composite_end_point)

paramx = rpart::rpart.control(minsplit = 4,
			      minbucket = 1,
			      maxdepth = 5,
			      maxcompete = 1000,
			      maxsurrogate = 100,
			      usesurrogate = 2,
			      cp = -Inf,
			      xval = 100)

paramy = rpart::rpart.control(minsplit = 5,
			      minbucket = 2,
			      maxdepth = 5,
			      maxcompete = 1000,
			      maxsurrogate = 100,
			      usesurrogate = 2,
			      cp = -Inf,
			      xval = 100)

'Class_NoCv' <- function(data) {
	fit = rpart(composite_end_point ~ ., data = data, method = "class", control = paramx)
	prd = predict(object = fit, newdata = data, type = "class")
	conf = table(data$composite_end_point, prd)
	mce = 100*(1-sum(diag(conf))/sum(conf))
	return(invisible(mce))
}

'Class_10Cv' <- function(data) {
	tidy_data = data %>% tidyr::drop_na()
	mce = list()
	for (j  in 1:100) {
		folds = split(tidy_data, cut(sample(1:nrow(tidy_data)),10))
		mce[[j]] = rep(NA, length(folds))
		for (i in 1:length(folds)) {
			test = ldply(folds[i], data.frame, .id = NULL)
			train = ldply(folds[-i], data.frame, .id = NULL)
			tmp.fit = rpart(composite_end_point ~ ., data = train, method = "class", control = paramy)
			tmp.prd = predict(tmp.fit, newdata = test, type = "class")
			conf = table(test$composite_end_point, tmp.prd)
			mce[[j]][i] <- 1-sum(diag(conf))/sum(conf)
		}
	}
	mce = 100*unlist(lapply(mce, mean))
	return(invisible(mce))
}

set.seed(12345)
Training_Error = list()
Test_Error = list()

#==================================================
# Relative HPV
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, rel_HPV_wk1, rel_HPV_wk2, rel_HPV_wk3) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point))

Training_Error[[1]] = Class_NoCv(data = data)
Test_Error[[1]] = Class_10Cv(data = data)

#==================================================
# Absolute HPV
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, abs_HPV_wk0, abs_HPV_wk1, abs_HPV_wk2, abs_HPV_wk3) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point))

Training_Error[[2]] = Class_NoCv(data = data)
Test_Error[[2]] = Class_10Cv(data = data)

#==================================================
# Relative AF
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, rel_AF_wk1, rel_AF_wk2, rel_AF_wk3) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point))

Training_Error[[3]] = Class_NoCv(data = data)
Test_Error[[3]] = Class_10Cv(data = data)

#==================================================
# Absolute AF
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, abs_AF_wk0, abs_AF_wk1, abs_AF_wk2, abs_AF_wk3) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point))

Training_Error[[4]] = Class_NoCv(data = data)
Test_Error[[4]] = Class_10Cv(data = data)

#==================================================
# Relative MRI
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, rel_MRI_wk1, rel_MRI_wk2, rel_MRI_wk3) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point))

Training_Error[[5]] = Class_NoCv(data = data)
Test_Error[[5]] = Class_10Cv(data = data)

#==================================================
# Absolute MRI
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, abs_MRI_wk0, abs_MRI_wk1, abs_MRI_wk2, abs_MRI_wk3) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point))

Training_Error[[6]] = Class_NoCv(data = data)
Test_Error[[6]] = Class_10Cv(data = data)

#==================================================
# Absolute AF + Absolute MRI
#==================================================
data = data_ %>%
       dplyr::select(composite_end_point, abs_AF_wk0, abs_AF_wk1, abs_AF_wk2, abs_AF_wk3, abs_MRI_wk0, abs_MRI_wk1, abs_MRI_wk2, abs_MRI_wk3) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point))

Training_Error[[7]] = Class_NoCv(data = data)
Test_Error[[7]] = Class_10Cv(data = data)

names(Training_Error) = c("Relative cfDNA HPV\nReads", "Absolute cfDNA HPV\nReads",
			  "Relative ctDNA\nFraction", "Absolute ctDNA\nFraction",
			  "Absolute MRI\nVolume", "Relative MRI\nVolume",
			  "ctDNA Fraction\n+ MRI Volume")

names(Test_Error) = c("Relative cfDNA HPV\nReads", "Absolute cfDNA HPV\nReads",
		      "Relative ctDNA\nFraction", "Absolute ctDNA\nFraction",
		      "Absolute MRI\nVolume", "Relative MRI\nVolume",
		      "ctDNA Fraction\n+ MRI Volume")


p1 = wilcox.test(x = Test_Error$`ctDNA Fraction\n+ MRI Volume`,
		 y = Test_Error$`Absolute ctDNA\nFraction`)$p.value

p2 = wilcox.test(x = Test_Error$`ctDNA Fraction\n+ MRI Volume`,
		 y = Test_Error$`Absolute MRI\nVolume`)$p.value

plot_ = Test_Error %>%
	dplyr::as_tibble() %>%
	reshape2::melt() %>%
	dplyr::mutate(category = case_when(
		grepl("Relative", variable) ~ "Relative",
		grepl("Absolute", variable) ~ "Absolute",
		TRUE ~ "Absolute"
	)) %>%
	dplyr::mutate(variable = gsub(pattern = "Relative |Absolute ", replacement = "", x = variable, perl = TRUE, fixed = FALSE)) %>%
	dplyr::mutate(category = factor(category, levels = c("Relative", "Absolute"), ordered = TRUE)) %>%
	dplyr::mutate(variable = factor(variable, levels = c("cfDNA HPV\nReads", "ctDNA\nFraction", "MRI\nVolume", "ctDNA Fraction\n+ MRI Volume"), ordered = TRUE)) %>%
	ggplot(aes(x = variable, y = value, color = category, group = variable:category)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA) +
	geom_jitter(stat = "identity", position = position_jitterdodge(jitter.width = .2), shape = 21, size = 2, fill = "white", alpha = .55) +
	geom_hline(yintercept = Training_Error$`ctDNA Fraction\n+ MRI Volume`, color = "red", linetype = 3) +
	geom_jitter(data = Training_Error %>%
		   	  dplyr::as_tibble() %>%
		   	  reshape2::melt() %>%
			  dplyr::mutate(category = case_when(
				  grepl("Relative", variable) ~ "Relative",
				  grepl("Absolute", variable) ~ "Absolute",
				  TRUE ~ "Absolute"
			  )) %>%
		   	  dplyr::mutate(variable = gsub(pattern = "Relative |Absolute ", replacement = "", x = variable, perl = TRUE, fixed = FALSE)) %>%
		   	  dplyr::mutate(category = factor(category, levels = c("Relative", "Absolute"), ordered = TRUE)) %>%
		   	  dplyr::mutate(variable = factor(variable, levels = c("cfDNA HPV\nReads", "ctDNA\nFraction", "MRI\nVolume", "ctDNA Fraction\n+ MRI Volume"), ordered = TRUE)),
		   mapping = aes(x = variable, y = value, fill = category, group = variable:category),
		   stat = "identity", position = position_jitterdodge(jitter.width = 0), shape = 8, inherit.aes = FALSE, show.legend = FALSE) +
	
	scale_color_brewer(type = "qual", palette = 7) +
	scale_x_discrete() +
	scale_y_continuous(limits = c(NA,100)) +
	xlab("") +
	ylab("10-Fold Misclass-\nification Error (%)") +
	geom_signif(annotation = formatC(p1, digits = 1),
		    y_position = 85, xmin = 2.3, xmax = 4, 
		    tip_length = c(0.03, 0.03), color = "black") +
	geom_signif(annotation = formatC(p2, digits = 1),
		    y_position = 75, xmin = 3.3, xmax = 4, 
		    tip_length = c(0.03, 0.03), color = "black") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
	      axis.text.y = element_text(size = 12)) +
	guides(color = guide_legend(title = "Relative or\nAbsolute\nvalue?"))

pdf(file = "../res/Decesion_Tree_Misclassification_Error.pdf", width = 4.25, height = 4.5)
print(plot_)
dev.off()
