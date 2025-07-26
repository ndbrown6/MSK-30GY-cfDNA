#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert()

mrd_smry = readr::read_tsv(file = url_mrd_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert()

bed_file = readr::read_tsv(file = url_bed_file, col_names = FALSE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::rename(chromosome = X1,
			 start = X2,
			 end = X3,
			 gene_name = X4) %>%
	   dplyr::select(-X5) %>%
	   dplyr::filter(gene_name != "E1^E4")

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

insert_size_metrics = readr::read_tsv(file = url_insert_metrics_by_gene, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert() %>%
		      dplyr::rename(sample_name = SAMPLE_NAME) %>%
		      dplyr::filter(TARGET_REGION != "NC001526.4:3372-4794")

mutation_smry = readr::read_tsv(file = url_mutation_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		readr::type_convert()

target_names = c("E1_1:1950", "E2_1891:2989", "E5_2985:3237", "E6_7124:7601", "E7_7603:7900", "L1_4774:6292")

#==================================================
# (+) Samples with mean AF > 5% from MRD assay
# (+) Samples with HPV subtype in assay
# (+) Samples with +ve MRD assay
#==================================================
smry_t_pos = mutation_smry %>%
	     dplyr::filter(FILTER == "PASS") %>%
	     dplyr::group_by(Tumor_Sample_Barcode) %>%
	     dplyr::summarize(mean_af = mean(t_maf)) %>%
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
			      dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
	     dplyr::filter(`MRD-Landmark_Result` == "PRESENT") %>%
	     dplyr::mutate(Is_ctDNA = "+ve")

#==================================================
# (-) Samples with mean AF < 0.01% from MRD assay
# (-) Samples with max AF < 0.1% from MRD assay
# (-) Patients with no nodal dissection ≥ 2 years
# (-) No duplicate patients
# (+) Samples with -ve MRD assay
#==================================================
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
			      dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
	     dplyr::filter(`MRD-Landmark_Result` == "ABSENT") %>%
	     dplyr::mutate(Is_ctDNA = "-ve")

smry_ft = dplyr::bind_rows(smry_t_pos %>%
			   dplyr::select(all_of(intersect(colnames(smry_t_pos), colnames(smry_t_neg)))),
			   smry_t_neg %>%
			   dplyr::select(all_of(intersect(colnames(smry_t_pos), colnames(smry_t_neg)))))

#==================================================
# Training set HPV-16 +ve
#==================================================
smry_ = insert_size_metrics %>%
	dplyr::left_join(mutation_smry %>%
			 dplyr::filter(FILTER == "PASS") %>%
			 dplyr::group_by(Tumor_Sample_Barcode) %>%
			 dplyr::summarize(mean_af = mean(t_maf),
					  median_af = median(t_maf),
					  max_af = max(t_maf)) %>%
			 dplyr::ungroup() %>%
			 dplyr::rename(sample_name = Tumor_Sample_Barcode), by = "sample_name") %>%
	dplyr::left_join(manifest, by = "sample_name") %>%
	dplyr::mutate(hpv_type_wes_wgs = ifelse(is.na(hpv_type_wes_wgs), "Unknown", hpv_type_wes_wgs)) %>%
	dplyr::mutate(chromosome = unlist(lapply(TARGET_REGION, function(x) { (strsplit(x = x, split = ":", fixed = TRUE)[[1]])[1]}))) %>%
	dplyr::mutate(start = unlist(lapply(TARGET_REGION, function(x) { (strsplit(x = x, split = ":|-", perl = TRUE)[[1]])[2]}))) %>%
	dplyr::mutate(end = unlist(lapply(TARGET_REGION, function(x) { (strsplit(x = x, split = ":|-", perl = TRUE)[[1]])[3]}))) %>%
	readr::type_convert() %>%
	dplyr::left_join(bed_file, by = c("chromosome", "start", "end")) %>%
	dplyr::mutate(chromosome = case_when(
			chromosome == "J04353.1" ~ "HPV-31",
			chromosome == "M12732.1" ~ "HPV-33",
			chromosome == "NC001357.1" ~ "HPV-18",
			chromosome == "NC001526.4" ~ "HPV-16",
			chromosome == "X74477.1" ~ "HPV-35"
	)) %>%
	dplyr::left_join(smry_ft %>%
			 dplyr::select(sample_name, Is_ctDNA), by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = ifelse(is.na(Is_ctDNA), "?", Is_ctDNA)) %>%
	dplyr::filter(Is_ctDNA != "?") %>%
	dplyr::mutate(mean_af = ifelse(mean_af == 0, .000001, mean_af)) %>%
	dplyr::mutate(median_af = ifelse(median_af == 0, .000001, median_af)) %>%
	dplyr::mutate(max_af = ifelse(max_af == 0, .000001, max_af)) %>%
	dplyr::mutate(start_end = paste0(start, ":", end)) %>%
	dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve"), ordered = TRUE)) %>%
	dplyr::filter(chromosome == "HPV-16" & hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::mutate(Is_ctDNA = case_when(
			hpv_type_wes_wgs == "HPV-16" & Is_ctDNA == "+ve" ~ "+ve",
			TRUE ~ "-ve"
	)) %>%
	dplyr::mutate(aligned_reads = log10(READ_PAIRS + 1)) %>%
	dplyr::mutate(insert_size = log10(MEAN_INSERT_SIZE + 1))

training_ = smry_ %>%
	    reshape2::dcast(formula = sample_name + Is_ctDNA ~ gene_name + start_end, value.var = "aligned_reads", fill = 0) %>%
	    dplyr::rename(`E1_1:1950_READ_COUNT` = `E1_1:1950`,
		          `E2_1891:2989_READ_COUNT` = `E2_1891:2989`,
		          `E5_2985:3237_READ_COUNT` = `E5_2985:3237`,
		          `E6_7124:7601_READ_COUNT` = `E6_7124:7601`,
		          `E7_7603:7900_READ_COUNT` = `E7_7603:7900`,
		          `L1_4774:6292_READ_COUNT` = `L1_4774:6292`) %>%
	    dplyr::left_join(smry_ %>%
			     reshape2::dcast(formula = sample_name + Is_ctDNA ~ gene_name + start_end, value.var = "insert_size", fill = 0) %>%
			     dplyr::rename(`E1_1:1950_INSERT_SIZE` = `E1_1:1950`,
				           `E2_1891:2989_INSERT_SIZE` = `E2_1891:2989`,
				           `E5_2985:3237_INSERT_SIZE` = `E5_2985:3237`,
				           `E6_7124:7601_INSERT_SIZE` = `E6_7124:7601`,
				           `E7_7603:7900_INSERT_SIZE` = `E7_7603:7900`,
				           `L1_4774:6292_INSERT_SIZE` = `L1_4774:6292`),
			     by = c("sample_name", "Is_ctDNA")) %>%
	    tibble::column_to_rownames(var = "sample_name") %>%
	    dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve"), ordered = TRUE))

'lda_cv' <- function(data, n = 100, f = 5) {
	mce = list()
	for (j  in 1:n) {
		folds = split(data, cut(sample(1:nrow(data)),f))
		mce[[j]] = foreach (i=1:length(folds)) %dopar% {
				test = ldply(folds[i], data.frame, .id = NULL)
				train = ldply(folds[-i], data.frame, .id = NULL)
				tmp.fit = MASS::lda(Is_ctDNA ~ ., data = train)
				tmp.prd = predict(tmp.fit, newdata = test)
				conf = table(test$Is_ctDNA, tmp.prd$class)
				1-sum(diag(conf))/sum(conf)
		}
		mce[[j]] = unlist(mce[[j]])
	}
	mce = 100*unlist(lapply(mce, mean))
	return(invisible(mce))
}

'qda_cv' <- function(data, n = 100, f = 5) {
	mce = list()
	for (j  in 1:n) {
		folds = split(data, cut(sample(1:nrow(data)),f))
		mce[[j]] = foreach (i=1:length(folds)) %dopar% {
			test = ldply(folds[i], data.frame, .id = NULL)
			train = ldply(folds[-i], data.frame, .id = NULL)
			tmp.fit = MASS::qda(Is_ctDNA ~ ., data = train)
			tmp.prd = predict(tmp.fit, newdata = test)
			conf = table(test$Is_ctDNA, tmp.prd$class)
			1-sum(diag(conf))/sum(conf)
		}
		mce[[j]] = unlist(mce[[j]])
	}
	mce = 100*unlist(lapply(mce, mean))
	return(invisible(mce))
}

'svm_cv' <- function(data, n = 100, f = 5) {
	mce = list()
	for (j  in 1:n) {
		folds = split(data, cut(sample(1:nrow(data)),f))
		mce[[j]] = foreach (i=1:length(folds)) %dopar% {
			test = ldply(folds[i], data.frame, .id = NULL)
			train = ldply(folds[-i], data.frame, .id = NULL)
			tmp.fit = e1071::svm(Is_ctDNA ~ ., data = train)
			tmp.prd = predict(tmp.fit, newdata = test)
			conf = table(test$Is_ctDNA, tmp.prd)
			1-sum(diag(conf))/sum(conf)
		}
		mce[[j]] = unlist(mce[[j]])
	}
	mce = 100*unlist(lapply(mce, mean))
	return(invisible(mce))
}

'nnet_cv' <- function(data, n = 100, f = 5) {
	mce = list()
	for (j  in 1:n) {
		folds = split(data, cut(sample(1:nrow(data)),f))
		mce[[j]] = foreach (i=1:length(folds)) %dopar% {
			test = ldply(folds[i], data.frame, .id = NULL)
			train = ldply(folds[-i], data.frame, .id = NULL)
			tmp.fit = nnet::nnet(Is_ctDNA ~ ., data = train, size = 10, trace = FALSE)
			tmp.prd = predict(tmp.fit, newdata = test, type = "class")
			conf = table(test$Is_ctDNA, tmp.prd)
			1-sum(diag(conf))/sum(conf)
		}
		mce[[j]] = unlist(mce[[j]])
	}
	mce = 100*unlist(lapply(mce, mean))
	return(invisible(mce))
}

'nb_cv' <- function(data, n = 100, f = 5) {
	mce = list()
	for (j  in 1:n) {
		folds = split(data, cut(sample(1:nrow(data)),f))
		mce[[j]] = foreach (i=1:length(folds)) %dopar% {
			test = ldply(folds[i], data.frame, .id = NULL)
			train = ldply(folds[-i], data.frame, .id = NULL)
			tmp.fit = e1071::naiveBayes(Is_ctDNA ~ ., data = train)
			tmp.prd = predict(tmp.fit, newdata = test, type = "class")
			conf = table(test$Is_ctDNA, tmp.prd)
			1-sum(diag(conf))/sum(conf)
		}
		mce[[j]] = unlist(mce[[j]])
	}
	mce = 100*unlist(lapply(mce, mean))
	return(invisible(mce))
}

'rpart_cv' <- function(data, n = 100, f = 5) {
	paramx = rpart::rpart.control(minsplit = 5,
			      minbucket = 2,
			      maxdepth = 5,
			      maxcompete = 1000,
			      maxsurrogate = 100,
			      usesurrogate = 2,
			      cp = -Inf,
			      xval = 100)
	
	mce = list()
	for (j  in 1:n) {
		folds = split(data, cut(sample(1:nrow(data)),f))
		mce[[j]] = foreach (i=1:length(folds)) %dopar% {
			test = ldply(folds[i], data.frame, .id = NULL)
			train = ldply(folds[-i], data.frame, .id = NULL)
			tmp.fit = rpart::rpart(Is_ctDNA ~ ., data = train, method = "class", control = paramx)
			tmp.prd = predict(tmp.fit, newdata = test, type = "class")
			conf = table(test$Is_ctDNA, tmp.prd)
			1-sum(diag(conf))/sum(conf)
		}
		mce[[j]] = unlist(mce[[j]])
	}
	mce = 100*unlist(lapply(mce, mean))
	return(invisible(mce))
}

x1 = lda_cv(data = training_, n = 50, f = 5)
x2 = qda_cv(data = training_, n = 50, f = 5)
x3 = svm_cv(data = training_, n = 50, f = 5)
x4 = nnet_cv(data = training_, n = 50, f = 5)
x5 = nb_cv(data = training_, n = 50, f = 5)
x6 = rpart_cv(data = training_, n = 50, f = 5)

plot_ = dplyr::tibble(
		`Linear\nDiscriminant\nAnalysis` = x1,
		`Support Vector\nMachines` = x3,
		`Quadratic\nDiscriminant\nAnalysis` = x2,
		`Neural\nNetwork` = x4,
		`Decision\ntree` = x6
	) %>%
	reshape2::melt(variable.name = "classifier", value.name = "mce") %>%
	dplyr::group_by(classifier) %>%
	dplyr::summarize(q1 = min(mce),
			 q2 = mean(mce),
			 q3 = ifelse(max(mce)>10, 10, max(mce))) %>%
	ggplot(aes(x = classifier, ymin = q1+1E-1, y = q2+1E-1, ymax=q3+1E-1)) +
	geom_pointrange(stat = "identity", shape = 21, color = "black", fill = "black", size = .35) +
	scale_x_discrete() +
	scale_y_log10(limits = c(.9E-1, 25),
		      breaks = c(.1, seq(from = .2, to = .8, by = .1), .9, seq(from = 1.9, to = 8.9, by = 1), 9.9),
		      labels = c(0, rep("", length(seq(from = .2, to = .8, by = .1))), 1, rep("", length(seq(from = 1.9, to = 8.9, by = 1))), 10)) +
	xlab("") +
	ylab("5-Fold Cross-validation\nMisclassification Error (%)") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
	      axis.text.y = element_text(size = 12)) +
	guides(color = FALSE)

pdf(file = "../res/Number_Read_Pairs_by_Insert_Size_Compare_Classifiers_Training_Set_HPV-16_Positive.pdf", width = 3.25, height = 3.25)
print(plot_)
dev.off()
