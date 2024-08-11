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

#########################################################
# (+) Samples with mean AF > 5% from MRD assay
# (+) Samples with HPV subtype in assay
# (+) Samples with +ve MRD assay
#########################################################
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

#########################################################
# (-) Samples with mean AF < 0.01% from MRD assay
# (-) Samples with max AF < 0.1% from MRD assay
# (-) Patients with no nodal dissection ≥ 2 years
# (-) No duplicate patients
# (+) Samples with -ve MRD assay
#########################################################
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

#########################################################
# Test sets HPV-16 -ve
#########################################################
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
	dplyr::mutate(mean_af = ifelse(mean_af == 0, .000001, mean_af)) %>%
	dplyr::mutate(median_af = ifelse(median_af == 0, .000001, median_af)) %>%
	dplyr::mutate(max_af = ifelse(max_af == 0, .000001, max_af)) %>%
	dplyr::mutate(start_end = paste0(start, ":", end)) %>%
	dplyr::filter(chromosome == "HPV-16" & hpv_type_wes_wgs != "HPV-16" & hpv_type_wes_wgs != "Unknown") %>%
	dplyr::mutate(aligned_reads = log10(READ_PAIRS + 1)) %>%
	dplyr::mutate(insert_size = log10(MEAN_INSERT_SIZE + 1)) %>%
	dplyr::filter(!(sample_name %in% smry_ft$sample_name))

test_ = smry_ %>%
	dplyr::left_join(mrd_smry %>%
			 dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			 dplyr::select(sample_name, `MRD-Landmark_Result`) %>%
			 dplyr::mutate(Is_ctDNA = case_when(
				 `MRD-Landmark_Result` == "ABSENT" ~ "-ve",
				 `MRD-Landmark_Result` == "PRESENT" ~ "+ve"
			 )), by = "sample_name") %>%
     	reshape2::dcast(formula = sample_name + Is_ctDNA ~ gene_name + start_end, value.var = "aligned_reads", fill = 0) %>%
	dplyr::rename(`E1_1:1950_READ_COUNT` = `E1_1:1950`,
		      `E2_1891:2989_READ_COUNT` = `E2_1891:2989`,
		      `E5_2985:3237_READ_COUNT` = `E5_2985:3237`,
		      `E6_7124:7601_READ_COUNT` = `E6_7124:7601`,
		      `E7_7603:7900_READ_COUNT` = `E7_7603:7900`,
		      `L1_4774:6292_READ_COUNT` = `L1_4774:6292`) %>%
	dplyr::left_join(smry_ %>%
			 dplyr::left_join(mrd_smry %>%
			 dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			 dplyr::select(sample_name, `MRD-Landmark_Result`) %>%
			 dplyr::mutate(Is_ctDNA = case_when(
				 `MRD-Landmark_Result` == "ABSENT" ~ "-ve",
				 `MRD-Landmark_Result` == "PRESENT" ~ "+ve"
			 )), by = "sample_name") %>%
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

fit_ = MASS::lda(Is_ctDNA ~ ., data = training_)
res_ = predict(object = fit_, newdata = test_)
	
plot_ = do.call(cbind, res_) %>%
	dplyr::as_tibble() %>%
	dplyr::mutate(set = "test") %>%
	dplyr::bind_rows(
		do.call(cbind, predict(object = fit_, newdata = training_)) %>%
		dplyr::as_tibble() %>%
		dplyr::mutate(set = "training")
	) %>%
	dplyr::mutate(class = factor(class, levels = c(1, 2), ordered = TRUE)) %>%
	dplyr::mutate(set = factor(set, levels = c("training", "test"), ordered = TRUE)) %>%
	ggplot(aes(x = class:set, y = log10(`+ve`))) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, fill = "white", show.legend = FALSE) +
	geom_jitter(data = do.call(cbind, res_) %>%
		    dplyr::as_tibble() %>%
		    dplyr::mutate(set = "test") %>%
		    dplyr::bind_rows(
			    do.call(cbind, predict(object = fit_, newdata = training_)) %>%
			    dplyr::as_tibble() %>%
			    dplyr::mutate(set = "training")
		    ) %>%
		    dplyr::mutate(class = factor(class, levels = c(1, 2), ordered = TRUE)) %>%
		    dplyr::mutate(set = factor(set, levels = c("training", "test"), ordered = TRUE)),
		    mapping = aes(x = class:set, y = log10(`+ve`), fill = set),
		    stat = "identity", width = .1, height = 0, shape = 21, alpha = .85, size = 3) +
	scale_fill_manual(values = c("training" = "#377eb8",
				     "test" = "#e41a1c")) +
	scale_x_discrete(breaks = c("1:training", "2:training", "2:test"),
			 labels = c("+ve", "-ve", "Unknown")) +
	scale_y_continuous(limits = c(-200, 20),
			   breaks = c(-200, -175, -150, -125, -100, -75, -50, -25, 0),
			   labels = c(-200, "", -150, "", -100, "", -50, "", 0)) +
	xlab("") +
	ylab(expression("Posterior Probability ("~-log[10]~")")) +
	geom_signif(stat = "signif",
		    comparisons = list(c("1:training", "2:training")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    y_position = 15,
		    tip_length = 0.01) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20), size = 14),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = FALSE)

pdf(file = "../res/Posterior_Probability_Test_HPV-16_Negative.pdf", width = 2.95, height = 3.25)
print(plot_)
dev.off()

plot_ = fit_$scaling %>%
	as.data.frame() %>%
	tibble::rownames_to_column("variable") %>%
	dplyr::as_tibble() %>%
	dplyr::rename(coefficient  = "LD1") %>%
	dplyr::mutate(is_insert_size = ifelse(grepl("INSERT_SIZE", variable, fixed = TRUE), "Mean Insert Size", "Aligned Read Pairs")) %>%
	dplyr::mutate(variable = unlist(lapply(variable, function(x) { strsplit(x, "_")[[1]][1] } ))) %>%
	dplyr::mutate(variable = gsub(pattern = "`", replacement = "", x = variable)) %>%
	dplyr::mutate(variable = factor(variable, levels = rev(c("E1", "E2", "E5", "L1", "E6", "E7")), ordered = TRUE)) %>%
	ggplot(aes(x = variable, ymin = 0, ymax = coefficient, y = coefficient, shape = is_insert_size)) +
	geom_linerange(stat = "identity", size = .5) +
	geom_hline(yintercept = 0, size = .5) +
	geom_point(stat = "identity", fill = "white", size = 1.5) +
	xlab("") +
	ylab("") +
	scale_shape_manual(values = c(21, 22)) +
	scale_y_continuous(limits = c(-40, 20),
			   breaks = c(-40, -30, -20, -10, 0, 10, 20),
			   labels = c(-40, "", -20, "", 0, "", 20)) +
	coord_flip() +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      panel.spacing = unit(1, "lines"),
	      strip.background = element_rect(colour="white", fill="white")) +
	facet_wrap(~is_insert_size, ncol = 1, scales = "free") +
	guides(shape = FALSE)

pdf(file = "../res/Coefficients_Linear_Discriminant_Analysis.pdf", width = 2, height = 3.25)
print(plot_)
dev.off()

fit_ = klaR::greedy.wilks(formula = Is_ctDNA ~ ., data = data.frame(training_), niveau = Inf)

plot_ = fit_ %>%
	.[["results"]] %>%
	dplyr::as_tibble() %>%
	dplyr::rename(variable  = vars) %>%
	dplyr::mutate(is_insert_size = ifelse(grepl("INSERT_SIZE", variable, fixed = TRUE), "Mean Insert Size", "Aligned Read Pairs")) %>%
	dplyr::mutate(variable = unlist(lapply(variable, function(x) { strsplit(x, "_")[[1]][1] } ))) %>%
	dplyr::mutate(variable = gsub(pattern = "`", replacement = "", x = variable)) %>%
	dplyr::mutate(index = 1:nrow(.)) %>%
	ggplot(aes(x = index, y = `Wilks.lambda`)) +
	geom_line(stat = "identity", size = .5) +
	geom_point(data = fit_ %>%
		   	  .[["results"]] %>%
		   	  dplyr::as_tibble() %>%
		   	  dplyr::rename(variable  = vars) %>%
		   	  dplyr::mutate(is_insert_size = ifelse(grepl("INSERT_SIZE", variable, fixed = TRUE), "Mean Insert Size", "Aligned Read Pairs")) %>%
		   	  dplyr::mutate(variable = unlist(lapply(variable, function(x) { strsplit(x, "_")[[1]][1] } ))) %>%
		   	  dplyr::mutate(variable = gsub(pattern = "`", replacement = "", x = variable)) %>%
		   	  dplyr::mutate(index = 1:nrow(.)),
		   mapping = aes(x = index, y = `Wilks.lambda`, shape = is_insert_size), stat = "identity", fill = "white", size = 1.5, inherit.aes = FALSE) +
	xlab("") +
	ylab(expression("Wilks'"~lambda)) +
	scale_shape_manual(values = c(21, 22)) +
	scale_x_continuous(limits = c(1, 12),
			   breaks = seq(from = 1, to = 12, by = 1),
			   labels =  fit_ %>%
				     .[["results"]] %>%
			   	     dplyr::as_tibble() %>%
			   	     dplyr::rename(variable  = vars) %>%
			   	     dplyr::mutate(is_insert_size = ifelse(grepl("INSERT_SIZE", variable, fixed = TRUE), "Mean Insert Size", "Aligned Read Pairs")) %>%
			   	     dplyr::mutate(variable = unlist(lapply(variable, function(x) { strsplit(x, "_")[[1]][1] } ))) %>%
			   	     dplyr::mutate(variable = gsub(pattern = "`", replacement = "", x = variable)) %>%
			   	     .[["variable"]]) +
	scale_y_continuous(breaks = c(0.006, 0.0065, 0.007, 0.0075, 0.008, 0.0085),
			   labels = c(".006", "", ".007", "", ".008", "")) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 7)) +
	guides(shape = guide_legend(title = ""))

pdf(file = "../res/Wilks_Lambda_Linear_Discriminant_Analysis.pdf", width = 4.25, height = 2.25)
print(plot_)
dev.off()
	