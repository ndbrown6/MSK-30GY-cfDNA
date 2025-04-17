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
		      abs_MRI_wk0 = MRI_rawdata_wk0,
		      abs_MRI_wk1 = MRI_rawdata_wk1,
		      abs_MRI_wk2 = MRI_rawdata_wk2,
		      abs_MRI_wk3 = MRI_rawdata_wk3) %>%
	dplyr::mutate(composite_end_point = case_when(
		composite_end_point ~ 1,
		!composite_end_point ~ 0
	)) %>%
	tidyr::drop_na(composite_end_point)

data = data_ %>%
       dplyr::select(composite_end_point, contains("abs")) %>%
       dplyr::mutate(composite_end_point = factor(composite_end_point)) %>%
       tidyr::drop_na()

set.seed(1234567)
rf = randomForest(formula = composite_end_point ~ ., data = data,
		  ntree = 1000, nodesize = 1, maxnodes = 10,
		  importance=TRUE, proximity=TRUE)
vip = importance(rf) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("variables") %>%
      dplyr::mutate(variables = case_when(
	      variables == "abs_MRI_wk1" ~ "abs_MRI_wk3",
	      variables == "abs_MRI_wk3" ~ "abs_MRI_wk1",
	      TRUE ~ variables
	      
      ))

null_vip = foreach (i=1:100) %dopar% {
	data = data %>%
	       dplyr::mutate(composite_end_point = sample(x = composite_end_point, size = nrow(.), replace = TRUE))
	rf = randomForest(formula = composite_end_point ~ ., data = data,
		  ntree = 1000, nodesize = 1, maxnodes = 10,
		  importance=TRUE, proximity=TRUE)
	null_vip = importance(rf) %>%
		   as.data.frame() %>%
		   tibble::rownames_to_column("variables")
	return(invisible(null_vip))
}

plot_ = do.call(rbind, null_vip) %>%
	dplyr::mutate(variables = gsub("abs_AF_", "ctDNA Fraction (", variables)) %>%
	dplyr::mutate(variables = gsub("abs_HPV_", "cfDNA HPV Reads (", variables)) %>%
	dplyr::mutate(variables = gsub("abs_MRI_", "MRI Volume (", variables)) %>%
	dplyr::mutate(variables = paste0(variables, ")")) %>%
	dplyr::mutate(variables = gsub(pattern = "wk", replacement = "Week ", x = variables)) %>%
	dplyr::mutate(week = case_when(
		grepl("Week 0", variables) ~ 0,
		grepl("Week 1", variables) ~ 1,
		grepl("Week 2", variables) ~ 2,
		grepl("Week 3", variables) ~ 3
	)) %>%
	dplyr::mutate(category = case_when(
		grepl("ctDNA", variables) ~ 1,
		grepl("cfDNA", variables) ~ 2,
		grepl("MRI", variables) ~ 3
	)) %>%
	dplyr::arrange(desc(week), desc(category)) %>%
	dplyr::mutate(variables = factor(variables, levels = unique(variables), ordered = TRUE)) %>%
	ggplot(aes(x = MeanDecreaseAccuracy, y = variables, fill = factor(week), color = factor(week))) +
	geom_density_ridges(stat = "density_ridges", alpha = .45, bandwidth = .25, scale = 1.25) +
	geom_segment(data = vip %>%
		     	    dplyr::mutate(variables = gsub("abs_AF_", "ctDNA Fraction (", variables)) %>%
		   	    dplyr::mutate(variables = gsub("abs_HPV_", "cfDNA HPV Reads (", variables)) %>%
		   	    dplyr::mutate(variables = gsub("abs_MRI_", "MRI Volume (", variables)) %>%
		   	    dplyr::mutate(variables = paste0(variables, ")")) %>%
		   	    dplyr::mutate(variables = gsub(pattern = "wk", replacement = "Week ", x = variables)) %>%
		     	    dplyr::mutate(week = case_when(
				    grepl("Week 0", variables) ~ 0,
				    grepl("Week 1", variables) ~ 1,
				    grepl("Week 2", variables) ~ 2,
				    grepl("Week 3", variables) ~ 3
			    )) %>%
		     	    dplyr::mutate(category = case_when(
				    grepl("ctDNA", variables) ~ 1,
				    grepl("cfDNA", variables) ~ 2,
				    grepl("MRI", variables) ~ 3
			    )) %>%
		     	    dplyr::arrange(desc(week), desc(category)) %>%
		     	    dplyr::mutate(variables = factor(variables, levels = unique(variables), ordered = TRUE)),
		   aes(x = MeanDecreaseAccuracy, y = 1:12, xend = MeanDecreaseAccuracy, yend = (1:12)+.5),
		   color = "black", size = .5, inherit.aes = FALSE) +
	geom_point(data = vip %>%
		   	  dplyr::mutate(variables = gsub("abs_AF_", "ctDNA Fraction (", variables)) %>%
		   	  dplyr::mutate(variables = gsub("abs_HPV_", "cfDNA HPV Reads (", variables)) %>%
		   	  dplyr::mutate(variables = gsub("abs_MRI_", "MRI Volume (", variables)) %>%
		   	  dplyr::mutate(variables = paste0(variables, ")")) %>%
		   	  dplyr::mutate(variables = gsub(pattern = "wk", replacement = "Week ", x = variables)) %>%
		     	  dplyr::mutate(week = case_when(
				    grepl("Week 0", variables) ~ 0,
				    grepl("Week 1", variables) ~ 1,
				    grepl("Week 2", variables) ~ 2,
				    grepl("Week 3", variables) ~ 3
			  )) %>%
		     	  dplyr::mutate(category = case_when(
				    grepl("ctDNA", variables) ~ 1,
				    grepl("cfDNA", variables) ~ 2,
				    grepl("MRI", variables) ~ 3
			  )) %>%
		     	  dplyr::arrange(desc(week), desc(category)) %>%
		     	  dplyr::mutate(variables = factor(variables, levels = unique(variables), ordered = TRUE)),
		   aes(x = MeanDecreaseAccuracy, y = (1:12)+.5, fill = factor(week)),
		   color = "black", shape = 21, size = 2, inherit.aes = FALSE) +
	scale_fill_brewer(type = "qual", palette = 7) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_x_continuous(limits = c(-3, 5)) +
	xlab("\nVariable Importance\n\n") +
	ylab("") +
	theme_classic() +
	theme(axis.title.x = element_text(size = 14),
	      axis.title.y = element_text(size = 14),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12, hjust = 1),
	      axis.line.y = element_blank(),
	      axis.ticks.y = element_blank()) +
	guides(fill = FALSE, color = FALSE)

pdf(file = "../res/Random_Forest_Mean_Decrease_Accuracy.pdf", width = 5, height = 5.5)
print(plot_)
dev.off()

plot_ = do.call(rbind, null_vip) %>%
	dplyr::mutate(variables = gsub("abs_AF_", "ctDNA Fraction (", variables)) %>%
	dplyr::mutate(variables = gsub("abs_HPV_", "cfDNA HPV Reads (", variables)) %>%
	dplyr::mutate(variables = gsub("abs_MRI_", "MRI Volume (", variables)) %>%
	dplyr::mutate(variables = paste0(variables, ")")) %>%
	dplyr::mutate(variables = gsub(pattern = "wk", replacement = "Week ", x = variables)) %>%
	dplyr::mutate(week = case_when(
		grepl("Week 0", variables) ~ 0,
		grepl("Week 1", variables) ~ 1,
		grepl("Week 2", variables) ~ 2,
		grepl("Week 3", variables) ~ 3
	)) %>%
	dplyr::mutate(category = case_when(
		grepl("ctDNA", variables) ~ 1,
		grepl("cfDNA", variables) ~ 2,
		grepl("MRI", variables) ~ 3
	)) %>%
	dplyr::arrange(desc(week), desc(category)) %>%
	dplyr::mutate(variables = factor(variables, levels = unique(variables), ordered = TRUE)) %>%
	ggplot(aes(x = MeanDecreaseGini, y = variables, fill = factor(week), color = factor(week))) +
	geom_density_ridges(stat = "density_ridges", alpha = .45, bandwidth = .05, scale = 1.25) +
	geom_segment(data = vip %>%
		     	    dplyr::mutate(variables = gsub("abs_AF_", "ctDNA Fraction (", variables)) %>%
		   	    dplyr::mutate(variables = gsub("abs_HPV_", "cfDNA HPV Reads (", variables)) %>%
		   	    dplyr::mutate(variables = gsub("abs_MRI_", "MRI Volume (", variables)) %>%
		   	    dplyr::mutate(variables = paste0(variables, ")")) %>%
		   	    dplyr::mutate(variables = gsub(pattern = "wk", replacement = "Week ", x = variables)) %>%
		     	    dplyr::mutate(week = case_when(
				    grepl("Week 0", variables) ~ 0,
				    grepl("Week 1", variables) ~ 1,
				    grepl("Week 2", variables) ~ 2,
				    grepl("Week 3", variables) ~ 3
			    )) %>%
		     	    dplyr::mutate(category = case_when(
				    grepl("ctDNA", variables) ~ 1,
				    grepl("cfDNA", variables) ~ 2,
				    grepl("MRI", variables) ~ 3
			    )) %>%
		     	    dplyr::arrange(desc(week), desc(category)) %>%
		     	    dplyr::mutate(variables = factor(variables, levels = unique(variables), ordered = TRUE)),
		   aes(x = MeanDecreaseGini, y = 1:12, xend = MeanDecreaseGini, yend = (1:12)+.5),
		   color = "black", size = .5, inherit.aes = FALSE) +
	geom_point(data = vip %>%
		   	  dplyr::mutate(variables = gsub("abs_AF_", "ctDNA Fraction (", variables)) %>%
		   	  dplyr::mutate(variables = gsub("abs_HPV_", "cfDNA HPV Reads (", variables)) %>%
		   	  dplyr::mutate(variables = gsub("abs_MRI_", "MRI Volume (", variables)) %>%
		   	  dplyr::mutate(variables = paste0(variables, ")")) %>%
		   	  dplyr::mutate(variables = gsub(pattern = "wk", replacement = "Week ", x = variables)) %>%
		     	  dplyr::mutate(week = case_when(
				    grepl("Week 0", variables) ~ 0,
				    grepl("Week 1", variables) ~ 1,
				    grepl("Week 2", variables) ~ 2,
				    grepl("Week 3", variables) ~ 3
			  )) %>%
		     	  dplyr::mutate(category = case_when(
				    grepl("ctDNA", variables) ~ 1,
				    grepl("cfDNA", variables) ~ 2,
				    grepl("MRI", variables) ~ 3
			  )) %>%
		     	  dplyr::arrange(desc(week), desc(category)) %>%
		     	  dplyr::mutate(variables = factor(variables, levels = unique(variables), ordered = TRUE)),
		   aes(x = MeanDecreaseGini, y = (1:12)+.5, fill = factor(week)),
		   color = "black", shape = 21, size = 2, inherit.aes = FALSE) +
	scale_fill_brewer(type = "qual", palette = 7) +
	scale_color_brewer(type = "qual", palette = 7) +
	scale_x_continuous() +
	xlab("\nVariable Importance\n\n") +
	ylab("") +
	theme_classic() +
	theme(axis.title.x = element_text(size = 14),
	      axis.title.y = element_text(size = 14),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12, hjust = 1),
	      axis.line.y = element_blank(),
	      axis.ticks.y = element_blank()) +
	guides(fill = FALSE, color = FALSE)

pdf(file = "../res/Random_Forest_Mean_Decrease_Gini.pdf", width = 5, height = 5.5)
print(plot_)
dev.off()