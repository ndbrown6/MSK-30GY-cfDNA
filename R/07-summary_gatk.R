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

posterior_smry = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		 readr::type_convert()

gatk_smry = readr::read_tsv(file = url_gatk_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	    readr::type_convert() %>%
	    dplyr::mutate(DP = T_AD) %>%
	    dplyr::rowwise() %>%
	    dplyr::mutate(AD = as.numeric(lapply(AS_D, str_split, split = ",", n = 2)[[1]][[1]])) %>%
	    dplyr::ungroup() %>%
	    dplyr::mutate(AF = AD/DP) %>%
	    dplyr::left_join(tibble::enframe(target_contigs, name = "HPV", value = "CHROM"), by = "CHROM") %>%
	    dplyr::mutate(CHROM = HPV) %>%
	    dplyr::select(-HPV) %>%
	    dplyr::rename(sample_uuid = sample_name)

plot_ = gatk_smry %>%
	dplyr::left_join(posterior_smry, by = "sample_uuid") %>%
	dplyr::filter(!is.na(sample_name)) %>%
	dplyr::select(CHROM, AD, DP, FS, MQ, QD, SOR, GQ, PL) %>%
	dplyr::mutate(FS = FS + .001) %>%
	reshape2::melt(id.vars = "CHROM", variable.name = "metric", value.name = "score") %>%
	dplyr::mutate(metric = case_when(
		metric == "AD" ~ "Allele depth (AD)",
		metric == "DP" ~ "Total depth (DP)",
		metric == "FS" ~ "P-value strand bias (FS)",
		metric == "MQ" ~ "Mapping quality (MQ)",
		metric == "QD" ~ "Variant confidence/quality by depth (QD)",
		metric == "SOR" ~ "Odds ratio strand bias (SOR)",
		metric == "GQ" ~ "Genotype quality (GQ)",
		metric == "PL" ~ "Genotype likelihoods (PL)"
	)) %>%
	dplyr::mutate(metric = factor(metric, levels = c("Allele depth (AD)", "Total depth (DP)", "Mapping quality (MQ)",
							 "Genotype quality (GQ)", "Genotype likelihoods (PL)", "Variant confidence/quality by depth (QD)",
							 "P-value strand bias (FS)", "Odds ratio strand bias (SOR)"), ordered = TRUE)) %>%
	ggplot(aes(x = CHROM, y = score, color = CHROM)) +
	geom_violin(stat = "ydensity", scale = "area", trim = FALSE, draw_quantiles = .5, fill = "white", size = .75, alpha = 1) +
	scale_color_manual(values = c("HPV-16" = "#e41a1c", "HPV-18" = "#377eb8", "HPV-31" = "#377eb8", "HPV-33" = "#377eb8", "HPV-35" = "#377eb8")) +
	scale_x_discrete() +
	scale_y_log10(labels = scientific_10) +
	xlab("HPV") +
	ylab("Score") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(color = FALSE) +
	facet_wrap(~metric, scales = "free_y", ncol = 2)

pdf(file = "../res/GATK_summary_metrics_by_HPV-16_Positive_All.pdf", width = 9, height = 7)
print(plot_)
dev.off()

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
	     dplyr::mutate(Is_ctDNA = "+ve") %>%
	     dplyr::select(sample_uuid = sample_name, Is_ctDNA, `MRD-Landmark_Result`)

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
	     dplyr::mutate(Is_ctDNA = "-ve") %>%
	     dplyr::select(sample_uuid = sample_name, Is_ctDNA, `MRD-Landmark_Result`)

plot_ = gatk_smry %>%
	dplyr::left_join(smry_t_pos %>%
			 dplyr::bind_rows(smry_t_neg), by = "sample_uuid") %>%
	dplyr::filter(!is.na(Is_ctDNA)) %>%
	dplyr::select(CHROM, Is_ctDNA, AD, DP, FS, MQ, QD, SOR, GQ, PL) %>%
	dplyr::mutate(FS = FS + .001) %>%
	reshape2::melt(id.vars = c("CHROM", "Is_ctDNA"), variable.name = "metric", value.name = "score") %>%
	dplyr::mutate(metric = case_when(
		metric == "AD" ~ "Allele depth (AD)",
		metric == "DP" ~ "Total depth (DP)",
		metric == "FS" ~ "P-value strand bias (FS)",
		metric == "MQ" ~ "Mapping quality (MQ)",
		metric == "QD" ~ "Variant confidence/quality by depth (QD)",
		metric == "SOR" ~ "Odds ratio strand bias (SOR)",
		metric == "GQ" ~ "Genotype quality (GQ)",
		metric == "PL" ~ "Genotype likelihoods (PL)"
	)) %>%
	dplyr::mutate(metric = factor(metric, levels = c("Allele depth (AD)", "Total depth (DP)", "Mapping quality (MQ)",
							 "Genotype quality (GQ)", "Genotype likelihoods (PL)", "Variant confidence/quality by depth (QD)",
							 "P-value strand bias (FS)", "Odds ratio strand bias (SOR)"), ordered = TRUE)) %>%
	ggplot(aes(x = CHROM, y = score, color = CHROM, fill = Is_ctDNA)) +
	geom_violin(stat = "ydensity", scale = "area", trim = FALSE, draw_quantiles = .5, size = .75) +
	scale_color_manual(values = c("HPV-16" = "#e41a1c", "HPV-18" = "#377eb8", "HPV-31" = "#377eb8", "HPV-33" = "#377eb8", "HPV-35" = "#377eb8")) +
	scale_fill_manual(values = c("+ve" = "#cccccc", "-ve" = "#f7f7f7")) +
	scale_x_discrete() +
	scale_y_log10(labels = scientific_10) +
	xlab("HPV") +
	ylab("Score") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(color = FALSE,
	       fill = FALSE) +
	facet_wrap(~metric, scales = "free_y", ncol = 2)

pdf(file = "../res/GATK_summary_metrics_by_HPV-16_Positive_Training_Set.pdf", width = 9, height = 7)
print(plot_)
dev.off()