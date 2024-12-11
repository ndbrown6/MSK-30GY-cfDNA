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

#==================================================
# Relative Mean AF
#==================================================
smry_ = readr::read_tsv(file = url_idx_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
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
	)) %>%
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(mean_af = mean(mean_af+1E-5, na.rm = TRUE)) %>%
	reshape2::dcast(formula = patient_id_mskcc ~ timepoint_weeks_since_start_of_RT,
			fun.aggregate = function(x) { mean(x, na.rm=TRUE) }, fill = NaN, value.var = "mean_af") %>%
	dplyr::select(patient_id_mskcc, `Pre-treatment`, wk1, wk2, wk3) %>%
	tidyr::drop_na(`Pre-treatment`, wk1, wk2, wk3) %>%
	dplyr::mutate(wk1 = log2(wk1/`Pre-treatment`),
		      wk2 = log2(wk2/`Pre-treatment`),
		      wk3 = log2(wk3/`Pre-treatment`),
		      `Pre-treatment` = 0)

set.seed(3)

tsne_ = smry_ %>%
        dplyr::select(wk1, wk2, wk3) %>%
        Rtsne::Rtsne(perplexity = 10, theta = 0, pca = TRUE, pca_center = FALSE, pca_scale = FALSE, normalize = FALSE, exaggeration_factor = 5)

pca_ = smry_ %>%
       dplyr::select(wk1, wk2, wk3) %>%
       stats::prcomp(center = FALSE, scale. = FALSE)

hcl_ = smry_ %>%
       dplyr::select(wk1, wk2, wk3) %>%
       dist(method = "euclidean") %>%
       hclust(method = "ward.D") %>%
       cutree(k = 4)

smry_ = smry_ %>%
	dplyr::bind_cols(tsne_$Y %>%
			 dplyr::as_tibble() %>%
			 dplyr::select(tSNE1 = V1, tSNE2 = V2)) %>%
	dplyr::bind_cols(pca_$x %>%
			 dplyr::as_tibble() %>%
			 dplyr::select(PC1, PC2)) %>%
	dplyr::mutate(cluster_id_af = as.character(hcl_)) %>%
	dplyr::mutate(cluster_id_af = case_when(
		cluster_id_af == "4" ~ "Fast",
		cluster_id_af == "3" ~ "Slow",
		cluster_id_af == "2" ~ "No clearance",
		cluster_id_af == "1" ~ "Intermediate"
	)) %>%
	dplyr::mutate(cluster_id_af = factor(cluster_id_af, levels = c("Fast", "Slow", "Intermediate", "No clearance"), ordered = TRUE))

set.seed(3)

tsne_ = clinical %>%
	dplyr::select(MRI_rawdata_wk0, MRI_rawdata_wk1, MRI_rawdata_wk2, MRI_rawdata_wk3) %>%
	tidyr::drop_na() %>%
	Rtsne::Rtsne(perplexity = 10, theta = 0, pca = TRUE, pca_center = TRUE, pca_scale = TRUE, normalize = TRUE, exaggeration_factor = 5)

pca_ = clinical %>%
       dplyr::select(MRI_rawdata_wk0, MRI_rawdata_wk1, MRI_rawdata_wk2, MRI_rawdata_wk3) %>%
       tidyr::drop_na() %>%
       stats::prcomp(center = TRUE, scale. = TRUE)

hcl_ = clinical %>%
       dplyr::select(MRI_rawdata_wk0, MRI_rawdata_wk1, MRI_rawdata_wk2, MRI_rawdata_wk3) %>%
       tidyr::drop_na() %>%
       dist(method = "euclidean") %>%
       hclust(method = "ward.D") %>%
       cutree(k = 4)

plot_ = clinical %>%
	tidyr::drop_na(MRI_rawdata_wk0, MRI_rawdata_wk1, MRI_rawdata_wk2, MRI_rawdata_wk3) %>%
	dplyr::bind_cols(tsne_$Y %>%
			 dplyr::as_tibble() %>%
			 dplyr::select(tSNE1 = V1, tSNE2 = V2)) %>%
	dplyr::bind_cols(pca_$x %>%
			 dplyr::as_tibble() %>%
			 dplyr::select(PC1, PC2)) %>%
	dplyr::mutate(cluster_id_vol = as.character(hcl_)) %>%
	ggplot(aes(x = tSNE1, y = tSNE2, color = cluster_id_vol)) +
	geom_point(stat = "identity")

plot_ = clinical %>%
	dplyr::select(patient_id_mskcc, MRI_rawdata_wk0, MRI_rawdata_wk1, MRI_rawdata_wk2, MRI_rawdata_wk3, MRI_rawdata_wk4) %>%
	tidyr::drop_na(MRI_rawdata_wk0, MRI_rawdata_wk1, MRI_rawdata_wk2, MRI_rawdata_wk3) %>%
	dplyr::mutate(cluster_id_vol = as.character(hcl_)) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "cluster_id_vol")) %>%
	ggplot(aes(x = variable, y = value, color = cluster_id_vol, shape = cluster_id_vol)) +
	geom_boxplot(stat = "boxplot", outlier.shape = NA, color = "black") +
	geom_jitter(stat = "identity", width = .1, height = 0, fill = "white", size = 3, alpha = .95) +
	scale_shape_manual(values = c("4" = 21, "1" = 22, "3" = 23, "2" = 24)) +
	scale_color_manual(values = c("4" = "#1b9e77", "1" = "#d95f02", "3" = "#7570b3", "2" = "#e7298a")) +
	scale_x_discrete(breaks = c("MRI_rawdata_wk0", "MRI_rawdata_wk1", "MRI_rawdata_wk2", "MRI_rawdata_wk3", "MRI_rawdata_wk4"),
			 labels = c(0, 1, 2, 3, 4)) +
	scale_y_continuous(limits = c(0, 70000),
			   labels = scientific_10) +
	xlab("") +
	ylab("MRI Volume") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 9.5),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	guides(color = FALSE, shape = FALSE) +
	facet_wrap(~cluster_id_vol, nrow = 1)

tt_ = smry_ %>%
      dplyr::left_join(clinical %>%
		       tidyr::drop_na(MRI_rawdata_wk0, MRI_rawdata_wk1, MRI_rawdata_wk2, MRI_rawdata_wk3) %>%
		       dplyr::mutate(cluster_id_vol = as.character(hcl_)),
		       by = "patient_id_mskcc") %>%
      dplyr::group_by(cluster_id_vol, cluster_id_af) %>%
      dplyr::summarize(n = n()) %>%
      dplyr::ungroup() %>%
      tidyr::drop_na() %>%
      reshape2::dcast(cluster_id_vol ~ cluster_id_af, fill = 0)

tt_ = clinical %>%
      tidyr::drop_na(MRI_rawdata_wk0, MRI_rawdata_wk1, MRI_rawdata_wk2, MRI_rawdata_wk3) %>%
      dplyr::mutate(cluster_id_vol = as.character(hcl_)) %>%
      dplyr::group_by(cluster_id_vol, simplified_hypoxia_group) %>%
      dplyr::summarize(n = n()) %>%
      dplyr::ungroup() %>%
      tidyr::drop_na() %>%
      reshape2::dcast(cluster_id_vol ~ simplified_hypoxia_group, fill = 0)


