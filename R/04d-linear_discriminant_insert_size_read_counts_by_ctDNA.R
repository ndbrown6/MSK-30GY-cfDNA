#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::mutate(sample_uuid = paste0(sample_id_mskcc, "-", sample_id_invitae))

bed_file = readr::read_tsv(file = url_bed_file, col_names = FALSE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::rename(chromosome = X1,
			 start = X2,
			 end = X3,
			 gene_name = X4) %>%
	   dplyr::select(-X5) %>%
	   dplyr::filter(gene_name != "E1^E4")

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
	     # Mean AF < 0.01%
	     dplyr::filter(mean_af < (0.01/100)) %>%
	     # Max AF < 0.1%
	     dplyr::filter(max_af < (0.1/100)) %>%
	     dplyr::left_join(manifest, by = "sample_name") %>%
	     dplyr::left_join(nodal_dissection_smry, by = "patient_id_mskcc") %>%
	     # No nodal dissection ≥ 2 years
	     dplyr::filter(nd_event == 0 & timepoint_days_since_start_of_RT > 730) %>%
	     # No duplicate patients
	     dplyr::left_join(mutation_smry %>%
	     		      dplyr::filter(FILTER == "PASS") %>%
			      dplyr::group_by(Tumor_Sample_Barcode) %>%
			      dplyr::summarize(mean_af = mean(t_maf),
					       max_af = max(t_maf)) %>%
			      dplyr::ungroup() %>%
			      dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
			      # Mean AF < 0.01%
			      dplyr::filter(mean_af < (0.01/100)) %>%
			      # Max AF < 0.1%
			      dplyr::filter(max_af < (0.1/100)) %>%
			      dplyr::left_join(manifest, by = "sample_name") %>%
			      dplyr::left_join(nodal_dissection_smry, by = "patient_id_mskcc") %>%
			      # No nodal dissection ≥ 2 years
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
# Baseline samples
#########################################################
smry_baseline__mrd_hpv = manifest %>%
	       	     	 dplyr::mutate(hpv_type_wes_wgs = case_when(
				 		is.na(hpv_type_wes_wgs) ~ "Unknown",
			     			TRUE ~ hpv_type_wes_wgs
			 )) %>%
			 dplyr::left_join(mrd_smry %>%
					  dplyr::mutate(sample_uuid = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)),
					  by = "sample_uuid") %>%
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
			 dplyr::mutate(Is_ctDNA = Is_ctDNA_MRD | Is_ctDNA_HPV)

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
	dplyr::filter(chromosome == "HPV-16" & hpv_type_wes_wgs == "HPV-16") %>%
	dplyr::mutate(aligned_reads = log10(READ_PAIRS + 1)) %>%
	dplyr::mutate(insert_size = log10(MEAN_INSERT_SIZE + 1)) %>%
	dplyr::left_join(smry_baseline__mrd_hpv %>%
			 dplyr::mutate(is_present = TRUE), by = "patient_name") %>%
	dplyr::filter(is_present) %>%
	dplyr::filter(timepoint_days_since_start_of_RT<=0)

test_ = smry_ %>%
	reshape2::dcast(formula = sample_name ~ gene_name + start_end, value.var = "aligned_reads", fill = 0) %>%
	dplyr::rename(`E1_1:1950_READ_COUNT` = `E1_1:1950`,
		      `E2_1891:2989_READ_COUNT` = `E2_1891:2989`,
		      `E5_2985:3237_READ_COUNT` = `E5_2985:3237`,
		      `E6_7124:7601_READ_COUNT` = `E6_7124:7601`,
		      `E7_7603:7900_READ_COUNT` = `E7_7603:7900`,
		      `L1_4774:6292_READ_COUNT` = `L1_4774:6292`) %>%
	dplyr::left_join(smry_ %>%
			 reshape2::dcast(formula = sample_name ~ gene_name + start_end, value.var = "insert_size", fill = 0) %>%
			 dplyr::rename(`E1_1:1950_INSERT_SIZE` = `E1_1:1950`,
				       `E2_1891:2989_INSERT_SIZE` = `E2_1891:2989`,
				       `E5_2985:3237_INSERT_SIZE` = `E5_2985:3237`,
				       `E6_7124:7601_INSERT_SIZE` = `E6_7124:7601`,
				       `E7_7603:7900_INSERT_SIZE` = `E7_7603:7900`,
				       `L1_4774:6292_INSERT_SIZE` = `L1_4774:6292`),
			 by = "sample_name") %>%
	tibble::column_to_rownames(var = "sample_name")


aligned_reads = seq(from = 0, to = 6, length = N_TILES)
insert_size = seq(from = 0, to = 3, length = N_TILES)
mat_ = dplyr::tibble(insert_size = rep(insert_size, each = N_TILES),
		     aligned_reads = rep(aligned_reads, times = N_TILES)) %>%
       dplyr::mutate(`E1_1:1950_READ_COUNT` = aligned_reads,
		     `E2_1891:2989_READ_COUNT` = aligned_reads,
		     `E5_2985:3237_READ_COUNT` = aligned_reads,
		     `E6_7124:7601_READ_COUNT` = aligned_reads,
		     `E7_7603:7900_READ_COUNT` = aligned_reads,
		     `L1_4774:6292_READ_COUNT` = aligned_reads,
		     `E1_1:1950_INSERT_SIZE` = insert_size,
		     `E2_1891:2989_INSERT_SIZE` = insert_size,
		     `E5_2985:3237_INSERT_SIZE` = insert_size,
		     `E6_7124:7601_INSERT_SIZE` = insert_size,
		     `E7_7603:7900_INSERT_SIZE` = insert_size,
		     `L1_4774:6292_INSERT_SIZE` = insert_size) %>%
      dplyr::select(-aligned_reads, -insert_size)


observed = simulated = list()

for (i in 1:length(target_names)) {
	fit_ = MASS::lda(Is_ctDNA ~ ., data = training_ %>% dplyr::select(c("Is_ctDNA", paste0(target_names[i], "_READ_COUNT"), paste0(target_names[i], "_INSERT_SIZE"))))
	res_ = predict(object = fit_, newdata = test_ %>% dplyr::select(c(paste0(target_names[i], "_READ_COUNT"), paste0(target_names[i], "_INSERT_SIZE"))))
	sim_ = predict(object = fit_, newdata = mat_ %>% dplyr::select(c(paste0(target_names[i], "_READ_COUNT"), paste0(target_names[i], "_INSERT_SIZE"))))
	
	observed[[i]] = test_ %>%
			dplyr::mutate(Pr = res_$posterior[,"+ve"]) %>%
			dplyr::mutate(Is_ctDNA = res_$class) %>%
			dplyr::select(c("Pr", paste0(target_names[i], "_READ_COUNT"), paste0(target_names[i], "_INSERT_SIZE"))) %>%
			dplyr::rename(read_count = paste0(target_names[i], "_READ_COUNT"),
				      insert_size = paste0(target_names[i], "_INSERT_SIZE")) %>%
			dplyr::mutate(target_name = target_names[i]) %>%
			dplyr::mutate(sample_name = rownames(test_))
	simulated[[i]] = mat_ %>%
			 dplyr::mutate(Pr = sim_$posterior[,"+ve"]) %>%
			 dplyr::select(c("Pr", paste0(target_names[i], "_READ_COUNT"), paste0(target_names[i], "_INSERT_SIZE"))) %>%
			 dplyr::rename(read_count = paste0(target_names[i], "_READ_COUNT"),
				       insert_size = paste0(target_names[i], "_INSERT_SIZE")) %>%
			 dplyr::mutate(target_name = target_names[i])
}

plot_ = do.call(rbind, observed) %>%
	dplyr::rename(sample_uuid = sample_name) %>%
	dplyr::left_join(smry_ %>%
			 dplyr::filter(!duplicated(sample_uuid)) %>%
			 dplyr::select(-insert_size), by = "sample_uuid") %>%
	dplyr::filter(!Is_ctDNA_MRD) %>%
	dplyr::mutate(end = gsub(pattern = "E1_1:|E2_1891:|E5_2985:|E6_7124:|E7_7603:|L1_4774:", replacement = "", x = target_name, perl = TRUE)) %>%
	readr::type_convert() %>%
	dplyr::arrange(end) %>%
	dplyr::mutate(target_name = factor(target_name, levels = unique(target_name), ordered = TRUE)) %>%
	ggplot(aes(x = read_count, y = insert_size, shape = patient_id_mskcc)) +
	geom_rect(data = do.call(rbind, simulated) %>%
		  	 dplyr::mutate(end = gsub(pattern = "E1_1:|E2_1891:|E5_2985:|E6_7124:|E7_7603:|L1_4774:", replacement = "", x = target_name, perl = TRUE)) %>%
		  	 readr::type_convert() %>%
		  	 dplyr::arrange(end) %>%
		  	 dplyr::mutate(target_name = factor(target_name, levels = unique(target_name), ordered = TRUE)),
		  mapping = aes(xmin = read_count - .025, xmax = read_count + .025, ymin = insert_size -.015, ymax = insert_size + .015, fill = Pr),
		  alpha = .25, inherit.aes = FALSE) +
	scale_fill_continuous(low = "#e7e1ef", high = "#c51b8a", breaks = c(.1, .3, .5, .7, .9)) +
	geom_point(stat = "identity", fill = NA, alpha = .55, size = 2) +
	scale_shape_manual(values = c(0:4,8,10)) +
	xlab("Number of Read Pairs Aligned") +
	ylab("Mean Insert size (bp)") +
	scale_x_continuous(breaks = log10(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6)+1),
			   labels = scientific_10(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6))) +
	scale_y_continuous(limits = log10(c(15, 150)+1),
			   breaks = log10(c(25, 50, 75, 100, 125)+1),
			   labels = c(25, 50, 75, 100, 125)) +
	facet_wrap(~target_name, ncol = 3) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      panel.spacing = unit(1, "lines")) +
	guides(shape = guide_legend(title = "ctDNA", order = 1),
	       fill = guide_colorbar(title = "Posterior\nProbability", order = 2))

pdf("../res/Number_Read_Pairs_by_Insert_Size_with_Posterior_Probability_Baseline_MRD_Negative.pdf", width = 15/(1.35*1.1), height = 9/(1.25*1.1))
print(plot_)
dev.off()

plot_ = do.call(rbind, observed) %>%
	dplyr::rename(sample_uuid = sample_name) %>%
	dplyr::left_join(smry_ %>%
			 dplyr::filter(!duplicated(sample_uuid)) %>%
			 dplyr::select(-insert_size), by = "sample_uuid") %>%
	dplyr::filter(!Is_ctDNA_HPV) %>%
	dplyr::mutate(end = gsub(pattern = "E1_1:|E2_1891:|E5_2985:|E6_7124:|E7_7603:|L1_4774:", replacement = "", x = target_name, perl = TRUE)) %>%
	readr::type_convert() %>%
	dplyr::arrange(end) %>%
	dplyr::mutate(target_name = factor(target_name, levels = unique(target_name), ordered = TRUE)) %>%
	ggplot(aes(x = read_count, y = insert_size, shape = patient_id_mskcc)) +
	geom_rect(data = do.call(rbind, simulated) %>%
		  	 dplyr::mutate(end = gsub(pattern = "E1_1:|E2_1891:|E5_2985:|E6_7124:|E7_7603:|L1_4774:", replacement = "", x = target_name, perl = TRUE)) %>%
		  	 readr::type_convert() %>%
		  	 dplyr::arrange(end) %>%
		  	 dplyr::mutate(target_name = factor(target_name, levels = unique(target_name), ordered = TRUE)),
		  mapping = aes(xmin = read_count - .025, xmax = read_count + .025, ymin = insert_size -.015, ymax = insert_size + .015, fill = Pr),
		  alpha = .25, inherit.aes = FALSE) +
	scale_fill_continuous(low = "#e7e1ef", high = "#c51b8a", breaks = c(.1, .3, .5, .7, .9)) +
	geom_point(stat = "identity", fill = NA, alpha = .55, size = 2) +
	scale_shape_manual(values = c(0:4,8,10:20)) +
	xlab("Number of Read Pairs Aligned") +
	ylab("Mean Insert size (bp)") +
	scale_x_continuous(breaks = log10(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6)+1),
			   labels = scientific_10(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6))) +
	scale_y_continuous(limits = log10(c(15, 150)+1),
			   breaks = log10(c(25, 50, 75, 100, 125)+1),
			   labels = c(25, 50, 75, 100, 125)) +
	facet_wrap(~target_name, ncol = 3) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      panel.spacing = unit(1, "lines")) +
	guides(shape = guide_legend(title = "ctDNA", order = 1),
	       fill = guide_colorbar(title = "Posterior\nProbability", order = 2))

pdf("../res/Number_Read_Pairs_by_Insert_Size_with_Posterior_Probability_Baseline_HPV_Negative.pdf", width = 15/(1.35*1.1), height = 9/(1.25*1.1))
print(plot_)
dev.off()

plot_ = do.call(rbind, observed) %>%
	dplyr::rename(sample_uuid = sample_name) %>%
	dplyr::left_join(smry_ %>%
			 dplyr::filter(!duplicated(sample_uuid)) %>%
			 dplyr::select(-insert_size), by = "sample_uuid") %>%
	dplyr::filter(!Is_ctDNA) %>%
	dplyr::mutate(end = gsub(pattern = "E1_1:|E2_1891:|E5_2985:|E6_7124:|E7_7603:|L1_4774:", replacement = "", x = target_name, perl = TRUE)) %>%
	readr::type_convert() %>%
	dplyr::arrange(end) %>%
	dplyr::mutate(target_name = factor(target_name, levels = unique(target_name), ordered = TRUE)) %>%
	ggplot(aes(x = read_count, y = insert_size, shape = patient_id_mskcc)) +
	geom_rect(data = do.call(rbind, simulated) %>%
		  	 dplyr::mutate(end = gsub(pattern = "E1_1:|E2_1891:|E5_2985:|E6_7124:|E7_7603:|L1_4774:", replacement = "", x = target_name, perl = TRUE)) %>%
		  	 readr::type_convert() %>%
		  	 dplyr::arrange(end) %>%
		  	 dplyr::mutate(target_name = factor(target_name, levels = unique(target_name), ordered = TRUE)),
		  mapping = aes(xmin = read_count - .025, xmax = read_count + .025, ymin = insert_size -.015, ymax = insert_size + .015, fill = Pr),
		  alpha = .25, inherit.aes = FALSE) +
	scale_fill_continuous(low = "#e7e1ef", high = "#c51b8a", breaks = c(.1, .3, .5, .7, .9)) +
	geom_point(stat = "identity", fill = NA, alpha = .55, size = 2) +
	scale_shape_manual(values = c(0:4,8,10:20)) +
	xlab("Number of Read Pairs Aligned") +
	ylab("Mean Insert size (bp)") +
	scale_x_continuous(breaks = log10(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6)+1),
			   labels = scientific_10(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6))) +
	scale_y_continuous(limits = log10(c(15, 150)+1),
			   breaks = log10(c(25, 50, 75, 100, 125)+1),
			   labels = c(25, 50, 75, 100, 125)) +
	facet_wrap(~target_name, ncol = 3) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      panel.spacing = unit(1, "lines")) +
	guides(shape = guide_legend(title = "ctDNA", order = 1),
	       fill = guide_colorbar(title = "Posterior\nProbability", order = 2))

pdf("../res/Number_Read_Pairs_by_Insert_Size_with_Posterior_Probability_Baseline_MRD_HPV_Negative.pdf", width = 15/(1.35*1.1), height = 9/(1.25*1.1))
print(plot_)
dev.off()
