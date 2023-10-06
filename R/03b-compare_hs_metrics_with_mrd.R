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

hs_metrics = readr::read_tsv(file = url_hs_target_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	     readr::type_convert() %>%
	     dplyr::select(sample_name = SAMPLE_NAME, chromosome = chrom, start, end, aligned_reads = read_count) %>%
	     dplyr::filter(chromosome %in% target_contigs) %>%
	     fuzzyjoin::genome_left_join(bed_file, by = c("chromosome", "start", "end")) %>%
	     dplyr::select(-chromosome.y, -start.y, -end.y) %>%
	     dplyr::rename(chromosome = chromosome.x,
			   start = start.x,
			   end = end.x)

hs_metrics_ft = readr::read_tsv(file = url_hs_target_metrics_ft, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      	readr::type_convert() %>%
		dplyr::select(sample_name = SAMPLE_NAME, chromosome = chrom, start, end, aligned_reads = read_count, fragment_length = FRAGMENT_LENGTH) %>%
	      	dplyr::filter(chromosome %in% target_contigs) %>%
	        fuzzyjoin::genome_left_join(bed_file, by = c("chromosome", "start", "end")) %>%
	        dplyr::select(-chromosome.y, -start.y, -end.y) %>%
	        dplyr::rename(chromosome = chromosome.x,
			      start = start.x,
			      end = end.x)

mutation_smry = readr::read_tsv(file = url_mutation_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		readr::type_convert()

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


for (i in 1:length(target_contigs)) {
	
	n = hs_metrics_ft %>%
	    dplyr::filter(chromosome == target_contigs[i]) %>%
	    dplyr::mutate(uuid = paste0(chromosome, "-", start, ":", end)) %>%
	    dplyr::filter(!duplicated(uuid)) %>%
	    nrow()
	
	if (names(target_contigs)[i] == "HPV-16") {
	
		smry_ = hs_metrics_ft %>%
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
			dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD) %>%
			dplyr::mutate(start_end = paste0(start, ":", end)) %>%
			dplyr::mutate(Is_ctDNA = factor(Is_ctDNA, levels = c("+ve", "-ve"), ordered = TRUE)) %>%
			dplyr::filter(chromosome == "HPV-16") %>%
			dplyr::mutate(Is_ctDNA = case_when(
				hpv_type_wes_wgs == "HPV-16" & Is_ctDNA == "+ve" ~ "+ve",
				TRUE ~ "-ve"
			)) %>%
			reshape2::dcast(formula = sample_name + Is_ctDNA ~ gene_name + start_end, value.var = "aligned_reads")

		fit = list()
		for (j in 3:ncol(smry_)) {
			fit[[j-2]] = smry_ %>%
				     dplyr::select(all_of(c(2,j)))
			colnames(fit[[j-2]]) = c("Is_ctDNA", "aligned_reads")
			fit[[j-2]] = fit[[j-2]] %>%
				     rpart::rpart(formula = Is_ctDNA ~ aligned_reads, data = ., method = "class")
			fit[[j-2]] = dplyr::tibble(hpv_type_wes_wgs = unique(smry_ft$hpv_type_wes_wgs),
						   chromosome = rep(names(target_contigs)[i], length(unique(smry_ft$hpv_type_wes_wgs))),
						   start_end = rep(strsplit(colnames(smry_)[j], "_", fixed = TRUE)[[1]][2], length(unique(smry_ft$hpv_type_wes_wgs))),
						   gene_name = rep(strsplit(colnames(smry_)[j], "_", fixed = TRUE)[[1]][1], length(unique(smry_ft$hpv_type_wes_wgs))),
						   xintercept = rep(fit[[j-2]]$splits[1,"index"], length(unique(smry_ft$hpv_type_wes_wgs))))
		}
		fit = do.call(rbind, fit)

		plot_ = hs_metrics_ft %>%
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
			dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD & chromosome == names(target_contigs)[i]) %>%
			dplyr::mutate(start_end = paste0(start, ":", end)) %>%
			dplyr::arrange(start) %>%
			dplyr::mutate(start_end = factor(start_end, levels = unique(start_end), ordered = TRUE)) %>%
			dplyr::arrange(desc(Is_ctDNA)) %>%
			ggplot(aes(x = aligned_reads+1, y = 100*mean_af, shape = hpv_type_wes_wgs, color = Is_ctDNA)) +
			geom_vline(data = fit, aes(xintercept = xintercept), alpha = 1, linetype = 3) +
			geom_point(stat = "identity", fill = NA, alpha = 1, size = 2.5) +
			scale_shape_manual(values = c(1, 2, 3, 4, 5, 6)) +
			scale_color_brewer(type = "qual", palette = 6) +
			xlab("Number of Read Pairs Aligned") +
			ylab("Mean Allele Fraction (%)") +
			scale_x_log10(limits = c(.9, 1e5),
				      labels = scientific_10) +
			scale_y_log10(limits = c(1e-4, 100),
				      breaks = c(1e-4, 1e-3, 1e-1, 1, 10, 100),
				      labels = c("ND", ".001", ".01", "1", "10", "100")) +
			facet_grid(hpv_type_wes_wgs ~ chromosome + gene_name + start_end) +
			theme(axis.title.x = element_text(margin = margin(t = 20)),
			      axis.title.y = element_text(margin = margin(r = 20)),
			      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			guides(shape = guide_legend(title = "Tumor\nHPV Type"),
			       color = guide_legend(title = "ctDNA", override.aes = list(alpha = 1)))

		pdf(file = paste0("../res/Number_Read_Pairs_Aligned_by_Mean_AF_by_target_", names(target_contigs)[i],".pdf"), width = (22.5/20) * (n) + 2.5, height = 15/2)
		print(plot_)
		dev.off()
		
	} else {
		
		plot_ = hs_metrics_ft %>%
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
			dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD & chromosome == names(target_contigs)[i]) %>%
			dplyr::mutate(start_end = paste0(start, ":", end)) %>%
			dplyr::arrange(start) %>%
			dplyr::mutate(start_end = factor(start_end, levels = unique(start_end), ordered = TRUE)) %>%
			dplyr::arrange(desc(Is_ctDNA)) %>%
			ggplot(aes(x = aligned_reads+1, y = 100*mean_af, shape = hpv_type_wes_wgs, color = Is_ctDNA)) +
			geom_point(stat = "identity", fill = NA, alpha = 1, size = 2.5) +
			scale_shape_manual(values = c(1, 2, 3, 4, 5, 6)) +
			scale_color_brewer(type = "qual", palette = 6) +
			xlab("Number of Read Pairs Aligned") +
			ylab("Mean Allele Fraction (%)") +
			scale_x_log10(limits = c(.9, 1e5),
				      labels = scientific_10) +
			scale_y_log10(limits = c(1e-4, 100),
				      breaks = c(1e-4, 1e-3, 1e-1, 1, 10, 100),
				      labels = c("ND", ".001", ".01", "1", "10", "100")) +
			facet_grid(hpv_type_wes_wgs ~ chromosome + gene_name + start_end) +
			theme(axis.title.x = element_text(margin = margin(t = 20)),
			      axis.title.y = element_text(margin = margin(r = 20)),
			      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
			guides(shape = guide_legend(title = "Tumor\nHPV Type"),
			       color = guide_legend(title = "ctDNA", override.aes = list(alpha = 1)))

		pdf(file = paste0("../res/Number_Read_Pairs_Aligned_by_Mean_AF_by_target_", names(target_contigs)[i],".pdf"), width = (22.5/20) * (n) + 2.5, height = 15/2)
		print(plot_)
		dev.off()
		
	}

}

smry_ = hs_metrics_ft %>%
	dplyr::mutate(chromosome = case_when(
			chromosome == "J04353.1" ~ "HPV-31",
			chromosome == "M12732.1" ~ "HPV-33",
			chromosome == "NC001357.1" ~ "HPV-18",
			chromosome == "NC001526.4" ~ "HPV-16",
			chromosome == "X74477.1" ~ "HPV-35"
	)) %>%
	dplyr::left_join(smry_ft %>%
			 dplyr::select(sample_name, Is_ctDNA, hpv_type_wes_wgs), by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = ifelse(is.na(Is_ctDNA), "?", Is_ctDNA)) %>%
	dplyr::filter(Is_ctDNA != "?") %>%
	dplyr::filter(!(Is_ctDNA=="+ve" & hpv_type_wes_wgs!="HPV-16")) %>%
	dplyr::filter(fragment_length == FRAGMENT_LENGTH_THRESHOLD & chromosome == "HPV-16") %>%
	reshape2::dcast(formula = sample_name + Is_ctDNA ~ gene_name + start + end, value.var = "aligned_reads")

target_contigs = colnames(smry_)[-c(1,2)]
target_contigs = target_contigs[order(unlist(lapply(target_contigs, str_split, split = "_", n = 2)))]

ii = 1
data_ = list()
for (i in 1:(length(target_contigs)-1)) {
	for (j in (i+1):length(target_contigs)) {
		data_[[ii]] = smry_ %>%
			      dplyr::select(all_of(c("sample_name", "Is_ctDNA", (target_contigs)[i], (target_contigs)[j]))) %>%
			      dplyr::rename(x = (target_contigs)[i],
					    y = (target_contigs)[j]) %>%
			      dplyr::mutate(xlab = (target_contigs)[i],
					    ylab = (target_contigs)[j])
		ii = ii + 1
	}
}

plot_ = do.call(bind_rows, data_) %>%
	dplyr::mutate(xlab = factor(xlab, levels = target_contigs, ordered = TRUE)) %>%
	dplyr::mutate(ylab = factor(ylab, levels = target_contigs, ordered = TRUE)) %>%
	ggplot(aes(x +1, y +1, color = Is_ctDNA)) +
	geom_abline(slope = 1, intercept = 0, linetype = 1, size = .75, alpha = .75, color = "goldenrod3") +
	geom_vline(xintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_hline(yintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_point(stat = "identity", shape = 21, alpha = .55, size = 1.5, fill = NA) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("Number of Read Pairs Aligned") +
	ylab("Number of Read Pairs Aligned") +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	facet_grid(ylab ~ xlab) +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	guides(color = guide_legend(title = "ctDNA"))

pdf(file = "../res/Number_Read_Pairs_Aligned_by_targets_HPV-16.pdf", width = 15, height = 15)
print(plot_)
dev.off()
