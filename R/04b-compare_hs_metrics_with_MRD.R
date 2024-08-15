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

#########################################################
# Training set only HPV-16 +ve
#########################################################
insert_size_metrics = readr::read_tsv(file = url_insert_metrics_by_gene, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert() %>%
		      dplyr::rename(sample_name = SAMPLE_NAME) %>%
		      dplyr::filter(TARGET_REGION != "NC001526.4:3372-4794")

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
	dplyr::mutate(aligned_reads = log10(READ_PAIRS+1)) %>%
	dplyr::mutate(insert_size = log10(MEAN_INSERT_SIZE+1)) %>%
	reshape2::dcast(formula = sample_name + Is_ctDNA ~ gene_name + start_end, value.var = "aligned_reads")

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
	ggplot(aes(x +1, y +1, color = Is_ctDNA, shape = Is_ctDNA)) +
	geom_abline(slope = 1, intercept = 0, linetype = 1, size = 1.5, alpha = .55, color = "goldenrod3") +
	geom_vline(xintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_hline(yintercept = 1, color = "black", linetype = 3, size = .5, alpha = .75) +
	geom_point(stat = "identity", alpha = .75, size = 1.95, fill = NA) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_shape_manual(values = c(21, 24)) +
	xlab("Number of Read Pairs Aligned") +
	ylab("Number of Read Pairs Aligned") +
	scale_x_continuous(limits = c(1,6),
			   breaks = log10(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6)+1),
			   labels = scientific_1e(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6))) +
	scale_y_continuous(limits = c(1,6),
			   breaks = log10(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6)+1),
			   labels = scientific_1e(c(0, 1E1, 1E2, 1E3, 1E4, 1E5, 1E6))) +
	stat_cor(data = do.call(bind_rows, data_) %>%
		 	dplyr::mutate(xlab = factor(xlab, levels = target_contigs, ordered = TRUE)) %>%
		 	dplyr::mutate(ylab = factor(ylab, levels = target_contigs, ordered = TRUE)),
		 mapping = aes(x +1, y +1),
		 method = "spearman", size = 3, inherit.aes = FALSE) +
	facet_wrap(ylab ~ xlab, scales = "free") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 14),
	      axis.title.y = element_text(margin = margin(r = 20), size = 14),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      panel.spacing = unit(1, "lines"),
	      strip.background = element_rect(colour="white", fill="white")) +
	guides(color = guide_legend(title = "ctDNA"),
	       shape = guide_legend(title = "ctDNA"))

pdf(file = "../res/Number_Read_Pairs_Aligned_by_targets_HPV-16.pdf", width = 3.25*2.75, height = 3.25*3)
print(plot_)
dev.off()

r = do.call(rbind, lapply(data_, function(x) { dplyr::tibble(rho = cor.test(x$x, x$y, method = "spearman", na.rm=TRUE)$estimate,
							     xlab = x$xlab[1],
							     ylab = x$ylab[1])} )) %>%
    as.data.frame()
m = matrix(1, nrow = 6, ncol = 6, dimnames = list(unique(c(r$xlab, r$ylab)),
						  unique(c(r$xlab, r$ylab))))
for (i in 1:nrow(r)) {
	m[r[i,"xlab"], r[i,"ylab"]] = r[i,"rho"]
	m[r[i,"ylab"], r[i,"xlab"]] = r[i,"rho"]
}

pdf(file = "../res/Cross_Correlation_.targets_HPV-16.pdf", width = 3.5, height = 3)
draw(Heatmap(matrix = m %>%
	     	      as.matrix(),
	     col = c(rev(viridis(n = 10)), rep("#440154FF", 3)),
	     name = "rho",
	     rect_gp = gpar(col = "white", lwd = 1),
	     border = NA,
	     
	     cluster_rows = TRUE,
	     clustering_distance_rows = "canberra",
	     clustering_method_rows = "ward.D",
	     show_row_names = TRUE,
	     row_names_side = "right",
	     row_names_gp = gpar(fontsize = 8),
	     
	     cluster_columns = TRUE,
	     clustering_distance_columns = "canberra",
	     clustering_method_columns = "ward.D",
	     show_column_names = TRUE,
	     column_names_side = "bottom",
	     column_names_gp = gpar(fontsize = 8),

	     use_raster = FALSE,
	     show_heatmap_legend = TRUE))
dev.off()