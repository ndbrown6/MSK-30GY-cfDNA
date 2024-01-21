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

manifest = manifest %>%
	   dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	   dplyr::mutate(patient_id_mskcc = case_when(
		   is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
		   TRUE ~ patient_id_mskcc
	   )) %>%
	   dplyr::mutate(SAMPLE_NAME = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	   dplyr::left_join(hpv_smry, by = "patient_id_mskcc")

aln_metrics = readr::read_tsv(file = url_aln_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::filter(CATEGORY == "PAIR")
	      
insert_size_metrics = readr::read_tsv(file = url_insert_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		      readr::type_convert()

	      
mutation_smry = readr::read_tsv(file = url_mutation_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		readr::type_convert()

smry_ = mutation_smry %>%
	dplyr::filter(FILTER == "PASS") %>%
	dplyr::group_by(Tumor_Sample_Barcode) %>%
	dplyr::summarize(mean_af = mean(t_maf),
		         max_af = max(t_maf)) %>%
	dplyr::ungroup() %>%
	dplyr::rename(SAMPLE_NAME = Tumor_Sample_Barcode) %>%
	dplyr::left_join(aln_metrics, by = "SAMPLE_NAME") %>%
	dplyr::left_join(insert_size_metrics, by = "SAMPLE_NAME") %>%
	dplyr::mutate(ctDNA = case_when(
		SAMPLE_NAME == "21-020-04087-GER2107031" ~ "+ve",
		SAMPLE_NAME %in% c("20-293-02575-GER2110407", "20-027-04124-GER2106983", "20-090-01921-GER2107058", "20-338-01963-GER2107081",
				   "21-336-01867-GER2213061", "20-098-01693-GER2106917", "20-106-00933-GER2110283", "21-063-02905-GER2110363",
				   "21-056-03771-GER2107040", "20-133-03273-GER2110315", "20-167-01502-GER2110203", "20-188-03485-GER2110226",
				   "20-216-02037-GER2107110", "20-234-03586-GER2107097", "20-232-02611-GER2106916", "20-344-03445-GER2107077",
				   "21-103-01853-GER2106755", "21-364-01538-GER2213112", "21-032-01135-GER2110367", "22-024-02646-GER2213114",
				   "20-198-02689-GER2110289", "20-226-01995-GER2106882", "18-316-02771-GER2110425", "18-340-02846-GER2110490",
				   "19-337-04606-GER2110241", "20-280-03572-GER2110172", "18-316-01478-GER2110415", "18-334-01797-GER2106689",
				   "18-355-03108-GER2106683", "19-053-03844-GER2106674", "20-311-03793-GER2106806", "20-195-05002-GER2110247",
				   "20-141-02240-GER2106953", "21-110-03786-GER2106754", "19-161-03236-GER2106824", "19-268-01264-GER2106820",
				   "19-280-02709-GER2107020", "21-195-02314-GER2110519", "20-220-03565-GER2107050", "19-234-03680-GER2107048",
				   "20-042-01779-GER2110291", "19-270-05067-GER2110377", "20-015-02809-GER2110359") ~ "-ve",
		TRUE ~ "Unknown"
	))

plot_ = smry_ %>%
	dplyr::mutate(mean_af = case_when(
		mean_af == 0 ~ 1E-6,
		TRUE ~ mean_af
	)) %>%
	ggplot(aes(x = 100*mean_af, y = TOTAL_READS/2, fill = ctDNA, color = ctDNA, alpha = ctDNA)) +
	geom_point(stat = "identity", shape = 21, size = 2) +
	scale_fill_manual(values = c("+ve" = "#e41a1c", "-ve" = "#377eb8", "Unknown" = "white")) +
	scale_color_manual(values = c("+ve" = "black", "-ve" = "black", "Unknown" = "grey50")) +
	scale_alpha_manual(values = c("+ve" = 1, "-ve" = 1, "Unknown" = .55)) +
	xlab("Mean AF (%)") +
	ylab("Total Read Pairs Aligned") +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(labels = scientific_10) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)))
	      
pdf(file = "../res/Number_of_Read_Pairs_by_AF.pdf", width = 6, height = 5)
print(plot_)
dev.off()

plot_ = smry_ %>%
	dplyr::mutate(mean_af = case_when(
		mean_af == 0 ~ 1E-6,
		TRUE ~ mean_af
	)) %>%
	ggplot(aes(x = 100*mean_af, y = MEAN_INSERT_SIZE, fill = ctDNA, color = ctDNA, alpha = ctDNA)) +
	geom_point(stat = "identity", shape = 21, size = 2) +
	scale_fill_manual(values = c("+ve" = "#e41a1c", "-ve" = "#377eb8", "Unknown" = "white")) +
	scale_color_manual(values = c("+ve" = "black", "-ve" = "black", "Unknown" = "grey50")) +
	scale_alpha_manual(values = c("+ve" = 1, "-ve" = 1, "Unknown" = .55)) +
	xlab("Mean AF (%)") +
	ylab("Mean Insert size") +
	scale_x_log10(labels = scientific_10) +
	scale_y_log10(limits = c(24, 150),
		      breaks = c(25, 50, 75, 100, 125),
		      labels = c(25, 50, 75, 100, 125)) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)))
	      
pdf(file = "../res/Mean_Insert_Size_by_AF.pdf", width = 6, height = 5)
print(plot_)
dev.off()
