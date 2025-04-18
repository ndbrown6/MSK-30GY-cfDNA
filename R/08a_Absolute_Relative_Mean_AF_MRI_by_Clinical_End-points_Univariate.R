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
	dplyr::group_by(patient_id_mskcc, timepoint_weeks_since_start_of_RT) %>%
	dplyr::summarize(mean_af = mean(mean_af+(1e-5), na.rm = TRUE)) %>%
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
		      neck_dissection_yes_no,
		      hypoxia_resolution,
		      abs_AF_wk0, abs_AF_wk1, abs_AF_wk2, abs_AF_wk3,
		      rel_AF_wk1, rel_AF_wk2, rel_AF_wk3,
		      abs_HPV_wk0, abs_HPV_wk1, abs_HPV_wk2, abs_HPV_wk3,
		      rel_HPV_wk1, rel_HPV_wk2, rel_HPV_wk3,
		      abs_MRI_wk0 = MRI_rawdata_wk0,
		      abs_MRI_wk1 = MRI_rawdata_wk2,
		      abs_MRI_wk2 = MRI_rawdata_wk3,
		      abs_MRI_wk3 = MRI_rawdata_wk4) %>%
	dplyr::mutate(rel_MRI_wk1 = (abs_MRI_wk1/abs_MRI_wk0),
		      rel_MRI_wk2 = (abs_MRI_wk2/abs_MRI_wk0),
		      rel_MRI_wk3 = (abs_MRI_wk3/abs_MRI_wk0)) %>%
	dplyr::mutate(composite_end_point = case_when(
		composite_end_point ~ 1,
		!composite_end_point ~ 0
	)) %>%
	dplyr::mutate(neck_dissection_yes_no = case_when(
		neck_dissection_yes_no == "yes" ~ "1",
		neck_dissection_yes_no == "no" ~ "0"
	)) %>%
	dplyr::mutate(hypoxia_resolution = case_when(
		hypoxia_resolution == "persistent" ~ "1",
		hypoxia_resolution == "resolved" ~ "0"
	)) %>%
	readr::type_convert() %>%
	dplyr::select(-patient_id_mskcc) %>%
	reshape2::melt(id.vars = c("composite_end_point", "neck_dissection_yes_no", "hypoxia_resolution")) %>%
	tidyr::drop_na()

uvars = unique(data_$variable)

p1 = vector(mode = "numeric", length = length(uvars))
for (i in 1:length(uvars)) {
	p1[i] = (data_ %>%
		dplyr::filter(variable == uvars[i]) %>%
		stats::glm(formula = composite_end_point ~ log2(value), data = ., family = binomial(link = "probit")) %>%
		summary() %>%
		.[["coefficients"]])[2,"Pr(>|z|)"]
}

p2 = vector(mode = "numeric", length = length(uvars))
for (i in 1:length(uvars)) {
	p2[i] = (data_ %>%
		dplyr::filter(variable == uvars[i]) %>%
		stats::glm(formula = neck_dissection_yes_no ~ log2(value), data = ., family = binomial(link = "probit")) %>%
		summary() %>%
		.[["coefficients"]])[2,"Pr(>|z|)"]
}

p3 = vector(mode = "numeric", length = length(uvars))
for (i in 1:length(uvars)) {
	p3[i] = (data_ %>%
		dplyr::filter(variable == uvars[i]) %>%
		stats::glm(formula = hypoxia_resolution ~ log2(value), data = ., family = binomial(link = "probit")) %>%
		summary() %>%
		.[["coefficients"]])[2,"Pr(>|z|)"]
}

q1 = vector(mode = "numeric", length = length(uvars))
for (i in 1:length(uvars)) {
	q1[i] = data_ %>%
		dplyr::filter(variable == uvars[i]) %>%
		wilcox.test(formula = value ~ composite_end_point, data = ., exact = FALSE) %>%
		.[["p.value"]]
}

q2 = vector(mode = "numeric", length = length(uvars))
for (i in 1:length(uvars)) {
	q2[i] = data_ %>%
		dplyr::filter(variable == uvars[i]) %>%
		wilcox.test(formula = value ~ neck_dissection_yes_no, data = ., exact = FALSE) %>%
		.[["p.value"]]
}

q3 = vector(mode = "numeric", length = length(uvars))
for (i in 1:length(uvars)) {
	q3[i] = data_ %>%
		dplyr::filter(variable == uvars[i]) %>%
		wilcox.test(formula = value ~ hypoxia_resolution, data = ., exact = FALSE) %>%
		.[["p.value"]]
}

res_ = dplyr::tibble(variables = uvars,
		     `De-escalation failure` = p2,
		     `PET response` = p3,
		     `Risk group` = p1) %>%
       dplyr::mutate(Method = "log") %>%
       dplyr::bind_rows(dplyr::tibble(variables = uvars,
				      `De-escalation failure` = q2,
				      `PET response` = q3,
				      `Risk group` = q1) %>%
			dplyr::mutate(Method = "u")) %>%
       dplyr::add_row(`variables` = "rel_AF_wk0",
		      `De-escalation failure` = 1,
		      `PET response` = 1,
		      `Risk group` = 1,
		      Method = "log") %>%
       dplyr::add_row(`variables` = "rel_MRI_wk0",
		      `De-escalation failure` = 1,
		      `PET response` = 1,
		      `Risk group` = 1,
		      Method = "log") %>%
       dplyr::add_row(`variables` = "rel_HPV_wk0",
		      `De-escalation failure` = 1,
		      `PET response` = 1,
		      `Risk group` = 1,
		      Method = "log") %>%
       dplyr::add_row(`variables` = "rel_AF_wk0",
		      `De-escalation failure` = 1,
		      `PET response` = 1,
		      `Risk group` = 1,
		      Method = "u") %>%
       dplyr::add_row(`variables` = "rel_MRI_wk0",
		      `De-escalation failure` = 1,
		      `PET response` = 1,
		      `Risk group` = 1,
		      Method = "u") %>%
       dplyr::add_row(`variables` = "rel_HPV_wk0",
		      `De-escalation failure` = 1,
		      `PET response` = 1,
		      `Risk group` = 1,
		      Method = "u")

plot_ = res_ %>%
	dplyr::filter(Method == "log") %>%
	dplyr::select(-Method) %>%
	reshape2::melt() %>%
	dplyr::mutate(value = -log10(value)) %>%
	dplyr::mutate(week = case_when(
		grepl("wk0", variables, fixed = TRUE) ~ "Pre-treat\n-ment",
		grepl("wk1", variables, fixed = TRUE) ~ "Week 1",
		grepl("wk2", variables, fixed = TRUE) ~ "Week 2",
		grepl("wk3", variables, fixed = TRUE) ~ "Week 3"
	)) %>%
	dplyr::mutate(assay = case_when(
		grepl("AF", variables, fixed = TRUE) ~ "PCM",
		grepl("HPV", variables, fixed = TRUE) ~ "HPV",
		grepl("MRI", variables, fixed = TRUE) ~ "MRI"
	)) %>%
	dplyr::mutate(mode = case_when(
		grepl("rel", variables, fixed = TRUE) ~ "Relative",
		grepl("abs", variables, fixed = TRUE) ~ "Absolute",
	)) %>%
	dplyr::select(-variables) %>%
	dplyr::mutate(mode = factor(mode, levels = c("Absolute", "Relative"), ordered = TRUE)) %>%
	dplyr::mutate(assay = factor(assay, levels = c("HPV", "PCM", "MRI"), ordered = TRUE)) %>%
	dplyr::mutate(variable = factor(variable, levels = c("Risk group", "PET response", "De-escalation failure"), ordered = TRUE)) %>%
	ggplot(aes(x = week, y = variable, fill = value)) +
	geom_tile(stat = "identity", color = "white", size = .5) +
	scale_fill_viridis(breaks = seq(from = 0, to = 3, by = .75)) +
	xlab("") +
	ylab("") +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	facet_grid(mode~assay) +
	guides(fill = guide_colourbar(title = expression(-Log[10](p))))

pdf(file = "../res/Univariate_AF_HPV_MRI_by_Clinical_End_Points.pdf", width = 10/1.2, height = 5/1.25)
print(plot_)
dev.off()

plot_ = data_ %>%
	dplyr::filter(variable == "abs_AF_wk2") %>%
	dplyr::mutate(composite_end_point = factor(composite_end_point, levels = c(0, 1), ordered = TRUE)) %>%
	ggplot(aes(x = composite_end_point, y = 100*value, fill = composite_end_point)) +
	geom_half_dotplot(binaxis = "y", stackdir = "up", dotsize = .65, right = TRUE, method = "dotdensity", fill = "white") +
	geom_half_violin(stat = "half_ydensity", side = "l", trim = FALSE, scale = "count", color = "black") +
	scale_fill_brewer(type = "qual", palette = 7) +
	scale_x_discrete(breaks = c(0, 1),
			 labels = c("Low", "High")) +
	scale_y_log10(limits = c(NA, 1000),
		      breaks = c(.001, .01, .1, 1, 10),
		      labels = scientific_10) +
	xlab("Risk group") +
	ylab("Mean AF (%)") +
	geom_signif(stat = "signif",
		    comparisons = list(c("0", "1")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    size = .55, textsize = 4, color = "grey25", vjust = -1, 
		    y_position = log10(500)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = FALSE, fill = FALSE)

pdf(file = "../res/Raw_Mean_AF_at_Wk2_by_Clinical_End_Points.pdf", width = 3, height = 4)
print(plot_)
dev.off()

plot_ = data_ %>%
	dplyr::filter(variable == "abs_MRI_wk2") %>%
	dplyr::mutate(composite_end_point = factor(composite_end_point, levels = c(0, 1), ordered = TRUE)) %>%
	ggplot(aes(x = composite_end_point, y = 100*value, fill = composite_end_point)) +
	geom_half_dotplot(binaxis = "y", stackdir = "up", dotsize = .65, right = TRUE, method = "dotdensity", fill = "white") +
	geom_half_violin(stat = "half_ydensity", side = "l", trim = FALSE, scale = "count", color = "black") +
	scale_fill_brewer(type = "qual", palette = 7) +
	scale_x_discrete(breaks = c(0, 1),
			 labels = c("Low", "High")) +
	scale_y_log10(limits = c(1E4, 1E8),
		      breaks = c(1E4, 1E5, 1E6, 1E7),
		      labels = scientific_10) +
	xlab("Risk group") +
	ylab("MRI Volume (mm3)") +
	geom_signif(stat = "signif",
		    comparisons = list(c("0", "1")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    size = .55, textsize = 4, color = "grey25", vjust = -1, 
		    y_position = log10(7E7)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = FALSE, fill = FALSE)

pdf(file = "../res/Raw_MRI_Volume_at_Wk2_by_Clinical_End_Points.pdf", width = 3, height = 4)
print(plot_)
dev.off()

plot_ = data_ %>%
	dplyr::filter(variable == "abs_HPV_wk2") %>%
	dplyr::mutate(composite_end_point = factor(composite_end_point, levels = c(0, 1), ordered = TRUE)) %>%
	ggplot(aes(x = composite_end_point, y = 100*value, fill = composite_end_point)) +
	geom_half_dotplot(binaxis = "y", stackdir = "up", dotsize = .65, right = TRUE, method = "dotdensity", fill = "white") +
	geom_half_violin(stat = "half_ydensity", side = "l", trim = FALSE, scale = "count", color = "black") +
	scale_fill_brewer(type = "qual", palette = 7) +
	scale_x_discrete(breaks = c(0, 1),
			 labels = c("Low", "High")) +
	scale_y_log10(limits = c(1E2, 1E12),
		      breaks = c(1E2, 1E4, 1E6, 1E8, 1E10),
		      labels = scientific_10) +
	xlab("Risk group") +
	ylab("cfDNA Aligned HPV Read Pairs") +
	geom_signif(stat = "signif",
		    comparisons = list(c("0", "1")),
		    test = "wilcox.test",
		    test.args = list(alternative = "two.sided", exact = FALSE),
		    size = .55, textsize = 4, color = "grey25", vjust = -1, 
		    y_position = log10(3E11)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = FALSE, fill = FALSE)

pdf(file = "../res/Raw_HPV_Reads_Pairs_at_Wk2_by_Clinical_End_Points.pdf", width = 3, height = 4)
print(plot_)
dev.off()
