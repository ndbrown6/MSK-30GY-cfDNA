#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::mutate(sample_uuid = paste0(sample_id_mskcc, "-", sample_id_invitae))

mrd_smry = readr::read_tsv(file = url_mrd_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert()

posterior_smry = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		 readr::type_convert()

plot_ = posterior_smry %>%
	dplyr::left_join(mrd_smry %>%
			 dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			 dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
	dplyr::filter(timepoint_days_since_start_of_RT>=0) %>%
	dplyr::mutate(Is_ctDNA = case_when(
		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
		`MRD-Landmark_Result` == "ABSENT" ~ "-ve"
	)) %>%
	dplyr::arrange(as.numeric(gsub("CTMS-", "", patient_id_mskcc))) %>%
	dplyr::mutate(patient_id_mskcc = factor(patient_id_mskcc, levels = unique(patient_id_mskcc), ordered = TRUE)) %>%
	ggplot(aes(x = timepoint_days_since_start_of_RT + 1, y = `Pr(x=1)`)) +
	geom_line(stat = "identity", size = .5, alpha = 1, color = "#e41a1c") +
	geom_jitter(data = posterior_smry %>%
			   dplyr::left_join(mrd_smry %>%
					    dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
					    dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
			   dplyr::filter(timepoint_days_since_start_of_RT>=0) %>%
			   dplyr::mutate(Is_ctDNA = case_when(
				   		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
				   		`MRD-Landmark_Result` == "ABSENT" ~ "-ve"
			   )) %>%
		    	   dplyr::arrange(as.numeric(gsub("CTMS-", "", patient_id_mskcc))) %>%
		    	   dplyr::mutate(patient_id_mskcc = factor(patient_id_mskcc, levels = unique(patient_id_mskcc), ordered = TRUE)),
		    mapping = aes(x = timepoint_days_since_start_of_RT + 1, y = `Pr(x=1)`, shape = Is_ctDNA),
		    width = 0, height = 0, inherit.aes = FALSE, color = "#e41a1c") +
	scale_shape_manual(values = c("-ve" = 1, "+ve" = 2)) +
	scale_x_log10(limits = c(1, 1096),
		      breaks = c(0, 7, 30, 182, 730)+1,
		      labels = c("RT", "7d", "1m", "6m", "2y")) +
	scale_y_continuous(limits = c(0, 1)) +
	xlab("Time since start of RT") +
	ylab("Posterior Probability") +
	facet_wrap(~patient_id_mskcc) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(shape = guide_legend(title = "ctDNA", order = 1))
	
pdf("../res/Time_Point_by_Posterior_Probability.pdf", width = 15, height = 10)
print(plot_)
dev.off()

plot_ = posterior_smry %>%
	dplyr::left_join(mrd_smry %>%
			 dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			 dplyr::select(sample_name, `MRD-Landmark_Result`, `MRD-Landmark_MRD_Quantification`), by = "sample_name") %>%
	dplyr::filter(timepoint_days_since_start_of_RT >= 0) %>%
	dplyr::mutate(Is_ctDNA = case_when(
		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
		`MRD-Landmark_Result` == "ABSENT" ~ "-ve"
	)) %>%
	dplyr::rename(frac_ctDNA = `MRD-Landmark_MRD_Quantification`) %>%
	dplyr::select(patient_id_mskcc, timepoint_days_since_start_of_RT, mean_af, Is_ctDNA, frac_ctDNA) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "timepoint_days_since_start_of_RT", "Is_ctDNA"), measure.vars = c("mean_af", "frac_ctDNA"),
		       variable.name = "summary", value.name = "AF") %>%
	dplyr::arrange(as.numeric(gsub("CTMS-", "", patient_id_mskcc))) %>%
	dplyr::mutate(patient_id_mskcc = factor(patient_id_mskcc, levels = unique(patient_id_mskcc), ordered = TRUE)) %>%
	dplyr::mutate(AF = ifelse(AF>0.05, 0.05, AF)) %>%
	dplyr::mutate(AF = ifelse(AF<0.00005, 0.00005, AF)) %>%
	dplyr::mutate(summary = case_when(
		summary == "mean_af" ~ "Mean",
		summary == "frac_ctDNA" ~ "ctDNA\nFraction"
	)) %>%
	ggplot(aes(x = timepoint_days_since_start_of_RT + 1, y = (100*AF), color = summary)) +
	geom_line(stat = "identity", size = .5, alpha = 1) +
	geom_jitter(data = posterior_smry %>%
		    	  dplyr::left_join(mrd_smry %>%
			 		   dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
					   dplyr::select(sample_name, `MRD-Landmark_Result`, `MRD-Landmark_MRD_Quantification`), by = "sample_name") %>%
		    	  dplyr::filter(timepoint_days_since_start_of_RT >= 0) %>%
		    	  dplyr::mutate(Is_ctDNA = case_when(
				  		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
				  		`MRD-Landmark_Result` == "ABSENT" ~ "-ve"
			  )) %>%
		    	  dplyr::rename(frac_ctDNA = `MRD-Landmark_MRD_Quantification`) %>%
		    	  dplyr::select(patient_id_mskcc, timepoint_days_since_start_of_RT, mean_af, Is_ctDNA, frac_ctDNA) %>%
		    	  reshape2::melt(id.vars = c("patient_id_mskcc", "timepoint_days_since_start_of_RT", "Is_ctDNA"), measure.vars = c("mean_af", "frac_ctDNA"),
					 variable.name = "summary", value.name = "AF") %>%
		    	  dplyr::arrange(as.numeric(gsub("CTMS-", "", patient_id_mskcc))) %>%
		    	  dplyr::mutate(patient_id_mskcc = factor(patient_id_mskcc, levels = unique(patient_id_mskcc), ordered = TRUE)) %>%
		    	  dplyr::mutate(AF = ifelse(AF>0.05, 0.05, AF)) %>%
			  dplyr::mutate(AF = ifelse(AF<0.00005, 0.00005, AF)) %>%
		    	  dplyr::mutate(summary = case_when(
				  		summary == "mean_af" ~ "Mean",
				  		summary == "frac_ctDNA" ~ "ctDNA\nFraction"
			  )),
		    mapping = aes(x = timepoint_days_since_start_of_RT + 1, y = (100*AF), color = summary, shape = Is_ctDNA),
		    width = 0, height = 0, inherit.aes = FALSE) +
	scale_shape_manual(values = c("-ve" = 1, "+ve" = 2)) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_x_log10(limits = c(1, 1096),
		      breaks = c(0, 7, 30, 182, 730)+1,
		      labels = c("RT", "7d", "1m", "6m", "2y")) +
	scale_y_log10(limits = c(0.001, 10),
		      breaks = c(1E-3, 1E-2, 1E-1, 1, 10),
		      labels = c("0.001", "0.01", "0.1", "1", "10")) +
	xlab("Time since start of RT") +
	ylab("Allele Fraction (%)") +
	facet_wrap(~patient_id_mskcc) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(color = guide_legend(title = "AF", override.aes = list(shape = NA), order = 1),
	       shape = guide_legend(title = "ctDNA", order = 2))
	
pdf("../res/Time_Point_by_AF.pdf", width = 15, height = 10)
print(plot_)
dev.off()
