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

posterior_smry = readr::read_tsv(file = "../res/posterior_probability_all.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
		 readr::type_convert()

Baseline_ctDNA_Neg = posterior_smry %>%
		     dplyr::left_join(mrd_smry %>%
				      dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
				      dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
		     dplyr::filter(timepoint_days_since_start_of_RT <= 0) %>%
		     dplyr::mutate(Is_ctDNA = case_when(
			     		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
			     		`MRD-Landmark_Result` == "ABSENT" ~ "-ve"
		     )) %>%
		     dplyr::group_by(patient_id_mskcc) %>%
		     dplyr::summarize(Is_ctDNA = any(Is_ctDNA=="+ve")) %>%
		     dplyr::filter(!Is_ctDNA) %>%
		     dplyr::mutate(Is_Baseline = TRUE) %>%
		     dplyr::select(patient_id_mskcc, Is_Baseline)

plot_ = posterior_smry %>%
	dplyr::left_join(mrd_smry %>%
			 dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			 dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%	
	dplyr::mutate(Is_ctDNA = case_when(
		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
		`MRD-Landmark_Result` == "ABSENT" ~ "-ve"
	)) %>%
	dplyr::left_join(Baseline_ctDNA_Neg, by = "patient_id_mskcc") %>%
	dplyr::filter(Is_Baseline) %>%
	dplyr::arrange(as.numeric(gsub("CTMS-", "", patient_id_mskcc))) %>%
	dplyr::mutate(patient_id_mskcc = factor(patient_id_mskcc, levels = unique(patient_id_mskcc), ordered = TRUE)) %>%
	dplyr::mutate(timepoint_days_since_start_of_RT = case_when(
		timepoint_days_since_start_of_RT > 200 ~ 200,
		TRUE ~ timepoint_days_since_start_of_RT
	)) %>%
	ggplot(aes(x = timepoint_days_since_start_of_RT, y = `Pr(x=1)`)) +
	geom_line(stat = "identity", size = .5, alpha = 1, color = "#e41a1c") +
	geom_jitter(data = posterior_smry %>%
			   dplyr::left_join(mrd_smry %>%
					    dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
					    dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
			   dplyr::mutate(Is_ctDNA = case_when(
				   		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
				   		`MRD-Landmark_Result` == "ABSENT" ~ "-ve"
			   )) %>%
		    	   dplyr::left_join(Baseline_ctDNA_Neg, by = "patient_id_mskcc") %>%
		    	   dplyr::filter(Is_Baseline) %>%
		    	   dplyr::arrange(as.numeric(gsub("CTMS-", "", patient_id_mskcc))) %>%
		    	   dplyr::mutate(patient_id_mskcc = factor(patient_id_mskcc, levels = unique(patient_id_mskcc), ordered = TRUE)) %>%
		    	   dplyr::mutate(timepoint_days_since_start_of_RT = case_when(
				   		timepoint_days_since_start_of_RT > 200 ~ 200,
				   		TRUE ~ timepoint_days_since_start_of_RT
			   )),
		    mapping = aes(x = timepoint_days_since_start_of_RT, y = `Pr(x=1)`, shape = Is_ctDNA),
		    width = 0, height = 0, inherit.aes = FALSE, color = "#e41a1c") +
	scale_shape_manual(values = c("-ve" = 1, "+ve" = 2)) +
	annotate(geom = "rect", xmin = -50, xmax = 0, ymin = 0, ymax = 1, fill = "#feb24c", color = NA, alpha = .25) +
	scale_x_continuous() +
	scale_y_continuous(limits = c(0, 1)) +
	xlab("Time since start of RT (Days)") +
	ylab("Posterior Probability") +
	facet_wrap(~patient_id_mskcc) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(shape = guide_legend(title = "ctDNA", order = 1))
	
pdf("../res/Time_Point_by_Posterior_Probability_Baseline_MRD_Negative.pdf", width = 15/2, height = 10/2)
print(plot_)
dev.off()

plot_ = posterior_smry %>%
	dplyr::left_join(mrd_smry %>%
			 dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			 dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
	dplyr::mutate(Is_ctDNA = case_when(
		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
		`MRD-Landmark_Result` == "ABSENT" ~ "-ve"
	)) %>%
	dplyr::left_join(Baseline_ctDNA_Neg, by = "patient_id_mskcc") %>%
	dplyr::filter(Is_Baseline) %>%
	dplyr::select(patient_id_mskcc, timepoint_days_since_start_of_RT, mean_af, max_af, median_af, Is_ctDNA) %>%
	reshape2::melt(id.vars = c("patient_id_mskcc", "timepoint_days_since_start_of_RT", "Is_ctDNA"), measure.vars = c("mean_af", "max_af", "median_af"),
		       variable.name = "summary", value.name = "AF") %>%
	dplyr::arrange(as.numeric(gsub("CTMS-", "", patient_id_mskcc))) %>%
	dplyr::mutate(patient_id_mskcc = factor(patient_id_mskcc, levels = unique(patient_id_mskcc), ordered = TRUE)) %>%
	dplyr::mutate(timepoint_days_since_start_of_RT = case_when(
		timepoint_days_since_start_of_RT > 200 ~ 200,
		TRUE ~ timepoint_days_since_start_of_RT
	)) %>%
	dplyr::mutate(AF = ifelse(AF>0.05, 0.05, AF)) %>%
	dplyr::mutate(AF = ifelse(AF<0.00005, 0.00005, AF)) %>%
	dplyr::mutate(summary = case_when(
		summary == "mean_af" ~ "Mean",
		summary == "median_af" ~ "Median",
		summary == "max_af" ~ "Max"
	)) %>%
	ggplot(aes(x = timepoint_days_since_start_of_RT, y = (100*AF), color = summary)) +
	geom_line(stat = "identity", size = .5, alpha = 1) +
	geom_jitter(data = posterior_smry %>%
		    	  dplyr::left_join(mrd_smry %>%
			 		   dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
					   dplyr::select(sample_name, `MRD-Landmark_Result`), by = "sample_name") %>%
		    	  dplyr::mutate(Is_ctDNA = case_when(
				  		`MRD-Landmark_Result` == "PRESENT" ~ "+ve",
				  		`MRD-Landmark_Result` == "ABSENT" ~ "-ve"
			  )) %>%
		    	  dplyr::left_join(Baseline_ctDNA_Neg, by = "patient_id_mskcc") %>%
		    	  dplyr::filter(Is_Baseline) %>%
		    	  dplyr::select(patient_id_mskcc, timepoint_days_since_start_of_RT, mean_af, max_af, median_af, Is_ctDNA) %>%
		    	  reshape2::melt(id.vars = c("patient_id_mskcc", "timepoint_days_since_start_of_RT", "Is_ctDNA"), measure.vars = c("mean_af", "max_af", "median_af"),
					 variable.name = "summary", value.name = "AF") %>%
		    	  dplyr::arrange(as.numeric(gsub("CTMS-", "", patient_id_mskcc))) %>%
		    	  dplyr::mutate(patient_id_mskcc = factor(patient_id_mskcc, levels = unique(patient_id_mskcc), ordered = TRUE)) %>%
		    	  dplyr::mutate(timepoint_days_since_start_of_RT = case_when(
						timepoint_days_since_start_of_RT > 200 ~ 200,
				  		TRUE ~ timepoint_days_since_start_of_RT
			  )) %>%
		    	  dplyr::mutate(AF = ifelse(AF>0.05, 0.05, AF)) %>%
			  dplyr::mutate(AF = ifelse(AF<0.00005, 0.00005, AF)) %>%
		    	  dplyr::mutate(summary = case_when(
				  		summary == "mean_af" ~ "Mean",
				  		summary == "median_af" ~ "Median",
				  		summary == "max_af" ~ "Max"
			  )),
		    mapping = aes(x = timepoint_days_since_start_of_RT, y = (100*AF), color = summary, shape = Is_ctDNA),
		    width = 0, height = 0, inherit.aes = FALSE) +
	scale_shape_manual(values = c("-ve" = 1, "+ve" = 2)) +
	scale_color_brewer(type = "qual", palette = 6) +
	annotate(geom = "rect", xmin = -50, xmax = 0, ymin = 0.001, ymax = 10, fill = "#feb24c", color = NA, alpha = .25) +
	scale_x_continuous() +
	scale_y_log10(limits = c(0.001, 10),
		      breaks = c(1E-3, 1E-2, 1E-1, 1, 10),
		      labels = c("0.001", "0.01", "0.1", "1", "10")) +
	xlab("Time since start of RT (Days)") +
	ylab("Allele Fraction (%)") +
	facet_wrap(~patient_id_mskcc) +
	theme_minimal() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20))) +
	guides(color = guide_legend(title = "AF", override.aes = list(shape = NA), order = 1),
	       shape = guide_legend(title = "ctDNA", order = 2))
	
pdf("../res/Time_Point_by_AF_Baseline_MRD_Negative.pdf", width = 15/2, height = 10/2)
print(plot_)
dev.off()
