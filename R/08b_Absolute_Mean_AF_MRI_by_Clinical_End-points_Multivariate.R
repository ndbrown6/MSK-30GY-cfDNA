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
			 dplyr::rename(abs_HPV_wk0 = `Pre-treatment`,
				       abs_HPV_wk1 = wk1,
				       abs_HPV_wk2 = wk2,
				       abs_HPV_wk3 = wk3),
			 by = "patient_id_mskcc") %>%
	dplyr::full_join(clinical, by = "patient_id_mskcc") %>%
	dplyr::select(patient_id_mskcc,
		      composite_end_point,
		      abs_AF_wk0, abs_AF_wk1, abs_AF_wk2,
		      #abs_AF_wk3,
		      #abs_HPV_wk0, abs_HPV_wk1, abs_HPV_wk2, abs_HPV_wk3,
		      abs_MRI_wk0 = MRI_rawdata_wk0,
		      abs_MRI_wk1 = MRI_rawdata_wk2,
		      abs_MRI_wk2 = MRI_rawdata_wk3) %>%
		      #abs_MRI_wk3 = MRI_rawdata_wk4) %>%
	dplyr::mutate(composite_end_point = case_when(
		composite_end_point ~ 1,
		!composite_end_point ~ 0
	)) %>%
	readr::type_convert() %>%
	dplyr::select(-patient_id_mskcc)

params = rpart::rpart.control(minsplit = 20,
			      minbucket = floor(20/3),
			      maxdepth = 4,
			      maxcompete = 1000,
			      maxsurrogate = 100,
			      usesurrogate = 2,
			      cp = 0,
			      xval = 100)
tree_all = rpart::rpart(composite_end_point ~ ., data = data_, method = "class", control = params)

tree_pruned = snip.rpart(x = tree_all, toss = 12)

pdf(file = "../res/Decesion_Tree_Initial_Rpart.pdf", width = 7, height = 7)
rpart.plot(tree_pruned, extra = 2)
dev.off()

tree_party = partykit::as.party(tree_pruned)

dens_df <- data.frame(x_dens = numeric(), y_dens = numeric(), id = numeric(), breaks = character())

id = 1
x_dens <- density(tree_party[id]$data$abs_AF_wk2, na.rm=T)$x
y_dens <- density(tree_party[id]$data$abs_AF_wk2, na.rm=T)$y
breaks <- rep("0", length(x_dens))
breaks[x_dens >= 0.02731] <- "1"
dens_df <- rbind(dens_df, data.frame(x_dens, y_dens, id, breaks))

id = 3
x_dens <- density(tree_party[id]$data$abs_MRI_wk1, na.rm=T)$x
y_dens <- density(tree_party[id]$data$abs_MRI_wk1, na.rm=T)$y/1000
breaks <- rep("0", length(x_dens))
breaks[x_dens >= 27472] <- "1"
dens_df <- rbind(dens_df, data.frame(x_dens, y_dens, id, breaks))

id = 4
x_dens <- density(tree_party[id]$data$abs_AF_wk1, na.rm=T)$x
y_dens <- density(tree_party[id]$data$abs_AF_wk1, na.rm=T)$y
breaks <- rep("0", length(x_dens))
breaks[x_dens >= 0.17969] <- "1"
dens_df <- rbind(dens_df, data.frame(x_dens, y_dens, id, breaks))


plot_ = ggparty(tree_party) +
	geom_edge() +
	geom_edge_label() +
	geom_node_plot(ids = 1,
		       gglist = list(geom_line(data = dens_df %>% dplyr::filter(id==1),
					       aes(x = x_dens, y = y_dens),
					       show.legend = FALSE,
					       alpha = 0.8),
				     geom_ribbon(data = dens_df %>% dplyr::filter(id==1),
						 aes(x = x_dens, ymin = 0, ymax = y_dens, fill = breaks),
						 show.legend = FALSE,
						 alpha = 0.8),
				     scale_x_continuous(labels = scientific_10),
				     scale_y_continuous(),
				     scale_fill_brewer(type = "qual", palette = 7),
				     xlab(""),
				     theme_classic(),
				     theme(axis.title.y = element_blank())),
		       size = 1.5,
		       height = 0.75) +
	geom_node_plot(ids = 3,
		       gglist = list(geom_line(data = dens_df %>% dplyr::filter(id==3),
					       aes(x = x_dens, y = y_dens),
					       show.legend = FALSE,
					       alpha = 0.8),
				     geom_ribbon(data = dens_df %>% dplyr::filter(id==3),
						 aes(x = x_dens, ymin = 0, ymax = y_dens, fill = breaks),
						 show.legend = FALSE,
						 alpha = 0.8),
				     scale_x_continuous(labels = scientific_10),
				     scale_y_continuous(),
				     scale_fill_brewer(type = "qual", palette = 7),
				     xlab(""),
				     theme_classic(),
				     theme(axis.title.y = element_blank())),
		       size = 1.5,
		       height = 0.75) +
	geom_node_plot(ids = 4,
		       gglist = list(geom_line(data = dens_df %>% dplyr::filter(id==4),
					       aes(x = x_dens, y = y_dens),
					       show.legend = FALSE,
					       alpha = 0.8),
				     geom_ribbon(data = dens_df %>% dplyr::filter(id==4),
						 aes(x = x_dens, ymin = 0, ymax = y_dens, fill = breaks),
						 show.legend = FALSE,
						 alpha = 0.8),
				     scale_x_continuous(labels = scientific_10),
				     scale_y_continuous(),
				     scale_fill_brewer(type = "qual", palette = 7),
				     xlab(""),
				     theme_classic(),
				     theme(axis.title.y = element_blank())),
		       size = 1.5,
		       height = 0.75) +
	geom_node_splitvar() +
	geom_node_plot(gglist = list(geom_bar(aes(x = "", fill = factor(composite_end_point)),
					      position = position_fill(), color = "black"),
				     scale_fill_brewer(type = "qual", palette = 7),
				     coord_polar("y", start = 0),
				     theme_void(),
				     guides(fill = guide_legend(title = "Risk group"))),
		       size = "log(nodesize)", shared_axis_labels = TRUE)

pdf(file = "../res/Decesion_Tree_Initial_PartyKit.pdf", width = 7, height = 7)
print(plot_)
dev.off()

dens_df <- data.frame(x_dens = numeric(), y_dens = numeric(), id = numeric(), breaks = character())

id = 1
x_dens <- density(log10(tree_party[id]$data$abs_AF_wk2), na.rm=T)$x
y_dens <- density(log10(tree_party[id]$data$abs_AF_wk2), na.rm=T)$y*10
breaks <- rep("0", length(x_dens))
breaks[x_dens >= log10(0.02731)] <- "1"
dens_df <- rbind(dens_df, data.frame(x_dens, y_dens, id, breaks))

id = 3
x_dens <- density(tree_party[id]$data$abs_MRI_wk1, na.rm=T)$x
y_dens <- density(tree_party[id]$data$abs_MRI_wk1, na.rm=T)$y*1E5
breaks <- rep("0", length(x_dens))
breaks[x_dens >= 27472] <- "1"
dens_df <- rbind(dens_df, data.frame(x_dens, y_dens, id, breaks))

id = 4
x_dens <- density(log10(tree_party[id]$data$abs_AF_wk1), na.rm=T)$x
y_dens <- density(log10(tree_party[id]$data$abs_AF_wk1), na.rm=T)$y*10
breaks <- rep("0", length(x_dens))
breaks[x_dens >= log10(0.17969)] <- "1"
dens_df <- rbind(dens_df, data.frame(x_dens, y_dens, id, breaks))

plot_ = dens_df %>%
	dplyr::mutate(x_dens = case_when(
		id == 3 ~ x_dens * 0.001,
		TRUE ~ x_dens
	)) %>%
	ggplot(aes(x = x_dens, y = y_dens)) +
	geom_line(stat = "identity", alpha = 0.8) +
	geom_ribbon(mapping = aes(x = x_dens, ymin = 0, ymax = y_dens, fill = breaks),
		    inherit.aes = TRUE, alpha = 0.8) +
	scale_x_continuous() +
	scale_y_continuous(limits = c(0, 5.5),
			   breaks = 0:5,
			   labels = 0:5)+
	scale_fill_brewer(type = "qual", palette = 7) +
	xlab("") +
	ylab("") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20)),
	      axis.title.y = element_text(margin = margin(r = 20)),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12),
	      strip.background = element_blank()) +
	guides(fill = FALSE) +
	facet_wrap(~id, scales = "free", ncol = 1)

pdf(file = "../res/Decesion_Tree_Initial_Densities.pdf", width = 3, height = 6)
print(plot_)
dev.off()

