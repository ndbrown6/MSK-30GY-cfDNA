#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

posterior_probability = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::filter(timepoint_days_since_start_of_RT<0) %>%
			dplyr::select(patient_id_mskcc, read_count = `E6_7124:7601_READ_COUNT`) %>%
			dplyr::group_by(patient_id_mskcc) %>%
			dplyr::summarize(read_count = mean(read_count)) %>%
			dplyr::ungroup() %>%
			dplyr::mutate(Is_ctDNA = case_when(
				read_count >= quantile(read_count, .5) ~ "+ve",
				read_count < quantile(read_count, .5) ~ "-ve",
				TRUE ~ "NA"
			)) %>%
			readr::type_convert() %>%
			tidyr::drop_na()

tpm_by_gene = readr::read_tsv(file = url_tpm_by_gene, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert()

posterior_probability = posterior_probability %>%
			dplyr::filter(patient_id_mskcc %in% colnames(tpm_by_gene))

tpm_by_gene = tpm_by_gene %>%
	      dplyr::select(all_of(c("Hugo_Symbol", posterior_probability$patient_id_mskcc))) %>%
       	      dplyr::arrange(Hugo_Symbol) %>%
	      dplyr::mutate(Entrez_ID = as.character(mapIds(org.Hs.eg.db, Hugo_Symbol, 'ENTREZID', 'SYMBOL'))) %>%
	      dplyr::filter(!is.na(Entrez_ID))

expr_set = tpm_by_gene %>%
	   dplyr::select(posterior_probability %>% .[["patient_id_mskcc"]]) %>%
	   as.matrix()
rownames(expr_set) = tpm_by_gene %>% .[["Entrez_ID"]]

clin_data = posterior_probability %>%
	    dplyr::mutate(SAMPLE_ID = patient_id_mskcc) %>%
	    as.data.frame() %>%
	    tibble::column_to_rownames("SAMPLE_ID")
clin_data = new("AnnotatedDataFrame", data = clin_data)

expr_eset = ExpressionSet(assayData = expr_set,
			  phenoData = clin_data,
			  annotation = "org.Hs.eg.db")

ft_eset = nsFilter(expr_eset,
		   require.entrez = TRUE,
		   remove.dupEntrez = TRUE,
		   var.filter = TRUE)

BroadGeneSet = getGmt(con = "../data/gmt_by_entrez_id/c2.5.6.8.h.v2023.2.Hs.entrez.gmt")

BroadGeneSet_SSGSEA = gsva(expr = ft_eset$eset,
			   gset.idx.list = BroadGeneSet,
			   annotation = org.Hs.eg.db,
			   method = "ssgsea",
			   ssgsea.norm = FALSE,
			   min.sz = 10,
			   max.sz = 500,
			   mx.diff = FALSE,
			   parallel.sz = 8,
			   verbose = TRUE)

# ctDNA +ve | ctDNA -ve
response = ft_eset$eset$Is_ctDNA
status = ifelse(response=="+ve", 1, NA)
status = ifelse(response=="-ve", 0, status)
design = model.matrix(~ factor(status))
colnames(design) = c("intercept", "status")

fit = lmFit(object = exprs(BroadGeneSet_SSGSEA)[,!is.na(status),drop=FALSE], design = design)
fiteb = eBayes(fit, trend = FALSE)
setof_allsignatures = topTable(fiteb, coef = "status", number = Inf) %>%
		      dplyr::mutate(signature_name = rownames(.)) %>%
		      dplyr::as_tibble()

n = unlist(lapply(BroadGeneSet, function(x) {length(x@geneIds)}))
signature_name = names(BroadGeneSet)

setof_allsignatures = dplyr::left_join(setof_allsignatures,
				       dplyr::tibble(signature_name = signature_name,
						     n = n),
				       by = "signature_name") %>%
		      dplyr::mutate(category = case_when(
			      		   grepl("HALLMARK", signature_name) ~ "Hallmark gene sets",
			      		   grepl("KEGG_MEDICUS", signature_name) ~ "c2 KEGG Medicus",
			      		   grepl("REACTOME", signature_name) ~ "c2 Reactome",
			      		   grepl("GOBP", signature_name) ~ "c5 GO Biological Process",
			      		   grepl("GOCC", signature_name) ~ "c5 GO Cellular Component",
			      		   grepl("GOMF", signature_name) ~ "c5 GO Molecular Function",
			      		   TRUE ~ "Other (c6-c8)")) %>%
		      dplyr::mutate(logFC = logFC/20,
				    P.Value  = ifelse(P.Value<0.001, 0.001, P.Value))

plot_ = setof_allsignatures %>%
	dplyr::filter(grepl("Hallmark|Reactome|KEGG|GO Biological", category, fixed = FALSE, perl = TRUE)) %>%
	dplyr::mutate(is_apoptosis = case_when(
		grepl("APOP|NECRO|DEATH", signature_name, perl = TRUE, fixed = FALSE, ignore.case = TRUE) ~ "Yes",
		TRUE ~ "No"
	)) %>%
	dplyr::mutate(is_mitosis = case_when(
		grepl("MITOS|MITOT|G1|G2|CELL_CYCLE", signature_name, perl = TRUE, fixed = FALSE, ignore.case = TRUE) ~ "Yes",
		TRUE ~ "No"
	)) %>%
	dplyr::mutate(is_apoptosis = factor(is_apoptosis, levels = c("Yes", "No"), ordered = TRUE)) %>%
	dplyr::mutate(is_mitosis = factor(is_mitosis, levels = c("Yes", "No"), ordered = TRUE)) %>%
	dplyr::arrange(desc(is_apoptosis), desc(is_mitosis)) %>%
	dplyr::mutate(is_pathway = case_when(
		is_apoptosis == "Yes" ~ "Apoptosis",
		is_mitosis == "Yes" ~ "Mitosis",
		TRUE ~ "Other"
	)) %>%
	dplyr::mutate(is_pathway = factor(is_pathway, levels = c("Apoptosis", "Mitosis", "Other"), ordered = TRUE)) %>%
	ggplot(aes(x = logFC, y = (1/P.Value), fill = is_pathway, color = is_pathway, alpha = is_pathway, shape = category)) +
	geom_point(stat = "identity", size = 2) +
	geom_hline(yintercept=1/.1, color = "goldenrod3", linetype = 2, size = .55) +
	geom_vline(xintercept=c(-2.5, 2.5), color = "goldenrod3", linetype = 2, size = .55) +
	scale_fill_manual(values = c("Apoptosis" = "#fc9272",
				     "Mitosis" = "#1d91c0",
				     "Other" = "#ffffd9")) +
	scale_color_manual(values = c("Apoptosis" = "#000000",
				      "Mitosis" = "#000000",
				      "Other" = "#c6c1bb"),
			   guide = FALSE) +
	scale_alpha_manual(values = c("Apoptosis" = 1,
				      "Mitosis" = 1,
				      "Other" = .55),
			   guide = FALSE) +
	scale_shape_manual(values = c("Hallmark gene sets" = 21,
				      "c2 KEGG Medicus" = 22,
				      "c2 Reactome" = 23,
				      "c5 GO Biological Process" = 24,
				      "c5 GO Cellular Component" = 24,
				      "c5 GO Molecular Function" = 24,
				      "Other (c6-c8)" = 24)) +
	scale_x_continuous() +
	scale_y_log10(limits = c(1, 1000),
		      labels = log_10) +
	xlab(expression(Log[2]~"FC")) +
	ylab("P-value") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 14),
	      axis.title.y = element_text(margin = margin(r = 20), size = 14),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = guide_legend(title="Pathway\nRelatedness", override.aes = list(size = 3, shape = 21), order = 1),
	       shape = guide_legend(title="MsigDB\nCollection", override.aes = list(size = 3), order = 2)) +
	geom_text_repel(data = subset(setof_allsignatures %>%
				      dplyr::filter(grepl("Hallmark|Reactome|KEGG|GO Biological", category, fixed = FALSE, perl = TRUE)) %>%
				      dplyr::mutate(P.Value = ifelse(P.Value<(1e-5), 1e-5, P.Value)) %>%
				      dplyr::mutate(is_pathway = case_when(
					      grepl("APOP|NECRO|DEATH|MITOS|MITOT|G1|G2|CELL_CYCLE", signature_name, perl = TRUE, fixed = FALSE, ignore.case = TRUE) ~ "Yes",
					      TRUE ~ "No"
				      )) %>%
				      dplyr::filter(is_pathway=="Yes") %>%
				      dplyr::filter(grepl("Hallmark|Reactome", category, fixed = FALSE, perl = TRUE)),
				      (logFC>2.5|logFC<(-2.5)) & (1/P.Value)>(1/.1)),
			mapping = aes(x = logFC, y = (1/P.Value), label = signature_name),
			size = 2.5, segment.size = 0.35, force = 1, show.legend = FALSE, inherit.aes = FALSE)

pdf(file = "../res/Differential_Expression_by_Pathway_ctDNA.pdf", width = 5*1.25, height = 4*1.25)
print(plot_)
dev.off()

data_ = exprs(BroadGeneSet_SSGSEA) %>%
	as.data.frame() %>%
	tibble::rownames_to_column("signature_name") %>%
	dplyr::as_tibble() %>%
	dplyr::filter(signature_name %in% (setof_allsignatures %>%
					   dplyr::filter(P.Value<.1) %>%
					   dplyr::filter(grepl("APOP|NECRO|DEATH|MITOS|MITOT|G1|G2|CELL_CYCLE", signature_name, perl = TRUE, fixed = FALSE)) %>%
					   dplyr::filter(grepl("HALLMARK|REACTOME|KEGG", signature_name, perl = TRUE, fixed = FALSE)) %>%
					   .[["signature_name"]])) %>%
	reshape2::melt() %>%
	reshape2::dcast(variable ~ signature_name, value.var = "value") %>%
	dplyr::as_tibble() %>%
	dplyr::rename(patient_id_mskcc = variable) %>%
        dplyr::left_join(readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt",
					 col_names = TRUE, col_types = cols(.default = col_character())) %>%
			 readr::type_convert() %>%
			 dplyr::filter(timepoint_days_since_start_of_RT<0) %>%
			 dplyr::group_by(patient_id_mskcc) %>%
			 dplyr::summarize(Is_ctDNA = case_when(
				 			any(Is_ctDNA == "+ve") ~ "+ve",
				 			all(Is_ctDNA == "-ve") ~ "-ve"),
					  mean_af = mean(mean_af),
					  `E6_7124:7601_INSERT_SIZE` = mean(`E6_7124:7601_INSERT_SIZE`)) %>%
			 dplyr::ungroup() %>% 
			 dplyr::select(patient_id_mskcc, Is_ctDNA, mean_af,
				       'e6_insert_size' = `E6_7124:7601_INSERT_SIZE`),
			 by = "patient_id_mskcc") %>%
	dplyr::mutate(across(where(is.numeric), scale, center = TRUE, scale = TRUE)) %>%
	dplyr::arrange(Is_ctDNA, mean_af, e6_insert_size)

sample_order_ctdna_neg = data_ %>%
			 dplyr::filter(Is_ctDNA == "-ve") %>%
			 dplyr::select(-patient_id_mskcc, -Is_ctDNA, -mean_af, -e6_insert_size) %>%
			 dist() %>%
			 hclust(method = "ward.D2") %>%
			 .[["order"]]

sample_order_ctdna_pos = data_ %>%
			 dplyr::filter(Is_ctDNA == "+ve") %>%
			 dplyr::select(-patient_id_mskcc, -Is_ctDNA, -mean_af, -e6_insert_size) %>%
			 dist() %>%
			 hclust(method = "ward.D2") %>%
			 .[["order"]]

sample_order = dplyr::tibble(patient_id_mskcc = c((data_ %>%
						   dplyr::filter(Is_ctDNA == "-ve") %>%
						   .[["patient_id_mskcc"]])[sample_order_ctdna_neg],
						  (data_ %>%
						   dplyr::filter(Is_ctDNA == "+ve") %>%
						   .[["patient_id_mskcc"]])[sample_order_ctdna_pos])) %>%
	       dplyr::mutate(index = 1:n())

data_ = data_ %>%
	dplyr::left_join(sample_order, by = "patient_id_mskcc") %>%
	dplyr::arrange(index) %>%
	dplyr::select(-index)

color_palette = c("#f0f9e8", "#bae4bc", "#7bccc4", "#43a2ca", "#0868ac", "#045a8d")

color_palette_ctDNA = colorRampPalette(color_palette)(length(unique(data_$Is_ctDNA)))
names(color_palette_ctDNA) = unique(data_$Is_ctDNA)

color_palette_af = colorRampPalette(color_palette)(length(unique(data_$mean_af)))
names(color_palette_af) = sort(unique(data_$mean_af))

color_palette_e6 = colorRampPalette(color_palette)(length(unique(data_$e6_insert_size)))
names(color_palette_e6) = sort(unique(data_$e6_insert_size))

ha = HeatmapAnnotation(
	df = data_ %>%
	     dplyr::select(Is_ctDNA, mean_af, e6_insert_size) %>%
	     as.data.frame(),
	col = list(Is_ctDNA = color_palette_ctDNA,
		   mean_af = color_palette_af,
		   e6_insert_size = color_palette_e6),
	show_legend = FALSE,
	simple_anno_size = unit(3, "mm"),
	annotation_name_gp = gpar(fontsize = 9)
)

pdf(file = "../res/Heatmap_by_Pathway_ctDNA.pdf", width = 10, height = 4)
draw(Heatmap(matrix = data_ %>%
	     	      dplyr::select(-patient_id_mskcc, -Is_ctDNA, -mean_af, -e6_insert_size) %>%
	     	      t(),
	     col = viridis(n = 15, direction = -1),
	     name = " ",
	     na_col = "#f0f0f0",
	     rect_gp = gpar(col = "white", lwd = .2),
	     border = NA,

	     cluster_rows = TRUE,
	     show_row_names = TRUE,
	     row_names_side = "left",
	     row_names_gp = gpar(fontsize = 8),
	     show_row_dend = FALSE,
	     
	     cluster_columns = FALSE,
	     show_column_names = FALSE,
	     
	     
	     top_annotation = ha,

	     use_raster = FALSE,
	     show_heatmap_legend = TRUE,
	     heatmap_legend_param = list(legend_height = unit(.5, "cm"), legend_width = unit(1.5, "cm"))))
dev.off()

exprs(BroadGeneSet_SSGSEA) %>%
as.data.frame() %>%
tibble::rownames_to_column("signature_name") %>%
dplyr::as_tibble() %>%
dplyr::filter(signature_name %in% (setof_allsignatures %>%
				   dplyr::filter(P.Value<.1) %>%
				   .[["signature_name"]])) %>%
reshape2::melt() %>%
reshape2::dcast(variable ~ signature_name, value.var = "value") %>%
dplyr::as_tibble() %>%
dplyr::rename(patient_id_mskcc = variable) %>%
readr::write_tsv(path = "../res/GSEA.txt", col_names = TRUE, append = FALSE)
