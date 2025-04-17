#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))

'log_10' <- function(x) {
	parse(text=gsub("+", "", gsub("e", " %.% 10^-", scales::scientific_format()(x)), fixed=TRUE))
}


if ("clusterProfiler" %in% rownames(installed.packages())) {

	library('dplyr')
	library('readr')
	library('magrittr')
	library('clusterProfiler')
	library('AnnotationDbi')
	library('org.Hs.eg.db')
	library('genefilter')
	library('limma')
	library('ggplot2')
	library('ggrepel')

	# HPV read counts upper versus lower quartile
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

	tpm_by_gene = readr::read_tsv(file = "../data/tpm_by_gene.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
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

	response = ft_eset$eset$Is_ctDNA
	status = ifelse(response=="+ve", 1, NA)
	status = ifelse(response=="-ve", 0, status)
	design = model.matrix(~ factor(status))
	colnames(design) = c("intercept", "status")
	
	logFC = log2(apply(exprs(ft_eset$eset)[,design[,"status"]==1,drop=FALSE], 1, mean) /
		     apply(exprs(ft_eset$eset)[,design[,"status"]==0,drop=FALSE], 1, mean))
	setof_allgenes = sort(logFC, decreasing = TRUE)
	BroadGeneSet = read.gmt("../data/gmt_by_entrez_id/c2.5.6.8.h.v2023.2.Hs.entrez.gmt")
	setof_allsignatures = GSEA(geneList = setof_allgenes,
				   minGSSize = 10,
				   maxGSSize = 500,
				   eps = 1e-10,
				   pvalueCutoff = 1,
				   pAdjustMethod = "fdr",
				   gson = NULL,
				   TERM2GENE = BroadGeneSet,
				   verbose = TRUE,
				   seed = FALSE,
				   by = "fgsea") %>%
			      .@result %>%
			      dplyr::as_tibble() %>%
			      dplyr::select(signature_name = ID,
					    n = setSize,
					    NES = NES,
					    pvalue = pvalue,
					    p.adjust = p.adjust,
					    qvalue = qvalue)
	
	plot_ = setof_allsignatures %>%
	dplyr::mutate(category = case_when(
			      		   grepl("HALLMARK", signature_name) ~ "Hallmark gene sets",
			      		   grepl("KEGG_MEDICUS", signature_name) ~ "c2 KEGG Medicus",
			      		   grepl("REACTOME", signature_name) ~ "c2 Reactome",
			      		   grepl("GOBP", signature_name) ~ "c5 GO Biological Process",
			      		   grepl("GOCC", signature_name) ~ "c5 GO Cellular Component",
			      		   grepl("GOMF", signature_name) ~ "c5 GO Molecular Function",
			      		   TRUE ~ "Other (c6-c8)")) %>%
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
	ggplot(aes(x = NES, y = (1/p.adjust), fill = is_pathway, color = is_pathway, alpha = is_pathway, shape = category)) +
	geom_point(stat = "identity", size = 2) +
	geom_hline(yintercept=1/.1, color = "goldenrod3", linetype = 2, size = .55) +
	geom_vline(xintercept=c(-1.5, 1.5), color = "goldenrod3", linetype = 2, size = .55) +
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
	scale_y_log10(limits = c(1, NA),
		      labels = log_10) +
	xlab("NES") +
	ylab("FDR Adjusted P-value") +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 14),
	      axis.title.y = element_text(margin = margin(r = 20), size = 14),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(fill = guide_legend(title="Pathway\nRelatedness", override.aes = list(size = 3, shape = 21), order = 1),
	       shape = guide_legend(title="MsigDB\nCollection", override.aes = list(size = 3), order = 2)) +
	geom_text_repel(data = subset(setof_allsignatures %>%
				      dplyr::mutate(category = case_when(
			      		   grepl("HALLMARK", signature_name) ~ "Hallmark gene sets",
			      		   grepl("KEGG_MEDICUS", signature_name) ~ "c2 KEGG Medicus",
			      		   grepl("REACTOME", signature_name) ~ "c2 Reactome",
			      		   grepl("GOBP", signature_name) ~ "c5 GO Biological Process",
			      		   grepl("GOCC", signature_name) ~ "c5 GO Cellular Component",
			      		   grepl("GOMF", signature_name) ~ "c5 GO Molecular Function",
			      		   TRUE ~ "Other (c6-c8)")) %>%
				      dplyr::filter(grepl("Hallmark|Reactome|KEGG|GO Biological", category, fixed = FALSE, perl = TRUE)) %>%
				      dplyr::mutate(is_pathway = case_when(
					      grepl("APOP|NECRO|DEATH|MITOS|MITOT|G1|G2|CELL_CYCLE", signature_name, perl = TRUE, fixed = FALSE, ignore.case = TRUE) ~ "Yes",
					      TRUE ~ "No"
				      )) %>%
				      dplyr::filter(is_pathway=="Yes") %>%
				      dplyr::filter(grepl("Hallmark|Reactome", category, fixed = FALSE, perl = TRUE)),
				      (NES>1.5|NES<(-1.5)) & (1/p.adjust)>(1/.1)),
			mapping = aes(x = NES, y = (1/p.adjust), label = signature_name),
			size = 2.5, segment.size = 0.35, force = 1, show.legend = FALSE, inherit.aes = FALSE)
	
}
