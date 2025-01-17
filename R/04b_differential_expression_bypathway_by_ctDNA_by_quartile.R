#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

# HPV read counts upper versus lower quartile
posterior_probability = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::filter(timepoint_days_since_start_of_RT<0) %>%
			dplyr::select(patient_id_mskcc, read_count = `E6_7124:7601_READ_COUNT`) %>%
			dplyr::group_by(patient_id_mskcc) %>%
			dplyr::summarize(read_count = mean(read_count)) %>%
			dplyr::ungroup() %>%
			dplyr::mutate(Is_ctDNA = case_when(
				read_count >= quantile(read_count, .75) ~ "+ve",
				read_count < quantile(read_count, .25) ~ "-ve",
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
	       shape = guide_legend(title="MsigDB\nCollection", override.aes = list(size = 3), order = 2))

pdf(file = "../res/Differential_Expression_by_Pathway_ctDNA_outer_quartile.pdf", width = 5*1.25, height = 4*1.25)
print(plot_)
dev.off()

# ctDNA positive versus negative
posterior_probability = readr::read_tsv(file = "../res/Posterior_Probability_ALL.txt", col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::filter(timepoint_days_since_start_of_RT<0) %>%
			dplyr::select(patient_id_mskcc, Is_ctDNA) %>%
			dplyr::group_by(patient_id_mskcc) %>%
			dplyr::summarize(Is_ctDNA = case_when(
				any(Is_ctDNA=="+ve") ~ "+ve",
				all(Is_ctDNA=="-ve") ~ "-ve"
			)) %>%
			dplyr::ungroup() %>%
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
	       shape = guide_legend(title="MsigDB\nCollection", override.aes = list(size = 3), order = 2))

pdf(file = "../res/Differential_Expression_by_Pathway_ctDNA_postive_negative.pdf", width = 5*1.25, height = 4*1.25)
print(plot_)
dev.off()
