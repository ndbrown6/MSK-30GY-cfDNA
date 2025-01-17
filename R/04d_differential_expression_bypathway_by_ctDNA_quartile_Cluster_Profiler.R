#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))

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
	
}
