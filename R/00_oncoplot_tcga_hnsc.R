#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")


Hugo_Symbols = c("TP53", "FAT1", "CDKN2A", "PIK3CA",
		 "NOTCH1", "KMT2D", "NSD1", "CASP8",
		 "NFE2L2", "AJUBA", "FBXW7", "TGFBR2",
		 "HRAS", "CUL3", "RB1", "PTEN", "PIK3R1")

vc_cols = RColorBrewer::brewer.pal(n = 8, name = "Set1")
names(vc_cols) = c('Frame_Shift_Del', 'Missense_Mutation',
		   'Nonsense_Mutation', 'Multi_Hit',
		   'Frame_Shift_Ins', 'In_Frame_Ins',
		   'Splice_Site', 'In_Frame_Del')

#==================================================
# TCGA oncoprint
#==================================================
TCGA = TCGAmutations::tcga_load(study = "HNSC", source = "MC3")

TCGA@data = TCGA@data %>%
	    dplyr::select(which(!duplicated(colnames(.))))

TCGA@clinical.data = TCGA@clinical.data %>%
		     dplyr::left_join(readr::read_tsv(file = url_tcga_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
				      readr::type_convert() %>%
				      dplyr::rename(Tumor_Sample_Barcode_min = Barcode,
						    HPV_Status = Final_HPV_Status) %>%
				      dplyr::mutate(is_selected = 1), by = "Tumor_Sample_Barcode_min")

TCGA = subsetMaf(TCGA, clinQuery = "!is.na(is_selected)")

color_palette_hpv = c("#7bccc4", "#0868ac")
names(color_palette_hpv) = unique(TCGA@clinical.data$HPV_Status)

pdf(file = "../res/OncoPrint_Cancer_Genes_TCGA_HPV_Negative_Positive.pdf", width = 4.75, height = 4)
oncoplot(maf = TCGA,
	 genes = Hugo_Symbols,
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 sortByMutation = TRUE,
	 sortByAnnotation = TRUE,
	 colors = vc_cols,
	 bgCol = "#EEEDEB",
	 drawRowBar = FALSE,
	 drawColBar = FALSE,
	 showTitle = FALSE,
	 gene_mar = 7,
	 showTumorSampleBarcodes = FALSE,
	 clinicalFeatures = c("HPV_Status"),
	 cohortSize = nrow(TCGA@clinical.data),
	 removeNonMutated = FALSE,
	 fill = TRUE,
	 annotationColor = list(HPV_Status = color_palette_hpv),
	 anno_height = .75)
dev.off()

#==================================================
# 30Gy oncoprint
#==================================================
clinical = readr::read_tsv(file = url_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::rename(Tumor_Sample_Barcode = tumor_id_wes,
			 T_Stage = t_stage,
			 N_Stage = n_stage,
			 Smoking_Status = smoking_category_yes_never) %>%
	   tidyr::drop_na(Tumor_Sample_Barcode) %>%
	   dplyr::select(-additional_crt_post_neck_dissection) %>%
	   dplyr::mutate(Final_HPV_Status = "Positive") %>%
	   dplyr::mutate(is_selected = 1)

mutation_smry = maftools::read.maf(maf = url_tumor_variants, clinical = clinical)

mutation_smry = subsetMaf(mutation_smry, clinQuery = "!is.na(is_selected)")

color_palette_t = c("#f0f9e8", "#bae4bc", "#2b8cbe")
names(color_palette_t) = unique(mutation_smry@clinical.data$T_Stage)

color_palette_n = c("#f0f9e8", "#bae4bc", "#7bccc4", "#2b8cbe")
names(color_palette_n) = unique(mutation_smry@clinical.data$N_Stage)

color_palette_s = c("#f0f9e8", "#bae4bc")
names(color_palette_s) = unique(mutation_smry@clinical.data$Smoking_Status)

pdf(file = "../res/OncoPrint_Cancer_Genes_30Gy_HPV_Positive.pdf", width = 2.75, height = 4.5)
oncoplot(maf = mutation_smry,
	 genes = Hugo_Symbols,
	 keepGeneOrder = TRUE,
	 GeneOrderSort = FALSE,
	 sortByMutation = FALSE,
	 sortByAnnotation = FALSE,
	 colors = vc_cols,
	 bgCol = "#EEEDEB",
	 drawRowBar = FALSE,
	 drawColBar = FALSE,
	 logColBar = FALSE,
	 showTitle = FALSE,
	 gene_mar = 7,
	 showTumorSampleBarcodes = FALSE,
	 clinicalFeatures = c("T_Stage", "N_Stage", "Smoking_Status"),
	 cohortSize = nrow(mutation_smry@clinical.data),
	 removeNonMutated = FALSE,
	 fill = TRUE,
	 annotationColor = list(T_Stage = color_palette_t,
			        N_Stage = color_palette_n,
			        Smoking_Status = color_palette_s),
	anno_height = 1.5)
dev.off()

#==================================================
# TCGA fraction of patients
#==================================================
TCGA = TCGAmutations::tcga_load(study = "HNSC", source = "MC3")
TCGA@data = TCGA@data %>%
	    dplyr::select(which(!duplicated(colnames(.))))
TCGA@clinical.data = TCGA@clinical.data %>%
		     dplyr::left_join(readr::read_tsv(file = url_tcga_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
				      readr::type_convert() %>%
				      dplyr::rename(Tumor_Sample_Barcode_min = Barcode) %>%
				      dplyr::mutate(is_selected = 1), by = "Tumor_Sample_Barcode_min")
HPV_Neg = subsetMaf(TCGA, clinQuery = "!is.na(is_selected) & Final_HPV_Status=='Negative'")
HPV_Pos = subsetMaf(TCGA, clinQuery = "!is.na(is_selected) & Final_HPV_Status=='Positive'")

# FoundationOne
Percent_Neg = HPV_Neg@data %>%
	      dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
	      dplyr::summarize(n = 1) %>%
	      dplyr::ungroup() %>%
	      dplyr::left_join(readr::read_tsv(file = url_gene_list[["foundation_one"]], col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(Is_Present = 1), by = "Hugo_Symbol") %>%
	      dplyr::filter(!is.na(Is_Present)) %>%
	      dplyr::group_by(Tumor_Sample_Barcode) %>%
	      dplyr::summarize(n = n()) %>%
	      {table(.$n)} %>%
	      tibble::enframe(name = "number_variants", value = "number_samples") %>%
	      dplyr::mutate(fraction_samples = 0)

Percent_Pos = HPV_Pos@data %>%
	      dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
	      dplyr::summarize(n = 1) %>%
	      dplyr::ungroup() %>%
	      dplyr::left_join(readr::read_tsv(file = url_gene_list[["foundation_one"]], col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(Is_Present = 1), by = "Hugo_Symbol") %>%
	      dplyr::filter(!is.na(Is_Present)) %>%
	      dplyr::group_by(Tumor_Sample_Barcode) %>%
	      dplyr::summarize(n = n()) %>%
	      {table(.$n)} %>%
	      tibble::enframe(name = "number_variants", value = "number_samples") %>%
	      dplyr::mutate(fraction_samples = 0)

for (i in 1:nrow(Percent_Neg)) {
	Percent_Neg$fraction_samples[i] = sum(Percent_Neg$number_samples[i:nrow(Percent_Neg)])/241 * 100
}

for (i in 1:nrow(Percent_Pos)) {
	Percent_Pos$fraction_samples[i] = sum(Percent_Pos$number_samples[i:nrow(Percent_Pos)])/36 * 100
}

plot_foundation = Percent_Neg %>%
		  dplyr::mutate(HPV = "-ve") %>%
		  dplyr::bind_rows(Percent_Pos %>%
				   dplyr::mutate(HPV = "+ve")) %>%
		  dplyr::mutate(number_variants = as.numeric(number_variants)) %>%
		  dplyr::filter(number_variants<10) %>%
		  dplyr::mutate(Assay = "FoudationOne CDx")
		
# Guardant360
Percent_Neg = HPV_Neg@data %>%
	      dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
	      dplyr::summarize(n = 1) %>%
	      dplyr::ungroup() %>%
	      dplyr::left_join(readr::read_tsv(file = url_gene_list[["guardant_360"]], col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(Is_Present = 1), by = "Hugo_Symbol") %>%
	      dplyr::filter(!is.na(Is_Present)) %>%
	      dplyr::group_by(Tumor_Sample_Barcode) %>%
	      dplyr::summarize(n = n()) %>%
	      {table(.$n)} %>%
	      tibble::enframe(name = "number_variants", value = "number_samples") %>%
	      dplyr::mutate(fraction_samples = 0)

Percent_Pos = HPV_Pos@data %>%
	      dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
	      dplyr::summarize(n = 1) %>%
	      dplyr::ungroup() %>%
	      dplyr::left_join(readr::read_tsv(file = url_gene_list[["guardant_360"]], col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(Is_Present = 1), by = "Hugo_Symbol") %>%
	      dplyr::filter(!is.na(Is_Present)) %>%
	      dplyr::group_by(Tumor_Sample_Barcode) %>%
	      dplyr::summarize(n = n()) %>%
	      {table(.$n)} %>%
	      tibble::enframe(name = "number_variants", value = "number_samples") %>%
	      dplyr::mutate(fraction_samples = 0)

for (i in 1:nrow(Percent_Neg)) {
	Percent_Neg$fraction_samples[i] = sum(Percent_Neg$number_samples[i:nrow(Percent_Neg)])/241 * 100
}

for (i in 1:nrow(Percent_Pos)) {
	Percent_Pos$fraction_samples[i] = sum(Percent_Pos$number_samples[i:nrow(Percent_Pos)])/36 * 100
}

plot_guardant = Percent_Neg %>%
		dplyr::mutate(HPV = "-ve") %>%
		dplyr::bind_rows(Percent_Pos %>%
				 dplyr::mutate(HPV = "+ve")) %>%
		dplyr::mutate(number_variants = as.numeric(number_variants)) %>%
		dplyr::filter(number_variants<10) %>%
		dplyr::mutate(Assay = "Guardant360 CDx")

# MSK-ACCESS v1
Percent_Neg = HPV_Neg@data %>%
	      dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
	      dplyr::summarize(n = 1) %>%
	      dplyr::ungroup() %>%
	      dplyr::left_join(readr::read_tsv(file = url_gene_list[["msk_accessv1"]], col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(Is_Present = 1), by = "Hugo_Symbol") %>%
	      dplyr::filter(!is.na(Is_Present)) %>%
	      dplyr::group_by(Tumor_Sample_Barcode) %>%
	      dplyr::summarize(n = n()) %>%
	      {table(.$n)} %>%
	      tibble::enframe(name = "number_variants", value = "number_samples") %>%
	      dplyr::mutate(fraction_samples = 0)

Percent_Pos = HPV_Pos@data %>%
	      dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
	      dplyr::summarize(n = 1) %>%
	      dplyr::ungroup() %>%
	      dplyr::left_join(readr::read_tsv(file = url_gene_list[["msk_accessv1"]], col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(Is_Present = 1), by = "Hugo_Symbol") %>%
	      dplyr::filter(!is.na(Is_Present)) %>%
	      dplyr::group_by(Tumor_Sample_Barcode) %>%
	      dplyr::summarize(n = n()) %>%
	      {table(.$n)} %>%
	      tibble::enframe(name = "number_variants", value = "number_samples") %>%
	      dplyr::mutate(fraction_samples = 0)

for (i in 1:nrow(Percent_Neg)) {
	Percent_Neg$fraction_samples[i] = sum(Percent_Neg$number_samples[i:nrow(Percent_Neg)])/241 * 100
}

for (i in 1:nrow(Percent_Pos)) {
	Percent_Pos$fraction_samples[i] = sum(Percent_Pos$number_samples[i:nrow(Percent_Pos)])/36 * 100
}

plot_ = Percent_Neg %>%
	dplyr::mutate(HPV = "-ve") %>%
	dplyr::bind_rows(Percent_Pos %>%
			 dplyr::mutate(HPV = "+ve")) %>%
	dplyr::mutate(number_variants = as.numeric(number_variants)) %>%
	dplyr::filter(number_variants<10) %>%
	dplyr::mutate(Assay = "MSK-ACCESS") %>%
	dplyr::bind_rows(plot_foundation) %>%
	dplyr::bind_rows(plot_guardant) %>%
	ggplot(aes(x = number_variants, y = fraction_samples, color = Assay, linetype = HPV, shape = HPV)) +
	geom_line(stat = "identity", size = .75) +
	geom_vline(xintercept = 3, linetype = 2, size = .5, color = "goldenrod3") +
	geom_point(stat = "identity", fill = "white", size = 2.5) +
	scale_shape_manual(values = c(21, NA)) +
	scale_color_brewer(type = "qual", palette = 6) +
	scale_linetype_manual(values = c(1, 3)) +
	xlab("Number of Variants") +
	ylab("Fraction of Patients (%)") +
	scale_x_continuous(limits = c(1, 10),
			   breaks = 1:10,
			   labels = 1:10) +
	scale_y_continuous() +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12)) +
	guides(color = guide_legend(override.aes = list(shape = 21, fill = "white")))

pdf(file = "../res/Fraction_Samples_by_Variant_TCGA.pdf", width = 5, height = 3)
print(plot_)
dev.off()

#==================================================
# 30Gy fraction of patients
#==================================================
clinical = readr::read_tsv(file = url_clinical, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::rename(Tumor_Sample_Barcode = tumor_id_wes) %>%
	   tidyr::drop_na(Tumor_Sample_Barcode) %>%
	   dplyr::select(-additional_crt_post_neck_dissection) %>%
	   dplyr::mutate(Final_HPV_Status = "Positive") %>%
	   dplyr::mutate(is_selected = 1)
HPV_Pos = maftools::read.maf(maf = url_tumor_variants, clinical = clinical)
HPV_Pos = subsetMaf(HPV_Pos, clinQuery = "!is.na(is_selected)")

# FoundationOne
Percent_Pos = HPV_Pos@data %>%
	      dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
	      dplyr::summarize(n = 1) %>%
	      dplyr::ungroup() %>%
	      dplyr::left_join(readr::read_tsv(file = url_gene_list[["foundation_one"]], col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(Is_Present = 1), by = "Hugo_Symbol") %>%
	      dplyr::filter(!is.na(Is_Present)) %>%
	      dplyr::group_by(Tumor_Sample_Barcode) %>%
	      dplyr::summarize(n = n()) %>%
	      {table(.$n)} %>%
	      tibble::enframe(name = "number_variants", value = "number_samples") %>%
	      dplyr::mutate(fraction_samples = 0)

for (i in 1:nrow(Percent_Pos)) {
	Percent_Pos$fraction_samples[i] = sum(Percent_Pos$number_samples[i:nrow(Percent_Pos)])/123 * 100
}

plot_foundation = Percent_Pos %>%
		  dplyr::mutate(number_variants = as.numeric(number_variants)) %>%
		  dplyr::filter(number_variants<10) %>%
		  dplyr::mutate(Assay = "FoudationOne CDx")
		
# Guardant360
Percent_Pos = HPV_Pos@data %>%
	      dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
	      dplyr::summarize(n = 1) %>%
	      dplyr::ungroup() %>%
	      dplyr::left_join(readr::read_tsv(file = url_gene_list[["guardant_360"]], col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(Is_Present = 1), by = "Hugo_Symbol") %>%
	      dplyr::filter(!is.na(Is_Present)) %>%
	      dplyr::group_by(Tumor_Sample_Barcode) %>%
	      dplyr::summarize(n = n()) %>%
	      {table(.$n)} %>%
	      tibble::enframe(name = "number_variants", value = "number_samples") %>%
	      dplyr::mutate(fraction_samples = 0)

for (i in 1:nrow(Percent_Pos)) {
	Percent_Pos$fraction_samples[i] = sum(Percent_Pos$number_samples[i:nrow(Percent_Pos)])/123 * 100
}

plot_guardant = Percent_Pos %>%
		dplyr::mutate(number_variants = as.numeric(number_variants)) %>%
		dplyr::filter(number_variants<10) %>%
		dplyr::mutate(Assay = "Guardant360 CDx")

# MSK-ACCESS v1
Percent_Pos = HPV_Pos@data %>%
	      dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
	      dplyr::summarize(n = 1) %>%
	      dplyr::ungroup() %>%
	      dplyr::left_join(readr::read_tsv(file = url_gene_list[["msk_accessv1"]], col_names = TRUE, col_types = cols(.default = col_character())) %>%
			       readr::type_convert() %>%
			       dplyr::mutate(Is_Present = 1), by = "Hugo_Symbol") %>%
	      dplyr::filter(!is.na(Is_Present)) %>%
	      dplyr::group_by(Tumor_Sample_Barcode) %>%
	      dplyr::summarize(n = n()) %>%
	      {table(.$n)} %>%
	      tibble::enframe(name = "number_variants", value = "number_samples") %>%
	      dplyr::mutate(fraction_samples = 0)

for (i in 1:nrow(Percent_Pos)) {
	Percent_Pos$fraction_samples[i] = sum(Percent_Pos$number_samples[i:nrow(Percent_Pos)])/123 * 100
}

plot_ = Percent_Pos %>%
	dplyr::mutate(number_variants = as.numeric(number_variants)) %>%
	dplyr::filter(number_variants<10) %>%
	dplyr::mutate(Assay = "MSK-ACCESS") %>%
	dplyr::bind_rows(plot_foundation) %>%
	dplyr::bind_rows(plot_guardant) %>%
	ggplot(aes(x = number_variants, y = fraction_samples, color = Assay)) +
	geom_line(stat = "identity", size = .75) +
	geom_vline(xintercept = 3, linetype = 2, size = .5, color = "goldenrod3") +
	geom_point(stat = "identity", shape = 21, fill = "white", size = 2.5) +
	scale_color_brewer(type = "qual", palette = 6) +
	xlab("Number of Variants") +
	ylab("Fraction of Patients (%)") +
	scale_x_continuous(limits = c(1, 10),
			   breaks = 1:10,
			   labels = 1:10) +
	scale_y_continuous(limits = c(1, 99)) +
	theme_classic() +
	theme(axis.title.x = element_text(margin = margin(t = 20), size = 12),
	      axis.title.y = element_text(margin = margin(r = 20), size = 12),
	      axis.text.x = element_text(size = 12),
	      axis.text.y = element_text(size = 12))

pdf(file = "../res/Fraction_Samples_by_Variant_30Gy.pdf", width = 5, height = 3)
print(plot_)
dev.off()

