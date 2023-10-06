#==================================================
# David Brown
# brownd7@mskcc.org
#==================================================
rm(list=ls(all=TRUE))
source("config.R")

manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(!is.na(bam_file_name_hpv))

preanalytical_conditions = readr::read_tsv(file = url_preanalytical_conidtions, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			   readr::type_convert()

manifest = manifest %>%
	   dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	   dplyr::mutate(patient_id_mskcc = case_when(
		   is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
		   TRUE ~ patient_id_mskcc
	   ))

#####################################################
## hpv bam files
#####################################################
patient_uuid = manifest %>%
	       dplyr::filter(!is.na(bam_file_name_hpv)) %>%
	       dplyr::arrange(patient_id_mskcc) %>%
	       .[["patient_id_mskcc"]] %>%
	       unique()

is_file = file.create("../res/samples__HPV.yaml", showWarnings = TRUE)
if (is_file) {
	for (i in 1:length(patient_uuid)) {
		cat(paste0("- name: ", patient_uuid[i], "\n"), file = "../res/samples__HPV.yaml", append = TRUE)
		cat("  tumor: [", file = "../res/samples__HPV.yaml", append = TRUE)
		msk_id = manifest %>%
			 dplyr::filter(patient_id_mskcc == patient_uuid[i]) %>%
			 .[["sample_id_mskcc"]]
		invitae_id = manifest %>%
			     dplyr::filter(patient_id_mskcc == patient_uuid[i]) %>%
			     .[["sample_id_invitae"]]
		for (j in 1:length(msk_id)) {
			cat(paste0(msk_id[j], "-", invitae_id[j]), file = "../res/samples__HPV.yaml", append = TRUE)
			if (j == length(msk_id)) {
				cat("]\n", file = "../res/samples__HPV.yaml", append = TRUE)
			} else {
				cat(",", file = "../res/samples__HPV.yaml", append = TRUE)
			}
		}
	}
}

#####################################################
## mrd vcf files
#####################################################
patient_uuid = manifest %>%
	       dplyr::filter(!is.na(vcf_file_name)) %>%
	       dplyr::arrange(patient_id_mskcc) %>%
	       .[["patient_id_mskcc"]] %>%
	       unique()

is_file = file.create("../res/samples__MRD.yaml", showWarnings = TRUE)
if (is_file) {
	for (i in 1:length(patient_uuid)) {
		cat(paste0("- name: ", patient_uuid[i], "\n"), file = "../res/samples__MRD.yaml", append = TRUE)
		cat("  tumor: [", file = "../res/samples__MRD.yaml", append = TRUE)
		msk_id = manifest %>%
			 dplyr::filter(patient_id_mskcc == patient_uuid[i]) %>%
			 .[["sample_id_mskcc"]]
		invitae_id = manifest %>%
			     dplyr::filter(patient_id_mskcc == patient_uuid[i]) %>%
			     .[["sample_id_invitae"]]
		for (j in 1:length(msk_id)) {
			cat(paste0(msk_id[j], "-", invitae_id[j]), file = "../res/samples__MRD.yaml", append = TRUE)
			if (j == length(msk_id)) {
				cat("]\n", file = "../res/samples__MRD.yaml", append = TRUE)
			} else {
				cat(",", file = "../res/samples__MRD.yaml", append = TRUE)
			}
		}
	}
}

#########################################################
# (+) Samples with mean AF > 1% from MRD assay
# (+) Samples with HPV subtype in assay
# (+) Samples with +ve MRD assay
# (-) Samples with mean AF = 0% from MRD assay
# (-) Samples with max AF = 0% from MRD assay
# (-) Patients with no nodal dissection ≥ 2 years
# (-) No duplicate patients
# (+) Samples with -ve MRD assay
#########################################################
manifest = readr::read_tsv(file = url_manifest, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::filter(!is.na(bam_file_name_hpv)) %>%
	   dplyr::mutate(sample_uuid = paste0(sample_id_mskcc, "-", sample_id_invitae))

preanalytical_conditions = readr::read_tsv(file = url_preanalytical_conidtions, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			   readr::type_convert()

mrd_smry = readr::read_tsv(file = url_mrd_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	   readr::type_convert() %>%
	   dplyr::mutate(sample_uuid = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer))

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

mutation_smry = readr::read_tsv(file = url_mutation_summary, col_names = TRUE, col_types = cols(.default = col_character())) %>%
		readr::type_convert()

nodal_dissection_smry = readr::read_tsv(file = url_no_node_dissection, col_names = TRUE, col_types = cols(.default = col_character())) %>%
			readr::type_convert() %>%
			dplyr::rename(patient_id_mskcc = sample)

aln_metrics = readr::read_tsv(file = url_aln_metrics, col_names = TRUE, col_types = cols(.default = col_character())) %>%
	      readr::type_convert() %>%
	      dplyr::filter(CATEGORY == "PAIR") %>%
	      dplyr::mutate(TOTAL_READS = TOTAL_READS/2)

manifest = manifest %>%
	   dplyr::left_join(preanalytical_conditions, by = "sample_id_mskcc") %>%
	   dplyr::mutate(patient_id_mskcc = case_when(
		   is.na(patient_id_mskcc) & sample_id_mskcc=="21-144-03654" ~ "CTMS-164",
		   TRUE ~ patient_id_mskcc
	   )) %>%
	   dplyr::mutate(sample_name = paste0(sample_id_mskcc, "-", sample_id_invitae)) %>%
	   dplyr::left_join(hpv_smry, by = "patient_id_mskcc")

smry_t_pos = mutation_smry %>%
	     dplyr::filter(FILTER == "PASS") %>%
	     dplyr::group_by(Tumor_Sample_Barcode) %>%
	     dplyr::summarize(mean_af = mean(t_maf),
			      max_af = max(t_maf)) %>%
	     dplyr::ungroup() %>%
	     dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
	     ## Mean AF > 1%
	     dplyr::filter(mean_af > (1/100)) %>%
	     dplyr::left_join(manifest, by = "sample_name") %>%
	     ## Known HPV
	     dplyr::filter(hpv_type_wes_wgs %in% names(target_contigs)) %>%
	     ## Positive MRD assay
	     dplyr::left_join(mrd_smry %>%
			      dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			      dplyr::select(sample_name, `MRD-Landmark_Result`, `MRD-Monitoring_MRD_Quantification`), by = "sample_name") %>%
	     dplyr::filter(`MRD-Landmark_Result` == "PRESENT") %>%
	     dplyr::mutate(Is_ctDNA = "+ve")

smry_t_neg = mutation_smry %>%
	     dplyr::filter(FILTER == "PASS") %>%
	     dplyr::group_by(Tumor_Sample_Barcode) %>%
	     dplyr::summarize(mean_af = mean(t_maf),
			      max_af = max(t_maf)) %>%
	     dplyr::ungroup() %>%
	     dplyr::rename(sample_name = Tumor_Sample_Barcode) %>%
	     ## Mean AF == 0
	     dplyr::filter(mean_af == 0) %>%
	     ## Max AF == 0
	     dplyr::filter(max_af == 0) %>%
	     dplyr::left_join(manifest, by = "sample_name") %>%
	     dplyr::left_join(nodal_dissection_smry, by = "patient_id_mskcc") %>%
	     ## No nodal dissection ≥ 0 years
	     dplyr::filter(nd_event == 0 & timepoint_days_since_start_of_RT > 0) %>%
	     ## Negative MRD assay
	     dplyr::left_join(mrd_smry %>%
			      dplyr::mutate(sample_name = gsub(pattern = "EP-D1-D1", replacement = "", x = Sample_ID_Archer, fixed = TRUE)) %>%
			      dplyr::select(sample_name, `MRD-Landmark_Result`, `MRD-Monitoring_MRD_Quantification`), by = "sample_name") %>%
	     dplyr::filter(`MRD-Landmark_Result` == "ABSENT") %>%
	     dplyr::mutate(Is_ctDNA = "-ve")

smry_ft = dplyr::bind_rows(smry_t_pos %>%
			   dplyr::select(all_of(intersect(colnames(smry_t_pos), colnames(smry_t_neg)))),
			   smry_t_neg %>%
			   dplyr::select(all_of(intersect(colnames(smry_t_pos), colnames(smry_t_neg))))) %>%
	  dplyr::filter(hpv_type_wes_wgs == "HPV-16")

smry_ft %>%
dplyr::left_join(aln_metrics %>%
		 dplyr::rename(sample_name = SAMPLE_NAME),
		 by = "sample_name") %>%
dplyr::arrange(Is_ctDNA, patient_name, mean_af) %>%
dplyr::select(patient_name, sample_name, mean_af, Is_ctDNA, TOTAL_READS) %>%
dplyr::mutate(mean_af = mean_af * 100) %>%
data.frame() %>%
pander::pander()

smry_t_pos = smry_ft %>%
	     dplyr::filter(sample_name == "21-020-04087-GER2107031")

smry_t_neg = smry_ft %>%
	     dplyr::filter(sample_name != "20-114-01597-GER2106900") %>%
	     dplyr::filter(sample_name != "18-239-01669-GER2213139") %>%
	     dplyr::filter(Is_ctDNA == "-ve") %>%
	     dplyr::arrange(patient_name)

is_file = file.create("../res/samples__LoD.yaml", showWarnings = TRUE)
if (is_file) {
	for (i in 1:nrow(smry_t_pos)) {
		cat(paste0("- name: ", smry_t_pos$patient_id_mskcc[i], "\n"), file = "../res/samples__LoD.yaml", append = TRUE)
		cat("  tumor: [", file = "../res/samples__LoD.yaml", append = TRUE)
		msk_id = smry_t_pos$sample_id_mskcc[i]
		invitae_id = smry_t_pos$sample_id_invitae[i]
		for (j in 1:length(msk_id)) {
			cat(paste0(msk_id[j], "-", invitae_id[j]), file = "../res/samples__LoD.yaml", append = TRUE)
			if (j == length(msk_id)) {
				cat("]\n", file = "../res/samples__LoD.yaml", append = TRUE)
			} else {
				cat(",", file = "../res/samples__LoD.yaml", append = TRUE)
			}
		}
	}
}


for (i in 1:nrow(smry_t_neg)) {
	cat(paste0("- name: ", smry_t_neg$patient_id_mskcc[i], "\n"), file = "../res/samples__LoD.yaml", append = TRUE)
	cat("  normal: ", file = "../res/samples__LoD.yaml", append = TRUE)
	msk_id = smry_t_neg$sample_id_mskcc[i]
	invitae_id = smry_t_neg$sample_id_invitae[i]
	for (j in 1:length(msk_id)) {
		cat(paste0(msk_id[j], "-", invitae_id[j]), file = "../res/samples__LoD.yaml", append = TRUE)
		if (j == length(msk_id)) {
			cat("\n", file = "../res/samples__LoD.yaml", append = TRUE)
		}
	}
}

smry_t_neg = smry_ft %>%
	     dplyr::left_join(aln_metrics %>%
		 	      dplyr::rename(sample_name = SAMPLE_NAME),
		 	      by = "sample_name") %>%
	     dplyr::filter(sample_name != "20-114-01597-GER2106900") %>%
	     dplyr::filter(sample_name != "18-239-01669-GER2213139") %>%
	     dplyr::filter(Is_ctDNA == "-ve") %>%
	     dplyr::filter(TOTAL_READS>=1000) %>%
	     dplyr::arrange(patient_name)

is_file = file.create("../res/samples__LoB.yaml", showWarnings = TRUE)
if (is_file) {
	for (i in 1:nrow(smry_t_neg)) {
		cat(paste0("- name: ", smry_t_neg$patient_id_mskcc[i], "\n"), file = "../res/samples__LoB.yaml", append = TRUE)
		cat("  normal: ", file = "../res/samples__LoB.yaml", append = TRUE)
		msk_id = smry_t_neg$sample_id_mskcc[i]
		invitae_id = smry_t_neg$sample_id_invitae[i]
		for (j in 1:length(msk_id)) {
			cat(paste0(msk_id[j], "-", invitae_id[j]), file = "../res/samples__LoB.yaml", append = TRUE)
			if (j == length(msk_id)) {
				cat("\n", file = "../res/samples__LoB.yaml", append = TRUE)
			}
		}
	}
}
