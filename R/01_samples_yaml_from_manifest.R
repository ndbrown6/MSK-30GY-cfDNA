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

#==================================================
# HPV BAM files
#==================================================
patient_uuid = manifest %>%
	       dplyr::filter(!is.na(bam_file_name_hpv)) %>%
	       dplyr::arrange(patient_id_mskcc) %>%
	       .[["patient_id_mskcc"]] %>%
	       unique()

is_file = file.create("../res/HPV.yaml", showWarnings = TRUE)
if (is_file) {
	for (i in 1:length(patient_uuid)) {
		cat(paste0("- name: ", patient_uuid[i], "\n"), file = "../res/HPV.yaml", append = TRUE)
		cat("  tumor: [", file = "../res/HPV.yaml", append = TRUE)
		msk_id = manifest %>%
			 dplyr::filter(patient_id_mskcc == patient_uuid[i]) %>%
			 .[["sample_id_mskcc"]]
		invitae_id = manifest %>%
			     dplyr::filter(patient_id_mskcc == patient_uuid[i]) %>%
			     .[["sample_id_invitae"]]
		for (j in 1:length(msk_id)) {
			cat(paste0(msk_id[j], "-", invitae_id[j]), file = "../res/HPV.yaml", append = TRUE)
			if (j == length(msk_id)) {
				cat("]\n", file = "../res/HPV.yaml", append = TRUE)
			} else {
				cat(",", file = "../res/HPV.yaml", append = TRUE)
			}
		}
	}
}

#==================================================
# MRD VCF files
#==================================================
patient_uuid = manifest %>%
	       dplyr::filter(!is.na(vcf_file_name)) %>%
	       dplyr::arrange(patient_id_mskcc) %>%
	       .[["patient_id_mskcc"]] %>%
	       unique()

is_file = file.create("../res/MRD.yaml", showWarnings = TRUE)
if (is_file) {
	for (i in 1:length(patient_uuid)) {
		cat(paste0("- name: ", patient_uuid[i], "\n"), file = "../res/MRD.yaml", append = TRUE)
		cat("  tumor: [", file = "../res/MRD.yaml", append = TRUE)
		msk_id = manifest %>%
			 dplyr::filter(patient_id_mskcc == patient_uuid[i]) %>%
			 .[["sample_id_mskcc"]]
		invitae_id = manifest %>%
			     dplyr::filter(patient_id_mskcc == patient_uuid[i]) %>%
			     .[["sample_id_invitae"]]
		for (j in 1:length(msk_id)) {
			cat(paste0(msk_id[j], "-", invitae_id[j]), file = "../res/MRD.yaml", append = TRUE)
			if (j == length(msk_id)) {
				cat("]\n", file = "../res/MRD.yaml", append = TRUE)
			} else {
				cat(",", file = "../res/MRD.yaml", append = TRUE)
			}
		}
	}
}
