#!/home/bioinformatics/software/R/R-4.2.2-gcc-8.5.0/bin/Rscript

# Bundle up QCs figures to ppptx slides, one slide per sample. Currently supports only 'hg_mm' in hg_mm_mode 
require(metacell, quietly = T)
require(officer, quietly = T)

args = commandArgs(trailingOnly = T)
if (length(args) == 0) {
  message("Usage: scRNA_qcs_to_ppt.R <base dir> <hg_mm_mode>\n")
  stop()
}

base_dir = args[1]
hg_mm_mode = args[2]

message(sprintf("Creating pptx slides with QC figs for %s (%s mode)", base_dir, hg_mm_mode))

figs_dir = paste0(base_dir, "/figs")
scfigs_init(figs_dir)

meta = read.csv(list.files(sprintf("%s/metadata", base_dir), pattern="*.csv", full.names = T))

bw = 1.875
sp = 0.25
pres = read_pptx()

for (sample_name in meta$Sample.name) {  
  pres = add_slide(pres)
  pres = ph_with(pres, value=fpar(ftext(sample_name, fp_text(font.size = 20))), location=ph_location(2*bw + sp, 0, bw*3, 0.5))

  hg_exists = file.exists(scfigs_fn(sample_name, "hg_all_UMIs"))
  mm_exists = file.exists(scfigs_fn(sample_name, "mm_all_UMIs"))
  hgmm_exists = file.exists(scfigs_fn(sample_name, "all_UMIs_by_specie"))
  
  if (hg_mm_mode == "hg_mm" && hgmm_exists) {
    pres = ph_with(pres, external_img(scfigs_fn(sample_name, "all_UMIs_by_specie"), bw, bw), location = ph_location(0, 0, bw, bw))
    pres = ph_with(pres, external_img(scfigs_fn(sample_name, "nonMT_UMIs_by_specie"), bw, bw), location = ph_location(bw, 0, bw, bw))
    pres = ph_with(pres, external_img(scfigs_fn(sample_name, "fMito_vs_logUMIs"), 2.5*bw, bw), location = ph_location(bw*2 + sp, bw, 2.5 * bw, bw))
    pres = ph_with(pres, external_img(scfigs_fn(sample_name, "hg_vs_mm_UMIs"), 2*bw, bw), location = ph_location(0, bw, 2*bw, bw))
    
    if (file.exists(scfigs_fn(sample_name, "fAmb_by_specie"))) {
      pres = ph_with(pres, external_img(scfigs_fn(sample_name, "fAmb_by_specie"), bw, bw), location = ph_location(0, 2*bw, bw, bw))
      pres = ph_with(pres, external_img(scfigs_fn(sample_name, "fAmb_vs_umis_by_specie"), 2*bw, bw), location = ph_location(bw + sp, 2*bw, 2*bw, bw))
      pres = ph_with(pres, external_img(scfigs_fn(sample_name, "fMT_vs_Amb_by_specie"), 2*bw, bw), location = ph_location(3*bw + 2*sp, 2*bw, 2*bw, bw))
    }
    
    if (file.exists(scfigs_fn(sample_name, "soc_hg_vs_mm_umis"))) {
      pres = ph_with(pres, external_img(scfigs_fn(sample_name, "soc_hg_vs_mm_umis"), 2*bw, bw), location = ph_location(0, 3*bw, 2*bw, bw))
      pres = ph_with(pres, external_img(scfigs_fn(sample_name, "soc_hg_vs_mm_umis_by_cl"), 2*bw, bw), location = ph_location(bw*2 + sp, 3*bw, 2*bw, bw))
      
    }
  } 
  else {
    stopifnot(hg_mm_mode != "hg_mm" || !all(c(hg_exists, mm_exists)))

    selected_specie = hg_mm_mode
    
    # When hg_mm mode but only hg/mm cells to figs are hg/mm only
    if (hg_mm_mode == "hg_mm") { 
      selected_specie = ifelse(hg_exists, "hg", "mm")
    }
    
    pres = ph_with(pres, external_img(scfigs_fn(sample_name, paste0(selected_specie, "_all_UMIs")), bw, bw), location = ph_location(0, 0, bw, bw))
    pres = ph_with(pres, external_img(scfigs_fn(sample_name, paste0(selected_specie, "_nonMT_UMIs")), bw, bw), location = ph_location(bw, 0, bw, bw))
    pres = ph_with(pres, external_img(scfigs_fn(sample_name, paste0(selected_specie, "_fMito_vs_logUMIs")), bw, bw), location = ph_location(bw*2 + sp, bw, bw, bw))
    
    if (file.exists(scfigs_fn(sample_name, paste0(selected_specie, "_fAmb")))) {
      pres = ph_with(pres, external_img(scfigs_fn(sample_name, paste0(selected_specie, "_fAmb")), bw, bw), location = ph_location(0, 2*bw, bw, bw))
      pres = ph_with(pres, external_img(scfigs_fn(sample_name, paste0(selected_specie, "_fAmb_vs_umis")), bw, bw), location = ph_location(bw + sp, 2*bw, bw, bw))
      pres = ph_with(pres, external_img(scfigs_fn(sample_name, paste0(selected_specie, "_fMT_vs_Amb")), bw, bw), location = ph_location(2*bw + 2*sp, 2*bw, bw, bw))
    }
    
  }
}
print(pres, target=sprintf("%s/report/%s_QCs.pptx", base_dir, basename(base_dir)))
 