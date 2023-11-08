#!/home/bioinformatics/software/R/R-4.2.2-gcc-8.5.0/bin/Rscript

# Demultiplex Metacell matrix to human and mouse submatrices (or just add the SoC output to the mat) by souporcell.

require(data.table, quietly = T)
require(metacell, quietly = T)
require(dplyr, quietly = T)
require(tibble, quietly = T)
require(MASS, quietly = T)
require(ggplot2, quietly = T)
require(class, quietly = T)

args = commandArgs(trailingOnly = T)

base_dir = args[1]
sample_name = args[2]
hg_pref = args[3]
mm_pref = args[4]

hg_mm_mode = hg_pref != "NA" && mm_pref != "NA"

samp_dir = sprintf("%s/%s/outs", base_dir, sample_name)

scdb_dir = paste0(base_dir, "/scdb") 
figs_dir = paste0(base_dir, "/figs")

scdb_init(scdb_dir)
scfigs_init(figs_dir)

mat = scdb_mat(sample_name)
stopifnot(!is.null(mat))

md = mat@cell_metadata

qc_ifn = sprintf("%s/report/tmp/%s_qcs.csv", base_dir, sample_name)
stopifnot(file.exists(qc_ifn))
qc_csv = read.csv(qc_ifn, header = T) %>%
  mutate(soc_tot_cells = NA,
         soc_singlet = NA,
         soc_singlet_human = NA,
         soc_singlet_mouse = NA,
         soc_unassigned = NA,
         soc_doublet = NA,
         soc_f_doublet = NA)

if (mm_pref == "NA") {
  qc_csv = qc_csv %>% 
    mutate(n_valid_mm = NA, f_valid_mm = NA, valid_mm_med_tot_umis = NA, valid_mm_med_nonMT_umis = NA, valid_mm_med_genes = NA, mm_min_nonMT_umi_cutoff = NA, mm_max_f_MT_cutoff = NA, mm_f_high_MT = NA, mm_high_MT_med_fMT = NA)
}
if (hg_pref == "NA") {
  qc_csv = qc_csv %>% 
    mutate(n_valid_hg = NA, f_valid_hg = NA, valid_hg_med_tot_umis = NA, valid_hg_med_nonMT_umis = NA, valid_hg_med_genes = NA, hg_min_nonMT_umi_cutoff = NA, hg_max_f_MT_cutoff = NA, hg_f_high_MT = NA, hg_high_MT_med_fMT = NA)
}

soc_clust_fn = sprintf("%s/souporcell/clusters.tsv", samp_dir)
if (file.exists(soc_clust_fn) && ! "soc_fine_assign" %in% colnames(md)) {
  message("adding souporcell clustering info...")
  
  soc = read.table(soc_clust_fn, header=T, stringsAsFactors = F)
  mds = md[soc$barcode, ]
  
  if (hg_mm_mode) {
    soc$is_hg = mds$hg_tot_umis > mds$mm_tot_umis & mds$valid_by_umis
    cl_t = table(soc[soc$status == "singlet", "assignment"], ifelse(soc[soc$status == "singlet", "is_hg"], 'Human', 'Mouse'))
    message("Human cells by souporcell clusters:")
    print(cl_t)
    
    cl_dict = apply(cl_t, 1, function(v) colnames(cl_t)[which.max(v)])
    soc$specie = cl_dict[as.character(soc$assignment)]
    bad_assign_ind = soc$is_hg & soc$specie == 'Mouse' | !soc$is_hg & soc$specie == 'Human'
    soc$specie[bad_assign_ind] = 'unassigned'
    soc$fine_assign = ifelse(soc$status == "singlet", paste(soc$status, soc$specie), soc$status)
    
    # deal with missed doublets (by soc hard cutoff of Prob[doub] > Prob[singlet]) - kmeans hg and mm umis, mark kmeans doublet cluster
    if (sum(soc$status == 'doublet') > 0) {
      # refine doublet assignment by clustering cell assignments on hg/mm umis (kmeans and knn)
      x = log2(1+mat@cell_metadata[soc$barcode, 'hg_tot_umis'])
      y = log2(1+mat@cell_metadata[soc$barcode, 'mm_tot_umis'])
      xy = cbind(x,y)
        
      # Kmeans
      # Select initial km centroids in the heart of the human, mouse and doublet clouds.
      centers = matrix(c(median(x[soc$is_hg]),  median(y[soc$is_hg]), 
                         median(x[!soc$is_hg]), median(y[!soc$is_hg]), 
                         median(x[soc$is_hg]),  median(y[!soc$is_hg])), nrow=3, byrow=T)
      km = try(expr=kmeans(xy, centers=centers))
      if (class(km) == 'try-error') {
        soc$doublet_by_kmeans = NA
      } else {
        d_cl = as.numeric(names(which.max(table(km$cluster[soc$status == 'doublet']))))
        soc$doublet_by_kmeans = km$cluster == d_cl
        if (d_cl != 3) {
          message(sprintf("Warning: Kmeans cluster of doublets is expected to be #3 and not #%d", d_cl))
        }
      }
      
      # Knn
      cl = soc$status 
      max_knn_iters = 40
      k = 30
      set.seed(42) # class::knn function is not deterministic, set seed for consistency
      for (i in 1:max_knn_iters) { 
        kcl = knn(train=xy, test=xy, cl=cl, k=k)
        cl = ifelse(cl == 'doublet' | kcl == 'doublet', 'doublet', 'singlet')
        if (sum(kcl == 'doublet') == sum(cl == 'doublet')) {
          break
        }
      }
      message(sprintf("After %d Knn iters: %d by SoC, %d kcl, %d both", i, sum(soc$status == 'doublet'), sum(kcl == 'doublet'), sum(soc$status == 'doublet' | kcl == 'doublet')))
      soc$doublet_by_knn = kcl == 'doublet'
    } else {
      soc$doublet_by_kmeans = F
      soc$doublet_by_knn = F
    }    
    
    qc_csv = qc_csv %>% 
      mutate(soc_tot_cells = nrow(soc),
             soc_singlet = sum(soc$fine_assign %in% c('singlet Human', 'singlet Mouse') & !soc$doublet_by_knn),
             soc_singlet_human = sum(soc$fine_assign == 'singlet Human' & !soc$doublet_by_knn),
             soc_singlet_mouse = sum(soc$fine_assign == 'singlet Mouse' & !soc$doublet_by_knn),
             soc_unassigned = sum(soc$fine_assign == "unassigned" & !soc$doublet_by_knn),
             soc_doublet = sum(soc$fine_assign == 'doublet' | soc$doublet_by_knn),
             soc_f_doublet = mean(soc$fine_assign == 'doublet' | soc$doublet_by_knn))
    
    bc_c = grepl('barcode', colnames(soc))
    colnames(soc)[!bc_c] = paste("soc", colnames(soc)[!bc_c], sep="_")
    
    md = mat@cell_metadata %>% 
      tibble::rownames_to_column(var='barcode') %>% 
      left_join(soc, by=c('barcode')) %>% 
      tibble::column_to_rownames(var="barcode")
    
    mat@cell_metadata = md
    scdb_add_mat(sample_name, mat)
    
    hg_cells = rownames(md)[!is.na(md$soc_fine_assign) & md$soc_fine_assign == "singlet Human"]
    mm_cells = rownames(md)[!is.na(md$soc_fine_assign) & md$soc_fine_assign == "singlet Mouse"]
    
    # selecting genes this way so eGFP will be in both if exists
    hg_genes = !grepl(sprintf("^%s_", mm_pref), mat@genes)
    mm_genes = !grepl(sprintf("^%s_", hg_pref), mat@genes)
    
    hg_m = mat@mat[hg_genes, hg_cells, drop=F]
    rownames(hg_m) = gsub(sprintf("%s_", hg_pref), "", rownames(hg_m))
    
    mm_m = mat@mat[mm_genes, mm_cells, drop=F]
    rownames(mm_m) = gsub(sprintf("%s_", mm_pref), "", rownames(mm_m))
    
    if (length(hg_cells) > 1) {
      scdb_add_mat(paste0(sample_name, "_hg"), tgScMat(hg_m, 'umi', md[hg_cells, , drop=F]))
      fwrite(as.data.frame(hg_cells), file=sprintf("%s/qc_hg_singlet_filt_barcodes.tsv.gz", samp_dir), quote=F, sep="\t", col.names = F, compress="gzip")
    }
    if (length(mm_cells) > 1) { 
      scdb_add_mat(paste0(sample_name, "_mm"), tgScMat(mm_m, 'umi', md[mm_cells, , drop=F]))
      fwrite(as.data.frame(mm_cells), file=sprintf("%s/qc_mm_singlet_filt_barcodes.tsv.gz", samp_dir), quote=F, sep="\t", col.names = F, compress="gzip")
    }
  } else {
    
    qc_csv = qc_csv %>% 
      mutate(soc_tot_cells = nrow(soc),
             soc_singlet = sum(soc$status == 'singlet'),
             soc_singlet_human = "NA",
             soc_singlet_mouse = "NA", 
             soc_unassigned = sum(soc$status == "unassigned"),
             soc_doublet = sum(soc$status == 'doublet'),
             soc_f_doublet = mean(soc$status == 'doublet'))
    
    bc_c = grepl('barcode', colnames(soc))
    colnames(soc)[!bc_c] = paste("soc", colnames(soc)[!bc_c], sep="_")
    
    mat@cell_metadata = mat@cell_metadata %>% 
      tibble::rownames_to_column(var='barcode') %>% 
      left_join(soc, by=c('barcode')) %>% 
      tibble::column_to_rownames(var="barcode")
    
    scdb_add_mat(sample_name, mat)
    
  }
}

write.csv(qc_csv, qc_ifn, quote=F, row.names = F)

if (hg_mm_mode) {
  mat = scdb_mat(sample_name)
  md = mat@cell_metadata
  
  
  # some plots
  
  status_nm = "soc_status"
  
  if (status_nm %in% colnames(md)) {
    #doub_by_km_nm = "soc_doublet_by_kmeans"
    doub_by_knn_nm = "soc_doublet_by_knn"
    
    mdf = md[!is.na(md[, status_nm]), ]
    status_cols = c(singlet="darkgrey", unassigned="blue", doublet="red")
  
    hg_x = log2(1+mdf$hg_tot_umis)
    mm_y = log2(1+mdf$mm_tot_umis)
    d = MASS::kde2d(hg_x, mm_y, n=50)
    
    png(scfigs_fn(sample_name, "soc_hg_vs_mm_umis"), 800, 400)
    par(mfrow=c(1,2))
    plot(hg_x, mm_y, pch=19, col=status_cols[mdf[, status_nm]], cex=0.7, xlab="Human UMIs (log2)", ylab="Mouse UMIs (log2)", main=paste(sample_name, "by SoC"))
    contour(d$x, d$y, d$z, drawlabels = F, col='black', add=T)
    legend("bottomleft", legend=names(status_cols), col=status_cols, pch=19, bty='n', cex=0.8)
    plot(hg_x, mm_y, pch=19, col=ifelse(mdf[, doub_by_knn_nm], 'purple', 'darkgrey'), cex=0.7, xlab="Human UMIs (log2)", ylab="Mouse UMIs (log2)", main=paste(sample_name, "by SoC"))
    legend("bottomleft", legend='Doublet by knn', col='purple', pch=19, bty='n', cex=0.8)
    dev.off()
    
    assign_nm = "soc_assignment"
    md_s = mdf[mdf[, status_nm] == 'singlet', ]
    cls = unique(md_s[, assign_nm])
    
    png(scfigs_fn(sample_name, "soc_hg_vs_mm_umis_by_cl"), length(cls) * 400, 400)
    par(mfrow=c(1,length(cls)))
    for (cl in cls) {
      ind = md_s[, assign_nm] == cl
      plot(log2(md_s$hg_tot_umis), log2(md_s$mm_tot_umis), pch=19, col="lightgray", cex=0.7, xlab="Human UMIs (log2)", ylab="Mouse UMIs (log2)", main=sprintf("%s by SoC: cl %s", sample_name, cl))
      points(log2(md_s$hg_tot_umis[ind]), log2(md_s$mm_tot_umis[ind]), pch=19, col="red", cex=0.7)
    }
    dev.off()
  } else {
    
  }
}

file.create(sprintf("%s/scrnaseq.done", samp_dir))
