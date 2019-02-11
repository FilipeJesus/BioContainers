#!/bigseq/analysis/images/r_3.4.3_v2.simg

# Rscript ./GSEA.R --de_obj_file /bigseq/analysis/20180603_IL13_DbDbx4.rnaseq/Group_IL13_vs_NIP228.de-lst.rds --gene_set_list "/bigseq/datasets/20180603_Reactome_db/Reactome.Ensembl.Mus musculus.db-version-64.RDS" --gene_set_name reactome_ensembl --p.val 1

library(reshape2)

args <- commandArgs( trailingOnly=TRUE )
args

if (! "--p.val" %in% args) {
  cat("\n")
  stop("'Species Name' flag [--species] absent")
} else {
  pval  <- args[grep("--p.val", args)+1]
}

if ("--gene_set_list" %in% args) {
  gene_set.file.idx  <- grep("--gene_set_list", args)
  gene_set.file.path <- args[ gene_set.file.idx+1 ]
  if (file.exists(gene_set.file.path)) {
    cat("\nLoading _serialised_ 'Gene Set' object from \"",gene_set.file.path,"\" ... ", sep="")
    gs_Object <- readRDS(gene_set.file.path)
    cat("done\n")
  } else {
    cat("\nLocation of file containing 'Gene Set' object : \"",gene_set.file.path,"\"\n", sep="")
    stop("File **DOES NOT EXIST**")
  }
  if("--gene_set_name" %in% args){
    gene_set_name <- grep("--gene_set_name", args)
    gene_set_name <- args[ gene_set_name+1 ]
    gsc <- gs_Object[[gene_set_name]]
  } else {
    cat("\n")
    stop("'gene set name' flag [--gene_set_name] absent")
  }
}

## Read file containing the **serialised** 'SummarizedExperiment' object & check its class is correct ##
if (! "--de_obj_file" %in% args) {
  cat("\n")
  stop("'Differential Expression' flag [--de_obj_file] absent")
} else {
  de.file.path <- args[ grep("--de_obj_file", args)+1 ]
  
  if (file.exists(de.file.path)) {
    cat("\nLoading _serialised_ 'differential expression' object from \"", de.file.path ,"\" ... ", sep="")
    Deseq.results <- readRDS(de.file.path)
    cat("done\n")
  }
  else {
    cat("\nLocation of file containing 'differential expression' object : \"", de.file.path,"\"\n", sep="")
    stop("File **DOES NOT EXIST**")
  }
}

Deseq.results <- Deseq.results$Test.Res

gsc <- gsc[gsc[,1] %in% rownames(Deseq.results),]
Gene_Sets <- lapply(unique(gsc[,2]), function(x) gsc[gsc[,2] == x,1])
names(Gene_Sets) <- unique(gsc[,2])

stats <- Deseq.results[!is.na(Deseq.results$prob),]
stats <- stats[stats$fail.detect.cnt <= 1,]
tmp <- stats[,"log2FoldChange"]
names(tmp) <- rownames(stats)
gsea.Results <- fgsea::fgsea(Gene_Sets,tmp, 10000, minSize = 10)
Results <- data.frame(gsea.Results[,2:7], row.names = gsea.Results$pathway)
GSEA_Object <- list("type" = args[ grep("--gene_set_name", args)+1 ], "gsea_output" = Results)
saveRDS(GSEA_Object, file = paste0(gsub("de-lst.rds","",basename(de.file.path)), args[ grep("--gene_set_name", args)+1 ],"_gsea.RDS"))

