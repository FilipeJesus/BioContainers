#!/bin/bash

# ./GSEA-init.sh /bigseq/analysis/20180613_CVRM_DbDbx4_timecourse.New.Group_week_20_vs_week_13.GSEA/*.de-lst.rds "/bigseq/datasets/20180603_Reactome_db/mmu_reactome_annotation.v92.rds" reactome_ensembl 1

echo -ne "\nBaseSpace directory='$1'\n\n"

for dir in $1; do
  echo -ne "\nDir='$dir'\n\n"
  Rscript ./pathway_analysis/GSEA.R --de_obj_file $dir --gene_set_list $2 --gene_set_name $3 --p.val $4

done