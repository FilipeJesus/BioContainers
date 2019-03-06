#!/usr/bin/env sh

# sh complete_run.sh cores lib_type index bam read_len out_dir annotation_file exon_size

sed -i -e 's/4/'$1'/g' /usr/local/lib/python2.7/site-packages/misopy/settings/miso_settings.txt
samtools index -b $4 $4'.bai'

if [ $2 = "Paired_end" ]
then
  exon_utils --get-const-exons $7 --min-exon-size $8 --output-dir "exon"

  pe_utils --compute-insert-len $4 ./exon/*const_exons.gff --output-dir insert_dist/ | tee log.txt

  grep -A 1 'mean' log.txt > log

  mean=$(awk 'NR == 4 {print $1}' log)
  sd=$(awk 'NR == 4 {print $2}' log)
  echo $mean
  echo $sd

  miso --run $3 $4 --output-dir miso_out/ --read-len $5 --paired-end $mean $sd
fi
if [ $2 = "Single" ]
then
  miso --run $3 $4 --output-dir $6 --read-len $5
fi

summarize_miso --summarize-samples $6 $6'_summary'
