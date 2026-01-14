#!/bin/bash
source /disk5/luosg/DIPseq20251215/workflow/scripts/download/ascp.sh
source /disk5/luosg/DIPseq20251215/workflow/scripts/download/metadata.sh

outdir=/disk5/luosg/DIPseq20251215/data/fq
log=/disk5/luosg/DIPseq20251215/log/PRJ2SAMPLESingle_PRJNA292544.log
# PRJ2SAMPLESingle PRJNA292544 ${outdir} ${log}
log=/disk5/luosg/DIPseq20251215/log/ENAdownload_6mA.log
infile=/disk5/luosg/DIPseq20251215/data/download.tsv
outdir=/disk5/luosg/DIPseq20251215/data/fq
# key=/home/luosg/miniconda3/envs/RNA_SNP/etc/asperaweb_id_dsa.openssh
# awk -F"\t" '{print $4}' ${infile} | while read -r SRA;do
#     ENAdownload ${SRA} SINGLE ${log} ${outdir} ${key}
# done
# ENAdownload SRR5815057 SINGLE ${log} ${outdir} ${key}
source /disk5/luosg/DIPseq20251215/workflow/scripts/download/checkFq.sh
fqDir=/disk5/luosg/DIPseq20251215/data/fq
log=/disk5/luosg/DIPseq20251215/log/gzipTestForSra.log
awk -F"\t" '{print $4}' ${infile} | while read -r SRA;do
    gzipTestForSra ${fqDir} ${SRA} ${log}
done
