#!/bin/bash
function macs3_callpeak(){
    local bam_treatment=$1
    local bam_control=$2
    local outdir=$3
    local name=$4
    macs3 callpeak \
    -t ${bam_treatment} \
    -c ${bam_control} \
    --bw 200 \
    -p 1e-5 \
    -g mm \
    --outdir ${outdir} \
    --name ${name} \
    --seed 2346
}
export -f macs3_callpeak
bam_Alkbh1_KO=/disk5/luosg/DIPseq20251215/output/Align/bam/mouse/GSM1847699.bam
bam_Alkbh1_wild=/disk5/luosg/DIPseq20251215/output/Align/bam/mouse/GSM1847697.bam
bam_control=/disk5/luosg/DIPseq20251215/output/Align/bam/mouse/GSM1847701.bam
macs3_callpeak ${bam_Alkbh1_KO} ${bam_control} /disk5/luosg/DIPseq20251215/output/macs3 Alkbh1_KO
macs3_callpeak ${bam_Alkbh1_wild} ${bam_control} /disk5/luosg/DIPseq20251215/output/macs3 Alkbh1_wild

