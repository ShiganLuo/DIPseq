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

function bam2bw(){
    # 6mA峰一般不宽，binSize 25
    local bam=$1
    local bigwig=$2
    bamCoverage -b ${bam} \
        -o ${bigwig} \
        --binSize 25 \
        --normalizeUsing CPM \
        --extendReads 200 \
        --ignoreDuplicates 
}
# macs3_callpeak ${bam_Alkbh1_KO} ${bam_control} /disk5/luosg/DIPseq20251215/output/macs3 Alkbh1_KO
# macs3_callpeak ${bam_Alkbh1_wild} ${bam_control} /disk5/luosg/DIPseq20251215/output/macs3 Alkbh1_wild
# samtools index ${bam_Alkbh1_KO}
# samtools index ${bam_Alkbh1_wild}
# samtools index  ${bam_control}
bam2bw ${bam_Alkbh1_KO} /disk5/luosg/DIPseq20251215/output/Align/bam/mouse/bigwig/GSM1847699_25_cpm.bigwig
bam2bw ${bam_Alkbh1_wild} /disk5/luosg/DIPseq20251215/output/Align/bam/mouse/bigwig/GSM1847697_25_cpm.bigwig
bam2bw ${bam_control} /disk5/luosg/DIPseq20251215/output/Align/bam/mouse/bigwig/GSM1847701_25_cpm.bigwig

# python /disk5/luosg/DIPseq20251215/workflow/DIPseq/scripts/peak/diff_peak_gene.py \
#     --ko-peaks /disk5/luosg/DIPseq20251215/output/macs3/Alkbh1_KO_peaks.narrowPeak \
#     --wt-peaks /disk5/luosg/DIPseq20251215/output/macs3/Alkbh1_wild_peaks.narrowPeak \
#     --ko-ip /disk5/luosg/DIPseq20251215/output/Align/bam/mouse/GSM1847699.bam \
#     --wt-ip /disk5/luosg/DIPseq20251215/output/Align/bam/mouse/GSM1847697.bam \
#     --outdir /disk5/luosg/DIPseq20251215/output/Diff/gene \
#     -g /disk5/luosg/Reference/GENCODE/mouse/GRCm39/gencode.vM38.primary_assembly.basic.annotation.gtf