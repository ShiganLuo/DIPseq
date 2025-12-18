from pathlib import Path
SNAKEFILE_FULL_PATH_Align = workflow.snakefile
SNAKEFILE_DIR_Align = os.path.dirname(SNAKEFILE_FULL_PATH_Align)
def get_yaml_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the workflow/RNA-SNP/snakemake/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR_Align ,f"{module_name}.yaml")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module configfile {module_name}.yaml not found at {module_path}")
    return module_path
alignYaml = get_yaml_path("Align")
configfile: alignYaml
logging.info(f"Include Align config: {alignYaml}")
rule trimming_Paired:
    input:
        fastq1 = indir + "/{sample_id}_1.fastq.gz",
        fastq2 = indir + "/{sample_id}_2.fastq.gz"
    output:
        fastq1 = temp(outdir + "/cutadapt/{sample_id}_1.fq.gz"),
        fastq2 = temp(outdir + "/cutadapt/{sample_id}_2.fq.gz"),
        report1 = outdir + "/log/Align/trimming/{sample_id}/trimming_statistics_1.txt",
        report2 = outdir + "/log/Align/trimming/{sample_id}/trimming_statistics_2.txt"
    params:
        outdir = outdir + "/cutadapt",
        quality = 30,
        trim_galore = config['tools']['trim_galore']
    threads: 6
    log:
        log = outdir + "/log/Align/{sample_id}/trimming.txt"
    shell:
        """
        # trim_galore can automatically judge the fq quality scoring system,it's no need to add such as --phred33 --phred64
        {params.trim_galore} --paired  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq1} {input.fastq2} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_val_1.fq.gz {output.fastq1}
        mv {params.outdir}/{wildcards.sample_id}_val_2.fq.gz {output.fastq2}
        mv {params.outdir}/{wildcards.sample_id}_1.fastq.gz_trimming_report.txt {output.report1}
        mv {params.outdir}/{wildcards.sample_id}_2.fastq.gz_trimming_report.txt {output.report2}
        """
rule trimming_Single:
    input:
        fastq = indir + "/{sample_id}.fastq.gz"
    output:
        fastq = temp(outdir + "/cutadapt/{sample_id}Single.fq.gz"),
        report = outdir + "/log/Align/trimming/{sample_id}/trimming_statistics.txt"
    params:
        outdir = outdir + "/cutadapt",
        quality = 30,
        trim_galore = config['tools']['trim_galore']
    threads: 6
    log:
        log = outdir + "/log/Align/{sample_id}/trimming.txt"
    shell:
        """
        {params.trim_galore} --phred33  --cores {threads} --quality {params.quality} \
            -o {params.outdir} --basename {wildcards.sample_id} {input.fastq} > {log.log} 2>&1
        mv {params.outdir}/{wildcards.sample_id}_trimmed.fq.gz {output.fastq}
        mv {params.outdir}/{wildcards.sample_id}.fastq.gz_trimming_report.txt {output.report}
        """

def get_alignment_input(wildcards):
    """
    function: Dynamically determines the input file type: paired-end or single-end sequencing.
    Based on the paired_samples and single_samples lists.This function is called in the star_align rule.

    param: 
        wildcards: Snakemake wildcards object containing the sample_id.
        paired_samples = ['sample1', 'sample2', ...]
        single_samples = ['sample3', 'sample4', ...]
    These lists must be defined in the Snakefile or config file.

    return: A list of input file paths for the STAR alignment step. 
    """
    logging.info(f"[get_alignment_input] called with wildcards: {wildcards}")
    # 构造可能的输入路径
    paired_r1 = f"{outdir}/cutadapt/{wildcards.sample_id}_1.fq.gz"
    paired_r2 = f"{outdir}/cutadapt/{wildcards.sample_id}_2.fq.gz"
    single = f"{outdir}/cutadapt/{wildcards.sample_id}Single.fq.gz"
    
    # 检查文件实际存在情况
    if wildcards.sample_id in paired_samples:
        logging.info(f"双端测序：{[paired_r1, paired_r2]}")
        return [paired_r1, paired_r2]
    elif wildcards.sample_id in single_samples:
        logging.info(f"单端测序：{[single]}")
        return [single]
    else:
        raise FileNotFoundError(
            f"Missing input files for sample {wildcards.sample_id}\n"
            f"Checked paths:\n- {paired_r1}\n- {paired_r2}\n- {single}"
        )


rule bowtie2_index_small:
    """
    fasta < 4GB
    """
    input:
        fasta = lambda wildcards: config['genome'][wildcards.genome]['fasta']
    output:
        bw2_1   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.1.bt2",
        bw2_2   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.2.bt2",
        bw2_3   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.3.bt2",
        bw2_4   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.4.bt2",
        bw2_rv1 = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.rev.1.bt2",
        bw2_rv2 = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.rev.2.bt2",
    log:
        outdir + "/log/Align/bowtie2/{genome}/bowtie2_index.log"
    params:
        bowtie2_build_cmd = config.get('tools', {}).get('bowtie2-build', 'bowtie2-build'),
        index_prefix = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome"
    conda:
        config['conda']['run']
    threads:
        8
    shell:
        """
        {params.bowtie2_build_cmd} --threads {threads} {input.fasta}  {params.index_prefix} > {log} 2>&1
        """

rule bowtie2_align:
    input:
        fastqs = get_alignment_input,
        bw2_1   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.1.bt2",
        bw2_2   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.2.bt2",
        bw2_3   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.3.bt2",
        bw2_4   = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.4.bt2",
        bw2_rv1 = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.rev.1.bt2",
        bw2_rv2 = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome.rev.2.bt2"
    output:
        bam = outdir + "/Align/bam/{genome}/{sample_id}.bam"
    log:
        outdir + "/log/Align/bowtie2/{genome}/{sample_id}/bowtie2_align.log"
    params:
        bowtie2_cmd = config.get('tools',{}).get('bowtie2','bowtie2'),
        samtools_cmd = config.get('tools',{}).get('samtools','samtools'),
        index_prefix = outdir + "/genome/{genome}/index/GRCm39.primary_assembly.genome"
    conda:
        config['conda']['run']
    threads:
        10
    shell:
        """
        # 将 input.fastqs 转换为 bash 数组
        fq_array=({input.fastqs})
        fq_count=${{#fq_array[@]}}

        if [ "$fq_count" -eq 2 ]; then
            # 双端逻辑
            {params.bowtie2_cmd} -x {params.index_prefix} \
                -1 "${{fq_array[0]}}" -2 "${{fq_array[1]}}" \
                -N 1 -L 30 \
                --threads {threads} 2> {log} \
                | {params.samtools_cmd} sort -@ {threads} -o {output.bam} - >> {log} 2>&1

        elif [ "$fq_count" -eq 1 ]; then
            # 单端逻辑
            {params.bowtie2_cmd} -x {params.index_prefix} \
                -U "${{fq_array[0]}}" \
                -N 1 -L 30 \
                --threads {threads} 2> {log} \
                | {params.samtools_cmd} sort -@ {threads} -o {output.bam} - >> {log} 2>&1
        else
            echo "Error: Expected 1 or 2 fastq files, got $fq_count" > {log}
            exit 1
        fi
        """
    
rule macs2:
    input:
        