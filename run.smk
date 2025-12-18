shell.prefix("set -x; set -e;")
import logging
import os
from itertools import chain
import sys
from snakemake.io import glob_wildcards
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    stream=sys.stdout,
	datefmt='%Y-%m-%d %H:%M:%S'
)
# containerize: "quay.nju.edu.cn"
EXECUTION_DIR = os.getcwd()
SNAKEFILE_FULL_PATH = workflow.snakefile
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE_FULL_PATH)
indir = config.get('indir', '../data')
outdir = config.get('outdir', '../output')
metadata = config.get('metadata')
logging.info("Workflow started.")
logging.info(f"metadata file path: {metadata}")
logging.info(f"Input directory: {indir}")
logging.info(f"Output directory: {outdir}")
logging.info(f"Snakefile path: {SNAKEFILE_FULL_PATH}")
logging.info(f"Execution directory: {EXECUTION_DIR}")
sys.path.append(f"{SNAKEFILE_DIR}/utils")
from fastq_utils import SNPMetadata
snpMetadata = SNPMetadata(metadata,indir,f"{outdir}/log/utils/fastq_utils.log")
groups = snpMetadata.run()

def get_output_files(groups):
    outfiles = []
    paired_samples = []
    single_samples = []
    all_samples = []
    genomes = []
    single_sample_genome_pairs = []
    paired_sample_genome_pairs = []
    for organism, types in groups.items():
        genomes.append(organism)
        for TYPE, samples in types.items():
            if TYPE == "PAIRED":
                for sample in samples:
                    paired_samples.append(sample)
                    outfiles.append(f"{outdir}/Align/bam/{organism}/{sample}.bam")
            elif TYPE == "SINGLE":
                for sample in samples:
                    single_samples.append(sample)
                    outfiles.append(f"{outdir}/Align/bam/{organism}/{sample}.bam")
            else:
                continue
    all_samples = paired_samples + single_samples
    return outfiles, paired_samples, single_samples, all_samples, genomes, single_sample_genome_pairs, paired_sample_genome_pairs
outfiles, paired_samples, single_samples, all_samples, genomes, single_sample_genome_pairs, paired_sample_genome_pairs = get_output_files(groups)

logging.info(f"genomes: {genomes}\npaired_samples: {paired_samples}\nsingle_samples: {single_samples}\nall_samples: {all_samples}\n\
single_sample_genome_pairs: {single_sample_genome_pairs}\npaired_sample_genome_pairs: {paired_sample_genome_pairs}\n\
outfiles: {outfiles}")

configfilePath = os.path.join(SNAKEFILE_DIR,"config","run.yaml")
configfile: configfilePath
logging.info(f"add cofigfile {configfilePath}")

def get_snakefile_path(module_name:str)->str:
    """
    function: Get the absolute path of a module in the core_snakefile_path/subworkflow/ directory.

    param: 
        module_name: Name of the module (without .smk extension).

    return: Absolute path of the module file.
    """
    module_path = os.path.join(SNAKEFILE_DIR, "subworkflow",module_name, f"{module_name}.smk")
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"Module snakefile {module_name}.smk not found at {module_path}")
    return module_path

alignSmk = get_snakefile_path("Align")
include: alignSmk
logging.info(f"Include Align workflow: {alignSmk}")

rule all:
    input:
        outfiles
        # outdir + "/multiqc/multiqc_report.html",




