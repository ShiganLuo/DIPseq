import subprocess
import pandas as pd
import re
import os
from typing import List
import multiprocessing
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Literal

def computeMatrixMultiple(bws:List[str],regionBed:str,outdir:str,prefix:str, mode:str = "region"):
    """
    costruct computeMatrix command
    mode: region or point
    """
    suffix = os.path.basename(regionBed).replace(".bed", ".gz")
    outName = f"{prefix}_{suffix}"
    outfile = os.path.join(outdir, outName)
    if mode == "region":
        cmd = ["computeMatrix",
            "scale-regions",
            "-S"] + bws + [
                "-R",regionBed,
                "-a", "3000",
                "-b", "3000",
                "--skipZeros",
                "-p","6",
                "-o",outfile
            ]
    elif mode == "point":
            cmd = ["computeMatrix","reference-point",
                "--referencePoint","center",
                "-S"] + bws + [
                "-R",regionBed,
                "-a", "3000",
                "-b", "3000",
                "--skipZeros",
                "-p","6",
                "-o",outfile
            ]
    else:
        raise ValueError(f"Invalid mode: {mode}. Must be 'region' or 'point'.")
    return cmd

def computeMatrixSingle(bw,regionBed:str,outdir:str,mode:str,prefix:str):
    """
    costruct computeMatrix command
    """
    sample = os.path.basename(bw).replace(".bigWig","")
    suffix = os.path.basename(regionBed).replace(".bed", ".gz")
    outName = f"{sample}_{prefix}_{suffix}"
    outfile = os.path.join(outdir, outName)
    if mode == "region":
        cmd = ["computeMatrix",
            "scale-regions",
            "-S",bw,
                "-R",regionBed,
                "-a", "3000",
                "-b", "3000",
                "--skipZeros",
                "-p","2",
                "-o",outfile
            ]
    elif mode == "point":
            cmd = ["computeMatrix","reference-point",
            "--referencePoint","center",
            "-S",bw,
            "-R",regionBed,
            "-a", "3000",
            "-b", "3000",
            "--skipZeros",
            "-p","2",
            "-o",outfile
            ]

    return cmd

def plotHeatmap(matrix:str, outdir:str):
    """
    costruct plotHeatmap command
    """
    outName = os.path.basename(matrix).replace(".gz", ".png")
    outfile = os.path.join(outdir, outName)
    cmd = ["plotHeatmap",
           "-m", matrix,
           "-out", outfile,
           "--dpi", "300",
           "--heatmapWidth","12"
           ]
    return cmd

def run_command(cmd:list):
    print("Executing command:", " ".join(cmd))
    subprocess.run(cmd)
    
    
if __name__ == "__main__":
    bws = ["/disk5/luosg/DIPseq20251215/output/Align/bam/mouse/bigwig/GSM1847697_25_cpm.bigwig",
            "/disk5/luosg/DIPseq20251215/output/Align/bam/mouse/bigwig/GSM1847699_25_cpm.bigwig",
            "/disk5/luosg/DIPseq20251215/output/Align/bam/mouse/bigwig/GSM1847701_25_cpm.bigwig"]
    cmd1 = computeMatrixMultiple(bws,regionBed="/disk5/luosg/DIPseq20251215/output/Diff/gene/L3mbtl1.bed",outdir="/disk5/luosg/DIPseq20251215/output/deeptools/L3mbtl1",
                           prefix="L3mbtl1")
    # print(cmd) 
    cmd2 = plotHeatmap("/disk5/luosg/DIPseq20251215/output/deeptools/L3mbtl1/L3mbtl1_L3mbtl1.gz","/disk5/luosg/DIPseq20251215/output/deeptools/L3mbtl1")
    with multiprocessing.Pool(processes=1) as pool:
        pool.map(run_command, [cmd2])