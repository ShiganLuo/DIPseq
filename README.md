# DIPseq

DNA免疫沉淀测序

1. 创建envs中的conda环境，对于各个subworflow，配置好同名yaml内容

注意:
    - trim_galore只是打包程序，需要确保cutadapt存在
    - 物种名是动态变化的，在yaml中，相关文件的键需要根据metadata改变，如果物种名为有空格隔开，需要改空格为_，如Mus musculus -> Mus_musculus

2. 运行

```sh
snakemake -s workflow/RNA-SNP/run.smk --config indir=data/fq outdir=output metadata=data/target_fq.tsv --cores 45 --use-conda
```
