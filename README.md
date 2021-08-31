# covid-pipeline
---
## *How to use*:

1- Use conda to create an environment using Covid19.yaml
ex: `conda env create --name NAME --file Covid19.yaml`

`conda activate NAME`

2- Go to the directory you want to run the pipeline in
ex: cd /user/home/myfastq

3- Copy the files that were supported by Vipin to the previous directory from step 2 (Snakefile and config.yaml)  

4- Run the following:
`snakemake  -s  ./Snakefile ~/E81581/result/1.txt   --latency-wait 120 -j all -p`
E81581:- Directory contains 2 fastq files; like H3GLWDRXY-1-IDUSU178_S63_L001_R1_001.fastq.gz H3GLWDRXY-1-IDUSU178_S63_L001_R2_001.fastq.gz
Naming is important, it should be anything_R1_001.fastq.gz and anything_R1_001.fastq.gz
result:- is the name of the output directory (feel free to use other names)
1.txt:- anchor file (feel free to use other names)

5- To change any parameter refer to config.yaml

