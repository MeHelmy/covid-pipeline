# import Lib
############
from snakemake.utils import min_version
import os
import glob
############################

# Snake Version
###############
min_version("5.7.1")


# Config File
#############
if os.path.isfile("config.yaml"):

    configfile: "config.yaml"


else:
    sys.exit(
        "Looks like there is no config.yaml file in "
        + os.getcwd()
        + " make sure there is one or at least specify one with the --configfile commandline parameter."
    )




# Functions
############
def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    return sorted(
        glob.glob(
            os.path.join(wildcards.mydir, wildcards.sampledir) + "/*_R1_001.fastq.gz"
        )
    )


def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    return sorted(
        glob.glob(
            os.path.join(wildcards.mydir, wildcards.sampledir) + "/*_R2_001.fastq.gz"
        )
    )



def get_reads_for_bam(wildcards):
    r1 = os.path.join(
        wildcards.mydir,
        wildcards.sampledir,
        wildcards.output,
        config["trim_dir"] + "/paired_R1.fastq.gz",
    )
    r2 = os.path.join(
        wildcards.mydir,
        wildcards.sampledir,
        wildcards.output,
        config["trim_dir"] + "/paired_R2.fastq.gz",
    )
    return {"r1": r1, "r2": r2}



onstart:
    shell("cat etc/start.txt")


rule all:
    input:
        "{mydir}/{sampledir}/{output}/{sample}.txt",


rule collect_analysis:
    input:
        qc="{mydir}/{sampledir}/{output}/.{sample}_QC.done",
        # trim = "{mydir}/{sampledir}/{output}/.{sample}.Trim.done",
        variant="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}_{pre}.tsv".format(
            pre=config["var_prefix"]
        ),
        variant_vcf="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}_{pre}.tsv".format(
            pre=config["var_prefix"]
        ),
        # con = "{{mydir}}/{{sampledir}}/{{output}}/{pre}.fa".format(pre=config['con_prefix']),
        compare="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}_{pre}.json".format(
            pre=config["next_prefix"]
        ),
        pan="{{mydir}}/{{sampledir}}/{{output}}/{outdir}/{{sample}}_{pre}.csv".format(
            outdir=config["pan_dir_prefix"], pre=config["pan_prefix"]
        ),
        # align_stat="{mydir}/{sampledir}/{output}/{{sample}}_{pre}.csv",
        align_human_matrix="{mydir}/{sampledir}/{output}/{sample}.human_metrics.tsv",
    output:
        "{mydir}/{sampledir}/{output}/{sample}.txt",
    shell:
        """
        touch {output}
        """


rule QC:
    """
    Input: pair-end Reads
    output: QC for reads in directory [QC] by default, can be changed from config.yaml file
    Method: uses fastqc
    """
    input:
        read1=get_fq1,
        read2=get_fq2,
    output:
        "{mydir}/{sampledir}/{output}/.{sample}_QC.done",
    message:
        "Running QC for input data {input.read1} and {input.read2}"
    threads: config["fastqc_threads"]
    params:
        dir_name=config["qc_dir"],
    shell:
        """
        mkdir -p {wildcards.mydir}/{wildcards.sampledir}/{wildcards.output}/{params.dir_name} && fastqc -t {threads} {input.read1} {input.read2} -o {wildcards.mydir}/{wildcards.sampledir}/{wildcards.output}/{params.dir_name} &&\
        touch {output}
        """


rule TrimData:
    """
    Input: reads infered from the directory.
    output: trimmed data in directory specified by trim_dir in config file
    """
    input:
        read1=get_fq1,
        read2=get_fq2,
    # output:"{mydir}/{sampledir}/{output}/.{sample}.Trim.done"
    output:
        trim1=(
            "{{mydir}}/{{sampledir}}/{{output}}/{trim_dir}/paired_R1.fastq.gz".format(
                trim_dir=config["trim_dir"]
            )
        ),
        trim2=(
            "{{mydir}}/{{sampledir}}/{{output}}/{trim_dir}/paired_R2.fastq.gz".format(
                trim_dir=config["trim_dir"]
            )
        ),
        trim1_unpaired=(
            "{{mydir}}/{{sampledir}}/{{output}}/{trim_dir}/unpaired_R1.fastq.gz".format(
                trim_dir=config["trim_dir"]
            )
        ),
        trim2_unpaired=(
            "{{mydir}}/{{sampledir}}/{{output}}/{trim_dir}/unpaired_R2.fastq.gz".format(
                trim_dir=config["trim_dir"]
            )
        ),
    # log: "{mydir}/{sampledir}/{output}/{sample}.Trim.log"
    message:
        "Triming {input.read1} and {input.read2}"
    threads: config["trim_threads"]
    params:
        dir_name=config["trim_dir"],
        summary="{mydir}/{sampledir}/{output}/Trim.summary",
        trim_option=config["trim_option"],
    shell:
        """
        trimmomatic PE -threads {threads} -trimlog {wildcards.mydir}/{wildcards.sampledir}/{wildcards.output}/Trim.log -summary {params.summary} -quiet \
        {input.read1} {input.read2} \
        {output.trim1} {output.trim1_unpaired} \
        {output.trim2} {output.trim2_unpaired} \
        {params.trim_option}
        """


rule AlignReads:
    input:
        unpack(get_reads_for_bam),
    output:
        "{mydir}/{sampledir}/{output}/{sample}.align.bam",
    message:
        "Align {input} data"
    log:
        "{mydir}/{sampledir}/{output}/{sample}.align.log",
    threads: config["align_threads"]
    params:
        other_options=config["bam_other_options"],
        ref=config["virus_ref"],
        bam_other_options=(
            config["bam_other_options"] if config["bam_other_options"] else ""
        ),
    shell:
        """
        bwa mem -U 17 -M -t {threads} {params.ref} {input.r1} {input.r2} | samtools sort -@ {threads}  -  > {output} 2> {log}
        """


rule IndexBam:
    input:
        "{dir}/{sample}.bam",
    output:
        "{dir}/{sample}.bam.bai",
    message:
        "Index bam file"
    shell:
        """
        samtools index {input}
        """


# rule VirusAlign:
#     """
#     Get virus reads from aligned bam
#     """
#     input:
#         bam="{mydir}/{sampledir}/{output}/align.bam",
#         bam_index="{mydir}/{sampledir}/{output}/align.bam.bai",
#     output:
#         bam="{mydir}/{sampledir}/{output}/align.virus.bam",
#     message:
#         "Select Virus reads from {input}"
#     threads: config["align_threads"]
#     params:
#         virus_contig=config["covid_contig"],
#     shell:
#         """
#         samtools view -h {input.bam} '{params.virus_contig}' | awk '!/^@/ {{print}} $0 ~ /SN:{params.virus_contig}/ {{print}} ($0 ~/^@/ && !/^@SQ/) {{print}} /$7="="/ {{print}}' | samtools view -F 4 -bS - > {output.bam} 2> /dev/null
#         """
# samtools view -h {input.bam} '{params.virus_contig}' | awk '!/^@/ {{print}} !/^@SQ/ {{print}} $0 ~ /"{params.virus_contig}"/ {{print}} /$7="="/ {{print}}' | samtools view -bS - > {output.bam} 2> /dev/null


rule ExtractNonVirusReads:
    """
    Aligning non virus reads to human genome
    """
    input:
        # bam="{{mydir}}/{{sampledir}}/{{output}}/{pre}.bam".format(pre=config["prefix"]),
        bam="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.sorted.bam".format(pre=config["prefix"]),
        bam_index="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.sorted.bam.bai".format(pre=config["prefix"]),
    output:
        bam=temp("{mydir}/{sampledir}/{output}/{sample}.non.virus.bam"),
    message:
        "Extract orphen reads and unaligned virus reads"
    log:
        "{mydir}/{sampledir}/{output}/{sample}.align.human.txt",
    threads: config["align_threads"]
    shell:
        """
        samtools  view -@ {threads} -h -f 4 -O BAM -o {output.bam}  {input.bam} > {log}
        """


rule ExtractHumanReads:
    """
    Extract non virus reads
    """
    input:
        bam="{mydir}/{sampledir}/{output}/{sample}.non.virus.bam",
        bam_index="{mydir}/{sampledir}/{output}/{sample}.non.virus.bam.bai",
    output:
        read1="{mydir}/{sampledir}/{output}/{sample}.non.virus.R1.001.fastq",
        read2="{mydir}/{sampledir}/{output}/{sample}.non.virus.R2.001.fastq",
        unpaired="{mydir}/{sampledir}/{output}/{sample}.non.virus.unpaired.fastq",
    message:
        "Convert reads from {input.bam} to fastq"
    log:
        "{mydir}/{sampledir}/{output}/{sample}.non.virus.reads.txt",
    threads: config["align_threads"]
    shell:
        """
        picard SamToFastq\
        I={input.bam}\
        F={output.read1}\
        F2={output.read2}\
        UNPAIRED_FASTQ={output.unpaired}\
        VALIDATION_STRINGENCY=LENIENT  > {log}
        """


rule Align2Human:
    input:
        R1="{mydir}/{sampledir}/{output}/{sample}.non.virus.R1.001.fastq",
        R2="{mydir}/{sampledir}/{output}/{sample}.non.virus.R2.001.fastq",
    output:
        "{mydir}/{sampledir}/{output}/{sample}.human.bam",
    message:
        "Aligning non virus reads to human"
    log:
        "{mydir}/{sampledir}/{output}/{sample}.human.align.txt",
    threads: config["align_threads"]
    params:
        ref=config["human_ref"],
    shell:
        """
        bwa mem -U 17 -M -t {threads} {params.ref} {input.R1} {input.R2} | samtools sort -@ {threads}  -  > {output} 2> {log}
        """


rule CollectHumanStatMetrics:
    """
    Collect Human aligned reads metrics.
    """
    input:
        bam="{mydir}/{sampledir}/{output}/{sample}.human.bam",
        bam_index="{mydir}/{sampledir}/{output}/{sample}.human.bam.bai",
    output:
        "{mydir}/{sampledir}/{output}/{sample}.human_metrics.tsv",
    message:
        "Getting aligning metrics for {input}"
    log:
        "{mydir}/{sampledir}/{output}/{sample}.human_metrics.log",
    params:
        ref=config["human_ref"],
    shell:
        """
        picard CollectAlignmentSummaryMetrics \
        R={params.ref} \
        I={input.bam}\
        O={output} > {log}
        """


rule TrimBam:
    input:
        bam="{mydir}/{sampledir}/{output}/{sample}.align.bam",
    output:
        bam="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.sorted.bam".format(pre=config["prefix"]),
    message:
        "Trimming Bam file {input}"
    log:
        "{mydir}/{sampledir}/{output}/{sample}.align.virus.trim.log",
    params:
        bed=config["trim_bed"],
        pre=config["prefix"],
        m=config["min_length"],
        q=config["min_quality"],
        s=config["slid"],
        x=config["offset_trimming"],
    shell:
        """
        ivar trim -x {params.x} -e -i {input.bam} -b {params.bed} -p {params.pre} -m {params.m} -q {params.q} -s {params.s} 2> {log} && samtools sort  {params.pre}.bam  -o {output.bam} && rm {params.pre}.bam
        """


# rule sortBam:
#     input:
#         "{anything}.bam",
#     output:
#         "{anything}.sorted.bam",
#     message:
#         "Sorting {input}"
#     shell:
#         """
#         samtools sort {input} -o {output}
#         """

rule CallVariant:
    input:
        bam="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.sorted.bam".format(
            pre=config["prefix"]
        ),
        bam_index="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.sorted.bam.bai".format(
            pre=config["prefix"]
        ),
    output:
        variant="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}_{pre}.tsv".format(
            pre=config["var_prefix"]
        ),
    message:
        "Calling varaint from {input}"
    params:
        ref=config["virus_ref"],
        pre=config["var_prefix"],
        q=config["var_min_quality"],
        t=config["var_min_freq"],
    shell:
        """
        samtools mpileup --reference {params.ref} -A -d 600000 -F 0 -B -Q 0 {input.bam} |  ivar variants -p {params.pre} -q {params.q} -t {params.t} && mv {params.pre}.tsv {output.variant}
        """


rule Consensus:
    input:
        "{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.sorted.bam".format(
            pre=config["prefix"]
        ),
    output:
        fa="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.fa".format(
            pre=config["con_prefix"]
        ),
        quality="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.qual.txt".format(
            pre=config["con_prefix"]
        ),
    message:
        "calling consensus for {input}"
    params:
        pre=config["con_prefix"],
        q=config["con_min_quality"],
        t=config["con_min_freq"],
    shell:
        """
        samtools mpileup -A -d 300000 -Q 0 -F 0 {input} | ivar consensus -p {params.pre} -q 20 -t 0 && mv {params.pre}.fa {output.fa} &&  mv  {params.pre}.qual.txt {output.quality}
        """


rule NextClade:
    input:
        "{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.fa".format(pre=config["con_prefix"]),
    output:
        compare="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}_{pre}.json".format(
            pre=config["next_prefix"]
        ),
        tree="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.tree".format(
            pre=config["next_prefix"]
        ),
    message:
        "Running NextClade for {input}"
    shell:
        """
        nextclade --input-fasta {input} -o {output.compare}  -T  {output.tree}
        """


rule Pangolin:
    input:
        "{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.fa".format(pre=config["con_prefix"]),
    output:
        "{{mydir}}/{{sampledir}}/{{output}}/{outdir}/{{sample}}_{pre}.csv".format(
            outdir=config["pan_dir_prefix"], pre=config["pan_prefix"]
        ),
    message:
        "running Pangolin for {input}"
    params:
        mydir=config["pan_dir_prefix"],
        outfile=config["pan_prefix"],
    threads: config["pan_threads"]
    shell:
        """
        pangolin {input} -o  {wildcards.mydir}/{wildcards.sampledir}/{wildcards.output}/{params.mydir} --outfile {wildcards.sample}_{params.outfile}.csv -t {threads}
        """


rule TSV2VCF:
    input:
        variant="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}_{pre}.tsv".format(
            pre=config["var_prefix"]
        ),
    output:
        variant="{{mydir}}/{{sampledir}}/{{output}}/{{sample}}_{pre}.vcf".format(
            pre=config["var_prefix"]
        ),
    message:
        "Converting variant to VCF"
    log:
        "{{mydir}}/{{sampledir}}/{{output}}/{{sample}}.{pre}.log".format(pre=config["var_prefix"]),
    params:
        pass_only="--pass_only" if config["pass_only"] else "",
        af=config["af"],
    shell:
        """
        scripts/ivar_variants_to_vcf.py -af {params.af} {params.pass_only} {input.variant} {output.variant} > {log}
        """


rule AlignStat:
    """
    Get result statitics for aligned bam
    """
    input:
        "{mydir}/{sampledir}/{output}/{sample}.align.bam",
    output:
        "{mydir}/{sampledir}/{output}/.{sample}.align.txt",
    message:
        "Calculate statistics for {input}"
    threads: config["align_stat_threads"]
    params:
        file_name=config["stat_output"],
    shell:
        """
        samtools stats -@ {threads} {input} > {wildcards.mydir}/{wildcards.sampledir}/{wildcards.output}/{params.file_name}.txt && touch {output}
        """


onsuccess:
    shell("cat etc/end.txt")
