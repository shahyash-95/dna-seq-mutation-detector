configfile : "config.yaml"
from glob import glob
import pandas as pd
import os
import boto3
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()
FILE_PATH = config["SAMPLE_DETAILS"]
df = pd.read_csv(FILE_PATH, index_col=None, header=None)
SAMPLES = list(df[0])
FILE_PATH = FILE_PATH.strip('sample_details.csv')
os.chdir(FILE_PATH)
##### code to download sequences file from s3
S3PATH=config["S3PATH"]
print('Working paths:',FILE_PATH, S3PATH) 
S3PATH=S3PATH.replace('s3://','')
s3list = S3PATH.split('/')
BUCKET_NAME = s3list[0]
PREFIX_PATH = '/'.join(s3list[1:])+'sequences'
if not os.path.exists('sequences'):
    os.makedirs('sequences')
s3_resource = boto3.resource('s3')
my_bucket = s3_resource.Bucket(BUCKET_NAME)
objects = my_bucket.objects.filter(Prefix=PREFIX_PATH)

for obj in objects:
    path, filename = os.path.split(obj.key)
    if not (os.path.isfile('sequences/'+filename)):
        if(filename.split('_')[0] in SAMPLES):
            print('Downloading file:', filename)
            my_bucket.download_file(obj.key, 'sequences/'+filename)

    #print('Path:'+path, 'filename:'+filename)

#variants=['mutect2','freebayes','platypus','vardict','varscan', 'lofreq', 'strelka']
#variants=['somaticseq', 'combined', 'freebayes','platypus']

rule all:
    input:
        S3.remote(expand(S3PATH+'Final-Output/{sample}/{sample}.Low_Coverage.png', sample=SAMPLES)),
        S3.remote(expand(S3PATH+'Final-Output/{sample}/{sample}.qc_report.pdf',sample=SAMPLES)),
        S3.remote(expand(S3PATH+'Final-Output/{sample}/{sample}.xlsx',sample=SAMPLES)),
        S3.remote(expand(S3PATH+'Final-Output/{sample}/{sample}.final.bam',sample=SAMPLES)),
        S3.remote(expand(S3PATH+'Final-Output/{sample}/{sample}.final.bam.bai',sample=SAMPLES)),
        S3.remote(expand(S3PATH+'Final-Output/{sample}/{sample}_getitd.zip',sample=SAMPLES)),
        S3.remote(S3PATH+'Final-Output/Coverview-Report.xlsx'),
        'DECON/bamfiles.txt'
        #expand('{sample}.analysis/variants/{sample}.strelka-somatic.vcf', sample=SAMPLES),

rule preprocess_reads:
    input:
       r1 = lambda wildcards: glob('sequences/{sample}_S*_L001_R1_001.fastq.gz'.format(sample=wildcards.sample)),
       r2 = lambda wildcards: glob('sequences/{sample}_S*_L001_R2_001.fastq.gz'.format(sample=wildcards.sample)),
       smmip_adaptors = config["adaptors"]
    output:
        r1_trim=temp('{sample}.analysis/processed_reads/{sample}.R1.trimmed.fastq'),
        r2_trim=temp('{sample}.analysis/processed_reads/{sample}.R2.trimmed.fastq')
    log:
        "logs/{sample}/process_reads.log"
    message:
        "Preprocessing reads for {wildcards.sample}"
    shell:
          '{config[EA_UTILS_PATH]}/fastq-mcf -o {output.r1_trim} -o {output.r2_trim} -l 53 -k 0 -q 0 {input.smmip_adaptors} {input.r1} {input.r2}  &>> {log}'

rule gzip:
        input:
          r1_zip='{sample}.analysis/processed_reads/{sample}.R1.trimmed.fastq',
          r2_zip='{sample}.analysis/processed_reads/{sample}.R2.trimmed.fastq'

        output:
            temp('{sample}.analysis/processed_reads/{sample}.R1.trimmed.fastq.gz'),
            temp('{sample}.analysis/processed_reads/{sample}.R2.trimmed.fastq.gz')
        message:
          "Zipping preprocessed reads for {wildcards.sample}"
        run:
          shell('gzip -f {input.r1_zip}')
          shell('gzip -f {input.r2_zip}')
          shell('wait')

rule pair_assembly:
            input:
                forward = '{sample}.analysis/processed_reads/{sample}.R1.trimmed.fastq.gz',
                reverse = '{sample}.analysis/processed_reads/{sample}.R2.trimmed.fastq.gz',
            output:
                temp('{sample}.analysis/assembled_reads/{sample}.assembled.fastq.gz')
            message:
                "Assemble preprocessed reads for {wildcards.sample}"
            threads: 15
            params:
                minlength = '53',
                interopt = '{sample}.analysis/assembled_reads/{sample}'
            log:
                "logs/{sample}/assemble_reads.log"
            run:
                shell('{config[PEAR_PATH]} -f {input.forward} -r {input.reverse} -o {params.interopt} -n {params.minlength} -j {threads} &>> {log}')
                shell('gzip {wildcards.sample}.analysis/assembled_reads/{wildcards.sample}.*.fastq')
                shell('wait')


rule rename_smips_reads:
        input:
            '{sample}.analysis/assembled_reads/{sample}.assembled.fastq.gz'
        output:
            temp('{sample}.analysis/assembled_reads/{sample}.assembled.fastq')
        message:
            "Rename assembled reads for {wildcards.sample} using perl script"
        log:
            "logs/{sample}/assemble_reads.log"
        run:
            shell('gunzip {input}')
            shell('perl ~/hemoseq_v2/mutation_detector/snakemake-pipeline/scripts/rename_smips_reads.pl {wildcards.sample}.analysis/assembled_reads/{wildcards.sample}.assembled.fastq &>> {log}')
            shell('mv {wildcards.sample}.analysis/assembled_reads/{wildcards.sample}.assembled.fastq.tags {output}')

rule mapping_reads:
        input:
            assembledfastq = '{sample}.analysis/assembled_reads/{sample}.assembled.fastq',
            genome = config["reference_genome"]
        output:
            temp('{sample}.analysis/mapped_reads/{sample}.sam')
        params:
            rg = "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:{sample}\\tPI:200"
        threads:
            config["cores"]
        benchmark:
            "benchmarks/{sample}.bwa.benchmark.txt"
        message:
            "Align reads for {wildcards.sample} using bwa"
        log:
            "logs/{sample}/bwa.log"
        run:
            shell('gzip {input.assembledfastq}')
            shell("bwa mem -R '{params.rg}' -t {threads} {input.genome} {input.assembledfastq}.gz > {output} 2>> {log}")

rule sam_conversion:
                input:
                    sam = '{sample}.analysis/mapped_reads/{sample}.sam',
                    genome = config["reference_genome"]
                output:
                    temp('{sample}.analysis/aligned_reads/{sample}.fxd_sorted.bam')
                message:'Sam tools processing for {wildcards.sample}'
                log:
                    "logs/{sample}/samtools.log"
                run:
                    shell('java -Xmx10G -jar $PICARD FixMateInformation I= {input.sam} O= {wildcards.sample}.analysis/mapped_reads/{wildcards.sample}.fxd.sam VALIDATION_STRINGENCY=SILENT TMP_DIR=tmp &>> {log}')
                    shell('samtools view -bT {input.genome} {wildcards.sample}.analysis/mapped_reads/{wildcards.sample}.fxd.sam > {wildcards.sample}.analysis/aligned_reads/{wildcards.sample}.fxd.bam')
                    shell('samtools sort {wildcards.sample}.analysis/aligned_reads/{wildcards.sample}.fxd.bam > {output}')
                    shell('samtools index {output}')

rule RealignerTargetCreator:
            input:
                bam = '{sample}.analysis/aligned_reads/{sample}.fxd_sorted.bam',
                refgenome = config["reference_genome"],
                site1 = config["known-variants"]["site1"]
            output:
                temp('{sample}.analysis/gatk38_processing/{sample}.intervals')
            threads:
                    config["cores"]
            benchmark:
                "benchmarks/{sample}.gatk-process.txt"
            log:
                "logs/{sample}/gatk-process.log"
            message:
                "Gatk Preprocessing for {wildcards.sample} realign targets"
            shell:
                'java -Xmx10G -jar {config[GATK38_PATH]} -T RealignerTargetCreator -R {input.refgenome} -nt {threads} -I {input.bam} --known {input.site1} -o {output} &>> {log}'

rule IndelRealigner:
            input:
                bam = '{sample}.analysis/aligned_reads/{sample}.fxd_sorted.bam',
                refgenome = config["reference_genome"],
                site1 = config["known-variants"]["site1"],
                targetIntervals = '{sample}.analysis/gatk38_processing/{sample}.intervals'
            output:
                temp('{sample}.analysis/gatk38_processing/{sample}.realigned.bam')
            log:
                "logs/{sample}/gatk-process.log"
            benchmark:
                "benchmarks/{sample}.gatk-process.txt"
            message:
                "Gatk Preprocessing for {wildcards.sample} realign indels"
            shell:
                'java -Xmx10G -jar {config[GATK38_PATH]} -T IndelRealigner -R {input.refgenome} -I {input.bam} -known {input.site1} --targetIntervals {input.targetIntervals} -o {output} &>> {log}'

rule BaseRecalibrator:
                input:
                    bam = '{sample}.analysis/gatk38_processing/{sample}.realigned.bam',
                    refgenome = config["reference_genome"],
                    site1 = config["known-variants"]["site2"],
                    site2 = config["known-variants"]["site3"]
                log:
                    "logs/{sample}/gatk-process.log"
                output:
                    temp('{sample}.analysis/gatk38_processing/{sample}.recal_data.table')
                benchmark:
                    "benchmarks/{sample}.gatk-process.txt"
                message:
                    "Gatk Preprocessing for {wildcards.sample} base recalibration"
                shell:
                    "java -Xmx10G -jar {config[GATK38_PATH]} -T BaseRecalibrator -R {input.refgenome}"
                    " -I {input.bam} -knownSites {input.site1} -knownSites {input.site2}"
                    " -o {output} &>> {log}"

rule PrintReads:
                input:
                    bam = '{sample}.analysis/gatk38_processing/{sample}.realigned.bam',
                    bqsr_table = '{sample}.analysis/gatk38_processing/{sample}.recal_data.table',
                    refgenome = config["reference_genome"]
                output:
                    temp('{sample}.analysis/gatk38_processing/{sample}.aligned.recalibrated.bam')
                benchmark:
                    "benchmarks/{sample}.gatk-process.txt"
                log:
                    "logs/{sample}/gatk-process.log"
                message:
                    "Gatk Preprocessing for {wildcards.sample} print reads"
                shell:
                    "java -Xmx10G -jar {config[GATK38_PATH]} -T PrintReads -R {input.refgenome} -I {input.bam}"
                    " --BQSR {input.bqsr_table} -o {output} &>> {log}"

rule generatefinalbam:
                input:
                    '{sample}.analysis/gatk38_processing/{sample}.aligned.recalibrated.bam'
                output:
                    bam = '{sample}.analysis/gatk38_processing/{sample}.final.bam',
                    bai = '{sample}.analysis/gatk38_processing/{sample}.final.bam.bai'
                log:
                    "logs/{sample}/gatk-process.log"
                message:
                    "Gatk Preprocessing for {wildcards.sample} final bam"
                run:
                    shell('samtools sort {input} > {output.bam}')
                    shell('samtools index {output.bam}')

#Qualimap
rule quality_check:
    input:
        '{sample}.analysis/gatk38_processing/{sample}.final.bam'
    params:
        outdir = '{sample}.analysis/qualimap_results',
        bedfile = config["bedfile"].replace('_file.bed', '_qualimap.bed')
    output:
        '{sample}.analysis/qualimap_results/{sample}.qc_report.pdf'
    log:
        "logs/{sample}/quality-check.log"
    message:
        "bam Quality check for {wildcards.sample}"
    shell:
        "{config[QUALIMAP_PATH]}/qualimap bamqc -bam {input}"
        " -gff {params.bedfile} -c"
        " -outdir {params.outdir} -outfile {wildcards.sample}.qc_report.pdf"

include: "rules/mutect2.smk"
include: "rules/freebayes.smk"
include: "rules/platypus.smk"
include: "rules/vardict.smk"
include: "rules/varscan.smk"
include: "rules/getITD.smk"
include: "rules/strelka.smk"
include: "rules/lofreq.smk"
include: "rules/coverage.smk"
include: "rules/somaticseq.smk"
include: "rules/annotation.smk"
include: "rules/formatannotation.smk"
include: "rules/finaloutput.smk"
include: "rules/combinevariants.smk"
include: "rules/coverview.smk"
include: "rules/cava.smk"
