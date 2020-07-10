rule somaticseq:
    input:
        bam='{sample}.analysis/gatk38_processing/{sample}.final.bam',
        mutect2='{sample}.analysis/variants/{sample}.mutect2.vcf',
        vardict='{sample}.analysis/variants/{sample}.vardict.vcf',
        varscan='{sample}.analysis/variants/{sample}.varscan.vcf',
        lofreq='{sample}.analysis/variants/{sample}.lofreq.vcf',
        strelka='{sample}.analysis/variants/{sample}.strelka.vcf'
    output:
        '{sample}.analysis/variants/{sample}.somaticseq/Consensus.sSNV.vcf',
        '{sample}.analysis/variants/{sample}.somaticseq/Consensus.sINDEL.vcf'
    message: " Call somaticseq variants for {wildcards.sample}"
    conda:
        "../envs/somaticseq.yaml"
    log:
        "logs/{sample}/somaticseq.log"
    params:
        outdir='{sample}.analysis/variants/{sample}.somaticseq'
    shell:
        ' somaticseq_parallel.py --output-directory {params.outdir} --genome-reference {config[reference_genome]} --inclusion-region {config[bedfile]} --threads {config[cores]} --algorithm xgboost  --dbsnp-vcf  ~/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf'
        ' single --bam-file {input.bam} --mutect2-vcf {input.mutect2}'
        ' --vardict-vcf {input.vardict} --varscan-vcf {input.varscan}'
        ' --lofreq-vcf {input.lofreq} --strelka-vcf {input.strelka}  --sample-name {wildcards.sample} 2> {log}'

rule somaticseq_combinevcf:
    input:
        snv='{sample}.analysis/variants/{sample}.somaticseq/Consensus.sSNV.vcf',
        indel='{sample}.analysis/variants/{sample}.somaticseq/Consensus.sINDEL.vcf'
    output:
        '{sample}.analysis/variants/{sample}.somaticseq.vcf'
    params:
        tempsnv='{sample}.analysis/variants/{sample}.somaticseq/somaticseq_snv.vcf',
        tempindel='{sample}.analysis/variants/{sample}.somaticseq/somaticseq_indel.vcf'
    run:
        shell('grep "^#" {input.snv} > {params.tempsnv}')
        shell('grep -v "^#" {input.snv} | sort -k1,1V -k2,2g >> {params.tempsnv}')
        shell('bgzip -c {params.tempsnv} > {params.tempsnv}.gz')
        shell('bcftools index -t {params.tempsnv}.gz')
        shell('grep "^#" {input.indel} > {params.tempindel}')
        shell('grep -v "^#" {input.indel} | sort -k1,1V -k2,2g >> {params.tempindel}')
        shell('bgzip -c {params.tempindel} > {params.tempindel}.gz')
        shell('bcftools index -t {params.tempindel}.gz')
        shell('bcftools concat -a {params.tempsnv}.gz {params.tempindel}.gz -o {output}')
