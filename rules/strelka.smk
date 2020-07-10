rule strelka2_germline:
    input:
        '{sample}.analysis/gatk38_processing/{sample}.final.bam'
    output:
        '{sample}.analysis/variants/{sample}.strelka.vcf'
    message: " Call strelka germline variants for {wildcards.sample}"
    log:
        "logs/{sample}/strelka.log"
    params:
        outdir='{sample}.analysis/variants/strelka'
    run:
        shell('~/programs/strelka2/bin/configureStrelkaGermlineWorkflow.py --bam {input} --referenceFasta {config[reference_genome]} --callRegions  {config[bedfile]}.gz --targeted --runDir {params.outdir}')
        shell('{params.outdir}/runWorkflow.py -m local -j 20 2> {log}')
        shell('gunzip -f {params.outdir}/results/variants/variants.vcf.gz')
        shell('mv {params.outdir}/results/variants/variants.vcf {output}')

rule strelka_somatic:
    input:
        '{sample}.analysis/gatk38_processing/{sample}.final.bam'
    output:
        '{sample}.analysis/variants/strelka-somatic/results/variants/somatic.indels.vcf.gz',
        '{sample}.analysis/variants/strelka-somatic/results/variants/somatic.snvs.vcf.gz'
    message: " Call strelka somatic variants for {wildcards.sample}"
    log:
        "logs/{sample}/strelka-somatic.log"
    params:
        outdir='{sample}.analysis/variants/strelka-somatic'
    run:
        shell('~/programs/strelka2/bin/configureStrelkaSomaticWorkflow.py --normalBam {config[NA12878_BAM]} --tumorBam {input} --referenceFasta {config[reference_genome]} --callRegions  {config[bedfile]}.gz --targeted --runDir {params.outdir}')
        shell('{params.outdir}/runWorkflow.py -m local -j 20 2> {log}')


rule strelkasom_variant:
    input:
        '{sample}.analysis/variants/strelka-somatic/results/variants/somatic.indels.vcf.gz',
        '{sample}.analysis/variants/strelka-somatic/results/variants/somatic.snvs.vcf.gz'
    output:
        '{sample}.analysis/variants/{sample}.strelka-somatic.vcf'
    shell:
        'bcftools concat -a {input} -o {output}'
