
rule freebayes_variant:
    input:
        '{sample}.analysis/gatk38_processing/{sample}.final.bam'
    output:
        '{sample}.analysis/variants/{sample}.freebayes.vcf'
    message: " Call freebayes variants for {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}.variantcalls.benchmark.txt"
    log:
        "logs/{sample}/freebayes.log"
    shell:
        '~/programs/freebayes -f {config[reference_genome]} -b {input} -t {config[bedfile]} > {output} 2>> {log}'
