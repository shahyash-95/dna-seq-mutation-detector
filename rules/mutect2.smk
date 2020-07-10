rule mutect2_variantcall:
    input:
        '{sample}.analysis/gatk38_processing/{sample}.final.bam'
    output:
        '{sample}.analysis/variants/{sample}.mutect2.vcf'
    message: " Call MUTECT2 variants for {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}.mutect2.benchmark.txt"
    log:
        "logs/{sample}/MuTect2.log"
    shell:
        "java -Xmx10G -jar {config[GATK38_PATH]} -T MuTect2 -R {config[reference_genome]}"
         " -I:tumor {input} -o {output}"
         " --dbsnp {config[known-variants][site2]} -L {config[bedfile]}"
         " -nct 30 -contamination 0.02 -mbq 30 &>> {log}"
