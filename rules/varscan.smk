rule mpileup:
    input:
        bam = '{sample}.analysis/gatk38_processing/{sample}.final.bam',
        refgenome = config["reference_genome"]
    output:
        '{sample}.analysis/variants/{sample}.mpileup'
    message: "Generate Samtools mpileup  file for {wildcards.sample}"
    log:
        "logs/{sample}/mpielup.log"
    shell:
        'samtools mpileup -f {input.refgenome} {input.bam} > {output} 2>> {log}'

rule varscan_variantcall:
    input:
        '{sample}.analysis/variants/{sample}.mpileup'
    output:
        '{sample}.analysis/variants/{sample}.varscan_snp.vcf',
        '{sample}.analysis/variants/{sample}.varscan_indel.vcf'
    message: "Call varscan variants for {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}.variantcalls.benchmark.txt"
    log:
        "logs/{sample}/varscan.log"
    run:
        shell('java -jar {config[VARSCAN_PATH]} mpileup2snp {input} --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > {output[0]} 2>> {log}')
        shell('java -jar {config[VARSCAN_PATH]} mpileup2indel {input} --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > {output[1]} 2>> {log}')

rule combine_vcf:
    input:
        '{sample}.analysis/variants/{sample}.varscan_snp.vcf',
        '{sample}.analysis/variants/{sample}.varscan_indel.vcf'
    output:
        '{sample}.analysis/variants/{sample}.varscan.vcf'
    run:
        shell('bgzip -c {input[0]} > {input[0]}.gz')
        shell('bgzip -c {input[1]} > {input[1]}.gz')
        shell('bcftools index -t {input[0]}.gz')
        shell('bcftools index -t {input[1]}.gz')
        shell('bcftools concat -a {input[0]}.gz {input[1]}.gz -o {output}')
