#bed tools
rule calculate_coverage:
    input:
        '{sample}.analysis/gatk38_processing/{sample}.final.bam'
    output:
        '{sample}.analysis/coverage/{sample}.counts.bed'
    message: " Calculate coverage for {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}.coverage.benchmark.txt"
    log:
        "logs/{sample}/coverage.log"
    run:
        shell('bedtools bamtobed -i {input} > {wildcards.sample}.analysis/coverage/{wildcards.sample}.bed')
        shell('bedtools coverage -counts -a {config[bedfile]} -b {wildcards.sample}.analysis/coverage/{wildcards.sample}.bed > {output} 2>> {log}')

rule coverage_plot:
    input:
        '{sample}.analysis/coverage/{sample}.counts.bed'
    output:
        '{sample}.analysis/coverage/{sample}.Low_Coverage.png'
    message: " Plot coverage for {wildcards.sample}"
    script:
        '../scripts/coverageplot.py'
