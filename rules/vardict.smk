rule vardict_variantcall:
    input:
        '{sample}.analysis/gatk38_processing/{sample}.final.bam'
    output:
        '{sample}.analysis/variants/{sample}.vardict.vcf'
    message: " Call vardict variants for {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}.variantcalls.benchmark.txt"
    log:
        "logs/{sample}/vardict.log"
    shell:
        'VarDict -G {config[reference_genome]} -f 0.03 -N {wildcards.sample} -b {input}'
        ' -c 1 -S 2 -E 3 -g 4 {config[bedfile]} | '
        " sed '1d' | teststrandbias.R | var2vcf_valid.pl -N {wildcards.sample} -E -f 0.03 > {output} 2>> {log}"
