
rule lofreq_variant:
    input:
        '{sample}.analysis/gatk38_processing/{sample}.final.bam'
    output:
        '{sample}.analysis/variants/{sample}.lofreq.vcf'
    message: " Call lofreq variants for {wildcards.sample}"
    log:
        "logs/{sample}/lofreq.log"
    shell:
        'lofreq call -f {config[reference_genome]} -o  {output} {input} -l {config[bedfile]} --call-indels -s -S {config[known-variants][site2]}.gz 2> {log}'
