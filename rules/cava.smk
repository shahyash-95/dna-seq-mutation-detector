rule annotate_CAVA:
    input:
        '{sample}.analysis/variants/{sample}.{variant}.vcf'
    output:
        '{sample}.analysis/CAVA/{sample}.{variant}.txt'
    params:
        outprefix = '{sample}.analysis/CAVA/{sample}.{variant}'
    message: " Annotate variants for sample : {wildcards.sample} variant: {wildcards.variant} using CAVA"
    shell:
          "~/programs/CAVA-1.2.3/cava  -c ~/programs/CAVA-1.2.3/config.txt -t 10 -i {input} -o {params.outprefix}"

rule cava_output:
    input:
        '{sample}.analysis/CAVA/{sample}.combined.txt',
        '{sample}.analysis/CAVA/{sample}.somaticseq.txt'
    output:
        '{sample}.analysis/CAVA/{sample}_CAVA.csv'
    message: " Generate CAVA csv output file for {wildcards.sample}"
    script:
        '../scripts/cava.py'
