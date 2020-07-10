rule sort_variants:
    input:
        '{sample}.analysis/variants/{sample}.freebayes.vcf',
        '{sample}.analysis/variants/{sample}.platypus.vcf'
    output:
        '{sample}.analysis/variants/{sample}.freebayes.sorted.vcf',
        '{sample}.analysis/variants/{sample}.platypus.sorted.vcf'
    run:
        shell('grep "^#" {input[0]} > {output[0]}')
        shell('grep -v "^#" {input[0]} | sort -k1,1V -k2,2g >> {output[0]}')
        shell('grep "^#" {input[1]} > {output[1]}')
        shell('grep -v "^#" {input[1]} | sort -k1,1V -k2,2g >> {output[1]}')


rule combined_variants:
    input:
        '{sample}.analysis/variants/{sample}.freebayes.sorted.vcf',
        '{sample}.analysis/variants/{sample}.platypus.sorted.vcf'
    output:
        '{sample}.analysis/variants/{sample}.combined.vcf'
    shell:
        " java -jar ~/programs/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T CombineVariants "
        " -R {config[reference_genome]} --variant {input[0]} --variant {input[1]}"
        " -o {output} -genotypeMergeOptions UNIQUIFY "
