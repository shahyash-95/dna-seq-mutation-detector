rule coverview:
    input:
        '{sample}.analysis/gatk38_processing/{sample}.final.bam'
    output:
        'Coverview/{sample}/{sample}.coverview_regions.txt'
    message: " Generate coverview for {wildcards.sample}"
    log:
        "logs/{sample}/coverview.log"
    params:
        outprefix = 'Coverview/{sample}/{sample}.coverview'

    run:
        shell("~/programs/CoverView-1.4.4/coverview -i {input} -b {config[bedfile]} -c ~/programs/CoverView-1.4.4/config/config.txt -o {params.outprefix} 2> {log}")

rule regions_csv:
    input:
        'Coverview/{sample}/{sample}.coverview_regions.txt'
    output:
        'Coverview/{sample}/{sample}.coverview_regions.csv'
    message: " Get coverview regions file for {wildcards.sample}"
    script:
        '../scripts/coverview.py'

rule coverview_report:
    input:
        expand('Coverview/{sample}/{sample}.coverview_regions.csv', sample=SAMPLES)
    output:
        report('Final-Output/Coverview-Report.xlsx', category='Coverview Report')
    script:
        '../scripts/coverview-report.py'

rule s3_coverview:
    input:
        'Final-Output/Coverview-Report.xlsx'
    output:
        S3.remote(S3PATH+'Final-Output/Coverview-Report.xlsx')
    shell:
        "cp {input} {output}"
"""
rule all:
    input:
        directory(expand('{sample}.analysis/coverview',sample=samples))
rule move:
        input:
            '{sample}.coverview_profiles.txt',
            '{sample}.coverview_summary.txt',
            '{sample}.coverview_regions.txt',
            '{sample}.coverview_meta.json'
        output:
            directory('{sample}.analysis/coverview')
        params:
            outdir = '{sample}.analysis/coverview'
        run:

            shell('mv {input} {output}/')
"""
