rule format_platypus:
    input:
        '{sample}.analysis/ANNOVAR/{sample}.platypus.hg19_multianno.csv'
    output:
        '{sample}.analysis/Annovar_Modified/{sample}.platypus.csv'
    message: " Modify annotation format for {wildcards.sample} platypus"
    script:
        '../scripts/platypusoutput-format.py'


rule format_freebayes:
    input:
        '{sample}.analysis/ANNOVAR/{sample}.freebayes.hg19_multianno.csv'
    output:
        '{sample}.analysis/Annovar_Modified/{sample}.freebayes.csv'
    message: " Modify annotation format for {wildcards.sample} freebayes"
    script:
        '../scripts/freebayesoutput-format.py'

rule format_somaticseq:
    input:
        '{sample}.analysis/ANNOVAR/{sample}.somaticseq.hg19_multianno.csv'
    output:
        '{sample}.analysis/Annovar_Modified/{sample}.somaticseq.csv'
    message: " Modify annotation format for {wildcards.sample} somaticseq"
    script:
        '../scripts/somaticseqoutput-format.py'

rule format_vardict:
    input:
        '{sample}.analysis/ANNOVAR/{sample}.vardict.hg19_multianno.csv'
    output:
        '{sample}.analysis/Annovar_Modified/{sample}.vardict.csv'
    message: " Modify annotation format for {wildcards.sample} vardict"
    script:
        '../scripts/vardictoutput-format.py'

#freebayes and platypus variants combined --merged
rule format_combined:
    input:
        '{sample}.analysis/ANNOVAR/{sample}.combined.hg19_multianno.csv'
    output:
        '{sample}.analysis/Annovar_Modified/{sample}.combined.csv'
    message: " Modify annotation format for {wildcards.sample} platypus and freebayes combined output"
    script:
        '../scripts/combineoutput-format.py'
