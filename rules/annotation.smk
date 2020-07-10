rule convert2annovar:
    input:
        '{sample}.analysis/variants/{sample}.{variant}.vcf'
    output:
        '{sample}.analysis/ANNOVAR/{sample}.{variant}.avinput'
    shell:
        'perl ~/programs/annovar/convert2annovar.pl -format vcf4 {input} --outfile {output}  --withzyg --includeinfo'


rule annotate_annovar:
    input:
        '{sample}.analysis/ANNOVAR/{sample}.{variant}.avinput'
    output:
        '{sample}.analysis/ANNOVAR/{sample}.{variant}.hg19_multianno.csv'
    params:
        outprefix = '{sample}.analysis/ANNOVAR/{sample}.{variant}'
    message: " Annotate variants for sample : {wildcards.sample} variant: {wildcards.variant} using ANNOVAR"
    log:
        "logs/{sample}/{variant}.annovar.log"
    shell:
        "perl ~/programs/annovar/table_annovar.pl {input} ~/programs/annovar/humandb/ \
   --out {params.outprefix} \
   --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all \
   --operation gx,r,f,f,f,f,f \
   --buildver hg19 \
   --nastring '.' \
   --otherinfo \
   --csvout --thread {config[cores]} --xreffile ~/programs/annovar/example/gene_fullxref.txt &>>  {log}"
