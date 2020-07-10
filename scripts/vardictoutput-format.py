# -*- coding: utf-8 -*-
"""
Created on Sat May  9 00:35:17 2020

@author: YASH
"""

import pandas as pd
filename = snakemake.input[0]
#filename = 'E:\\Bionfotmatics Projects\\NGS\\OCIAMl\\OCIAMl.annovar\\OCIAMl3.vardict.hg19_multianno.csv'
df = pd.read_csv(filename)
x = df['Otherinfo1']
discarded_column=df.columns.get_loc('Otherinfo2')
data = dict()

data.setdefault('REF_COUNT', [])
data.setdefault('ALT_COUNT', [])
data.setdefault('VAF', [])
for row in x:
    rowitems=row.split('\t')
    formatval=rowitems[-1].split(':')
    readdepth=formatval[3]
    vaf=float(formatval[4])
    readslist=list(map(int,readdepth.split(',')))
    data['REF_COUNT'].append(readslist[0])
    data['ALT_COUNT'].append(readslist[1])
    data['VAF'].append("{:.2%}".format(vaf))


df1=df.iloc[:,:5]
df2=pd.DataFrame(data, columns=data.keys())
df3=df.iloc[:,5:discarded_column]

horizontal_stack = pd.concat([df1, df2, df3], axis=1)
horizontal_stack.rename(columns = {'Func.refGene':'Variant Site', 'ExonicFunc.refGene':'Variant Function'}, inplace = True) 
#horizontal_stack.to_csv('E:\\Bionfotmatics Projects\\NGS\\OCIAMl\\OCIAMl.annovar\\OCIAMl3.vardict.modified.csv', index=False)
horizontal_stack.to_csv(snakemake.output[0], index=False)

