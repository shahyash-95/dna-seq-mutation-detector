# -*- coding: utf-8 -*-
"""
Created on Sat May  9 00:35:17 2020

@author: YASH
"""

import pandas as pd
filename = snakemake.input[0]
df = pd.read_csv(filename)
x = df['Otherinfo1']
discarded_column=df.columns.get_loc('Otherinfo2')
data = dict()
data.setdefault('Variant_Callers', [])
data.setdefault('REF_COUNT', [])
data.setdefault('ALT_COUNT', [])
data.setdefault('VAF', [])
for row in x:
    rowitems=row.split('\t')
    formatval=rowitems[-1].split(':')
    refcount=formatval[-2]
    altcount=formatval[-1]
    if ',' in refcount:
        refcount=refcount.split(',')[0]
    if ',' in altcount:
        altcount=altcount.split(',')[0]
    data['REF_COUNT'].append(refcount)
    data['ALT_COUNT'].append(altcount)
    allele_fraction=float(int(altcount)/int(refcount))    
    vaf="{:.2%}".format(allele_fraction)
    data['VAF'].append(vaf)
    data['Variant_Callers'].append('Platypus')

df1=df.iloc[:,:5]
df2=pd.DataFrame(data, columns=data.keys())
df3=df.iloc[:,5:discarded_column]

horizontal_stack = pd.concat([df1, df2, df3], axis=1)
horizontal_stack.to_csv(snakemake.output[0], index=False)
print(horizontal_stack.shape)
