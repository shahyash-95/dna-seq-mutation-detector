# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 16:47:19 2020

@author: YASH
"""

import pandas as pd

file1 = snakemake.input[0]
file2 = snakemake.input[1]
outfile = snakemake.output[0]

df1 = pd.read_csv(file1, delimiter="\t")

df2= pd.read_csv(file2, delimiter="\t")

concat= pd.concat([df1,df2],axis=0)

concat.drop(columns=['ID'], inplace=True)
concat.replace(to_replace='.', value='-1', inplace=True)
concat.to_csv(outfile, index=False)