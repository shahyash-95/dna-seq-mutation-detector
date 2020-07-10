# -*- coding: utf-8 -*-
"""
Created on Mon May 18 13:38:48 2020

@author: YASH
"""

import pandas as pd
import os
#filepath='E:\\Bionfotmatics Projects\\NGS\\OCIAMl\\'
csvfilenames=snakemake.input
#csvfilenames=[filepath+'OCIAMl3.somaticseq.modified.csv',filepath+'OCIAMl3.combined.csv']
#writer = pd.ExcelWriter(filepath+'OCIAMl3.variants.xlsx')
writer = pd.ExcelWriter(snakemake.output[0]) # output name
for csvfilename in csvfilenames:
    sheetname=os.path.split(csvfilename)[1]
    df = pd.read_csv(csvfilename)
    print('process file:', csvfilename, 'shape:', df.shape)
    df.to_excel(writer,sheet_name=os.path.splitext(sheetname)[0], index=False)
writer.save()