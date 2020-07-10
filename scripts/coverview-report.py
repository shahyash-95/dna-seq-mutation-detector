# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 16:34:57 2020

@author: YASH
"""

import pandas as pd
import os

#file1= 'C:\\Users\\YASH\\Desktop\\Coverview\\20NGS345.coverview\\20NGS345.coverview_regions.csv'
#file2= 'C:\\Users\\YASH\\Desktop\\Coverview\\20NGS123.coverview\\20NGS123.coverview_regions.csv'
#file3= 'C:\\Users\\YASH\\Desktop\\Coverview\\19NGS1434.coverview\\19NGS1434.coverview_regions.csv'
writer = pd.ExcelWriter(snakemake.output[0])
dataframes=[]
commoncols=[]
#files = [file1,file2,file3]
files=snakemake.input
for file in files:    
    path, filename = os.path.split(file)
    samplename=os.path.splitext(filename)[0].split('.')[0]
    df= pd.read_csv(file, index_col=0)
    df=df.iloc[:,:6]
    df.rename(columns = {'Read count':samplename+'_Read count', 'Median coverage':samplename+'_Median coverage', 'Pass_or_flag':samplename+'_FILTER'},inplace=True)
    dataframes.append(df)
    
commoncols = list(dataframes[0].columns)[:3]
for i in range(1,len(dataframes)):
    dataframes[i].drop(commoncols, axis=1, inplace=True)
        
concat=pd.concat(dataframes, axis=1)
concat.to_excel(writer)
writer.save()

#merged = pd.merge(dataframes[0], dataframes[1], on=['#Region','Chromosome', 'Start_position', 'End_position'])

#joined = dataframes[0].join(other=dataframes[1:], how='inner')
#joined.to_csv('Coverview-Report.csv', index=None)                 
