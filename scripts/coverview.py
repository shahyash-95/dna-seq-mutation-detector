# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 00:52:55 2020

@author: YASH
"""

import pandas as pd

#file = 'C:\\Users\\YASH\\Desktop\\Coverview\\19NGS1434.coverview\\19NGS1434.coverview_regions.txt'
#outfile = file.replace('.txt','.csv')

file = snakemake.input[0]
outfile = snakemake.output[0]

regions = pd.read_csv(file, delimiter="\t")
						
regions.rename(columns = {'RC':'Read count', 'MEDCOV':'Median coverage', 'MINCOV':'Minimum coverage',
                          'MEDQCOV':'Median quality coverage','MINQCOV':'Minimum quality coverage',
                          'MAXFLMQ':'Maximum fraction of low mapping quality',
                          'MAXFLBQ':'Maximum fraction of low base quality'}, inplace = True) 

regions.to_csv(outfile, index=False)


"""
### Read Summary File -- generate csv file and bar plot reads in and reads out and put info in plot as text
file1 = 'C:\\Users\\YASH\\Desktop\\19NGS1434.coverview\\19NGS1434\\19NGS1434.coverview_summary.txt'
summary = pd.read_csv(file1, delimiter="\t")
info = summary[:3]
chrom_plt = summary[3:]
chromosomes = list(chrom_plt.iloc[:,0])

## code to generate plot from chrom dataframe
#chrom_plt = chrom_plt.sort_values('RCIN')
reads_list = list(map(int,chrom_plt['RCIN'].tolist()))
samplename='19NGS1434'
fig, ax = plt.subplots(figsize=(20,10))
fig.suptitle(f"Chromosome Level Summary of Mapped Read Counts for {samplename}", fontsize=20, fontweight='bold')
ax.bar(chrom_plt['#CHROM'],reads_list)
ax.set_xticklabels(chrom_plt['#CHROM'], rotation=90, horizontalalignment='right',fontsize='12')
#ax.set_title(f'Number_of_MIPS ={df.shape[0]}')
ax.set_ylabel('Read Counts')
plt.savefig('C:\\Users\\YASH\\Desktop\\19NGS1434.coverview\\19NGS1434\\19NGS1434.coverview.png', bbox_inches='tight')
plt.show()


"""
# Read regions and generate flagged file and flagged summary
"""
flagged = regions.loc[(regions['Pass_or_flag']=='FLAG')]
flag_summary = dict()
flag_summary.setdefault('Chromosome', chromosomes)
flag_summary.setdefault('Total Regions', [])
flag_summary.setdefault('Flagged Regions', [])

for chrom in chromosomes:
    flag_summary['Total Regions'].append(len(regions.loc[(regions['Chromosome']==chrom)]))
    flag_summary['Flagged Regions'].append(len(flagged.loc[(flagged['Chromosome']==chrom)]))

region_sum=pd.DataFrame(flag_summary, columns=flag_summary.keys())
flag_sum = region_sum.loc[(region_sum['Flagged Regions']>0)]

##merge chrom summary, regions, flagged and flagged summary into one excel
writer = pd.ExcelWriter(outfile) # output name
sheet1='Read-Count Summary'
summary.to_excel(writer,sheet_name=sheet1, index=False)
sheet2='Regions Summary'
flag_sum.to_excel(writer,sheet_name=sheet2, index=False)
sheet3='Regions'
regions.to_excel(writer,sheet_name=sheet3, index=False)
sheet4='Flagged Regions'
flagged.to_excel(writer,sheet_name=sheet4, index=False)

writer.save()
"""