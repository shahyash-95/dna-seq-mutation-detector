import pandas as pd
import matplotlib.pyplot as plt

samplename = snakemake.wildcards.sample
filename = snakemake.input[0]
df = pd.read_csv(filename, sep='\t', header=None)
df = df.loc[(df[4] < 100)]
df = df.sort_values(4, ascending=False)

fig, ax = plt.subplots(figsize=(20,10))
fig.suptitle(f"Less than 100X coverage regions for {samplename} Myeloid smMIPS assay", fontsize=20, fontweight='bold')
ax.bar(df[3],df[4])
ax.set_xticklabels(df[3], rotation=90, horizontalalignment='right',fontsize='12')
ax.set_title(f'Number_of_MIPS ={df.shape[0]}')
ax.set_ylabel('Coverage')
plt.savefig(snakemake.output[0], bbox_inches='tight')
plt.show()
