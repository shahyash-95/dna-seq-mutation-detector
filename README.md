Steps to install dependencies
1. Miniconda installation
	Downlaod  Miniconda python 3 version for linux
	 bash Miniconda3-latest-Linux-x86_64.sh
1 a. Check for insatllation and config
	conda version
	conda list
	conda config --set auto_activate_base false 

2. Install snakemake using conda
	conda create -c conda-forge -c bioconda -n snakemake snakemake
	conda activate snakemake
	snakemake --help
	
3. Install other tools dependencies 
a. Coverview
b. CAVA
c. ANNOVAR
d. conda install -c bioconda bedtools
e. All variant callers ( use source code or binary installations -- do not use conda installation for lofreq)
f. gatk 3.8 jar
g. qualimap - check website http://qualimap.bioinfo.cipf.es/


