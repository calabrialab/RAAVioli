# RAAVioli
Bioinformatics pipeline for AAV integration sites and viral rearrangements
### Install
Run ./setup.sh it will check if all requirements are installed otherwise it will try to install the missing ones in the bin directory.
A *config.txt* with all needed paths will be created. We suggest to add the installation dir to the system path so you can run the pipeline from anywhere on your system. You can do this using:  
```
export PATH:$PATH:/path/to/RAAVioli
```
(note that this is not permanently)
For the same reason we suggest to specify all paths to parameters.
## Usage
**The Viral genome must be a single sequence named chrV. This is a must for later steps so please do it.**
Sample usage
```
./pipeline.sh -i sample_label.tsv -t threads -v viral_genome.fa -r reference.fa
               -R 1 -A annotation.bed -o output_dir -m mixed_genome.fa


Sample usage: ./pipeline.sh -i sample_label.tsv -t threads -v viral_genome.fa -r reference.fa\
                            -R 1 -A annotation.bed -o output_dir -m mixed_genome.fa


	-i the .tsv file with the paths to fastq files as last column.

	-t max threads to be used.

	-v (optional) the fasta file with viral genome (e.g. AAV).
	   The bwa-index will be created in the same directory. 
	   If you have already an index please see -V.
	   You must specify -V if you don't specify -v.

	-r (optional) the fasta file with the reference genome (e.g. hg19).
	   The bwa-index will be created in the same directory. 
	   If you have already a bwa-index please see -R. N.B.
	   You must specify -R if you don't specify -r.

	-V (optional) path to the viral bwa-index with basename
	   (e.g. if you have the index in /home/resources/genome/index
	   directory and it has as basename aav.fa
	   you have to specify home/resources/genome/index/aav.fa ).
	   If specified the index of the viral genome will not be made.
	   If you don't specify -V you must specify -v.

	-R (optional) path to the reference bwa-index with basename
	   (e.g. if you have the index in /home/resources/genome/index
	   directory and it has as basename hg19.fa
	   you have to specify home/resources/genome/index/hg19.fa ).
	   If specified the index of the reference genome will not be made.
	   If you don't specify -R you must specify -r.

	-m (optional) the fasta file with the mixed genome
	   N.B. viral genome must be appended at the end of reference genome
	   with the sequence name chrV.
	   Please note that if not specified it will be created and 
	   you must specify -v and -r 
	   (since index could be located in a different dir
	   and to create the mixed genome both genomes are needed). 
	   In this case if you already have
	   the viral index and/or the reference index 
	   in the same directory you can specify -V 1 and/or -R 1 instead 
	   of specifying twice the same path for -v and -V (or -r and -R).

	-M (optional) bwa-index of the mixed_genome.
	   If specified you can omit -m.

	-a the bed file with the custom annotation.

	-o path to the output directory.

	-b (optional) any value. If specified also 2820 will be made.
```
### Options
**-i**: The file .tsv with the paths to fastq files as last column  
 This file must have the same format of the sample_label.tsv. Specify in the last column the path to the already trimmed fastq file(s).  
**-t**: max threads to be used for file. NB since the pipeline is runned in parallel for the different fastq files specified with the -i option, each subprocess will use -t THREADS.  
**-v**: (optional, see -V) path to the fasta file with viral genome (e.g. AAV). **The viral genome must be a single sequence named chrV. This is essential for later steps.** If you have already made the index of the fasta file you can use the -V option otherwise the index will be created in the same directory of the fasta file. You can specify only one between -v or -V or both. -v is required only if you do not specify -m or -M since both viral genome and reference genome are needed.  
**-V**: (optional) path to the viral bwa-index with basename. If specified the index of the viral genome will not be made. If you have already created the mixed genome (see -m and -M) you can specify only -V (e.g. -V /path/to/index.fa -M /path/to/mixed.index.fa). If the index is the same directory of the viral genome fasta file you can use -V 1 e.g. `-v /path/to/viral.fa -V 1`. More example in the **example** section.  
**-r**: (optional, see -R) the fasta file with the reference genome (e.g. hg19). If you have already made the index of the fasta file you can use the -R option. Otherwise the index will be created in the same directory of the fasta file. You can specify only one between -r or -R or both. -r is required only if you do not specify -m or -M since both viral genome and reference genome are needed.  
**-R**: (optional) path to the reference bwa-index with basename. If specified the index of the reference genome will not be made (see -V for more info).  
**-m**: (optional) the fasta file with the mixed genome (viral genome must be appended at the end with the sequence name chrV).
 You can create this file running the following commands: 
```
cp reference.fa mixed.fa
cat viral.fa >> reference.fa
```
If the option is not specified the mixed genome will be created and indexed by the program.  
**-M**: (optional) path to bwa-index of the mixed genome. If specified the index of the mixed genome will not be made. If you specify -M you do not need to specify -m.  
**-a** :the bed file with the custom annotation  
**-o**: path to the output directory  
**-b**: (optional) any value. If specified also 2820 will be made.  
### Other parameters
The first step of the pipeline is aligning the reads to the viral genome. This is done through bwa mem. You can change the bwa mem parameters and the samtools filter parameters changing the values in the viral_variables file. In the second step reads will be aligned to the mixed genome and optionally (if -b is used) to the reference genome. Even in this case you can change aligning and filtering parameters in the mixed_variables file or reference_variables file. You can also modify the species (it will be used in the files names and tags).
## Examples
Suppose you have the **viral genome** (in */resources/viral/aav.fa*) and the **reference genome** (in */resources/human/hg19.fa*) but **no indexes and no mixed genome**.  
You can run the following:  
```
pipeline.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.bed -v /resources/viral/aav.fa\ 
-r /resources/human/hg19.fa
```
If you have the **bwa-index of the hg19.fa** in *the same directory of the genome* you can run:  
```
pipeline.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.bed -v /resources/viral/aav.fa\ 
             -r /resources/human/hg19.fa -R 1
```
or  
```
pipeline.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.bed -v /resources/viral/aav.fa\ 
               -r /resources/human/hg19.fa -R /resources/human/hg19.fa
```
Now suppose you have also the **bwa-index of the aav.fa** but in a different directory (*/resources/viral/index/aav.fa*). You can run:  
```
pipeline.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.bed -v /resources/viral/aav.fa\ 
             -V /resources/viral/index/aav.fa -r /resources/human/hg19.fa -R 1
```
Please note that if you did not create the mixed genome before you have to specify both the viral genome and the reference genome.  
Instead if you have already created the **mixed genome** you do not need to specify the -v viral genome and the -r reference genome.
E.g. you have the mixed genome in */resources/mixed/mixed.fa* and you have both the **viral index** (in */resources/viral/index/aav.fa*)  and the 
**reference index** (in **/resources/human/hg19.fa**). You can use:  
```
pipeline.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.bed -V /resources/viral/index/aav.fa\ 
             -R /resources/human/hg19.fa -m /resources/mixed/mixed.fa
```
If you also have the **bwa-index of the mixed genome** (in */resources/mixed/index/mixed.fa*) you can use -M instead of -m (in fact if -M is specified, -m will not be considered):  
```
pipeline.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.bed -V /resources/viral/index/aav.fa\ 
                -R /resources/human/hg19.fa -M /resources/mixed/index/mixed.fa
```
If instead you do have the mixed_genome and its index but not the Viral index you naturally can use:  
```
pipeline.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.bed -v /resources/viral/aav.fa\ 
                -R /resources/human/hg19.fa -M /resources/mixed/index/mixed.fa
```
The same thing applies for the reference genome.
