# RAAVioli
_Recombinant Adeno-Associated Viral IntegratiOn anaLysIs (RAAVIoli)_

Bioinformatics pipeline for the identification and characterization of AAV integration sites and viral rearrangements.

### Install
```
chmod +x mamba_setup.sh
./mamba_setup.sh
```
it will install all the required tools in the bin sub-directory.
A *config.txt* with all needed paths will be created. 
In case of errors please check the *error.log files to see whetever an error has occured during the installation of the tools needed.

## Usage
Run the pipeline from the installation dir.  
**The Viral genome must be a single sequence named chrV. This is a must for later steps.**  
The entry point is a .tsv file containing info about the pools and as last column the path to the fastq.gz file. It can contain multiple files. E.g.
```
sample	sample_label	sample_name	sample_source	path
sample1	InVivo1	B01	B01	path/to/sample1.fastq.gz
sample2	InVivo2	B02	B02	path/to/sample2.fastq.gz
```
This file is included into the repository. You can modify it or copy it elsewhere.  

Sample usage:
```
./RAAVioli_Long.sh -i sample_label.tsv -t threads -v viral_genome.fa -r reference.fa
               -R 1 -a annotation.gtf -o output_dir -m mixed_genome.fa -c variables_mixed -w variables_viral -y variables_rscript 


Sample usage: ./RAAVioli.sh -i sample_label.tsv -t threads -v viral_genome.fa -r reference.fa\
                            -R 1 -a annotation.gtf -o output_dir -m mixed_genome.fa -c variables_mixed -w variables_viral -y variables_rscript 


	-i the .tsv file with the paths to fastq files as last column.

	-t max threads to be used.
    -c path to variables_mixed
    -w path to variables_viral
    -y path to variables_rscript 
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

	-a the gtf file with the custom annotation. 

	-o path to the output directory.


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
**-a** :the gtf file with the custom annotation (containing both reference and viral annotations). It must be sorted  
You can sort it using the command 
```
sort -k1,1 -k4,4n -k5,5n
```
**-o**: path to the output directory  

### Other parameters
The first step of the pipeline is aligning the reads to the viral genome. This is done through bwa mem. You can change the bwa mem parameters and the samtools filter parameters changing the values in the viral_variables file. In the second step the reads will be aligned to the mixed genome. Even in this case you can change aligning and filtering parameters in the mixed_variables file or reference_variables file. You can also modify the species (it will be used in the files names and tags).
## Examples
Suppose you have the **viral genome** (in */resources/viral/aav.fa*) and the **reference genome** (in */resources/human/hg19.fa*) but **no indexes and no mixed genome**.  
You can run the following:  
```
RAAVioli.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.gtf -v /resources/viral/aav.fa\ 
-r /resources/human/hg19.fa
```
If you have the **bwa-index of the hg19.fa** in *the same directory of the genome* you can run:  
```
RAAVioli.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.gtf -v /resources/viral/aav.fa\ 
             -r /resources/human/hg19.fa -R 1
```
or  
```
RAAVioli.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.gtf -v /resources/viral/aav.fa\ 
               -r /resources/human/hg19.fa -R /resources/human/hg19.fa
```
Now suppose you have also the **bwa-index of the aav.fa** but in a different directory (*/resources/viral/index/aav.fa*). You can run:  
```
RAAVioli.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.gtf -v /resources/viral/aav.fa\ 
             -V /resources/viral/index/aav.fa -r /resources/human/hg19.fa -R 1
```
Please note that if you did not create the mixed genome before you have to specify both the viral genome and the reference genome.  
Instead if you have already created the **mixed genome** you do not need to specify the -v viral genome and the -r reference genome.
E.g. you have the mixed genome in */resources/mixed/mixed.fa* and you have both the **viral index** (in */resources/viral/index/aav.fa*)  and the 
**reference index** (in **/resources/human/hg19.fa**). You can use:  
```
RAAVioli.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.gtf -V /resources/viral/index/aav.fa\ 
             -R /resources/human/hg19.fa -m /resources/mixed/mixed.fa
```
If you also have the **bwa-index of the mixed genome** (in */resources/mixed/index/mixed.fa*) you can use -M instead of -m (in fact if -M is specified, -m will not be considered):  
```
RAAVioli.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.gtf -V /resources/viral/index/aav.fa\ 
                -R /resources/human/hg19.fa -M /resources/mixed/index/mixed.fa
```
If instead you do have the mixed_genome and its index but not the Viral index you can use:  
```
RAAVioli.sh -i sample_label.tsv -t 4 -o outputDir -a annotation.gtf -v /resources/viral/aav.fa\ 
                -R /resources/human/hg19.fa -M /resources/mixed/index/mixed.fa
```
The same thing applies for the reference genome.
