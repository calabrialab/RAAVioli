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

The simplest way to run RAAVioli_long is:
```
./RAAVioli_long.sh -i sample_label.tsv -t threads -V viral_genome.fa -M mixed_genome.fa \
                              -a annotation.gtf -o output_dir -c variables_mixed -w variables_viral -y variables_rscript 
    
    -i the .tsv file with the paths to fastq files as last column.
    -t max threads to be used.
    -c path to variables_mixed
    -w path to variables_viral
    -y path to variables_rscript 

    -V path to the viral bwa-index with basename
       (e.g. if you have the index in /home/resources/genome/index
       directory and it has as basename aav.fa
       you have to specify home/resources/genome/index/aav.fa ).
       If specified the index of the viral genome will not be made.
       If you don't specify -V you must specify -v.
       
    -M path to bwa-index of the mixed_genome.
       The mixed genome is obtained by appending the viral genome at 
       the end of reference genome with the sequence name chrV.
    
    -a the gtf file with the custom annotation. 
    
    -o path to the output directory.


```
___If you don't have indexes or mixed genome, they can be created using other options (-r, -v -m), see full options below.___



### Other parameters
The first step of the pipeline is aligning the reads to the viral genome. This is done through bwa mem. You can change
the bwa mem parameters and the samtools filter parameters changing the values in the viral_variables file. 
In the second step the reads will be aligned to the mixed genome. 
Even in this case you can change aligning and filtering parameters in the mixed_variables file or reference_variables 
file. You can also modify the species (it will be used in the files names and tags).
#### viral_variables
The following variables control the behavior of BWA-MEM during the initial viral-genome alignment step.
___This step is intentionally permissive: we aim to retain any read showing even a weak signal of alignment
to the viral genome before stricter downstream filtering.___

```
bwa_mem_k=18       # minimum seed length (k-mer size)
bwa_mem_r=2        # seed occurrence threshold
bwa_mem_A=1        # matching score
bwa_mem_T=10       # minimum score to output alignment
bwa_mem_d=10000    # maximum occurrences for a seed
bwa_mem_B=2        # mismatch penalty
bwa_mem_O=1        # gap open penalty
bwa_mem_E=1        # gap extension penalty
bwa_mem_L=0        # clipping penalty
sam_view_q=6       # MAPQ threshold when filtering SAM output
bam_filter_AS=50   # minimum alignment score (AS tag) to retain a read
```

All parameters correspond directly to BWA-MEM command-line options (expect for sam_view_q referring to samtools
and bam_filter_AS referring to bamtools).
Users may adjust them depending on read length, expected divergence, or stringency requirements, but I suggest
to operate more on the variables_mixed or variables_rscript.
The official BWA-MEM documentation describing each option is available at:

https://github.com/lh3/bwa/blob/master/bwa.1

#### mixed_variables
he following parameters control BWA-MEM during alignment against the mixed host + vector genome.
This step is intentionally more stringent than the initial viral-only alignment.
At this stage, the goal is to ensure that reads retained for downstream classification map unambiguously and with high confidence to either the host or vector genome.
```
bwa_mem_k=18       # seed length
bwa_mem_r=2        # seed occurrence threshold
bwa_mem_A=1        # match score
bwa_mem_T=15       # minimum score to output alignment
bwa_mem_d=100      # drop chains with excessive hits
bwa_mem_B=4        # mismatch penalty
bwa_mem_O=5        # gap open penalty
bwa_mem_E=5        # gap extension penalty
bwa_mem_L=0        # clipping penalty
sam_view_q=16      # MAPQ threshold for filtering
bam_filter_AS=800  # minimum alignment score required THIS CAN BE VERY STRINGENT
bam_filter_NM=800  # maximum allowed edit distance (NM tag)
SPEC="mixed"       # classification label used in the output files
```
These settings are tuned to:

- increase stringency relative to the viral-genome prescreening step, favoring high-quality, full-length alignments

- reduce ambiguous or low-confidence placements in the mixed genome and enforce strong penalties for mismatches and gaps (B, O, E)

- apply stricter MAPQ and alignment-score filters  and reduce spurious hits in repetitive regions (lower d)

Practically, this improves the specificity of integration-site calling.
Depending on the source pfo your reads you can change the parameters as suggested by the
-x parameter of bwa mem.

___AS score is setted very high since we worked with PacBio high fidelity reads___
you can lower it to get more "results" but expect them to be noiser.

## OUTPUT
The main output is a _Summary file that reports detailed information for every read.
Each read may appear multiple times, depending on the number of alignments it produces. For example, 
if a read contains both a vector rearrangement and an integration-site (IS) locus, it will appear
three times: twice for the alignments on chrV and once for the alignment on the host chromosome associated 
with the IS.

For each alignment, the file includes:

* the genomic start and end coordinates

* query_start and query_end, indicating where the aligned segment lies within the read

* (these fields allow reconstruction of the alignment structure)

* gene annotation

* additional variables used during analysis (like cigar string ecc)

Chimeric reads include, in the corresponding alignment rows, the inferred vector junction locus and the integration locus.

### Options
**-i**: The file .tsv with the paths to fastq files as last column  
 This file must have the same format of the sample_label.tsv. Specify in the last column the path to the already 
trimmed fastq file(s).  
**-t**: max threads to be used for file. NB since the pipeline is runned in parallel for the different fastq files 
specified with the -i option, each subprocess will use -t THREADS.  
**-v**: (optional, see -V) path to the fasta file with viral genome (e.g. AAV). **The viral genome must be a single 
sequence named chrV. This is essential for later steps.** If you have already made the index of the fasta file you 
can use the -V option otherwise the index will be created in the same directory of the fasta file. You can specify 
only one between -v or -V or both. -v is required only if you do not specify -m or -M since both viral genome and 
reference genome are needed.  
**-V**: (optional) path to the viral bwa-index with basename. If specified the index of the viral genome will not 
be made. If you have already created the mixed genome (see -m and -M) you can specify only 
-V (e.g. -V /path/to/index.fa -M /path/to/mixed.index.fa). If the index is the same directory of the viral genome 
fasta file you can use -V 1 e.g. `-v /path/to/viral.fa -V 1`. More example in the **example** section.  
**-r**: (optional, see -R) the fasta file with the reference genome (e.g. hg19). If you have already made the index
of the fasta file you can use the -R option. Otherwise the index will be created in the same directory of the fasta 
file. You can specify only one between -r or -R or both. -r is required only if you do not specify -m or -M since 
both viral genome and reference genome are needed.  
**-R**: (optional) path to the reference bwa-index with basename. If specified the index of the reference genome will 
not be made (see -V for more info).  
**-m**: (optional) the fasta file with the mixed genome (viral genome must be appended at the end with the sequence 
name chrV).
You can create this file running the following commands: 
```
cp reference.fa mixed.fa
cat viral.fa >> mixed.fa
```
If the option is not specified the mixed genome will be created and indexed by the program.  
**-M**: (optional) path to bwa-index of the mixed genome. If specified the index of the mixed genome will not be made. If you specify -M you do not need to specify -m.  
**-a** :the gtf file with the custom annotation (containing both reference and viral annotations). It must be sorted  
You can sort it using the command 
```
sort -k1,1 -k4,4n -k5,5n
```
**-o**: path to the output directory  

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
