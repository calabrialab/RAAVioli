# RAAVioli
_Recombinant Adeno-Associated Viral IntegratiOn anaLysIs (RAAVIoli)_

Bioinformatics pipeline for the identification and characterization of AAV integration sites and viral rearrangements.

# RAAVioli Short-Read Pipeline  
### Install
To install RAAVioli Short Pipeline run:
```
chmod +x setup_short.sh
./setup_short.sh
```
### Configuration and Variable Documentation

This document describes all configuration variables used in the **RAAVioli Short-Read** pipeline.  
They are organized into three configuration files:

1. `mandatory_vars.txt` – core workflow and directory setup  
2. `step_vars.txt` – alignment and fusion detection parameters  
3. `isr_vars.txt` – integration site reconstruction and annotation parameters  

---

## `mandatory_vars.txt`
**NGSWORKINGPATH:**  
Root working directory where all analysis outputs files will be stored.

**DISEASE:**  
Name of the experimental condition or study group. 

**PATIENT:**  
Identifier for the patient or donor. Defines the sample-level subdirectory.

**POOL:**  
Sequencing pool identifier. 

Output files will be under `${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/${POOL}`

**TMPDIR:**  
Directory for temporary files. Should have sufficient disk space.

**REMOVE_TMP_DIR:**  
Flag controlling removal of temporary files after completion.  
Accepts `remove_tmp_yes` or `remove_tmp_no`.
Since, tmp files can be very large, you should always delete this dir., unless you want to try different 
IS identification parameters. See below (1) for more information.


**R1_FASTQ:**  
Path to the folder containing R1 demultiplexed FASTQ Gzipped files. R1 must contain the primer.

**R2_FASTQ:**  
Path to the folder containing R2 demultiplexed FASTQ Gzipped files.

**ASSOCIATIONFILE:**  
Path to a TSV file mapping barcodes or identifiers to sample IDs.  
Required for demultiplexing and sample tracking. See below (2) for more information.

**MAXTHREADS:**  
Maximum number of threads to use for parallel processing (alignment, sorting, etc.).

**OUTPUT_NAME:**  
Custom tag or version label appended to all output files (e.g. `v1`, `test1`).

**alignment_vars_file:**  
Path to the file containing alignment-related variables (`alignment_vars.txt`).

**isr_vars_file:**  
Path to the file defining integration site reconstruction parameters (`isr_vars.txt`).

---

### `alignment_vars.txt`

**FUSIORERRORRATE:**  
Maximum allowed error rate for detecting primer.

**FUSION_PRIMERS:**  
FASTA file containing primer sequences used during sequencing.

**GENOME:**  
Path to the mixed genome FASTA used for alignment. See below (3) on how to create mixed genome.

**VECTORGENOME:**  
Path to the vector or viral genome FASTA used for identifying vector-derived reads.

**BWA_MIN_ALN_LEN:**  
Minimum alignment length (in base pairs) required for a valid `bwa mem` alignment.

**minmapQ:**  
Minimum mapping quality threshold for host alignments. Reads below this are discarded.

**mapQvec:**  
Minimum mapping quality threshold for vector alignments. 
Read do not containing at least a vector alignment are discarded.

---

### `isr_vars.txt`

**SUBOPTH:**  
Removes integration locus falling in repetitve regions of the genome as explained in Vispa2. 
Although default values is equal to 40 you can lower this threshold.

**MINAAVMATCHES:**  
Minimum number of matched bases on the AAV (vector) side required for a valid integration read.

**MAXCLUSTERD:**  
Maximum intracluster distance (in base pairs) when merging nearby integration sites.

**MERGECOL:**  
Column name used for merging integration site to estimate abundance.  
Must be one of `"CompleteAmplificationID"` or `"AddedField1"` when you want to merge IS from different replicates.

**ANNOTATIONGTF:**  
Path to the GTF annotation file of the target genome.  
Used to assign integration sites to genomic features and genes.

**ITR_DF:**  
Path to a TSV file containing Inverted Terminal Repeat (ITR) coordinates or sequences.  
Leave empty (`""`) if unused. See below (4).

**system_sequences_df:**  
Optional TSV listing system-specific sequences in nested PCR. Used to remove PCR artifacts.  
Leave empty (`""`) if unused. See below (5).

---

### Additional Notes
**1 REMOVE_TMP_DIR**

Since tmp files can be very large, you should always delete them by using remove_tmp_yes, however, 
if you want to try different variables/threshold in the IS identification step, set remove_tmp_no, and then re-run 
RAAVioli without specifying the alignment_vars_file variable. This will avoid re-mapping the reads, speeding up the process.
---
**2 ASSOCIATIONFILE**

A tab-separated file that must contain the following columns:

`TagID`: 

The tags of the demultiplexd fastq. RAAVioli will look for `${R1_FASTQ}/${tag}.r1.fastq.gz` and 
`${R1_FASTQ}/${tag}.r2.fastq.gz` files.

`CompleteAmplificationID`:

The sample name, must be unique. Here you can insert all information about the sample.

`AddedField1`:

Insert here replicates ID. You can specify this field to merge IS in the latest step. In this way, IS coming from different samples (different CompleteAmplificationID) but having same AddedField1 will be joined together.
If you don't have replicates, paste the CompleteAmplificationID

`concatenatePoolIDSeqRun`:

The pool id. Must be equal to the POOL variable.

***examples of association_file cam be found in `RAAVioli/Data/Short/Runs` directory.***

---

**3 GENOME**

This is the path to the mixed genome, made of a host genome and the vector genome with a single chromosome named chrV, to create a mixed genome run:
`RAAVioli_Short/utils_scripts/create_mixed_genome.sh -g <host_genome.fa> -v <vector_genome.fa> -o <output_dir> -n <output_name.fa>`
This will create also the bwa index, and can take some time. Once created, you can reuse it for all the experiments using the same host/vector genomes.

**4 ITR_DF**

**5 system_sequences_df**

