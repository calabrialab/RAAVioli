# RAAVioli  
_Recombinant Adeno-Associated Viral IntegratiOn anaLysIs (RAAVIoli)_

RAAVioli is a bioinformatics pipeline for the identification and characterization of AAV integration sites and viral genome rearrangements.

---

# RAAVioli Short-Read Pipeline  

## Installation
To install the RAAVioli Short-Read pipeline, run:
```
chmod +x setup_short.sh
./setup_short.sh
```

---

## Output

RAAVioli generates two primary output matrices: **SeqCount** and **ShsCount**.  
Each row corresponds to a unique integration site (IS), while each column represents a sample.  
Values represent either:
- **Sequence count:** number of distinct reads supporting that IS  
- **Shear site estimate:** inferred number of unique shear sites for that IS  

For each sample, RAAVioli also produces a `results_only_aav` file containing the characterization of reads fully derived from the AAV vector.  

> **Note:**  
> RAAVioli detects chimeric reads only in **R1** (the read containing the primer).  
> Reads without chimeras in R1 are classified as *AAV-only*.  
> R2 is used to validate chimeras found in R1.  
> Although it is technically possible to search for integration sites in R2, this would likely produce a high number of false positives due to missing contextual information.

In addition, a compressed file `results.tsv.gz` is generated concatenating each sample result.  
It provides per-read details such as integration locus, rearrangements, and alignment features — useful for plotting, debugging, or downstream analysis.

---

## Configuration and Variable Documentation

The RAAVioli Short-Read pipeline uses three configuration files:

1. **`mandatory_vars.txt`** – core workflow and directory setup  
2. **`step_vars.txt`** – alignment and fusion detection parameters  
3. **`isr_vars.txt`** – integration site reconstruction and annotation parameters  

---

### `mandatory_vars.txt`

**NGSWORKINGPATH:**  
Root directory where all analysis output files will be stored.

**DISEASE:**  
Experimental condition or study group name.

**PATIENT:**  
Patient or donor identifier. Defines the sample-level subdirectory.

**POOL:**  
Sequencing pool identifier.  
Output files will be located under `${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/${POOL}`.

**TMPDIR:**  
Directory for temporary files.  
RAAVioli writes temp files under `${TMPDIR}/${POOL}`.

**REMOVE_TMP_DIR:**  
Controls removal of temporary files after execution.  
Accepts `remove_tmp_yes` or `remove_tmp_no`.  
Temporary files can be large; they should normally be deleted unless you plan to rerun IS identification with modified parameters (see Note 1 below).

**R1_FASTQ:**  
Path to the folder containing R1 demultiplexed FASTQ (gzipped) files.  
R1 must contain the primer sequence.

**R2_FASTQ:**  
Path to the folder containing R2 demultiplexed FASTQ (gzipped) files.

**ASSOCIATIONFILE:**  
Path to a TSV file mapping identifiers to sample IDs
(see Note 2 below).

**MAXTHREADS:**  
Maximum number of threads to use for parallel operations such as alignment and sorting.

**OUTPUT_NAME:**  
Custom tag or version label appended to all output files (e.g., `v1`, `test1`).

**GENOME:**  
Path to the mixed genome FASTA used for alignment (see Note 3 below).

**VECTORGENOME:**  
Path to the vector or viral genome FASTA used for identifying vector-derived reads.
The index must be in the same directory. If you do not have it you can create it by running:
```
conda activate RAAVioliShort_env
bwa index /path/to/vectorgenome.fa
```

**alignment_vars_file:**  
Path to the file defining alignment-related variables (`alignment_step_vars.txt`).

**isr_vars_file:**  
Path to the file defining integration site reconstruction parameters (`isr_step_vars.txt`).

---

### `alignment_step_vars.txt`

**FUSIORERRORRATE:**  
Maximum allowed error rate when detecting primers.

**FUSION_PRIMERS:**  
FASTA file containing the primer sequences used for library preparation.

**BWA_MIN_ALN_LEN:**  
Minimum alignment length (in base pairs) required for a valid `bwa mem` alignment.

**minmapQ:**  
Minimum mapping quality threshold for host alignments. Reads below this threshold are discarded.

**mapQvec:**  
Minimum mapping quality threshold for vector alignments.  
Reads without at least one vector alignment are discarded.

---

### `isr_step_vars.txt`

**SUBOPTH:**  
Filters out integration loci falling within repetitive genomic regions, as described in VISPA2.  
Default value is `40`, but it can be lowered if needed.

**MINAAVMATCHES:**  
Minimum number of matched bases on the AAV (vector) side required for a valid integration read.

**MAXCLUSTERD:**  
Maximum distance (in base pairs) for merging nearby integration sites into a single cluster.

**MERGECOL:**  
Column used to merge integration sites when estimating abundance.  
Must be `"CompleteAmplificationID"` or `"AddedField1"` if merging IS from replicates.

**ANNOTATIONGTF:**  
Path to the genome annotation (GTF) file.  
Used to assign integration sites to genes or genomic features.

**ITR_DF:**  
Path to a TSV file containing Inverted Terminal Repeat (ITR) coordinates or sequences.  
Leave empty (`""`) if unused (see Note 4).

**system_sequences_df:**  
Optional TSV file listing system-specific sequences (e.g., nested PCR primers).  
Used to filter PCR artifacts.  
Leave empty (`""`) if unused (see Note 5).

---

## Additional Notes

**1. REMOVE_TMP_DIR**  
Temporary files can be very large.  
Normally, set `remove_tmp_yes` to automatically delete them.  
However, if you want to rerun the IS identification step with different thresholds, set `remove_tmp_no` and rerun RAAVioli **without** specifying `alignment_vars_file`.  
This reuses existing mappings and speeds up the process.

---

**2. ASSOCIATIONFILE**  
A tab-separated file with the following columns:

- **TagID:** the FASTQ tag prefix (RAAVioli expects `${R1_FASTQ}/${tag}.r1.fastq.gz` and `${R2_FASTQ}/${tag}.r2.fastq.gz`).  
- **CompleteAmplificationID:** unique sample name; can include metadata.  
- **AddedField1:** replicate ID. Integration sites from samples with different `CompleteAmplificationID` but the same `AddedField1` can be merged if the latter is specified in the isr_vars file. If no replicates are present, copy `CompleteAmplificationID`.  
- **concatenatePoolIDSeqRun:** pool ID (must match the `POOL` variable).

Example association files are available in  
`RAAVioli/Data/Short/files/`.

---

**3. GENOME**  
The mixed genome combines the host genome with the vector genome, which must have a single chromosome named `chrV`.  
To create it, run:
```
RAAVioli_Short/utils_scripts/create_mixed_genome.sh     -g <host_genome.fa>     -v <vector_genome.fa>     -o <output_dir>     -n <output_name.fa>
```
This script also generates the BWA index.  
The process can take some time but only needs to be done once per host/vector combination.

---

**4. ITR_DF**  
File listing ITR coordinates or reference sequences. Used to detect rearrangements involving ITR regions.

**5. system_sequences_df**  
File listing system-specific sequences (e.g., primer or linker artifacts) to be excluded during analysis.
