# RAAVioli  
_Recombinant Adeno-Associated Viral IntegratiOn anaLysIs (RAAVIoli)_

A bioinformatics pipeline for the identification and characterization of AAV integration sites and viral genome rearrangements.

---

## Run Examples

This directory contains example FASTQ files to test the RAAVioli pipeline using an AAV Tomato vector and the **hg19** reference genome.

Before running any analysis, you must first create a **mixed genome** that includes both the host (hg19) and vector (AAV) sequences.

### Step 1 — Create the mixed genome
You need the `hg19.fa` FASTA file. Then, run the following command (assuming you are in the `Data` directory):

```
../RAAVioli_Short/utils_scripts/create_mixed_genome.sh \
    -g <path/to/hg19.fa> \
    -v <PATH_TO_RAAVIOLI_DIR>/Data/genome/vector.fa \
    -o <PATH_TO_RAAVIOLI_DIR>/Data/genome \
    -n mixed.fa
```

---

### Step 2 — Short-read examples

This repository provides two short-read example runs:
- One containing **chimeric** reads (AAV-host integrations)
- One containing **AAV-only** reads

Before running RAAVioli, open both `mandatory_vars` and `alignment_vars` and replace all instances of `${RAAVIOLI_DATA_DIR}` with the actual path to your `Data` directory.  

For example, if RAAVioli is installed in `/opt/applications/`, use:
```
/opt/applications/RAAVioli/Data/
```

Then, run one of the example pipelines from within the `RAAVioli_short` directory:

```
./RAAVioli_short.sh path/to/RAAVioli/Data/Short/RUNS/mandatory_vars_chimeras.txt
```

or

```
./RAAVioli_short.sh path/to/RAAVioli/Data/Short/RUNS/mandatory_vars_aavonly.txt
```

---

### Step 3 — Long-read examples

There are three example runs available for long-read sequencing data.

Before running any of them, edit the corresponding `sample_labels.tsv` file located in `Data/Long/RUNS/` to update the input data paths.

Once this is done and the mixed genome has been generated, you can launch the pipeline from the `RAAVioli_long` directory with:

```
nohup  RAAVioli_long.sh \
    -i <PATH_TO_RAAVIOLI_DIR>/Data/Long/RUNS/sample_labels_hg19AAV.tsv \
    -t 10 \
    -V <PATH_TO_RAAVIOLI_DIR>/Data/genome/vector.fa \
    -a <PATH_TO_RAAVIOLI_DIR>/Data/Long/RUNS/sampleannot.gtf \
    -o <PATH_TO_RAAVIOLI_DIR>/OUTPUT/RESULTS-LONG/ \
    -R <PATH_TO_RAAVIOLI_DIR>/Data/genome/mixed.fa \
    -M <PATH_TO_RAAVIOLI_DIR>/Data/genome/mixed.fa \
    -c <PATH_TO_RAAVIOLI_DIR>/Data/Long/RUNS/variables_mixed \
    -w <PATH_TO_RAAVIOLI_DIR>/Data/Long/RUNS/variables_viral \
    -y <PATH_TO_RAAVIOLI_DIR>/Data/Long/RUNS/variables_rscript \
    > log.out 2>&1 &
```
