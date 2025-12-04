# Docker Jobs
This directory contains the 3 Docker jobs you can run inside docker.

- [Create Mixed Genome](#create-mixed-genome)
- [Run RAAVioli Short]((#running-raavioli-on-short-read-data))
- [Run RAAVioli Long](#running-raavioli-on-long-read-data)

_When you run RAAVioli trough docker, you need to specify in the config
paths docker paths and then mount them when running the job (see below)_

## Creating the Mixed Genome

RAAVioli requires a mixed genome that contains the host reference genome (e.g. hg19)
concatenated with the AAV vector genome.
This job provides a generic, Docker-compatible wrapper around the RAAVioli script:

`RAAVioli_short/utils_scripts/create_mixed_genome.sh`

to create and index the mixed genome.

---

### Input Requirements

You must provide:

1. Host genome FASTA (e.g. hg19.fa, mm10.fa, or custom reference)
2. Vector genome FASTA (typically vector.fa)
3. Output directory for writing the mixed genome
4. (Optional) Output file name (default: mixed.fa)

---

### Arguments

- host <path>:      Path to the host reference FASTA
- vector <path>:    Path to the vector FASTA
- outdir <path>:   Directory where the mixed genome will be created
- name <filename>:  Optional output FASTA name (default: mixed.fa)

---

### Docker Usage
```
docker run -d -name job_name \
    -v /path/to/host_genome:/opt/host \
    -v /path/to/vector_genome:/opt/vector \
    -v /path/to/output_dir:/opt/output \
    raavioli_docker \
    /opt/Raavioli/DockerJobs/setup_genome.sh \
        --host /opt/host/hg19.fa \
        --vector /opt/vector/vector.fa \
        --outdir /opt/output \
        --name mixed.fa
```
---

### Output

Output files will appear in the specified --outdir:

```
mixed.fa
mixed.fa.amb
mixed.fa.ann
mixed.fa.bwt
mixed.fa.pac
mixed.fa.sa
```

These files constitute the fully indexed mixed genome ready for alignment.

----

---
## Running RAAVioli on Short-Read Data

RAAVioliShort requires the configuration files (`mandatory_vars.txt`) that defines all
the paths needed by the pipeline. When running inside Docker, **every path inside
the config file must point to a path *inside the container***.

This means you must:

1. Mount your local directories into fixed locations inside Docker.
2. Write the config file using those *container paths*.

---

### 1. Create configuration files

Create a file such as:

`/path/to/YourProject/configs/mandatory_vars.txt`

Inside it, define all required variables as explained in the RAAVioli_short dir
**using container paths**:

### Example mandatory_vars.txt
```
VECTOR_GENOME=/opt/genomes/vector/vector.fa
MIXED_GENOME=/opt/genomes/mixed/mixed.fa
...
OUTPUT_DIR=/opt/output/sample1
isr_vars_file=/opt/configs/isr_step.txt
```
Important notes:

- /opt/raw, /opt/genomes, and /opt/output must be mounted when running Docker.
- Do **not** put host paths (like /home/user/...) inside the config.  
  Only use the container paths shown above.

---

### 2. Run Docker Jobs run_short.sh

Every directory referenced in your config must be mounted into Docker, for example:

- The directory with FASTQs → /opt/raw
- The directory with genomes → /opt/genomes
- The output directory → /opt/output
- The directory with your config files → /opt/config

Example (remember to assign a job name to show logs later on) :
```
docker run -d --name job_short \
    -v /path/to/raw_fastq:/opt/raw \
    -v /path/to/genomes:/opt/genomes \
    -v /path/to/output:/opt/output \
    -v /path/to/configs:/opt/config \
    raavioli_docker \
    /opt/RAAVioli/DockerJobs/run_short.sh \
    /opt/config/mandatory_vars.txt
```
---

---
## Running RAAVioli on Long-Read Data (RAAVioli_long)

RAAVioli_long requires several input files (FASTQs, genome indices, annotation, and
variables files) whose locations must be provided to `RAAVioli_long.sh`.

When running inside Docker, every path provided to RAAVioli_long must refer to a
path inside the container, not the host filesystem.

This means you must:

1. Mount your local directories to fixed locations inside Docker.
2. Write the TSV and variables files using container paths.

---

### 1. Prepare the input files

RAAVioli_long requires:

- A sample label file (sample_label.tsv)
  containing sample metadata and the FASTQ paths in the last column
- The mixed-genome BWA index
- The viral-genome BWA index
- The gene annotation GTF
- Three variables files:
  - variables_mixed
  - variables_viral
  - variables_rscript

All paths inside these files must be container paths, such as:
```
/opt/raw/sample1.fastq.gz
/opt/genomes/mixed/mixed.fa
/opt/genomes/vector/vector.fa
/opt/annot/custom.gtf
```
---

### 2. Example sample_label.tsv

The last column must contain FASTQ paths inside the container:
```
sample1    conditionA    /opt/raw/sample1.fastq.gz
sample2    conditionA    /opt/raw/sample2.fastq.gz
```
---
### 3. Running the Docker job
```
docker run -d --name job_long \
    -v /path/to/raw_fastq:/opt/raw \
    -v /path/to/genomes:/opt/genomes \
    -v /path/to/annotation:/opt/annot \
    -v /path/to/configs:/opt/config \
    -v /path/to/output:/opt/output \
    raavioli_docker \
    /opt/Raavioli/DockerJobs/run_long.sh \
        -i /opt/config/sample_label.tsv \
        -t 8 \
        -V /opt/genomes/vector/vector.fa \
        -M /opt/genomes/mixed/mixed.fa \
        -a /opt/annot/custom.gtf \
        -o /opt/output \
        -c /opt/config/variables_mixed \
        -w /opt/config/variables_viral \
        -y /opt/config/variables_rscript
```
Monitor progress using:
```
docker logs job_long
```

---

