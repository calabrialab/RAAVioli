# RAAVioli
_Recombinant Adeno-Associated Viral IntegratiOn anaLysIs (RAAVIoli)_

Bioinformatics pipeline for the identification and characterization of AAV integration sites and viral rearrangements.

---

## Install

You can use RAAVioli through Docker or install it locally on your server.

### Docker Installation

To build the Docker image manually:
```
git clone https://github.com/calabrialab/RAAVioli.git
docker build -t raavioli_docker .
```
Alternatively, you can download the pre-built `raavioli_docker` image from Docker Hub (to be added).

All Docker-ready jobs are available in the `DockerJobs` directory, showing how to run RAAVioli using containerized workflows.

---

### Local Installation

A working conda or mamba installation is required.

Clone the repository:
```
git clone https://github.com/calabrialab/RAAVioli.git
```
Then install one or both pipelines:

#### RAAVioliShort
```
cd RAAVioli_short
chmod +x setup_short.sh
./setup_short.sh
```
#### RAAVioliLong
```
cd RAAVioli_long
chmod +x mamba_setup.sh
./mamba_setup.sh
```
### Examples
Simulated data, example configurations, and scripts used in the MTM publication are available in the repository: MTMRaavioli
https://github.com/calabrialab/MTMRAAVioli.git.

---

## Q & A

**Q: What is the difference between the Short and Long pipelines?**  
A: RAAVioli supports both short- and long-read sequencing technologies.  
The short-read pipeline is optimized for experiments using vector-specific primers, where reads begin directly on the vector genome. It is designed for high-throughput PCR-based datasets.  
The long-read pipeline has no such constraints and works on a wider range of experimental setups.

**Q: Can I use RAAVioli with WGS or TES experiments?**  
A: Yes, the long-read pipeline can process WGS or TES by treating paired-end FASTQs as single-end reads.  
The short-read pipeline is not compatible with WGS/TES due to its vector-anchored primer assumptions.

**Q: Can I use RAAVioli with double-AAV vector experiments?**  
A: ___Although not formally tested___, two approaches are possible:

1. Run RAAVioli twice, once per vector.
2. Merge the two vector genomes by concatenation (removing duplicated ITRs if needed).  
   This approach also captures hybrid rearrangements between vectors.

_Adding two vectors as separate chromosomes is not supported in the current version and would require pipeline modifications._

**Q: Does RAAVioli_long support aligners other than BWA?**  
A: Not at the moment. Although minimap2 can be faster, internal tests showed limited differences for this specific problem. Support for additional aligners may be added in future updates.

---

## What’s Next?

Development is ongoing:

- Improved memory-efficient short-read processing
- Additional visualization examples from published datasets
- Feature refinements and usability improvements

Check the repository periodically for updates.

### Problems and Resolution
For any question about RAAVioli please open an issue here.