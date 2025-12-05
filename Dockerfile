FROM debian:12

ENV DEBIAN_FRONTEND=noninteractive

# ------------------------------------------------------------
# 1. Basic tools
# ------------------------------------------------------------
RUN apt-get update && apt-get install -y \
    wget curl git bzip2 ca-certificates nano \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# 2. Install Miniconda
# ------------------------------------------------------------
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Use bash login-style for build steps
SHELL ["/bin/bash", "-lc"]

# ------------------------------------------------------------
# 3. Configure conda
# ------------------------------------------------------------
RUN source /opt/conda/etc/profile.d/conda.sh && \
    conda config --system --set subdir linux-64 && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r


# ------------------------------------------------------------
# 4. Install RAAVioli
# ------------------------------------------------------------
WORKDIR /opt
RUN git clone https://github.com/calabrialab/RAAVioli.git

WORKDIR /opt/RAAVioli/RAAVioli_short
RUN chmod +x setup_short.sh && \
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate base && \
    ./setup_short.sh

WORKDIR /opt/RAAVioli/RAAVioli_long
RUN chmod +x mamba_setup.sh && \
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate base && \
    ./mamba_setup.sh

# ------------------------------------------------------------
# 5. Ensure conda loads in all shells
# ------------------------------------------------------------
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> /etc/bash.bashrc
# Optional auto-activation of base environment:
# RUN echo "conda activate base" >> /etc/bash.bashrc

# ------------------------------------------------------------
# 6. Entry point logic
# ------------------------------------------------------------
SHELL ["/bin/bash", "-lc"]
# Make container start with a login shell
ENTRYPOINT ["/bin/bash", "-l"]
