FROM nfcore/base:1.13.1
LABEL authors="Sarai Varona and Sara Monzon" \
      description="Docker image containing all software requirements for the nf-core/viralrecon pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-viralrecon-1.2.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-viralrecon-1.2.0dev > nf-core-viralrecon-1.2.0dev.yml
