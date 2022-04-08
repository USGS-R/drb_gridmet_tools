 FROM condaforge/mambaforge

 COPY environment.yml /tmp/environment.yml

 RUN mamba env create -f /tmp/environment.yml && \
     mamba clean --all --yes
