FROM nfcore/base
LABEL authors="urmovosa@ut.ee" \
      description="Docker image containing requirements for heritability estimation of eQTLGen results"

COPY environment.yml /
RUN conda update conda
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/eQTLGenHeritabilityEstimation/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
RUN apt-get update && apt-get install -y gcc libcurl4
COPY data/SumTool-master.zip /
RUN R -e "devtools::install_local('SumTool-master.zip', dependencies = TRUE)"
RUN R -e "remotes::install_version('arrow', version = '11.0.0.2', dependencies = TRUE, repos = 'http://cran.rstudio.com/', upgrade = 'never')"
RUN R -e "remotes::install_version('bigsnpr', version = '1.10.8', dependencies = TRUE, repos = 'http://cran.rstudio.com/', upgrade = 'never')"

