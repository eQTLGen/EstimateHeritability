FROM nfcore/base
LABEL authors="c.a.warmerdam@umcg.nl" \
      description="Docker image containing requirements for running LDSC"

COPY environment.yml /
COPY requirements.txt /
RUN conda update conda
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ldsc/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
RUN apt-get update && apt-get install -y gcc libcurl4 g++
RUN pip install -r requirements.txt