FROM centos:8

MAINTAINER Sean Houghton <sbn.houghton@gmail.com>

RUN yum -y update && yum clean all

RUN yum install -y wget

RUN mkdir main && \
	cd main && \
	mkdir fastq && \
	mkdir output


CMD ['wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./miniconda.sh']

CMD ['bash ./miniconda.sh -b -p ./miniconda']

CMD ['PATH=$PATH:/main/miniconda/bin']

CMD ['conda config --add channels defaults']
CMD ['conda config --add channels bioconda']
CMD ['conda config --add channels conda-forge']

COPY environment.yml /
RUN conda env create -f ./environment.yml # && conda clean -ay

CMD ['source activate regenRNA-1.0']
