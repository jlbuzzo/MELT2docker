###############################################################################
#
# Dockerfile to make a java container for MELT.
# Image: MELT2docker
#
###############################################################################



# Import base image: Debian:strech-slim with openjdk.
FROM galantelab/biotools

# Some metadata
MAINTAINER Jose Leonel L. Buzzo



# System wide installations: apt-get.
RUN apt-get -y update && apt-get -y upgrade
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



# Begin installations.
WORKDIR /build

#> melt
COPY ./deps/MELTv2.1.4 MELTv2.1.4
#-

#> bowtie2
COPY ./deps/bowtie2-2.3.4.1 bowtie2-2.3.4.1
RUN cd bowtie2-2.3.4.1 && make
RUN echo "export PATH=/build/bowtie2-2.3.4.1:$PATH" >> /home/lion/.bashrc
ENV PATH=/build/bowtie2-2.3.4.1:$PATH
#-

# Finish installations.
WORKDIR /home/lion



CMD ["bash"]
