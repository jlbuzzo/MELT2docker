###############################################################################
#
# Dockerfile to make a java container for Melt.
#
###############################################################################



# Import base image.
FROM openjdk:8-jdk-slim

# Some metadata
MAINTAINER Jose Leonel L. Buzzo




WORKDIR /

# All installations
RUN apt-get update && \
	apt-get install -y \
		wget \
		gcc \
		g++ \
		make \
		parallel \
		autoconf \
		automake \
		libtool \
		libtbb-dev \
		libncurses5-dev \
		libbz2-dev \
		liblzma-dev \
		zlib1g-dev \
		pigz \
		bzip2 \
		gzip && \
	apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*




# Create inputs, outputs and assets directories.
RUN mkdir \
	/home/config \
	/home/inputs \
	/home/outputs \
	/home/reference \
	/home/annotation \
	/home/assets \
	/home/tmp

# Copy MELT into it's directory.
#RUN mkdir /home/MELT
#RUN wget URL-MELT -o /home
#RUN tar xzvf MELTv2.1.4.tar.gz -o /home/MELT
COPY ./MELTv2.1.4 /home/MELTv2.1.4

# Copy Bowtie2 into it's directory and compile it.
#RUN mkdir /home/bowtie2
#RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip -o /home
#RUN unzip bowtie2-2.3.4.1-source.zip -o /home/MELT
COPY ./bowtie2-2.3.4.1 /home/bowtie2-2.3.4.1
RUN make -f /home/bowtie2-2.3.4.1/Makefile -C /home/bowtie2-2.3.4.1
RUN echo "export PATH=/home/bowtie2-2.3.4.1:$PATH" > /home/.bashrc
ENV PATH=/home/bowtie2-2.3.4.1:$PATH

# Create the work directory and adentrate it.
WORKDIR /home

# Add the main script.
ADD ./Makefile .
ADD ./config.mk .

# An important switch!
RUN touch ponga_switch

CMD ["bash"]
