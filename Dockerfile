###############################################################################
#
# Dockerfile to make a java container for MELT.
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
		tar \
		unzip \
		bzip2 \
		gzip && \
	apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*




# Create inputs, outputs and assets directories.
RUN mkdir \
	/home/config \
	/home/inputs \
	/home/outputs \
	/home/assets \
	/home/reference \
	/home/annotation \
	/home/tmp \
	/home/extra

### Copy MELT into it's container directory.
# Copy from internet.
#RUN mkdir /home/MELT
#RUN wget URL-MELT -o /home
#RUN tar xzvf MELTv2.1.4.tar.gz -o /home/MELT

# Offline copy.
COPY ./MELTv2.1.4 /home/MELTv2.1.4


### Copy Bowtie2 into it's directory and compile it.
# Copy from internet.
#RUN mkdir /home/bowtie2
#RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip -o /home
#RUN unzip bowtie2-2.3.4.1-source.zip -o /home/MELT

# Offline copy.
COPY ./bowtie2-2.3.4.1 /home/bowtie2-2.3.4.1
RUN make -f /home/bowtie2-2.3.4.1/Makefile -C /home/bowtie2-2.3.4.1
RUN echo "export PATH=/home/bowtie2-2.3.4.1:$PATH" > /home/.bashrc
ENV PATH=/home/bowtie2-2.3.4.1:$PATH


# Install Samtools from internet.
ENV SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 \
	SAMTOOLS_PATH=/samtools-1.7
RUN wget $SAMTOOLS_URL && tar xjf samtools-1.7.tar.bz2 && (cd $SAMTOOLS_PATH && make)
ENV PATH=$SAMTOOLS_PATH:$PATH



# Create the work directory and adentrate it.
WORKDIR /home

# Add the main script.
ADD ./Makefile .

# An important switch!
RUN touch .ponga_switch

CMD ["bash"]
