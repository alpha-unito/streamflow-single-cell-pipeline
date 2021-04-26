FROM centos:7 AS builder
ARG CELLRANGER_VERSION="3.1.0"

COPY cellranger-${CELLRANGER_VERSION}.tar.gz /tmp/

RUN cd /tmp \
    && tar -xzf cellranger-${CELLRANGER_VERSION}.tar.gz \
    && mv cellranger-${CELLRANGER_VERSION} /opt/cellranger


FROM centos:7
LABEL maintainer="iacopo.colonnelli@unito.it"

ENV PATH=/opt/cellranger:${PATH}

COPY --from=builder /opt/cellranger /opt/cellranger
COPY bcl2fastq2-v2-20-0-linux-x86-64.zip /tmp/

RUN cd /tmp \
    && yum install -y \
        unzip \
        which \
    && unzip bcl2fastq2-v2-20-0-linux-x86-64.zip \
    && yum install -y bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm \
    && yum clean all \
    && rm -f bcl2fastq2-*
