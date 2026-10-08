FROM python:3.14-slim

ARG BLAST_VERSION=2.17.0
ARG R_VERSION=4.6.1

ENV DEBIAN_FRONTEND=noninteractive \
    PATH=/opt/blast/bin:/usr/local/bin:${PATH} \
    PYTHONPATH=/opt/fungani:/usr/local/lib/python3.14/site-packages:${PYTHONPATH}

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    file \
    gfortran \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libpcre2-dev \
    libreadline-dev \
    libssl-dev \
    libxml2-dev \
    make \
    wget \
    xz-utils \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/blast && \
    wget -q "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz" -O /tmp/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
    tar -xzf /tmp/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz -C /opt/blast --strip-components=1 && \
    rm -f /tmp/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz

RUN wget -q "https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz" -O /tmp/R-${R_VERSION}.tar.gz && \
    tar -xzf /tmp/R-${R_VERSION}.tar.gz -C /tmp && \
    cd /tmp/R-${R_VERSION} && \
    ./configure --prefix=/usr/local --enable-R-shlib --without-x && \
    make -j"$(nproc)" && \
    make install && \
    rm -rf /tmp/R-${R_VERSION} /tmp/R-${R_VERSION}.tar.gz

RUN git clone https://github.com/gopalpeddinti/fungani.git /opt/fungani

WORKDIR /opt/fungani

RUN python -m pip install --upgrade pip && \
    python -m pip install --no-cache-dir -r requirements.txt && \
    python -m pip install --no-cache-dir . && \
    Rscript -e "install.packages(c('ggplot2', 'patchwork'), repos='https://cloud.r-project.org')"

CMD ["fungani", "-h"]
