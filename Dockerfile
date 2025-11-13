FROM fedora:42
RUN dnf update
RUN dnf install -y samtools pigz R-core-devel openssl-devel libcurl-devel libxml2-devel
RUN dnf clean all
RUN Rscript -e 'utils::chooseCRANmirror(F,1); \
  install.packages(c("optparse", "readr", "BiocManager", "tictoc", "processx","httr2", "jsonlite", "curl"));\
  BiocManager::install(c("GenomicRanges", "plyranges"))' 

RUN mkdir -p /work/in && mkdir /work/out
WORKDIR /work

RUN curl -L https://github.com/nanoporetech/modkit/releases/download/v0.5.1-rc1/modkit_v0.5.1rc1_u16_x86_64.tar.gz | tar xvzf -
ENV PATH="$PATH:/work/dist_modkit_v0.5.1_8fa79e3"

RUN curl -O --insecure https://mnp-flex.org/MNP-flex-methylationepic-v-1-0-b5-manifest-file-MGMT-complete.sorted.bed
COPY make_mnpflex.R .
COPY entrypoint.R .

ENTRYPOINT ["/usr/bin/Rscript", "/work/entrypoint.R"]

