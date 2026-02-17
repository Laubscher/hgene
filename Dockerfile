FROM bioconductor/bioconductor_docker:3.19 AS base

RUN R -q -e "install.packages(c('tidyverse','rmarkdown','knitr','Biostrings','GenomicAlignments','Rsamtools','kableExtra','officer'), repos='https://cloud.r-project.org')"
RUN R -q -e "BiocManager::install(c('Biostrings','Rsamtools','GenomicRanges','GenomeInfoDb','GenomicAlignments','kableExtra','officer'), ask=FALSE, update=FALSE)"

# Install packages dependencies
RUN apt-get update && apt-get install -y \
      libz-dev \
      liblzma-dev \
      libbz2-dev \
      libnetcdf-dev \
      bzip2 \
      pigz \
      curl \
      git \
      make \
    && rm -rf /var/cache/apt/* /var/lib/apt/lists/*;

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# samtools
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
FROM base AS samtools
RUN mkdir -p /app \
   && curl -kL https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 \
   | tar -C /app --strip-components=1 -jxf -
RUN cd /app && ./configure && make -j4

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# bcftools
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
FROM base AS bcftools
RUN mkdir -p /app \
   && curl -kL https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 \
   | tar -C /app --strip-components=1 -jxf -
RUN cd /app && ./configure && make -j4

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# htslib
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
FROM base AS htslib
RUN mkdir -p /app \
   && curl -kL https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 \
   | tar -C /app --strip-components=1 -jxf -
RUN cd /app && ./configure && make -j4 \
  && sed -e 's#@-includedir@#/usr/include#g;s#@-libdir@#/usr/lib#g;s#@-PACKAGE_VERSION@#1.19.1#g' htslib.pc.tmp > htslib.pc


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# minimap2
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
FROM base AS minimap2
RUN mkdir -p /app \
   && curl -kL https://github.com/lh3/minimap2/releases/download/v2.27/minimap2-2.27.tar.bz2 \
   | tar -C /app --strip-components=1 -jxf -
RUN cd /app \
  && (if [ "$TARGETPLATFORM" = "linux/arm64" ]; then \
      make arm_neon=1 aarch64=1; \
  elif [ "$TARGETPLATFORM" = "linux/amd64" ]; then \
      make; \
  fi)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Main app
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
FROM base AS main

RUN apt-get update && apt-get install -y \
      seqtk \
      fastp \
      bwa \
      kmc \
      jellyfish \
      libkmc-dev \
      bcalm \
      rna-star \
      hisat2 \
      stringtie \
      ncbi-blast+ \
      bowtie \
      bowtie2 \
      python3 \
      python3-pip \
      fastqc \
      cutadapt \
      subread \
      salmon \
      porechop \
&& rm -rf /var/cache/apt/* /var/lib/apt/lists/*; ln -s /usr/bin/python3 /usr/bin/python

RUN pip3 install --no-cache-dir pysam==0.22.0

# minimap + (k8/paftools.js)
#COPY --from=k8 /app/k8/k8 /usr/local/bin/
COPY --from=minimap2 /app/minimap2 /app/misc/paftools.js /usr/local/bin/

# samtools 
COPY --from=samtools \
  /app/samtools \
  /app/misc/ace2sam \
  /app/misc/maq2sam-long \
  /app/misc/maq2sam-short \
  /app/misc/md5fa \
  /app/misc/md5sum-lite \
  /app/misc/wgsim \
  /app/misc/plot-* \
  /app/misc/*.pl \
  /usr/bin/
  
# bcftools
COPY --from=bcftools \
  /app/bcftools \
  /app/misc/*.py \
  /app/misc/run-roh.pl \
  /app/misc/vcfutils.pl \
  /usr/bin/
COPY --from=bcftools /app/plugins/*.so /usr/local/libexec/bcftools/

# htslib
COPY --from=htslib /app/annot-tsv /app/bgzip /app/htsfile /app/tabix /usr/bin
COPY --from=htslib /app/htslib/*.h /usr/include/htslib/
COPY --from=htslib /app/libhts.so /usr/lib/libhts.so.1.19.1
COPY --from=htslib /app/libhts.a /usr/lib/
COPY --from=htslib /app/htslib.pc /usr/lib/pkgconfig/
RUN ln -sf libhts.so.1.19.1 /usr/lib/libhts.so \
 && ln -sf libhts.so.1.19.1 /usr/lib/libhts.so.3


COPY --chmod=755 bin/* /usr/local/bin/
COPY --chmod=755 src/* /usr/local/bin/
COPY --chmod=755 db/* /usr/local/db/
COPY --chmod=755 notebooks/* /usr/local/notebooks/
COPY data/db/ /usr/local/data/db/
COPY LICENSE.txt /usr/local/share/licenses/hgene/LICENSE

   
ENV PATH="/usr/local/bin:${PATH}"

# --- Release metadata (OCI labels) + default entrypoint (for Singularity/Apptainer "run") ---
ARG VERSION=dev
ARG VCS_REF=unknown
ARG BUILD_DATE=unknown

LABEL org.opencontainers.image.title="hsv-hgene" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}" \
      org.opencontainers.image.source="https://github.com/Laubscher/hsv-hgene.git" \
      org.opencontainers.image.licenses="MIT"

# Makes `singularity run <image>.sif ...` work (runs hgene by default)
ENTRYPOINT ["hgene"]
CMD ["--help"]
