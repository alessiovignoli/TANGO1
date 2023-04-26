#wise2 container wise2 version 2.4.1 
#from bioconda
FROM quay.io/biocontainers/wise2:2.4.1--h40d77a6_1
LABEL maintainer="Alessio Vignoli" \
        name="alessiovignoli3/tango-project:selenoprofiles" \
        description="Docker image with all dependencies for selenoprofiles program"


#python:3.9.16-slim
FROM python@sha256:389f4d4e7180c8e41d89ff48e0d0522bb314c956787454b969aaef7d36d7d885

# Install apt get related packages
RUN apt-get update \
    && apt-get install wget \
    && apt-get install gawk \
    && apt-get install gettext -y \
    && apt-get install sqlite3 \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/*

# Install exonerate v2.4.0
# form gitrepo of nathanweeks
RUN git clone https://github.com/nathanweeks/exonerate.git \
    && cd exonerate \
    && git checkout v2.4.0 \
    && ./configure \
    && make AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: \
    && make check AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=: \
    && make install AUTOCONF=: AUTOHEADER=: AUTOMAKE=: ACLOCAL=:

# Install mafft
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.511-1_amd64.deb \
    && dpkg -i mafft_7.511-1_amd64.deb

# Install blast 
RUN mkdir -p /usr/local/blast_tmp/ \
    && cd /usr/local/blast_tmp/ \
    && wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz  \
    && tar -xf blast-2.2.26-x64-linux.tar.gz \
    && rm blast-2.2.26-x64-linux.tar.gz

# Install python related packages
# and selenoprofiles
RUN pip install \
    biopython==1.78 \
    brotlipy==0.7.0 \
    certifi==2022.6.15 \
    charset_normalizer==2.0.4 \
    cryptography==37.0.1 \
    easyterm==1.0.0 \
    idna==3.3 \
    obonet==0.2.3 \
    pyopenssl==22.0.0 \
    PySocks==1.7.1 \
    requests==2.27.1 \
    urllib3==1.26.9 \
    selenoprofiles4==4.4.5 \
    setuptools==61.2.0 \
    six==1.16.0

# Make selenoprofiles command line executable and changing all lines neccessary to fix:
# relative python imports of selenoprofiles 
# opt.resolve_links problem of easyterm package 
# -download -y option  
# translatable to a singularity image.
# however selenoprofiles -setup and selenoprofiles -download  has to be run by the user on execution
RUN sed -i 's/config_filename = os.path.expanduser(\"\~\") + \"\/.selenoprofiles_config.txt\"/config_filename = os.path.abspath(".\/.selenoprofiles_config.txt")/g' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/selenoprofiles_data_dir = ~\/selenoprofiles_data\//selenoprofiles_data_dir = selenoprofiles_data/" /usr/local/lib/python3.9/site-packages/selenoprofiles4/load_config.py \
	&& sed -i 's/if not opt\["y"\]/if not sys.argv[2]=="-y"/' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py 

# Creation of the file in bin that imports the selenoprofiles4 main script
RUN echo '#!/usr/local/bin/python3.9' > /usr/bin/selenoprofiles \
	&& echo '# -*- coding: utf-8 -*-' >> /usr/bin/selenoprofiles \
	&& echo 'import re' >> /usr/bin/selenoprofiles \
	&& echo 'import sys' >> /usr/bin/selenoprofiles \
	&& echo 'from selenoprofiles4.selenoprofiles4 import main' >> /usr/bin/selenoprofiles \
	&& echo "if __name__ == '__main__':" >> /usr/bin/selenoprofiles \
	&& echo "BUBBAsys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])" >> /usr/bin/selenoprofiles \
	&& echo 'BUBBAsys.exit(main())' >> /usr/bin/selenoprofiles \
	&& sed -i "s/BUBBA/    /g" /usr/bin/selenoprofiles \
	&& chmod 775 /usr/bin/selenoprofiles
	
# Adding also to the path the directory where blast is so that is callable from command line
# and it knows where the data/ dir is for the BLOSUM matrices and so on
ENV PATH "${PATH}:/usr/local/blast_tmp/blast-2.2.26/bin"

# Very important line for installing genewise wise2 package from biocontainer
COPY --from=0 /usr/local/bin/gene* /usr/bin
