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

# Install blast 
#RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.26/ncbi-blast-2.2.26+-x64-linux.tar.gz \
#    && tar -xf ncbi-blast-2.2.26+-x64-linux.tar.gz \
#    && mv ncbi-blast-2.2.26+/bin/* /usr/bin/ \
#    && rm -r ncbi-blast-2.2.26+

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz  \
    && tar -xf blast-2.2.26-x64-linux.tar.gz \
    && mv blast-2.2.26/bin/* /usr/bin/ \
    && rm -r blast-2.2.26 \
    && rm blast-2.2.26-x64-linux.tar.gz



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

# Make selenoprofiles command line executable
# and fixing python imports of selenoprofiles
# fixing also -download option  
RUN chmod +x /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/from ._version import __version__/from _version import __version__/g"  /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/from .load_config import selenoprofiles_config_content/from load_config import selenoprofiles_config_content/g" /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/from .MMlib3 import \*/from MMlib3 import \*/g" /usr/local/lib/python3.9/site-packages/selenoprofiles4/*.py \
	&& sed -i "s/from .selenoprofiles/from selenoprofiles/g" /usr/local/lib/python3.9/site-packages/selenoprofiles4/*.py \
	&& sed -i '1314i BUBBAopt["selenoprofiles_install_dir"] = selenoprofiles_install_dir' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/BUBBAopt/        opt/" /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/selenoprofiles_data_dir = ~\/selenoprofiles_data\//selenoprofiles_data_dir = ~\/selenoprofiles_data/" /usr/local/lib/python3.9/site-packages/selenoprofiles4/load_config.py \
	&& ln -s /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py /usr/bin/selenoprofiles \
	&& selenoprofiles -setup \
	&& sed -i 's/if not opt\["y"\]/if not sys.argv[2]=="-y"/' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/selenoprofiles_data_dir = ~\/selenoprofiles_data/selenoprofiles_data_dir = selenoprofiles_data/" ~/.selenoprofiles_config.txt
#	&& selenoprofiles -download -y 


# & sed -i "s/selenoprofiles_data_dir = ~\/selenoprofiles_data/selenoprofiles_data_dir = selenoprofiles_data/" ~/.selenoprofiles_config.txt \

ENV PYTHONPATH "${PYTHONPATH}:/usr/local/lib/python3.9/site-packages/selenoprofiles4"


# Very important line for installing genewise wise2 package from biocontainer
COPY --from=0 /usr/local/bin/gene* /usr/bin


