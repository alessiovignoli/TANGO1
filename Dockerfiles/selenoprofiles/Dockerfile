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


#    && echo '[NCBI]' > ~/.ncbirc \
#    && echo '' >> ~/.ncbirc \
#    && echo 'data=/usr/local/blast_tmp/blast-2.2.26/data' >> ~/.ncbirc \
#    && rm blast-2.2.26-x64-linux.tar.gz


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
RUN chmod +x /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/from ._version import __version__/from _version import __version__/g"  /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/from .load_config import selenoprofiles_config_content/from load_config import selenoprofiles_config_content/g" /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/from .MMlib3 import \*/from MMlib3 import \*/g" /usr/local/lib/python3.9/site-packages/selenoprofiles4/*.py \
	&& sed -i "s/from .selenoprofiles/from selenoprofiles/g" /usr/local/lib/python3.9/site-packages/selenoprofiles4/*.py \
	&& sed -i 's/home_config_filename = os.path.expanduser(\"\~\") + \"\//home_config_filename = \"/' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i 's/config_filename = os.path.expanduser(\"\~\") + \"\/.selenoprofiles_config.txt\"/config_filename = home_config_filename/' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& sed -i "s/selenoprofiles_data_dir = ~\/selenoprofiles_data\//selenoprofiles_data_dir = selenoprofiles_data/" /usr/local/lib/python3.9/site-packages/selenoprofiles4/load_config.py \
	&& sed -i 's/if not opt\["y"\]/if not sys.argv[2]=="-y"/' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
	&& ln -s /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py /usr/bin/selenoprofiles 


#&& sed -i '1281i BUBBAprint("inside if -setup config file location :", config_filename)' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
#&& sed -i '1317i BUBBAprint("after the resolve_link stuff" )' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
#&& sed -i "s/BUBBAprint/        print/" /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py 

#&& sed -i 's/home_config_filename = os.path.expanduser(\"\~\") + \"\//home_config_filename = \"\/usr\/local\/selenoprofiles_workdir\//' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
# Here it is doing so that it also work for singularity where the ~ could create problems

#&& selenoprofiles -setup \
#&& sed -i 's/if not opt\["y"\]/if not sys.argv[2]=="-y"/' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
#&& sed -i "s/selenoprofiles_data_dir = ~\/selenoprofiles_data/selenoprofiles_data_dir = selenoprofiles_data/" ~/.selenoprofiles_config.txt \
#&& cd ~ \
#&& selenoprofiles -download -y \

#&& sed -i "s/selenoprofiles_data_dir = ~\/selenoprofiles_data\//selenoprofiles_data_dir = ~\/selenoprofiles_data/" /usr/local/lib/python3.9/site-packages/selenoprofiles4/load_config.py \
#&& sed -i '25i print("home config dir :", home_config_filename)' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
#&& sed -i '26i print("install dir ;", selenoprofiles_install_dir)' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py \
#&& sed -i '27i print("dependencies :", os.path.abspath(obonet.__file__),  os.path.abspath(easyterm.__file__),  os.path.abspath(gzip.__file__))' /usr/local/lib/python3.9/site-packages/selenoprofiles4/selenoprofiles4.py 



# this is the last piece to fix relative import in python
ENV PYTHONPATH "${PYTHONPATH}:/usr/local/lib/python3.9/site-packages/selenoprofiles4"

# Adding also to the path the directory where blast is so that is callable from command line
# and it knows where the data/ dir is for the BLOSUM matrices and so on
ENV PATH "${PATH}:/usr/local/blast_tmp/blast-2.2.26/bin"

# Very important line for installing genewise wise2 package from biocontainer
COPY --from=0 /usr/local/bin/gene* /usr/bin


# Creating script to set some things that will only be done at lounch of the container
# this is to solve difference in dir architectures and sytems
#RUN mkdir -p /usr/local/start_dir/ \
#    && echo '#!/bin/bash' > /usr/local/start_dir/start.sh\
#    && echo 'selenoprofiles -setup' >> /usr/local/start_dir/start.sh \
#    && echo 'selenoprofiles -download -y' >> /usr/local/start_dir/start.sh \
#    && echo "echo '[NCBI]' > ~/.ncbirc" >> /usr/local/start_dir/start.sh \
#    && echo "echo '' >> ~/.ncbirc" >> /usr/local/start_dir/start.sh \
#    && echo "echo 'data=/usr/local/blast_tmp/blast-2.2.26/data' >> ~/.ncbirc" >> /usr/local/start_dir/start.sh \
#    && chmod +x /usr/local/start_dir/start.sh


#ENTRYPOINT ["/usr/local/start_dir/start.sh"]


# Here it is doing so that it also work for singularity where the ~ could create problems
# for this an artificial directory is created and set as entrypoint at executio
#RUN mkdir -p /usr/local/selenoprofiles_workdir \
#	&& cd /usr/local/selenoprofiles_workdir \
#	&& selenoprofiles -setup \
#	&& selenoprofiles -download -y
 
#WORKDIR /usr/local/selenoprofiles_workdir
