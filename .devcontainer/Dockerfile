# debian buster
FROM maxh/base_conda:0.1

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# update apt-get and install needed packages
RUN sudo apt-get update --fix-missing \
    && sudo apt-get install -y --no-install-recommends \
        build-essential \
        libcairo2-dev \
        # install git iproute2, process tools, lsb-release (common in install
        # instructions for CLIs) and libzip for R Tools extension
        iproute2 \
        procps \
        lsb-release \
        libzip-dev \
    # clean up
    && sudo apt-get autoremove -y \
    && sudo apt-get clean -y \
    && sudo rm -rf /var/lib/apt/lists/*

# install .Net runtime for R Tools extension
# this doesn't work when attached to above (something to do with wget)
RUN wget --quiet https://packages.microsoft.com/config/debian/10/packages-microsoft-prod.deb \
    && sudo dpkg -i packages-microsoft-prod.deb \
    && sudo apt-get install apt-transport-https \
    && sudo apt-get update \
    && sudo apt-get -y install dotnet-runtime-2.2 \
    # clean up again
    && sudo apt-get autoremove -y \
    && sudo apt-get clean -y \
    && sudo rm -rf /var/lib/apt/lists/*

# install R through conda
RUN conda update conda \
    && conda install -c conda-forge r

# set default R mirror
RUN echo '# Default repo \n\
local({ \n\
    r <- getOption("repos") \n\
    r["CRAN"] <- "https://cloud.r-project.org" \n\
    options(repos=r) \n\
})' >> /home/maxh/.Rprofile

# install R packages available on conda-forge
RUN conda install -c conda-forge \
        r-import \
        r-roxygen2 \
        r-scrime \
        r-tidyverse \
        r-igraph \
    && conda install -c bioconda \
        bioconductor-seqarray \
        bioconductor-snprelate

RUN sudo ln -s /bin/tar /bin/gtar \
    && Rscript -e 'install.packages(c("stringi"))'
