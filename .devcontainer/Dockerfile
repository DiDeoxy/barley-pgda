#-------------------------------------------------------------------------------------------------------------
# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT License. See https://go.microsoft.com/fwlink/?linkid=2090316 for license information.
#-------------------------------------------------------------------------------------------------------------

FROM rocker/r-apt:bionic

# The r-apt image comes with a base non-root 'docker' user which this Dockerfile
# gives sudo access. However, for Linux, this user's GID/UID must match your local
# user UID/GID to avoid permission issues with bind mounts. Update USER_UID / USER_GID 
# if yours is not 1000. See https://aka.ms/vscode-remote/containers/non-root-user.
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Configure apt and install packages
RUN apt-get update \
    && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends apt-utils dialog 2>&1 \
    #
    # install git iproute2, process tools, lsb-release (common in install instructions for CLIs) and libzip for R Tools extension
    && apt-get -y install git openssh-client less iproute2 procps lsb-release libzip-dev \
    #
    # Register Microsoft key and feed
    && wget -q https://packages.microsoft.com/config/ubuntu/18.04/packages-microsoft-prod.deb \
    && dpkg -i packages-microsoft-prod.deb \
    #
    # Install the .NET Runtime
    && add-apt-repository universe \
    && apt-get install apt-transport-https \
    && apt-get update \
    && apt-get -y install dotnet-runtime-3.1 \
    #
    # [Optional] Update a non-root user to match UID/GID - see https://aka.ms/vscode-remote/containers/non-root-user.
    && if [ "$USER_GID" != "1000" ]; then groupmod docker --gid $USER_GID; fi \
    && if [ "$USER_UID" != "1000" ]; then usermod --uid $USER_UID docker; fi \
    # [Optional] Add sudo support for the non-root user
    && apt-get install -y sudo \
    && echo docker ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/docker\
    && chmod 0440 /etc/sudoers.d/docker \
    #
    # Clean up
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# setup ssh
RUN echo 'eval $(ssh-agent -s)' >> /home/docker/.bashrc \
    && echo 'ssh-add ~/.ssh/id_rsa' >> /home/docker/.bashrc

# set default R mirror
RUN echo '# Default repo \n\
local({ \n\
    r <- getOption("repos") \n\
    r["CRAN"] <- "https://cloud.r-project.org" \n\
    options(repos=r) \n\
})' >> /home/docker/.Rprofile

RUN apt-get update \
    && apt-get install -y \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages( \
    c( \
        "BiocManager", "countrycode", "dbscan", "htmlwidgets", "igraph", \
        "import", "magrittr", "parallel", "plotly", "plyr", "readr", "stringr", \
        "tibble" \
    ) \
)'

RUN Rscript -e 'BiocManager::install( \
    c( \
        "SNPRelate", "SeqArray" \
    ) \
)'