# dockerCluster
#
# VERSION	0.0.1

FROM ubuntu:14.04
MAINTAINER Louis Gioia <lhgioia@scripps.edu>

# make sure the package repository is up to date
RUN apt-get update -y

# install dependencies and vim editor
RUN apt-get install -y python-setuptools
RUN apt-get install -y autoconf
RUN apt-get install -y python-dev
RUN apt-get install -y vim
RUN apt-get install -y nano
RUN apt-get install -y wget

# create StarCluster ssh key directory
RUN mkdir ~/.ssh

# install StarCluster
RUN easy_install StarCluster
RUN mkdir ~/.starcluster
RUN mkdir ~/.starcluster/plugins

# download config file (used to operate StarCluster)
RUN wget https://bitbucket.org/sulab/omics_pipe/downloads/starcluster_config_custombuild -O ~/.starcluster/config

# download plugin python script (used to generate custom StarCluster nodes)
RUN wget https://bitbucket.org/sulab/omics_pipe/downloads/omicspipe_config.py -O ~/.starcluster/plugins/omicspipe_config.py