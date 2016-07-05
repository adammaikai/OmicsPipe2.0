#!/bin/bash

mkdir .local/easybuild/modules/all/snpeff
mkdir -p .local/easybuild/software/snpeff/3.6
wget -c -P .local/easybuild/software/snpeff/3.6 http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download
unzip .local/easybuild/software/snpeff/3.6/snpEff_latest_core.zip -d .local/easybuild/software/snpeff/3.6/
cp /omics_pipe/dist/modulefiles/snpeff .local/easybuild/modules/all/snpeff/3.6