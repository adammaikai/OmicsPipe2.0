#!/bin/bash

mkdir .local/easybuild/modules/all/cufflinks
mkdir -p .local/easybuild/software/cufflinks/2.2.1
wget -c -P .local/easybuild/software/cufflinks/2.2.1 http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar xvzf .local/easybuild/software/cufflinks/2.2.1/cufflinks-2.2.1.Linux_x86_64.tar.gz -C .local/easybuild/software/cufflinks/2.2.1/
cp /omics_pipe/dist/modulefiles/cufflinks .local/easybuild/modules/all/cufflinks/2.2.1