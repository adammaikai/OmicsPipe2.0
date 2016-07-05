#!/usr/bin/perl
package SNPiR::config;

use lib qw(../);

use base 'Exporter';


################################################################################
################################################################################
# $Revision: $
# Authors: Robert Piskol ( piskol@stanford.edu )
# Last modification $Author: piskol $
#
# configuration file

our $BLATEXE = '/opt/applications/ucsc_tools/273/gnu/bin/blat'; #this is the path to the blat executable
our $SAMTOOLSEXE = '/opt/applications/samtools/0.1.19/gnu/bin/samtools'; #this is the path to the samtools executable
our $FASTAFROMBED = '/opt/applications/bedtools/2.17.0/gnu/bin/fastaFromBed'; #this is the path to the fastaFromBed executable in the bedtools package 

our @EXPORT = (
        '$BLATEXE',
        '$SAMTOOLSEXE',
        '$FASTAFROMBED',
);

1;