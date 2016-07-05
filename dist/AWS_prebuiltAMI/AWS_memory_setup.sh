#!/bin/bash

# SGE Memory Limit Setup for AWS EC2
# Kenneth Ban 30 July 2014
# To run once on headnode prior to running jobs

# Set memory limits for each node (complex values field: h_vmem=<memory of node>
# - retrieve nodename and total memory per node as key:value string
# - loop through each key:value and set h_vmem for memory limit per node

list=`qhost | tail -n+4 | awk '{print $1,":", $5}'` 

while read -r line; do
     node=$(echo $line | cut -d':' -f1 | tr -d ' ')
     mem=$(echo $line | cut -d':' -f2 | tr -d ' ')
     qconf -rattr exechost complex_values h_vmem=$mem $node
done <<< "$list"

# Set h_vmem as a consumable resource
# - find the h_vmem line and replace consumable flag from NO to YES
# - write temp file and use it to replace SGE configuration

qconf -sc | sed '/h_vmem              h_vmem     MEMORY      <=    YES         NO         0        0/ c\h_vmem              h_vmem     MEMORY      <=    YES         YES        0        0' > sge.conf
qconf -Mc sge.conf && rm sge.conf