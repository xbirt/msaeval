#!/bin/bash

# Download the aligned SILVA databases from the Mothur website
# https://mothur.org/wiki/silva_reference_files/

wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_2.tgz
wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_2.tgz

tar xf silva.nr_v138_2.tgz
tar xf silva.seed_v138_2.tgz

rm -f silva.nr_v138_2.tgz silva.seed_v138_2.tgz