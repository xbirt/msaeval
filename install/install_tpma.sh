#!/bin/bash

set -e

sudo apt update
sudo apt install -y --no-install-recommends git build-essential gcc g++ make time screen
git clone https://github.com/malabz/TPMA.git
cd TPMA
make
cd ..
