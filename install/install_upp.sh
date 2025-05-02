#!/bin/bash

set -e

sudo apt update && sudo apt install -y --no-install-recommends git python3 python3-pip python3-venv hmmer openjdk-17-jdk mafft fasttree
mkdir -p ~/venv/upp
python3 -m venv ~/venv/upp
source ~/venv/upp/bin/activate
git clone https://github.com/smirarab/sate-tools-linux.git
export PASTA_TOOLS_DEVDIR="$HOME/sate-tools-linux"
echo "export PASTA_TOOLS_DEVDIR=\"$PASTA_TOOLS_DEVDIR\"" >> ~/.bashrc
git clone https://github.com/smirarab/PASTA.git
cd PASTA
pip install -e . --config-settings="--sate-tools-dir=$PASTA_TOOLS_DEVDIR"
cd ~
git clone https://github.com/smirarab/sepp.git
cd sepp
pip install -e .
python3 setup.py config
python3 setup.py install
python3 setup.py upp
