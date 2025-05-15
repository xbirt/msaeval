#!/bin/bash

set -e

sudo apt update
sudo apt install -y --no-install-recommends git python3 python3-pip python3-venv hmmer openjdk-17-jdk mafft fasttree time screen
mkdir -p ~/venv/pasta
python3 -m venv ~/venv/pasta
source ~/venv/pasta/bin/activate
git clone https://github.com/smirarab/sate-tools-linux.git
export PASTA_TOOLS_DEVDIR="$HOME/sate-tools-linux"
echo "export PASTA_TOOLS_DEVDIR=\"$PASTA_TOOLS_DEVDIR\"" >> ~/.bashrc
git clone https://github.com/smirarab/PASTA.git
cd PASTA
pip install -e . --config-settings="--sate-tools-dir=$PASTA_TOOLS_DEVDIR"
cd ~
