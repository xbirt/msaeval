#!/bin/bash

set -e

#sudo apt update && sudo apt install -y --no-install-recommends mafft

#wget https://mafft.cbrc.jp/alignment/software/mafft_7.526-1_amd64.deb
#sudo dpkg -i mafft_7.526-1_amd64.deb
#rm -f mafft_7.526-1_amd64.deb

sudo apt update
# golang este necesar pentru DASH, time pentru benchmark
sudo apt install -y --no-install-recommends build-essential g++ make golang time screen

wget https://mafft.cbrc.jp/alignment/software/mafft-7.525-with-extensions-src.tgz
tar xzf mafft-7.525-with-extensions-src.tgz
cd mafft-7.525-with-extensions/core
# Activăm DASH
sed -i '/^#DASH_CLIENT = dash_client/s/^#//' Makefile
# Instalăm în /usr, nu /usr/local
sed -i 's|^PREFIX = /usr/local$|PREFIX = /usr|' Makefile
make
sudo make install
cd ../extensions
# Nu în /usr, trebuie neapărat în /usr/local; nu merge altfel, deci nu modificăm PREFIX
make
sudo make install
echo 'export PATH="$PATH:/usr/libexec/mafft/"' >> ~/.bashrc
rm -rf mafft-7.525-with-extensions mafft-7.525-with-extensions-src.tgz