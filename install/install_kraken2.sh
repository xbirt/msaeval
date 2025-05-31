#!/bin/bash

set -e

sudo apt update
sudo apt install -y --no-install-recommends build-essential g++ zlib1g-dev make time screen rsync

wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.5.tar.gz
tar xf v2.1.5.tar.gz
rm -f v2.1.5.tar.gz
mv kraken2-2.1.5 kraken2
cd kraken2

sudo ./install_kraken2.sh ./bin

echo 'export PATH="$PATH:~/kraken2/bin/"' >> ~/.bashrc
source ~/.bashrc