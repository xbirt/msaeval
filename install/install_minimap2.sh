#!/bin/bash

set -e

sudo apt update
sudo apt install -y --no-install-recommends time screen

wget https://github.com/lh3/minimap2/releases/download/v2.29/minimap2-2.29_x64-linux.tar.bz2
tar xf minimap2-2.29_x64-linux.tar.bz2
rm -f minimap2-2.29_x64-linux.tar.bz2
mv minimap2-2.29_x64-linux minimap2
