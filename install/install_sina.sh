#!/bin/bash

set -e

sudo apt update
sudo apt install -y --no-install-recommends time screen

wget https://github.com/epruesse/SINA/releases/download/v1.7.2/sina-1.7.2-linux.tar.gz
tar xf sina-1.7.2-linux.tar.gz
rm -f sina-1.7.2-linux.tar.gz
mv sina-1.7.2-linux sina
