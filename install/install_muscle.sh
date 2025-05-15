#!/bin/bash

set -e

sudo apt update
sudo apt install -y --no-install-recommends time screen
#sudo apt install -y --no-install-recommends muscle

wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3
sudo mv muscle-linux-x86.v5.3 /usr/bin/muscle