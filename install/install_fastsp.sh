#!/bin/bash

set -e

sudo apt update
sudo apt-get install -y --no-install-recommends default-jdk git time screen
git clone https://github.com/smirarab/FastSP.git