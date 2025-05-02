#!/bin/bash

sudo apt update && sudo apt install -y --no-install-recommends unzip libreadline-dev && \
sudo ln -s /lib/x86_64-linux-gnu/libreadline.so.8 /usr/lib/x86_64-linux-gnu/libreadline.so.7 && \
wget https://github.com/mothur/mothur/releases/download/v1.48.2/Mothur.linux_8.zip && \
unzip Mothur.linux_8.zip
