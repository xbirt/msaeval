#!/bin/bash

set -e

sudo apt update && sudo apt install -y --no-install-recommends build-essential libssl-dev libbz2-dev libreadline-dev libsqlite3-dev libopenmpi-dev zlib1g-dev mafft clustalw git curl wget
curl -sSL https://pyenv.run | bash
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(pyenv init -)"' >> ~/.bashrc
export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init - bash)"
pyenv install 2.7.18
pyenv global 2.7.18
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz
tar xzf blast-2.2.26-x64-linux.tar.gz
rm -f blast-2.2.26-x64-linux.tar.gz
sudo cp blast-2.2.26/bin/* /usr/local/bin/
sudo cp -r blast-2.2.26/data /usr/local/blast-data
python2 -m pip install numpy==1.16.6 scipy==1.2.3 cogent==1.5.3 mpi4py
git clone https://github.com/qiime/pynast.git pynast
cd pynast
python2 setup.py install
ln -s "$HOME/.pyenv/versions/2.7.18/bin/pynast" "$HOME/.pyenv/bin/pynast"
echo 'export BLASTMAT="/usr/local/blast-data"' >> ~/.bashrc
export BLASTMAT="/usr/local/blast-data"
