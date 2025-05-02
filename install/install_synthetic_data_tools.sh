#!/bin/bash

set -e

export PERL_MM_USE_DEFAULT=1
SEQKIT_VERSION="2.10.0"
SEQKIT_TAR="seqkit_linux_amd64.tar.gz"
SEQKIT_URL="https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/${SEQKIT_TAR}"

sudo apt update
sudo apt install -y --no-install-recommends git curl build-essential r-base r-cran-rcpp \
    r-cran-matrixstats python3-dev python3-pip python3-venv libboost-all-dev hmmer mafft \
    fasttree perl make cpanminus bioperl cmake default-jre libcurl4-openssl-dev libssl-dev \
    libxml2-dev libglpk-dev libgsl-dev libfftw3-dev libtiff-dev libjpeg-dev libpng-dev \
    libx11-dev libxt-dev libbz2-dev liblzma-dev libpcre3-dev libicu-dev libfontconfig1-dev \
    libfreetype6-dev libharfbuzz-dev libfribidi-dev libgmp-dev libmpfr-dev libblas-dev \
    liblapack-dev gfortran parallel porechop vsearch seqtk minimap2 racon bcftools
wget https://sourceforge.net/projects/biogrinder/files/biogrinder/Grinder-0.5.4/Grinder-0.5.4.tar.gz/download -O Grinder-0.5.4.tar.gz
tar xzf Grinder-0.5.4.tar.gz
rm -f Grinder-0.5.4.tar.gz
cd Grinder-0.5.4
cpanm --installdeps . 
sudo cpan -f Module::Install
perl Makefile.PL
make
sudo make install
cd ~
git clone https://github.com/caboche/CuReSim-LoRM

sudo R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
sudo R -e 'if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")'
sudo R -e 'devtools::install_github("biobakery/SparseDOSSA2")'
sudo R -e 'if (!require("SparseDOSSA2")) BiocManager::install("SparseDOSSA2")'
sudo R -e 'BiocManager::install(c("S4Vectors", "SummarizedExperiment"))'
sudo R -e 'install.packages(c("remotes", "Rcpp", "irlba", "igraph", "statmod", "dqrng"), repos="https://cloud.r-project.org/")'
sudo R -e 'BiocManager::install(c("scuttle", "BiocSingular", "scran", "bluster", "edgeR"))'
sudo R -e 'remotes::install_git("https://gitlab.com/sysbiobig/metasparsim.git", 
                              dependencies = TRUE, 
                              upgrade = "never",
                              build_vignettes = FALSE)'

wget -O "$SEQKIT_TAR" "$SEQKIT_URL"
tar -zxvf "$SEQKIT_TAR"
sudo mv seqkit /usr/local/bin/
sudo chmod +x /usr/local/bin/seqkit
rm -f "$SEQKIT_TAR"
