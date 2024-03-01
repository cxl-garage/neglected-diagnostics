#!/bin/bash

# Install essential linux packages
sudo apt update
sudo apt install build-essential
sudo apt install mafft

# Create and setup virtual environment for Seq2Dx
conda init bash

conda create -y -n negDia python=3.10 pip
conda activate negDia
pip install -e ".[dev,docs]"

# Start Seq2Dx
streamlit run src/app/Home.py --server.headless true
