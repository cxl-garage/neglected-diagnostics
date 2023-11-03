#!/bin/bash

rm -rf dist
python -m venv neglected-diagnostics-venv
source neglected-diagnostics-venv/bin/activate
pip install -r requirements.txt
pyinstaller ndiag.spec --clean
deactivate
rm -rf neglected-diagnostics-venv
rm -rf build
cp -r pages dist/ndiag/
cp -r assets dist/ndiag/
cp Home.py dist/ndiag/
cp -r ../.streamlit dist/ndiag/