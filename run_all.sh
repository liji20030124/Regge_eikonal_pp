#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 regge_eikonal_pp.py output
python3 full_comparison.py data/TOTEM13_dsdt.csv output
python3 extensions.py data/TOTEM13_dsdt.csv output
TFIT=4.1 python3 extensions.py data/TOTEM13_dsdt.csv output_t4
python3 goodwalker.py data/TOTEM13_dsdt.csv output
TFIT=4.1 python3 goodwalker.py data/TOTEM13_dsdt.csv output_t4
cd figures
python3 make_figures.py
