#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 regge_eikonal_pp.py output
python3 full_comparison.py data/TOTEM13_dsdt.csv output
python3 extensions.py data/TOTEM13_dsdt.csv output
TFIT=4.1 python3 extensions.py data/TOTEM13_dsdt.csv output_t4
python3 goodwalker.py data/TOTEM13_dsdt.csv output
TFIT=4.1 python3 goodwalker.py data/TOTEM13_dsdt.csv output_t4
python3 profile_inversion.py output
python3 energy_blackness.py output
MODELS=M1,M2,M3,M4 NSEED=4 python3 kmr_fits.py output
python3 universal_curves.py output
cd figures
python3 make_figures.py
python3 make_new_figures.py ../output
