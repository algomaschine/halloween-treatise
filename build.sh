#!/usr/bin/env bash
# Rebuild everything from the raw inputs. Requires: python3 (numpy, pandas, scipy, statsmodels, prophet,
# matplotlib, pillow), pandoc >= 3, a TeX Live with xelatex + latexmk, and the fonts Liberation Serif/Sans,
# DejaVu Sans Mono and DejaVu Math TeX Gyre.
set -euo pipefail
cd "$(dirname "$0")"

python3 code/parse_sw.py data/SW-All.txt                      # -> kp_ap_daily.csv
mv kp_ap_daily.csv data/
( cd figures/source && python3 ../../code/digitise.py )       # -> digitised_*.csv (reads the three source PNGs)
mv figures/source/digitised_*.csv data/
( cd data && python3 ../code/analysis.py )                    # -> results.json
( cd data && python3 ../code/figures.py && mv fig_S*.png ../figures/ )
( cd data && python3 ../code/listing1.py )                    # prints the headline numbers of Section 3

pandoc treatise.md -s -t latex -o treatise.tex \
  --pdf-engine=xelatex -V documentclass=article -V fontsize=10pt -V papersize=a4 -V geometry:margin=2.2cm \
  -V mainfont="Liberation Serif" -V sansfont="Liberation Sans" -V monofont="DejaVu Sans Mono" \
  -V monofontoptions="Scale=0.82" -V mathfont="DejaVu Math TeX Gyre" \
  -V linkcolor=black -V urlcolor=blue -V colorlinks=true --highlight-style=monochrome -H header.tex
sed -i '/^\\usepackage{lmodern}$/d' treatise.tex              # lmodern is not needed under xelatex/fontspec
latexmk -xelatex -interaction=nonstopmode -halt-on-error treatise.tex
pandoc treatise.md -o treatise.docx
