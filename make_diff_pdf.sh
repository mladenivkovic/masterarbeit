#!/bin/bash

olddir=../paper_submitted
# olddir=../paper_updated

if [ ! -d $olddir ]; then
    echo "Didn't find $olddir directory, cant work like this"
    exit 1
fi

latexdiff --flatten --allow-spaces --exclude-safecmd="citet,citep" $olddir/main.tex main.tex > diff_main.tex

pdflatex -jobname paper_diff diff_main.tex
bibtex diff_main
pdflatex -jobname paper_diff diff_main.tex
pdflatex -jobname paper_diff diff_main.tex



