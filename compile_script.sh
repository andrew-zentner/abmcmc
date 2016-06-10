perl -pi -w -e 's/haloe/halo/g;' ms.tex
perl -pi -w -e 's/colour/color/g;' ms.tex
git clean -dfx
pdflatex ms.tex
bibtex ms.bib
pdflatex ms.tex
pdflatex ms.tex
pdflatex ms.tex
