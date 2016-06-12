perl -pi -w -e 's/haloe/halo/g;' ms.tex
perl -pi -w -e 's/colour/color/g;' ms.tex
rm *.aux
rm *.log
rm *.synctex.gz
rm *.bbl
rm *.blg
rm *.out
pdflatex ms
bibtex ms
pdflatex ms
pdflatex ms
pdflatex ms
