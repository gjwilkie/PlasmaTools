pdflatex operator.tex
if [ -n $? ]
then
   bibtex operator
   pdflatex operator.tex
fi
