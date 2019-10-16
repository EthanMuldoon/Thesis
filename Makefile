all: thesis
	evince pdfs/thesis.pdf &

thesis:
	pdflatex thesis.tex
	mv thesis.pdf pdfs/thesis.pdf

%:: %/chapter.tex
	pdflatex -jobname=$(firstword $(subst /, ,$<)) "\includeonly{$(firstword $(subst ., ,$<))}\input{thesis}"
	mv $(firstword $(subst /, ,$<)).pdf pdfs/$(firstword $(subst /, ,$<)).pdf

clean:
	rm pdfs/*
	rm *.log
	rm *.aux
	rm *.toc
	rm *.lof
	rm c[0-9]*/*.aux
