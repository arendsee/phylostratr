TARGET=ml

${TARGET}.pdf: ${TARGET}.Rnw
	Rscript -e "library(knitr); knit('${TARGET}.Rnw')"
	latexmk -pdf ${TARGET}.tex

.PHONY: clean
clean:
	latexmk -f -c ${TARGET}.tex
	rm -f ${TARGET}.tex *pdf

.PHONY: deepclean
deepclean:
	latexmk -f -c ${TARGET}.tex
	rm -rf cache figure ${TARGET}.tex ${TARGET}.pdf
