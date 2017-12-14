# For generating PDFs from LaTeX source
#LATEX = latexmk -pdfxe  # for XeLaTeX
#LATEX = latexmk -pdflua # for LuaLaTeX
LATEX = latexmk -pdf # for pdflatex; probably safest

# For converting LaTeX to ODT, which can then be cleaned up and saved as
# DOCX for those who need Word-formatted documents.
MK4HT = mk4ht oolatex

# For converting ODT to DOCX
SOFFICE = /Applications/LibreOffice.app/Contents/MacOS/soffice

all: pdf odt docx

pdf: Act01.pdf Act03.pdf

odt: Act01.odt Act03.odt

docx: Act01.docx Act03.docx

%.pdf: %.tex mdpchem.sty
	$(LATEX) $<

%.odt: %.tex
	sed -e 's/tuftesque/usefonts=false/' < $< > $*.work.tex
	$(MK4HT) $*.work.tex
	mv $*.work.odt $*.odt
	rm -fv $**.tmp $*.{4ct,4tc,dvi,idv,lg,log,xref}
	rm -fv $*.work.*

%.docx: %.odt
	$(SOFFICE) --headless --convert-to docx $<

clean:
	rm -fv *.tmp *.aux *.fdb_latexmk *.fls *.log *.out *.toc *.synctex* *.xdv *.dvi
	rm -fv *.xref *.4ct *.4tc *.4o? *.idv *.lg *.html *.css Act*.png Act*.pdf
	rm -fv *.work.*
	rm -rf .auctex-auto

# So named because there's no guarantee that the odt and docx files
# haven't been hand-patched.
nukefromorbit: clean
	rm -fv *.docx *.odt

%PHONY: all pdf odt clean nukefromorbit
