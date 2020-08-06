export SHELL:=/bin/bash
.ONESHELL:
UNAME:=$(shell uname)

define help
Help message goes here
endef
export help
help:
	@printf "$$help"
.PHONY : help

# ~~~~~ Install Dependencies ~~~~~ #
# need to use Python 2.7 because Python 3 gives different results
export PATH:=$(CURDIR)/conda/bin:$(CURDIR)/bin:$(PATH)
unexport PYTHONPATH
unexport PYTHONHOME

ifeq ($(UNAME), Darwin)
CONDASH:=Miniconda3-4.7.12.1-MacOSX-x86_64.sh
endif

ifeq ($(UNAME), Linux)
CONDASH:=Miniconda3-4.7.12.1-Linux-x86_64.sh
endif

CONDAURL:=https://repo.anaconda.com/miniconda/$(CONDASH)

conda:
	echo ">>> Setting up conda..."
	wget "$(CONDAURL)"
	bash "$(CONDASH)" -b -p conda
	rm -f "$(CONDASH)"

install: conda
	conda install -y \
	r-base==3.6.1 \
	r::r-rmarkdown==1.12 \
	r::r-ggplot2==3.1.1 \
	r::r-knitr==1.22 \
	r::r-dt==0.5 \
	conda-forge::pandoc==2.10.1


bash:
	bash

run:
	R --vanilla <<<'rmarkdown::render(input="report.Rmd", params=list(foo="baz"), output_format="html_document", output_file="report.html")'

# R --vanilla <<<'print("foo")'
# R --vanilla <<E0F
# rmarkdown::render(
#     input = "compare.Rmd",
# params = list(
#     old_unpaired_annot = "old.unpaired.annot.tsv",
#     new_unpaired_annot = "new.unpaired.annot.tsv",
#     old_paired_annot = "old.paired.annot.tsv",
#     new_paired_annot = "new.paired.annot.tsv",
#     old_unpaired_annot_path = "${old_unpaired_annot_path}",
#     new_unpaired_annot_path = "${new_unpaired_annot_path}",
#     old_paired_annot_path = "${old_paired_annot_path}",
#     new_paired_annot_path = "${new_paired_annot_path}"
# ),
# output_format = "html_document",
# output_file = "${html_output}"
# )
# E0F
