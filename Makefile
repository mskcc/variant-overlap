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
	conda-forge::pandoc==2.10.1 \
	r::r-testthat==2.1.1


bash:
	bash

INPUT_FILES:=analyst_file.old.maf analyst_file.new.maf
OUTPUT_DIR:=output
$(OUTPUT_DIR):
	mkdir -p "$(OUTPUT_DIR)"
run: $(OUTPUT_DIR)
	R --vanilla <<<'input_files <- commandArgs()[3:length(commandArgs())]
	rmarkdown::render(
	input="report.Rmd",
	params=list(input_files=input_files, output_dir = "$(OUTPUT_DIR)"),
	output_format="html_document",
	output_file="report.html")
	' \
	$(INPUT_FILES)
# html_output <- rmarkdown::html_document(mathjax=NULL)
test:
	./test_overlap.R
