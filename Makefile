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
	r::r-testthat==2.1.1 \
	conda-forge::r-devtools \
	bioconda::bcftools==1.9 \
	bioconda::vcf2maf==1.6.19 \
	bioconda::gatk==3.8


# Commands for converting bcf and vcf into tsv for easier parsing in script
bcf2vcf:
	bcftools view input.bcf > input.vcf

BCF_DIR:=input
convert_all_bcfs:
	find "$(BCF_DIR)/" -type f -name "*.bcf" | \
	while read i; do \
	bcftools view "$${i}" > "$${i%%.bcf}.vcf"; \
	done

hg19_fa:=/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta
hg19_fai:=/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai
hg19_dict:=/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict
hg19_index:=/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.index
hg19_ms:=/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.microsatellites.list
hg19_excl:=/juno/work/taylorlab/cmopipeline/mskcc-igenomes/grch37/sv_calling/human.hg19.excl.tsv

# sed -i -e 's|9e0009d1-c993-4247-9706-88ee84591dec|TUMOR|g' DO220841.INV.samplenames.vcf
# sed -i -e 's|78960b37-34e1-4093-8cc1-d3fc639805e5|NORMAL|g' DO220841.INV.samplenames.vcf
# --tumor-id ${idTumor} \
# --normal-id ${idNormal} \
# --vcf-tumor-id ${idTumor} \
# --vcf-normal-id ${idNormal} \
# --retain-info ${infoCols} \
# --retain-fmt ${formatCols} \
# --custom-enst ${isoforms} \
# --filter-vcf 0
TUMOR_ID:=9e0009d1-c993-4247-9706-88ee84591dec
NORMAL_ID:=78960b37-34e1-4093-8cc1-d3fc639805e5
vcf2maf:
	vcf2maf.pl \
	--inhibit-vep \
	--input-vcf DO220841.INV.vcf \
	--output-maf DO220841.INV.maf \
	--tumor-id "$(TUMOR_ID)" \
	--normal-id "$(NORMAL_ID)" \
	--vcf-tumor-id "$(TUMOR_ID)" \
	--vcf-normal-id "$(NORMAL_ID)" \
	--retain-info CHROM,POS,ID,REF,ALT,QUAL,FILTER,CIEND,CIPOS,CHR2,MAPQ,SR,SVTYPE,INSLEN \
	--retain-fmt GT,GQ,RC,RCL,RCR,CN \
	--ref-fasta /juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta

VCF_DIR:=input
VCF_FILES:=$(shell find "$(VCF_DIR)/" -type f -name "*.vcf")
$(VCF_FILES):
	@vcf_file="$$(echo $@)"
	maf_file="$${vcf_file%%.vcf}.maf"
	echo "$$maf_file"
	vcf2maf.pl \
	--inhibit-vep \
	--input-vcf "$${vcf_file}" \
	--output-maf "$${maf_file}" \
	--tumor-id "$(TUMOR_ID)" \
	--normal-id "$(NORMAL_ID)" \
	--vcf-tumor-id "$(TUMOR_ID)" \
	--vcf-normal-id "$(NORMAL_ID)" \
	--retain-info CHROM,POS,ID,REF,ALT,QUAL,FILTER,CIEND,CIPOS,CHR2,MAPQ,SR,SVTYPE,INSLEN \
	--retain-fmt GT,GQ,RC,RCL,RCR,CN \
	--ref-fasta /juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta
vcf2maf_all: $(VCF_FILES)
.PHONY: $(VCF_FILES)

GATK_JAR:=conda/opt/gatk-3.8/GenomeAnalysisTK.jar
varaints_to_table:
	java -Xms8G -Xmx8G -jar "$(GATK_JAR)" -T VariantsToTable \
	-R "$(hg19_fa)" \
	-V DO220841.INV.vcf \
	-F CHROM \
	-F POS \
	-F ID \
	-F REF \
	-F ALT \
	-F QUAL \
	-F FILTER \
	-F CIEND \
	-F CIPOS \
	-F CHR2 \
	-F MAPQ \
	-F SR \
	-F SVTYPE \
	-F INSLEN \
	-GF GT \
	-GF GQ \
	-GF RC \
	-GF RCL \
	-GF RCR \
	-GF CN \
	-o "DO220841.INV.tsv"
# command for SNPs and variant calls
# java -Xms8G -Xmx8G -jar "$(GATK_JAR)" -T VariantsToTable \
# -R "$(hg19_fa)" \
# -V DO220841.INV.vcf \
# -F CHROM \
# -F POS \
# -F ID \
# -F REF \
# -F ALT \
# -F QUAL \
# -F FILTER \
# -GF DP \
# -GF AD \
# -GF RD \
# -GF FREQ \
# -GF RBQ \
# -GF ABQ \
# -o "DO220841.INV.tsv"


# run the report
OLD_MAF:=examples/Sample1.maf
NEW_MAF:=examples/Sample2.maf
INPUT_FILES:=$(OLD_MAF) $(NEW_MAF)
OUTPUT_DIR:=output
$(OUTPUT_DIR):
	mkdir -p "$(OUTPUT_DIR)"
run: $(OUTPUT_DIR)
	./run.R $(INPUT_FILES)
# R --vanilla <<<'input_files <- commandArgs()[3:length(commandArgs())]
# rmarkdown::render(
# input="report.Rmd",
# params=list(input_files=input_files, output_dir = "$(OUTPUT_DIR)"),
# output_format="html_document",
# output_file="report.html")
# ' \
# $(INPUT_FILES)
# html_output <- rmarkdown::html_document(mathjax=NULL)

# run the test cases
test:
	Rscript -e 'testthat::test_dir(".")'

# interactive session
bash:
	bash
