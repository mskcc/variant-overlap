#!/usr/bin/env Rscript
# wrapper script to compile the report with the given args

input_files <- commandArgs(T)
output_dir <- "output"
output_file <-"report.html"
output_path <- file.path(output_dir, output_file)

rmarkdown::render(
    input="report.Rmd",
    params=list(
        input_files = input_files, 
        output_dir = output_dir),
    output_format = "html_document",
    output_file = output_path)