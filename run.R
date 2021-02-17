#!/usr/bin/env Rscript
# script to run through the variant overlap workflow
library("ggplot2")
source("overlap.R")
input_files <- commandArgs(T)
output_dir <- "output"
id_colnames <- c(
    "Chromosome",
    "Start_Position",
    "End_Position", 
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", 
    "Matched_Norm_Sample_Barcode"
    )

# load datasets
input_list <- list()
for(i in input_files){
    # need to make a label for each input file; remove any slashes that might come from filepaths
    label <- gsub('/', '.', i)
    input_list[[label]] <- load_variants(i, id_colnames)
}
save.image()

# create variant lists
variant_list <- list()
for(name in names(input_list)){
    variant_list[[name]] <- input_list[[name]][["VariantID"]]
}
save.image()

# overlap the datasets
overlaps <- create_overlaps(variant_list, 
                            col_name = "VariantID", 
                            overlap_prefix = '',
                            comb_delim = '_')
save.image()


# aggregate results metrics
aggr_table <- aggregate_overlaps(overlaps)
output_path <- file.path(output_dir, "overlap.tsv")
write.table(x = aggr_table, file = output_path, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
save.image()

# subset the dataframes for the overlaps and save output files
overlap_variant_tables <- list()

# get the variants unique to each file; comb == name
for(name in names(input_list)){
    # output filepath
    output_filename <- sprintf("uniq.%s.tsv", basename(name))
    output_path <- file.path(output_dir, output_filename)
    # get the combination label
    uniq_variantID_df <- overlaps[overlaps[["comb"]] == name ,]
    # subset the dataset for only that label
    uniq_df <- merge(x = input_list[[name]], y = uniq_variantID_df, by = "VariantID")
    overlap_variant_tables[[name]] <- uniq_df[, c(id_colnames, "VariantID")]
    # save to file
    write.table(x = uniq_df, file = output_path, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}
save.image()


# NOTE: I think I only ever tested this with 2 input files, not sure if it works with multiple input files
# TODO: how to handle outputting the variants from the overlap combinations?
# get the variants from all other combinations
for(name in unique(overlaps[["comb"]])){
    if(! name %in% names(input_list)){
        # output filepath
        output_filename <- sprintf("common.%s.tsv", basename(name))
        output_path <- file.path(output_dir, output_filename)
        
        # get all the variants for this combination
        common_variant_IDs <- overlaps[overlaps[["comb"]] == name, "VariantID"]
        df <- Reduce(function(x, y){
            merge(x[x[["VariantID"]] %in% common_variant_IDs, c(id_colnames, "VariantID")], 
                  y[y[["VariantID"]] %in% common_variant_IDs, c(id_colnames, "VariantID")], 
                  all = TRUE) # by = "VariantID",
        }, input_list)
        overlap_variant_tables[[name]] <- df[, c(id_colnames, "VariantID")]
        write.table(x = df, 
                    file = output_path, 
                    quote = FALSE, 
                    sep = '\t', 
                    row.names = FALSE, 
                    col.names = TRUE)
    }
}
save.image()

output_pdf <- file.path(output_dir, "overlap.pdf")
pdf(file = output_pdf)
overlap_barplot(aggr_table, plot_title = "Number of overlaps per combination")
dev.off()


output_Rdata <- file.path(output_dir, ".RData")
save.image(output_Rdata)

