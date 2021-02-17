#!/usr/bin/env Rscript
# script to run through the variant overlap workflow
library("ggplot2")
source("overlap.R")
input_files <- commandArgs(T)
output_dir <- "output"
# columns in each table to concatenate to make a unique variant identifier
id_colnames <- c(
    "Chromosome",
    "Start_Position",
    "End_Position", 
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", 
    "Matched_Norm_Sample_Barcode"
    )

# columns with numberic values that should be plotted
metric_colnames <- c(
    't_depth',
    't_alt_count',
    'vcf_qual',
    't_QUAL'
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
    write.table(x = uniq_df, 
                file = output_path, 
                quote = FALSE, 
                sep = '\t', 
                row.names = FALSE,
                col.names = TRUE)
}
save.image()

# get the labels that are only for combination categories
combination_labels <- unique(overlaps[["comb"]])
combination_labels <- combination_labels[! combination_labels %in% names(input_list)]

# NOTE: I think I only ever tested this with 2 input files, not sure if it works with multiple input files
# TODO: how to handle outputting the variants from the overlap combinations?
# get the variants from all other combinations
for(name in combination_labels ){
    # output filepath
    output_filename <- sprintf("common.%s.tsv", basename(name))
    output_path <- file.path(output_dir, output_filename)
    
    # get all the variants for this combination
    common_variant_IDs <- overlaps[overlaps[["comb"]] == name, "VariantID"]
    df <- get_all_variants(
        variant_IDs = common_variant_IDs, 
        df_list = input_list, 
        merge_colnames = c(id_colnames, "VariantID"))
    overlap_variant_tables[[name]] <- df[, c(id_colnames, "VariantID")]
    
    # save output file
    write.table(x = df, 
                file = output_path, 
                quote = FALSE, 
                sep = '\t', 
                row.names = FALSE, 
                col.names = TRUE)
}
save.image()

# save a plot of the overlap counts; horizontal bar plot Venn
# NOTE: consider using Upset plot here in the future
output_pdf <- file.path(output_dir, "overlap.pdf")
pdf(file = output_pdf)
overlap_barplot(aggr_table, plot_title = "Number of overlaps per combination")
dev.off()


# make a differential plot on the metrics in the overlaps
# TODO: integrate this with the other loops
for (metric in metric_colnames){
    for(name in combination_labels ){
        common_variant_IDs <- overlaps[overlaps[["comb"]] == name, "VariantID"]
        df <- get_all_variants(
            variant_IDs = common_variant_IDs, 
            df_list = input_list, 
            value_var = metric)
        
        # these will be the new columns in the output df; t_QUALexamples.Sample1.maf, t_QUALexamples.Sample2.maf, etc
        labels <- paste0(metric, names(input_list))
        # subtract the values; new from old
        # NOTE: this only uses the first two items in theinputs; how to set this up for >2 inputs?
        diff_label <- paste0(metric, "_diff")
        df[[diff_label]] <- as.numeric(df[,labels[2]]) - as.numeric(df[,labels[1]])
        
        
        p <- ggplot(data = df, aes(x = VariantID, y = get(diff_label))) + 
            geom_point(aes(color = get(diff_label)))
        
        output_file <- file.path(output_dir, paste0(diff_label, ".pdf"))
        pdf(output_file)
        # NOTE: do not forget that you need to use print with a plot in a loop...
        print(p)
        dev.off()
        
    }
}


output_Rdata <- file.path(output_dir, ".RData")
save.image(output_Rdata)

