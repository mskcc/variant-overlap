#!/usr/bin/env Rscript
source("overlap.R")
library(testthat)


test_that("no overlap", {
    variant_list <- list(
        old = c("chr1:1:1", "chr1:2:1"),
        new = c("chr1:3:1", "chr1:4:1")
    )
    overlap <- create_overlaps(variant_list, col_name = "VariantID")
    # VariantID comb
    # 1  chr1:1:1  old
    # 2  chr1:2:1  old
    # 3  chr1:3:1  new
    # 4  chr1:4:1  new

    expect_equal(overlap[["VariantID"]], 
                 factor(c("chr1:1:1", "chr1:2:1", "chr1:3:1", "chr1:4:1"), 
                        levels = c("chr1:1:1", "chr1:2:1", "chr1:3:1", "chr1:4:1")))
    expect_equal(overlap[["comb"]], c("old", "old", "new", "new"))
    expect_equal(dim(overlap), c(4,2))
})



test_that("two variants overlap", {
    variant_list <- list(
        old = c("chr1:1:1", "chr1:2:1"),
        new = c("chr1:1:1", "chr1:2:1")
    )
    overlap <- create_overlaps(variant_list, col_name = "VariantID")
    # VariantID             comb
    # 1  chr1:1:1 overlap: old.new
    # 2  chr1:2:1 overlap: old.new

    expect_equal(overlap[["VariantID"]], 
                 factor(c("chr1:1:1", "chr1:2:1"), 
                        levels = c("chr1:1:1", "chr1:2:1")))
    expect_equal(overlap[["comb"]], c("overlap: old.new", "overlap: old.new"))
    expect_equal(dim(overlap), c(2,2))
})

test_that("two variants overlap, no overlap prefix", {
    variant_list <- list(
        old = c("chr1:1:1", "chr1:2:1"),
        new = c("chr1:1:1", "chr1:2:1")
    )
    overlap <- create_overlaps(variant_list, col_name = "VariantID", overlap_prefix = '')
    # VariantID             comb
    # 1  chr1:1:1 old.new
    # 2  chr1:2:1 old.new
    
    expect_equal(overlap[["VariantID"]], 
                 factor(c("chr1:1:1", "chr1:2:1"), 
                        levels = c("chr1:1:1", "chr1:2:1")))
    expect_equal(overlap[["comb"]], c("old.new", "old.new"))
    expect_equal(dim(overlap), c(2,2))
})


test_that("no overlap prefix and different combination delimeter", {
    variant_list <- list(
        old = c("chr1:1:1", "chr1:2:1"),
        new = c("chr1:1:1", "chr1:2:1")
    )
    overlap <- create_overlaps(variant_list, 
                               col_name = "VariantID", 
                               overlap_prefix = '',
                               comb_delim = '&')
    # VariantID             comb
    # 1  chr1:1:1 old&new
    # 2  chr1:2:1 old&new
    
    expect_equal(overlap[["VariantID"]], 
                 factor(c("chr1:1:1", "chr1:2:1"), 
                        levels = c("chr1:1:1", "chr1:2:1")))
    expect_equal(overlap[["comb"]], c("old&new", "old&new"))
    expect_equal(dim(overlap), c(2,2))
})




test_that("two variants overlap, one doesnt overlap", {
    variant_list <- list(
        old = c("chr1:1:1", "chr1:2:1", "chr1:3:1"),
        new = c("chr1:1:1", "chr1:2:1")
    )
    overlap <- create_overlaps(variant_list, col_name = "VariantID")
    # VariantID             comb
    # 1  chr1:3:1              old
    # 2  chr1:1:1 overlap: old.new
    # 3  chr1:2:1 overlap: old.new
    
    expect_equal(overlap[["VariantID"]], 
                 factor(c("chr1:3:1", "chr1:1:1", "chr1:2:1"),
                        levels = c("chr1:3:1", "chr1:1:1", "chr1:2:1")))
    expect_equal(overlap[["comb"]], c("old", "overlap: old.new", "overlap: old.new"))
    expect_equal(dim(overlap), c(3,2))
})


test_that("two variants overlap, two doesnt overlap", {
    variant_list <- list(
        old = c("chr1:1:1", "chr1:2:1", "chr1:3:1"),
        new = c("chr1:1:1", "chr1:2:1", "chr1:4:1")
    )
    overlap <- create_overlaps(variant_list, col_name = "VariantID")
    # VariantID             comb
    # 1  chr1:3:1              old
    # 2  chr1:4:1              new
    # 3  chr1:1:1 overlap: old.new
    # 4  chr1:2:1 overlap: old.new
    # 
    
    expect_equal(overlap[["VariantID"]], 
                 factor(c("chr1:3:1", "chr1:4:1", "chr1:1:1", "chr1:2:1"), 
                        levels = c("chr1:3:1", "chr1:4:1", "chr1:1:1", "chr1:2:1")))
    expect_equal(overlap[["comb"]], c("old", "new", "overlap: old.new", "overlap: old.new"))
    expect_equal(dim(overlap), c(4,2))
})


test_that("three sets of overlaps", {
    variant_list <- list(
        old1 = c("chr1:1:1", "chr1:2:1", "chr1:3:1", "chr1:3:2"),
        old2 = c("chr1:1:1", "chr1:2:1", "chr1:3:1", "chr1:3:3"),
        new1 = c("chr1:1:1", "chr1:2:1", "chr1:4:1")
    )
    overlap <- create_overlaps(variant_list, col_name = "VariantID")
    # VariantID                    comb
    # 1  chr1:3:2                    old1
    # 2  chr1:3:3                    old2
    # 3  chr1:4:1                    new1
    # 4  chr1:3:1      overlap: old1.old2
    # 5  chr1:1:1 overlap: old1.old2.new1
    # 6  chr1:2:1 overlap: old1.old2.new1
    
    variants <- c("chr1:3:2", "chr1:3:3", "chr1:4:1", "chr1:3:1", "chr1:1:1", "chr1:2:1")
    combs <- c("old1", "old2", "new1", "overlap: old1.old2", "overlap: old1.old2.new1", 
               "overlap: old1.old2.new1")
    expect_equal(overlap[["VariantID"]], 
                 factor(variants, levels = variants))
    expect_equal(overlap[["comb"]], combs)
    expect_equal(dim(overlap), c(6,2))
})


test_that("overlap aggregate table", {
    variant_list <- list(
        old = c("chr1:1:1", "chr1:2:1", "chr1:3:1"),
        new = c("chr1:1:1", "chr1:2:1", "chr1:4:1")
    )
    overlap <- create_overlaps(variant_list, col_name = "VariantID")
    aggr <- aggregate_overlaps(overlap)
    # comb n pcnt
    # 1              new 1   25 
    # 2              old 1   25 
    # 3 overlap: old.new 2   50 

    expect_equal(aggr[["comb"]], 
                 c("new", "old", "overlap: old.new"))
    expect_equal(aggr[["n"]], c(1, 1, 2))
    expect_equal(aggr[["pcnt"]], c(25, 25, 50))
    expect_equal(dim(aggr), c(3,3))
})


# keepcols<-c("Hugo_Symbol", "Entrez_Gene_Id", "Center",
#             "NCBI_Build", "Chromosome", "Start_Position",
#             "End_Position", "Strand", "Variant_Classification",
#             "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1",
#             "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
#             "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "t_depth", "n_depth")



# set up some fixture data for next tests
id_colnames <- c(
    "Chromosome",
    "Start_Position",
    "End_Position", 
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", 
    "Matched_Norm_Sample_Barcode"
)

df1 <- structure(list(
    Hugo_Symbol = c("HPSE2", "GOT1", "SUFU", "SUFU", "SUFU", "SUFU", "SUFU"), 
Entrez_Gene_Id = c("60495", "2805", "51684", "51684", "51684", "51684", "51684"), 
Center = c("mskcc.org", "mskcc.org", "mskcc.org", "mskcc.org", "mskcc.org", "mskcc.org", "mskcc.org"), 
NCBI_Build = c("GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh37"), 
Chromosome = c("10", "10", "10", "10", "10", "10", "10"), 
Start_Position = c("100995756", "101167315", "104309781", "104356995", "104375092", "104375107", "104375147"), 
End_Position = c("100995756", "101167315", "104309781", "104356995", "104375092", "104375107", "104375147"), 
Strand = c("+", "+", "+", "+", "+", "+", "+"), 
Variant_Classification = c("5'Flank", "Intron", "Silent", "Silent", "Missense_Mutation", "Missense_Mutation", "Missense_Mutation"), 
Variant_Type = c("SNP", "SNP", "SNP", "SNP", "SNP", "SNP", "SNP"), 
Reference_Allele = c("C", "G", "G", "C", "C", "G", "C"), 
Tumor_Seq_Allele1 = c("T", "G", "G", "C", "C", "G", "C"), 
Tumor_Seq_Allele2 = c("T", "T", "A", "T", "T", "A", "T"), 
Tumor_Sample_Barcode = c("Sample1", "Sample1", "Sample1", "Sample1", "Sample1", "Sample1", "Sample1"), 
Matched_Norm_Sample_Barcode = c("Sample3", "Sample3", "Sample3", "Sample3", "Sample3", "Sample3", "Sample3"), 
Match_Norm_Seq_Allele1 = c("C", "G", "G", "C", "C", "G", "C"), 
Match_Norm_Seq_Allele2 = c("C", "G", "G", "C", "C", "G", "C"), 
t_depth = c("4", "7", "145", "156", "178", "192", "179"), 
n_depth = c("2", "1", "1168", "1021", "1079", "1119", "1035")), 
row.names = c(1L, 3L, 6L, 7L, 8L, 9L, 10L), class = "data.frame")

df2 <- structure(list(
    Hugo_Symbol = c("HPSE2", "GOT1", "GOT1", "TCEB1P3", "SUFU", "SUFU", "SUFU", "SUFU"), 
Entrez_Gene_Id = c("60495", "2805", "2805", "0", "51684", "51684", "51684", "51684"), 
Center = c("mskcc.org", "mskcc.org", "mskcc.org", "mskcc.org", "mskcc.org", "mskcc.org", "mskcc.org", "mskcc.org"), 
NCBI_Build = c("GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh37", "GRCh37"), 
Chromosome = c("10", "10", "10", "10", "10", "10", "10", "10"), 
Start_Position = c("100995756", "101167306", "101167315", "10216215", "104301397", "104309781", "104356995", "104375107"), 
End_Position = c("100995756", "101167306", "101167315", "10216215", "104301397", "104309781", "104356995", "104375107"), 
Strand = c("+", "+", "+", "+", "+", "+", "+", "+"), 
Variant_Classification = c("5'Flank", "Intron", "Intron", "RNA", "Intron", "Silent", "Silent", "Missense_Mutation"), 
Variant_Type = c("SNP", "SNP", "SNP", "SNP", "SNP", "SNP", "SNP", "SNP"), 
Reference_Allele = c("C", "C", "G", "G", "C", "G", "C", "G"), 
Tumor_Seq_Allele1 = c("T", "C", "G", "G", "C", "G", "C", "G"), 
Tumor_Seq_Allele2 = c("T", "T", "T", "A", "T", "A", "T", "A"), 
Tumor_Sample_Barcode = c("Sample1", "Sample1", "Sample1", "Sample1", "Sample1", "Sample1", "Sample1", "Sample1"), 
Matched_Norm_Sample_Barcode = c("Sample3", "Sample3", "Sample3", "Sample3", "Sample3", "Sample3", "Sample3", "Sample3"), 
Match_Norm_Seq_Allele1 = c("C", "C", "G", "G", "C", "G", "C", "G"), 
Match_Norm_Seq_Allele2 = c("C", "C", "G", "G", "C", "G", "C", "G"), 
t_depth = c("5", "7", "6", "153", "9", "145", "150", "199"), 
n_depth = c("2", "1", "1", "849", "1", "1160", "1025", "1100")), 
row.names = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 9L), 
class = "data.frame")

test_that("get all the variants from a list of data frames merged together", {
    # make list of dataframes per sample
    df_list <- list(
        Sample1 = add_variantID(df1, id_colnames),
        Sample3 = add_variantID(df2, id_colnames))
    
    # make list of variant ID's per sample
    variant_list <- list()
    for(name in names(df_list)){
        variant_list[[name]] <- df_list[[name]][["VariantID"]]
    }
    
    # make list of overlaps of variants by ID
    overlaps <- create_overlaps(variant_list, 
                                col_name = "VariantID", 
                                overlap_prefix = '',
                                comb_delim = '_')
    
    # set of sample overlap combinations
    combinations <- unique(overlaps[["comb"]]) # "Sample1", "Sample3", "Sample1_Sample3"
    
    # get the list of variants for each combination
    sample1_variantIDs <- overlaps[overlaps[["comb"]] == "Sample1", "VariantID"]
    sample2_variantIDs <- overlaps[overlaps[["comb"]] == "Sample3", "VariantID"]
    common_variantIDs <- overlaps[overlaps[["comb"]] == "Sample1_Sample3", "VariantID"]
    
    # get a single df with just the variants from each combination
    sample1_df <- get_all_variants(sample1_variantIDs, df_list, merge_colnames = c("VariantID",id_colnames))
    sample2_df <- get_all_variants(sample2_variantIDs, df_list, merge_colnames = c("VariantID",id_colnames))
    # get a df here that includes one of the value vars we want
    common_df <- get_all_variants(common_variantIDs, df_list, value_var = "t_depth")

    expected_sample1_df <- structure(
        list(VariantID = c("10:104375092:104375092:C:T:Sample1:Sample3",
"10:104375147:104375147:C:T:Sample1:Sample3"), Chromosome = c("10",
"10"), Start_Position = c("104375092", "104375147"), End_Position = c("104375092",
"104375147"), Reference_Allele = c("C", "C"), Tumor_Seq_Allele2 = c("T",
"T"), Tumor_Sample_Barcode = c("Sample1", "Sample1"), Matched_Norm_Sample_Barcode = c("Sample3",
"Sample3")), row.names = c(NA, -2L), class = "data.frame")
    # VariantID Chromosome Start_Position
    # 1 10:104375092:104375092:C:T:Sample1:Sample3         10      104375092
    # 2 10:104375147:104375147:C:T:Sample1:Sample3         10      104375147
    # End_Position Reference_Allele Tumor_Seq_Allele2 Tumor_Sample_Barcode
    # 1    104375092                C                 T              Sample1
    # 2    104375147                C                 T              Sample1
    # Matched_Norm_Sample_Barcode
    # 1                     Sample3
    # 2                     Sample3
    expect_equal(sample1_df, expected_sample1_df)
    
    
    expected_sample2_df <- structure(
        list(VariantID = c("10:101167306:101167306:C:T:Sample1:Sample3",
"10:10216215:10216215:G:A:Sample1:Sample3", "10:104301397:104301397:C:T:Sample1:Sample3"
), Chromosome = c("10", "10", "10"), Start_Position = c("101167306",
"10216215", "104301397"), End_Position = c("101167306", "10216215",
"104301397"), Reference_Allele = c("C", "G", "C"), Tumor_Seq_Allele2 = c("T",
"A", "T"), Tumor_Sample_Barcode = c("Sample1", "Sample1", "Sample1"
), Matched_Norm_Sample_Barcode = c("Sample3", "Sample3", "Sample3"
)), row.names = c(NA, -3L), class = "data.frame")
    # VariantID Chromosome Start_Position
    # 1 10:101167306:101167306:C:T:Sample1:Sample3         10      101167306
    # 2   10:10216215:10216215:G:A:Sample1:Sample3         10       10216215
    # 3 10:104301397:104301397:C:T:Sample1:Sample3         10      104301397
    # End_Position Reference_Allele Tumor_Seq_Allele2 Tumor_Sample_Barcode
    # 1    101167306                C                 T              Sample1
    # 2     10216215                G                 A              Sample1
    # 3    104301397                C                 T              Sample1
    # Matched_Norm_Sample_Barcode
    # 1                     Sample3
    # 2                     Sample3
    # 3                     Sample3
    expect_equal(sample2_df, expected_sample2_df)
    
    expected_common_df <-structure(
        list(VariantID = c("10:100995756:100995756:C:T:Sample1:Sample3",
"10:101167315:101167315:G:T:Sample1:Sample3", "10:104309781:104309781:G:A:Sample1:Sample3",
"10:104356995:104356995:C:T:Sample1:Sample3", "10:104375107:104375107:G:A:Sample1:Sample3"
), t_depthSample1 = c("4", "7", "145", "156", "192"), t_depthSample3 = c("5",
"6", "145", "150", "199")), row.names = c(NA, -5L), class = "data.frame")
    # VariantID t_depthSample1 t_depthSample3
    # 1: 10:100995756:100995756:C:T:Sample1:Sample3              4              5
    # 2: 10:101167315:101167315:G:T:Sample1:Sample3              7              6
    # 3: 10:104309781:104309781:G:A:Sample1:Sample3            145            145
    # 4: 10:104356995:104356995:C:T:Sample1:Sample3            156            150
    # 5: 10:104375107:104375107:G:A:Sample1:Sample3            192            199
    expect_equal(common_df, expected_common_df)

})
