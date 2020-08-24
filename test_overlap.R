#!/usr/bin/env Rscript
source("overlap.R")
library(testthat)


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

test_that("two variants overlap, two doesnt overlap", {
    variant_list <- list(
        old = c("chr1:1:1", "chr1:2:1", "chr1:3:1"),
        new = c("chr1:1:1", "chr1:2:1", "chr1:4:1")
    )
    overlap <- create_overlaps(variant_list, col_name = "VariantID")
    aggr <- aggregate_overlaps(overlap)
    # comb n pcnt all
    # 1              new 1   25   .
    # 2              old 1   25   .
    # 3 overlap: old.new 2   50   .
    # 

    expect_equal(aggr[["comb"]], 
                 c("new", "old", "overlap: old.new"))
    expect_equal(aggr[["n"]], c(1, 1, 2))
    expect_equal(aggr[["pcnt"]], c(25, 25, 50))
    expect_equal(aggr[["all"]], c('.', '.', '.'))
    expect_equal(dim(aggr), c(3,4))
})


