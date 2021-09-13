## context("Test preprocessing of GWAS summary datasets")
library(GRAPPLE)
library(data.table)

## ## Obtain small datasets for testing

## files <- c("../data/BMI-giant17eu.csv", "../data/CAD-c4d11.csv",
##            "../data/T2D-diagram12-Im.csv")
## tmp <- get_snp_intersection(files)

## snps_kept <- unique(c(subset(tmp, pval < 1e-6)$SNP, sample(tmp$SNP, 2000)))

## data <- fread("../data/BMI-giant17eu.csv")
## data <- data[SNP %in% snps_kept]
## fwrite(data, "../data/BMI-giant17eu-subset.csv")

## data <- fread("../data/CAD-c4d11.csv")
## data <- data[SNP %in% snps_kept]
## fwrite(data, "../data/CAD-c4d11-subset.csv")

## data <- fread("../data/T2D-diagram12-Im.csv")
## data <- data[SNP %in% snps_kept]
## fwrite(data, "../data/T2D-diagram12-Im-subset.csv")

shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
if (shell == "Windows") {
    plink_exe <- "./plink"
} else {
    plink_exe <-"../../util/plink_mac/plink"
}


sel_files <- "../../data/BMI-giant17eu-subset.csv"
exp_files <- "../../data/CAD-c4d11-subset.csv"
out_files <- "../../data/T2D-diagram12-Im-subset.csv"
files <- c(sel_files, exp_files, out_files)

suppressMessages(SNP_inter <- get_SNP_intersection(files))

test_that("Get intersection of SNPs", {
    expect_equal(nrow(SNP_inter), 2644)
})

suppressMessages(
    snps_inference <-
        filter_SNPs(sel_files, exp_files, out_files,
                    SNP_inter,
                    utility = "inference",
                    plink_exe = plink_exe,
                    plink_refdat = "../../util/data_maf0.01_rs")
)

test_that("Filter SNPs using LD clumping", {
    expect_equal(length(snps_inference), 67)
})

suppressMessages(dat_list <- extract_data(files, snps_inference))

test_that("Data extraction", {
    expect_equal(length(dat_list), 3)
    expect_equal(sapply(dat_list, nrow), rep(67, 3))
})

suppressMessages(dat_list <- harmonise_data_list(dat_list))
## Check the data tables in the list have the same SNP, effect_allele, and reference_allele
check_same <- function(x) {
    prod(sapply(x, all.equal, x[[1]]))
}

test_that("Data harmonising", {
    expect_equal(check_same(lapply(dat_list, function(x) x$SNP)), 1)
    expect_equal(check_same(lapply(dat_list, function(x) x$effect_allele)), 1)
    expect_equal(check_same(lapply(dat_list, function(x) x$other_allele)), 1)
})

inference_data <- format_data(dat_list, exp_files = exp_files, out_files = out_files)

minimal_object_names <- c("meta_data", "beta_exp", "se_exp", "beta_out", "se_out")
test_that("Get inference data", {
    expect_equal(intersect(minimal_object_names, names(inference_data)), minimal_object_names)
    expect_equal(nrow(inference_data$meta_data), 67)
})

## data <- data.frame(gamma_exp = inference_data$beta_exp,
##                    gamma_out = inference_data$beta_out,
##                    se_exp = inference_data$se_exp,
##                    se_out = inference_data$se_out)
## grappleRobustEst(data)
