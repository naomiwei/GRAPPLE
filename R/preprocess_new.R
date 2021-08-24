loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


#' Read GWAS summary statistics data
#'
#' @importFrom tools file_ext
#' @importFrom data.table fread
#'
#'

read_gwas_summary <- function(file, message = TRUE) {

    if (message) {
        message(paste("Reading", file))
    }

    if (file_ext(file) %in% c("rda", ".rData")) {
        dat <- loadRData(file)
    } else {
        dat <- fread(file)
    }
    dat
}

#' Get intersection of the SNPs in GWAS summary datasets
#'
#' @import data.table
#'
#' @details The first file is the selection dataset.
#'
get_snp_intersection <- function(files, message = TRUE) {

    message("Finding intersection of SNPs from all datasets...")
    sel_dat <- read_gwas_summary(files[1], message)
    snp_intersection <- unique(sel_dat$SNP)

    if (length(files) > 1) {
        for (file in files[-1]) {

            dat_new <- read_gwas_summary(file, message)
            snp_intersection <- intersect(snp_intersection, dat_new$SNP)
        }
    }

    output <- sel_dat[SNP %in% snp_intersection, c("SNP", "pval")]
    output[!duplicated(SNP), ]

}

#' Extract selected SNPs from data
#'
#' @examples
#' files <- file.path("../data", list.files("../data"))
#' tmp <- get_snp_intersection(files)
#' tmp2 <- plink_clump(tmp, "../util/plink_mac/plink", refdat = "../util/data_maf0.01_rs")
#' dat <- extract_data(files, tmp2$SNP)
#' dat <- harmonise_data_list(dat)
#'
extract_data <- function(files, SNP_list = NULL, message = TRUE) {

    dat <- list()
    for (i in 1:length(files)) {
        dat[[i]] <- read_gwas_summary(files[i], message = message)
        dat[[i]] <- dat[[i]][SNP %in% SNP_list]
        dat[[i]] <- dat[[i]][!duplicated(SNP), ]
    }
    dat
}

#' Harmonise a list of GWAS summary data tables
#'
#' Input: dat, a list of dataframes
#' Output: dat_new, a list of (harmonised) dataframes
#'
#' mode = 1 for TwoSampleMR hormonise_data
#' mode = 2 for GRAPPLE harmonise_data
#'
#' @importFrom TwoSampleMR harmonise_data
#'
harmonise_data_list <- function(dat, fast = T) {
    message("Harmonising a list of data...")

    dat_new <- list()
    SNP_list <- dat[[1]]$SNP
    dat_new[[1]] <- dat[[1]]

    ## Reference dataset
    data1 <- formatData(dat[[1]], "exposure")

    stopifnot(length(dat) > 1)

    for (i in 2:length(dat)) {
        ## Format and harmonise with the reference dataset
        data2 <- formatData(dat[[i]], "outcome")
        if (!("outcome" %in% colnames(data2))) {
            data2[, "outcome"] <- rep("outcome", nrow(data2)) }

        if (fast == T) {
            data2 <- suppressMessages(harmonise_data1(data1, data2))
        } else {
            data2 <- suppressMessages(harmonise_data(data1, data2))
        }



        ## Extract the columns corresponding to the new dataset
        data2 <- data2[data2$mr_keep, ]

        data2$SNP <- as.character(data2$SNP)

        data2 <- data2[, c(1, grep(".outcome", colnames(data2))), drop = FALSE]

        colnames(data2) <- gsub(".outcome", "", colnames(data2))

        dat_new[[i]] <- data.table(data2)
        SNP_list <- intersect(SNP_list, data2$SNP)
    }

    dat_new <- lapply(dat_new, function(x) x[SNP %in% SNP_list][!duplicated(SNP), ])

    ## Check the data tables in the list have the same SNP, effect_allele, and reference_allele
    check_same <- function(x) {
        prod(sapply(x, all.equal, x[[1]])) == 1
    }
    stopifnot(check_same(lapply(dat_new, function(x) x$SNP)))
    stopifnot(check_same(lapply(dat_new, function(x) x$effect_allele)))
    stopifnot(check_same(lapply(dat_new, function(x) x$other_allele)))


    message('Done harmonising.')
    dat_new
}


## Left to be done:
## Generate Jingshu's list of (data, marker_data, cor_mat) using the functions above_
# select SNPs for main data and correlation

selectSNPs <- function(files, max_p_thres = 1, clump_r2 = 0.001, cal_cor = T,
                       p_thres_cor = 0.5, plink_exe = './plink',
                       plink_refdat){
    message("Selecting SNPs for inference and correlation estimation...")
    tmp <- get_snp_intersection(files)

    snp_inter <<- tmp

    # if calculating correlation, use non-significant SNPs (pval > 0.5) without LD clumping
    if (cal_cor) {
        cor_SNPs <- as.character(tmp$SNP[tmp$pval > p_thres_cor]) # (no bonferroni correction here)
    } else {
        cor_SNPs <- NULL
    }


    message("Start clumping using PLINK ...")
    tmp2 <- plink_clump(tmp, plink_exe, refdat = plink_refdat, clump_r2 = clump_r2,
                        clump_p1 = max_p_thres)
    sel_SNPs <- as.character(tmp2$SNP) # selected SNPs


    list(sel_SNPs, cor_SNPs)
}


selectMarkerSNPs <- function(snp_intersection, max_marker_p_thres = 1, clump_r2_formarkers = 0.001, plink_exe = './plink',
                             plink_refdat){
    message("Start clumping using PLINK (markers)")
    mar_SNPs <- plink_clump(snp_intersection, plink_exe, refdat = plink_refdat,
                            clump_r2 = clump_r2_formarkers, clump_p1 = max_marker_p_thres)

    mar_SNPs <- as.character(mar_SNPs$SNP)
    mar.SNPs
}


getInput <- function(sel_files, exp_files, out_files, sel_SNPs = NULL,
                     cor_SNPs = NULL, p_thres_cor = (1-1e-3), cal_cor = TRUE,
                     get_marker_candidates, marker_p_source = "exposure",
                     marker_p_thres = 1e-5,
                     clump_r2_formarkers = 0.05,
                     plink_exe = "./plink",
                     plink_refdat = "./util/data_maf0.01_rs_ref",
                     max_p_thres = 0.01, clump_r2 = 0.001) {

    num_sel <- length(sel_files)
    num_exp <- length(exp_files)
    num_out <- length(out_files)

    files = c(sel_files, exp_files, out_files) # list of all files to be processed

    if (is.null(sel_SNPs)){
        temp <- selectSNPs(files, max_p_thres, clump_r2, p_thres_cor, cal_cor, plink_exe,
                           plink_refdat)
        sel_SNPs <- temp[[1]]
        cor_SNPs <- temp[[2]]

        selectedSNPs <- sel_SNPs
        selectedcorSNPs <- cor_SNPs
        message('Selected SNPs for inference data and estimating correlation, extracting data...')
    }


    # get list of SNPs for data and harmonise

    dat <- extract_data(files, sel_SNPs)
    dat <- harmonise_data_list(dat, fast = F)  # list of harmonised data sets
    dat$num_sel <- num_sel
    dat$num_exp <- num_exp
    dat$num_out <- num_out

    data <- get_data(dat)
    message('Inference data completed')
    head(data)

    if (cal_cor) {
        cor_dat <- extract_data(files, cor_SNPs)
        cor_dat <- harmonise_data_list(dat, fast = T)

        cor_dat <- get_data(dat, num_sel, num_exp, num_out)
        corr <- calcCorr(cor_dat)
    }

}

#' get data file with meta_data, beta_exp, se_exp, beta_out, se_out, selection_pvals)
#' dat is a list of sel, exp, out dataframes, harmonised.
#'
get_data <- function(dat_list) {

    list2env(dat_list[c("num_sel", "num_exp", "num_out")], environment())

    tryCatch(num_files <- num_sel + num_exp + num_out,
             error = function(e) {
                 stop("At least one of num_sel, num_exp, and num_out are not available in dat_list.")
             })

    meta_data <- dat_list[[1]][, c("SNP", "effect_allele", "other_allele")]
    beta_exp <- lapply(dat_list[(num_sel+1):(num_sel+num_exp)], function(x) x$beta)
    beta_exp <- do.call("rbind", beta_exp)

    se_exp <- lapply(dat_list[(num_sel+1):(num_sel+num_exp)], function(x) x$se)
    se_exp <- do.call("rbind", se_exp)


    beta_out <- lapply(dat_list[(num_sel+num_exp+1):num_files], function(x) x$beta)
    beta_out <- do.call("rbind", beta_out)

    se_out <- lapply(dat_list[(num_sel+num_exp+1):num_files], function(x) x$se)
    se_out <- do.call("rbind", se_out)

    pval <- do.call("rbind",(lapply(dat_list[1:num_sel], function(x) x$pval)))
    pval <- apply(pval*num_sel, 2, min)
    pval <- pmin(pval, 1)

    list(meta_data = meta_data,
         beta_exp = beta_exp,
         se_exp = se_exp,
         beta_out = beta_out,
         se_out = se_out,
         selection_pvals = pval)

}

#  expecting cor_dat to be a dataframe with columns:
# SNP, beta_exp, se_exp, beta_out, se_out
calcCorr <- function(cor_data, cor_SNPs){
    z_values <- cbind(cor_data$beta_exp[SNP %in% cor_SNPs]/cor_data$se_exp[SNP %in% cor_SNPs],
                      cor_data$beta_out[SNP %in% cor_SNPs]/cor_data$se_out[SNP %in% cor_SNPs])


    colnames(z_values) <- c(paste0("exposure", 1:length(exp_files)),
                            paste0("outcome", 1:length(out_files)))

    z_values <- as_matrix(z_values)
    z_values <- z_values[rowSums(is.na(z_values)) == 0,  , drop = F]
    covv <- t(z_values) %*% z_values / nrow(z_values)
    varr <- colMeans(z_values^2, na.rm = T)
    corr <- t(covv / sqrt(varr))/sqrt(varr)
    corr
}
