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

    print("Finding intersection of SNPs from all datasets...")
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
    print("Harmonising a list of data...")

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
    stopifnot(check_same(lapply(dat_new, function(x) x$reference_allele)))


    print('Done harmonising.')
    dat_new
}


## Left to be done:
## Generate Jingshu's list of (data, marker.data, cor.mat) using the functions above.
# select SNPs for main data and correlation

selectSNPs <- function(files, max.p.thres = 1, clump_r2 = 0.001, cal.cor = T,
                       p.thres.cor = 0.5, plink_exe = './plink',
                       plink_refdat){
    print("Selecting SNPs for inference and correlation estimation...")
    tmp <- get_snp_intersection(files)

    snp_inter <<- tmp

    # if calculating correlation, use non-significant SNPs (pval > 0.5) without LD clumping
    if (cal.cor) {
        cor.SNPs <- as.character(tmp$SNP[tmp$pval > p.thres.cor]) # (no bonferroni correction here)
    } else {
        cor.SNPs <- NULL
    }


    print("Start clumping using PLINK ...")
    tmp2 <- plink_clump(tmp, plink_exe, refdat = plink_refdat, clump_r2 = clump_r2,
                        clump_p1 = max.p.thres)
    sel.SNPs <- as.character(tmp2$SNP) # selected SNPs


    list(sel.SNPs, cor.SNPs)
}


selectMarkerSNPs <- function(snp_intersection, max.marker.p.thres = 1, clump_r2_formarkers = 0.001, plink_exe = './plink',
                             plink_refdat){
    print("Start clumping using PLINK (markers)")
    mar.SNPs <- plink_clump(snp_intersection, plink_exe, refdat = plink_refdat,
                            clump_r2 = clump_r2_formarkers, clump_p1 = max.marker.p.thres)

    mar.SNPs <- as.character(mar.SNPs$SNP)
    mar.SNPs
}


getInput <- function(sel.files, exp.files, out.files, sel.SNPs = NULL,
                                cor.SNPs = NULL, p.thres.cor = (1-1e-3), cal.cor = TRUE,
                                get.marker.candidates, marker.p.source = "exposure",
                                marker.p.thres = 1e-5,
                                clump_r2_formarkers = 0.05,
                                plink_exe = "./plink",
                                plink_refdat = "./util/data_maf0.01_rs_ref",
                                max.p.thres = 0.01, clump_r2 = 0.001) {

    num_sel <- length(sel.files)
    num_exp <- length(exp.files)
    num_out <- length(out.files)

    files = c(sel.files, exp.files, out.files) # list of all files to be processed

    if (is.null(sel.SNPs)){
        temp <- selectSNPs(files, max.p.thres, clump_r2, p.thres.cor, cal.cor, plink_exe,
                           plink_refdat)
        sel.SNPs <- temp[[1]]
        cor.SNPs <- temp[[2]]

        selectedSNPs <<- sel.SNPs
        selectedcorSNPs <<- cor.SNPs
        print('Selected SNPs for inference data and estimating correlation, extracting data...')
    }


    # get list of SNPs for data and harmonise

    dat <- extract_data(files, sel.SNPs)
    dat <- harmonise_data_list(dat, fast = F)  # list of harmonised data sets

    data <- getData(dat, num_sel, num_exp, num_out)
    print('Inference data completed')
    head(data)

    if (cal.cor) {
        cor_dat <- extract_data(files, cor.SNPs)
        cor_dat <- harmonise_data_list(dat, fast = T)

        cor_dat <- getCorDat(dat, num_sel, num_exp, num_out)
        corr <- calcCorr(cor_dat)
    }

}

#' get data file with meta_data, beta_exp, se_exp, beta_out, se_out, selection_pvals)
#' dat is a list of sel, exp, out dataframes, harmonised.
#'
getData <- function(dat, num_sel, num_exp, num_out) {
    num_files <- num_sel + num_exp + num_out

    meta_data <- dat[[1]][, c("SNP", "effect_allele",
                                               "other_allele")]
    beta_exp <- lapply(dat[(num_sel+1):(num_sel+num_exp)], function(x) x$beta)
    beta_exp <- do.call("rbind", beta_exp)

    se_exp <- lapply(dat[(num_sel+1):(num_sel+num_exp)], function(x) x$se)
    se_exp <- do.call("rbind", se_exp)


    beta_out <- lapply(dat[(num_sel+num_exp+1):num_files], function(x) x$beta)
    beta_out <- do.call("rbind", beta_out)

    se_out <- lapply(dat[(num_sel+num_exp+1):num_files], function(x) x$se)
    se_out <- do.call("rbind", se_out)

    pval <- do.call("rbind",(lapply(dat[1:num_sel], function(x) x$pval)))
    pval <- apply(pval*num_sel, 2, min)
    pval <- pmin(pval, 1)

    # cor_data <- cbind(meta_data$SNP, beta_exp, se_exp, beta_out, se_out)[SNP %in% cor.SNPs]
    #
    #
    # #calculate correlation matrix
    #
    # if (cal.cor) {
    #     z.values <- cbind(cor_data$beta_exp[SNP %in% cor.SNPs]/cor_data$se_exp[SNP %in% cor.SNPs],
    #                       cor_data$beta_out[SNP %in% cor.SNPs]/cor_data$se_out[SNP %in% cor.SNPs])
    #
    #
    #     colnames(z.values) <- c(paste0("exposure", 1:length(exp.files)),
    #                             paste0("outcome", 1:length(out.files)))
    #
    #     z.values <- as.matrix(z.values)
    #     z.values <- z.values[rowSums(is.na(z.values)) == 0,  , drop = F]
    #     covv <- t(z.values) %*% z.values / nrow(z.values)
    #     varr <- colMeans(z.values^2, na.rm = T)
    #     corr <- t(covv / sqrt(varr))/sqrt(varr)
    # } else
    #     corr <- NULL

    #
    data <- cbind(meta_data, beta_exp, se_exp, beta_out, se_out, pval)
    # marker.data <- cbind(meta_data, beta_exp, se_exp, beta_out, se_out, pval)[SNP %in% marker.SNPs]


    data.list <- list(meta_data = meta_data, beta_exp = beta_exp,
                      se_exp = se_exp,
                   beta_out = beta_out, se_out = se_out, selection_pvals = pval)

}

getCorDat <- function(dat, num_sel, num_exp, num_out) {
    num_files <- num_sel + num_exp + num_out

    snp <- dat[[1]][,"SNP"]
    beta_exp <- lapply(dat[(num_sel+1):(num_sel+num_exp)], function(x) x$beta)
    beta_exp <- do.call("rbind", beta_exp)

    se_exp <- lapply(dat[(num_sel+1):(num_sel+num_exp)], function(x) x$se)
    se_exp <- do.call("rbind", se_exp)


    beta_out <- lapply(dat[(num_sel+num_exp+1):num_files], function(x) x$beta)
    beta_out <- do.call("rbind", beta_out)

    se_out <- lapply(dat[(num_sel+num_exp+1):num_files], function(x) x$se)
    se_out <- do.call("rbind", se_out)

    cor_dat <- cbind(snp, beta_exp, se_exp, beta_out, se_out)
    cor_dat
}


#  expecting cor_dat to be a dataframe with columns:
# SNP, beta_exp, se_exp, beta_out, se_out
calcCorr <- function(cor_dat){
    z.values <- cbind(cor_data$beta_exp[SNP %in% cor.SNPs]/cor_data$se_exp[SNP %in% cor.SNPs],
                      cor_data$beta_out[SNP %in% cor.SNPs]/cor_data$se_out[SNP %in% cor.SNPs])


    colnames(z.values) <- c(paste0("exposure", 1:length(exp.files)),
                            paste0("outcome", 1:length(out.files)))

    z.values <- as.matrix(z.values)
    z.values <- z.values[rowSums(is.na(z.values)) == 0,  , drop = F]
    covv <- t(z.values) %*% z.values / nrow(z.values)
    varr <- colMeans(z.values^2, na.rm = T)
    corr <- t(covv / sqrt(varr))/sqrt(varr)
    corr
}
