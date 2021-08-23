#' Read GWAS summary statistics data
#'
#' @importFrom tools file_ext
#' @importFrom data.table fread
#'
read_gwas_summary <- function(file, message = TRUE) {

    if (message) {
        message(paste("Reading", file))
    }

    if (file_ext(file) %in% c("rda", ".rData")) {
        load(file)
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
#' files is a list of the file names
#' 
#' Outputs df with columns SNP and pval
#' 
#' 
get_snp_intersection <- function(files, message = TRUE) {

    sel_dat <- read_gwas_summary(files[1], message)
    snp_intersection <- unique(sel_dat$SNP)

    if (length(files) > 1) {
        for (file in files[-1]) {

            dat_new <- read_gwas_summary(file, message)
            snp_intersection <- intersect(snp_intersection, dat_new$SNP)
            
        }
    }

    output <- sel_dat[SNP %in% snp_intersection, c("SNP", "pval")] # store pvals for all files
    output[!duplicated(SNP), ]

}

#' Extract selected SNPs from data
#' Input: list of datafile names, list of SNPs to be extracted
#' 
#' Outputs a list of length(files) dataframes with only the SNPs in SNP_list
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
#' @importFrom TwoSampleMR harmonise_data
#'
harmonise_data_list <- function(dat) {
    
    
    dat_new <- list()
    dat[[1]] <- dat[[1]][order(SNP)]
    SNP_list <- dat[[1]]$SNP
    dat_new[[1]] <- dat[[1]]

    ## Reference dataset
    
    data1 <- formatData(dat[[1]], "exposure")
    if (!("exposure" %in% colnames(data1))) {
            data1[, "exposure"] <- rep("exposure", nrow(data1)) }

    
    stopifnot(length(dat) > 1)

    for (i in 2:length(dat)) {
        ## Format and harmonise with the reference dataset
        data2 <- formatData(dat[[i]], "outcome")
        if (!("outcome" %in% colnames(data2))) {
            data2[, "outcome"] <- rep("outcome", nrow(data2)) }
        
        
        data2 <- suppressMessages(harmonise_data(data1, data2))
        
        

        ## Extract the columns corresponding to the new dataset
        data2 <- data2[data2$mr_keep, ]
        
        data2$SNP <- as.character(data2$SNP)
        
        data2 <- data2[, c(1, grep(".outcome", colnames(data2))), drop = FALSE]
        colnames(data2) <- gsub(".outcome", "", colnames(data2))
        
        # change types back to character after harmonising
        data2$effect_allele <- as.character(data2$effect_allele)
        data2$other_allele <- as.character(data2$other_allele)

        dat_new[[i]] <- data.table(data2)
        SNP_list <- intersect(SNP_list, data2$SNP)
    }

    dat_new <- lapply(dat_new, function(x) x[SNP %in% SNP_list][!duplicated(SNP), ])
    dat_new1 <<- dat_new
    
    ## Check the data tables in the list have the same SNP, effect_allele, and reference_allele
    check_same <- function(x) {
        prod(sapply(x, identical, x[[1]])) == 1
        
    }
    
    
    stopifnot(check_same(lapply(dat_new, function(x) x$SNP)))
    stopifnot(check_same(lapply(dat_new, function(x) x$effect_allele)))
    stopifnot(check_same(lapply(dat_new, function(x) x$reference_allele)))

    dat_new
}


## Left to be done:
## Generate Jingshu's list of (data, marker.data, cor.mat) using the functions above.

# select SNPs for main data and correlation

selectSNPs <- function(files, max.p.thres = 1, clump_r2 = 0.001, p.thres.cor = 0.05, cal.cor = TRUE, plink_exe = './plink',
                       plink_refdat, get.marker.candidates = TRUE, max.marker.p.thres = 1, clump_r2_formarkers = 0.001){
    
    tmp <- get_snp_intersection(files)
    
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
    
    if (get.marker.candidates == TRUE) {
        mar.SNPs <- selectMarkerSNPs(tmp, max.marker.p.thres, clump_r2_formarkers, plink_exe = plink_exe, refdat = plink_refdat)
    }
    
    list(sel.SNPs, cor.SNPs, mar.SNPs)
}


selectMarkerSNPs <- function(snp_intersection, max.marker.p.thres = 1, clump_r2_formarkers = 0.001, plink_exe = './plink',
                             plink_refdat){
    print("Start clumping using PLINK (markers)")
    mar.SNPs <- plink_clump(snp_intersection, plink_exe, refdat = plink_refdat,
                            clump_r2 = clump_r2_formarkers, clump_p1 = max.marker.p.thres)
    
    mar.SNPs <- as.character(mar.SNPs$SNP)
    mar.SNPs
}


#' get data file with meta_data, beta_exp, se_exp, beta_out, se_out, selection_pvals)
#' 
getData <- function(sel.files, exp.files, out.files, sel.SNPs = NULL, 
                    cor.SNPs = NULL, p.thres.cor = (1-1e-3), cal.cor = TRUE,
                    get.marker.candidates, marker.p.source = "exposure",
                    marker.p.thres = 1e-5,
                    clump_r2_formarkers = 0.05
                    plink_exe = "./plink",
                    plink_refdat = "./util/data_maf0.01_rs_ref",
                    max.p.thres = 0.01, clump_r2 = 0.001) {
   

    num_sel <- length(sel.files)
    num_exp <- length(exp.files)
    num_out <- length(out.files)
    num_files <- num_sel + num_exp + num_out
    
    files = c(sel.files, exp.files, out.files) # list of all files to be processed
    
    
     # get intersection of SNPs, do LD clumping and harmonise, get cor SNPs
    if (is.null(sel.SNPs)){
        temp <- selectSNPs(files, max.p.thres,  clump_r2, p.thres.cor, cal.cor, plink_exe,
        plink_refdat)
        sel.SNPs <- temp[[1]]
        cor.SNPs <- temp[[2]]
        marker.SNPs <- temp[[3]]
        print('Selected SNPs for data, marker data and correlation, extracting data...')
    }
    
    # get list of SNPs for data, cor and markers and then harmonise them all together
    
    dat <- extract_data(files, union(sel.SNPs, cor.SNPs))
    print()
    dat <- harmonise_data_list(dat)  # list of harmonised data sets
    
    

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

    cor_data <- cbind(meta_data$SNP, beta_exp, se_exp, beta_out, se_out)[SNP %in% cor.SNPs]
    
    
    #calculate correlation matrix

    if (cal.cor) {
        z.values <- cbind(cor_data$beta_exp[SNP %in% cor.SNPs]/cor_data$se_exp[SNP %in% cor.SNPs],
                          cor_data$beta_out[SNP %in% cor.SNPs]/cor_data$se_out[SNP %in% cor.SNPs])


        colnames(z.values) <- c(paste0("exposure", 1:length(exp.files)),
                                paste0("outcome", 1:length(out.files)))

        z.values <- as.matrix(z.values)
        z.values <- z.values[rowSums(is.na(z.values)) == 0,  , drop = F]
        covv <- t(z.values) %*% z.values / nrow(z.values)
        varr <- colMeans(z.values^2, na.rm = T)
        corr <- t(covv / sqrt(varr))/sqrt(varr)
    } else
        corr <- NULL

    #
    data <- cbind(meta_data, beta_exp, se_exp, beta_out, se_out, pval)[SNP %in% sel.SNPs]
    marker.data <- cbind(meta_data, beta_exp, se_exp, beta_out, se_out, pval)[SNP %in% marker.SNPs]
    
    return(list(data = data, marker.data = marker.data, cor.mat = corr))

    
    
    # data.list <- list(meta_data = meta_data, beta_exp = beta_exp , se_exp = se_exp,
    #              beta_out = beta_out, se_out = se_out, selection_pvals = pval)
    
    
    
    
    
}