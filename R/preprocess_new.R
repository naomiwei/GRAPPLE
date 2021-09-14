#' Loads an R data file, and returns it
#'
#' @param file Name of the R data file (.rda or .RData).
#' @return First object in the R data file.
#'
#' @keywords internal
#'
load_R_data <- function(file){
    #
    load(file_name)
    object_names <- ls()[ls() != "file_name"]
    if (length(object_names) > 1) {
        warning(paste0("More than one object in",
                       file,
                       ". Only loading the first object called ",
                       object_names[1]))
    }
    get(object_names[1])
}

#' Read one GWAS summary statistics dataset
#'
#' @param file Name of the GWAS summary statistics file
#' @param message Logical indicator for message printing
#'
#' @return A \code{data.table} object.
#'
#' @details The dataset must contain the following columns: SNP, effect_allele, other_allele, beta, se. When there are duplicated SNPs, only the first entry is kept.
#'
#' @importFrom tools file_ext
#' @importFrom data.table fread data.table
#'
#' @keywords internal
#'
read_gwas_summary <- function(file, message = TRUE) {

    if (message) {
        message(paste("Reading", file))
    }

    if (file_ext(file) %in% c("rda", ".rData")) {
        dat <- data.table(load_R_data(file))
    } else {
        dat <- fread(file)
    }

    ## Check necessary columns
    necessary_columns <- c("SNP", "effect_allele", "other_allele", "beta", "se")
    missing_columns <- setdiff(necessary_columns, names(dat))
    if (length(missing_columns) > 1){
        stop(paste0("The following columns are missing in ", file, ":", missing_columns, "."))
    }

    ## Compute p-values if not already present
    if (!("pval" %in% names(dat))) {
        dat$pval <- pnorm(- abs(dat$beta / dat$se)) * 2
    }

    ## Only return data columns that GRAPPLE is going to use
    output_columns <- c(necessary_columns, "pval")
    dat[!duplicated(dat$SNP), ..output_columns]
}

#' Get the intersection of the SNPs in several GWAS summary datasets
#'
#' @param files A vector of file names
#' @param file_labels A vector of labels the same length as \code{files} used to label the different p-value columns. If null, the original file names are used.
#' @param message Logical indicator for message printing
#'
#' @import data.table
#'
#' @return A \code{data.table} object with the SNP name and p-values in all summary datasets.
#'
#' @export
#'
get_SNP_intersection <- function(files, file_labels = NULL, message = TRUE) {

    message("Finding intersection of SNPs from all datasets...")
    dat <- read_gwas_summary(files[1], message)
    dat <- dat[, c("SNP", "pval")]

    if (is.null(file_labels)) {
        labels <- files
    } else {
        labels <- file_labels
    }
    
    names(labels) <- files
    
    if (length(files) > 1) {
        for (file in files[-1]) {

            dat_new <- read_gwas_summary(file, message)

            dat <- merge(dat, dat_new[, c("SNP", "pval")],
                         by = "SNP",
                         suffixes = c("", paste0("_", labels[file])))
        }
    }

    names(dat)[names(dat) == "pval"] <- paste0("pval_", labels[1])
    dat

}

#' Extract selected SNPs from data
#'
#' @param files A vector of file names.
#' @param SNP_list If not null, only returns SNPs in this list; otherwise returns the original dataset.
#' @param message Logical indicator for message printing.
#'
#' @return A list of data tables, each corresponding to one file.
#'
#' @export
#'
extract_data <- function(files, SNP_list = NULL, message = TRUE) {

    dat <- list()
    for (i in 1:length(files)) {
        dat[[i]] <- read_gwas_summary(files[i], message = message)
        if (!is.null(SNP_list)) {
            dat[[i]] <- dat[[i]][SNP %in% SNP_list]
        }
    }
    names(dat) <- files
    dat
}

#' Harmonise a list of GWAS summary data tables
#'
#' @param dat A list of \code{data.table} objects containing GWAS summary data (typically returned by \code{extract_data})
#' @param fast Whether to use a fast implementation of GWAS data harmonisation.
#'
#' @details When \code{fast} is \code{TRUE}, this function uses a modification of the \emph{TwoSampleMR::harmonise_data} that skips certain checks.
#'
#' @importFrom TwoSampleMR harmonise_data
#'
#' @export
#'
harmonise_data_list <- function(dat_list, fast = T) {
    message("Harmonising a list of data...")

    dat_new <- list()
    
    SNP_list <- dat_list[[1]]$SNP
    dat_new[[1]] <- dat_list[[1]]

    ## Reference dataset
    data1 <- formatData(dat_list[[1]], "exposure")

    stopifnot(length(dat_list) > 1)

    for (i in 2:length(dat_list)) {

        ## Format and harmonise with the reference dataset
        data2 <- formatData(dat_list[[i]], "outcome")

        if (fast == T) {
            data2 <- suppressMessages(harmonise_data_fast(data1, data2))
        } else {
            if (!("outcome" %in% colnames(data2))) {
                data2[, "outcome"] <- rep("outcome", nrow(data2))
            }
            data2 <- suppressMessages(harmonise_data(data1, data2))
        }

        ## Extract the columns corresponding to the new dataset
        data2 <- data2[data2$mr_keep, ]

        data2$SNP <- as.character(data2$SNP)

        data2 <- data2[, c(1, grep(".outcome", colnames(data2))), drop = FALSE]

        colnames(data2) <- gsub(".outcome", "", colnames(data2))
        data2 <- data.table(data2)
        dat_new[[i]] <- data2
        
        SNP_list <- intersect(SNP_list, data2$SNP)
    }

    dat_new <- lapply(dat_new, function(x) x[SNP %in% SNP_list][!duplicated(SNP), ])
    
    # Order SNPs consistently and make sure effect_allele and other_alle
    dat_new <- lapply(dat_new, function(x) x[order(SNP)])
    dat_new <- lapply(dat_new, 
                      function(x) x[,(c("effect_allele", "other_allele")):= lapply(.SD, as.character), .SDcols = c("effect_allele", "other_allele")])

    ## Check the data tables in the list have the same SNP, effect_allele, and reference_allele
    check_same <- function(x) {
        prod(sapply(x, all.equal, x[[1]])) == 1
    }
    
    
    stopifnot(check_same(lapply(dat_new, function(x) x$SNP)))
    stopifnot(check_same(lapply(dat_new, function(x) x$effect_allele)))
    stopifnot(check_same(lapply(dat_new, function(x) x$other_allele)))

    message('Done harmonising.')
    names(dat_new) <- names(dat_list)
    dat_new
}


## Left to be done:
## Generate Jingshu's list of (data, marker_data, cor_mat) using the functions above_

# select SNPs for main data and correlation


#' Select SNPs by p-value filtering and LD clumping
#'
#' @param files A vector of file names
#' @param max_pval_files Files used for upper-bounding p-values
#' @param max_pval_thres Threshold to upper-bound p-values.
#' @param min_pval_files Files used for lower-bounding p-values
#' @param min_pval_thres Threshold to lower-bound p-values.
#' @param clump Whether LD clumping should be performed
#' @param clump_pval_files Files to be used to compute minimum p-value for clumping.
#'
#' @inheritParams plink_clump
#' @inheritParams filter_SNPs
#'
#' @details When upper-bounding the p-values, only one p-value needs to be below \code{max_pval_thres}. When lower-bounding the p-values, all p-values need to be above \code{min_pval_thres}.
#'
#' @keywords internal
#'
filter_SNPs_basic <- function(files = NULL,
                              SNP_inter = NULL,
                              max_pval_files = files,
                              max_pval_thres = 1,
                              min_pval_files = files,
                              min_pval_thres = 0,
                              clump = FALSE,
                              clump_pval_files = files,
                              clump_r2 = 0.001,
                              plink_exe = NULL,
                              plink_refdat = NULL) {

    if (clump && (is.null(plink_exe) || is.null(plink_refdat))) {
        stop("Executable or reference data file for PLINK clumping is not provided.")
    }

    if (is.null(SNP_inter)) {
        SNP_inter <- get_SNP_intersection(files)
    }

    SNP_inter <- data.table(SNP_inter)

    max_pval_names <- paste0("pval_", max_pval_files)
    min_pval_in_max_pval_files <- apply(SNP_inter[, ..max_pval_names], 1, min) # e.g. at least one pval to be below threshold

    min_pval_names <- paste0("pval_", min_pval_files)
    min_pval_in_min_pval_files <- apply(SNP_inter[, ..min_pval_names], 1, min) # e.g. need all pvals > 0.5

    SNP_inter <-  SNP_inter[min_pval_in_max_pval_files < max_pval_thres &
                            min_pval_in_min_pval_files > min_pval_thres, ]
    selected_SNP <- SNP_inter$SNP

    if (clump) {
        message("Start clumping using PLINK ...")
        clump_pval_names <- paste0("pval_", clump_pval_files)
        clump_pval <- apply(SNP_inter[, ..clump_pval_names], 1, min)

        SNP_inter$pval <- clump_pval

        tmp2 <- plink_clump(SNP_inter,
                            plink_exe = plink_exe,
                            plink_refdat = plink_refdat,
                            clump_r2 = clump_r2,
                            clump_p1 = max_pval_thres)
        selected_SNP <- tmp2$SNP
    }

    selected_SNP <- as.character(selected_SNP) # selected SNPs

    selected_SNP
}

#' Filter SNPs for different purporses
#'
#' @param sel_files A vector of GWAS summary data file names for selection
#' @param exp_files A vector of GWAS summary data file names for exposures
#' @param out_files A vector of GWAS summary data file names for outcomes
#' @param SNP_inter "Intersection" SNPs that are in all files returned by \code{get_SNP_intersection}.
#' @param utility Purpose of filtering SNPs.
#' @param pval_thres P-value threshold used for filtering SNPs (lower bound for \code{utility=correlation} and upper bound otherwise).
#' @inheritParams plink_clump
#'
#' @return a vector of SNP names.
#'
#' @export
#'
filter_SNPs <- function(sel_files = NULL, exp_files = NULL, out_files = NULL,
                        SNP_inter = NULL,
                        utility = c("inference", "marker_sel", "marker_exp", "correlation"),
                        pval_thres = 5e-8,
                        clump_r2 = switch(utility,
                                          inference = 0.01,
                                          marker_sel = 0.3, # Could discuss with Jingshu
                                          marker_exp = 0.3),
                        plink_exe = NULL,
                        plink_refdat = NULL) {

    utility <- match.arg(utility)
    files  <-  c(sel_files, exp_files, out_files)

    if (utility == "inference") {

        filter_SNPs_basic(files, SNP_inter,
                          max_pval_files = sel_files, max_pval_thres = pval_thres,
                          clump = TRUE, clump_pval_files = sel_files, clump_r2 = clump_r2,
                          plink_exe = plink_exe, plink_refdat = plink_refdat)

    } else if (utility == "marker-exp") {

        filter_SNPs_basic(files, SNP_inter,
                          max_pval_files = exp_files, max_pval_thres = pval_thres,
                          clump = TRUE, clump_pval_files = exp_files, clump_r2 = clump_r2,
                          plink_exe = plink_exe, plink_refdat = plink_refdat)

    } else if (utility == "marker-sel") {

        filter_SNPs_basic(files, SNP_inter,
                          max_pval_files = sel_files, max_pval_thres = pval_thres,
                          clump = TRUE, clump_pval_files = sel_files, clump_r2 = clump_r2,
                          plink_exe = plink_exe, plink_refdat = plink_refdat)

    } else if (utility == "correlation") {

        filter_SNPs_basic(files, SNP_inter,
                          min_pval_files = sel_files, min_pval_thres = pval_thres,
                          clump = FALSE)

    }

}

## #' Main function for data preprocessing
## #'
## #'
## getInput <- function(sel_files, exp_files, out_files, sel_SNPs = NULL,
##                      cor_SNPs = NULL, pval_thres_cor = (1-1e-3), cal_cor = TRUE,
##                      mar_SNPs = NULL,
##                      get_marker_candidates = T, marker_p_source = "exposure",
##                      marker_pval_thres = 1e-5,
##                      clump_r2_formarkers = 0.05,
##                      plink_exe = "./plink",
##                      plink_refdat = "./util/data_maf0.01_rs_ref",
##                      max_pval_thres = 0.01, clump_r2 = 0.001) {

##     num_sel <- length(sel_files)
##     num_exp <- length(exp_files)
##     num_out <- length(out_files)

##     files  <-  c(sel_files, exp_files, out_files) # list of all files to be processed

##     if (num_exp > 1) {
##         if (get.marker.candidates)
##             message("Marker candidates will not be obtained as number of risk factors k > 1")
##         get.marker.candidates <- F
##     }


##     if (is.null(sel_SNPs)){
##         temp <- selectSNPs(files, max_pval_thres, clump_r2, pval_thres_cor, cal_cor, plink_exe,
##                            plink_refdat)

##         if (get_marker_candidates & (marker_p_source == "selection")) {
##             SNP_inter <- temp[[3]] # used later for getting markers
##         }

##         sel_SNPs <- temp[[1]]
##         cor_SNPs <- temp[[2]]

##         selectedSNPs <- sel_SNPs
##         selectedcorSNPs <- cor_SNPs
##         message('Selected SNPs for inference data and estimating correlation, extracting data...')
##     }


##     # get list of SNPs for data and harmonise

##     dat <- extract_data(files, sel_SNPs)
##     dat <- harmonise_data_list(dat, fast = F)  # list of harmonised data sets
##     dat$num_sel <- num_sel
##     dat$num_exp <- num_exp
##     dat$num_out <- num_out

##     data <- get_data(dat)
##     message('Inference data completed.')
##     head(data)

##     if (cal_cor) {
##         cor_dat <- extract_data(files, cor_SNPs)
##         cor_dat <- harmonise_data_list(dat, fast = T)

##         cor_dat <- get_data(dat)
##         corr <- calc_corr(cor_dat)
##         message('Correlation matrix computed.')
##         head(corr)
##     }

##     if (get_marker_candidates & is.null(mar_SNPs)) {
##         if (marker_p_source == "selection") {
##             message('Extracting marker data, source = selection...')
##             tmp <- selectSNPs(files = NULL, SNP_inter = SNP_inter, max_pval_thres = 1, clump_r2 = clump_r2_formarkers, cal_cor = F,
##             pval_thres_cor = 0.5, plink_exe = './plink', plink_refdat)

##             mar_dat <- extract_data(files, mar_SNPs)
##             mar_dat <- harmonise_data_list(dat, fast = T)

##             marker_data <- get_data(mar_dat)
##             marker_data$SNP <- marker_data$SNP[pval < marker_pval_thres]

##             message('Extracted marker data.')
##             head(mar_dat)
##         } else {
##             # If marker_p_source is exposure, only take intersection of exposure files and harmonise them
##             message('Extracting Marker data, source = exposure...')
##             exp_files <- files[1:num_exp]
##             mar_SNP_inter <- get_SNP_intersection(exp_files, pval_thres = marker_pval_thres/num_sel) # bonferroni correction

##             mar_SNPs <- selectSNPs(files = NULL, SNP_inter = mar_SNP_inter, max_pval_thres = 1, clump_r2 = clump_r2_formarkers, cal_cor = F,
##                               pval_thres_cor = 0.5, plink_exe = './plink', plink_refdat)[[1]]

##             mar_dat <- extract_data(exp_files, mar_SNPs)
##             mar_dat <- harmonise_data_list(exp_files, fast = T)


##             marker_data <- get_data(mar_dat)
##             message('Marker data extracted')
##         }
##     }
## }

#' Format GWAS data
#'
#' Group the same columns in a list of GWAS data tables together.
#'
#' @param dat_list A list of GWAS data tables (typically returned by \code{harmonise_data_list}).
#' @inheritParams filter_SNPs
#'
#' @return A list with the following variables:
#' \itemize{
#' \item \code{meta_data} (a data frame);
#' \item \code{beta_sel}, \code{se_sel}, \code{pval_sel} (returned if \code{sel_files} is not NULL);
#' \item \code{beta_exp}, \code{se_exp}, \code{pval_exp} (returned if \code{exp_files} is not NULL);
#' \item \code{beta_out}, \code{se_out}, \code{pval_out} (returned if \code{out_files} is not NULL);
#' \item \code{beta}, \code{se}, \code{pval} (returned if all files are NULL);
#' }
#'
#' @export
#'
format_data <- function(dat_list,
                        sel_files = NULL,
                        exp_files = NULL,
                        out_files = NULL) {

    ## Check the data tables in the list have the same SNP, effect_allele, and reference_allele
    check_same <- function(x) {
        prod(sapply(x, all.equal, x[[1]])) == 1
    }
    stopifnot(check_same(lapply(dat_list, function(x) x$SNP)))
    stopifnot(check_same(lapply(dat_list, function(x) x$effect_allele)))
    stopifnot(check_same(lapply(dat_list, function(x) x$other_allele)))

    meta_data <- dat_list[[1]][, c("SNP", "effect_allele", "other_allele")]

    beta <- lapply(dat_list, function(x) x$beta)
    beta <- do.call("cbind", beta)

    se <- lapply(dat_list, function(x) x$se)
    se <- do.call("cbind", se)

    pval <- beta
    for (i in 1:length(dat_list)) {
        if ("pval" %in% colnames(dat_list[i])) {
            pval[, i] <- dat_list[i]$pval
        } else {
            pval[, i] <- pnorm(- abs(beta[, i] / se[, i])) * 2
        }
    }

    colnames(beta) <- names(dat_list)
    colnames(se) <- names(dat_list)
    colnames(pval) <- names(dat_list)

    output <- list(meta_data = meta_data,
                   sel_files = sel_files,
                   exp_files = exp_files,
                   out_files = out_files)

    if (!is.null(sel_files)) {
        output$beta_sel <- beta[, sel_files]
        output$se_sel <- se[, sel_files]
        output$pval_sel <- pval[, sel_files]
    }

    if (!is.null(exp_files)) {
        output$beta_exp <- beta[, exp_files]
        output$se_exp <- se[, exp_files]
        output$pval_exp <- pval[, exp_files]
    }

    if (!is.null(out_files)) {
        output$beta_out <- beta[, out_files]
        output$se_out <- se[, out_files]
        output$pval_out <- pval[, out_files]
    }

    if (is.null(sel_files) && is.null(exp_files) && is.null(out_files)) {
        output$beta <- beta
        output$se <- se
        output$pval <- pval
    }

    output
}

#' Estimate the sample correlation between GWAS datasets
#'
#' @param data A list with following elements: SNP, beta_exp, se_exp, beta_out, se_out (typically returned by \code{format_data}).
#' @param cor_SNPs A vector of SNPs used to compute the correlation.
#'
#' @return Estimated sample correlation matrix
#'
#' @export
#'
calc_corr <- function(data, cor_SNPs){

    z_values <- cbind(data$beta_exp[SNP %in% cor_SNPs, ]/data$se_exp[SNP %in% cor_SNPs, ],
                      data$beta_out[SNP %in% cor_SNPs, ]/data$se_out[SNP %in% cor_SNPs, ])

    colnames(z_values) <- c(data$exp_files, data$out_files)

    z_values <- as_matrix(z_values)
    z_values <- z_values[rowSums(is.na(z_values)) == 0,  , drop = F]
    covv <- t(z_values) %*% z_values / nrow(z_values)
    varr <- colMeans(z_values^2, na.rm = T)
    corr <- t(covv / sqrt(varr))/sqrt(varr)
    corr
}
