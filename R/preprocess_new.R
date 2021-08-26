#' Loads an R data file, and returns it
#'
#' @param file Name of the R data file (.rda or .RData).
#' @return First object in the R data file.
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
    dat[!duplicated(dat$SNP), c(necessary_columns, "pval")]
}

#' Get intersection of the SNPs in several GWAS summary datasets
#'
#' @param files A vector of file names
#' @param message Logical indicator for message printing
#'
#' @importFrom data.table merge
#'
#' @return A \code{data.table} object with the SNP name and p-values in all summary datasets.
#'
get_snp_intersection <- function(files, message = TRUE) {

    message("Finding intersection of SNPs from all datasets...")
    dat <- read_gwas_summary(files[1], message)
    dat <- dat[, c("SNP", "pval")]

    if (length(files) > 1) {
        for (file in files[-1]) {

            dat_new <- read_gwas_summary(file, message)

            dat <- merge(dat, dat_new[, c("SNP", "pval")],
                         by = "SNP",
                         suffixes = c("", paste0("_", file)))
        }
    }

    names(dat)[names(dat) == "pval"] <- paste0("pval_", files[1])
    dat

}

#' Extract selected SNPs from data
#'
#' @param files A vector of file names
#' @param SNP_list If not null, only returns SNPs in this list; otherwise returns the original dataset.
#' @param message Logical indicator for message printing
#'
#' @return A list of data tables, each corresponding to one file.
#'
extract_data <- function(files, SNP_list = NULL, message = TRUE) {

    dat <- list()
    for (i in 1:length(files)) {
        dat[[i]] <- read_gwas_summary(files[i], message = message)
        if (!is.null(SNP_list)) {
            dat[[i]] <- dat[[i]][SNP %in% SNP_list]
        }
    }
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
harmonise_data_list <- function(dat, fast = T) {
    message("Harmonising a list of data...")

    dat_new <- list()
    SNP_list <- dat[[1]]$SNP
    dat_new[[1]] <- dat[[1]]

    ## Reference dataset
    data1 <- formatData(dat[[1]], "exposure")

    stopifnot(length(dat) > 1)

    for (i in 2:length(dat)) {

        if (fast == T) {
            data2 <- suppressMessages(harmonise_data_fast(data1, data2))
        } else {
            ## Format and harmonise with the reference dataset
            data2 <- formatData(dat[[i]], "outcome")
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

selectSNPs <- function(files = NULL,
                       snp_inter = NULL,
                       max_pval_files = files,
                       max_pval_thres = 1,
                       min_pval_files = files,
                       min_pval_thres = 0,
                       clump = FALSE,
                       clump_r2 = 0.001,
                       plink_exe = './plink',
                       plink_refdat) {

    message("Selecting SNPs for inference and correlation estimation...")

    if (is.null(snp_inter)) {
        snp_inter <- get_snp_intersection(files)
    }

    ## TODO: Select SNPs based on p-values


    message("Start clumping using PLINK ...")
    tmp2 <- plink_clump(snp_inter, plink_exe, refdat = plink_refdat, clump_r2 = clump_r2,
                        clump_p1 = max_pval_thres)

    sel_SNPs <- as.character(tmp2$SNP) # selected SNPs
    sel_SNPs_pvals <- tmp2$pval # keep pvals for marker SNPs pval threshold


    list(sel_SNPs, cor_SNPs, snp_inter, sel_SNPs_pvals)
}


getInput <- function(sel_files, exp_files, out_files, sel_SNPs = NULL,
                     cor_SNPs = NULL, pval_thres_cor = (1-1e-3), cal_cor = TRUE,
                     mar_SNPs = NULL,
                     get_marker_candidates = T, marker_p_source = "exposure",
                     marker_pval_thres = 1e-5,
                     clump_r2_formarkers = 0.05,
                     plink_exe = "./plink",
                     plink_refdat = "./util/data_maf0.01_rs_ref",
                     max_pval_thres = 0.01, clump_r2 = 0.001) {

    num_sel <- length(sel_files)
    num_exp <- length(exp_files)
    num_out <- length(out_files)

    files = c(sel_files, exp_files, out_files) # list of all files to be processed

    if (num_exp > 1) {
        if (get.marker.candidates)
            message("Marker candidates will not be obtained as number of risk factors k > 1")
        get.marker.candidates <- F
    }


    if (is.null(sel_SNPs)){
        temp <- selectSNPs(files, max_pval_thres, clump_r2, pval_thres_cor, cal_cor, plink_exe,
                           plink_refdat)

        if (get_marker_candidates & (marker_p_source == "selection")) {
            snp_inter <- temp[[3]] # used later for getting markers
        }

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
    message('Inference data completed.')
    head(data)

    if (cal_cor) {
        cor_dat <- extract_data(files, cor_SNPs)
        cor_dat <- harmonise_data_list(dat, fast = T)

        cor_dat <- get_data(dat)
        corr <- calc_corr(cor_dat)
        message('Correlation matrix computed.')
        head(corr)
    }

    if (get_marker_candidates & is.null(mar_SNPs)) {
        if (marker_p_source == "selection") {
            message('Extracting marker data, source = selection...')
            tmp <- selectSNPs(files = NULL, snp_inter = snp_inter, max_pval_thres = 1, clump_r2 = clump_r2_formarkers, cal_cor = F,
            pval_thres_cor = 0.5, plink_exe = './plink', plink_refdat)

            mar_dat <- extract_data(files, mar_SNPs)
            mar_dat <- harmonise_data_list(dat, fast = T)

            marker_data <- get_data(mar_dat)
            marker_data$SNP <- marker_data$SNP[pval < marker_pval_thres]

            message('Extracted marker data.')
            head(mar_dat)
        } else {
            # If marker_p_source is exposure, only take intersection of exposure files and harmonise them
            message('Extracting Marker data, source = exposure...')
            exp_files <- files[1:num_exp]
            mar_snp_inter <- get_snp_intersection(exp_files, pval_thres = marker_pval_thres/num_sel) # bonferroni correction

            mar_SNPs <- selectSNPs(files = NULL, snp_inter = mar_snp_inter, max_pval_thres = 1, clump_r2 = clump_r2_formarkers, cal_cor = F,
                              pval_thres_cor = 0.5, plink_exe = './plink', plink_refdat)[[1]]

            mar_dat <- extract_data(exp_files, mar_SNPs)
            mar_dat <- harmonise_data_list(exp_files, fast = T)


            marker_data <- get_data(mar_dat)
            message('Marker data extracted')
        }


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
calc_corr <- function(cor_data, cor_SNPs){
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
