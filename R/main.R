#' Mutational signature fitting with mmsig
#'
#' @param muts.input mutation input data frame as specified under input.format
#' @param sig.input mutational signature reference with mutational classes as rows and signature exposures as columns: Substitution.Type, Trinucleotide, signature.1, ... signature.n
#' @param input.format sigs.input: pre-computed mutation count matrix in mut.to.sigs.input() output format (samples as rows, 96 mutational classes as columns). classes: samples as columns and the 96 mutational classes as rows
#' @param sample.sigt.profs NULL = use all signatures provided in the reference. Optionally provide list with signatures to consider for each sample
#' @param bootstrap TRUE/FALSE for whether bootstrapping is to be performed
#' @param iterations number of bootstrapping iterations to perform (only if bootstrap == TRUE)
#' @param strandbias TRUE/FALSE for whether transcriptional strand bias should be tested for (only for vcf-like input format)
#' @param cos_sim_threshold cosine similarity threshold below which signatures are removed from the final profile
#' @param force_include vector with the names of signatures to always keep in the final profile of every sample
#' @param dbg FALSE = silent; TRUE = verbose
#' @importFrom dplyr left_join
#'
#' @return mutational signature fitting results for all samples
#' @export
#'
mm_fit_signatures = function(muts.input,
                             sig.input,
                             input.format = "sigs.input",
                             sample.sigt.profs=NULL,
                             bootstrap=FALSE,
                             iterations=1000,
                             strandbias=FALSE,
                             cos_sim_threshold=0.01,
                             force_include=c("SBS1", "SBS5"),
                             dbg=FALSE) {

  options(scipen = 999)

  #####################################################
  ########        Read/assign input data       ########
  #####################################################

  # Mutational signature reference
  consigts.defn <- sig.input

  # Input mutation data
  if(input.format == "sigs.input") {

    # Validate that input looks like mut.to.sigs.input() output:
    # rows = samples, cols = 96 trinucleotide mutation classes
    if(!is.data.frame(muts.input) && !is.matrix(muts.input)){
      stop("ERROR: muts.input must be a data frame or matrix in mut.to.sigs.input() output format
           (samples as rows, 96 mutational classes as columns)")
    }

    if(ncol(muts.input) < 96){
      stop("ERROR: muts.input must have at least 96 columns corresponding to trinucleotide mutation classes.
           Please provide input in mut.to.sigs.input() output format.")
    }

    # Transpose: convert from (samples x classes) to (classes x samples)
    # to match the internal format expected downstream
    samples.muts <- as.data.frame(t(muts.input))
    samples.muts <- samples.muts[1:96, , drop=FALSE]   # keep only the 96 mutation class rows

  } else if(input.format == "classes"){

    # classes format: rows = 96 mutational classes, columns = samples
    samples.muts <- muts.input
    samples.muts <- samples.muts[names(samples.muts) != "Total"]        # remove totals column if present
    samples.muts <- samples.muts[1:96, , drop=FALSE]                    # keep only the 96 mutation class rows

  } else{
    stop("ERROR: Invalid input format. Please provide:
         'sigs.input' - mut.to.sigs.input() output format (samples as rows, 96 classes as columns), or
         'classes'    - 96 mutational classes as rows, samples as columns")
  }

  # Process signature reference
  tcons.defn <- t(consigts.defn[,3:ncol(consigts.defn)])
  mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
  consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric)
  rownames(consigts.defn) <- mutlist
  ref_signatures <- colnames(consigts.defn)  # names of signatures in reference

  # Process sample mutational profiles
  samples <- colnames(samples.muts)

  # Assign which signatures to fit to each sample
  if(is.list(sample.sigt.profs)){
    spit(dbg, "using mm signature profiles from input argument")
    sigt.profs <- sample.sigt.profs
  } else {
    spit(dbg, "defaulting to use all signatures in the provided reference")
    mm.sigts <- ref_signatures

    sigt.profs <- list()
    for (i in 1:length(samples)) {
      sigt.profs[[samples[i]]] <- mm.sigts
    }
  }

  # Define output variable
  output <- list()

  #####################################################
  ########       Fit mutational signatures     ########
  #####################################################

  sigfit <- fit_signatures(samples.muts=samples.muts,
                           consigts.defn=consigts.defn,
                           sigt.profs=sigt.profs,
                           cos_sim_threshold=cos_sim_threshold,
                           force_include=force_include,
                           dbg=dbg)

  output$estimate <- sigfit

  #####################################################
  ########             Bootstrapping           ########
  #####################################################

  if(bootstrap){

    # Point estimates
    sig_est <- sigfit
    sig_est$sample <- row.names(sig_est)
    sig_est <- melt(sig_est[names(sig_est) != "mutations"],
                    id.vars = "sample",
                    variable.name = 'signature',
                    value.name = 'estimate')

    sig_est$signature <- as.character(sig_est$signature)

    sigboot <- bootstrap_fit_signatures(samples.muts = samples.muts,
                                        consigts.defn = consigts.defn,
                                        sigt.profs = sigt.profs,
                                        iterations = iterations,
                                        cos_sim_threshold = cos_sim_threshold,
                                        force_include = force_include)

    sigboot <- left_join(sigboot, sig_est, by = c("sample", "signature"))

    output$bootstrap <- sigboot
  }

  #####################################################
  ########              Strand Bias            ########
  #####################################################

  if(strandbias){
    warning("Transcriptional strand bias cannot be estimated from pre-computed mutation class input.
             Please provide vcf-like input to a strand bias function separately.")
  }

  output$mutmatrix <- samples.muts

  return(output)
}
