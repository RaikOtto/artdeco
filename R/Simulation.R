#' Simulate expression profiles of specific cell-types.
#'
#' \code{simulateCellTypes} simulates \emph{in silico} expression profiles
#' of specific cell-types using a negative binomial distribution.
#' Simulation is based on biological data and marker genes
#' of these cell-types.
#'
#' @param referenceCellTypes Matrix of single-cell expression values
#'  of the cell-type to be simulated. Used to estimate parameters for
#'  negative binomial distribution. Samples as columns and genes as rows.
#'  Counts have to be normalized beforehand.
#' @param markerGenes Character vector containing marker genes which
#'  characterize the cell-type to be simulated. Gene identifiers need
#'  to be consistend with \code{referenceCellTypes}
#' @param numSamples A integer specifing the number of samples to be simulated.
#' @param seed A integer specifing seed for simulation. Default is \code{NULL}.
#'  Random seed will be generated if no seed is provided.
#' @param verbose Logical, indicating whether status updates will be printed.
#'  Default is \code{TRUE}.
#' @importFrom msir loess.sd
#' @importFrom stats rnbinom
#' @return Matrix with simulated count values. Samples as columns and genes
#'  as rows.
#' @usage
#' simulateCellTypes(
#'     referenceCellTypes,
#'     markerGenes,
#'     numSamples,
#'     seed,
#'     verbose
#' )
#' @examples
#' markerGenes = 1:10
#' referenceCellTypes = matrix(
#'     rbinom(20*20, c(100,500,1000),prob = .8),
#'     ncol = 20,
#'     dimnames = list(
#'         1:20,1:20
#'     )
#' )
#' simulateCellTypes(
#'     referenceCellTypes = referenceCellTypes,
#'     markerGenes = markerGenes,
#'     numSamples = 20,
#'     seed = 123,
#'     verbose = TRUE
#' )
#' @export
simulateCellTypes = function(
    referenceCellTypes,
    markerGenes,
    numSamples,
    seed = NULL,
    verbose = TRUE
){

    # Split reference into marker genes and non marker genes
    if(verbose){message("Processing input files ...")}
    markerCounts = referenceCellTypes[
        rownames(referenceCellTypes) %in% markerGenes,]
    
    if(nrow(markerCounts) == 0){
        stop("Non of the marker genes present in reference cell-types!")
    }
    otherCounts = referenceCellTypes[
        !rownames(referenceCellTypes) %in% markerGenes,]

    # Estimate parameters for marker genes and other genes separately
    if(verbose){message("Estimating parameters ...")}

        paramMarkers = estimateParameters(
            markerCounts,
            simMarkers = TRUE
        )
        paramOthers = estimateParameters(
            otherCounts,
            simMarkers = FALSE
        )

    # Simulate counts
    if(verbose){message("Simulating counts ...")}
    if(is.null(seed)) {
        seed = sample(1:1000000, size = 1)
    }
    #set.seed(seed)

    sim_Counts_Markers = simulateCounts(
        simParam = paramMarkers,
        nSamples = numSamples,
        nGenes = nrow(markerCounts),
        simSeed = seed,
        simMarkers = TRUE)

    sim_Counts_Others = simulateCounts(
        simParam = paramOthers,
        nSamples = numSamples,
        nGenes = nrow(otherCounts),
        simSeed = seed,
        simMarkers = FALSE
    )

    # Rename genes of simulated data, join and create sample IDs
    if(verbose){message("Preparing output ...")}
    rownames(sim_Counts_Markers) = rownames(markerCounts)
    rownames(sim_Counts_Others) = rownames(otherCounts)

    res = rbind(sim_Counts_Markers, sim_Counts_Others)
    colnames(res) = paste0("simu_", 1:numSamples)

    if(verbose){message("Done!")}
    return(res)
}

#' Estimate parameters for negative binomial distribution
#'
#' \code{estimateParameters} estimates the parameters for the negative binomial
#' distribution from the provided reference expression profiles. This function
#' is used by \code{\link{simulateCellTypes}}.
#'
#' @param countData Matrix with the expression values of the reference
#'  cell-types. Parameters are estimated based on this count data.
#' @param simMarkers Logical, indicating whether genes to be simulated are
#'  markers or not
#' @return Returns a list with the estimated parameters.
estimateParameters = function(
    countData,
    simMarkers = TRUE
){
    # Set parameters
    sigma = 1.96

    # Kick out empty samples and keep only expressed genes
    totalS = ncol(countData)
    totalG = nrow(countData)
    fullS = colSums(countData, na.rm = TRUE) > 0
    detectG = rowMeans(countData, na.rm = TRUE) > 0
    countData = countData[detectG, fullS]

    nsamples = dim(countData)[2]
    counts0 = countData == 0
    nn0 = rowSums(!counts0)

    # the negative binomial
    mu = rowSums(countData) / ncol(countData)
    s2 = rowSums((countData - mu) ^ 2) / ncol(countData)
    size = mu ^ 2 / (s2 - mu + 1e-04)
    size = ifelse(size > 0, size, NA)
    p0 = (nsamples - nn0) / nsamples
    mu = mu[!is.na(size)]
    p0 = p0[!is.na(size)]
    remove = rownames(countData)[is.na(size)]
    detectG[names(detectG) %in% remove] = FALSE
    size = size[!is.na(size)]
    phi.g = 1 / size
    phi.c = mean(phi.g)
    ldisp = log2(phi.g)
    lsize = log2(size)
    lmu = log2(mu + 1)

    estG = length(mu)
    estS = length(ncol(countData))

    # meansizefit
    meansizefit = loess.sd(lsize ~ lmu, nsigma = sigma)

    # meandispfit
    meandispfit = loess.sd(ldisp ~ lmu, nsigma = sigma)

    # return object
    paramData = list(means = mu,
                     dispersion = phi.g,
                     common.dispersion = phi.c,
                     size = size,
                     p0 = p0,
                     meansizefit = meansizefit,
                     meandispfit = meandispfit,
                     estS = estS,
                     estG = estG,
                     totalS = totalS,
                     totalG = totalG,
                     detectG = detectG,
                     sigma = sigma)

    return(paramData)
}

#' Simulate count data based on negative binomial distribution
#'
#' \code{simulateCounts} simulates the expression values using
#' a negative binomial distribution with parameters estimated by
#' \code{\link{estimateParameters}}. This function is used by
#' \code{\link{simulateCellTypes}}.
#' @param simParam List of parameters estimated by
#'  \code{\link{estimateParameters}}.
#' @param nSamples A integer specifing the number of samples to be simulated.
#' @param nGenes A integer specifing the number of genes to be simulated.
#' @param simSeed A integer specifing seed for simulation
#' @param simMarkers Logical, indicating whether genes to be simulated are
#'  markers or not
#' @return Matrix with simulated count values.
simulateCounts = function(
    simParam,
    nSamples,
    nGenes,
    simSeed = NULL,
    simMarkers = TRUE
){
    if(simMarkers){
        lfcs = as.matrix(rep(0, nGenes))
    } else {
        lfcs = as.matrix(rep(0, nGenes))
    }

    # define NB params
    mu = simParam$means
    meansizefit = simParam$meansizefit

    # For markers use observed mean parameters, for other genes sample
    if(simMarkers) {
        present = simParam$detectG
        if(sum(present) < nGenes){
            warning("Detected one or more marker genes with no expression
                    in reference samples!")
            true.means = rep(0, nGenes)
            true.means[present] = mu
        } else {
            true.means = mu
        }
    } else {
        index = sample(1:length(mu), size = nGenes, replace = TRUE)
        true.means = mu[index]
    }

    # estimate size parameter associated with true mean values
    lmu = log2(true.means + 1)
    predsize.mean = approx(meansizefit$x, meansizefit$y, xout = lmu, rule = 2)$y
    predsize.sd = approx(meansizefit$x, meansizefit$sd, xout = lmu, rule = 2)$y
    sizevec = rnorm(n = length(lmu), mean = predsize.mean, sd = predsize.sd)

    # size factor
    #all.facs = rep(1, nSamples)
    all.facs = sample(seq(1, 2.5, by = 0.1), size = nSamples, replace = TRUE)

    # effective means
    effective.means = outer(true.means, all.facs, "*")

    mod = as.matrix(rep(1, nSamples))

    # make mean expression with beta coefficients added as defined
    # by model matrix
    mumat = log2(effective.means + 1) + lfcs %*% t(mod)
    mumat[mumat < 0] = min(log2(effective.means + 1))

    # result count matrix
    counts = matrix(
        rnbinom(nSamples * nGenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec),
        ncol = nSamples,
        nrow = nGenes,
        dimnames = list(paste0(rownames(mumat),"_", seq_len(nGenes)),
                        NULL))

    return(counts)
}

#' Simulate expression data which does not follow a cell-type specific profile.
#'
#' \code{simulateNegativeControls} simulates random expression profiles
#' following a negative binomial distribution. Parameters of the negative
#' binomial are randomly sampled from a normal distribution.
#' Can be used as negative controls for benchmarking purposes.
#' @param nGenes A integer specifing the number of genes to be simulated.
#' @param numSamples A integer specifing the number of samples to be simulated.
#' @param normMean A integer specifing the mean parameter of the normal distribution
#'  which is used to generate mean expression values for the simulation.
#' @param normSD A integer specifing the standard deviation parameter of the
#'  normal distribution which is used to to generate mean expression values for
#'  the simulation.
#' @param seed A integer specifing seed for simulation. Default is \code{NULL}.
#'  Random seed will be generated if no seed is provided.
#' @param verbose Logical, indicating whether status updates will be printed.
#'  Default is \code{TRUE}.
#' @importFrom stats rnbinom
#' @return Matrix with simulated count values. Samples as columns and genes
#'  as rows.
#' @usage
#' simulateNegativeControls(
#'     nGenes,
#'     numSamples,
#'     normMean,
#'     normSD,
#'     seed,
#'     verbose
#' )
#' @examples
#' simulateNegativeControls(
#'     nGenes = 20000,
#'     numSamples = 10,
#'     normMean = 50,
#'     normSD = 500,
#'     seed = 123,
#'     verbose = TRUE
#' )
#' @export
simulateNegativeControls = function(
    nGenes,
    numSamples,
    normMean = 50,
    normSD = 500,
    seed = NULL,
    verbose = TRUE
){
    if(is.null(seed)) {
        seed = sample(1:1000000, size = 1)
    }
    #set.seed(seed)

    # Simulate counts based on randomly sampled mean and size values
    # for negative binomial distribution
    means = rnorm(nGenes, mean = normMean, sd = normSD)
    means[means < 0] = 0
    all.facs = sample(seq(1, 3, by = 0.1), size = numSamples, replace = TRUE)
    effective.means = outer(means, all.facs, "*")
    mumat = log2(effective.means + 1)
    mumat[mumat < 0] = min(log2(effective.means + 1))
    sizevec = rnorm(nGenes, mean = 0, sd = 4)

    counts = matrix(
        rnbinom(numSamples * nGenes, mu = 2 ^ mumat - 1, size = 2 ^ sizevec),
        ncol = numSamples,
        nrow = nGenes
    )
    colnames(counts) = paste0("Negative_Control_", c(1:numSamples))

    return(counts)
}
