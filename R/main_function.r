#' Categorical Clustering - Unified Interface
#'
#' Run any of the supported categorical clustering methods and evaluate
#' the partition using the Quality Index (QI).
#'
#' @param data A data.frame or matrix with categorical variables
#' @param method Character string specifying the method. One of:
#'   \itemize{
#'     \item \code{"aac_sequential"} - AAC algorithm (pure R)
#'     \item \code{"aac_parallel"} - AAC with mclapply (macOS native)
#'     \item \code{"aac_rcpp"} - AAC with Rcpp + mclapply (fastest)
#'     \item \code{"hierarchical"} - Ward hierarchical with Gower
#'     \item \code{"diana"} - DIANA divisive clustering
#'     \item \code{"dbscan"} - DBSCAN density-based
#'     \item \code{"hdbscan"} - HDBSCAN hierarchical density-based
#'     \item \code{"spectral"} - Spectral clustering
#'     \item \code{"pop"} - Partitioning Around PAMoids
#'   }
#' @param k_range Range of k values to evaluate (method-dependent)
#' @param n_cores Number of cores for parallel methods
#' @param seed Random seed for reproducibility
#' @param verbose Logical, print progress
#' @param ... Additional arguments passed to the specific method
#'
#' @return A list with components: method, k_opt, QI, groups, time
#'
#' @examples
#' \dontrun{
#' data(Felinos)
#' res <- categorical_clustering(Felinos, method = "aac_rcpp")
#' print(res)
#' }
#'
#' @export
categorical_clustering <- function(data, method = "aac_rcpp",
                                   k_range = NULL, n_cores = NULL,
                                   seed = 42, verbose = FALSE, ...) {
  
  set.seed(seed)
  method <- tolower(method)
  
  valid_methods <- c("aac_sequential", "aac_parallel", "aac_rcpp",
                     "hierarchical", "diana", "dbscan", "hdbscan",
                     "spectral", "pop")
  
  if (!method %in% valid_methods) {
    stop("Invalid method. Choose one of: ",
         paste(valid_methods, collapse = ", "))
  }
  
  t0 <- proc.time()
  
  result <- switch(method,
    "aac_sequential" = aac_sequential(data, k_range = k_range, verbose = verbose),
    "aac_parallel"   = aac_parallel(data, k_range = k_range, n_cores = n_cores, verbose = verbose),
    "aac_rcpp"       = aac_rcpp_parallel(data, k_range = k_range, n_cores = n_cores, verbose = verbose),
    "hierarchical"   = hierarchical_clustering(data, k_range = k_range, ...),
    "diana"          = diana_clustering(data, k_range = k_range),
    "dbscan"         = dbscan_clustering(data, ...),
    "hdbscan"        = hdbscan_clustering(data, seed = seed, ...),
    "spectral"       = spectral_clustering(data, k_range = k_range, seed = seed),
    "pop"            = pop_clustering(data, ...)
  )
  
  elapsed <- (proc.time() - t0)[["elapsed"]]
  result$time <- elapsed
  result$data <- data
  result$seed <- seed
  class(result) <- c("AACclustMac_result", "list")
  return(result)
}

#' Search for Optimal Number of Clusters
#'
#' Runs all available methods and returns the partition with the highest QI.
#'
#' @param data A data.frame or matrix with categorical variables
#' @param k_range Range of k values to evaluate
#' @param n_cores Number of cores for parallel methods
#' @return List with the best method and its results
#' @export
Find_k_optimo <- function(data, k_range = NULL, n_cores = NULL) {
  methods <- c("aac_rcpp", "hierarchical", "diana", "dbscan",
               "hdbscan", "spectral", "pop")
  
  results <- list()
  for (m in methods) {
    cat(sprintf(">>> Running method: %s\n", m))
    res <- tryCatch({
      categorical_clustering(data, method = m, k_range = k_range, n_cores = n_cores)
    }, error = function(e) {
      warning(sprintf("Method %s failed: %s", m, e$message))
      NULL
    })
    if (!is.null(res)) {
      results[[m]] <- res
    }
  }
  
  # Find best QI
  qis <- sapply(results, function(x) x$QI)
  best_method <- names(which.max(qis))
  
  list(
    best_method = best_method,
    best_result = results[[best_method]],
    all_results = results,
    comparison = data.frame(
      Method = names(qis),
      k = sapply(results, function(x) x$k_opt),
      QI = qis,
      Time = sapply(results, function(x) x$time),
      stringsAsFactors = FALSE
    )
  )
}

#' Compare Multiple Clustering Methods
#'
#' @param data A data.frame or matrix
#' @param methods Character vector of method names
#' @param k_range Range of k values
#' @param n_cores Number of cores for parallel methods
#' @return A data.frame with comparison results
#' @export
compare_methods <- function(data,
                            methods = c("aac_rcpp", "hierarchical", "diana",
                                        "dbscan", "hdbscan", "spectral", "pop"),
                            k_range = NULL, n_cores = NULL) {
  
  results <- list()
  for (m in methods) {
    cat(sprintf("\n>>> Running method: %s\n", m))
    res <- tryCatch({
      categorical_clustering(data, method = m, k_range = k_range, n_cores = n_cores)
    }, error = function(e) {
      warning(sprintf("Method %s failed: %s", m, e$message))
      NULL
    })
    if (!is.null(res)) {
      results[[m]] <- data.frame(
        Method = res$method,
        k = res$k_opt,
        Time = res$time,
        QI = res$QI,
        stringsAsFactors = FALSE
      )
    }
  }
  
  df <- do.call(rbind, results)
  df <- df[order(-df$QI), ]
  rownames(df) <- NULL
  class(df) <- c("AACclustMac_comparison", "data.frame")
  return(df)
}

#' Print method for AACclustMac_result
#' @param x An AACclustMac_result object
#' @param ... Additional arguments
#' @export
print.AACclustMac_result <- function(x, ...) {
  cat(sprintf("AACclustMac Result\n"))
  cat(sprintf("  Method : %s\n", x$method))
  cat(sprintf("  k*     : %s\n", x$k_opt))
  cat(sprintf("  QI     : %.4f\n", x$QI))
  cat(sprintf("  Time   : %.3f s\n", x$time))
  cat(sprintf("  Clusters: %d\n", length(x$groups)))
  invisible(x)
}

#' Print method for AACclustMac_comparison
#' @param x An AACclustMac_comparison object
#' @param ... Additional arguments
#' @export
print.AACclustMac_comparison <- function(x, ...) {
  cat("\n=== Method Comparison ===\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}