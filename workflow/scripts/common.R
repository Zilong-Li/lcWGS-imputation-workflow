library(data.table)


modified_calculate_pse <- function(
    test,
    truth,
    which_snps,
    seed = NULL) {
  ## for testing
  if (is.null(seed) == FALSE) {
    set.seed(seed)
  }
  truth <- truth[which_snps, ]
  test <- test[which_snps, ]
  ## rownames(test) <- 1:nrow(test)
  which_sites <-
    rowSums(truth == 0 | truth == 1, na.rm = TRUE) == 2 &
      rowSums(truth, na.rm = TRUE) == 1 &
      rowSums(is.na(truth)) == 0
  truth <- truth[which_sites, ]
  test <- test[which_sites, ]
  snps <- which_snps[which_sites,]
  if (nrow(test) == 0) {
    return(NA)
  }
  ## as these sites, discrepency as well
  disc <- sum(rowSums(test) != 1)
  ## round test data for now. choose hets at random
  ## for homs
  testO <- test
  truthO <- truth
  for (i_option in 1:2) {
    test <- testO
    truth <- truthO
    ## specifically remove from consideration double phase switch errors
    w <- rowSums(test) == 1
    w2 <- which(diff(abs(test[w, 1] - truth[w, 1])) != 0)
    to_remove <- NULL
    if (length(w2) > 0) {
      w3 <- which(diff(w2) == 1)
      if (length(w3) > 0) {
        for (a in w3) {
          c <- w2[c(a, a + 1)]
          to_remove <- c(to_remove, which(w)[c])
        }
      }
    }
    ## double pse are two consecutive
    if (length(to_remove) > 0) {
      test <- test[-to_remove, ]
      truth <- truth[-to_remove, ]
      snps <- snps[-to_remove]
    }
    ##
    if (i_option == 1) {
      ## only consider non-discrepent sites
      ## chose best start
      if (test[1, 1] != truth[1, 1]) {
        test <- test[, c(2, 1)]
      }
      ## calculate number of differences
      w <- rowSums(test) == 1
      snps <- snps[w]
      if (sum(w) == 0) {
        print("Test has no hets! possibly an error or homo over region, possibly no record dosages turned on in impute_all")
        switches1 <- cbind(i1 = NA, i2 = NA, l1 = NA, l2 = NA)
        phase_errors_def1 <- 0
        phase_sites_def1 <- 0
      } else {
        y <- diff(abs(test[w, 1] - truth[w, 1])) != 0
        phase_errors_def1 <- sum(y)
        phase_sites_def1 <- sum(w) - 1
        ## we need rownames(test) <- 1:nrow(test) beforehand
        ## s <- as.integer(rownames(test[w, , drop = FALSE][c(as.logical(y), FALSE), , drop = FALSE]))
        ## e <- as.integer(rownames(test[w, , drop = FALSE][c(FALSE, as.logical(y)), , drop = FALSE]))
        ## switches1 <- cbind(i1 = s, i2 = e)
      }
    }
    if (i_option == 2) {
      choose_at_random <- which(rowSums(test) != 1)
      if (length(choose_at_random) > 0) {
        test[choose_at_random, ] <- 0
        r <- sample(
          c(1, 2),
          length(choose_at_random),
          replace = TRUE
        )
        test[cbind(choose_at_random, r)] <- 1
      }
      ## chose best start
      if (test[1, 1] != truth[1, 1]) {
        test <- test[, c(2, 1)]
      }
      ## calculate number of differences
      phase_errors_def2 <- sum(diff(abs(test[, 1] - truth[, 1])) != 0)
      phase_sites_def2 <- nrow(test) - 1
    }
  }
  ##
  return(
    list(
      values = c(
        het_sites = snps,
        phase_errors_def1 = phase_errors_def1,
        phase_sites_def1 = phase_sites_def1,
        phase_errors_def2 = phase_errors_def2,
        phase_sites_def2 = phase_sites_def2,
        disc_errors = disc,
        dist_n = nrow(test)
      ),
      switches1 = NULL
    )
  )
}


acc_phasing_single_matrix <- function(test,truth, id) {
  n <- ncol(truth)/2
  a <- lapply(1:n,function(i) {
    t1 <- as.matrix(test[id,1:2+(i-1)*2])
    t0 <- as.matrix(truth[id,1:2+(i-1)*2])
    values <- modified_calculate_pse(t1, t0, id)$values
    sites <- values["het_sites"]
    pse <- values["phase_errors_def1"] / values["phase_sites_def1"]
    pse <- round(100 * pse, 1)
    disc <- round(100 * values["disc_errors"] / values["dist_n"], 1)
    list(pse = pse, sites = sites)
  })
  a
}

# d0:truth, d1: quilt2, d2:glimpse2, d3:quilt1, d4:glimpse1
acc_phasing <- function(d0, d1, d2, d3, d4) {
  id <- intersect(intersect(intersect(intersect(rownames(d0), rownames(d1)), rownames(d2)), rownames(d3)), rownames(d4))
  quilt2 <- acc_phasing_single_matrix(d1, d0, id)
  glimpse2 <- acc_phasing_single_matrix(d2, d0, id)
  quilt1 <- acc_phasing_single_matrix(d3, d0, id)
  glimpse1 <- acc_phasing_single_matrix(d4, d0, id)
  res <- list(QUILT2 = quilt2, GLIMPSE2 = glimpse2, QUILT1 = quilt1, GLIMPSE1 = glimpse1)
  res
}

## input is matrix
r2_by_freq <- function(breaks, af, truthG, testDS, which_snps = NULL, flip = FALSE, per_snp = FALSE) {
  if (!is.null(which_snps)) {
    af <- af[which_snps]
    truthG <- truthG[which_snps, ]
    testDS <- testDS[which_snps, ]
  }
  truthG <- as.matrix(truthG)
  testDS <- as.matrix(testDS)
  af <- as.numeric(af)
  if (flip) {
    w <- af > 0.5
    af[w] <- 1 - af[w]
    truthG[w, ] <- 2 - truthG[w, ]
    testDS[w, ] <- 2 - testDS[w, ]
  }
  x <- cut(af, breaks = breaks)
  if (ncol(truthG) > 1 && per_snp) {
    # for multiple sample, calculate r2 per snp then average them
    cors_per_af <- tapply(1:length(x), x, function(w) {
      c(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        simple = mean(sapply(w, function(ww) {
          cor(truthG[ww, ], testDS[ww, ], use = "pairwise.complete")**2
        }), na.rm = TRUE),
        norm = mean(sapply(w, function(ww) {
          cor(truthG[ww, ] - 2 * af[ww], testDS[ww, ] - 2 * af[ww], use = "pairwise.complete")**2
        }), na.rm = TRUE)
      )
    })
  } else {
    cors_per_af <- tapply(1:length(x), x, function(w) {
      c(
        n = length(w),
        nA = sum(truthG[w, ], na.rm = TRUE),
        simple = cor(as.vector(truthG[w, ]), as.vector(testDS[w, ]), use = "pairwise.complete")**2,
        norm = cor(as.vector(truthG[w, ] - 2 * af[w]), as.vector(testDS[w, ] - 2 * af[w]), use = "pairwise.complete")**2
      )
    })
  }
  # fill with NA for AF bins without SNPs
  cors_per_af <- t(sapply(cors_per_af, function(a) {
    if (is.null(a[1])) {
      return(c(n = NA, nA = NA, simple = NA, norm = NA))
    }
    a
  }))
  return(cors_per_af)
}

## acc_r2_by_af <- function(d0, d1, af, bins) {
##   id <- intersect(intersect(rownames(d0), rownames(d1)), names(af))
##   res <- r2_by_freq(breaks = bins, af[id], truthG = d0[id,], testDS = d1[id,])
##   as.data.frame(cbind(bin = bins[-1], single = res[, "simple"], orphan = res[, "simple"]))
## }

# d0:truth, d1: quilt2, d2:glimpse2, d3:quilt1, d4:glimpse1
acc_r2_by_af <- function(d0, d1, d2, d3, d4, af, bins, flip = TRUE, per_snp = FALSE) {
  id <- intersect(intersect(intersect(intersect(rownames(d0), rownames(d1)), rownames(d2)), rownames(d3)), rownames(d4))
  id <- intersect(id, names(af))
  res1 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d1, which_snps = id, flip = flip, per_snp = per_snp)
  res2 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d2, which_snps = id, flip = flip, per_snp = per_snp)
  res3 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d3, which_snps = id, flip = flip, per_snp = per_snp)
  res4 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d4, which_snps = id, flip = flip, per_snp = per_snp)
  as.data.frame(cbind(bin = bins[-1], nsnps = res1[, "n"], quilt2 = res1[, "simple"], glimpse2 = res2[, "simple"], quilt1 = res3[, "simple"], glimpse1 = res4[, "simple"]))
}

quilt_r2_by_freq <- function(breaks, af, truthG, testDS, which_snps = NULL, flip = FALSE) {
  if (flip) {
    w <- af > 0.5
    af[w] <- 1 - af[w]
    truthG[w] <- 2 - truthG[w]
    testDS[w] <- 2 - testDS[w]
  }
  if (!is.null(which_snps)) {
    af <- af[which_snps]
    truthG <- truthG[which_snps]
    testDS <- testDS[which_snps]
  }
  x <- cut(af, breaks = breaks) # <NA> is x
  cors_per_af <- tapply(1:length(x), x, function(w) {
    c(
      n = length(w),
      nA = sum(truthG[w], na.rm = TRUE),
      simple = cor(truthG[w], testDS[w], use = "pairwise.complete")**2,
      norm = cor(truthG[w] - 2 * af[w], testDS[w] - 2 * af[w], use = "pairwise.complete")**2
    )
  })
  cors_per_af <- t(sapply(cors_per_af, function(a) {
    if (is.null(a[1])) {
      return(c(n = NA, nA = NA, simple = NA, norm = NA))
    }
    a
  }))
  return(cors_per_af)
}

gettimes <- function(ss) {
  sapply(strsplit(ss, ":"), function(s) {
    s <- as.numeric(s)
    n <- length(s)
    sum(sapply(1:n, function(i) {
      s[i] * 60^(n - i)
    }))
  })
}

rmnull <- function(l) {
  l <- l[!sapply(l, is.null)]
  unlist(l)
}

rmna <- function(l) {
  l <- l[!sapply(l, is.na)]
  unlist(l)
}

parse.quilt.gts <- function(fn) {
  d1 <- fread(cmd = paste("awk '{for(i=1;i<=NF;i=i+3) printf $i\" \"; print \"\"}'", fn), data.table = F)
  id <- d1[, 1]
  d1 <- as.matrix(sapply(d1[, -1], as.numeric))
  rownames(d1) <- id
  return(d1)
}

## https://stackoverflow.com/questions/14984989/how-to-avoid-warning-when-introducing-nas-by-coercion
## as.num(c("1", "2", "X"), na.strings="X")
as.num <- function(x, na.strings = "NA") {
        # stopifnot(is.character(x))
        na = x %in% na.strings
        x[na] = "0"
        x = as.numeric(x)
        x[na] = NA_real_
        x
}

## (gt0, gt1, ds) x nsamples
parse.imputed.gts2 <- function(fn) {
  d1 <- fread(fn, data.table = F)
  id <- d1[, 1]
  d1 <- as.matrix(suppressWarnings(sapply(d1[, -1], as.numeric)))
  rownames(d1) <- id
  return(d1)
}

mycols <- c(QUILT2 = "#e69f00", GLIMPSE2 = "#d55e00", QUILT1 = "#56b4e9", GLIMPSE1 = "#cc79a7")
