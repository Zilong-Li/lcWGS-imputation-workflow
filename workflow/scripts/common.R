
library(data.table)

## input is matrix
r2_by_freq <- function(breaks, af, truthG, testDS, which_snps = NULL, flip = FALSE, per_snp = FALSE) {
    if (flip) {
        w <- af > 0.5
        af[w] <- 1 - af[w]
        truthG[w,] <- 2 - truthG[w,]
        testDS[w,] <- 2 - testDS[w,]
    }
    if (!is.null(which_snps)) {
        af <- af[which_snps]
        truthG <- truthG[which_snps,]
        testDS <- testDS[which_snps,]
    }
    truthG <- as.matrix(truthG)
    testDS <- as.matrix(testDS)
    x <- cut(af, breaks = breaks)
    if (ncol(truthG) > 1 && per_snp) {
        # for multiple sample, calculate r2 per snp then average them
        cors_per_af <- tapply(1:length(x), x, function(w) {
            c(
                n = length(w),
                nA = sum(truthG[w,], na.rm = TRUE),
                simple = mean(sapply(w, function(ww){ cor(truthG[ww,], testDS[ww,], use = 'pairwise.complete') ** 2 }), na.rm = TRUE),
                norm = mean(sapply(w, function(ww){ cor(truthG[ww,] - 2 * af[ww], testDS[ww,] - 2 * af[ww], use = 'pairwise.complete') ** 2 }), na.rm = TRUE)
            )
        })
    } else {
        cors_per_af <- tapply(1:length(x), x, function(w) {
            c(
                n = length(w),
                nA = sum(truthG[w,], na.rm = TRUE),
                simple = cor(as.vector(truthG[w,]), as.vector(testDS[w,]), use = 'pairwise.complete') ** 2,
                norm = cor(as.vector(truthG[w,] - 2 * af[w]), as.vector(testDS[w,] - 2 * af[w]), use = 'pairwise.complete') ** 2
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

# d1: quilt2, d2:glimpse2, d3:quilt1, d4:glimpse1
acc_r2_by_af <- function(d0, d1, d2, d3, d4, af, bins) {
  id <- intersect(rownames(d0), rownames(d1))
  res1 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d1, which_snps = id)
  res2 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d2, which_snps = id)
  res3 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d3, which_snps = id)
  res4 <- r2_by_freq(breaks = bins, af, truthG = d0, testDS = d4, which_snps = id)
  as.data.frame(cbind(bin = bins[-1], quilt2 = res1[, "simple"], glimpse2 = res2[, "simple"], quilt1 = res3[, "simple"], glimpse1 = res4[, "simple"]))
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
            simple = cor(truthG[w], testDS[w], use = 'pairwise.complete') ** 2,
            norm = cor(truthG[w] - 2 * af[w], testDS[w] - 2 * af[w], use = 'pairwise.complete') ** 2
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
  id <- d1[,1]
  d1 <- as.matrix(sapply(d1[, -1], as.numeric))
  rownames(d1) <- id
  return(d1)
}
