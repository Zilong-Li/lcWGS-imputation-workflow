
library(data.table)

## input is matrix
r2_by_freq <- function(breaks, af, truthG, testDS, which_snps = NULL, flip = FALSE) {
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
    x <- cut(af, breaks = breaks)
    if (ncol(truthG) > 1) {
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
