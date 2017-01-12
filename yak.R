library("Rcpp")

library("microbenchmark")

sizes <- c(500, 1000, 5000)
crits <- c(10, 20, 30)

scales <- c(1, 2, 4, 8, 16, 32)

sourceCpp("yak_1.cpp")

scaleTime <-
    function(target, scales, prob, sizes, crits) {

        vapply(scales,
               function(ind, prob, sizes, crits) {
                   summary(microbenchmark(target(prob,
                                                 ind * sizes,
                                                 ind * crits)),
                           unit="ms")$min
               }, 0, prob, sizes, crits)
    }



simple <- scaleTime(gsProbs_1, scales, 0.01, sizes, crits)


sourceCpp("yak_2.cpp")

better <- scaleTime(gsProbs_2, scales, 0.01, sizes, crits)

Sys.setenv(PKG_LIBS="-lfftw3")
sourceCpp("yak_3.cpp", rebuild=TRUE)

fft <- scaleTime(gsProbs_3, scales, 0.01, sizes, crits)

Sys.setenv(PKG_LIBS="-lfftw3")
sourceCpp("yak_4.cpp", rebuild=TRUE)

fft.2 <- scaleTime(gsProbs_4, scales, 0.01, sizes, crits)

matplot(scales,
        cbind(simple,
              better,
              fft,
              fft.2), log='xy')


gsProbs_2(0.01, sizes, crits)
