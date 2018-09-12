#
# this is the generic genetic power calculator to compute the NCP
#
# fA        = risk allele frequency
# k         = prevalence
# rAa       = genotype relative risk (Aa)
# rAA       = genotype relative risk (AA)
# n_case    = number of cases
# n_control = number of controls
# verbose   = optional flag
#

cc_gpc <- function( fA, k, rAa, rAA, n_case, n_control, verbose )
{
        # frequency of the protective allele (non-risk)
        fa <- 1 - fA

        # calculate genotype freqs
        fAA <- fA * fA
        fAa <- 2 * fA * fa
        faa <- fa * fa

        # baseline risk in aa, i.e. an indiv. may have disease even if phenotype 'aa'; there is some baseline risk of the 'aa' genotype
        raa <- k / ( fA*fA*rAA + 2*fA*fa*rAa + fa*fa )

        # actual (not relative) risk in genotypes
        rrAA <- rAA * raa
        rrAa <- rAa * raa
        rraa <- raa

        # non-relative risk (i.e. how much is disease reduced) in genotypes
        nrrAA <- 1 - rrAA
        nrrAa <- 1 - rrAa
        nrraa <- 1 - rraa

        # odds ratio
        orAa <- ( rrAa / nrrAa ) / ( rraa / nrraa )
        orAA <- ( rrAA / nrrAA ) / ( rraa / nrraa )

        v1_AA <- fAA * rrAA
        v1_Aa <- fAa * rrAa
        v1_aa <- faa * rraa
        v1_sum <- v1_AA + v1_Aa + v1_aa

        v2_AA <- fAA * nrrAA
        v2_Aa <- fAa * nrrAa
        v2_aa <- faa * nrraa
        v2_sum <- v2_AA + v2_Aa + v2_aa

        v1_AA <- v1_AA / v1_sum
        v1_Aa <- v1_Aa / v1_sum
        v1_aa <- v1_aa / v1_sum

        v2_AA <- v2_AA / v2_sum
        v2_Aa <- v2_Aa / v2_sum
        v2_aa <- v2_aa / v2_sum
        case_A <- v1_AA + v1_Aa / 2
        control_A <- v2_AA + v2_Aa / 2

        # make into case-control units
        n_case_A <- 2 * n_case * case_A
        n_case_a <- 2 * n_case * (1 - case_A)
        n_control_A <- 2 * n_control * control_A
        n_control_a <- 2 * n_control * (1 - control_A)

        # complete the chi-square test; the chi-square statistic is equivalent to the ncp
        x2 <- chisq.test(matrix(c(n_case_A, n_control_A, n_case_a, n_control_a), nrow=2, ncol=2), correct=F)
        ncp <- x2$statistic

        if ( verbose ) {
                cat("high risk allele frequency = ", fA, "\n")
                cat("prevalence                 = ", k, "\n")
                cat("genotype relative risk Aa  = ", rAa, "\n")
                cat("genotype relative risk AA  = ", rAA, "\n")
                cat("penetrance aa              = ", rraa, "\n")
                cat("penetrance Aa              = ", rrAa, "\n")
                cat("penetrance AA              = ", rrAA, "\n")
                cat("genotypic odds ratio Aa    = ", orAa, "\n")
                cat("genotypic odds ratio AA    = ", orAA, "\n")
                cat("frequency A in cases       = ", case_A, "\n")
                cat("genotypic odds ratio AA    = ", orAA, "\n")
                cat("frequency A in cases       = ", case_A, "\n")
                cat("frequency A in controls    = ", control_A, "\n")
                cat("NCP (1df)                  = ", ncp, "\n")
        }

        ncp
}


#
# to combine NCP into z-scores, calculate power for stage 1 and stage 2
#
# e1 =          ncp stage 1
# e2 =          ncp stage 2
# pval =        alpha
# r =           fraction of samples (case + control) in stage 1
#

power.meta <- function(e1, e2, pval)
{
        crit.val <- qnorm(pval/2,0,1)

        power1 <- 0
        power2 <- 0
        power.joint <- 0

        # calculate power for stage 1, given ncp > 0
        if ( e1 > 0 ) {

                power1 <- 1-pnorm( abs(crit.val), mean=sqrt(e1), sd=1 )

                #z1 <- rnorm( N, sqrt(e1), 1 )
                #pass1 <- z1[z1 >= abs(crit.val)]
                #power1 <- length(pass1) / N
        }

        # calculate power for stage 2, given ncp > 0
        if ( e2 > 0 ) {

                power2 <- 1-pnorm( abs(crit.val), mean=sqrt(e2), sd=1 )

                #z2 <- rnorm( N, sqrt(e2), 1 )
                #pass2 <- z2[z2 >= abs(crit.val)]
                #power2 <- length(pass2) / N
        }

        # calculate power if either ncp of stage 1 or ncp of stage 2 == 0
        if ( e1 == 0 ) {
                power.joint <- power2
        } else if ( e2 == 0 ) {
                power.joint <- power1
        } else {

                power.joint <- 1-pnorm( abs(crit.val), mean=sqrt(e1+e2), sd=1 )


                        # this is equivalent to
                        #z.joint <- rnorm( N, sqrt(e1+e2), 1 )
                #z.joint <- ( z1 * sqrt(r) + z2 * sqrt(1-r) )
                #pass.joint <- z.joint[z.joint >= abs(crit.val)]
                #power.joint <- length(pass.joint) / N
        }

        p <- c( power1, power2, power.joint )

        p
}

#
# to calculate power for the study
#
# grr =         genotype relative risk
# pval =        pvalue for the study
# maf1 =        minor allele frequency (1)
# maf2 =        minor allele frequency (2)
# flip =        T or F, used to determine which is the risk allele
#

ccPowerMultiPop <- function(grr, pval, maf1, maf2, flip)
{

        # set number of cases/controls in stages 1 and 2
        n_case1 <- 2500
        n_control1 <- 2500
        n_case2 <- 2500
        n_control2 <- 2500

        n_samples1 <- n_case1 + n_control1
        n_samples2 <- n_case2 + n_control2
        n_joint <- n_samples1 + n_samples2

        #R <- 100000
        k <- 0.01

        # flip a coin to determine which is the risk allele
        if ( flip ) {
                raf1 <- maf1
                raf2 <- maf2
        } else {
                raf1 <- 1 - maf1
                raf2 <- 1 - maf2
        }

        ncp1 <- 0
        ncp2 <- 0

        # calculate ncp for stage 1, calling cc_gpc function
        if ( raf1 > 0 && raf1 < 1 ) {
                ncp1 <- cc_gpc( raf1, k, grr, grr^2, n_case1, n_control1, F )
        }
        # calculate ncp for stage 2, calling cc_gpc function
        if ( raf2 > 0 && raf2 < 1 ) {
                ncp2 <- cc_gpc( raf2, k, grr, grr^2, n_case2, n_control2, F )
        }

        # call power.meta to calculate power for stage 1, stage 2, and joint power (if ncp for either stage == 0)
        power <- power.meta( ncp1, ncp2, pval )

}

#
# Use data from inputted values (commands for input below)
# ccPowerMultiPopBigLoop calls all of the above functions to compute power.
# afyri, afceu, afchb are allele frequencies inputted (below) using data from HapMap
#

ccPowerMultiPopBigLoop <- function(afceu, afpop2)
{
        # each iteration (e.g. results.pop.number) yields three numbers - power 1, power 2, joint.power

        results.flip.F <- ccPowerMultiPop( grr, 5e-8, afceu, afpop2, F )
        results.flip.T <- ccPowerMultiPop( grr, 5e-8, afceu, afpop2, T )

        # set up an array to store the calculated powers
        power.array <- rep(NA, 9)

        # depending on the allele frequency in CEU, determine a place to store the calculated power in power.array

        array.index <- 0

        if ( afceu == 0.00 ) {
                array.index <- 1
        }

        if ( afceu > 0.00 && afceu <= 0.01 ) {
                array.index <- 2
        }

        if ( afceu > 0.01 && afceu <= 0.02 ) {
                array.index <- 3
        }

        if ( afceu > 0.02 && afceu <= 0.03 ) {
                array.index <- 4
        }

        if ( afceu > 0.03 && afceu <= 0.04 ) {
                array.index <- 5
        }

        if ( afceu > 0.04 && afceu <= 0.05 ) {
                array.index <- 6
        }
        if ( afceu > 0.05 && afceu <= 0.10 ) {
                array.index <- 7
        }

        if ( afceu > 0.10 && afceu <= 0.25 ) {
                array.index <- 8
        }

        if ( afceu > 0.25 && afceu <= 0.50 ) {
                array.index <- 9
        }

        # calculate joint power for the snp, adding flip=T and flip=F values and averaging them
        power.joint.avg <- ( results.flip.F[3] + results.flip.T[3] ) / 2

        # store the averaged joint power in the array

        power.array[array.index] <- power.joint.avg

        # return the array of calculated power
        power.array

}

# scan commandargs for datatype
x <- 0
repeat {
        x <- x+1
        if (commandArgs()[x] == "-CL") {
        grr <- commandArgs()[x+1]; grr <- as.numeric(substr(grr, 2, nchar(grr)))
        pop2 <- commandArgs()[x+2]; pop2 <- substr(pop2, 2, nchar(pop2))
        letter <- commandArgs()[x+3]; letter <- substr(letter, 2, nchar(letter))
        break
        }
        if (x == length(commandArgs())) {
                print("remember the -CL command!")
                break}
        }

rm(x)

print(grr); print(pop2)

inputfile <- paste("1kg.frequencies.", letter, sep="")
outputfile <- paste("power.10K.", pop2, ".", grr, ".", letter, ".1KG.txt", sep="")

print(inputfile)

# input allele frequencies
freqs <- read.table(inputfile, header=T, stringsAsFactors=F, as.is=T)

# run power calculations
meta <- mapply(ccPowerMultiPopBigLoop, freqs$AFCEU, freqs[,pop2])
#print(meta)

# sum the powers in the rows (will take the average later)
sums <- apply(meta, 1, sum, na.rm=T)
sums <- as.matrix(sums)

# output how many values per column are NOT NA, so that the average can be taken later
not.na <- rowSums(!is.na(meta))

# make "not.na" a 1x9 matrix to output with the sums
not.na <- as.matrix(not.na)

print(not.na)
print(sums)

# create matrix in which to store sums and "not.na.t" values, place sums and not.na.t values into matrix
results <- matrix(data=NA, nrow=1, ncol=18)

# store sums and not.na counts in matrix, alternating between sums and counts
results[1,1] <- sums[1,1]; results[1,2] <- not.na[1,1]
results[1,3] <- sums[2,1]; results[1,4] <- not.na[2,1]
results[1,5] <- sums[3,1]; results[1,6] <- not.na[3,1]
results[1,7] <- sums[4,1]; results[1,8] <- not.na[4,1]
results[1,9] <- sums[5,1]; results[1,10] <- not.na[5,1]
results[1,11] <- sums[6,1]; results[1,12] <- not.na[6,1]
results[1,13] <- sums[7,1]; results[1,14] <- not.na[7,1]
results[1,15] <- sums[8,1]; results[1,16] <- not.na[8,1]
results[1,17] <- sums[9,1]; results[1,18] <- not.na[9,1]

print(results)

#create the output file
write.table(results, outputfile, row.names=F, col.names=c("zero", "zero.count", "one", "one.count", "two", "two.count", "three", "three.count", "four", "four.count", "five", "five.count", "ten", "ten.count", "twentyfive", "twentyfive.count", "fifty", "fifty.count"))

