erf  <- function(x) 2* pnorm(x * sqrt(2)) - 1   # erf(x) = 2/sqrt(pi) int_0^x exp(-t^2) dt
erf2 <- function(x) 2* pnorm(x) - 1             # erf2(x) = erf( x * sqrt(2) )


probLevel <- function( thresh ){
    # thresh : vector of threshold voltages normalized by sigma (0 or positive values)
    numThresh <- length(thresh)
    prob <- numeric(numThresh)
    for(indexLevel in 1:(numThresh - 1)){
        prob[indexLevel] <- 0.5* (erf2(thresh[indexLevel + 1]) - erf2(thresh[indexLevel]))
    }
    prob[numThresh] <- 0.5 - 0.5*erf2(thresh[numThresh])
    return( prob )
}

varDigital <- function( prob, weights ){
    # prob   : vector of probabilities in each level
    # weights: vector of weights
    return( as.numeric(prob %*% weights^2) )
}

expectDigital <- function(thresh){  # expectation between thresholds
    numThresh <- length(thresh)
    expect <- numeric( numThresh )
        for(indexLevel in 1:(numThresh - 1)){
        expect[indexLevel] <-  exp(-0.5*thresh[indexLevel]^2) - exp(-0.5*thresh[indexLevel+1]^2) 
    }
    expect[numThresh] <- exp(-0.5*thresh[numThresh]^2) 
    return(as.numeric(expect))
}

corrDigital <- function(expect, weights){   # calculate <x hat[x]>^2
    # expect : vector of <x> between thresholds
    # weights: vector of weights
    return( as.numeric(weights %*% expect)^2 / pi )
}

digitalEff <- function(thresh){
    prob    <- probLevel( thresh )
    expect  <- expectDigital(thresh)
    weights <- expect / prob / sqrt(2*pi)
    return( corrDigital(expect, weights)  / varDigital(prob, weights) )
}
#-------- Simulated annealing
thresh128 <- c(0, 32, 64, 96, 128, 160, 192, 256)   # Initial thresholds (uniformly spaced)
for(index_gen in 1:500000){        # generation
    thresh_matrix <- matrix(rep(thresh128, 9), ncol=9)
    index_mutation <- sample(2:8,2,replace=F)           # Random selection of two threshold levels to mutate
    step_mutation  <- sample(1:5,2,replace=T)           # Random step size
    thresh_matrix[index_mutation[1], seq(2,8,by=3)] <- thresh_matrix[index_mutation[1], seq(2,8,by=3)] + step_mutation[1]
    thresh_matrix[index_mutation[1], seq(3,9,by=3)] <- thresh_matrix[index_mutation[1], seq(3,9,by=3)] - step_mutation[1]
    thresh_matrix[index_mutation[2], 4:6]           <- thresh_matrix[index_mutation[2], 4:6]           + step_mutation[2]
    thresh_matrix[index_mutation[2], 7:9]           <- thresh_matrix[index_mutation[2], 7:9]           - step_mutation[2]
    effGen <- apply(thresh_matrix/128, 2, digitalEff)   # Calculate quantization efficiency 
    index_selection <- which.max(effGen)                # Select the best efficiency
    if(index_selection > 1){                            # In the case of evolution
        thresh128 <- thresh_matrix[,index_selection]
        cat(sprintf('Gen %d :', index_gen)); cat(thresh128); cat(sprintf(' : eff=%.7f\n', effGen[index_selection]))
    }
}
# Champion Result
# Gen 304 :0 33 67 102 140 183 235 306 : eff=0.9904974
#-------- case study
#
Label <- 'Current ALMA         '
thresh  <- 0.586* (seq(4)-1)
weights <- seq(4) - 0.5
prob    <- probLevel( thresh )
expect  <- expectDigital(thresh)
etaALMA <- corrDigital(expect, weights)  / varDigital(prob, weights)
cat(sprintf('%s : ', Label)); cat(thresh); cat(' / '); cat(weights); cat(sprintf(': efficiency = %.8f\n', etaALMA))
#
Label <- 'Modified weights     '
thresh  <- 0.586* (seq(4)-1)
prob    <- probLevel( thresh )
expect  <- expectDigital(thresh)
weights <- expect / prob / sqrt(2*pi)
etaALMA <- corrDigital(expect, weights)  / varDigital(prob, weights)
cat(sprintf('%s : ', Label)); cat(thresh); cat(' / '); cat(weights); cat(sprintf(': efficiency = %.8f\n', etaALMA))
#
#
Label <- 'Powers of sqrt(2)    '
thresh  <- 0.284* c(0.0, sqrt(2), 2+sqrt(2), 2+3*sqrt(2))
prob    <- probLevel( thresh )
expect  <- expectDigital(thresh)
weights <- expect / prob / sqrt(2*pi)
etaALMA <- corrDigital(expect, weights)  / varDigital(prob, weights)
cat(sprintf('%s : ', Label)); cat(thresh); cat(' / '); cat(weights); cat(sprintf(': efficiency = %.8f\n', etaALMA))
#
#
Label <- 'Equal probability    '
thresh  <- 1.53* qnorm(0.5 + 0.125*(seq(4)-1))
prob    <- probLevel( thresh )
expect  <- expectDigital(thresh)
weights <- expect / prob / sqrt(2*pi)
etaALMA <- corrDigital(expect, weights)  / varDigital(prob, weights)
cat(sprintf('%s : ', Label)); cat(thresh); cat(' / '); cat(weights); cat(sprintf(': efficiency = %.8f\n', etaALMA))
#
#
Label <- 'Uniform 16 levels    '
thresh  <- 0.335* (seq(8)-1)
prob    <- probLevel( thresh )
expect  <- expectDigital(thresh)
weights <- seq(8) - 0.5
etaALMA <- corrDigital(expect, weights)  / varDigital(prob, weights)
cat(sprintf('%s : ', Label)); cat(thresh); cat(' / '); cat(weights); cat(sprintf(': efficiency = %.8f\n', etaALMA))
#
#
Label <- 'Uniform 16 levels, modified weights'
thresh  <- 0.335* (seq(8)-1)
prob    <- probLevel( thresh )
expect  <- expectDigital(thresh)
weights <- expect / prob / sqrt(2*pi)
etaALMA <- corrDigital(expect, weights)  / varDigital(prob, weights)
cat(sprintf('%s : ', Label)); cat(thresh); cat(' / '); cat(weights); cat(sprintf(': efficiency = %.8f\n', etaALMA))
#
#
Label <- 'Equal probability 16 '
thresh  <- 1.589* qnorm(0.5 + 0.0625*(seq(8)-1))
prob    <- probLevel( thresh )
expect  <- expectDigital(thresh)
weights <- expect / prob / sqrt(2*pi)
etaALMA <- corrDigital(expect, weights)  / varDigital(prob, weights)
cat(sprintf('%s : ', Label)); cat(thresh); cat(' / '); cat(weights); cat(sprintf(': efficiency = %.8f\n', etaALMA))
#
#
