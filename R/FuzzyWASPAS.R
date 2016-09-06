#' Implementation of Fuzzy WASPAS Method for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{FuzzyWASPAS} function implements the Fuzzy Weighted Aggregated Sum Product ASsessment (Fuzzy WASPAS) Method.
#' @param decision The decision matrix (\emph{m} x (\emph{n}*3)) with the values of the \emph{m} alternatives, for the \emph{n} criteria, and multiplied by 3 since they are triangular fuzzy numbers.
#' @param weights A vector of length \emph{n}*3, containing the fuzzy weights for the criteria.
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @param lambda A value in [0,1]. It is used in the calculation of the W index.
#' @return \code{FuzzyWASPAS} returns a data frame which contains the score of the W index and the ranking of the alternatives.
#' @references Turskis, Z. and Zavadskas, E. K. and Antucheviciene, J. and Kosareva, N. A Hybrid Model Based on Fuzzy AHP and Fuzzy WASPAS for Construction Site Selection. International Journal of Computers Communications & Control, 10(6), 873-888, 2015.
#' @examples
#'
#'  d <- matrix(c(0.5,0.6,0.6,0.6,0.6,0.7,0.7,0.7,0.7,0.8,0.8,0.8,0.6,0.6,0.8,0.5,0.7,0.7,
#'  0.9,0.6,0.8,0.8,1,0.7,0.8,0.5,0.6,0.6,0.9,0.6,0.7,0.7,1,0.7,0.8,0.8,0.5,0.6,0.5,0.4,0.6,
#'  0.7,0.6,0.5,0.7,0.8,0.7,0.6,0.8,0.7,0.6,0.5,0.9,0.8,0.7,0.6,1,0.9,0.8,0.7,0.5,0.8,0.6,
#'  0.8,0.6,0.9,0.7,0.9,0.7,1,0.8,1,0.4,0.5,0.8,0.7,0.5,0.6,0.9,0.8,0.6,0.7,1,0.9,0.5,0.4,
#'  0.4,0.5,0.6,0.5,0.5,0.6,0.7,0.6,0.6,0.7),nrow=4,ncol=24)
#'  w <- c(0.21,0.28,0.35,0.16,0.20,0.23,0.14,0.16,0.17,0.09,0.12,0.17,0.07,0.08,0.12,0.05,
#'  0.06,0.09,0.03,0.05,0.07,0.01,0.03,0.06)
#'  cb <- c('max','max','max','max','max','max','max','max')
#'  lambda <- 0.49
#'  FuzzyWASPAS(d,w,cb,lambda)

FuzzyWASPAS <- function(decision, #matrix with all the alternatives
                   weights,  #vector with the numeric values of the weights
                   cb,       #vector with the "type" of the criteria (benefit = "max", cost = "min")
                   lambda    #value with the real number of the 'lambda' parameter to calculate W
)
{
  #Checking the arguments
  if(! is.matrix(decision))
    stop("'decision' must be a matrix with the values of the alternatives")
  if(missing(weights))
    stop("a vector containing n weigths should be provided")
#   if(sum(weights[seq(2, length(weights), 3)]) != 1)
#     stop("The sum of 'weights' is not equal to 1")
  if(! is.character(cb))
    stop("'cb' must be a character vector with the type of the criteria")
  if(! all(cb == "max" | cb == "min"))
    stop("'cb' should contain only 'max' or 'min'")
  if(length(weights) != ncol(decision))
    stop("length of 'weights' does not match the number of the criteria")
  if(length(cb) != (ncol(decision)/3))
    stop("length of 'cb' does not match the number of the criteria")
  if(missing(lambda))
    stop("a value for 'lambda' in [0,1] should be provided")


  #Fuzzy WASPAS method

  # Conversion of cb in "fuzzy" values
  new_cb <- c(1:ncol(decision))
  k=1
  for(j in seq(1, ncol(decision), 3)){
    if (cb[k] == 'max'){
      new_cb[j] <- 'max'
      new_cb[j+1] <- 'max'
      new_cb[j+2] <- 'max'
    }
    else{
      new_cb[j] <- 'min'
      new_cb[j+1] <- 'min'
      new_cb[j+2] <- 'min'
    }
    k=k+1
  }

  #1. Normalization

  Norm <- as.integer(new_cb == "max") * apply(decision, 2, max) +
    as.integer(new_cb == "min") * apply(decision, 2, min)

  N <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for(j in seq(1, ncol(decision), 3)){
    if (new_cb[j] == 'max'){
      N[,j] <- decision[,j] / Norm[j+2]
      N[,j+1] <- decision[,j+1] / Norm[j+2]
      N[,j+2] <- decision[,j+2] / Norm[j+2]
    }
    else{
      N[,j] <- Norm[j] / decision[,j+2]
      N[,j+1] <- Norm[j] / decision[,j+1]
      N[,j+2] <- Norm[j] / decision[,j]
    }
  }


  #2. WSM
  W <- diag(weights)
  NW <- N%*%W
  WSM <- matrix(nrow = nrow(decision), ncol = 3)

  WSM[,1] <- apply(NW[,seq(1, ncol(decision), 3)], 1, sum)
  WSM[,2] <- apply(NW[,seq(2, ncol(decision), 3)], 1, sum)
  WSM[,3] <- apply(NW[,seq(3, ncol(decision), 3)], 1, sum)


  #3. WPM
  NW2 <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for(j in seq(1, ncol(decision), 3)){
    NW2[,j] <- N[,j]^weights[j+2]
    NW2[,j+1] <- N[,j+1]^weights[j+1]
    NW2[,j+2] <- N[,j+2]^weights[j]
  }
  WPM <- matrix(nrow = nrow(decision), ncol = 3)

  WPM[,1] <- apply(NW2[,seq(1, ncol(decision), 3)], 1, prod)
  WPM[,2] <- apply(NW2[,seq(2, ncol(decision), 3)], 1, prod)
  WPM[,3] <- apply(NW2[,seq(3, ncol(decision), 3)], 1, prod)


  #4. Q index
  # Defuzzification
  Def_WSM <- (apply(WSM[,1:3], 1, sum))/3
  Def_WPM <- (apply(WPM[,1:3], 1, sum))/3

  Q <- (lambda*Def_WSM) + ((1-lambda)*Def_WPM)

  #5. Ranking the alternatives
  return(data.frame(Alternatives = 1:nrow(decision), WSM = Def_WSM, WPM = Def_WPM, W = Q, Ranking = rank(-Q, ties.method= "first")))

}
