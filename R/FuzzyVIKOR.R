#' Implementation of Fuzzy VIKOR Method for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{FuzzyVIKOR} function implements the Fuzzy "VIseKriterijumska Optimizacija I Kompromisno Resenje" (Fuzzy VIKOR) Method.
#' @param decision The decision matrix (\emph{m} x (\emph{n}*3)) with the values of the \emph{m} alternatives, for the \emph{n} criteria, and multiplied by 3 since they are triangular fuzzy numbers.
#' @param weights A vector of length \emph{n}*3, containing the fuzzy weights for the criteria.
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @param v A value in [0,1]. It is used in the calculation of the Q index.
#' @return \code{FuzzyVIKOR} returns a data frame which contains the score of the S, R and Q indixes and the ranking of the alternatives according to Q index.
#' @references Opricovic, S. Fuzzy VIKOR with an application to water resources planning. Expert Systems with Applications, 38(10), 12983-12990, 2011.
#' @examples
#'
#'  d <- matrix(c(38,20,24.58,44.54,33.33,33.86,40.01,21.06,25.87,46.89,33.33,33.86,48,24,
#'  29.75,56.27,43.33,42.32,3.26,2.57,2.82,2.46,2.25,2.47,4.08,2.87,2.97,2.73,2.5,2.74,4.08,
#'  2.87,2.97,2.73,2.62,2.85,43,6,38,60,6,6,47,6,42,62,6,6,48,6,50,68,6,6,10,10,1,0,2,3,10,
#'  10,1,0,2,3,10,10,1,0,2,3),nrow=6,ncol=12)
#'  w <- c(1/12,1/12,1/12,1/12,1/12,1/12,1/12,1/12,1/12,1/12,1/12,1/12)
#'  cb <- c('min','max','min','min')
#'  v <- 0.625
#'  FuzzyVIKOR(d,w,cb,v)

FuzzyVIKOR <- function(decision, #matrix with all the alternatives
                  weights,  #vector with the numeric values of the weights
                  cb,       #vector with the "type" of the criteria (benefit = "max", cost = "min")
                  v         #value with the real number of the 'v' parameter to calculate Q
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
  if(missing(v))
    stop("a value for 'v' in [0,1] should be provided")


  #VIKOR method

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

  #1. Ideal solutions
  posI <- as.integer(new_cb == "max") * apply(decision, 2, max) +
    as.integer(new_cb == "min") * apply(decision, 2, min)
  negI <- as.integer(new_cb == "min") * apply(decision, 2, max) +
    as.integer(new_cb == "max") * apply(decision, 2, min)

  #2. S and R index
  #Normalization
  d <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for(i in seq(1, ncol(decision), 3)){
    if (new_cb[i] == 'max'){
      denominator = posI[i+2]-negI[i]
      d[,i] = (posI[i] - decision[,i+2])/denominator
      d[,i+1] = (posI[i+1] - decision[,i+1])/denominator
      d[,i+2] = (posI[i+2] - decision[,i])/denominator
    }
    else{
      denominator = negI[i+2]-posI[i]
      d[,i] = (decision[,i] - posI[i+2])/denominator
      d[,i+1] = (decision[,i+1] - posI[i+1])/denominator
      d[,i+2] = (decision[,i+2] - posI[i])/denominator
    }
  }

  W <- diag(weights)
  NW <- d%*%W

  S <- matrix(nrow = nrow(decision), ncol = 3)
  R <- matrix(nrow = nrow(decision), ncol = 3)

  S[,1] <- apply(NW[,seq(1, ncol(decision), 3)], 1, sum)
  S[,2] <- apply(NW[,seq(2, ncol(decision), 3)], 1, sum)
  S[,3] <- apply(NW[,seq(3, ncol(decision), 3)], 1, sum)

  R[,1] <- apply(NW[,seq(1, ncol(decision), 3)], 1, max)
  R[,2] <- apply(NW[,seq(2, ncol(decision), 3)], 1, max)
  R[,3] <- apply(NW[,seq(3, ncol(decision), 3)], 1, max)


  #3. Q index
  Q1 <- matrix(nrow = nrow(decision), ncol = 3)
  Q2 <- matrix(nrow = nrow(decision), ncol = 3)

  denominatorS = max(S[,3])-min(S[,1])
  Q1[,1] = (S[,1] - min(S[,3]))/denominatorS
  Q1[,2] = (S[,2] - min(S[,2]))/denominatorS
  Q1[,3] = (S[,3] - min(S[,1]))/denominatorS

  denominatorR = max(R[,3])-min(R[,1])
  Q2[,1] = (R[,1] - min(R[,3]))/denominatorR
  Q2[,2] = (R[,2] - min(R[,2]))/denominatorR
  Q2[,3] = (R[,3] - min(R[,1]))/denominatorR


  if (v==1) {
    Q <- Q1
  } else if (v==0) {
    Q <- Q2
  } else
    Q <- v*Q1 + (1-v)*Q2

  Def_S <- (S[,1]+S[,2]*2+S[,3])/4
  Def_R <- (R[,1]+R[,2]*2+R[,3])/4
  Def_Q <- (Q[,1]+Q[,2]*2+Q[,3])/4

  #4. Ranking the alternatives
  return(data.frame(Alternatives = 1:nrow(decision), S = S, Def_S = Def_S, R = R, Def_R = Def_R, Q = Q, Def_Q = Def_Q, Ranking = rank(Def_Q, ties.method= "first")))

}
