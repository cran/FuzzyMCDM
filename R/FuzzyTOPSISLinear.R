#' Implementation of Fuzzy TOPSIS Method for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{FuzzyTOPSISLinear} function implements the Fuzzy Technique for Order of Preference by Similarity to Ideal Solution (Fuzzy TOPSIS) Method with de linear transformation (maximum) as normalization method.
#' @param decision The decision matrix (\emph{m} x (\emph{n}*3)) with the values of the \emph{m} alternatives, for the \emph{n} criteria, and multiplied by 3 since they are triangular fuzzy numbers.
#' @param weights A vector of length \emph{n}*3, containing the fuzzy weights for the criteria.
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @return \code{FuzzyTOPSISLinear} returns a data frame which contains the score of the R index and the ranking of the alternatives.
#' @references Chen, C.T. Extensions of the TOPSIS for group decision-manking under fuzzy environment. Fuzzy Sets and Systems, 114, 1-9, 2000.
#' @examples
#'
#'  d <- matrix(c(5.7,6.3,6.3,7.7,8.3,8,9.3,9.7,9,5,9,7,7,10,9,9,10,10,5.7,8.3,7,7.7,9.7,9,
#'  9,10,10,8.33,9,7,9.67,10,9,10,10,10,3,7,6.3,5,9,8.3,7,10,9.7),nrow=3,ncol=15)
#'  w <- c(0.7,0.9,1,0.9,1,1,0.77,0.93,1,0.9,1,1,0.43,0.63,0.83)
#'  cb <- c('max','max','max','max','max')
#'  FuzzyTOPSISLinear(d,w,cb)

FuzzyTOPSISLinear <- function(decision, #matrix with all the alternatives
                   weights,  #vector with the numeric values of the weights
                   cb        #vector with the "type" of the criteria (benefit = "max", cost = "min")
)
{
  #Checking the arguments
  if(! is.matrix(decision))
    stop("'decision' must be a matrix with the values of the alternatives")
  if(missing(weights))
    stop("a vector containing n weigths, adding up to 1, should be provided")
  #   if(sum(weights[seq(2, length(weights), 3)]) != 1)
  #     stop("The sum of 'weights' is not equal to 1")
  if(! is.character(cb))
    stop("'cb' must be a character vector with the type of the criteria")
  if(! all(cb == "max" | cb == "min"))
    stop("'cb' should contain only 'max' or 'min'")
  if(length(weights) != ncol(decision))
    stop("length of 'weights' does not match the number of the criteria")
  if(length(cb) != ncol(decision)/3)
    stop("length of 'cb' does not match the number of the criteria")


  #TOPSIS method

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

  #1. Normalization and weighting
  N <- matrix(nrow = nrow(decision), ncol = ncol(decision))

  for(i in seq(1, ncol(decision), 3)){
    if (new_cb[i] == 'max'){
      denominator = max(decision[,i+2])
      N[,i] = decision[,i]/denominator
      N[,i+1] = decision[,i+1]/denominator
      N[,i+2] = decision[,i+2]/denominator
    }
    else{
      denominator = min(decision[,i])
      N[,i] = denominator/decision[,i+2]
      N[,i+1] = denominator/decision[,i+1]
      N[,i+2] = denominator/decision[,i]
    }
  }

  W <- diag(weights)
  NW <- N%*%W

  #2. Ideal solutions
  posI <- rep(1,ncol(decision))
  negI <- rep(0,ncol(decision))

  #3. Distances to the ideal solutions
  distance_pos = matrix(0,nrow = nrow(decision), ncol = ncol(decision))
  distance_neg = matrix(0,nrow = nrow(decision), ncol = ncol(decision))
  for(j in seq(1, ncol(decision), 3)){
    distance_pos[,j] = ((1/3)*((NW[,j]-posI[j])^2+(NW[,j+1]-posI[j+1])^2+(NW[,j+2]-posI[j+2])^2))^(1/2)
    distance_neg[,j] = ((1/3)*((NW[,j]-negI[j])^2+(NW[,j+1]-negI[j+1])^2+(NW[,j+2]-negI[j+2])^2))^(1/2)
  }

  posDis <- apply(distance_pos, 1, sum)
  negDis <- apply(distance_neg, 1, sum)

  #4. R index
  R <- negDis/(negDis+posDis)

  #5. Rank the alternatives
  return(data.frame(Alternatives = 1:nrow(decision), R = R, Ranking = rank(-R, ties.method= "first")))

}
