#' Implementation of Fuzzy TOPSIS Method for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{FuzzyTOPSISVector} function implements the Fuzzy Technique for Order of Preference by Similarity to Ideal Solution (Fuzzy TOPSIS) Method with the vector normalization procedure.
#' @param decision The decision matrix (\emph{m} x (\emph{n}*3)) with the values of the \emph{m} alternatives, for the \emph{n} criteria, and multiplied by 3 since they are triangular fuzzy numbers.
#' @param weights A vector of length \emph{n}*3, containing the fuzzy weights for the criteria.
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @return \code{FuzzyTOPSISVector} returns a data frame which contains the score of the R index and the ranking of the alternatives.
#' @references Garcia-Cascales, M. S.; Lamata, M. T. and Sanchez-Lozano, J. M. Evaluation of photovoltaic cells in a multi-criteria decision making process. Annals of Operations Research, 199(1), 373-391, 2012.
#' @examples
#'
#'  d <- matrix(c(0.68,0.4,0.6,0.2,0.4,1.44,0.67,0.9,0.45,0.6,2.2,0.95,1.2,0.7,0.8,18,8,8,
#'  25,6,21,11.5,11.5,32.5,9,24,15,15,40,12,9,0.66,0.66,0,0,10,2.33,2.33,0.66,0.33,10,4.33,
#'  4.33,2.33,1.66,5,1.33,1.33,5.66,1,7,3,3,7.66,2,8.66,5,5,9.33,3.66,2.33,0.66,0.33,1.33,
#'  1.66,4.33,2,1.33,3,2.66,6.33,3.66,3,5,4.33),nrow=5,ncol=15)
#'  w <- c(0.189,0.214,0.243,0.397,0.432,0.462,0.065,0.078,0.096,0.068,0.084,0.106,0.174,
#'  0.190,0.207)
#'  cb <- c('min','max','max','min','min')
#'  FuzzyTOPSISVector(d,w,cb)

FuzzyTOPSISVector <- function(decision, #matrix with all the alternatives
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

  denominator = (apply(decision^2,2,sum))^(1/2)
  for(i in seq(1, ncol(decision), 3)){
    N[,i] = decision[,i]/denominator[i+2]
    N[,i+1] = decision[,i+1]/denominator[i+1]
    N[,i+2] = decision[,i+2]/denominator[i]
  }

  W <- diag(weights)
  NW <- N%*%W

  #2. Ideal solutions
  posI <- as.integer(new_cb == "max") * apply(NW, 2, max) +
    as.integer(new_cb == "min") * apply(NW, 2, min)
  negI <- as.integer(new_cb == "min") * apply(NW, 2, max) +
    as.integer(new_cb == "max") * apply(NW, 2, min)

  #3. Distances to the ideal solutions
  distance_pos = matrix(0,nrow = nrow(decision), ncol = ncol(decision))
  distance_neg = matrix(0,nrow = nrow(decision), ncol = ncol(decision))
  for(j in seq(1, ncol(decision), 3)){
    distance_pos[,j] = (NW[,j]-posI[j])^2
    distance_pos[,j+1] = (NW[,j+1]-posI[j+1])^2
    distance_pos[,j+2] = (NW[,j+2]-posI[j+2])^2
    distance_neg[,j] = (NW[,j]-negI[j])^2
    distance_neg[,j+1] = (NW[,j+1]-negI[j+1])^2
    distance_neg[,j+2] = (NW[,j+2]-negI[j+2])^2
  }

  posDis <- matrix(nrow = nrow(decision), ncol = 3)
  negDis <- matrix(nrow = nrow(decision), ncol = 3)
  posDis[,1] <- apply(distance_pos[,seq(1, ncol(decision), 3)], 1, sum)^(1/2)
  posDis[,2] <- apply(distance_pos[,seq(2, ncol(decision), 3)], 1, sum)^(1/2)
  posDis[,3] <- apply(distance_pos[,seq(3, ncol(decision), 3)], 1, sum)^(1/2)
  negDis[,1] <- apply(distance_neg[,seq(1, ncol(decision), 3)], 1, sum)^(1/2)
  negDis[,2] <- apply(distance_neg[,seq(2, ncol(decision), 3)], 1, sum)^(1/2)
  negDis[,3] <- apply(distance_neg[,seq(3, ncol(decision), 3)], 1, sum)^(1/2)

  #4. R index
  R <- matrix(nrow = nrow(decision), ncol = 3)
  denominator = negDis+posDis
  R[,1] <- negDis[,1]/denominator[,3]
  R[,2] <- negDis[,2]/denominator[,2]
  R[,3] <- negDis[,3]/denominator[,1]

  Def_R <- (1/3)*((R[,1]+R[,2]*4+R[,3])/2)

  #5. Rank the alternatives
  return(data.frame(Alternatives = 1:nrow(decision), R = R, Def_R = Def_R, Ranking = rank(-Def_R, ties.method= "first")))

}
