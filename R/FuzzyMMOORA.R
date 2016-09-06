#' Implementation of Fuzzy MULTIMOORA Method for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{FuzzyMMOORA} function implements both the Fuzzy Multi-Objetive Optimization by Ration Analysis (MOORA) and the Fuzzy "Full Multiplicative Form" (Fuzzy MULTIMOORA).
#' @param decision The decision matrix (\emph{m} x (\emph{n}*3)) with the values of the \emph{m} alternatives, for the \emph{n} criteria, and multiplied by 3 since they are triangular fuzzy numbers.
#' @param weights A vector of length \emph{n}*3, containing the fuzzy weights for the criteria.
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @return \code{FuzzyMMOORA} returns a data frame which contains the scores and the four rankings calculated (Ratio System, Reference Point, Multiplicative Form and Multi-MOORA ranking).
#' @references Balezentis, T. and Balezentis, A. A Survey on Development and Applications of the Multi-criteria Decision Making Method MULTIMOORA. Journal of Multi-Criteria Decision Analysis, 21(3-4), 209-222, 2014.
#' @examples
#'
#'  d <- matrix(c(0.63,0.42,0.63,0.67,0.8,0.59,0.8,0.84,0.92,0.75,0.92,0.92,0.29,0.71,0.75,
#'  0.42,0.46,0.88,0.92,0.59,0.63,1,1,0.71,0.75,0.59,0.42,0.42,0.92,0.75,0.58,0.59,1,0.88,
#'  0.76,0.75,0.59,0.71,0.42,0.33,0.75,0.88,0.58,0.51,0.88,0.96,0.71,0.67,0.5,0.67,0.67,
#'  0.67,0.67,0.84,0.84,0.84,0.84,0.92,0.96,0.96,0.67,0.54,0.54,0.25,0.84,0.71,0.71,0.42,
#'  0.96,0.88,0.88,0.59,0.67,0.71,0.42,0.25,0.84,0.88,0.59,0.42,0.96,0.96,0.75,0.58,0.54,
#'  0.625,0.625,0.295,0.705,0.79,0.795,0.46,0.88,0.92,0.875,0.62),nrow=4,ncol=24)
#'  w <- c(1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24,
#'  1/24,1/24,1/24,1/24,1/24,1/24,1/24,1/24)
#'  cb <- c('max','max','max','max','max','max','max','max')
#'  FuzzyMMOORA(d,w,cb)

FuzzyMMOORA <- function(decision, #matrix with all the alternatives
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


  #MMOORA method

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
    denominator = (sum(decision[,i]^2+decision[,i+1]^2+decision[,i+2]^2))^(1/2)
    N[,i] = decision[,i]/denominator
    N[,i+1] = decision[,i+1]/denominator
    N[,i+2] = decision[,i+2]/denominator
  }

  W <- diag(weights)
  NW <- N%*%W

  #2. Ration system
  NR <- NW
  for(i in seq(1, ncol(decision), 3)){
    if (new_cb[i] == 'min'){
      NR[,i] <- NW[,i]*(-1)
      NR[,i+1] <- NW[,i+1]*(-1)
      NR[,i+2] <- NW[,i+2]*(-1)
    }
  }

  RS <- matrix(nrow = nrow(decision), ncol = 3)

  RS[,1] <- apply(NR[,seq(1, ncol(decision), 3)], 1, sum)
  RS[,2] <- apply(NR[,seq(2, ncol(decision), 3)], 1, sum)
  RS[,3] <- apply(NR[,seq(3, ncol(decision), 3)], 1, sum)

  Def_RS <- ((RS[,3]-RS[,1])+(RS[,2]-RS[,1]))/3+RS[,1]

  #3. Reference point
  Ref <- as.integer(new_cb == "max") * apply(NW, 2, max) +
    as.integer(new_cb == "min") * apply(NW, 2, min)

  RefP <- matrix(nrow = nrow(decision), ncol = (ncol(decision))/3)

  k=1
  for(j in seq(1, ncol(decision), 3)){
    RefP[,k] <- (((NW[,j]-Ref[j])^2+(NW[,j+1]-Ref[j+1])^2+(NW[,j+2]-Ref[j+2])^2))^(1/2)
    k=k+1
  }
  RP <- apply(RefP, 1, max)

  #4. Multiplicative form

  max <- decision
  min <- decision
  for(j in seq(1, ncol(decision), 3)){
    if (new_cb[j] == 'max'){
      min[,j] <- 1
      min[,j+1] <- 1
      min[,j+2] <- 1
    }else{
      max[,j] <- 1
      max[,j+1] <- 1
      max[,j+2] <- 1
    }
  }

  A <- matrix(nrow = nrow(decision), ncol = 3)
  B <- matrix(nrow = nrow(decision), ncol = 3)
  M <- matrix(nrow = nrow(decision), ncol = 3)

  A[,1] <- apply(max[,seq(1, ncol(decision), 3)], 1, prod)
  A[,2] <- apply(max[,seq(2, ncol(decision), 3)], 1, prod)
  A[,3] <- apply(max[,seq(3, ncol(decision), 3)], 1, prod)
  B[,1] <- apply(min[,seq(1, ncol(decision), 3)], 1, prod)
  B[,2] <- apply(min[,seq(2, ncol(decision), 3)], 1, prod)
  B[,3] <- apply(min[,seq(3, ncol(decision), 3)], 1, prod)

  M[,1] <- A[,1]/B[,3]
  M[,2] <- A[,2]/B[,2]
  M[,3] <- A[,3]/B[,1]

  Def_M <- ((M[,3]-M[,1])+(M[,2]-M[,1]))/3+M[,1]

  #5. Ranking the alternatives
  Rrs <- rank(-Def_RS, ties.method= "first")
  Rrp <- rank(RP, ties.method= "first")
  Rm <- rank(-Def_M, ties.method= "first")

  MMRanking = TheoryOfDominance(Rrs,Rrp,Rm,decision)

  return(data.frame(Alternatives = 1:nrow(decision), RatioSystem = Def_RS, Ranking = Rrs, ReferencePoint = RP, Ranking = Rrp, MultiplicativeForm = Def_M, Ranking = Rm, MultiMooraRanking = MMRanking))

}
