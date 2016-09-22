#' Implementation of MetaRanking function for Multi-Criteria Decision Making Problems.
#'
#' @description The \code{MetaRanking} function internally calls functions  \code{FuzzyMMOORA}, \code{FuzzyTOPSISLinear}, \code{FuzzyTOPSISVector}, \code{FuzzyVIKOR} and \code{FuzzyWASPAS} and then calculates a sum of the their rankings and an aggregated ranking by applying the \code{RankAggreg} package.
#' @param decision The decision matrix (\emph{m} x \emph{n}) with the values of the \emph{m} alternatives, for the \emph{n} criteria.
#' @param weights A vector of length \emph{n}, containing the weights for the criteria. The sum of the weights has to be 1.
#' @param cb A vector of length \emph{n}. Each component is either \code{cb(i)='max'} if the \emph{i-th} criterion is benefit or \code{cb(i)='min'} if the \emph{i-th} criterion is a cost.
#' @param lambda A value in [0,1]. It is used in the calculation of the W index for WASPAS method.
#' @param v A value in [0,1]. It is used in the calculation of the Q index for VIKOR method.
#' @return \code{MetaRanking} returns a data frame which contains the rankings of the Fuzzy Multi-MOORA, Fuzzy TOPSIS (linear transformation and vectorial normalization), Fuzzy VIKOR, Fuzzy WASPAS Methods and the MetaRankings of the alternatives.
#' @examples
#'
#'  d <- matrix(c(0.68,0.4,0.6,0.2,0.4,1.44,0.67,0.9,0.45,0.6,2.2,
#'  0.95,1.2,0.7,0.8,18,8,8,25,6,21,11.5,11.5,32.5,9,24,15,15,40,
#'  12,9,0.66,0.66,0,0,10,2.33,2.33,0.66,0.33,10,4.33,4.33,2.33,
#'  1.66,5,1.33,1.33,5.66,1,7,3,3,7.66,2,8.66,5,5,9.33,3.66,2.33,
#'  0.66,0.33,1.33,1.66,4.33,2,1.33,3,2.66,6.33,3.66,3,5,4.33),
#'  nrow=5,ncol=15)
#'  w <- c(0.189,0.214,0.243,0.397,0.432,0.462,0.065,0.078,0.096,
#'  0.068,0.084,0.106,0.174,0.190,0.207)
#'  cb <- c('min','max','max','min','min')
#'  lambda <- 0.5
#'  v <- 0.5
#'  MetaRanking(d,w,cb,lambda,v)

MetaRanking <- function(decision, #matrix with all the alternatives
                        weights,  #vector with the numeric values of the weights
                        cb,       #vector with the "type" of the criteria (benefit = "max", cost = "min")
                        lambda,   #value with the real number of the 'lambda' parameter to calculate W
                        v         #value with the real number of the 'v' parameter to calculate Q
)
{

  #Multi-MOORA method
  MMoora = FuzzyMMOORA(decision,weights,cb)

  #TOPSIS method
  TopsisV = FuzzyTOPSISVector(decision,weights,cb)
  TopsisL = FuzzyTOPSISLinear(decision,weights,cb)

  #VIKOR method
  Vikor = FuzzyVIKOR(decision,weights,cb,v)

  #WASPAS method
  Waspas = FuzzyWASPAS(decision,weights,cb,lambda)

  #Meta-Ranking
  MetaR = MMoora[,8]+TopsisV[,6]+TopsisL[,3]+Vikor[,14]+Waspas[,5]

  #Ranking Aggregated
  ra = rbind(MMoora[,8],TopsisV[,6],TopsisL[,3],Vikor[,14],Waspas[,5])
  if(nrow(decision)<=10)
    RA = RankAggreg::BruteAggreg(ra, nrow(decision), distance="Spearman")
  else
    RA = RankAggreg::RankAggreg(ra, nrow(decision), method = "GA", distance = "Spearman", verbose=FALSE)

  return(data.frame(Alternatives = 1:nrow(decision), MMOORA = MMoora[,8], TOPSISVector = TopsisV[,6],
                    TOPSISLinear = TopsisL[,3], VIKOR = Vikor[,14], WASPAS = Waspas[,5],
                    MetaRanking_Sum = rank(MetaR, ties.method= "first"), MetaRanking_Aggreg = RA$top.list))

}
