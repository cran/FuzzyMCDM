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
#'  d <- matrix(rpois(36, 9), nrow = 4)
#'  w <- c(0.19, 0.2, 0.21, 0.19, 0.2, 0.21, 0.56, 0.6, 0.61)
#'  cb <- c('max','min','max')
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
