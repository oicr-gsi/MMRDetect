#' Default MMRDclassifier.
#'
#' MMRDclassifier trained on ~270 colorectal cancers.
#'
#' @format R object list of 30
#' 
"MMRDclassifier"




#' MMRDetect classifier
#'
#' @param mutationVariable A list of input variables,"Del_rep_mean","RepIndel_num","MMR_sum","maxcossim"
#' @param classifier provided classifier
#' @return classification result
#' @export
MMRDetect.classify <- function(mutationVariable, classifier = MMRDclassifier) {
  classifyset = mutationVariable[,c("Del_rep_mean","RepIndel_num","MMR_sum","maxcossim")]
  
#  classifyset$RepIndel_num <- classifyset$RepIndel_num/291592 #max(classifyset$RepIndel_num)
  #classifyset$MMR_sum <- classifyset$MMR_sum/1011838.248 #max(classifyset$MMR_sum)
  
  predictedGLM <- stats::predict.glm(classifier, newdata=classifyset, type="response",se.fit=T)
  
  mutationVariable$glm_prob = predictedGLM$fit
  mutationVariable$se.fit = predictedGLM$se.fit
  
  return(mutationVariable)
} 
