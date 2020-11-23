pamCluster <- function(data=NULL,...)
{
  pm <- cluster::pam(data,...);
  
  result <- list(classification = pm$clustering,pam = pm);
  class(result) <- "pamCluster"
  return(result);
}

predict.pamCluster <- function(object,...)
{
  parameters <- list(...);
  testData <- parameters[[1]];
  class <- nearestCentroid(testData,object$pam$medoids)
  result <- list(classification=class)
  return(result);
}