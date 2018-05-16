pday_touse <- function(pday){
  
  if(FALSE){
    if(pday == 1){
      res <- c(pday, pday+1)
    }else if(pday == 48){
      res <- c(47, 48)
    }else{
      res <- c(pday - 1, pday)
    }
    res
  }
  
  pday
}