###                    
### Calculates the Pseudo Bayes Factor.
###
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 09-05-2006.

"DPpsBF"<-
function(...) 
{
	model.list <- list(...)
	M <- length(model.list)  
	this.call <- match.call()
	this.call.string <- deparse(this.call)
	model.names <- paste("Model", 1:M, sep=" ")
	
	for (i in 1:M)
	{
	    for(j in 1:M)
	    {
		if (class(model.list[[i]])!=class(model.list[[j]]))
		{
			stop("arguments are not of the same class\n")
  	  	}
            }
  	}
  
  	PsBF.mat <- matrix(0, M, M)
  	PsBF2l.mat <- matrix(0, M, M)
  	rownames(PsBF.mat) <- colnames(PsBF.mat) <-
  	rownames(PsBF2l.mat) <- colnames(PsBF2l.mat) <- model.names

  	PsBF.call <- NULL
  
  	for (i in 1:M)
  	{
		PsBF.call <- c(PsBF.call, model.list[[i]]$call)
		for (j in 1:M)
		{
			if (identical(attr(model.list[[i]], "y"), attr(model.list[[j]], "y")))
			{
				PsBF.mat[i,j] <- exp(sum(
				                 model.list[[i]]$cpo -
				                 model.list[[j]]$cpo
				                 ))
				PsBF2l.mat[i,j] <- 2*log(PsBF.mat[i,j])
      			}	
    		}
  	}
  
  	return(structure(list(PsBF=PsBF.mat, PsBF2l=PsBF2l.mat,
                         PsBF.call=PsBF.call,model.names=model.names),
                   class="DPpsBF"))
}


"print.DPpsBF" <- function(x, ...){

  cat("Models:\n")
  M <- length(x$PsBF.call)
  for (i in 1:M){
    cat("\n")
    cat(rownames(x$PsBF)[i], "<-\n")
    print(x$PsBF.call[[i]])
  }
  cat("\n\n")
  cat("Pseudo Bayes Factors:\n")
  print(x$PsBF, digits=3)
  cat("\n2*log(Pseudo Bayes Factors):\n")
  print(x$PsBF2l, digits=3)
}

