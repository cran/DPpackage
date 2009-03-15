### DPpsBF.R                   
### Calculates the Pseudo Bayes Factor.
###
### Copyright: Alejandro Jara, 2006-2009.
###
### Last modification: 01-07-2006.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### The author's contact information:
###
###      Alejandro Jara
###      Department of Statistics
###      Facultad de Ciencias Físicas y Matemáticas
###      Universidad de Concepción
###      Avenida Esteban Iturra S/N
###      Barrio Universitario
###      Concepción
###      Chile
###      Voice: +56-41-2203163  URL  : http://www2.udec.cl/~ajarav
###      Fax  : +56-41-2251529  Email: ajarav@udec.cl
###

"DPpsBF"<-
function(...) 
{
	model.list <- list(...)
	M <- length(model.list)  
	this.call <- match.call()
	this.call.string <- deparse(this.call)
	model.names <- paste("Model", 1:M, sep=" ")
	
	#for (i in 1:M)
	#{
	#    for(j in 1:M)
	#    {
	#	if (class(model.list[[i]])!=class(model.list[[j]]))
	#	{
	#		stop("arguments are not of the same class\n")
  	#  	}
        #    }
  	#}
  
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
				                 log(model.list[[i]]$cpo) -
				                 log(model.list[[j]]$cpo)
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

