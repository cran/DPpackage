### HDPdensity.R                   
### Fit a hierarchical DP mixture of normals model
###
### Copyright: Peter Mueller
### Modified by: Alejandro Jara
### Last modification: 07-06-2008.
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
###      Peter Mueller
###      Department of Biostatistics
###      The University of Texas MD Anderson Cancer Center
###      1515 Holcombe Blvd, Unit 447 
###      Houston TX 77030-4009, USA
###      Voice: (713) 563-4296  URL  : http://www.mdanderson.org/departments/biostats
###      Fax  : (713) 563-4243  Email: pmueller@mdanderson.org
###

"HDPdensity"<-
function(formula,study,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail,work.dir=NULL)
UseMethod("HDPdensity")

"HDPdensity.default"<-
function(formula,
         study,
         prior,
         mcmc,
         state,
         status,
         data=sys.frame(sys.parent()),
         na.action=na.fail,
         work.dir=NULL)  
{
       #########################################################################################
       # call parameters
       #########################################################################################

	 cl <- match.call()
	 mf <- match.call(expand.dots = FALSE)
	 m <- match(c("formula", "data","na.action"), names(mf), 0)
	 mf <- mf[c(1, m)]
	 mf$drop.unused.levels <- TRUE
	 mf[[1]] <- as.name("model.frame")
	 mf <- eval(mf, parent.frame())

  	 y <- as.matrix(model.response(model.frame(formula, data=data)))
	 nrec <- nrow(y)
         nvar <- ncol(y)
         x <- as.matrix(model.matrix(formula, data=data))
	 p <- ncol(x)
         p <- p-1

         m <- mcall <- cl <- match.call()
         nm <- names(m)[-1]
         keep <- is.element(nm, c("data", "na.action"))
         for (i in nm[!keep]) m[[i]] <- NULL
         
         allvars <- c(all.vars(formula), all.vars(study))

         Terms <- if (missing(data)) 
              terms(formula)
         else terms(formula, data = data)
         
         cl$fixed <- eval(formula)
         cl$study <- eval(study)
         m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
         environment(m$formula) <- environment(formula)
         m$drop.unused.levels <- TRUE
         m[[1]] <- as.name("model.frame")
         mf <- eval.parent(m)
        
       #########################################################################################
       # data structure
       #########################################################################################
         nrec <- dim(mf)[1]
         y <- mf[,1:nvar]
         x <- mf[,(nvar+1):(nvar+p)]
         study <- mf[,(nvar+p+1)]
  	 namesxm <- colnames(x)

         pz <- nvar + p
         px <- p
         Z <- cbind(y,x)

         ## moments of initial r.e. for initializations..
         mhat <- apply(Z,2,mean)
         V <- var(Z)
         v <- diag(V)

         if (sum(abs(study-round(study)))>0)
         { # not integers
             cat("\n *** Error: study indicators need to be integers 1...nstudies.\n")
             return(-1)
         }
         if (length(unique(study)) != max(study))
         {
            cat("\n *** Error: studies need to be indexed 1...nstudies.\n")
            return(-1)
         }

         nstudies <- max(study)+1
         npatients <- nrec 

       #########################################################################################
       # prior specification
       #########################################################################################

         # prob idiosyncratic cluster
         # p(eps) ~ pe0*I(0)+pe1*I(1)+(1-pe0-pe1)*Be(ae,be)

         if(is.null(prior$pe1))
         { 
            pe1 <- 0.1
         }
         else
         {
            pe1 <- prior$pe1
         } 

         if(is.null(prior$pe0))
         { 
            pe0 <- 0.1
         }
         else
         {
            pe0 <- prior$pe0
         } 

         if (pe1+pe0 >= 1.0)
         {
             cat("\n *** Error: need pe0+pe1 < 1.\n")
             return(-1)
         }

         if(is.null(prior$eps))
         { 
            mcmc.eps <- 1
            eps <- runif(1)
         }
         else
         {
            mcmc.eps <- 0
            eps <- prior$eps
         } 

         if(is.null(prior$ae))
         { 
            ae <- 1
         }
         else
         {
            ae <- prior$ae
         } 
         
         if(is.null(prior$be))
         { 
            be <- 1
         }
         else
         {
            be <- prior$be
         } 

         if(is.null(prior$q))
         { 
            q <- 10
            R <- 0.25*diag(v)
            S.init <- prior$S
            S.prior <- 0
         }
         else
         {
            q <- prior$q
            R <- prior$R 
            S.init <- 0.25*V
            S.prior <- 1
            if(is.null(R))
            {
               R <- 0.25*diag(v)
               cat(" setting default for R = 1/4*diag(var(Z)).\n")
            }
         } 

         if(is.null(prior$a))
         { 
            m.prior <- 0
            m.int <- prior$m
            if (length(m.init) != pz)
            {
                cat("\n *** Error: length(m) != p.\n")
                cat(" Note: dim(m) includes the covariates!.\n")
                return(-1)
            }
            a <- rep(0,pz)
            A <- matrix(0,nrow=pz,ncol=pz)  
         }
         else
         {
            m.prior <- 1
            m.init <- mhat
            a <- prior$a
            A <- prior$A    

            if (length(a) != pz)
            {
                cat("\n *** Error: length(a) != p.\n")
                cat(" Note: dim(a) includes the covariates!.\n")
                return(-1)
            }

            if ( (nrow(A) != pz) | (ncol(A) != pz))
            {
                cat("\n *** Error: dim(A) != (p x p).\n")
                cat(" Note: dim(A) includes the covariates!.\n")
                return(-1)
            }
         } 

         if(is.null(prior$a0))
         { 
            alpha.prior <- 0
            a0 <- 1
            b0 <- 1
            alpha <- prior$alpha
         }
         else
         {
            alpha.prior <- 1
            a0 <- prior$a0
            b0 <- prior$b0
            alpha <- rgamma(1,shape=a0,rate=b0)
         } 


         if(is.null(prior$cc))
         { 
            cc <- 10
            C <- 0.25*diag(v)
            B.init <- prior$B
            B.prior <- 0
         }
         else
         {
            cc <- prior$cc
            C <- prior$C 
            B.init <- diag(v)
            B.prior <- 1
            if ( (nrow(C) != pz) | (ncol(C) != pz))
            {
               cat("\n *** Error: dim(C) != (p x p).\n")
               cat(" Note: dim(C) includes the covariates!.\n")
               return(-1)
            }
         } 

       #########################################################################################
       # warnings
       #########################################################################################
 
         if (cc < pz+2)
             cat(" *** Warning: should use c > p+2 for good mixing MCMC.\n")
         if (q < pz+2)
             cat(" *** Warning: should use q > p+2 for good mixing MCMC.\n")

       #########################################################################################
       # mcmc specification
       #########################################################################################
         n.discard <- mcmc$nburn
         n.batch <- mcmc$nskip
         n.iter <- n.discard + (n.batch+1)*mcmc$nsave
         n.batch <- n.batch + 1
         n.verbose <- mcmc$ndisplay 
         verbose <- 0
         if(n.verbose>0) verbose <- 3

         if(is.null(mcmc$npredupdate))
         {
            n.predupdate <- 100
         }
         else
         {
            n.predupdate <- mcmc$npredupdate
         }     

         if (n.batch >= n.iter/2){
             cat("\n *** Error: n.iter/2 <= n.batch. \n")
             return(-1)
         }

       #########################################################################################
       # parameters depending on status
       #########################################################################################
        
    	 if(status==TRUE)
	 {
            k0 <- npatients
   	 }
	 
      	 if(status==FALSE)
	 {
            alpha <- state$alpha
            m.int <- state$m
            B.int <- state$B
            eps <- state$eps
            k0 <- state$ncluster
	 }    

       #########################################################################################
       # working space
       #########################################################################################
         seed1 <- sample(1:29000,1)
         seed2 <- sample(1:29000,1)

       #########################################################################################
       # change working directory (if requested..)
       #########################################################################################
         if(!is.null(work.dir))
         {
            cat("\n Changing working directory to ",work.dir,"\n")
            old.dir <- getwd()  # by default work in current working directory
            setwd(work.dir)
         }

       #########################################################################################
       # creating files for Peter's functions
       #########################################################################################

         if(is.null(pz))pz <- ncol(Z)
         if(pz != ncol(Z))
         {
            cat("\n *** Error: ncol(Z) != pz. Need to match.\n")
            return(-1)
         }

         out <- file("init.giis",open="w")
         cat("data-npatients ",npatients,"\n",
             "data-nstudies ", nstudies, "\n",
             "model-p1 ", pz, "\n",
             "data-study ", format(study), "\n",
             "mcmc-niter", n.iter, "\n",
             "mcmc-npi ", n.batch, "\n",
             "mcmc-nbatch 5 \n",
             "mcmc-ndiscard ",n.discard, "\n",
             "mcmc-verbose ", verbose, "\n",
             "mcmc-seed1 ", seed1, "\n",
             "mcmc-seed2 ", seed2, "\n",
             "mcmc-eps ", mcmc.eps, "\n",
             "model-eps ", eps, "\n",
             "model-ae ", ae,   "\n",
             "model-be ", be,   "\n",
             "model-pe1 ", pe1, "\n",
             "model-pe0 ", pe0, "\n",
             file=out)
         close(out)

         out <- file("init.mdp",open="w")
         cat(" data-p ",pz,
             "\n model-m_prior ", m.prior,
             "\n model-B_prior ", B.prior,
             "\n model-S_prior ", S.prior,
             "\n model-alpha_prior ", alpha.prior,
             "\n mcmc-npredupdate ", n.predupdate,
             "\n par-S \n", rbind( format(S.init),"\n"),
             "\n par-q ", q,
             "\n par-R  \n", rbind( format(R),"\n"),
             "\n par-B \n", rbind(format(B.init),"\n"),
             "\n par-c ", cc,
             "\n par-C \n ", rbind(format(C),"\n"),
             "\n par-m ", format(m.init),
             "\n par-aa ", format(a),
             "\n par-Ainv \n", rbind( format(solve(A)),"\n"),
             "\n par-alpha ", alpha,
             "\n par-aalpha ", a0, 
             "\n par-balpha ", b0, 
             "\n init-k0 ", k0, "\n",
             file=out)
         close(out)

         out <- file("data-z.giis",open="w")
         cat("data-z \n",file=out)
         write(format(t(Z)),ncolumns=pz,file=out)
         close(out)

       #########################################################################################
       # calling the c code
       #########################################################################################

         if (verbose>0)
           cat(" Printing trace of form:   iteration:k=k0+k1+k2..+kJ\n",
               " where k  = total # clusters and k0 = # clusters in common measure,",
               "       kj = # clusters in idiosyncratic measure j.\n\n")
         .C("hdpmn",
            package="DPpackage")
         if (verbose>1)
         {
            cat("\n\n")
            cat(" Results are saved in files '*.mdp' and '*.giis'. \n")
            cat(" in the working directory \n")
         }

       #########################################################################################
       # save state
       #########################################################################################

         model.name<-"Hierarchical DP mixture of normals model"		


       # return error code
         if(!is.null(work.dir))
         {
           cat("\n Changing working directory back to ",old.dir,"\n")
           setwd(old.dir)
         }


	 z<-list(modelname=model.name,
	         call=cl,
                 prior=prior,
                 mcmc=mcmc,
                 state=state,
                 work.dir=work.dir,
                 pz=pz,
                 px=px)
                 
         cat("\n\n")        

         class(z)<-c("HDPdensity")
         return(z) 
}


###                    
### Predictive information.
###

"predict.HDPdensity"<-
function(object, data.pred=NULL, j=1, r=0, nsim=100, idx.x=NULL, ...)
{
   if(is(object, "HDPdensity"))
   {
         work.dir <- object$work.dir
         p <- object$pz
         px <- object$px 
         X <- data.pred
  
         ## 1. read in covariates for desired posterior predictive
         if(is.null(X))
         {
            cat("\n *** Error: need X, covariates for future patients.")
            cat("\n     specify as file name or matrix.\n")
            return(-1)
         }
         npa <- nrow(X)

         ## 2. change working directory (if requested..)
         if(!is.null(work.dir))
         {
            cat("\n Changing working directory to ",work.dir,"\n")
            old.dir <- getwd()  # by default work in current working directory
            setwd(work.dir)
         }

         ## 3. check paramters
         if (!is.null(p))
         {
             if (p != ncol(X))
             {
                 cat("\n *** Error: need p=ncol(X).")
                 cat("   Note: X should include dummy random effects (use 0, e.g.).\n")
                 return(-1)
             }
         } 

         if (is.null(idx.x))
         {
             cat("\n *** Warning: did not specify the indices of the covariates.\n",
                   "    assuming the last ",px," coordinates are covariates.\n")
             idx.x <- (p-px+1):p
         }
         if (j != abs(round(j))){ # not an integer
             cat("\n *** Error: j needs to 0 or a study index 1..(J+1),\n",
                   "     where J = number of studies.\n")
             return(-1)
         }
         if (!is.element(r,0:1))
         {
             cat("\n *** Error: r=0 or 1 required.\n")
             return(-1)
         }
  
         out <- file("init.p",open="w")
         cat(" n ",nsim,"\n",
             " p ",p,"\n",
             "npa ", npa, "\n",
             "px ", px,"\n",
             "idx_c", idx.x-1,"\n",  # NOTE: -1 for C indexing strting at 0
             "data-z \n",
             rbind(format(t(X)),"\n"),
             file=out)
         close(out)

         cat("\n Computing predictive for (j,r)=(",c(j,r),
                 ") i.e., study ",j)
         if(r==0)        cat(" with ")
         else            cat(" without ")
         cat(" common measure.\n")
         .C("predictN",
            jpredp=as.integer(j),
            rpredp=as.integer(r),
            package="DPpackage")

         fn <- paste("z-",j,r,".p",sep="") # output file
         zout <- as.matrix(read.table(fn))

         ## change working directory back (if changed earlier..)
         if(!is.null(work.dir))
         {
            cat("\n Changing working directory back to ",old.dir,"\n")
            setwd(old.dir)
         }
         return(zout)
   }
}

