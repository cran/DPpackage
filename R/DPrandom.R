"DPrandom"<-
function(object,centered=FALSE,predictive=FALSE)
UseMethod("DPrandom")

"DPrandom.default"<-
function(object,centered=FALSE,predictive=FALSE)
{
   if(is(object, "DPlmm"))
   {
      random<-matrix(0,nrow=object$nsubject,ncol=object$nrandom)
      predp<-rep(0,object$nrandom)
      predsd<-rep(0,object$nrandom)
      predse<-rep(0,object$nrandom)
      predl<-rep(0,object$nrandom)
      predu<-rep(0,object$nrandom)
      predm<-rep(0,object$nrandom)
      
      randommat<-matrix(object$save.state$randsave,
                 nrow=object$mcmc$nsave,ncol=object$nrandom*(object$nsubject+1))
      
      dimnames(randommat)<-dimnames(object$save.state$randsave)
      
      thetamat<-matrix(object$save.state$thetasave,nrow=object$mcmc$nsave, 
                       ncol=object$dimen)
      
      counter<-0
      for(i in 1:object$nsubject){
          for(j in 1:object$nrandom){
              counter<-counter+1
              if(centered)
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter]-
                                     object$save.state$thetasave[,j])
              }
              else
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter])              
              }
          }
      }
      
      dimnames(random)<-list(object$namesre1,object$namesre2)
      z<-list(randomm=random,randommat=randommat,thetamat=thetamat,centered=centered,
              predictive=predictive,nsubject=object$nsubject,nrandom=object$nrandom,
              modelname=object$modelname,call=object$call)
      
      if(predictive==TRUE)
      {
      	 for(i in 1:object$nrandom)
      	 { 		
      	     counter<-counter+1	
      	     if(centered)
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predm[i]<-median(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i]))      	     	

                vec<-object$save.state$randsave[,counter]-object$save.state$thetasave[,i]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	     else
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter])      	     	

                predm[i]<-median(object$save.state$randsave[,counter])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]))      	     	

                vec<-object$save.state$randsave[,counter]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	 }
      	 
      	 predtable <- cbind(predp, predm, predsd, predse , predl , predu)
         dimnames(predtable) <- list(object$namesre2, c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	 z$prediction<-predtable
      }
      class(z)<-c("DPrandom") 
      return(z)
     
   }


   if(is(object, "DPglmm")){
      random<-matrix(0,nrow=object$nsubject,ncol=object$nrandom)
      predp<-rep(0,object$nrandom)
      predsd<-rep(0,object$nrandom)
      predse<-rep(0,object$nrandom)
      predl<-rep(0,object$nrandom)
      predu<-rep(0,object$nrandom)
      predm<-rep(0,object$nrandom)
      
      randommat<-matrix(object$save.state$randsave,
                 nrow=object$mcmc$nsave,ncol=object$nrandom*(object$nsubject+1))
      
      dimnames(randommat)<-dimnames(object$save.state$randsave)
      
      thetamat<-matrix(object$save.state$thetasave,nrow=object$mcmc$nsave, 
                       ncol=object$dimen)
      
      counter<-0
      for(i in 1:object$nsubject){
          for(j in 1:object$nrandom){
              counter<-counter+1
              if(centered)
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter]-
                                     object$save.state$thetasave[,j])
              }
              else
              {
                   random[i,j]<-mean(object$save.state$randsave[,counter])              
              }
          }
      }
      
      dimnames(random)<-list(object$namesre1,object$namesre2)
      z<-list(randomm=random,randommat=randommat,thetamat=thetamat,centered=centered,
              predictive=predictive,nsubject=object$nsubject,nrandom=object$nrandom,
              modelname=object$modelname,call=object$call)
      
      if(predictive==TRUE)
      {
      	 for(i in 1:object$nrandom)
      	 { 		
      	     counter<-counter+1	
      	     if(centered)
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predm[i]<-median(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]-
                                object$save.state$thetasave[,i]))      	     	

                vec<-object$save.state$randsave[,counter]-object$save.state$thetasave[,i]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	     else
      	     {
                predp[i]<-mean(object$save.state$randsave[,counter])      	     	

                predm[i]<-median(object$save.state$randsave[,counter])      	     	

                predsd[i]<-sqrt(var(object$save.state$randsave[,counter]))      	     	

                vec<-object$save.state$randsave[,counter]
                
                n<-length(vec)
                
                alpha<-0.05
                
                alow<-rep(0,2)
                
                aupp<-rep(0,2)
                
       
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                
                predse[i]<-predsd[i]/sqrt(n)
      	     }
      	 }
      	 
      	 predtable <- cbind(predp, predm, predsd, predse , predl , predu)
         dimnames(predtable) <- list(object$namesre2, c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
      	 
      	 z$prediction<-predtable
      }
      class(z)<-c("DPrandom") 
      return(z)
     
   }



   if(is(object, "DPdensity"))
   {

       if (centered) 
       { 
	   stop("This option is not implemented  for DPdensity.\n")
       }


       counter<-0 
       dimen<-object$nvar+object$nvar*(object$nvar+1)/2
       random<-matrix(0,nrow=object$nrec,ncol=dimen)
       randommat=object$save.state$randsave

       for(i in 1:object$nrec)
       {
           for(j in 1:dimen)
           {
                counter<-counter+1
                random[i,j]<-mean(object$save.state$randsave[,counter])       
           }
       }
       
       if(is.null(object$varnames))
       {
          object$varnames<-all.vars(object$call)[1:object$nvar]
       }
       
       namesre<-NULL
       for(i in 1:object$nvar)
       {
           namesre<-c(namesre,paste("mu",object$varnames[i],sep="-"))
       }

       for(i in 1:object$nvar)
       {
           for(j in i:object$nvar)
           {
               if(i==j)namesre<-c(namesre,paste("var",object$varnames[i],sep="-"))
               if(i!=j)
               {
                   tempname<-paste(object$varnames[i],object$varnames[j],sep="-")
                   namesre<-c(namesre,paste("sigma",tempname,sep="-"))
               }    
           }
       }
       
       predtable=NULL
       
       if(predictive==TRUE)
       {       
           start<-object$nrec*(object$nvar+object$nvar*(object$nvar+1)/2)
           predp<-rep(0,dimen)
           predm<-rep(0,dimen)
           predsd<-rep(0,dimen)
           predse<-rep(0,dimen)
           predl<-rep(0,dimen)
           predu<-rep(0,dimen)
           for(i in 1:dimen)
           {
                predp[i]<-mean(object$save.state$randsave[,(start+i)])      	     	
                predm[i]<-median(object$save.state$randsave[,(start+i)])      	     	
                predsd[i]<-sqrt(var(object$save.state$randsave[,(start+i)]))      	     	
                vec<-object$save.state$randsave[,(start+i)]
                n<-length(vec)
                alpha<-0.05
                alow<-rep(0,2)
                aupp<-rep(0,2)
                
                a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                            alow=as.double(alow),aupp=as.double(aupp),PACKAGE="DPpackage")
                predl[i]<-a$alow[1]            
                predu[i]<-a$aupp[1]
                predse[i]<-predsd[i]/sqrt(n)
           }
           predtable <- cbind(predp, predm, predsd, predse , predl , predu)
           dimnames(predtable) <- list(namesre, c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Low","95%HPD-Upp"))
       }

       dimnames(random)<-list(seq(1,object$nrec),namesre)
       
       z<-list(randomm=random,modelname=object$modelname,call=object$call,predictive=predictive,
               nsubject=object$nrec,nrandom=dimen,centered=centered,randommat=randommat,
               prediction=predtable)
               
       class(z)<-c("DPrandom") 
       return(z)
   }

}

