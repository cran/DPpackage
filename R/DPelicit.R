###                    
### Function to Prior elicitation of the precision parameter of a 
### Dirichlet process prior.
### 
### Copyright: Alejandro Jara Vallejos, 2006
### Last modification: 24-03-2006.


DPelicit<-function(n,method='JGL',a0=NULL,b0=NULL,mean=NULL,std=NULL)
UseMethod("DPelicit")

DPelicit.default<-function(n,method='JGL',a0=NULL,b0=NULL,mean=NULL,std=NULL)
{
        cl <- match.call()
        
	if(method=='JGL')
	{
		if(is.null(a0))
		{
			foo<-newton(c(2,2),n,mean,std**2,fn1)
			inp<-matrix(c(foo[1],foo[2]),ncol=2)
			out<-matrix(c(mean,std),ncol=2)
			colnames(inp)<-c("a0","b0")
			colnames(out)<-c("Mean","Std.")
			rownames(inp)<-c(" ")
			rownames(out)<-c(" ")

		}
		else
		{
			mean<-(a0/b0)*(digamma((a0+n*b0)/n)-digamma(a0/b0))		
			var<-(a0/b0)*log(1+n*b0/a0)-a0/b0+(((a0/b0)*(trigamma((a0+n*b0)/n)-
			     trigamma(a0/b0))+digamma((a0+n*b0)/n)-digamma(a0/b0))**2)*(a0/b0**2)
			
			inp<-matrix(c(a0,b0),ncol=2)
			out<-matrix(c(mean,sqrt(var)),ncol=2)
			colnames(inp)<-c("a0","b0")
			colnames(out)<-c("Mean","Std.")
			rownames(inp)<-c(" ")
			rownames(out)<-c(" ")
		}
	
	}

	if(method=='KMQ')
	{
		if(is.null(a0))
		{
			foo<-newton(c(2,2),n,mean,std**2,fn2)

			inp<-matrix(c(foo[1],foo[2]),ncol=2)
			out<-matrix(c(mean,std),ncol=2)
			colnames(inp)<-c("a0","b0")
			colnames(out)<-c("Mean","Std.")
			rownames(inp)<-c(" ")
			rownames(out)<-c(" ")
		}
		else
		{
			mean<-(a0/b0)*log(1+n*b0/a0)
   			var<-(a0/b0)*log(1+n*b0/a0)-a0/b0+((log(1+n*b0/a0)-n*b0/
			     (a0+n*b0))**2)*(a0/b0**2)

			inp<-matrix(c(a0,b0),ncol=2)
			out<-matrix(c(mean,sqrt(var)),ncol=2)
			colnames(inp)<-c("a0","b0")
			colnames(out)<-c("Mean","Std.")
			rownames(inp)<-c(" ")
			rownames(out)<-c(" ")
		}	
	
	}
	
	output<-list(inp=inp,out=out,call=cl)
	
	class(output)<-"DPelicit"
	return(output)
}



fn1<-function(y,n,mean,var)
# evaluate the function at a0 and b0
# JGL
{
	x<-exp(y)
	fvec<-rep(0,2)
        fvec[1]<-mean-((x[1]/x[2])*(digamma((x[1]+n*x[2])/n)-digamma(x[1]/x[2])))
        fvec[2]<-var-((x[1]/x[2])*log(1+n*x[2]/x[1])-x[1]/x[2]+(((x[1]/x[2])*(trigamma((x[1]+n*x[2])/n)-
		 trigamma(x[1]/x[2]))+digamma((x[1]+n*x[2])/n)-digamma(x[1]/x[2]))**2)*(x[1]/x[2]**2))
	return(fvec)         
}



fn2<-function(y,n,mean,var)
# evaluate the function at a0 and b0
# KMQ
{
	x<-exp(y)
	fvec<-rep(0,2)
	fvec[1]<-mean-(x[1]/x[2])*log(1+n*x[2]/x[1])
	fvec[2]<-var-((x[1]/x[2])*log(1+n*x[2]/x[1])-x[1]/x[2]+
	         ((log(1+n*x[2]/x[1])-n*x[2]/(x[1]+n*x[2]))**2)*(x[1]/x[2]**2))
	return(fvec)         
}


fjacb<-function(x,fvec,n,mean,var,fn)
#   This function computes forward-difference approximation to Jacobian
#   of the function fn.
#   On input x is the point at which the Jacobian es to be evaluated,
#   fvec is the vector of function values at the point.
#   eps is the approximate square root of the machine precision
#   A.J.V., 2006.
{
	eps<-1e-06
	fjac<-matrix(0,ncol=2,nrow=2)
	for(j in 1:2)
	{
		temp<-x[j]
		h<-eps*abs(temp)
		if(h==0)h<-eps
		x[j]<-temp+h
		h<-x[j]-temp
		f<-fn(x,n,mean,var)
		x[j]<-temp
		for(i in 1:2)
		{
			fjac[i,j]<-(f[i]-fvec[i])/h
		}
	}
	return(fjac)
}



newton<-function(x0,n,mean,var,fn,ntrial=200,tolx=0.00001,tolf=0.00001)
#   Given a initial guess x0 for a root, take ntrial Newton-Raphson steps
#   to improve the root. Stop if the root converges in either summed absolute
#   variable increments tolx or summed absolute function values tolf.
#   A.J.V., 2006.
{
	x<-log(x0)
	p<-rep(0,2)
	for(k in 1:ntrial)
	{
		fvec<-fn(x,n,mean,var)
		fjac<-fjacb(x,fvec,n,mean,var,fn)
		errf<-0
		for(i in 1:2)
		{
			errf<-errf+abs(fvec[i])
		}
		if(errf<=tolf)break

		for(i in 1:2)
		{
			p[i]<--fvec[i]
		}
		p<-solve(fjac, p)
		errx<-0
		for(i in 1:2)
		{
			errx<-errx+abs(p[i])
			x[i]<-x[i]+p[i]
		}
		if(errx<=tolx)break
	}
	return(exp(x))
}

