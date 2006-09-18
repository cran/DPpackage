    
c=======================================================================                      
      subroutine ptlmp(maxm,ngrid,nrec,p,x,y,
     &                 a0b0,betapm,betapv,tau,  
     &                 mcmc,nsave,propv,tune2,tune3,
     &                 seed,
     &                 acrate,randsave,thetasave,cpo,f,
     &                 alpha,beta,mu,sigma2,v,
     &                 betac,iflag,vc,workm1,workm2,
     &                 workmh1,workv1,workv2,grid)
c=======================================================================                      
c
c     Version 1.0: 
c     Last modification: 08-10-2006.
c
c     Subroutine `ptlmp' to run a Markov chain in the  
c     semiparametric linear regression model using a partially 
c     specified Mixture of Polya trees. 
c
c     Copyright: Alejandro Jara Vallejos, 2006
c
c     This program is free software; you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published by
c     the Free Software Foundation; either version 2 of the License, or (at
c     your option) any later version.
c
c     This program is distributed in the hope that it will be useful, but
c     WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     General Public License for more details.
c
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c
c     The author's contact information:
c
c     Alejandro Jara Vallejos
c     Biostatistical Centre
c     Katholieke Universiteit Leuven
c     U.Z. Sint-Rafaël
c     Kapucijnenvoer 35
c     B-3000 Leuven
c     Voice: +32 (0)16 336892 
c     Fax  : +32 (0)16 337015 
c     URL  : http://student.kuleuven.be/~s0166452/
c     Email: Alejandro.JaraVallejos@med.kuleuven.be
c
c---- Data -------------------------------------------------------------
c
c        maxm        :  integer giving the number of binary partitions
c                       in the Finite Polya tree prior.
c        ngrid       :  integer giving the size of the grid where
c                       the density estimate is evaluated.
c        nrec        :  integer giving the number of observations.
c        p           :  integer giving the number of fixed coefficients.
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        y           :  real vector giving the limits of the responses,
c                       y(nrec).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        aa0, ab0    :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(aa0,ab0). If aa0<0 the precision 
c                       parameter is considered as a constant.
c        betapm      :  real vector giving the prior mean of regression
c                       coefficients, betapm(p).
c        betapv      :  real matrix giving the prior covariance of 
c                       regression coefficients, betapv(p,p).
c        tau1, tau2  :  reals giving the hyperparameters of the prior 
c                       distribution for the inverse of the 
c                       variance of the normal baseline distribution,
c                       1/sigma ~ Gamma(tau1/2,tau2/2).
c
c-----------------------------------------------------------------------
c
c---- MCMC parameters --------------------------------------------------
c
c        nburn       :  integer giving the number of burn-in scans.
c        ndisplay    :  integer giving the number of saved scans to be
c                       displayed on screen.
c        nskip       :  integer giving the thinning interval.
c        nsave       :  integer giving the number of scans to be saved.
c        propv       :  real matrix giving the variance of the normal
c                       proposal for the mh algorithm, propv(p,p).
c        tune2       :  Metropolis tunnig parameter for the scale
c                       parameter.
c        tune3       :  Metropolis tunnig parameter for the precision
c                       parameter.
c        
c-----------------------------------------------------------------------
c
c---- Random numbers ---------------------------------------------------
c
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        seed3       :  seed for random number generation.
c        
c-----------------------------------------------------------------------
c
c---- Output ----------------------------------------------------------- 
c
c        acrate      :  real giving the MH acceptance rate. 
c        cpo         :  real giving the cpo, cpo(nrec).
c        f           :  real vector giving the density estimate at the
c                       grid, f(ngrid).
c        randsave    :  real matrix containing the mcmc samples for
c                       the errors and prediction,
c                       randsave(nsave,nrec+1).
c        thetasave   :  real vector containing the mcmc sample for the
c                       regression parameters, betsave(nsave,p+2). 
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Polya tree process.
c        beta        :  real vector giving the current value of the 
c                       regression coefficients, beta(p).
c        mu          :  real giving the mean of the baseline 
c                       distribution.
c        sigma2      :  real giving the current value of the variance
c                       of the baseline distribution.
c        v           :  real vector giving the current value of the 
c                       errors, v(nrec).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        ainf        :  logical working variable.
c        alphac      :  real giving the candidate value of the precision
c                       parameter.
c        betac       :  real vector giving the candidate value of the 
c                       candidate for regression parameters, betac(p).
c        binf        :  logical working variable.
c        countern    :  index.
c        countero    :  index.
c        denom       :  real working variable.
c        dispcount   :  index. 
c        dlnrm       :  density of a lognormal normal distribution.
c        dnrm        :  density of a normal distribution.
c        grid        :  real vector giving the grid where the density
c                       estimate is evaluated, grid(ngrid).
c        i           :  index.
c        iflag       :  integer vector used to evaluate the prior
c                       distribution for the regression coefficients, 
c                       iflag(p).
c        invcdfnorm  :  inverse of the cdf of a normal distribution.
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        je2         :  index. 
c        k           :  index. 
c        k1          :  index. 
c        k2          :  index. 
c        l           :  index. 
c        linf        :  real working variable.
c        logcgko     :  real working variable.
c        logcgkn     :  real working variable.
c        loglikec    :  real working variable.
c        loglikeo    :  real working variable.
c        logpriorc   :  real working variable.
c        logprioro   :  real working variable.
c        lsup        :  real working variable.
c        nint        :  integer indicator.
c        nscan       :  index.
c        ok          :  integer indicator.
c        prob        :  real working variable.
c        quan        :  real working variable.
c        ratio       :  real working variable.
c        rgamma      :  real gamma random number generator.
c        rtlnorm     :  real truncated log normal random number generator.
c        rtnorm      :  real truncated normal random number generator.
c        runif       :  real uniform random number generator.
c        sd          :  real working variable.
c        sdc         :  real working variable.
c        sisgma2c    :  real working variable.
c        sprint      :  integer function to print on screen.
c        skipcount   :  index. 
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        tmp3        :  real used to accumulate quantities.
c        vc          :  real vector giving the candidate value of the 
c                       errors, vc(nrec).
c        vpred       :  real working variable.
c        workm1      :  real matrix used to update the fixed effects,
c                       workm1(p,p).
c        workm2      :  real matrix used to update the fixed effects,
c                       workm2(p,p).
c        workmh1     :  real vector used to update the fixed effects,
c                       workmh1(p*(p+1)/2).
c        workv1      :  real vector used to update the fixed effects,
c                       workv1(p).
c        workv2      :  real vector used to update the fixed effects,
c                       workv2(p).
c
c=======================================================================                  
     
      implicit none 
      
c+++++Partially specified Polya Trees parameter      
      integer maxm
     
c+++++Observed variables
      integer ngrid,nrec,p  
      real*8 x(nrec,p),y(nrec)

c+++++Prior information
      real*8 a0b0(2),aa0,ab0,betapm(p),betapv(p,p)
      real*8 tau(2),tau1,tau2
      
c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      real*8 propv(p,p),tune2,tune3

c+++++Random numbers
      integer seed(3),seed1,seed2,seed3

c+++++Stored output
      real*8 acrate(3),randsave(nsave,nrec+1),thetasave(nsave,p+2)
      real*8 cpo(nrec),f(ngrid)
      
c+++++Current values of the parameters
      real*8 alpha,beta(p),mu,sigma2,v(nrec)
      
c+++++Working space
      integer countern,countero
      integer dispcount
      integer i
      integer iflag(p)
      integer iscan
      integer isave
      integer j
      integer je2
      integer k
      integer k1,k2
      integer l
      integer nint
      integer nscan
      integer ok
      integer sprint
      integer skipcount
      real*8 alphac
      real*8 betac(p)
      real*8 denom
      real*8 dlnrm
      real*8 dnrm
      real*8 grid(ngrid)
      real*8 invcdfnorm
      real*8 linf
      real*8 logcgkn,logcgko,loglikec,loglikeo,logpriorc,logprioro
      real*8 lsup
      real*8 prob
      real*8 quan
      real*8 ratio
      real*8 rtlnorm
      real*8 rtnorm
      real*8 sd,sdc
      real*8 sigma2c
      real*8 tmp1,tmp2,tmp3
      real*8 vc(nrec)
      real*8 vpred
      real*8 workm1(p,p),workm2(p,p),workmh1(p*(p+1)/2)
      real*8 workv1(p),workv2(p)
      
      real runif
      
      logical ainf,binf

c+++++CPU time
      real*8 sec00,sec0,sec1,sec

      
c++++ initialize variables

      aa0=a0b0(1)
      ab0=a0b0(2)

      mu=0.d0 
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

      tau1=tau(1)
      tau2=tau(2)
      
      sd=sqrt(sigma2)


c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      seed3=seed(3)
      
      call setall(seed1,seed2)


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)

      call cpu_time(sec0)
      sec00=0.d0

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ First computation of loglikelihood (to reduce CPU time)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      loglikeo=0.d0

      do i=1,nrec
         call rchkusr()   
         
         tmp1=0.d0
         do j=1,p
            tmp1=tmp1+x(i,j)*beta(j)
         end do
         v(i)=y(i)-tmp1
      

c+++++++ first observation
         if(i.eq.1)then

            loglikeo=dnrm(v(1) ,mu,sd ,1)

c+++++++ following observations
           else

              quan=0.d0
              countero=0
              
              if(v(i).le.quan) then
                 do l=1,i-1
                    if(v(l).le.quan)countero=countero+1
                 end do
               else
                 do l=1,i-1
                    if(v(l).gt.quan)countero=countero+1
                 end do
              end if  

              if(countero.eq.0)go to 100
              
              ok=1
              j=2
              do while(ok.eq.1.and.j.le.maxm)
                 nint=2**j
                 je2=j**2
                 prob=1.d0/dble(nint)
                 k=1
                 
                 quan=invcdfnorm(prob,mu,sd,1,0)
                 
                 do while(v(i).gt.quan.and.k.le.(nint-1))
                    k=k+1
                    if(k.lt.nint)then
                      quan=invcdfnorm(dble(k)*prob,mu,sd,1,0)
                    end if  
                 end do
                 
                 countern=0
                 
                 if(k.eq.1)then
                    do l=1,i-1
                       if(v(l).le.quan)countern=countern+1
                    end do
                  else if(k.eq.nint)then
                    quan=invcdfnorm(dble(k-1)*prob,mu,sd,1,0) 
                    do l=1,i-1
                       if(v(l).gt.quan)countern=countern+1
                    end do
                  else
                    tmp1=invcdfnorm(dble(k-1)*prob,mu,sd,1,0)
                    tmp2=invcdfnorm(dble(k  )*prob,mu,sd,1,0)

                    if(tmp1.ge.tmp2)then
                      call rexit("Error in the limits")
                    end if  
                  
                    do l=1,i-1
                       if(v(l).gt.tmp1.and.v(l).le.tmp2)then
                          countern=countern+1
                       end if   
                    end do
                 end if
                 
                 loglikeo=loglikeo+
     &                    log(2.d0*alpha*dble(je2)+dble(2*countern))-
     &                    log(2.d0*alpha*dble(je2)+dble(  countero))

                 if(countern.eq.0)then
                    ok=0
                  else  
                    countero=countern
                    j=j+1
                 end if   
              end do

100           continue
              
              loglikeo=loglikeo+dnrm(v(i),mu,sd,1)
         end if
      end do


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Scanning the posterior distribution
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()
	
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to update the regression coefficients                   +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ sample candidates

         call rmvnorm(p,beta,propv,workmh1,workv1,betac)
         
c+++++++ evaluate log-prior for candidate value of the parameters

         call dmvn(p,betac,betapm,betapv,logpriorc,workv1,workm1,
     &             workm2,workv2,iflag)  

c+++++++ evaluate log-prior for current value of parameters

         call dmvn(p,beta,betapm,betapv,logprioro,workv1,workm1,
     &             workm2,workv2,iflag)  


c+++++++ evaluate log-likelihood

         loglikec=0.d0

         do i=1,nrec
         
            call rchkusr()   
            
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+x(i,j)*betac(j)
            end do
            vc(i)=y(i)-tmp1

c++++++++++ first observation
            if(i.eq.1)then

                 loglikec=dnrm(vc(1),mu,sd ,1)

c++++++++++ following observations
              else

                 quan=0.d0
                 
                 countero=0
                 
                 if(vc(i).le.quan) then
                     do l=1,i-1
                        if(vc(l).le.quan)countero=countero+1
                     end do
                  else
                     do l=1,i-1
                        if(vc(l).gt.quan)countero=countero+1
                     end do
                 end if  

                 if(countero.eq.0)go to 101

                 ok=1
                 j=2
                 do while(ok.eq.1.and.j.le.maxm)
                 
                    nint=2**j
                    je2=j**2
                    prob=1.d0/dble(nint)
                    k=1
                    
                    quan=invcdfnorm(prob,mu,sd,1,0)
                    
                    do while(vc(i).gt.quan.and.k.le.(nint-1))
                       k=k+1
                       if(k.lt.nint)then
                          quan=invcdfnorm(dble(k)*prob,mu,sd,1,0)
                       end if  
                    end do
                    
                    countern=0
                    
                    if(k.eq.1)then
                       do l=1,i-1
                          if(vc(l).le.quan)countern=countern+1
                       end do
                     else if(k.eq.nint)then
                       quan=invcdfnorm(dble(k-1)*prob,mu,sd,1,0) 
                       do l=1,i-1
                          if(vc(l).gt.quan)countern=countern+1
                       end do
                     else
                       tmp1=invcdfnorm(dble(k-1)*prob,mu,sd,1,0)
                       tmp2=invcdfnorm(dble(k  )*prob,mu,sd,1,0)

                       if(tmp1.ge.tmp2)then
                         call rexit("Error in the limits")
                       end if  
                       
                       do l=1,i-1
                          if(vc(l).gt.tmp1.and.vc(l).le.tmp2)then
                             countern=countern+1
                          end if   
                       end do
                    end if
                    
                    loglikec=loglikec+
     &                       log(2.d0*alpha*dble(je2)+dble(2*countern))-
     &                       log(2.d0*alpha*dble(je2)+dble(  countero))

                    if(countern.eq.0)then
                       ok=0
                     else  
                       countero=countern
                       j=j+1
                    end if   
                 end do

101              continue
                 
                 loglikec=loglikec+dnrm(vc(i),mu,sd,1)

            end if
            
         end do


c+++++++ acceptance step

         ratio=dexp(loglikec+logpriorc-loglikeo-logprioro)

         if(dble(runif()).lt.ratio)then
            do j=1,p
               beta(j)=betac(j)
            end do
            
            do j=1,nrec
               v(j)=vc(j)
            end do
            
            loglikeo=loglikec
            acrate(1)=acrate(1)+1.d0
         end if


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to update the scale parameter                           +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()
	
c+++++++ sample candidates

         sdc=rtlnorm(log(sd),tune2*0.1d0,0,0,.true.,.true.)
         sigma2c=sdc**2

         logcgkn=dlnrm(sd,log(sdc),tune2*0.1d0,1) 
         logcgko=dlnrm(sdc,log(sd ),tune2*0.1d0,1) 


c+++++++ evaluate log-prior for candidate value of the parameters

         call dgamma2(1.d0/(sigma2c),0.5d0*tau1,0.5d0*tau2,logpriorc)  
       

c+++++++ evaluate log-prior for current value of parameters

         call dgamma2(1.d0/sigma2,0.5d0*tau1,0.5d0*tau2,logprioro)  
      

c+++++++ evaluate log-likelihood

         loglikec=0.d0

         do i=1,nrec
         
            call rchkusr()   
            
c++++++++++ first observation
            if(i.eq.1)then

                 loglikec=dnrm(v(1) ,mu,sdc,1)

c++++++++++ following observations
              else

                 quan=0.d0
                 
                 countero=0
                 
                 if(v(i).le.quan) then
                     do l=1,i-1
                        if(v(l).le.quan)countero=countero+1
                     end do
                  else
                     do l=1,i-1
                        if(v(l).gt.quan)countero=countero+1
                     end do
                 end if  

                 if(countero.eq.0)go to 102

                 ok=1
                 j=2
                 do while(ok.eq.1.and.j.le.maxm)
                 
                    nint=2**j
                    je2=j**2
                    prob=1.d0/dble(nint)
                    k=1
                    
                    quan=invcdfnorm(prob,mu,sdc,1,0)
                    
                    do while(v(i).gt.quan.and.k.le.(nint-1))
                       k=k+1
                       if(k.lt.nint)then
                          quan=invcdfnorm(dble(k)*prob,mu,sdc,1,0)
                       end if  
                    end do
                    
                    countern=0
                    
                    if(k.eq.1)then
                       do l=1,i-1
                          if(v(l).le.quan)countern=countern+1
                       end do
                     else if(k.eq.nint)then
                       quan=invcdfnorm(dble(k-1)*prob,mu,sdc,1,0) 
                       do l=1,i-1
                          if(v(l).gt.quan)countern=countern+1
                       end do
                     else
                       tmp1=invcdfnorm(dble(k-1)*prob,mu,sdc,1,0)
                       tmp2=invcdfnorm(dble(k  )*prob,mu,sdc,1,0)

                       if(tmp1.ge.tmp2)then
                         call rexit("Error in the limits")
                       end if  
                       
                       do l=1,i-1
                          if(v(l).gt.tmp1.and.v(l).le.tmp2)then
                             countern=countern+1
                          end if   
                       end do
                    end if
                    
                    loglikec=loglikec+
     &                       log(2.d0*alpha*dble(je2)+dble(2*countern))-
     &                       log(2.d0*alpha*dble(je2)+dble(  countero))

                    if(countern.eq.0)then
                       ok=0
                     else  
                       countero=countern
                       j=j+1
                    end if   
                 end do

102              continue
                 
                 loglikec=loglikec+dnrm(v(i),mu,sdc,1)

            end if
            
         end do


c+++++++ acceptance step

         ratio=dexp(loglikec+logpriorc-loglikeo-logprioro+
     &              logcgkn-logcgko)

         if(dble(runif()).lt.ratio)then
            sd=sdc
            sigma2=sd**2
            loglikeo=loglikec
            acrate(2)=acrate(2)+1.d0
         end if

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(aa0.gt.0.d0)then
	
c+++++++++++ sample candidates

             alphac=rtlnorm(log(alpha),tune3*0.1d0,0,0,.true.,.true.)
             logcgkn=dlnrm(alpha ,log(alphac),tune3*0.1d0,1) 
             logcgko=dlnrm(alphac,log(alpha ),tune3*0.1d0,1) 

c+++++++++++ evaluate log-prior for candidate value of the parameters

             call dgamma2(alphac,aa0,ab0,logpriorc)  

c+++++++++++ evaluate log-prior for current value of parameters

             call dgamma2(alpha,aa0,ab0,logprioro)

c+++++++++++ evaluate log-likelihood

             loglikec=0.d0

             do i=1,nrec
         
                call rchkusr()   
            
c++++++++++++++ first observation
                if(i.eq.1)then

                     loglikec=dnrm(v(1) ,mu,sd,1)

c++++++++++++++ following observations
                  else

                     quan=0.d0
                 
                     countero=0
                 
                     if(v(i).le.quan) then
                         do l=1,i-1
                            if(v(l).le.quan)countero=countero+1
                         end do
                      else
                         do l=1,i-1
                            if(v(l).gt.quan)countero=countero+1
                         end do
                     end if  

                     if(countero.eq.0)go to 103

                     ok=1
                     j=2
                     do while(ok.eq.1.and.j.le.maxm)
                 
                        nint=2**j
                        je2=j**2
                        prob=1.d0/dble(nint)
                        k=1
                    
                        quan=invcdfnorm(prob,mu,sd,1,0)
                    
                        do while(v(i).gt.quan.and.k.le.(nint-1))
                           k=k+1
                           if(k.lt.nint)then
                              quan=invcdfnorm(dble(k)*prob,mu,sd,1,0)
                           end if  
                        end do
                    
                        countern=0
                    
                        if(k.eq.1)then
                           do l=1,i-1
                              if(v(l).le.quan)countern=countern+1
                           end do
                         else if(k.eq.nint)then
                           quan=invcdfnorm(dble(k-1)*prob,mu,sd,1,0) 
                           do l=1,i-1
                              if(v(l).gt.quan)countern=countern+1
                           end do
                         else
                           tmp1=invcdfnorm(dble(k-1)*prob,mu,sd,1,0)
                           tmp2=invcdfnorm(dble(k  )*prob,mu,sd,1,0)

                           if(tmp1.ge.tmp2)then
                             call rexit("Error in the limits")
                           end if  
                       
                           do l=1,i-1
                              if(v(l).gt.tmp1.and.v(l).le.tmp2)then
                                 countern=countern+1
                              end if   
                           end do
                        end if
                    
                        loglikec=loglikec+
     &                     log(2.d0*alphac*dble(je2)+dble(2*countern))-
     &                     log(2.d0*alphac*dble(je2)+dble(  countero))

                        if(countern.eq.0)then
                           ok=0
                         else  
                           countero=countern
                           j=j+1
                        end if   
                     end do

103                  continue
                 
                     loglikec=loglikec+dnrm(v(i),mu,sdc,1)

                end if
            
             end do


c+++++++++++ acceptance step

             ratio=dexp(loglikec+logpriorc-loglikeo-logprioro+
     &              logcgkn-logcgko)

             if(dble(runif()).lt.ratio)then
                alpha=alphac
                acrate(3)=acrate(3)+1.d0
                loglikeo=loglikec
             end if
         end if 


c+++++++ save samples
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1
              
c+++++++++++++ regression coefficient information
               do j=1,p
                  thetasave(isave,j)=beta(j)
               end do

c+++++++++++++ base line information
               thetasave(isave,p+1)=sigma2

c+++++++++++++ precision parameter information
               thetasave(isave,p+2)=alpha


c+++++++++++++ errors and predictive information
               do j=1,nrec
                  randsave(isave,j)=v(j)
               end do

               quan=0.d0
               
               countero=0
               
               if(dble(runif()).le.0.5d0)then
                  k=1
                  do l=1,nrec
                     if(v(l).le.quan)countero=countero+1
                  end do
                else
                  k=2
                  do l=1,nrec
                     if(v(l).gt.quan)countero=countero+1
                  end do
               end if
               
               ok=1
               j=2
               do while(ok.eq.1.and.j.le.maxm)
                  
                  nint=2**j
                  je2=j**2
                  prob=1.d0/dble(nint)
                  
                  k1=2*(k-1)+1
                  k2=2*(k-1)+2
                  
                  countern=0
                    
                  if(k1.eq.1)then
                     quan=invcdfnorm(dble(k1)*prob,mu,sd,1,0) 
                     do l=1,nrec
                        if(v(l).le.quan)countern=countern+1
                     end do
                   else if(k1.eq.nint)then
                     quan=invcdfnorm(dble(k1-1)*prob,mu,sd,1,0) 
                     do l=1,nrec
                        if(v(l).gt.quan)countern=countern+1
                     end do
                   else
                     tmp1=invcdfnorm(dble(k1-1)*prob,mu,sd,1,0)
                     tmp2=invcdfnorm(dble(k1  )*prob,mu,sd,1,0)

                     if(tmp1.ge.tmp2)then
                        call rexit("Error in the limits")
                     end if  
                       
                     do l=1,nrec
                          if(v(l).gt.tmp1.and.v(l).le.tmp2)then
                             countern=countern+1
                          end if   
                     end do
                  end if
                  
                  tmp3=exp(
     &                 log(     alpha*dble(je2)+dble(countern))-
     &                 log(2.d0*alpha*dble(je2)+dble(countero)))                  


                  if(dble(runif()).le.tmp3)then
                      k=k1
                    else
                      k=k2
                      countern=countero-countern
                  end if

                  if(countern.eq.0)then
                     ok=0
                    else 
                     countero=countern
                     j=j+1
                   end if   
               end do
               
               if(j.gt.maxm)j=maxm


c+++++++++++++ Now j indicates the partition and k the interval

               nint=2**j
               prob=1.d0/dble(nint)

               if(k.eq.1)then
                  quan=invcdfnorm(dble(k)*prob,mu,sd,1,0) 

                  ainf=.true.
                  binf=.false.
                  linf=0.d0
                  lsup=quan
                  vpred=rtnorm(mu,sd,linf,lsup,ainf,binf)
                  
                else if(k.eq.nint)then
                  quan=invcdfnorm(dble(k-1)*prob,mu,sd,1,0) 

                  ainf=.false.
                  binf=.true.
                  linf=quan
                  lsup=0.d0
                  vpred=rtnorm(mu,sd,linf,lsup,ainf,binf)

                else
                  tmp1=invcdfnorm(dble(k-1)*prob,mu,sd,1,0)
                  tmp2=invcdfnorm(dble(k  )*prob,mu,sd,1,0)

                  if(tmp1.ge.tmp2)then
                     call rexit("Error in the limits")
                  end if  

                  ainf=.false.
                  binf=.false.
                  linf=tmp1
                  lsup=tmp2
                  vpred=rtnorm(mu,sd,linf,lsup,ainf,binf)
               end if
               
               randsave(isave,nrec+1)=vpred
               
c+++++++++++++ cpo

               do i=1,nrec
                  loglikec=0.d0

                  quan=0.d0
                  denom=1.0d0
                  
                  countero=0

                  if(v(i).le.quan) then
                      do l=1,nrec
                         if(v(l).le.quan.and.l.ne.i)countero=countero+1
                      end do
                    else
                      do l=1,nrec
                         if(v(l).gt.quan.and.l.ne.i)countero=countero+1
                      end do
                  end if  
                 
                  if(countero.eq.0)go to 104

                  ok=1
                  j=2
                  do while(ok.eq.1.and.j.le.maxm)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)
                     k=1

                     denom=prob
                     
                     quan=invcdfnorm(prob,mu,sd,1,0)
                    
                     do while(v(i).gt.quan.and.k.le.(nint-1))
                        k=k+1
                        if(k.lt.nint)then
                          quan=invcdfnorm(dble(k)/dble(nint),mu,sd,1,0)
                        end if  
                     end do
                    
                     countern=0

                     if(k.eq.1)then
                        do l=1,nrec
                           if(v(l).le.quan.and.l.ne.i)then
                              countern=countern+1
                           end if   
                        end do
                      else if(k.eq.nint)then
                        quan=invcdfnorm(dble(k-1)/dble(nint),mu,sd,1,0) 
                        do l=1,nrec
                           if(v(l).gt.quan.and.l.ne.i)then
                              countern=countern+1
                           end if   
                        end do
                      else
                        tmp1=invcdfnorm(dble(k-1)/dble(nint),mu,sd,1,0)
                        tmp2=invcdfnorm(dble(k  )/dble(nint),mu,sd,1,0)

                        if(tmp1.ge.tmp2)then
                          call rexit("Error in the limits")
                        end if  
                     
                        do l=1,nrec
                           if(l.ne.i)then
                           if(v(l).gt.tmp1.and.v(l).le.tmp2)then
                             countern=countern+1
                           end if
                           end if
                        end do
                     end if
                    
                     loglikec=loglikec+
     &                   log(2.d0)+                
     &                   log(     alpha*dble(je2)+dble(  countern))-
     &                   log(2.d0*alpha*dble(je2)+dble(  countero))


                     if(countern.eq.0)then
                         ok=0
                       else  
                         countero=countern
                         j=j+1
                     end if   
                  end do

104               continue
                 
                  
                  loglikec=loglikec+dnrm(v(i),mu,sd,1)

                  cpo(i)=cpo(i)+exp(loglikec)  
            
               end do


c+++++++++++++ density estimate

               do i=1,ngrid
               
                  loglikec=0.d0

                  quan=0.d0
                  denom=1.0d0
                  
                  countero=0
                 
                  if(grid(i).le.quan) then
                      do l=1,nrec
                         if(v(l).le.quan)countero=countero+1
                      end do
                    else
                      do l=1,nrec
                         if(v(l).gt.quan)countero=countero+1
                      end do
                  end if  

                  if(countero.eq.0)go to 105

                  ok=1
                  j=2
                  do while(ok.eq.1.and.j.le.maxm)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)
                     k=1

                     denom=prob
                     
                     quan=invcdfnorm(prob,mu,sd,1,0)
                    
                     do while(grid(i).gt.quan.and.k.le.(nint-1))
                        k=k+1
                        if(k.lt.nint)then
                          quan=invcdfnorm(dble(k)/dble(nint),mu,sd,1,0)
                        end if  
                     end do
                    
                     countern=0
                    
                     if(k.eq.1)then
                        do l=1,nrec
                           if(v(l).le.quan)countern=countern+1
                        end do
                      else if(k.eq.nint)then
                        quan=invcdfnorm(dble(k-1)/dble(nint),mu,sd,1,0) 
                        do l=1,nrec
                           if(v(l).gt.quan)countern=countern+1
                        end do
                      else
                        tmp1=invcdfnorm(dble(k-1)/dble(nint),mu,sd,1,0)
                        tmp2=invcdfnorm(dble(k  )/dble(nint),mu,sd,1,0)

                        if(tmp1.ge.tmp2)then
                          call rexit("Error in the limits")
                        end if  
                     
                        do l=1,nrec
                           if(v(l).gt.tmp1.and.v(l).le.tmp2)then
                             countern=countern+1
                           end if   
                        end do
                     end if
                    
                     loglikec=loglikec+
     &                   log(2.d0)+                
     &                   log(     alpha*dble(je2)+dble(  countern))-
     &                   log(2.d0*alpha*dble(je2)+dble(  countero))


                     if(countern.eq.0)then
                         ok=0
                       else  
                         countero=countern
                         j=j+1
                     end if   
                  end do

105               continue
                 
                  
                  loglikec=loglikec+dnrm(grid(i),mu,sd,1)

                  f(i)=f(i)+exp(loglikec)  
               end do

c+++++++++++++ print
               skipcount = 0
               if(dispcount.ge.ndisplay)then
                  call cpu_time(sec1)
                  sec00=sec00+(sec1-sec0)
                  sec=sec00
                  sec0=sec1
                  tmp1=sprint(isave,nsave,sec)
                  dispcount=0
               end if   
            end if         
         end if
      end do
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ post chain analysis                                +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      acrate(1)=acrate(1)/dble(nscan)        
      acrate(2)=acrate(2)/dble(nscan)        
      acrate(3)=acrate(3)/dble(nscan)        
      
      do i=1,nrec
         cpo(i)=cpo(i)/dble(nsave)
      end do

      do i=1,ngrid
         f(i)=f(i)/dble(nsave)
      end do

      
      return
      end


