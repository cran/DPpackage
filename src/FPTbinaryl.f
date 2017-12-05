c=======================================================================                      
      subroutine fptbinaryl(model,nrec,p,sens,spec,x,yobs,
     &                      nlink,xlink,      
     &                      a0b0,betapm,betapv,ninter,nlevel,
     &                      mcmcvec,nsave,tune2,propv,
     &                      acrate,fsave,ppsave,randsave,thetasave,cpo,
     &                      alpha,beta,v,y,
     &                      accums,assign,betac,counter,
     &                      endp,eta,etan,
     &                      iflag,intpn,intpo,
     &                      prob,rvecs,
     &                      seed,
     &                      workm1,workm2,workmh1,workv1,workv2) 
c=======================================================================                  
c
c     Subroutine `fptbinaryl' to run a Markov chain in the  
c     binary regression model using a Finite Polya Tree prior 
c     for the link, centered in a logistic distribution.
c
c     Copyright: Alejandro Jara and Tim Hanson, 2006-2010.
c
c     Version 1.0: 
c
c     Last modification: 16-08-2006.
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
c     The authors's contact information:
c
c      Alejandro Jara
c      Department of Statistics
c      Facultad de Matematicas
c      Pontificia Universidad Catolica de Chile
c      Casilla 306, Correo 22 
c      Santiago
c      Chile
c      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
c      Fax  : +56-2-3547729  Email: atjara@uc.cl
c
c      Tim Hanson
c      Division of Biostatistics
c      University of Minnesota
c      School of Public Health
c      A460 Mayo Building, 
c      MMC 303
c      420 Delaware St SE
c      Minneapolis, MN 55455
c      Voice: 612-626-7075  URL  : http://www.biostat.umn.edu/~hanson/
c      Fax  : 612-626-0660  Email: hanson@biostat.umn.edu
c
c----- Data ------------------------------------------------------------
c
c        model       :  integer indicating if the model correct for
c                       missclasification (1) or not (0).
c        nrec        :  integer giving the number of observations.
c        p           :  integer giving the number of fixed coefficients.
c        sens        :  real vector of sensitivity, sens(nrec).
c        spec        :  real vector of specificity, spec(nrec).
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        yobs        :  integer vector giving the oberved binary 
c                       response, yobs(nrec).
c
c-----------------------------------------------------------------------
c
c---- Link information -------------------------------------------------
c       
c        nlink       :  integer giving the number of grid points to 
c                       evaluate the link.
c        xlink       :  real vector giving the grid points, xlink(nlink)
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
c        ninter      :  integer giving the number of final intervals 
c                       in the Finite Polya tree prior.
c        nlevel      :  integer giving the number of binary partitions
c                       in the Finite Polya tree prior.
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
c        tune2       :  real giving the Metropolis tunig parameter
c                       for the precision parameter.
c        
c-----------------------------------------------------------------------
c
c---- Output ----------------------------------------------------------- 
c
c        acrate      :  real vector giving the MH acceptance rate,
c                       acrate(2). 
c        ppsave      :  real matrix containing the mcmc samples for
c                       the probability for each final interval,
c                       ppsave(nsave,ninter).
c        fsave       :  real matrix containing the mcmc samples for
c                       the link, fsave(nsave,nlink).
c        randsave    :  real matrix containing the mcmc samples for
c                       the errors and prediction,
c                       randsave(nsave,nrec+1).
c        thetasave   :  real vector containing the mcmc sample for the
c                       regression parameters, betsave(nsave,p+1). 
c        cpo         :  real giving the cpo, cpo(nrec).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the PT process.
c        beta        :  real vector giving the current value of the 
c                       regression coefficients, beta(p).
c        y           :  integer vector giving the current value of the
c                       true binary responses, y(nrec).
c        v           :  real vector giving the current value of the 
c                       errors, v(nrec).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        accums      :  real matrix giving the accumulated log 
c                       probabilities used in the computation of
c                       final interval probabilities,
c                       accums(nlevel,ninter).
c        alphac      :  real giving the current value of the candidate
c                       for the precision parameter.
c        assign      :  integer matrix giving the possition of each
c                       observation in each partition,
c                       assign(nrec,nlevel).
c        betac       :  real vector giving the current value of the 
c                       candidate for regression parameters, betac(p).
c        cdflogis    :  cdf of a logistic distribution.
c        counter     :  integer matrix giving the number of subjects
c                       in each binary partition, counter(nlevel,ninter).
c        dbet        :  density of a beta distribution.
c        dispcount   :  index. 
c        dlnrm       :  density of a log-normal distribution.
c        endp        :  real vector giving the end of the intervals,
c                       endp(ninter-1).
c        eta         :  real vector giving the linear predictor, 
c                       eta(nrec).
c        etan        :  real vector giving the linear predictor, 
c                       etan(nrec).
c        evali       :  index.
c        i           :  index.
c        iflag       :  integer vector used to evaluate the prior
c                       distribution for the regression coefficients, 
c                       iflag(p).
c        imax        :  index.
c        imin        :  index.
c        intlp       :  index.
c        intpn       :  integer vector giving the interval possition
c                       for the current value of the linear predictor,
c                       intpn(nrec).
c        intpo       :  integer vector giving the interval possition
c                       for the candidate value of the linear predictor,
c                       intpo(nrec).
c        invcdflogis :  quantile function for a logistic distribution.
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        je2         :  index.
c        k           :  index.   
c        k1          :  index.   
c        k2          :  index.  
c        logcgkn     :  real working variable.
c        logcgko     :  real working variable.
c        loglikn     :  real working variable.
c        logliko     :  real working variable.
c        logpriorn   :  real working variable.
c        logprioro   :  real working variable.
c        mu          :  real working variable.
c        nint        :  index
c        npoints     :  index.
c        nscan       :  index.
c        ok          :  integer indicator.
c        pprn        :  index.
c        prob        :  real vector giving the probability of 
c                       the intervals, prob(ninter).
c        quan        :  real working variable.
c        ratio       :  real working variable.
c        rbeta       :  real beta random number generator.
c        rtlnorm     :  real truncated log normal random number generator.
c        runif       :  real uniform random number generator.
c        rvecs       :  real matrix giving the random vectors for the
c                       Polya tree,  rvecs(nlevel,ninter).
c        sd          :  real working variable.
c        sec         :  cpu time working variable.
c        sec0        :  cpu time working variable.
c        sec00       :  cpu time working variable.
c        sec1        :  cpu time working variable.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        skipcount   :  index. 
c        sprint      :  integer function to print on screen.
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        tmp3        :  real used to accumulate quantities.
c        uni         :  real working variable.
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

c+++++Constants
      double precision zero,one
      parameter(zero=0.d0)
      parameter(one =1.d0)

c+++++Observed variables
      integer model,nrec,p,yobs(nrec)
      double precision sens(nrec),spec(nrec)      
      double precision x(nrec,p)

c+++++Link information
      integer nlink 
      double precision xlink(nlink)
      
c+++++Prior information
      integer ninter,nlevel
      double precision a0b0(2),aa0,ab0
      double precision betapm(p),betapv(p,p)

c+++++MCMC parameters
      integer mcmcvec(3),nburn,nskip,nsave,ndisplay
      double precision propv(p,p),tune2

c+++++Stored output
      double precision acrate(2),fsave(nsave,nlink)
      double precision ppsave(nsave,ninter)
      double precision randsave(nsave,nrec+1),thetasave(nsave,p+1)
      double precision cpo(nrec)

c+++++Current values of the parameters
      integer y(nrec)
      double precision alpha,beta(p),v(nrec)

c+++++Working space - CPU time
      double precision sec00,sec0,sec1,sec

c+++++Working space - Distributions
      double precision cdflogis
      double precision dbet,dlnrm
      double precision invcdflogis 

c+++++Working space - General
      integer i,j,je2,k,k1,k2
      integer nint,npoints,pprn,sprint
      double precision mu,quan,sd
      double precision tmp1,tmp2,tmp3
      double precision vpred

c+++++Working space - Linear predictor
      double precision eta(nrec),etan(nrec)

c+++++Working space - Link
      integer assign(nrec,nlevel)
      integer counter(nlevel,ninter)
      double precision accums(nlevel,ninter)
      double precision endp(ninter-1)
      double precision prob(ninter)
      double precision rvecs(nlevel,ninter)

c+++++Working space - MCMC
      integer dispcount,isave,iscan,nscan,skipcount

c+++++Working space - MH regression coefficients
      integer iflag(p)
      integer imax,imin
      integer intlp
      integer intpn(nrec),intpo(nrec)
      integer ok
      double precision betac(p)
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision ratio
      double precision workm1(p,p),workm2(p,p),workmh1(p*(p+1)/2)
      double precision workv1(p),workv2(p)

c+++++Working space - MH precision parameter
      double precision alphac,logcgkn,logcgko

c+++++Working space - Random number generator
      integer seed(2),seed1,seed2
      double precision rbeta,rtlnorm,uni
      real runif

c+++++Working space - True Binary Variable
      integer evali

c++++ initialize variables

      nburn=mcmcvec(1)
      nskip=mcmcvec(2)
      ndisplay=mcmcvec(3)
      
      aa0=a0b0(1)
      ab0=a0b0(2)
     
c++++ set random number generator

      seed1=seed(1)
      seed2=seed(2)

      call setall(seed1,seed2)
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Check consistency with the data before starting 
c++++ the algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c++++ check if the user has requested an interrupt
      call rchkusr()
 
      do i=1,nrec
         if(y(i).eq.1)then
           if(v(i).gt.eta(i))then
              call intpr("i",-1,i,1)
              call intpr("y",-1,y(i),1)
              call dblepr("eta",-1,eta(i),1)
              call dblepr("v",-1,v(i),1)
              call rexit("Errors not consistent with data in S0")
           end if   
         end if
         if(y(i).eq.0)then
           if(v(i).le.eta(i))then
              call intpr("i",-1,i,1)
              call intpr("y",-1,y(i),1)
              call dblepr("eta",-1,eta(i),1)
              call dblepr("v",-1,v(i),1)
              call rexit("Errors not consistent with data in S0")
           end if   
         end if
      end do


c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)
      
      npoints=ninter-1
      
      mu=0.d0
      sd=1.d0

      call cpu_time(sec0)
      sec00=0.d0

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Generating from the full conditionals for the Link function ++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         do i=1,nlevel
            nint=2**i
            do j=1,nint
               counter(i,j)=0
               accums(i,j)=0.d0
               rvecs(i,j)=0.d0
            end do
         end do   
         
         do i=1,nrec
            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdflogis(tmp1,mu,sd,1,0)
            if(v(i).le.quan)then
                assign(i,1)=1
                counter(1,1)=counter(1,1)+1
              else
                assign(i,1)=2
                counter(1,2)=counter(1,2)+1
            end if  
        
            do j=2,nlevel
               nint=2**j
               tmp1=1.d0/dble(nint)            
               k=assign(i,j-1)
               k1=2*(k-1)+1
               k2=2*(k-1)+2
               
               quan=invcdflogis(dble(k1)*tmp1,mu,sd,1,0) 
               
               if(v(i).le.quan)then
                  assign(i,j)=k1
                  counter(j,k1)=counter(j,k1)+1
                 else
                  assign(i,j)=k2
                  counter(j,k2)=counter(j,k2)+1
               end if  
            end do
         end do

c+++++++ generating (Y0,Y1),(Y00,Y01,Y10,Y11),.... 

         tmp3=0.5d0 
         rvecs(1,1)=tmp3
         rvecs(1,2)=1.d0-tmp3
         accums(1,1)=log(tmp3)
         accums(1,2)=log(1.d0-tmp3)

         rvecs(2,1)=tmp3
         rvecs(2,2)=1.d0-tmp3
         accums(2,1)=log(tmp3)+accums(1,1)
         accums(2,2)=log(1.d0-tmp3)+accums(1,1)

         je2=(2)**2
         tmp1=alpha*dble(je2)+counter(2,3)
         tmp2=alpha*dble(je2)+counter(3,4)
         tmp3=rbeta(tmp1,tmp2)

         rvecs(2,3)=tmp3
         rvecs(2,4)=(1.d0-tmp3)
         accums(2,3)=log(tmp3)+accums(1,2)
         accums(2,4)=log(1.d0-tmp3)+accums(1,2)

         do i=2,nlevel-1
            nint=2**i
            je2=(i+1)**2
            do j=1,nint
               k1=2*(j-1)+1
               k2=2*(j-1)+2            
               tmp1=alpha*dble(je2)+dble(counter(i+1,k1))
               tmp2=alpha*dble(je2)+dble(counter(i+1,k2))
               tmp3=rbeta(tmp1,tmp2)
               rvecs(i+1,k1)=tmp3
               rvecs(i+1,k2)=1.d0-tmp3
               accums(i+1,k1)=log(tmp3)+accums(i,j)
               accums(i+1,k2)=log(1.d0-tmp3)+accums(i,j)               
            end do
         end do   

         do i=1,ninter
            prob(i)=exp(accums(nlevel,i))
         end do

c+++++++ end points

         tmp1=1.d0/dble(ninter)  
         do i=1,npoints
            endp(i)=invcdflogis(dble(i)*tmp1,mu,sd,1,0)             
         end do

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to updating regression coefficients       +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ sample candidates

         call rmvnorm(p,beta,propv,workmh1,workv1,betac)

c+++++++ evaluate log-prior for candidate value of parameters

         call dmvn(p,betac,betapm,betapv,logpriorn,workv1,workm1,
     &             workm2,workv2,iflag)  

c+++++++ evaluate log-prior for current value of parameters

         call dmvn(p,beta,betapm,betapv,logprioro,workv1,workm1,
     &             workm2,workv2,iflag)  

      
c+++++++ evaluate log-likelihood for current and candidate value 
c+++++++ of parameters

         logliko=0.d0
         loglikn=0.d0

         do i=1,nrec

c++++++++++ lp candidate
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+x(i,j)*betac(j)                   
            end do
            etan(i)=tmp1            

c++++++++++ possition lp current

            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdflogis(tmp1,mu,sd,1,0)
            if(eta(i).le.quan)then
                intlp=1
              else
                intlp=2
            end if  
        
            do j=2,nlevel
               nint=2**j
               tmp1=1.d0/dble(nint)            
               k=intlp
               k1=2*(k-1)+1
               k2=2*(k-1)+2
               
               quan=invcdflogis(dble(k1)*tmp1,mu,sd,1,0) 
               
               if(eta(i).le.quan)then
                  intlp=k1
                 else
                  intlp=k2
               end if  
            end do
            
            intpo(i)=intlp
            
c++++++++++ possition lp candidate

            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdflogis(tmp1,mu,sd,1,0)
            if(etan(i).le.quan)then
                intlp=1
              else
                intlp=2
            end if  
        
            do j=2,nlevel
               nint=2**j
               tmp1=1.d0/dble(nint)            
               k=intlp
               k1=2*(k-1)+1
               k2=2*(k-1)+2
               
               quan=invcdflogis(dble(k1)*tmp1,mu,sd,1,0) 
               
               if(etan(i).le.quan)then
                  intlp=k1
                 else
                  intlp=k2
               end if  
            end do
            
            intpn(i)=intlp

c++++++++++ likelihood current

            if(intpo(i).eq.1)then
               tmp1=0.d0
               tmp2=prob(1)*
     &              ( dble(ninter)*
     &                cdflogis(eta(i),mu,sd,1,0)
     &              )
              else
               imin=1
               imax=intpo(i)-1
               tmp1=0.d0
               do j=imin,imax
                  tmp1=tmp1+prob(j)
               end do
               tmp2=prob(intpo(i))*
     &              ( dble(ninter)*
     &                cdflogis(eta(i),mu,sd,1,0)-
     &                dble(intpo(i)-1)
     &              )
               tmp2=tmp1+tmp2                    
            end if  
            tmp3=tmp2*sens(i)+(1.d0-spec(i))*(1.d0-tmp2) 
            if(tmp3.lt.zero)go to 100
            if(tmp3.gt.one )go to 100

            if(yobs(i).eq.1)then
               logliko=logliko+log(tmp3)
             else
               logliko=logliko+log(1.d0-tmp3)
            end if   


c++++++++++ likelihood candidate

            if(intpn(i).eq.1)then
               tmp1=0.d0
               tmp2=prob(1)*
     &              ( dble(ninter)*
     &                cdflogis(etan(i),mu,sd,1,0)
     &              )
              else
               imin=1
               imax=intpn(i)-1
               tmp1=0.d0
               do j=imin,imax
                  tmp1=tmp1+prob(j)
               end do
               tmp2=prob(intpn(i))*
     &              ( dble(ninter)*
     &                cdflogis(etan(i),mu,sd,1,0)-
     &                dble(intpn(i)-1)
     &              )
               tmp2=tmp1+tmp2                    
            end if  
            tmp3=tmp2*sens(i)+(1.d0-spec(i))*(1.d0-tmp2) 
            if(tmp3.lt.zero)go to 100
            if(tmp3.gt.one )go to 100

            if(yobs(i).eq.1)then
               loglikn=loglikn+log(tmp3)
             else
               loglikn=loglikn+log(1.d0-tmp3)
            end if   

         end do

c+++++++ acceptance step

         ok=0
         ratio=dexp(loglikn+logpriorn-logliko-logprioro)

         if(dble(runif()).lt.ratio)then
            do j=1,p
               beta(j)=betac(j)
            end do
            do i=1,nrec
               eta(i)=etan(i)
               intpo(i)=intpn(i)
            end do
            acrate(1)=acrate(1)+1.d0
         end if


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ updating the true binary variables               +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(model.eq.1)then
           do i=1,nrec
              if(intpo(i).eq.1)then
                 tmp1=0.d0
                 tmp2=prob(1)*
     &                ( dble(ninter)*
     &                  cdflogis(eta(i),mu,sd,1,0)
     &                )
                else
                 imin=1
                 imax=intpo(i)-1
                 tmp1=0.d0
                 do j=imin,imax
                    tmp1=tmp1+prob(j)
                 end do
                 tmp2=prob(intpo(i))*
     &                ( dble(ninter)*
     &                  cdflogis(eta(i),mu,sd,1,0)-
     &                  dble(intpo(i)-1)
     &                )
                 tmp2=tmp1+tmp2                    
              end if  

              if(yobs(i).eq.1)then
                 tmp3=sens(i)*tmp2/
     &                 (sens(i)*tmp2+(1.d0-spec(i))*(1.d0-tmp2))
               else 
                 tmp3=(1.d0-sens(i))*tmp2/
     &                 ((1.d0-sens(i))*tmp2+spec(i)*(1.d0-tmp2))
              end if

              call rbinom(1,tmp3,evali)
              
              if(evali.ne.0.and.evali.ne.1)then
                call dblepr("prob",-1,tmp3,1) 
                call intpr("evali",-1,evali,1) 
                call rexit("Error in the generation of y")
              end if
              y(i)=evali
           end do
           
         end if

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ updating the errors                              +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,nrec

            uni=dble(runif())
            if(intpo(i).eq.1)then
               tmp1=0.d0
               tmp2=prob(1)*
     &              ( dble(ninter)*
     &                cdflogis(eta(i),mu,sd,1,0)
     &              )
              else
               imin=1
               imax=intpo(i)-1
               tmp1=0.d0
               do j=imin,imax
                  tmp1=tmp1+prob(j)
               end do
               tmp2=prob(intpo(i))*
     &              ( dble(ninter)*
     &                cdflogis(eta(i),mu,sd,1,0)-
     &                dble(intpo(i)-1)
     &              )
               tmp2=tmp1+tmp2                    
            end if           
            
            if(y(i).eq.1)tmp3=uni*tmp2
            if(y(i).eq.0)tmp3=uni+(1-uni)*tmp2

            tmp1=prob(1)
            j=1
            do while(tmp3.gt.tmp1.and.j.lt.ninter)
               j=j+1
               tmp1=tmp1+prob(j)
            end do

            tmp2=(tmp3-tmp1+dble(j)*prob(j))/(dble(ninter)*prob(j))
            v(i)=invcdflogis(tmp2,mu,sd,1,0) 

c++++++++++ Check consistency with the data    

            if(y(i).eq.1)then
              if(v(i).gt.eta(i))then
                 call intpr("i",-1,i,1)
                 call intpr("y",-1,y(i),1)
                 call dblepr("eta",-1,eta(i),1)
                 call dblepr("v",-1,v(i),1)
                 call rexit("Errors not consistent with data in S1")
              end if   
            end if
            if(y(i).eq.0)then
              if(v(i).le.eta(i))then
                 call intpr("i",-1,i,1)
                 call intpr("y",-1,y(i),1)
                 call dblepr("eta",-1,eta(i),1)
                 call dblepr("v",-1,v(i),1)
                 call rexit("Errors not consistent with data in S1")
              end if   
            end if
            
         end do

100      continue

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++ Prediction                                      +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         tmp3=dble(runif())
         tmp1=prob(1)
         j=1
         do while(tmp3.gt.tmp1.and.j.lt.ninter)
            j=j+1
            tmp1=tmp1+prob(j)
         end do
         tmp2=(tmp3-tmp1+dble(j)*prob(j))/(dble(ninter)*prob(j))
         vpred=invcdflogis(tmp2,mu,sd,1,0) 


c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++
         if(aa0.gt.0.d0)then

c++++++++++ sample candidates

            alphac=rtlnorm(log(alpha),tune2*0.2d0,0,0,.true.,.true.)
            logcgkn=dlnrm(alpha ,log(alphac),tune2*0.2d0,1) 
            logcgko=dlnrm(alphac,log(alpha ),tune2*0.2d0,1) 

c++++++++++ evaluate the log-prior for candidate value of parameters

            call dgamma2(alphac,aa0,ab0,logpriorn)  

c++++++++++ evaluate the log-prior for current value of parameters

            call dgamma2(alpha,aa0,ab0,logprioro)

c++++++++++ evaluate the log-likelihood

            je2=(2)**2
            tmp1=alpha*dble(je2)
            tmp2=alpha*dble(je2)
            logliko=dbet(rvecs(2,3),tmp1,tmp2,1)

            tmp1=alphac*dble(je2)
            tmp2=alphac*dble(je2)
            loglikn=dbet(rvecs(2,3),tmp1,tmp2,1)

            do i=2,nlevel-1
               nint=2**i
               je2=(i+1)**2
               do j=1,nint
                  k1=2*(j-1)+1
                  k2=2*(j-1)+2            
                  tmp1=alpha*dble(je2)
                  tmp2=alpha*dble(je2)
                  logliko=logliko+dbet(rvecs(i+1,k1),tmp1,tmp2,1)

                  tmp1=alphac*dble(je2)
                  tmp2=alphac*dble(je2)
                  loglikn=loglikn+dbet(rvecs(i+1,k1),tmp1,tmp2,1)
               end do
            end do   

c++++++++++ acceptance step
            ratio=dexp(loglikn+logpriorn-logliko-logprioro+
     &                 logcgkn-logcgko)

            if(dble(runif()).lt.ratio)then
               alpha=alphac
               acrate(2)=acrate(2)+1.d0
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

c+++++++++++++ cluster information
               thetasave(isave,p+1)=alpha

c+++++++++++++ link information

               do i=1,ninter
                  ppsave(isave,i)=prob(i)
               end do

               do i=1,nlink
                  nint=2
                  tmp1=1.d0/dble(nint)
                  quan=invcdflogis(tmp1,mu,sd,1,0)
                  if(xlink(i).le.quan)then
                     intlp=1
                   else
                     intlp=2
                  end if  
        
                  do j=2,nlevel
                     nint=2**j
                     tmp1=1.d0/dble(nint)            
                     k=intlp
                     k1=2*(k-1)+1
                     k2=2*(k-1)+2
                     quan=invcdflogis(dble(k1)*tmp1,mu,sd,1,0) 
               
                     if(xlink(i).le.quan)then
                       intlp=k1
                     else
                       intlp=k2
                     end if  
                  end do
            
                  if(intlp.eq.1)then
                     tmp1=0.d0
                     tmp2=prob(1)*
     &                    ( dble(ninter)*
     &                      cdflogis(xlink(i),mu,sd,1,0)
     &                    )
                    else
                     imin=1
                     imax=intlp-1
                     tmp1=0.d0
                     do j=imin,imax
                        tmp1=tmp1+prob(j)
                     end do
                     tmp2=prob(intlp)*
     &                   ( dble(ninter)*
     &                     cdflogis(xlink(i),mu,sd,1,0)-
     &                     dble(intlp-1)
     &                   )
                     tmp2=tmp1+tmp2                    
                  end if  
                  fsave(isave,i)=tmp2
               end do

c+++++++++++++ cpo, errors and predictive information
               
               do i=1,nrec
                  randsave(isave,i)=v(i) 

                  if(intpo(i).eq.1)then
                     tmp1=0.d0
                     tmp2=prob(1)*
     &                   ( dble(ninter)*
     &                     cdflogis(eta(i),mu,sd,1,0)
     &                   )
                   else
                     imin=1
                     imax=intpo(i)-1
                     tmp1=0.d0
                     do j=imin,imax
                        tmp1=tmp1+prob(j)
                     end do
                     tmp2=prob(intpo(i))*
     &                   ( dble(ninter)*
     &                     cdflogis(eta(i),mu,sd,1,0)-
     &                     dble(intpo(i)-1)
     &                   )
                     tmp2=tmp1+tmp2                    
                  end if  
                  tmp3=tmp2*sens(i)+(1.d0-spec(i))*(1.d0-tmp2)

                  if(yobs(i).eq.1)then
                     cpo(i)=cpo(i)+1.0d0/tmp3 
                   else
                     tmp3=1.0d0-tmp3 
                     cpo(i)=cpo(i)+1.0d0/tmp3 
                  end if          
                  
               end do
               
               randsave(isave,nrec+1)=vpred

c+++++++++++++ print
               skipcount = 0
               if(dispcount.ge.ndisplay)then
                  call cpu_time(sec1)
                  sec00=sec00+(sec1-sec0)
                  sec=sec00
                  sec0=sec1
                  pprn=sprint(isave,nsave,sec)
                  dispcount=0
               end if   
            end if         
         end if
      end do
      
      acrate(1)=acrate(1)/dble(nscan)      
      acrate(2)=acrate(2)/dble(nscan)      
      
      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do

      return
      end
     

c=======================================================================     
      subroutine fptpredl(ninter,nlevel,npred,nsave,lp,ppsave,out)
c=======================================================================
      implicit none
      integer i,j,jj,k,k1,k2
      integer imin,imax
      integer ninter,nlevel,npred,nsave
      integer intlp,nint
      double precision lp(npred,nsave),out(npred,nsave),
     1  ppsave(nsave,ninter)
      double precision cdflogis,invcdflogis,quan,tmp1,tmp2
      
      do i=1,npred
         do j=1,nsave
            nint=2
            tmp1=1.d0/dble(nint)
            quan=invcdflogis(tmp1,0.d0,1.d0,1,0)
            if(lp(i,j).le.quan)then
               intlp=1
             else
               intlp=2
            end if  

            do jj=2,nlevel
               nint=2**jj
               tmp1=1.d0/dble(nint)            
               k=intlp
               k1=2*(k-1)+1
               k2=2*(k-1)+2
               quan=invcdflogis(dble(k1)*tmp1,0.d0,1.d0,1,0) 
            
               if(lp(i,j).le.quan)then
                 intlp=k1
               else
                 intlp=k2
               end if  
            end do

            if(intlp.eq.1)then
               tmp1=0.d0
               tmp2=ppsave(j,1)*
     &              ( dble(ninter)*
     &                cdflogis(lp(i,j),0.d0,1.d0,1,0)
     &              )
              else
               imin=1
               imax=intlp-1
               tmp1=0.d0
               do jj=imin,imax
                  tmp1=tmp1+ppsave(j,jj)
               end do
               tmp2=ppsave(j,intlp)*
     &             ( dble(ninter)*
     &               cdflogis(lp(i,j),0.d0,1.d0,1,0)-
     &               dble(intlp-1)
     &             )
               tmp2=tmp1+tmp2                    
            end if
            out(i,j)=tmp2
         end do
      end do
      
      return
      end
      
