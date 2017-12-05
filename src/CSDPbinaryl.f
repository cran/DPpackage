    
c=======================================================================                      
      subroutine csdpbinaryl(nrec,p,sens,spec,x,yobs,
     &                       nlink,xlink,
     &                       a0b0,betapm,betapv,
     &                       mcmcvec,nsave,propv,extra,      
     &                       acrate,fsave,randsave,thetasave,cpo, 
     &                       alpha,beta,ncluster,theta,y,v,
     &                       betac,
     &                       endp,endp2,eta,etan,
     &                       clusts,     
     &                       index,
     &                       intcount,intcount2,intind,intind2,iflag, 
     &                       intposso,intpossn, 
     &                       limr,lpsav,lpsavc,
     &                       maxint,maxend,maxm,mass,
     &                       massurn1,massurn2,massurn3,massurn4,
     &                       prob,
     &                       proburn1,s,
     &                       seed,urn,uvec,vvec,vnew,
     &                       workm1,workm2,workmh1,workv1,workv2,wvec)
c=======================================================================                  
c
c     Subroutine `cspdbinaryl' to run a Markov chain in the  
c     semiparametric logistic regression model using a Centrally
c     Standarized Dirichlet process prior for the link as in Newton et 
c     al., 1996. 
c
c     Copyright: Alejandro Jara, 2006 - 2011.
c
c     Version 4.0: 
c
c     Last modification: 22-12-2006.
c     
c     Changes and Bug fixes: 
c
c     Version 1.0 to Version 2.0:
c          - Uses qsort3 subroutine to sort errors and endpoints.
c          - Keeps the observations that belong to each interval
c            in memory to allow a faster assignment of the errors.
c          - Performs a binary search to find the endpoints corresponding
c            to each interval.
c          - Computation of G.
c          - Metropolis ratio for theta.
c
c     Version 2.0 to Version 3.0:
c          - Add the link points to generate the partitions. Therefore 
c            the link is partially sampled up to the specified partition.
c
c     Version 3.0 to Version 4.0:
c          - Computation of baseline.
c          - Added strategy for sampling the precision parameter alpha.
c          - Sampling from the partially sampled CSDP 
c          - Added an extra step which moves the clusters (optional). 
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
c---- Data -------------------------------------------------------------
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
c        ppar        :  real value giving the contribution of each
c                       Dirichlet process prior to the final
c                       distribution.
c        d           :  real value giving the distance between 
c                       F{-1}((1+ppar)/2) and F{-1}((1-ppar)/2.
c
c-----------------------------------------------------------------------
c
c---- MCMC parameters --------------------------------------------------
c
c        nburn       :  integer giving the number of burn-in scans.
c        ndisplay    :  integer giving the number of saved scans to be
c                       displayed on screen.
c        nskip       :  integer giving the thinning interval.
c        ntheta      :  integer giving the thinning interval for theta.
c        nsave       :  integer giving the number of scans to be saved.
c        propv       :  real matrix giving the variance of the normal
c                       proposal for the mh algorithm, propv(p,p).
c        
c-----------------------------------------------------------------------
c
c---- Output ----------------------------------------------------------- 
c
c        acrate      :  real vector giving the MH acceptance rate,
c                       acrate(2).
c        fsave       :  real matrix containing the mcmc samples for
c                       the link.
c        randsave    :  real matrix containing the mcmc samples for
c                       the errors and prediction,
c                       randsave(nsave,nrec+1).
c        thetasave   :  real vector containing the mcmc sample for the
c                       regression parameters, thetasave(nsave,p+3). 
c        cpo         :  real giving the cpo, cpo(nrec).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Dirichlet process.
c        beta        :  real vector giving the current value of the 
c                       regression coefficients, beta(p).
c        ncluster    :  integer giving the number of clusters in the
c                       random effects.
c        theta       :  real giving the cutoff determining 4 urns.
c        y           :  integer vector giving the current value of the
c                       true binary responses, y(nrec).
c        v           :  real vector giving the current value of the 
c                       errors, v(nrec).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        area        :  real giving the area of a region.
c        ac          :  real vector giving the mass assigned to each 
c                       urn, ac(4).
c        ac2         :  real vector giving the mass assigned to each 
c                       urn for candidate theta, ac2(4).
c        ainf        :  logical working variable.
c        binf        :  logical working variable.
c        betac       :  real vector giving the current value of the 
c                       candidate for regression parameters, betac(p).
c        cdflogis    :  cdf of a logistic distribution.
c        cdfbaselinel:  cdf of the centrally standarized logistic 
c                       distribution.
c        counter     :  index.
c        dgamlog     :  function to compute the log gamma.
c        dispcount   :  index. 
c        efind       :  integer function to find the intervals where the
c                       errors should be sampled from.
c        endp        :  real vector giving the end of the intervals,
c                       endp(maxend).
c        endp2       :  real vector giving the end of the intervals,
c                       endp2(maxend).
c        eta         :  real vector giving the linear predictor, 
c                       eta(nrec).
c        etan        :  real vector giving the linear predictor, 
c                       etan(nrec).
c        evali       :  integer indicator used in updating the state.
c        i           :  index.
c        iflag       :  integer vector used to evaluate the prior
c                       distribution for the regression coefficients, 
c                       iflag(p).
c        imax        :  index.
c        imin        :  index.
c        index       :  integer working vector, index(maxm)
c        indi        :  index.
c        intcount    :  integer vector used to count the number of 
c                       observations for each interval, 
c                       intcount(maxint).
c        intcount2   :  integer vector used to count the number of 
c                       observations for each interval, 
c                       intcount2(maxint).
c        intcounturn :  integer vector to save the number of 
c                       intervals within each urn, 
c                       intcounturn(4).
c        intind      :  integer vector giving the interval where 
c                       each obbservation belong, intind(nrec+1).
c        intind2     :  integer vector giving the interval where 
c                       each obbservation belong, intind2(nrec+1).
c        intposso    :  integer matrix giving the observations that
c                       belong to each interval, intposso(maxint,nrec+1).
c        intpossn    :  integer matrix giving the observations that
c                       belong to each interval, intpossn(maxint,nrec+1).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index.   
c        keepbeta    :  real working variable.
c        liminf      :  real working variable.
c        limsup      :  real working variable.
c        logliko     :  real working variable.
c        loglikn     :  real working variable.
c        logprioro   :  real working variable.
c        logpriorn   :  real working variable.
c        lpsav       :  integer vector used to store the position of 
c                       the endpoint for the linear predictor,
c                       lpsav(nrec). 
c        lpsavc      :  integer vector used to store the position of 
c                       the endpoint for the linear predictor,
c                       lpsavc(nrec). 
c        maxint      :  integer giving the maximum number of intervals.
c        maxend      :  integer giving the maximum number of intervals.
c        maxm        :  integer giving the maximum number of points in
c                       the finite approximation of the DP.
c        mass        :  real vector giving the mass of the DP,
c                       mass(maxint).
c        massurn1    :  real vector giving the mass for the urn 1 of 
c                       the CSDP, massurn1(maxint).
c        massurn2    :  real vector giving the mass for the urn 2 of 
c                       the CSDP, massurn2(maxint).
c        massurn3    :  real vector giving the mass for the urn 3 of 
c                       the CSDP, massurn3(maxint).
c        massurn4    :  real vector giving the mass for the urn 4 of 
c                       the CSDP, massurn4(maxint).
c        maxu        :  real working variable.
c        mrand       :  index.
c        nclusteru   :  integer giving the number of clusters in urns.
c        nclusteruc  :  integer giving the number of clusters in urns.
c        nclusw      :  index. 
c        npoints     :  index.
c        ns          :  index.
c        nscan       :  index.
c        nso         :  index.
c        numc        :  integer vector giving the number of v_i in 
c                       each urn, numc(4).
c        numc2       :  integer vector giving the number of v_i in 
c                       each urn for the candidate thetac, 
c                       numc2(4).
c        ok          :  integer indicator.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(maxint).
c        proburn1    :  real vector giving F_1 in the CSDP, 
c                       proburn1(maxint).
c        ratio       :  real working variable.
c        rbeta       :  real beta random number generator.
c        rtslogistic2:  real truncated logistic random number generator.
c        runif       :  real uniform random number generator.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        skipcount   :  index. 
c        sprint      :  integer function to print information on the
c                       screen.
c        thetac      :  real giving the value of the candidate of theta.
c        thetaskip   :  index
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        urn         :  integer vector indicating the urn
c                       where each interval belongs,
c                       urn(maxint).
c        uvec        :  real working vector, uvec(maxm).
c        vpred       :  real working variable. 
c        vnew        :  real working vector, vnew(nrec+1).
c        vvec        :  real working vector, vvec(maxm).
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
c        wvec        :  real working vector, wvec(maxm).
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
      double precision a0b0(4),aa0,ab0
      double precision betapm(p),betapv(p,p),d,ppar

c+++++MCMC parameters
      integer extra,mcmcvec(5),nburn,nskip,nsave,ntheta,ndisplay
      double precision propv(p,p)

c+++++Stored output
      double precision acrate(2),fsave(nsave,nlink)
      double precision randsave(nsave,nrec+1),thetasave(nsave,p+3)
      double precision cpo(nrec)

c+++++Current values of the parameters
      integer ncluster,y(nrec)
      double precision alpha,beta(p),theta,v(nrec)

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c+++++Working space - integers
      integer maxm,maxend,maxint
      integer clusts(nrec,nrec+1)
      integer counter,dispcount 
      integer efind
      integer evali
      integer i 
      integer index(maxm)
      integer iflag(p)
      integer imax,imin
      integer indi 
      integer intcount(maxint),intcount2(maxint)
      integer intcounturn(4)
      integer intposso(maxint,nrec+1),intpossn(maxint,nrec+1)
      integer intind(nrec+1),intind2(nrec+1)
      integer isave,iscan
      integer j
      integer k
      integer lpsav(nrec),lpsavc(nrec)
      integer mrand
      integer nclusteru(4)
      integer nclusteruc(4)
      integer npoints
      integer ns,nscan,nso
      integer ok
      integer s(nrec)
      integer seed(2),seed1,seed2
      integer skipcount
      integer sprint
      integer thetaskip
      integer urn(maxint)

c+++++Working space - double precision
      double precision area,ac(4),ac2(4) 
      double precision betac(p),cdflogis 
      double precision dgamlog
      double precision endp(maxend),endp2(maxend)
      double precision eta(nrec),etan(nrec)
      double precision keepbeta
      double precision liminf,limsup
      double precision limr(nrec,2)
      double precision logliko,loglikn,logprioro,logpriorn
      double precision mass(maxint)
      double precision massurn1(maxint),massurn2(maxint)
      double precision massurn3(maxint),massurn4(maxint)
      double precision maxu
      double precision nclusw
      double precision numc(4),numc2(4)
      double precision prob(maxint)
      double precision proburn1(maxint)
      double precision ratio
      double precision rbeta,rtslogistic2
      double precision thetac
      double precision tmp1,tmp2
      double precision uvec(maxm)
      double precision vpred
      double precision vnew(nrec+1)
      double precision vvec(maxm)
      double precision workm1(p,p),workm2(p,p)
      double precision workmh1(p*(p+1)/2)
      double precision workv1(p),workv2(p)
      double precision wvec(maxm)

c+++++Working space - single precision
      real runif

c+++++Working space - logical
      logical ainf,binf

c++++ initialize variables

      nburn=mcmcvec(1)
      nskip=mcmcvec(2)
      ndisplay=mcmcvec(3)
      ntheta=mcmcvec(4)
      model=mcmcvec(5)
      
      aa0=a0b0(1)
      ab0=a0b0(2)
      d=a0b0(3)
      ppar=a0b0(4)

c++++ set random number generator

      seed1=seed(1)
      seed2=seed(2)

      call setall(seed1,seed2)

      call rdisc(1,nrec,evali)
      vpred=v(evali)
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Check consistency with the data before starting 
c++++ the algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c++++ check if the user has requested an interrupt
      call rchkusr()
 
      do i=1,nrec
         
         if(model.eq.0)then
            sens(i)=1.d0
            spec(i)=1.d0
         end if
      
         limr(i,1)=-999.d0         
         limr(i,2)= 999.d0
      end do
 
      nclusw=1 
      do i=1,nrec
         if(i.eq.1)then
           s(i)=1
          else
           j=0 
           ok=0
           do while(ok.eq.0.and.j.lt.(i-1))          
              j=j+1
              if(v(i).eq.v(j))then
                ok=1
                s(i)=s(j)
              end if
           end do
           if(ok.eq.0)then
             nclusw=nclusw+1
             s(i)=nclusw
           end if  
         end if
         
         if(y(i).eq.1)then
            if(eta(i).lt.limr(s(i),2))limr(s(i),2)=eta(i)          
          else
            if(eta(i).gt.limr(s(i),1))limr(s(i),1)=eta(i)
         end if 

         clusts(s(i),1)=clusts(s(i),1)+1
         clusts(s(i),clusts(s(i),1)+1)=i      
      
         vnew(i)=v(i) 
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
              call intpr("int",-1,intind(i),1)
              call dblepr("eta",-1,eta(i),1)
              call dblepr("v",-1,v(i),1)
              call rexit("Errors not consistent with data in S0")
           end if   
         end if
      end do

      ncluster=nclusw

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      thetaskip=0
      nscan=nburn+(nskip+1)*(nsave)

      call cpu_time(sec0)
      sec00=0.d0


c++++ evaluate log-prior for current value of parameters

      call dmvn(p,beta,betapm,betapv,logprioro,workv1,workm1,
     &          workm2,workv2,iflag)  


      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to updating regression coefficients and v +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ sample candidates

         call rmvnorm(p,beta,propv,workmh1,workv1,betac)

c+++++++ evaluate log-prior for candidate value of parameters

         call dmvn(p,betac,betapm,betapv,logpriorn,workv1,workm1,
     &             workm2,workv2,iflag)  

c+++++++ define the end points
         
         counter=0
         do i=1,nrec
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+x(i,j)*betac(j)                   
            end do
            etan(i)=tmp1
            counter=counter+1
            endp(counter)=eta(i)
            counter=counter+1
            endp(counter)=etan(i)
         end do

         counter=counter+1
         endp(counter)=0.d0
         counter=counter+1
         endp(counter)=theta
         counter=counter+1
         endp(counter)=theta-d

         do i=1,nlink
            counter=counter+1
            endp(counter)=xlink(i)
         end do


c+++++++ sort end points

         call sortvec(maxend,endp,1,counter)
         
         npoints=1
         endp2(1)=endp(1)
         do i=2,counter
            if(endp(i).ne.endp2(npoints))then
               npoints=npoints+1
               endp2(npoints)=endp(i)
            end if
         end do

         if(maxint.lt.(npoints+1))then
            call rexit("Error in the maximum number of intervals")
         end if   

c+++++++ updating G base on the partition defined by beta and beta*

         do i=1,maxint
            mass(i)=0.d0
            massurn1(i)=0.d0
            massurn2(i)=0.d0
            massurn3(i)=0.d0
            massurn4(i)=0.d0
            prob(i)=0.d0
            proburn1(i)=0.d0
            intcount(i)=0
            intcount2(i)=0
            urn(i)=0
         end do

         do i=1,4
            intcounturn(i)=0
         end do

         do i=1,npoints
            evali=4 
            if(endp2(i).le. (theta)  )evali=3 
            if(endp2(i).le. 0.d0     )evali=2 
            if(endp2(i).le. (theta-d))evali=1 
            urn(i)=evali
            intcounturn(evali)=intcounturn(evali)+1  
         end do

         urn(npoints+1)=4
         intcounturn(4)=intcounturn(4)+1  

         do i=1,nrec
         
            imax=efind(eta(i),maxend,endp2,npoints)
            lpsav(i)=imax
         
            indi=npoints+1
            ok=0
            j=npoints
            do while(j.ge.1.and.ok.eq.0)
               if(v(i).le.endp2(j))then
                  indi=indi-1
                  j=j-1
                else
                  ok=1
               end if  
            end do
            
            mass(indi)=mass(indi)+1.d0
            intind2(i)=indi
            intcount2(indi)=intcount2(indi)+1
            
            intposso(indi,intcount2(indi))=i
         end do

         ac(1)=alpha*cdflogis(theta-d,0.d0,1.d0,1,0)
         ac(2)=alpha*(0.5d0-cdflogis(theta-d,0.d0,1.d0,1,0))
         ac(3)=alpha*(cdflogis(theta,0.d0,1.d0,1,0)-0.5d0)
         ac(4)=alpha*(1.d0-cdflogis(theta,0.d0,1.d0,1,0))
         
         do i=1,npoints+1
            area=0.d0
            if(i.eq.1)then
                area=alpha*cdflogis(endp2(1),0.d0,1.d0,1,0)
              else if(i.gt.npoints)then
                area=alpha*cdflogis(endp2(npoints),0.d0,1.d0,0,0)
              else
                area=alpha*(cdflogis(endp2(i  ),0.d0,1.d0,1,0)-  
     &                      cdflogis(endp2(i-1),0.d0,1.d0,1,0))
            end if  
            mass(i)=mass(i)+area
         end do
         
         do i=1,intcounturn(1) 
            massurn1(i)=mass(i)
         end do
          
         counter=0 
         do i=intcounturn(1)+1,intcounturn(1)+intcounturn(2) 
            counter=counter+1
            massurn2(counter)=mass(i)
         end do

         counter=0 
         do i=intcounturn(1)+intcounturn(2)+1,intcounturn(1)+
     &        intcounturn(2)+intcounturn(3) 
            counter=counter+1
            massurn3(counter)=mass(i)
         end do

         counter=0 
         do i=intcounturn(1)+intcounturn(2)+intcounturn(3)+1,
     &        intcounturn(1)+intcounturn(2)+intcounturn(3)+
     &        intcounturn(4)  
            counter=counter+1
            massurn4(counter)=mass(i)
         end do


         tmp1=0.d0
c+++++++ sampling G1

         if(intcounturn(1).gt.1)then
            call dirichlet(massurn1,maxint,intcounturn(1),proburn1)
           else
            proburn1(1)=1.d0
         end if   

         do i=1,intcounturn(1) 
            prob(i)=exp(log(proburn1(i))+log(1.d0-ppar)+log(0.5d0))
            tmp1=tmp1+prob(i)
         end do

c+++++++ sampling G2

         if(intcounturn(2).gt.1)then
            call dirichlet(massurn2,maxint,intcounturn(2),proburn1)
           else
            proburn1(1)=1.d0
         end if   

         counter=0 
         do i=intcounturn(1)+1,intcounturn(1)+intcounturn(2) 
            counter=counter+1
            prob(i)=exp(log(proburn1(counter))+log(ppar)+log(0.5d0))
            tmp1=tmp1+prob(i)
         end do

c+++++++ sampling G3

         if(intcounturn(3).gt.1)then
            call dirichlet(massurn3,maxint,intcounturn(3),proburn1)
           else
            proburn1(1)=1.d0
         end if   

         counter=0 
         do i=intcounturn(1)+intcounturn(2)+1,intcounturn(1)+
     &        intcounturn(2)+intcounturn(3) 
            counter=counter+1
            prob(i)=exp(log(proburn1(counter))+log(ppar)+log(0.5d0))
            tmp1=tmp1+prob(i)
         end do

c+++++++ sampling G4

         if(intcounturn(4).gt.1)then
            call dirichlet(massurn4,maxint,intcounturn(4),proburn1)
           else
            proburn1(1)=1.d0
         end if   

         counter=0 
         do i=intcounturn(1)+intcounturn(2)+intcounturn(3)+1,
     &        intcounturn(1)+intcounturn(2)+intcounturn(3)+
     &        intcounturn(4)
            counter=counter+1
            prob(i)=exp(log(proburn1(counter))+
     &                  log(1.d0-ppar)+log(0.5d0))
            tmp1=tmp1+prob(i)
         end do

c+++++++ standarize G         
c         do i=1,npoints+1
c            prob(i)=prob(i)/tmp1
c         end do


c+++++++ evaluate log-likelihood for current and candidate value 
c+++++++ of parameters

        
         logliko=0.d0
         loglikn=0.d0

         do i=1,nrec

            imin=1
            imax=lpsav(i)
            tmp1=0.d0
            do j=imin,imax
               tmp1=tmp1+prob(j)
            end do
         
            imin=1
            imax=efind(etan(i),maxend,endp2,npoints)
            tmp2=0.d0
            do j=imin,imax
               tmp2=tmp2+prob(j)
            end do
            lpsavc(i)=imax

            if(yobs(i).eq.1)then
               tmp1=tmp1*sens(i)+(1.d0-spec(i))*(1.d0-tmp1) 
               if(tmp1.lt.zero)go to 100
               if(tmp1.gt.one )go to 100
               tmp2=tmp2*sens(i)+(1.d0-spec(i))*(1.d0-tmp2) 
               if(tmp2.lt.zero)go to 100
               if(tmp2.gt.one )go to 100
               logliko=logliko+log(tmp1)
               loglikn=loglikn+log(tmp2)
            else
               tmp1=tmp1*sens(i)+(1.d0-spec(i))*(1.d0-tmp1) 
               if(tmp1.lt.zero)go to 100
               if(tmp1.gt.one )go to 100
               tmp2=tmp2*sens(i)+(1.d0-spec(i))*(1.d0-tmp2) 
               if(tmp2.lt.zero)go to 100
               if(tmp2.gt.one )go to 100
               logliko=logliko+log(1.d0-tmp1)
               loglikn=loglikn+log(1.d0-tmp2)
            end if   

         end do

c+++++++ acceptance step

         ok=0
         ratio=dexp(loglikn+logpriorn-logliko-logprioro)

         if(dble(runif()).lt.ratio)then
            do j=1,p
               beta(j)=betac(j)
            end do
            logprioro=logpriorn
            do i=1,nrec
               eta(i)=etan(i)
               lpsav(i)=lpsavc(i)
            end do
            acrate(1)=acrate(1)+1.d0
           else 
            ok=1            
         end if

         if(ok.eq.1)go to 100 

           if(model.eq.1)then
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++ updating the true binary variables               +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             do i=1,nrec
                
                imin=1
                imax=lpsav(i)
                tmp1=0.d0
                do j=imin,imax
                   tmp1=tmp1+prob(j)
                end do
                
                if(yobs(i).eq.1)then
                   tmp2=sens(i)*tmp1/
     &               (sens(i)*tmp1+(1.d0-spec(i))*(1.d0-tmp1))
                 else 
                   tmp2=(1.d0-sens(i))*tmp1/
     &               ((1.d0-sens(i))*tmp1+spec(i)*(1.d0-tmp1))
                end if

                call rbinom(1,tmp2,evali)
                
                if(evali.ne.0.and.evali.ne.1)then
                  call dblepr("prob",-1,tmp2,1) 
                  call intpr("evali",-1,evali,1) 
                  call rexit("Error in the generation of y")
                end if
                y(i)=evali
             end do
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++ end updating the true binary variables           +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           end if


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ updating errors using the Sethuraman-Tiwari (1982) +++
c+++++++ representation                                     +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c         do i=1,maxint
c            intcount(i)=0
c         end do
         
         do i=1,nrec

            limr(i,1)=-999.d0         
            limr(i,2)= 999.d0
            clusts(i,1)=0

            imin=1
            imax=lpsav(i)
            tmp1=0.d0
            do j=imin,imax
               tmp1=tmp1+prob(j)
            end do

            if(y(i).eq.0)then
              imin=imax+1
              imax=npoints+1
              tmp1=1.d0-tmp1
            end if
            
            if(tmp1.eq.0.d0)then
               call intpr("Imin",-1,imin,1)
               call intpr("Imax",-1,imax,1)
               call rexit("Zero Probability")
              else
               call simdiscint(prob,maxint,imin,imax,evali)
            end if   

            if(evali.lt.imin)then
               call intpr("Evali",-1,evali,1)
               call intpr("Imin",-1,imin,1)
               call intpr("Imax",-1,imax,1)
               call rexit("Error in the intervals")
            end if  
            if(evali.gt.imax)then
               call intpr("Evali",-1,evali,1)
               call intpr("Imin",-1,imin,1)
               call intpr("Imax",-1,imax,1)
               call rexit("Error in the intervals")
            end if  
            
            if(evali.lt.1.and.evali.gt.(npoints+1))then
               call intpr("Evali",-1,evali,1)
               call intpr("Imin",-1,imin,1)
               call intpr("Imax",-1,imax,1)
               call rexit("Error in the intervals")
            end if  

            intcount(evali)=intcount(evali)+1
            intind(i)=evali

            intpossn(evali,intcount(evali))=i

         end do

c+++++++ start predictions         

         call simdiscint(prob,maxint,1,npoints+1,evali)
         intcount(evali)=intcount(evali)+1
         intind(nrec+1)=evali
         intpossn(evali,intcount(evali))=nrec+1
         
c+++++++ end predictions         

         do i=1,npoints+1
            ns=intcount(i)
            nso=intcount2(i)
            area=0.d0
            if(i.eq.1)then
               area=alpha*cdflogis(endp2(1),0.d0,1.d0,1,0)
               ainf=.true.
               binf=.false.
               liminf=0.d0
               limsup=endp2(1)
             else if(i.gt.npoints)then
               area=alpha*cdflogis(endp2(npoints),0.d0,1.d0,0,0)
               ainf=.false.
               binf=.true.
               liminf=endp2(npoints)
               limsup=0.d0
             else
               area=alpha*(cdflogis(endp2(i  ),0.d0,1.d0,1,0)-  
     &                     cdflogis(endp2(i-1),0.d0,1.d0,1,0))
               ainf=.false.
               binf=.false.
               liminf=endp2(i-1)
               limsup=endp2(i)
            end if
            
            if(ns.gt.0)then
               maxu=0.d0
               do j=1,ns
                  uvec(j)=dble(runif())
                  if(uvec(j).gt.maxu)maxu=uvec(j)
               end do

               mrand=1
               keepbeta=rbeta(1.d0,mass(i))
               wvec(mrand)=keepbeta
               tmp1=1.d0-keepbeta
               tmp2=keepbeta
               
               do while(maxu.ge.tmp2.and.mrand.lt.maxm)
                  mrand=mrand+1
                  if(mrand.gt.maxm)then
                    call rexit("Increase maxn in R function")
                  end if  
                  keepbeta=rbeta(1.d0,mass(i))
                  wvec(mrand)=keepbeta*tmp1
                  tmp1=tmp1*(1.d0-keepbeta)
                  tmp2=tmp2+wvec(mrand)
               end do

               do j=1,mrand
                  index(j)=j 
                  if(dble(runif()).le.(1.d0-(area/mass(i))))then
                     call rdisc(1,nso,evali)
                     if(evali.lt.1.or.evali.gt.nso)return
                     vvec(j)=v(intposso(i,evali))
                   else
                     tmp1=rtslogistic2(ainf,binf,liminf,limsup)
                     vvec(j)=tmp1
                  end if
               end do
               
               if(mrand.gt.1)then  
                  call sortvecind(maxm,index,vvec,wvec,1,mrand)
               end if   

               do j=1,ns
                  tmp1=0.d0
                  k=1
                  do while(uvec(j).gt.tmp1.and.k.lt.mrand)
                     tmp1=tmp1+wvec(k)
                     k=k+1  
                  end do
                  vnew(intpossn(i,j))=vvec(k)
               end do
               
            end if 
         end do

         nclusw=1
         
         do i=1,nrec
            v(i)=vnew(i)

            if(i.eq.1)then
              s(i)=1
             else
              j=0 
              ok=0
              do while(ok.eq.0.and.j.lt.(i-1))          
                 j=j+1
                 if(v(i).eq.v(j))then
                   ok=1
                   s(i)=s(j)
                 end if
              end do
              if(ok.eq.0)then
                nclusw=nclusw+1
                s(i)=nclusw
              end if  
            end if     
            
            if(y(i).eq.1)then
               if(eta(i).lt.limr(s(i),2))limr(s(i),2)=eta(i)          
             else
               if(eta(i).gt.limr(s(i),1))limr(s(i),1)=eta(i)
            end if 

            clusts(s(i),1)=clusts(s(i),1)+1
            clusts(s(i),clusts(s(i),1)+1)=i                    
         end do
         
         ncluster=nclusw         
         vpred=vnew(nrec+1)

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Check consistency with the data    
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
c         call rchkusr()
c         do i=1,nrec
c            vnew(i)=v(i) 
c            if(y(i).eq.1)then
c              if(v(i).gt.eta(i))then
c                 call intpr("i",-1,i,1)
c                 call intpr("y",-1,y(i),1)
c                 call dblepr("eta",-1,eta(i),1)
c                 call dblepr("v",-1,v(i),1)
c                 call rexit("Errors not consistent with data in S1")
c              end if   
c            end if
c            if(y(i).eq.0)then
c              if(v(i).le.eta(i))then
c                 call intpr("i",-1,i,1)
c                 call intpr("y",-1,y(i),1)
c                 call intpr("int",-1,intind(i),1)
c                 call dblepr("eta",-1,eta(i),1)
c                 call dblepr("v",-1,v(i),1)
c                 call rexit("Errors not consistent with data in S2")
c              end if   
c            end if
c         end do

100      continue

c+++++++++++++++++++++++++++++++++++++++++++++         
c+++++++ Extra Step which moves the clusters
c+++++++++++++++++++++++++++++++++++++++++++++

        if(extra.eq.1)then
           do i=1,ncluster
              ainf=.false.
              binf=.false.
              liminf=limr(i,1)
              limsup=limr(i,2)         
              tmp1=rtslogistic2(ainf,binf,liminf,limsup)
              do j=1,clusts(i,1)
                 v(clusts(i,j+1))=tmp1
              end do
           end do
        end if

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Some computations necessary to update theta
c+++++++ are performed here to increase the efficiency.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         thetac=dble(runif())*d
         
         do i=1,4
            numc(i)=0.d0
            numc2(i)=0.d0
            nclusteru(i)=0
            nclusteruc(i)=0
         end do

         
         do i=1,ncluster
            if(v(clusts(i,2)).le.(theta-d))then
               nclusteru(1)=nclusteru(1)+1
               numc(1)=numc(1)+clusts(i,1)
             else if (v(clusts(i,2)).le.0.d0)then 
               nclusteru(2)=nclusteru(2)+1
               numc(2)=numc(2)+clusts(i,1)
             else if (v(clusts(i,2)).le.theta)then   
               nclusteru(3)=nclusteru(3)+1
               numc(3)=numc(3)+clusts(i,1)               
             else 
               nclusteru(4)=nclusteru(4)+1
               numc(4)=numc(4)+clusts(i,1)               
            end if   

            if(v(clusts(i,2)).le.(thetac-d))then
               nclusteruc(1)=nclusteruc(1)+1
               numc2(1)=numc2(1)+clusts(i,1)
             else if (v(clusts(i,2)).le.0.d0)then 
               nclusteruc(2)=nclusteruc(2)+1
               numc2(2)=numc2(2)+clusts(i,1)
             else if (v(clusts(i,2)).le.thetac)then   
               nclusteruc(3)=nclusteruc(3)+1
               numc2(3)=numc2(3)+clusts(i,1)               
             else 
               nclusteruc(4)=nclusteruc(4)+1
               numc2(4)=numc2(4)+clusts(i,1)               
            end if   
         end do

c         call intpr("nclus",-1,nclusteru,4)
c         call dblepr("num",-1,numc,4)

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ updating theta using an independent sampler        +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         thetaskip = thetaskip + 1
         if(thetaskip.ge.ntheta)then

            ac(1)=alpha*cdflogis(theta-d,0.d0,1.d0,1,0)
            ac(2)=alpha*(0.5d0-cdflogis(theta-d,0.d0,1.d0,1,0))
            ac(3)=alpha*(cdflogis(theta,0.d0,1.d0,1,0)-0.5d0)
            ac(4)=alpha*(1.d0-cdflogis(theta,0.d0,1.d0,1,0))

            ac2(1)=alpha*cdflogis(thetac-d,0.d0,1.d0,1,0)
            ac2(2)=alpha*(0.5d0-cdflogis(thetac-d,0.d0,1.d0,1,0))
            ac2(3)=alpha*(cdflogis(thetac,0.d0,1.d0,1,0)-0.5d0)
            ac2(4)=alpha*(1.d0-cdflogis(thetac,0.d0,1.d0,1,0))

            ratio=(numc2(2)+numc2(3)-numc(2)-numc(3))*log(ppar)
            ratio=ratio+(numc2(1)+numc2(4)-numc(1)-numc(4))*
     &                  log(1.d0-ppar) 

            do j=1,4
               ratio=ratio+dgamlog(dble(numc(j) +ac(j)))
     &                    -dgamlog(dble(numc2(j)+ac2(j)))
               ratio=ratio+dgamlog(ac2(j))-dgamlog(ac(j))
            end do

            if(log(runif()).lt.ratio)then
               theta=thetac
               acrate(2)=acrate(2)+1.d0
               do i=1,4
                  ac(i)=ac2(i)
                  numc(i)=numc2(i)
               end do
               
            end if
            thetaskip=0
         end if
             
c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++
         if(aa0.gt.0.d0)then
             call samalphcsdp(alpha,aa0,ab0,ncluster,numc,ac,theta)
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

c+++++++++++++ theta
               thetasave(isave,p+1)=theta

c+++++++++++++ cluster information
               thetasave(isave,p+2)=ncluster
               thetasave(isave,p+3)=alpha
               
               
c+++++++++++++ link information,
c+++++++++++++ cpo, errors and predictive information

               do i=1,nrec

                  randsave(isave,i)=v(i)

                  imin=1
                  imax=lpsav(i)
                  tmp1=0.d0
                  do j=imin,imax
                     tmp1=tmp1+prob(j)
                  end do
                  
                  tmp2=sens(i)*tmp1+(1.d0-spec(i))*(1.d0-tmp1)

                  if(yobs(i).eq.1)then
                    cpo(i)=cpo(i)+1.d0/tmp2 
                  else
                    cpo(i)=cpo(i)+1.d0/(1.d0-tmp2)
                  end if                          
               
               end do
               
               randsave(isave,nrec+1)=vpred

               do i=1,nlink
                  imin=1
                  imax=efind(xlink(i),maxend,endp2,npoints)
                  tmp1=0.d0
                  do j=imin,imax
                    tmp1=tmp1+prob(j)
                  end do
                  fsave(isave,i)=tmp1
               end do
               
c               if(isave.eq.3602)then
c                 call dblepr("prob",-1,prob,npoints+1)
c               end if  

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

      acrate(1)=acrate(1)/dble(nscan)      
      acrate(2)=acrate(2)*dble(ntheta)/dble(nscan)      

      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do

      return
      end


c=======================================================================                  
      subroutine csdppredl(nsave,nrec,npred,alpha,theta,v,d,ppar,lp,
     &                     out)
c=======================================================================                  
c     Draws samples from F(t), where F follows a CSDP prior.
c     This function uses the fact that F(t) follows a modified beta 
c     distribution under the CSDP prior.
c     A.J.V., 2006
c=======================================================================                  
      implicit none
      integer nsave,nrec,npred
      double precision alpha(nsave),theta(nsave),v(nsave,nrec),d,ppar
      double precision lp(npred,nsave),out(npred,nsave)
      integer i,j,k,count,urn,ntotal
      double precision area,mass
      double precision tmp1,tmp2,tmp3,tmp4,x
      double precision cdflogis,rbeta

      do i=1,npred
         do j=1,nsave

c++++++++++ check if the user has requested an interrupt
            call rchkusr()
    
            x=lp(i,j)
            if(x.le.theta(j)-d)then
               mass=alpha(j)*cdflogis(theta(j)-d,0.d0,1.d0,1,0)
               area=alpha(j)*cdflogis(x,0.d0,1.d0,1,0)
               urn=1
             else if(x.le.0.d0)then
               mass=alpha(j)*(0.5d0-cdflogis(theta(j)-d,0.d0,1.d0,1,0))
               area=alpha(j)*(cdflogis(x,0.d0,1.d0,1,0)-
     &                        cdflogis(theta(j)-d,0.d0,1.d0,1,0))
               urn=2
             else if(x.le.theta(j))then   
               mass=alpha(j)*(cdflogis(theta(j),0.d0,1.d0,1,0)-0.5d0)
               area=alpha(j)*(cdflogis(x,0.d0,1.d0,1,0)-0.5d0)
               urn=3
             else  
               mass=alpha(j)*(1.d0-cdflogis(theta(j),0.d0,1.d0,1,0))
               area=alpha(j)*(cdflogis(x,0.d0,1.d0,1,0)-
     &                        cdflogis(theta(j),0.d0,1.d0,1,0))
               urn=4
            end if   
            
            count=0
            ntotal=0
            do k=1,nrec
               if(urn.eq.1)then
                 if(v(j,k).le.(theta(j)-d))ntotal=ntotal+1         
                 if(v(j,k).le.lp(i,j))count=count+1
                else if(urn.eq.2)then
                 if(v(j,k).gt.(theta(j)-d).and.v(j,k).le.0.d0)then
                    ntotal=ntotal+1      
                 end if   
                 if(v(j,k).gt.(theta(j)-d).and.v(j,k).le.lp(i,j))then
                    count=count+1
                 end if   
                else if(urn.eq.3)then
                 if(v(j,k).gt.0.d0.and.v(j,k).le.theta(j))then
                    ntotal=ntotal+1      
                 end if   
                 if(v(j,k).gt.0.d0.and.v(j,k).le.lp(i,j))then
                    count=count+1
                 end if   
                else 
                 if(v(j,k).gt.theta(j))then
                    ntotal=ntotal+1      
                 end if   
                 if(v(j,k).gt.theta(j).and.v(j,k).le.lp(i,j))then
                    count=count+1
                 end if   
               end if                
            end do

            tmp1=(area+dble(count))
            tmp2=(mass-area+dble(ntotal-count))
            tmp3=rbeta(tmp1,tmp2)
            
            if(urn.eq.1)then
              tmp4=0.5d0*(1.d0-ppar)*tmp3
             else if(urn.eq.2)then 
              tmp4=0.5d0*(1.d0-ppar)+0.5d0*ppar*tmp3
             else if(urn.eq.3)then 
              tmp4=0.5d0+0.5d0*ppar*tmp3
             else 
              tmp4=0.5d0*(1.d0+ppar)+0.5d0*(1.d0-ppar)*tmp3
            end if 
            out(i,j)=tmp4
            
         end do
      end do
      
      return
      end
      

c=======================================================================                  
      subroutine samalphcsdp(alpha,aa0,ab0,ncluster,numc,ac,theta)
c=======================================================================                  
c     Draws a sample for the precision parameter of a CSDP prior
c     when a Gamma(aa0,ab0) prior is used for it.
c     The strategy is similar to the one of Escobar and West (1995)
c     for DP priors.
c     A.J.V., 2006
c=======================================================================                  
      implicit none
      integer i,evali
      integer ncluster
      double precision ac(4),numc(4),theta
      double precision alpha,aa0,ab0
      double precision ac2(4)
      double precision dgamlog,rbeta,rgamma,eta(4),ww(5)
      double precision tmp1,tmp2,tmp3
      double precision total

      do i=1,5
         ww(i)=0.d0
      end do
      
      tmp1=1.d0
      tmp2=1.d0
      tmp3=0.d0
      do i=1,4
         ac2(i)=ac(i)/alpha
         eta(i)=rbeta(ac(i)+1.d0,numc(i))
         tmp1=tmp1*ac2(i)
         tmp2=tmp2*numc(i)
         tmp3=tmp3+ac2(i)*log(eta(i))
      end do  
      
      ww(1)=tmp1

      ww(2)=       numc(1)*ac2(2)*ac2(3)*ac2(4)
      ww(2)=ww(2)+ numc(2)*ac2(1)*ac2(3)*ac2(4)
      ww(2)=ww(2)+ numc(3)*ac2(1)*ac2(2)*ac2(4)
      ww(2)=ww(2)+ numc(4)*ac2(1)*ac2(2)*ac2(3)
      
      ww(3)=       numc(1)*numc(2)*ac2(3)*ac2(4)
      ww(3)=ww(3)+ numc(1)*numc(3)*ac2(2)*ac2(4)
      ww(3)=ww(3)+ numc(1)*numc(4)*ac2(2)*ac2(3)
      ww(3)=ww(3)+ numc(2)*numc(3)*ac2(1)*ac2(4)
      ww(3)=ww(3)+ numc(2)*numc(4)*ac2(1)*ac2(3)
      ww(3)=ww(3)+ numc(3)*numc(4)*ac2(1)*ac2(2)

      ww(4)=       numc(1)*numc(2)*numc(3)*ac2(4)
      ww(4)=ww(4)+ numc(1)*numc(2)*numc(4)*ac2(3)
      ww(4)=ww(4)+ numc(1)*numc(3)*numc(4)*ac2(2)
      ww(4)=ww(4)+ numc(2)*numc(3)*numc(4)*ac2(1)

      ww(5)=tmp2

      total=0.d0
      
      
      ww(1)=ww(1)*exp(dgamlog(aa0+dble(ncluster))-
     &               (aa0+dble(ncluster))*log(theta))
  
      ww(2)=ww(2)*exp(dgamlog(aa0+dble(ncluster-1))-
     &               (aa0+dble(ncluster-1))*log(theta))

      ww(3)=ww(3)*exp(dgamlog(aa0+dble(ncluster-2))-
     &               (aa0+dble(ncluster-2))*log(theta))

      ww(4)=ww(4)*exp(dgamlog(aa0+dble(ncluster-3))-
     &               (aa0+dble(ncluster-3))*log(theta))

      ww(5)=ww(5)*exp(dgamlog(aa0+dble(ncluster-4))-
     &               (aa0+dble(ncluster-4))*log(theta))

  
      call simdiscint(ww,5,1,5,evali)               
  
      tmp3=ab0-tmp3
  
      if(evali.eq.1)then
         tmp1=aa0+dble(ncluster)
         alpha=rgamma(tmp1,tmp3)
       else if(evali.eq.2)then
         tmp1=aa0+dble(ncluster-1)
         alpha=rgamma(tmp1,tmp3)
       else if(evali.eq.3)then
         tmp1=aa0+dble(ncluster-2)
         alpha=rgamma(tmp1,tmp3)
       else if(evali.eq.4)then       
         tmp1=aa0+dble(ncluster-3)
         alpha=rgamma(tmp1,tmp3)
       else
         tmp1=aa0+dble(ncluster-4)
         alpha=rgamma(tmp1,tmp3)
      end if 
  
  
      return
      end

