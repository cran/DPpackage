    
c=======================================================================                      
      subroutine psgampois
     &                (nrec,nfixed,p,nsmooth,q,starts,ends,nsmooths,    
     &                 maxq,x,z,y,roffset,ngrid,znew,xreal,possfp,
     &                 bet0,prec,sb,taub,iapm,japm,apm,nznh,pordv,      
     &                 beta,b,sigmab,                                   
     &                 mcmc,nsave,                                      
     &                 acrate,cpo,randsave,thetasave,pssave,            
     &                 seed,                                            
     &                 iflagp,betac,workmp1,workmp2,workmhp1,workvp1,         
     &                 xtx,xty,                                         
     &                 iflagq,bc,workmq1,workmq2,workmhq1,workvq1,            
     &                 ztz,zty,ztzinv,theta,                            
     &                 mc,betasave,bsave,workvps)                       
c=======================================================================                      
c     # of arguments = 58.
c
c     Subroutine `psgampois' to run a Markov chain in a semiparametric 
c     poisson model, using a B-splines and penalties.
c
c     Copyright: Alejandro Jara, 2007-2010.
c
c     Version 1.0:
c
c     Last modification: 15-09-2010.
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
c        nrec        :  integer giving the number of observations.
c        nfixed      :  integer giving the number of parametric effects,
c                       if nfixed is 0 then p=1.
c        p           :  integer giving the number of parametric
c                       coefficients.
c        nsmooth     :  integer giving the number of smoothers.
c        nsmooths    :  integer giving the number of smoothers.
c        q           :  integer giving the dimension of the desing
c                       matrix in the mixed model representation of
c                       splines.
c        starts      :  integer vector with the starting locations
c                       for each smoother, starts(nsmooths)
c        ends        :  integer vector with the final locations
c                       for each smoother, starts(nsmooths)
c        maxq        :  integer used to define the maximum number
c                       of elements in the smooth representation.
c        x           :  real matrix giving the design matrix for the 
c                       parametric effects, x(nrec,p). 
c        y           :  integer vector giving the response variable,
c                       y(nrec).
c        z           :  real matrix giving the design matrix for the 
c                       smoothers, z(nrec,q). 
c        roffset     :  real vector of offsets, roffset(nrec).
c        ngrid       :  integer giving the number of grid points
c                       where the curves are evaluated. 
c        znew        :  real matrix giving the design for a grid
c                       of values of the covariates, znew(ngrid,nsmooths). 
c        xreal       :  real matrix giving the grid
c                       of values of the covariates, xreal(ngrid,nsmooths). 
c        possfp      :  integer vector giving the location of the 
c                       smoothers in the parametric part of the model,
c                       possfp(nsmooths).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        bet0        :  real vector giving the prior mean for the fixed 
c                       effects, bet0(p).
c        prec        :  real matrix giving the prior precision matrix
c                       for the fixed effects, prec(p,p).
c        sb          :  real vector giving the product of the prior 
c                       precision and the prior mean for the fixed 
c                       effects, sb(p).
c        taub1, taub2:  reals giving the hyperparameters of the prior 
c                       distribution for the inverse of the
c                       variances, 1/sigmab ~ Gamma(taub1/2,taub2/2).
c        iapm,japm,apm: linked list representation of the penalty matrix.
c        nznh         : integer giving the number of elements in the
c                       linked list representation of the penalty matrix.
c        pordv        : integer vector giving the penalty type for
c                       each smoother, pordv(nsmooths). 
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        beta        :  real vector giving the current value of the 
c                       fixed effects, beta(p).
c        b           :  real vector giving the random effects,
c                       b(q) .
c        sigmab      :  real vector giving the current value of the
c                       variances .
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
c
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        acrate      :  real vector giving the MH acceptance rate, 
c                       acrate(1+nsmooths). 
c        cpo         :  real giving the cpo's. 
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,q).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the fixed effects and variances, 
c                       thetsave(nsave,p+nsmooths).
c        pssave       : real matrix giving the curve and variances,
c                       pssave(ngrid*nsmooths,2).
c
c
c-----------------------------------------------------------------------
c
c
c---- Working space ----------------------------------------------------
c
c        bc          :  real vector giving the candidate value of the 
c                       smoothers, bc(maxq).
c        betac       :  real vector giving the candidate value of the 
c                       fixed effects, betac(p).
c        betasave    :  real vector used to save the posterior mean
c                       of fixed effects, betasave(p).
c        bsave       :  real vector used to save the posterior mean
c                       of random effects, bsave(q).
c        iflagp      :  integer vector used to invert the of the lhs
c                       least square solution for the fixed effects,
c                       iflagp(p).
c        iflagq      :  integer vector used to invert the of the lhs
c                       least square solution for the smoothers,
c                       iflagp(maxq).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index. 
c        l           :  index.
c        logcgko     :  real working variable.
c        logcgkn     :  real working variable.
c        loglikn     :  real working variable.
c        logliko     :  real working variable.
c        logpriorn   :  real working variable.
c        logprioro   :  real working variable.
c        mc          :  real working vactor used for model comparison,
c                       mc(5). 
c        ratio       :  real working variable.
c        sec         :  cpu time working variable.
c        sec0        :  cpu time working variable.
c        sec00       :  cpu time working variable.
c        sec1        :  cpu time working variable.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        skipcount   :  index. 
c        sprint      :  integer function to print on screen.
c        theta       :  real working vector, theta(maxq).
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        tmp3        :  real used to accumulate quantities.
c        workmp1     :  real matrix used to update the fixed effects,
c                       workvmp1(p,p).
c        workmp2     :  real matrix used to update the fixed effects,
c                       workvmp2(p,p).
c        workmq1     :  real matrix used to update the smoothers,
c                       workvmq1(maxq,maxq).
c        workmq2     :  real matrix used to update the smoothers,
c                       workvmq2(maxq,maxq).
c        workmhp1    :  real vector used to update the fixed effects,
c                       workmhp1(p*(p+1)/2).
c        workmhq1    :  real vector used to update the smoothers,
c                       workmhp1(maxq*(maxq+1)/2).
c        workvp1     :  real vector used to update the fixed effects,
c                       workvp1(p).
c        workvq1     :  real vector used to update the smoothers,
c                       workvq1(maxq).
c        workvps     :  real working vector, workvps(ngrid*nsmooths). 
c        xtx         :  real matrix giving the product X^tX, xtx(p,p).
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p).
c        ztz         :  real matrix giving the product Z^tZ, ztz(maxq,
c                       maxq).
c        zty         :  real vector used to save the product 
c                       Zt(Y-Xbeta), zty(q).
c
c=======================================================================                  
      implicit none 

c+++++Data
      integer nrec,nfixed,p,nsmooth,nsmooths,q,ngrid,maxq
      integer y(nrec)
      integer starts(nsmooths),ends(nsmooths)
      integer possfp(nsmooths)
      real*8 roffset(nrec),x(nrec,p),z(nrec,q),znew(ngrid,q)
      real*8 xreal(ngrid,nsmooths)

c+++++Prior 
      integer nznh,iapm(q+1),japm(nznh)
      real*8 bet0(p),prec(p,p),sb(p)
      real*8 taub(2),taub1,taub2
      real*8 apm(nznh)
      real*8 pordv(nsmooths)
      
c+++++Current values of the parameters
      real*8 beta(p)
      real*8 b(q)
      real*8 sigmab(nsmooths)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      real*8 acrate(1+nsmooth) 
      real*8 cpo(nrec,2)
      real*8 randsave(nsave,q)
      real*8 thetasave(nsave,p+nsmooth)
      real*8 pssave(ngrid*nsmooths,2)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++Fixed effects
      integer iflagp(p) 
      real*8 betac(p)
      real*8 workmp1(p,p),workmp2(p,p)
      real*8 workmhp1(p*(p+1)/2)
      real*8 workvp1(p)
      real*8 xtx(p,p),xty(p)

c+++++Smoothers
      integer iflagq(maxq) 
      real*8 bc(maxq)
      real*8 theta(maxq)
      real*8 workmq1(maxq,maxq),workmq2(maxq,maxq)
      real*8 workmhq1(maxq*(maxq+1)/2)
      real*8 workvq1(maxq)
      real*8 ztz(maxq,maxq),zty(maxq)
      real*8 ztzinv(maxq,maxq)
      real*8 workvps(ngrid*nsmooths)

c++++ model's performance
      real*8 mc(5)
      real*8 betasave(p),bsave(q)
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer i,ii,j,jj,l,ll
      integer sprint
      integer yij     
      real*8 detlog
      real*8 eta,gprime
      real*8 mean
      real*8 logcgkn,logcgko
      real*8 loglikn,logliko
      real*8 logpriorn,logprioro
      real*8 offset
      real*8 ratio
      real*8 tmp1,tmp2
      real*8 ytilde

c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++CPU time
      real*8 sec00,sec0,sec1,sec

c+++++RNG and distributions
      real*8 dpoiss,rgamma
      real runif

c++++ model�s performance
      real*8 dbarc,dbar,dhat,pd,lpml

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      
c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)

c++++ unzip
      taub1=taub(1)
      taub2=taub(2)

c++++ opening temporary file

      open(unit=1,file='temp_dppackage.txt',status='unknown',
     &     form='unformatted')

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      dbar=0.d0
      isave=0
      skipcount=0
      dispcount=0
      nburn=nburn

      nscan=nburn+(nskip+1)*nsave

      call cpu_time(sec0)
      sec00=0.d0
      
      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++
c+++++++ fixed effects        +++
c++++++++++++++++++++++++++++++++

         if(nfixed.gt.0)then
            do i=1,p
               do j=1,p
                  xtx(i,j)=prec(i,j)
               end do
               xty(i)=sb(i)
            end do

            logliko=0.d0

            do i=1,nrec
               eta=0.d0
               offset=0.d0
               mean=0.d0
               gprime=0.d0
               yij=y(i)
            
               do j=1,p
                  eta=eta+x(i,j)*beta(j)
               end do
               
               do j=1,q
                  eta=eta+z(i,j)*b(j)
                  offset=offset+z(i,j)*b(j)
               end do
               offset=offset+roffset(i)

               mean=exp(eta)
               gprime=exp(-eta)
               ytilde=eta+(dble(yij)-mean)*gprime-offset
 
               logliko=logliko+dpoiss(dble(yij),mean,1)
              
               do j=1,p
                  do l=1,p
                     xtx(j,l)=xtx(j,l)+x(i,j)*x(i,l)*mean
                  end do
                  xty(j)=xty(j)+x(i,j)*ytilde*mean
               end do
            
            end do

            call inverse(xtx,p,iflagp)      

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+xtx(i,j)*xty(j) 
               end do
               workvp1(i)=tmp1
            end do

            call rmvnorm(p,workvp1,xtx,workmhp1,xty,betac)

c++++++++++ evaluating the candidate generating kernel

            call dmvn2(p,betac,workvp1,xtx,logcgko,
     &                 xty,workmp1,workmp2,iflagp)                 
            
c++++++++++ evaluating the likelihood

            do i=1,p
               do j=1,p
                  xtx(i,j)=prec(i,j)
               end do
               xty(i)=sb(i)
            end do

            loglikn=0.d0

            do i=1,nrec
               eta=0.d0
               offset=0.d0
               mean=0.d0
               gprime=0.d0
               yij=y(i)

               do j=1,p
                  eta=eta+x(i,j)*betac(j)
               end do
               
               do j=1,q
                  eta=eta+z(i,j)*b(j)
                  offset=offset+z(i,j)*b(j)
               end do
               offset=offset+roffset(i)

               mean=exp(eta)
               gprime=exp(-eta)
               ytilde=eta+(dble(yij)-mean)*gprime-offset

               loglikn=loglikn+dpoiss(dble(yij),mean,1)
               
               do j=1,p
                  do l=1,p
                     xtx(j,l)=xtx(j,l)+x(i,j)*x(i,l)*mean
                  end do
                  xty(j)=xty(j)+x(i,j)*ytilde*mean
               end do
            end do

            call inverse(xtx,p,iflagp)      

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+xtx(i,j)*xty(j) 
               end do
               workvp1(i)=tmp1
            end do

c++++++++++ evaluating the candidate generating kernel

            call dmvn2(p,beta,workvp1,xtx,logcgkn,
     &                 xty,workmp1,workmp2,iflagp)                 

c++++++++++ prior ratio
            logprioro=0.d0
            logpriorn=0.d0
         
            do i=1,p
               do j=1,p
                  logpriorn=logpriorn+(betac(i)-bet0(i))* 
     &                       prec(i,j)       *
     &                      (betac(j)-bet0(j))

                  logprioro=logprioro+(beta(i)-bet0(i))* 
     &                      prec(i,j)      *
     &                      (beta(j)-bet0(j))

               end do
            end do
      
            logpriorn=-0.5d0*logpriorn
            logprioro=-0.5d0*logprioro

c++++++++++ mh step

            ratio=loglikn-logliko+logcgkn-logcgko+
     &            logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               acrate(1)=acrate(1)+1.d0
               do i=1,p
                  beta(i)=betac(i) 
               end do
            end if
         end if          


c+++++++++++++++++++++++++++++++++
c+++++++ smoothers             +++ 
c+++++++++++++++++++++++++++++++++

        if(nsmooth.gt.0)then

           do ii=1,nsmooth
              jj=ends(ii)-starts(ii)+1
              
              do i=1,jj
                 theta(i)=b(starts(ii)+i-1)
                 do j=1,jj
                    ztz(i,j)=0.d0
                 end do
                 zty(i)=0.d0  
              end do
              
              do i=starts(ii),ends(ii)
                 do j=iapm(i),iapm(i+1)-1
                    ztz(i-starts(ii)+1,japm(j)-starts(ii)+1)=apm(j)*
     &                   (1.d0/sigmab(ii))                       
                 end do
              end do

              logliko=0.d0

              do i=1,nrec
                 eta=0.d0
                 offset=0.d0
                 mean=0.d0
                 gprime=0.d0
                 yij=y(i)
            
                 do j=1,p
                    eta=eta+x(i,j)*beta(j)
                    offset=offset+x(i,j)*beta(j)
                 end do

                 do j=1,q
                    eta=eta+z(i,j)*b(j)
                    if(j.lt.starts(ii).or.j.gt.ends(ii))then
                       offset=offset+z(i,j)*b(j)
                    end if   
                 end do
               
                 offset=offset+roffset(i)

                 mean=exp(eta)
                 gprime=exp(-eta)
                 ytilde=eta+(dble(yij)-mean)*gprime-offset
                 
                 logliko=logliko+dpoiss(dble(yij),mean,1)
               
                 do j=1,jj
                    do l=1,jj
                       ztz(j,l)=ztz(j,l)+z(i,starts(ii)+j-1)*
     &                          z(i,starts(ii)+l-1)*mean
                    end do
                    zty(j)=zty(j)+z(i,starts(ii)+j-1)*ytilde*mean
                 end do
              end do

              call invdet2(ztz,maxq,jj,ztzinv,detlog,iflagq,workvq1)

              do i=1,jj
                 tmp1=0.d0
                 do j=1,jj
                    tmp1=tmp1+ztzinv(i,j)*zty(j) 
                 end do
                 workvq1(i)=tmp1
              end do

              call rmvnorm2(maxq,jj,workvq1,ztzinv,workmhq1,zty,bc)
              
c++++++++++++ evaluating the candidate generating kernel

              call dmvn3(maxq,jj,bc,workvq1,ztzinv,logcgko,
     &                   zty,workmq1,workmq2,iflagq)

c++++++++++++ likelihood evaluation

              do i=1,jj
                 do j=1,jj
                    ztz(i,j)=0.d0
                 end do
                 zty(i)=0.d0  
              end do
              
              do i=starts(ii),ends(ii)
                 do j=iapm(i),iapm(i+1)-1
                    ztz(i-starts(ii)+1,japm(j)-starts(ii)+1)=apm(j)*
     &                   (1.d0/sigmab(ii))                       
                 end do
              end do

              loglikn=0.d0

              do i=1,nrec
                 eta=0.d0
                 offset=0.d0
                 mean=0.d0
                 gprime=0.d0
                 yij=y(i)
           
                 do j=1,p
                    eta=eta+x(i,j)*beta(j)
                    offset=offset+x(i,j)*beta(j)
                 end do

                 do j=1,q
                    if(j.lt.starts(ii).or.j.gt.ends(ii))then
                       eta=eta+z(i,j)*b(j)
                       offset=offset+z(i,j)*b(j)
                    end if   
                 end do

                 ll=0
                 do j=starts(ii),ends(ii)
                    ll=ll+1
                    eta=eta+z(i,j)*bc(ll)
                 end do

                 offset=offset+roffset(i)

                 mean=exp(eta)
                 gprime=exp(-eta)
                 ytilde=eta+(dble(yij)-mean)*gprime-offset
                 
                 loglikn=loglikn+dpoiss(dble(yij),mean,1)
 
                 do j=1,jj
                    do l=1,jj
                       ztz(j,l)=ztz(j,l)+z(i,starts(ii)+j-1)*
     &                          z(i,starts(ii)+l-1)*mean
                    end do
                    zty(j)=zty(j)+z(i,starts(ii)+j-1)*ytilde*mean
                 end do
              end do

              call invdet2(ztz,maxq,jj,ztzinv,detlog,iflagq,workvq1)

              do i=1,jj
                 tmp1=0.d0
                 do j=1,jj
                    tmp1=tmp1+ztzinv(i,j)*zty(j) 
                 end do
                 workvq1(i)=tmp1
              end do

c++++++++++++ evaluating the candidate generating kernel

              call dmvn3(maxq,jj,theta,workvq1,ztzinv,logcgkn,
     &                   zty,workmq1,workmq2,iflagq)

c++++++++++++ prior ratio
              logprioro=0.d0
              logpriorn=0.d0
            
              do i=starts(ii),ends(ii)
                 do j=iapm(i),iapm(i+1)-1
                    logpriorn=logpriorn+
     &              apm(j)*bc(i-starts(ii)+1)*
     &                     bc(japm(j)-starts(ii)+1)
                 end do
              end do

              do i=starts(ii),ends(ii)
                 do j=iapm(i),iapm(i+1)-1
                    logprioro=logprioro+
     &              apm(j)*theta(i-starts(ii)+1)*
     &                     theta(japm(j)-starts(ii)+1)
                 end do
              end do

              logpriorn=-0.5d0*logpriorn/sigmab(ii)
              logprioro=-0.5d0*logprioro/sigmab(ii)

c++++++++++++ mh step

              ratio=loglikn-logliko+logcgkn-logcgko+
     &              logpriorn-logprioro

              if(log(dble(runif())).lt.ratio)then
                 acrate(1+ii)=acrate(1+ii)+1.d0
                 do i=1,jj
                    b(starts(ii)+i-1)=bc(i) 
                 end do
              end if
           end do 

         end if          

c+++++++++++++++++++++++++++++++++
c+++++++ penalty parameters    +++ 
c+++++++++++++++++++++++++++++++++
 
         if(nsmooth.gt.0)then
            do ii=1,nsmooth
               tmp1=0.d0

               do i=starts(ii),ends(ii)
                  do j=iapm(i),iapm(i+1)-1
                     tmp1=tmp1+
     &               apm(j)*b(i)*
     &                      b(japm(j))
                  end do
               end do
               
               j=ends(ii)-starts(ii)+1-pordv(ii)
               
               sigmab(ii)=1.d0/   
     &           rgamma(0.5d0*(dble(j)+taub1),0.5d0*(tmp1+taub2))
            end do
            
         end if

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         

         if(iscan.gt.nburn)then
            skipcount=skipcount+1

c++++++++++ defining the grid

            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ regression coefficients

               tmp1=0.d0
               do i=1,nsmooth
                  tmp2=0.d0
                  do j=starts(i),ends(i)
                     tmp2=tmp2+b(j)
                  end do
                  tmp2=tmp2/dble(ends(i)-starts(i)+1)  
                  tmp1=tmp1+tmp2
               end do

               if(nfixed.gt.0)then
                  do i=1,p
                     if(i.eq.1)then
                        thetasave(isave,i)=beta(i)+tmp1
                       else
                        thetasave(isave,i)=beta(i)
                     end if   
                     betasave(i)=betasave(i)+beta(i)
                  end do
               end if   

c+++++++++++++ smoothers

               if(nsmooth.gt.0)then
                  if(nfixed.gt.0)then
                     do i=1,nsmooth
                        thetasave(isave,p+i)=sigmab(i)
                     end do
                    else
                     do i=1,nsmooth
                        thetasave(isave,i)=sigmab(i)
                     end do
                  end if

                  do i=1,nsmooth
                     tmp1=0.d0
                     do j=starts(i),ends(i)
                        tmp1=tmp1+b(j)
                     end do
                     tmp1=tmp1/dble(ends(i)-starts(i)+1)
                     do j=starts(i),ends(i)
                        randsave(isave,j)=b(j)-tmp1
                        bsave(j)=bsave(j)+b(j)
                     end do
                  end do
               end if   

               if(nsmooth.gt.0)then
                  ll=0
                  do i=1,nsmooth
                     do j=1,ngrid
                        ll=ll+1
                        tmp1=0.d0
                        tmp2=0.d0
                        do l=starts(i),ends(i)
                           tmp1=tmp1+znew(j,l)*b(l)
                           tmp2=tmp2+b(l)
                        end do
                        tmp2=tmp2/dble(ends(i)-starts(i)+1)  

                        pssave(ll,1)=pssave(ll,1)+tmp1-tmp2 
                        pssave(ll,2)=tmp1-tmp2
                     end do
                  end do
                  write(unit=1)(pssave(j,2),j=1,nsmooth*ngrid)
               end if

c+++++++++++++ cpo
               dbarc=0.d0
               do i=1,nrec
                  yij=y(i)
                  eta=0.d0
                  do j=1,p
                     eta=eta+x(i,j)*beta(j)
                  end do
                  do j=1,q
                     eta=eta+z(i,j)*b(j)
                  end do
                  eta=eta+roffset(i)

                  mean=exp(eta)
                  
                  tmp1=dpoiss(dble(yij),mean,0)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp1
                  cpo(i,2)=cpo(i,2)+tmp1                  
                  tmp1=dpoiss(dble(yij),mean,1)
                  dbarc=dbarc+tmp1
               end do

c+++++++++++++ dic
               dbar=dbar-2.d0*dbarc

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
      
      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do

      do i=1,1+nsmooth
         acrate(i)=acrate(i)/dble(nscan)
      end do

      do i=1,p
         betasave(i)=betasave(i)/dble(nsave)
      end do

      do i=1,q
         bsave(i)=bsave(i)/dble(nsave)
      end do   

      do i=1,nsmooth*ngrid
         pssave(i,1)=pssave(i,1)/dble(nsave)
         pssave(i,2)=0.d0
      end do

      rewind(unit=1)
      do ii=1,nsave
         read(1)(workvps(j),j=1,nsmooth*ngrid)
         
         do i=1,nsmooth*ngrid
            pssave(i,2)=pssave(i,2)+(workvps(i)-pssave(i,1))*
     &                              (workvps(i)-pssave(i,1))
         end do
      end do
      close(unit=1,status='delete')

      do i=1,nsmooth*ngrid
         pssave(i,2)=pssave(i,2)/dble(nsave)
      end do
 
      dhat=0.d0
      lpml=0.d0
      do i=1,nrec
         yij=y(i)
         eta=0.d0
         do j=1,p
            eta=eta+x(i,j)*betasave(j)
         end do
         do j=1,q
            eta=eta+z(i,j)*bsave(j)
         end do
         eta=eta+roffset(i)
         mean=exp(eta)
         dhat=dhat+dpoiss(dble(yij),mean,1)
         lpml=lpml+log(cpo(i,1))
      end do
      dhat=-2.d0*dhat      

      dbar=dbar/dble(nsave)
      pd=dbar-dhat
      
      mc(1)=dbar
      mc(2)=dhat
      mc(3)=pd
      mc(4)=dbar+pd
      mc(5)=lpml

      return
      end
      
     
