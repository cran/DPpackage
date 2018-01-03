    
c=======================================================================                      
      subroutine psgamprob
     &                (nrec,nfixed,p,nsmooth,q,starts,ends,nsmooths,     
     &                 maxq,x,z,yr,roffset,ngrid,znew,xreal,      
     &                 possfp,                                          
     &                 prec,sb,taub,iapm,japm,apm,nznh,pordv,      
     &                 beta,b,sigmab,y,                                 
     &                 mcmc,nsave,                                      
     &                 cpo,randsave,thetasave,pssave,                   
     &                 seed,                                            
     &                 iflagp,workmhp1,workvp1,          
     &                 xtx,xty,                                          
     &                 iflagq,bc,workmhq1,workvq1,             
     &                 ztz,zty,ztzinv,theta,                             
     &                 mc,betasave,bsave,workvps)                       
c=======================================================================                      
c     # of arguments = 53.
c
c     Subroutine `psgamprob' to run a Markov chain in a semiparametric 
c     probit model, using a B-splines and penalties.
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
c        y           :  real vector giving the response variable,
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
      integer starts(nsmooths),ends(nsmooths)
      integer possfp(nsmooths)
      integer yr(nrec)
      double precision roffset(nrec),x(nrec,p),z(nrec,q),znew(ngrid,q)
      double precision xreal(ngrid,nsmooths)
      
c+++++Prior 
      integer nznh,iapm(q+1),japm(nznh)
      double precision prec(p,p),sb(p)
      double precision taub(2),taub1,taub2
      double precision apm(nznh)
      double precision pordv(nsmooths)
      
c+++++Current values of the parameters
      double precision y(nrec) 
      double precision beta(p)
      double precision b(q)
      double precision sigmab(nsmooths)

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      double precision cpo(nrec,2)
      double precision randsave(nsave,q)
      double precision thetasave(nsave,p+nsmooth)
      double precision pssave(ngrid*nsmooths,2)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++Fixed effects
      integer iflagp(p) 
      double precision workmhp1(p*(p+1)/2)
      double precision workvp1(p)
      double precision xtx(p,p),xty(p)

c+++++Smoothers
      integer iflagq(maxq) 
      double precision bc(maxq)
      double precision theta(maxq)
      double precision workmhq1(maxq*(maxq+1)/2)
      double precision workvq1(maxq)
      double precision ztz(maxq,maxq),zty(maxq)
      double precision ztzinv(maxq,maxq)
      double precision workvps(ngrid*nsmooths)

c++++ model's performance
      double precision mc(5)
      double precision betasave(p),bsave(q)
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer i,ii,j,jj,l,ll
      integer sprint
      double precision detlog
      double precision disp      
      double precision eta,gprime
      double precision mean
      double precision offset
      double precision tmp1,tmp2
      double precision ytilde
      double precision yij
      logical ainf,asup
      
c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c+++++RNG and distributions
      double precision rgamma,cdfnorm,dbin,rtnorm

c++++ model's performance
      double precision dbarc,dbar,dhat,pd,lpml

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

      disp=1.d0
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

c++++++++++++++++++++++++++++++++++
c+++++++ latent variable
c++++++++++++++++++++++++++++++++++

         do i=1,nrec
            tmp1=0.d0
            if(nfixed.gt.0)then
              do j=1,p
                 tmp1=tmp1+x(i,j)*beta(j)
              end do
            end if
            
            do j=1,q
               tmp1=tmp1+z(i,j)*b(j) 
            end do
            
            if(yr(i).eq.1)then
              ainf=.false.
              asup=.true.
              y(i)=rtnorm(tmp1,1.d0,0.d0,0.d0,ainf,asup) 
            end if
            
            if(yr(i).eq.0)then
              ainf=.true.
              asup=.false.
              y(i)=rtnorm(tmp1,1.d0,0.d0,0.d0,ainf,asup) 
            end if
         end do

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

               mean=eta
               gprime=disp
               ytilde=(yij-offset)*gprime
 
               do j=1,p
                  do l=1,p
                     xtx(j,l)=xtx(j,l)+x(i,j)*x(i,l)*gprime
                  end do
                  xty(j)=xty(j)+x(i,j)*ytilde
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

            call rmvnorm(p,workvp1,xtx,workmhp1,xty,beta)

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

                 mean=eta
                 gprime=disp
                 ytilde=(yij-offset)*gprime

                 do j=1,jj
                    do l=1,jj
                       ztz(j,l)=ztz(j,l)+z(i,starts(ii)+j-1)*
     &                          z(i,starts(ii)+l-1)*disp
                    end do
                    zty(j)=zty(j)+z(i,starts(ii)+j-1)*ytilde
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
              
              do i=1,jj
                 b(starts(ii)+i-1)=bc(i) 
              end do

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

                  mean=cdfnorm(eta,0.d0,1.d0,1,0)

                  tmp2=dbin(dble(yr(i)),1.d0,mean,0)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp2
                  cpo(i,2)=cpo(i,2)+tmp2

                  tmp2=dbin(dble(yr(i)),1.d0,mean,1)
                  dbarc=dbarc+tmp2
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
         eta=0.d0
         do j=1,p
            eta=eta+x(i,j)*betasave(j)
         end do
         do j=1,q
            eta=eta+z(i,j)*bsave(j)
         end do
         eta=eta+roffset(i)

         mean=cdfnorm(eta,0.d0,1.d0,1,0)

         tmp1=dbin(dble(yr(i)),1.d0,mean,1)
         dhat=dhat+tmp1   
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
      
     
