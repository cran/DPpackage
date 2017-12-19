c=======================================================================                      
      subroutine ptmdensity(ngrid,nrec,nvar,y,
     &                      ab,murand,sigmarand,jfr,m0,s0,nu0,tinv,
     &                      mcmcvec,nsave,tune1,tune2,tune3,
     &                      acrate,cpo,f,randsave,thetasave,
     &                      cpar,mu,sigma,
     &                      grid1,grid2,
     &                      iflag,whicho,whichn,limw,
     &                      linf,lsup,narea,mass,massi,
     &                      parti,pattern,patterns,
     &                      s,sigmainv,sigmainvc,
     &                      ybar,z,zc,zwork,vv,
     &                      workmh,workh1,workh2,workm1,workm2,
     &                      muc,sigmac,propv,propv1,propv2,
     &                      seed)
c=======================================================================                  
c
c     Subroutine `ptmdensity' to run a Markov chain for Density
c     estimation using a Multivariate Mixture of Polya Trees prior. 
c     The Multivariate Finite Polya Tree is centered in a multivariate 
c     N(mu,sigma) distribution. The Jeffery's prior is used for mu and
c     sigma.
c
c     Copyright: Alejandro Jara, 2006-2010.
c
c     Version 2.0: 
c
c     Last modification: 25-01-2010.
c
c     Changes and Bug fixes: 
c
c     Version 1.0 to Version 2.0:
c          - Centering parameters can be fixed.
c          - Proper prior can be used on the centering parameters.
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
c        ngrid       :  integer giving the size of the grid where
c                       the density estimate is evaluated.
c        nrec        :  integer giving the number of observations.
c        nvar        :  integer giving the number of variables.
c        y           :  real matrix giving the response variables,
c                       y(nrec,nvar).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        ca, cb      :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       c ~ Gamma(ca,cb). If ca<0 the precision 
c                       parameter is considered as a constant.
c        jfr         :  integer vector indicating whether Jeffery's
c                       prior is used for the centering parameters.
c        m0          :  real vector giving the mean of the normal prior
c                       for the centering mean.
c        s0          :  real matrix giving the precision of the normal
c                       prior for the centering mean.
c        nu0         :  integer giving the prior degrees of freedom
c                       for the centering covariance matrix.
c        tinv        :  real matrix giving the scale matrix for the
c                       inv-wishart prior for the centering covariance
c                       matrix.
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
c        tune1       :  real giving the tuning parameter for MH of
c                       mean baseline.
c        tune2       :  real giving the tuning parameter for MH of
c                       covariance matrix baseline.
c        tune3       :  real giving the tuning parameter for MH of
c                       precision parameter.
c        
c-----------------------------------------------------------------------
c
c---- Output ----------------------------------------------------------- 
c
c        acrate      :  real vector giving the MH acceptance rate, 
c                       acrate(3). 
c        f           :  real matrix giving the density estimate at the
c                       grid, f(ngrid,ngrid) (only in bivariate case).
c        cpo         :  real vector giving the cpos, cpo(nrec). 
c        randsave    :  real matrix containing the mcmc samples from
c                       the posterior predictive distribution, 
c                       randsave(nsave,nvar).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, 
c                       thetasave(nsave,nvar+nvar*(nvar+1)/2+1).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        cpar        :  real giving the current value of the precision
c                       parameter of the Polya Tree.
c        mu          :  real vector giving the current value of the 
c                       baseline means, mu(nvar)
c        sigma       :  real matrix giving the he current value of the
c                       baseline covariance matrix, sigma(nvar,nvar).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        binaryrep   :  integer function to compute a number based 
c                       on its binary representation. 
c        countern    :  index.
c        countero    :  index.
c        cparc       :  real giving the value of the candidate
c                       for the precision parameter.
c        detlogl     :  real giving the log of the determinat of a
c                       matrix.
c        detloglc    :  real giving the log of the determinat of a
c                       matrix.
c        detlog1     :  real giving the log of the determinat of a
c                       matrix.
c        detlog2     :  real giving the log of the determinat of a
c                       matrix.
c        detlog3     :  real giving the log of the determinat of a
c                       matrix.
c        dispcount   :  index. 
c        dlnrm       :  density of a log-normal distribution.
c        dnrm        :  density of a normal distribution.
c        evali       :  integer indicator.
c        evali2      :  integer indicator.
c        final       :  integer indicator.
c        grid1       :  real vector giving the grid where the density
c                       estimate is evaluated for first coordinate,
c                       grid1(ngrid).
c        grid2       :  real vector giving the grid where the density
c                       estimate is evaluated for second coordinate,
c                       grid2(ngrid).
c        i           :  index.
c        i1          :  index.
c        iflag       :  integer vector used to invert the baseline 
c                       covariance matrix, iflag(nvar).
c        ihmssf      :  integer function to evaluate the position of an
c                       element in a matrix based on a half-stored 
c                       version.
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index.
c        j1          :  index.
c        je2         :  index. 
c        k           :  index. 
c        k1          :  index. 
c        k2          :  index. 
c        l           :  index. 
c        limw        :  real vector giving the limits of partitions, 
c                       limw(nvar).
c        linf        :  real vector giving the limits of partitions, 
c                       linf(nvar).
c        logcgko     :  real working variable.
c        logcgkn     :  real working variable.
c        loglikec    :  real working variable.
c        loglikeo    :  real working variable.
c        logpriorc   :  real working variable.
c        logprioro   :  real working variable.
c        lsup        :  real vector giving the limits of partitions, 
c                       lsup(nvar).
c        mass        :  real vector giving the probabilities in each 
c                       partitio area, mass(narea).
c        massi       :  integer vector giving the number of observation
c                       in each element of the partition, massi(narea).
c        muc         :  real vector giving the candidate value of the 
c                       baseline means, muc(nvar)
c        narea       :  integer giving the total number of areas per 
c                       partitio, narea=2**nvar.
c        nint        :  integer indicator.
c        nscan       :  index.
c        nu          :  real working variable.
c        ok          :  integer indicator.
c        parti       :  integer vector giving the partition,
c                       parti(nvar). 
c        pattern     :  integer vector giving the pattern of an observation,
c                       pattern(nvar). 
c        patterns    :  integer vector giving the pattern of an observation,
c                       patterns(nvar). 
c        pprn        :  index.
c        prob        :  real working variable.
c        propv       :  real working matrix, propv(nvar,nvar).
c        propv1      :  real working matrix, propv1(nvar,nvar).
c        propv2      :  real working matrix, propv2(nvar,nvar).
c        quan        :  real working variable.
c        ratio       :  real working variable.
c        rgamma      :  real gamma random number generator.
c        rtlnorm     :  real truncated log normal random number generator.
c        rtnorm      :  real truncated normal random number generator.
c        runif       :  real uniform random number generator.
c        s           :  real matrix giving the sample variance, s(nvar,nvar). 
c        sec         :  cpu time working variable.
c        sec0        :  cpu time working variable.
c        sec00       :  cpu time working variable.
c        sec1        :  cpu time working variable.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        sigmac      :  real matrix giving the candidate value of the
c                       baseline covariance matrix, sigmac(nvar,nvar).
c        sigmainv    :  real matrix giving the inverse of the current
c                       value of the baseline covariance matrix, 
c                       sigmainv(nvar,nvar).
c        sigmainvc   :  real matrix giving the inverse of the candidate
c                       value of the baseline covariance matrix, 
c                       sigmainvc(nvar,nvar).
c        skipcount   :  index. 
c        sprint      :  integer function to print on screen.
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        vv          :  real working vector, vv(nvar).
c        whicho      :  integer vector giving the observation in each
c                       partition, whicho(nrec).
c        whichn      :  integer vector giving the observation in each
c                       partition, whichn(nrec).
c        workh1      :  real working vector, workh1(nvar*(nvar+1)/2)
c        workh2      :  real working vector, workh2(nvar*(nvar+1)/2)
c        workm1      :  real working matrix, workm1(nvar,nvar)
c        workm2      :  real working matrix, workm2(nvar,nvar)
c        workmh      :  real working vector, workmh(nvar*(nvar+1)/2)
c        ybar        :  real vector giving the sample mean, ybar(nvar).
c        z           :  real matrix giving the current value of the
c                       standarized observations, z(nrec,nvar).
c        zc          :  real matrix giving the candidate value of the
c                       standarized observations, zc(nrec,nvar).
c        zwork       :  real working vector, zwork(nvar).
c
c=======================================================================

      implicit none 

c+++++Data
      integer ngrid,nrec,nvar
      double precision y(nrec,nvar)

c+++++Prior information
      double precision ab(2),ca,cb
      integer murand,sigmarand,jfr(2),nu0
      double precision m0(nvar),s0(nvar,nvar)
      double precision tinv(nvar,nvar)

c+++++MCMC parameters
      integer mcmcvec(3),nburn,nskip,nsave,ndisplay
      double precision tune1,tune2,tune3

c+++++Stored output
      double precision acrate(3)
      double precision cpo(nrec)
      double precision f(ngrid,ngrid)
      double precision randsave(nsave,nvar)
      double precision thetasave(nsave,nvar+nvar*(nvar+1)/2+1)

c+++++Current values of the parameters
      double precision cpar,mu(nvar),sigma(nvar,nvar)

c+++++Working space - CPU time
      double precision sec00,sec0,sec1,sec

c+++++Working space - Density
      double precision grid1(ngrid),grid2(ngrid)

c+++++Working space - Distributions
      double precision dnrm,dlnrm
      double precision invcdfnorm

c+++++Working space - General
      integer binaryrep
      integer countero,countern
      integer evali,evali2
      integer final
      integer i,i1,j,j1,je2,k,k1,k2,l
      integer iflag(nvar)
      integer ihmssf
      integer narea,nint,ok,pattern(nvar),patterns(nvar)
      integer massi(narea)   
      integer parti(nvar)
      integer pprn,sprint
      integer whicho(nrec),whichn(nrec)
      double precision detlogl,detloglc,detlog1,detlog2
      double precision limw(nvar),linf(nvar),lsup(nvar)
      double precision mass(narea)
      double precision prob
      double precision quan
      double precision s(nvar,nvar)
      double precision sigmainv(nvar,nvar),sigmainvc(nvar,nvar)
      double precision tmp1,tmp2
      double precision ybar(nvar)
      double precision z(nrec,nvar),zc(nrec,nvar)
      double precision zwork(nvar)
      double precision vv(nvar)
      double precision workmh(nvar*(nvar+1)/2)
      double precision workh1(nvar*(nvar+1)/2),workh2(nvar*(nvar+1)/2)
      double precision workm1(nvar,nvar)
      double precision workm2(nvar,nvar)
      
c+++++Working space - MCMC scans
      integer dispcount,isave,iscan,nscan,skipcount

c+++++Working space - MH steps
      double precision cparc
      double precision logcgkn,logcgko
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision muc(nvar),sigmac(nvar,nvar)
      double precision propv(nvar,nvar),propv1(nvar,nvar),
     1  propv2(nvar,nvar)
      double precision ratio,nu

c+++++Working space - Random number generator
      integer seed(2),seed1,seed2
      double precision rtnorm,rnorm
      real runif

c++++ initialize variables
      nburn=mcmcvec(1)
      nskip=mcmcvec(2)
      ndisplay=mcmcvec(3)
      
      ca=ab(1)
      cb=ab(2)
      
c++++ set random number generator

      seed1=seed(1)
      seed2=seed(2)

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
  
      do i=1,nvar 
         parti(i)=0
         ybar(i)=0.d0
         do j=1,nvar
            s(i,j)=0.d0
         end do
      end do
      
      do i=1,nrec
         do j=1,nvar
            ybar(j)=ybar(j)+y(i,j) 
         end do   
      end do
      
      do i=1,nvar
         ybar(i)=ybar(i)/dble(nrec)
      end do   
      
      do i=1,nrec
         do j=1,nvar
            do k=1,nvar
               s(j,k)=s(j,k)+(y(i,j)-ybar(j))*(y(i,k)-ybar(k))
            end do    
         end do    
      end do
      
      do i=1,nvar
         do j=1,nvar
            s(i,j)=s(i,j)/dble(nrec)
         end do
      end do
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ First computation of loglikelihood (to reduce CPU time)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      logliko=0.d0
      
      do i=1,nvar
         do j=1,nvar
            propv(i,j)=sigma(i,j)
         end do
      end do
      call invdet(propv,nvar,sigmainv,detlogl,iflag,vv)

      do i=1,nvar
         do j=1,nvar
            propv(i,j)=0.d0            
         end do
      end do
      
      call cholesky(nvar,sigma,workmh)
      
      do i=1,nvar
         do j=1,i
            sigmainv(i,j)=workmh(ihmssf(i,j,nvar))
         end do
      end do
      call inverse(sigmainv,nvar,iflag)
      
      do i=1,nrec

         do j=1,nvar
            z(i,j)=y(i,j)-mu(j)
         end do
         
         do j=1,nvar
            tmp1=0.d0
            do k=1,nvar
               tmp1=tmp1+sigmainv(j,k)*z(i,k)   
            end do
            zwork(j)=tmp1
         end do
         
         do j=1,nvar
            z(i,j)=zwork(j)
         end do

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first observation
         if(i.eq.1)then
            logliko=-0.5d0*detlogl
            do j=1,nvar
               logliko=logliko+dnrm(z(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following observations
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nvar
               if(z(i,j).le.quan)then
                  linf(j)=-999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)=999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nvar
                  if(z(l,j).gt.lsup(j).or.z(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do
            
            logliko=logliko+
     &              log((2.d0**nvar)*cpar+dble((2.d0**nvar)*countero))-
     &              log((2.d0**nvar)*cpar+dble(i-1))

            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.20000)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nvar
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(z(i,k).le.quan)then
                    parti(k)=k1 
                    lsup(k)=quan
                   else 
                    parti(k)=k2
                    linf(k)=quan
                  end if
               end do                 
               
               countern=0
               do l=1,countero
                  final=1
                  do k=1,nvar
                     if(z(whicho(l),k).gt.lsup(k).or.
     &                  z(whicho(l),k).lt.linf(k)    )then
                        
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               logliko=logliko+
     &              log((2.d0**nvar)*cpar*dble(je2)+
     &                  dble((2.d0**nvar)*countern))-
     &              log((2.d0**nvar)*cpar*dble(je2)+dble(countero))

               if(countern.eq.0)then
                  ok=0
                else  
                  countero=countern
                  do l=1,countern
                     whicho(l)=whichn(l)
                  end do
                  j=j+1
               end if   
            end do

1           continue

            logliko=logliko-0.5d0*detlogl
            do j=1,nvar
               logliko=logliko+dnrm(z(i,j),0.d0, 1.d0, 1)
            end do   
            
         end if
        
      end do

c      call dblepr("loglik",-1,logliko,1)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Scanning the posterior distribution
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do iscan=1,nscan
  
c+++++++ check if the user has requested an interrupt
         call rchkusr()
 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating mu using a MH step                  +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(murand.eq.1)then

            do i=1,nvar
               do j=1,nvar
                  propv(i,j)=tune1*s(i,j)/dble(nrec)
               end do
            end do
         
            call rmvnorm(nvar,mu,propv,workmh,vv,muc)

            loglikn=0.d0
      
            do i=1,nrec

               do j=1,nvar
                  zc(i,j)=y(i,j)-muc(j)
               end do
         
               do j=1,nvar
                  tmp1=0.d0
                  do k=1,nvar
                     tmp1=tmp1+sigmainv(j,k)*zc(i,k)   
                  end do
                  zwork(j)=tmp1
               end do
         
               do j=1,nvar
                  zc(i,j)=zwork(j)
               end do

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

c+++++++++++++ first observation
               if(i.eq.1)then
                  loglikn=-0.5d0*detlogl
                  do j=1,nvar
                     loglikn=loglikn+dnrm(zc(i,j),0.d0, 1.d0, 1)
                  end do   

c+++++++++++++ following observations
                else

                  nint=2
                  prob=1.d0/dble(nint)
                  quan=invcdfnorm(prob,0.d0,1.d0,1,0)

                  countero=0
            
                  do j=1,nvar
                     if(zc(i,j).le.quan)then
                        linf(j)=-999.d0
                        lsup(j)=quan
                        parti(j)=1
                      else
                        linf(j)=quan
                        lsup(j)=999.d0
                        parti(j)=2
                     end if
                  end do
           
                  do l=1,i-1
                     final=1
                     do j=1,nvar
                        if(zc(l,j).gt.lsup(j).or.zc(l,j).lt.linf(j))then
                          final=0
                        end if
                     end do
               
                     if(final.eq.1)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                  end do
            
                  loglikn=loglikn+
     &              log((2.d0**nvar)*cpar+dble((2.d0**nvar)*countero))-
     &              log((2.d0**nvar)*cpar+dble(i-1))

                  if(countero.eq.0) go to 2

                  ok=1
                  j=2
                  do while(ok.eq.1.and.j.le.20000)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)

                     do k=1,nvar
                        k1=2*(parti(k)-1)+1
                        k2=2*(parti(k)-1)+2
                        quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                        if(zc(i,k).le.quan)then
                          parti(k)=k1 
                          lsup(k)=quan
                         else 
                          parti(k)=k2
                          linf(k)=quan
                        end if
                     end do                 
               
                     countern=0
                     do l=1,countero
                        final=1
                        do k=1,nvar
                           if(zc(whicho(l),k).gt.lsup(k).or.
     &                        zc(whicho(l),k).lt.linf(k)    )then
                        
                              final=0 
                          end if   
                        end do
                  
                        if(final.eq.1)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                        end if
                     end do

                     loglikn=loglikn+
     &                log((2.d0**nvar)*cpar*dble(je2)+
     &                  dble((2.d0**nvar)*countern))-
     &                log((2.d0**nvar)*cpar*dble(je2)+dble(countero))

                     if(countern.eq.0)then
                        ok=0
                      else  
                        countero=countern
                        do l=1,countern
                           whicho(l)=whichn(l)
                        end do
                        j=j+1
                     end if   
                  end do

2                 continue

                  loglikn=loglikn-0.5d0*detlogl
                  do j=1,nvar
                     loglikn=loglikn+dnrm(zc(i,j),0.d0, 1.d0, 1)
                  end do   
            
               end if
        
            end do

c++++++++++ acceptance step

            logpriorn=0.d0
            logprioro=0.d0

            if(jfr(1).eq.0)then
               do i=1,nvar
                  do j=1,nvar
                     logpriorn=logpriorn+(muc(i)-m0(i))* 
     &                                    s0(i,j)    *
     &                                   (muc(j)-m0(j))

                     logprioro=logprioro+(mu(i)-m0(i))* 
     &                                    s0(i,j)    *
     &                                   (mu(j)-m0(j))
                  end do
               end do
         
               logpriorn=-0.5d0*logpriorn
               logprioro=-0.5d0*logprioro
            end if

            ratio=loglikn-logliko+logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               do i=1,nvar
                  mu(i)=muc(i)
               end do
               do i=1,nrec
                  do j=1,nvar
                     z(i,j)=zc(i,j)
                  end do   
               end do
               logliko=loglikn
               acrate(1)=acrate(1)+1.d0
            end if

         end if

         
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating sigma using a MH step               +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(sigmarand.eq.1)then

c++++++++++ Candidate generating kernel

            nu=dble(nrec)*tune2
         
            do i=1,nvar
               do j=1,nvar
                  sigmac(i,j)=(nu-dble(nvar+1))*sigma(i,j)
               end do
            end do

            call riwishart(nvar,int(nu),sigmac,sigmainvc,workm2,vv,
     &                     workh1,workh2,iflag)

            do i=1,nvar
               do j=1,nvar
                  propv1(i,j)=sigma(i,j)
               end do
            end do
            call invdet(propv1,nvar,sigmainv,detlog1,iflag,vv)

            do i=1,nvar
               do j=1,nvar
                  propv2(i,j)=sigmac(i,j)
               end do
            end do
            call invdet(propv2,nvar,propv,detlog2,iflag,vv)

            logcgkn=0.d0
            logcgkn=-(nu+dble(nvar+1))*detlog1
            logcgkn=logcgkn+nu*detlog2

            tmp2=0.d0
            do i=1,nvar
               do j=1,nvar
                  tmp1=0.d0 
                  do k=1,nvar
                     tmp1=tmp1+(nu-dble(nvar+1))*
     &                          sigmac(i,k)*sigmainv(k,j)        
                  end do
                  if(i.eq.j)tmp2=tmp2+tmp1
               end do
            end do   
            logcgkn=logcgkn-tmp2
            logcgkn=logcgkn/dble(2)

            logcgko=0.d0
            logcgko=-(nu+dble(nvar+1))*detlog2
            logcgko=logcgko+nu*detlog1
         
            tmp2=0.d0
            do i=1,nvar
               do j=1,nvar
                  tmp1=0.d0 
                  do k=1,nvar
                     tmp1=tmp1+(nu-dble(nvar+1))*
     &                          sigma(i,k)*sigmainvc(k,j)        
                  end do
                  if(i.eq.j)tmp2=tmp2+tmp1
               end do
            end do   
            logcgko=logcgko-tmp2
            logcgko=logcgko/dble(2)

c++++++++++ End candidate generating kernel

            loglikn=0.d0
      
            do i=1,nvar
               do j=1,nvar
                  propv(i,j)=sigmac(i,j)
               end do
            end do
            call invdet(propv,nvar,sigmainvc,detloglc,iflag,vv)

            do i=1,nvar
               do j=1,nvar
                   propv(i,j)=0.d0 
               end do
            end do
            call cholesky(nvar,sigmac,workmh)
      
            do i=1,nvar
               do j=1,i
                  sigmainvc(i,j)=workmh(ihmssf(i,j,nvar))
               end do
            end do
            call inverse(sigmainvc,nvar,iflag)

            do i=1,nrec

               do j=1,nvar
                  zc(i,j)=y(i,j)-mu(j)
               end do
         
               do j=1,nvar
                  tmp1=0.d0
                  do k=1,nvar
                     tmp1=tmp1+sigmainvc(j,k)*zc(i,k)   
                  end do
                  zwork(j)=tmp1
               end do
         
               do j=1,nvar
                  zc(i,j)=zwork(j)
               end do

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

c+++++++++++++ first observation
               if(i.eq.1)then
                  loglikn=-0.5d0*detloglc
                  do j=1,nvar
                     loglikn=loglikn+dnrm(zc(i,j),0.d0, 1.d0, 1)
                  end do   

c+++++++++++++ following observations
                else

                  nint=2
                  prob=1.d0/dble(nint)
                  quan=invcdfnorm(prob,0.d0,1.d0,1,0)

                  countero=0
            
                  do j=1,nvar
                     if(zc(i,j).le.quan)then
                        linf(j)=-999.d0
                        lsup(j)=quan
                        parti(j)=1
                      else
                        linf(j)=quan
                        lsup(j)=999.d0
                        parti(j)=2
                     end if
                  end do
           
                  do l=1,i-1
                     final=1
                     do j=1,nvar
                        if(zc(l,j).gt.lsup(j).or.zc(l,j).lt.linf(j))then
                          final=0
                        end if
                     end do
               
                     if(final.eq.1)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                  end do
            
                  loglikn=loglikn+
     &              log((2.d0**nvar)*cpar+dble((2.d0**nvar)*countero))-
     &              log((2.d0**nvar)*cpar+dble(i-1))

                  if(countero.eq.0) go to 3

                  ok=1
                  j=2
                  do while(ok.eq.1.and.j.le.20000)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)

                     do k=1,nvar
                        k1=2*(parti(k)-1)+1
                        k2=2*(parti(k)-1)+2
                        quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                        if(zc(i,k).le.quan)then
                          parti(k)=k1 
                          lsup(k)=quan
                         else 
                          parti(k)=k2
                          linf(k)=quan
                        end if
                     end do                 
               
                     countern=0
                     do l=1,countero
                        final=1
                        do k=1,nvar
                           if(zc(whicho(l),k).gt.lsup(k).or.
     &                        zc(whicho(l),k).lt.linf(k)    )then
                        
                              final=0 
                           end if   
                        end do
                  
                        if(final.eq.1)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                        end if
                     end do

                     loglikn=loglikn+
     &                 log((2.d0**nvar)*cpar*dble(je2)+
     &                  dble((2.d0**nvar)*countern))-
     &                 log((2.d0**nvar)*cpar*dble(je2)+dble(countero))

                     if(countern.eq.0)then
                        ok=0
                      else  
                        countero=countern
                        do l=1,countern
                           whicho(l)=whichn(l)
                        end do
                        j=j+1
                     end if   
                  end do

3                 continue

                  loglikn=loglikn-0.5d0*detloglc
                  do j=1,nvar
                     loglikn=loglikn+dnrm(zc(i,j),0.d0, 1.d0, 1)
                  end do   
               
               end if
        
            end do

c++++++++++ acceptance step

            logprioro= -(dble(nvar+1)/dble(2))*(detlog2 - detlog1)
            logprioro=0.d0

            if(jfr(2).eq.0)then
               call diwishart(nvar,nu0,sigmac,tinv,workm1,workm2,
     &                        vv,iflag,logpriorn)        

               call diwishart(nvar,nu0,sigma,tinv,workm1,workm2,
     &                        vv,iflag,logprioro)        

            end if

            ratio=loglikn-logliko+logcgkn-logcgko+
     &            logpriorn-logprioro 

            if(log(dble(runif())).lt.ratio)then
               do i=1,nvar
                  do j=1,nvar
                     sigma(i,j)=sigmac(i,j)
                     sigmainv(i,j)=sigmainvc(i,j)
                  end do
               end do
               do i=1,nrec
                  do j=1,nvar
                     z(i,j)=zc(i,j)
                  end do   
               end do 
               detlogl=detloglc
               logliko=loglikn
               acrate(2)=acrate(2)+1.d0
             else
               do i=1,nvar
                  do j=1,nvar
                     propv(i,j)=sigma(i,j)
                  end do
               end do
               call invdet(propv,nvar,sigmainv,detlogl,iflag,vv)
               do i=1,nvar
                  do j=1,nvar
                     propv(i,j)=sigma(i,j)
                  end do
               end do
               call cholesky(nvar,sigma,workmh)
      
               do i=1,nvar
                  do j=1,i
                     sigmainv(i,j)=workmh(ihmssf(i,j,nvar))
                  end do
               end do
               call inverse(sigmainv,nvar,iflag)
            end if

         end if

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ MH to update the c parameter                 +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(ca.gt.0.d0)then

c++++++++++ sample candidates

            cparc=exp(rnorm(log(cpar),tune3*1.d0))
            logcgkn=dlnrm(cpar ,log(cparc),tune3*1.0,1) 
            logcgko=dlnrm(cparc,log(cpar ),tune3*1.0,1) 

c++++++++++ evaluate log-prior for candidate value of the parameters

            call dgamma2(cparc,ca,cb,logpriorn)  

c++++++++++ evaluate log-prior for current value of parameters

            call dgamma2(cpar ,ca,cb,logprioro)

            loglikn=0.d0
      
            do i=1,nrec

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

c+++++++++++++ first observation
               if(i.eq.1)then
                  loglikn=-0.5d0*detlogl
                  do j=1,nvar
                     loglikn=loglikn+dnrm(z(i,j),0.d0, 1.d0, 1)
                  end do   

c+++++++++++++ following observations
                else

                  nint=2
                  prob=1.d0/dble(nint)
                  quan=invcdfnorm(prob,0.d0,1.d0,1,0)

                  countero=0
            
                  do j=1,nvar
                     if(z(i,j).le.quan)then
                        linf(j)=-999.d0
                        lsup(j)=quan
                        parti(j)=1
                      else
                        linf(j)=quan
                        lsup(j)=999.d0
                        parti(j)=2
                     end if
                  end do
           
                  do l=1,i-1
                     final=1
                     do j=1,nvar
                        if(z(l,j).gt.lsup(j).or.z(l,j).lt.linf(j))then
                          final=0
                        end if
                     end do
               
                     if(final.eq.1)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                  end do
            
                  loglikn=loglikn+
     &             log((2.d0**nvar)*cparc+dble((2.d0**nvar)*countero))-
     &             log((2.d0**nvar)*cparc+dble(i-1))

                  if(countero.eq.0) go to 4

                  ok=1
                  j=2
                  do while(ok.eq.1.and.j.le.20000)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)

                     do k=1,nvar
                        k1=2*(parti(k)-1)+1
                        k2=2*(parti(k)-1)+2
                        quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                        if(z(i,k).le.quan)then
                          parti(k)=k1 
                          lsup(k)=quan
                         else 
                          parti(k)=k2
                          linf(k)=quan
                        end if
                     end do                 
               
                     countern=0
                     do l=1,countero
                        final=1
                        do k=1,nvar
                           if(z(whicho(l),k).gt.lsup(k).or.
     &                        z(whicho(l),k).lt.linf(k)    )then
                        
                              final=0 
                           end if   
                        end do
                  
                        if(final.eq.1)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                        end if
                     end do

                     loglikn=loglikn+
     &                 log((2.d0**nvar)*cparc*dble(je2)+
     &                     dble((2.d0**nvar)*countern))-
     &                 log((2.d0**nvar)*cparc*dble(je2)+dble(countero))

                     if(countern.eq.0)then
                        ok=0
                      else  
                        countero=countern
                        do l=1,countern
                           whicho(l)=whichn(l)
                        end do
                        j=j+1
                     end if   
                  end do

4                 continue
  
                  loglikn=loglikn-0.5d0*detlogl
                  do j=1,nvar
                     loglikn=loglikn+dnrm(z(i,j),0.d0, 1.d0, 1)
                  end do   
            
               end if
        
            end do

c++++++++++ acceptance step
            ratio=loglikn+logpriorn-logliko-logprioro+
     &                 logcgkn-logcgko

            if(log(dble(runif())).lt.ratio)then
               cpar=cparc
               acrate(3)=acrate(3)+1.d0
               logliko=loglikn
            end if            

         end if 


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Save samples                                 +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1
              
c+++++++++++++ mu
               do i=1,nvar
                  thetasave(isave,i)=mu(i)
               end do   

c+++++++++++++ sigma
               k1=0
               do i=1,nvar
                  do j=i,nvar
                     k1=k1+1
                     thetasave(isave,nvar+k1)=sigma(i,j)
                  end do
               end do

c+++++++++++++ c parameter
               
               thetasave(isave,nvar+k1+1)=cpar

c+++++++++++++ cpo

               do i=1,nrec

                  loglikn=0.d0

c++++++++++++++++ check if the user has requested an interrupt
                  call rchkusr()

                  nint=2
                  prob=1.d0/dble(nint)
                  quan=invcdfnorm(prob,0.d0,1.d0,1,0)

                  countero=0
             
                  do j=1,nvar
                     if(z(i,j).le.quan)then
                        linf(j)=-999.d0
                        lsup(j)=quan
                        parti(j)=1
                      else
                        linf(j)=quan
                        lsup(j)=999.d0
                        parti(j)=2
                     end if
                  end do
           
                  do l=1,nrec
                     final=1
                     if(l.ne.i)then
                     do j=1,nvar
                        if(z(l,j).gt.lsup(j).or.z(l,j).lt.linf(j))then
                          final=0
                        end if
                     end do
               
                     if(final.eq.1)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                     end if
                  end do
            
                  loglikn=loglikn+
     &              log((2.d0**nvar)*cpar+dble((2.d0**nvar)*countero))-
     &              log((2.d0**nvar)*cpar+dble(nrec-1))

                  if(countero.eq.0) go to 5

                  ok=1
                  j=2
                  do while(ok.eq.1.and.j.le.20000)
                     nint=2**j
                     je2=j**2
                     prob=1.d0/dble(nint)

                     do k=1,nvar
                        k1=2*(parti(k)-1)+1
                        k2=2*(parti(k)-1)+2
                        quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                        if(z(i,k).le.quan)then
                          parti(k)=k1 
                          lsup(k)=quan
                         else 
                          parti(k)=k2
                          linf(k)=quan
                        end if
                     end do                 
               
                     countern=0
                     do l=1,countero
                        final=1
                        do k=1,nvar
                           if(z(whicho(l),k).gt.lsup(k).or.
     &                        z(whicho(l),k).lt.linf(k)    )then
                        
                              final=0 
                           end if   
                        end do
                  
                        if(final.eq.1)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                        end if
                     end do

                     loglikn=loglikn+
     &                log((2.d0**nvar)*cpar*dble(je2)+
     &                    dble((2.d0**nvar)*countern))-
     &                log((2.d0**nvar)*cpar*dble(je2)+dble(countero))

                     if(countern.eq.0)then
                        ok=0
                      else  
                        countero=countern
                        do l=1,countern
                           whicho(l)=whichn(l)
                        end do
                        j=j+1
                     end if   
                  end do

5                 continue

                  loglikn=loglikn-0.5d0*detlogl
                  do j=1,nvar
                     loglikn=loglikn+dnrm(z(i,j),0.d0, 1.d0, 1)
                  end do   
            
                  cpo(i)=cpo(i)+1.0d0/exp(loglikn)
               end do

c+++++++++++++ sample from the predictive distribution

               nint=2
               prob=1.d0/dble(nint)
               quan=invcdfnorm(prob,0.d0,1.d0,1,0)

               countero=0
               
               do i=1,narea
                  massi(i)=0
                  mass(i)=0.d0
               end do
               
               do l=1,nrec
                  do j=1,nvar
                     evali=1 
                     if(z(l,j).le.quan)evali=0  
                     pattern(j)=evali
                     zwork(j)=z(l,j)
                  end do
                  evali=binaryrep(nvar,pattern)
                  massi(evali)=massi(evali)+1
               end do   
      
               do i=1,narea
                  mass(i)=(cpar+dble(massi(i)))/
     &                    ((2**nvar)*cpar+dble(nrec))
               end do
               
               call simdisc(mass,narea,narea,evali)  
               evali2=evali
               call binaryrepinv(nvar,evali2,patterns)
                
               do l=1,nrec
                  final=1
                  do j=1,nvar
                     evali=1 
                     if(z(l,j).le.quan)evali=0  
                     pattern(j)=evali
                     if(pattern(j).ne.patterns(j))final=0
                  end do
               
                  if(final.eq.1)then
                     countero=countero+1
                     whicho(countero)=l
                  end if   
               end do 
               
               do i=1,nvar
                  if(patterns(i).eq.0)then
                    linf(i)=-999.d0
                    lsup(i)=quan
                    parti(i)=1
                   else
                    linf(i)=quan
                    lsup(i)=999.d0
                    parti(i)=2
                  end if 
               end do

               if(countero.eq.0) go to 6

               ok=1
               j=2
               countern=0
               
               do while(ok.eq.1.and.j.le.20000)
                  nint=2**j
                  je2=j**2
                  prob=1.d0/dble(nint)
                  
                  do k=1,nvar
                     
                     k1=2*(parti(k)-1)+1
                     k2=2*(parti(k)-1)+2
                     quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)

                     limw(k)=quan
                     
                     if(quan.gt.lsup(k).or.quan.lt.linf(k))then
                        call intpr("j",-1,j,1)
                        call dblepr("linf",-1,linf,nvar)
                        call dblepr("lsup",-1,lsup,nvar)
                        call dblepr("limw",-1,limw,nvar)
                        call rexit("Errors in limits")
                     end if
                  end do   
                     
                  do k=1,narea
                     massi(k)=0
                     mass(k)=0.d0
                  end do
               
                  do l=1,countero
                     do k=1,nvar
                        evali=1 
                        if(z(whicho(l),k).le.limw(k))evali=0  
                        pattern(k)=evali
                     end do
                     evali=binaryrep(nvar,pattern)
                     massi(evali)=massi(evali)+1
                  end do                      

                  do l=1,narea
                     mass(l)=(cpar*dble(je2)+dble(massi(l)))/
     &                   ((2**nvar)*cpar*dble(je2)+dble(countero))
                  end do

                  call simdisc(mass,narea,narea,evali)  
                  evali2=evali
                  call binaryrepinv(nvar,evali2,patterns)

                  countern=0
                  do l=1,countero
                     final=1
                     do k=1,nvar
                        evali=1 
                        if(z(whicho(l),k).le.limw(k))evali=0  
                        pattern(k)=evali
                        if(pattern(k).ne.patterns(k))final=0
                     end do
                     if(final.eq.1)then
                        countern=countern+1
                        whichn(countern)=whicho(l)
                     end if   
                  end do                      

                  do k=1,nvar
                     if(patterns(k).eq.0)then
                       parti(k)=2*(parti(k)-1)+1
                       lsup(k)=limw(k)
                      else
                       parti(k)=2*(parti(k)-1)+2
                       linf(k)=limw(k)
                     end if 
                  end do

                  if(countern.eq.0)then
                     ok=0
                   else  
                     countero=countern
                     do l=1,countern
                        whicho(l)=whichn(l)
                     end do
                     j=j+1
                  end if   
               end do

6              continue 
               
               do k=1,nvar
                  if(linf(k).ge.lsup(k))then
                     call intpr("j",-1,j,1)
                     call dblepr("linf",-1,linf,nvar)
                     call dblepr("lsup",-1,lsup,nvar)
                     call rexit("Errors in limits")
                  end if
               end do   

               do i=1,nvar
                  zwork(i)=rtnorm(0.d0,1.d0,linf(i),
     &                            lsup(i),.false.,.false.)
               end do

               do i=1,nvar
                  do j=1,nvar
                     propv(i,j)=sigma(i,j)
                  end do
               end do
        
               call cholesky(nvar,propv,workmh)

               do i=1,nvar
                  do j=1,nvar
                     propv(i,j)=0.d0
                  end do
               end do
       
               do i=1,nvar
                  do j=1,i
                     propv(i,j)=workmh(ihmssf(i,j,nvar))
                  end do
               end do

               do j=1,nvar
                  tmp1=0.d0
                  do k=1,nvar
                     tmp1=tmp1+propv(j,k)*zwork(k)   
                  end do
                  randsave(isave,j)=mu(j)+tmp1
               end do

c+++++++++++++ density 

               if(nvar.eq.2)then
                  do i1=1,ngrid
                     do j1=1,ngrid
            
c++++++++++++++++++++++ check if the user has requested an interrupt
                        call rchkusr()

                        loglikn=0.d0
      
                        vv(1)=grid1(i1)
                        vv(2)=grid2(j1)
                        
                        do j=1,nvar
                           vv(j)=vv(j)-mu(j)
                        end do
         
                        do j=1,nvar
                           tmp1=0.d0
                           do k=1,nvar
                              tmp1=tmp1+sigmainv(j,k)*vv(k)   
                           end do
                           zwork(j)=tmp1
                        end do
         
                        do j=1,nvar
                           vv(j)=zwork(j)
                        end do

                        nint=2
                        prob=1.d0/dble(nint)
                        quan=invcdfnorm(prob,0.d0,1.d0,1,0)

                        countero=0
            
                        do j=1,nvar
                           if(vv(j).le.quan)then
                              linf(j)=-999.d0
                              lsup(j)=quan
                              parti(j)=1
                            else
                              linf(j)=quan
                              lsup(j)=999.d0
                              parti(j)=2
                           end if
                        end do
           
                        do l=1,nrec
                           final=1
                           do j=1,nvar
                              if(z(l,j).gt.lsup(j).or.
     &                           z(l,j).lt.linf(j))then
                                 final=0
                              end if
                           end do
               
                           if(final.eq.1)then
                              countero=countero+1
                              whicho(countero)=l
                           end if   
                        end do
            
                        loglikn=loglikn+
     &                      log((2.d0**nvar)*cpar+
     &                      dble((2.d0**nvar)*countero))-
     &                      log((2.d0**nvar)*cpar+dble(nrec))

                        if(countero.eq.0) go to 7

                        ok=1
                        j=2
                        do while(ok.eq.1.and.j.le.20000)
                           nint=2**j
                           je2=j**2
                           prob=1.d0/dble(nint)

                           do k=1,nvar
                              k1=2*(parti(k)-1)+1
                              k2=2*(parti(k)-1)+2
                              quan=invcdfnorm(dble(k1)*prob,
     &                             0.d0,1.d0,1,0)
               
                              if(vv(k).le.quan)then
                                parti(k)=k1 
                                lsup(k)=quan
                               else 
                                parti(k)=k2
                                linf(k)=quan
                              end if
                           end do                 
               
                           countern=0
                           do l=1,countero
                              final=1
                              do k=1,nvar
                                 if(z(whicho(l),k).gt.lsup(k).or.
     &                              z(whicho(l),k).lt.linf(k)    )then
                        
                                    final=0 
                                 end if   
                              end do
                  
                              if(final.eq.1)then
                                countern=countern+1
                                whichn(countern)=whicho(l)
                              end if
                           end do

                           loglikn=loglikn+
     &                       log((2.d0**nvar)*cpar*dble(je2)+
     &                           dble((2.d0**nvar)*countern))-
     &                       log((2.d0**nvar)*cpar*dble(je2)+
     &                            dble(countero))

                           if(countern.eq.0)then
                              ok=0
                            else  
                              countero=countern
                              do l=1,countern
                                 whicho(l)=whichn(l)
                              end do
                              j=j+1
                           end if   
                        end do

7                       continue

                        loglikn=loglikn-0.5d0*detlogl
                        do j=1,nvar
                           loglikn=loglikn+dnrm(vv(j),0.d0, 1.d0, 1)
                        end do   
                        
                        f(i1,j1)=f(i1,j1)+exp(loglikn)

                     end do
                  end do   
               end if

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

      do i=1,3
         acrate(i)=acrate(i)/dble(nscan)      
      end do   
     
      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do

      do i=1,ngrid
         do j=1,ngrid
            f(i,j)=f(i,j)/dble(nsave)
         end do
      end do

      return
      end
      
