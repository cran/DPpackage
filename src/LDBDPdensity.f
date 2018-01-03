c=======================================================================                      
      subroutine ldbdpdensity(y,x,nrec,p,
     &                        npred,ngrid,grid,xpred,
     &                        maxn,nu,alpha,lambda,tau1,tau2,psiinv,
     &                        s0invm,s0inv,
     &                        kk,gp,beta,mub,sb,v,
     &                        mcmc,nsave,slice,
     &                        acrate,thetasave,randsave,
     &                        fmean,flow,fupp,meanfpm,meanfpl,meanfph, 
     &                        cpo,seed,
     &                        iflag,sbinv,workm1,workv1,workv2,
     &                        workmh1,workmh2,betal,betar,beta1,
     &                        vl,vr,v1,workdpw,weight,fw,fw2,
     &                        fs,fm,worksam)

c=======================================================================
c     # 57 arguments
c
c     Subroutine `ldbdpdensity' to run a Markov chain for a 
c     Linear Dependent Bernstein-Dirichlet Process prior for 
c     bounded conditional density estimation.
c
c     Copyright: Felipe Barrientos and Alejandro Jara, 2010.
c
c     Version 1.0: 
c
c     Last modification: 30-08-2010.
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
c      Felipe Barrientos
c      Department of Statistics
c      Facultad de Matematicas
c      Pontificia Universidad Catolica de Chile
c      Casilla 306, Correo 22 
c      Santiago
c      Chile
c      Email: afbarrie@uc.cl
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
c        nrec        :  integer giving the number of data points. 
c        p           :  integer giving the number of predictors.
c        y           :  real vector giving the responses, y(nrec). 
c        x           :  real matrix giving the design matrix, x(nrec,p).
c
c-----------------------------------------------------------------------
c
c---- Prediction -------------------------------------------------------
c 
c        cband       :  integer value indicating whether the 
c                       credible bands need to be computed or not.
c        ngrid       :  integer giving the number of grid points where 
c                       the density estimates are evaluated.. 
c        npred       :  integer giving the number of predictions.
c        grid        :  real vector giving the grid, grid(ngrid). 
c        tband       :  integer indicating the type of credible 
c                       band that need to be computed.
c        xpred       :  real matrix giving the design matrix of the 
c                       predictions, zpred(npred,p).
c
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        alpha       :  real giving the value of the DP total mass 
c                       parameter.
c        lambda      :  real giving the parameter of the truncated
c                       Poisson prior for the degree of the Bernstein
c                       polynomial.
c        maxn        :  integer giving the truncation of the DP.
c        tau1        :  real giving the value of the tau1 hyperparameter
c                       in the Gamma component of the centering 
c                       distribution. 
c        tau1,tau2   :  real giving the values of hyperparameters of the
c                       Gamma prior for the g-prior (only internal use).
c        s0invm      :  real vector giving the product of the prior
c                       precision and the prior mean of the normal
c                       normal prior for mub, s0ivnm(p).
c        s0inv       :  real vector giving the inverse of the covariance
c                       matrix of the normal prior for mub, 
c                       s0inv(p,p).
c        nu          :  integer giving the degres of freedom parameter
c                       of the inverted-Wishart prior for sb.
c        psiinv      :  real matrix giving the scale matrix of the
c                       inverted-Wishart prior for sd, psiinv(p,p).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        kk          :  integer giving the degree of the Bernstein 
c                       polynomial.
c        gp          :  real giving the value of the G-prior 
c                       parameter (only for internal use).
c        beta        :  real matrix giving the value of the regression  
c                       coefficients, beta(maxn,p).
c        mub         :  real vector giving the mean of the normal 
c                       centering distribution, mub(p).
c        sb          :  real matrix giving the variance of the normal
c                       centering distribution, sb(p,p).
c        v           :  real vector giving the stick-breaking
c                       parameters, v(maxm).
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
c        betaw       :  real giving the Slice sampling parameter for
c                       the regression coefficients.
c        gw          :  real giving the Slice sampling parameter for
c                       the g-prior parameter.
c        vw          :  real giving the Slice sampling parameter for
c                       the stick-breaking parameters.
c        
c-----------------------------------------------------------------------
c
c---- Output -----------------------------------------------------------
c
c        acrate      :  real vector giving the average number of
c                       evaluation of the posterior in the Slice
c                       sampling step and the acceptance rate 
c                       of the MH step for the degree of the 
c                       polynomial, acrate(2).
c        cpo         :  real giving the cpos and fsos, cpo(nrec,2). 
c        randsave    :  real matrix containing the mcmc samples for
c                       the regression coeff and stick-breaking
c                       parameters, randsave(nsave,p*maxm,maxm).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the parameters, thetasave(nsave,2+p+p*(p+1)/2)
c        fmean       :  real matrix giving the posterior mean of the 
c                       density, fmean(npred,ngrid).
c        flow        :  real matrix giving the lower limit of the 
c                       HPD of the density, flow(npred,ngrid).
c        fupp        :  real matrix giving the upper limit of the  
c                       HPD of the density, fupp(npred,ngrid).
c        meanfpm     :  real vector giving the posterior mean of the 
c                       mean function, meanfpm(npred).
c        meanfpl     :  real vector giving the lower limit of the 
c                       HPD of the mean function meanfpl(npred).
c        meanfph     :  real vector giving the upper limit of the  
c                       HPD of the mean function, meanfph(npred).
c
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        betal       :  real matrix used for the slice-sampling of the
c                       regression coefficients, betal(maxn,p).
c        betar       :  real matrix used for the slice-sampling of the
c                       regression coefficients, betar(maxn,p).
c        beta1       :  real matrix used for the slice-sampling of the
c                       regression coefficients, beta1(maxn,p).
c        fs          :  real vector used to evaluate the conditional
c                       densities, fs(ngrid).
c        fm          :  real vector used to evaluate the conditional
c                       means, fs(npred).
c        fw          :  real vector used to compute cpos, 
c                       fw(nrec).
c        fw2         :  real matrix used to evaluate the predictions,
c                       fw2(npred,ngrid).
c        iflag       :  integer vector used to evaluate reg. coeff.,
c                       iflag(p).
c        sbinv       :  real matrix used to keep the inverse of sb,
c                       sbinv(p,p).
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        vl          :  real vector used for the slice sampling
c                       of stick-breaking parameters, vl(maxn).
c        vr          :  real vector used for the slice sampling
c                       of stick-breaking parameters, vr(maxn).
c        v1          :  real vector used for the slice sampling
c                       of stick-breaking parameters, v1(maxn).
c        weight      :  real vector used to compute the DP weights,
c                       weight(maxn).
c        workdpw     :  real vector used to compute the DP weights,
c                       workdpw(maxn+1).                       

c        workm1      :  real matrix used to update centering param.,
c                       workm1(p,p).
c        workmh1     :  real vector used to update centering param.,
c                       workmh1(p*(p+1)/2).
c        workmh2     :  real vector used to update centering param.,
c                       workmh2(p*(p+1)/2).
c        workv1      :  real vector used to update centering param.,
c                       workv1(p).
c        workv2      :  real vector used to update centering param.,
c                       workv2(p).
c        worksam     :  real vector used to comput HPD bands,
c                       worksam(nsave).
c
c=======================================================================                  

      implicit none 

c+++++Data
      integer nrec,p
      double precision x(nrec,p)
      double precision y(nrec)

c+++++Predictions
      integer cband,npred,ngrid,tband
      double precision  grid(ngrid),xpred(npred,p)

c+++++Prior 
      integer maxn,nu
      double precision alpha
      double precision lambda
      double precision tau1,tau2
      double precision psiinv(p,p)
      double precision s0invm(p)
      double precision s0inv(p,p)
      
c+++++Current values of the parameters
      integer kk
      double precision gp
      double precision beta(maxn,p)
      double precision mub(p)
      double precision sb(p,p)
      double precision v(maxn)

c+++++MCMC parameters
      integer mcmc(5),nburn,nskip,nsave,ndisplay
      double precision slice(3),betaw,gw,vw

c+++++Output
      double precision acrate(2)
      double precision thetasave(nsave,1+p+p*(p+1)/2+1)
      double precision randsave(nsave,p*maxn+maxn)
      double precision fmean(npred,ngrid)
      double precision flow(npred,ngrid)
      double precision fupp(npred,ngrid)
      double precision meanfpm(npred)
      double precision meanfpl(npred)
      double precision meanfph(npred)
      double precision cpo(nrec,2)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c++++ Inversion
      integer iflag(p)
      double precision sbinv(p,p)
      double precision workm1(p,p)
      double precision workv1(p)
      double precision workv2(p)
      double precision workmh1(p*(p+1)/2)
      double precision workmh2(p*(p+1)/2)

c++++ Slice sampler
      double precision betal(maxn,p)
      double precision betar(maxn,p)
      double precision beta1(maxn,p)
      double precision vl(maxn)
      double precision vr(maxn)
      double precision v1(maxn) 

c++++ DP
      double precision workdpw(maxn+1)
      double precision weight(maxn)

c++++ CPO
      double precision fw(nrec)

c++++ predictions
      double precision fw2(npred,ngrid)

c++++ credible bands
      double precision fs(ngrid)
      double precision fm(npred) 
      double precision worksam(nsave)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer count,i,i1,ii,j,j1,k
      integer sprint
      double precision tmp1

c+++++Distributions
      double precision dbet,dpoiss
      
c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++Slice sampler
      integer kk1
      double precision evaly,evalf,aux
      double precision gl
      double precision gr
      double precision g1

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c+++++RNG and distributions
      real runif

c++++ opening files

      open(unit=1,file='dppackage1.out',status='unknown',
     &     form='unformatted')
      open(unit=2,file='dppackage2.out',status='unknown',
     &     form='unformatted')


c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      cband=mcmc(4)
      tband=mcmc(5)

      betaw=slice(1)
      gw=slice(2)
      vw=slice(3)
      
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

      do i=1,p
         do j=1,p
            sbinv(i,j)=sb(i,j)
         end do
      end do
      call inverse(sbinv,p,iflag)      

      
      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()
      
c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Slice sampler step for beta, v, and g.
c+++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ Step a  

         call dbdplogposteri(kk,maxn,nrec,p,
     &                       alpha,beta,gp,v,lambda,
     &                       tau1,tau2,mub,sbinv,x,y,
     &                       workdpw,weight,aux)

         evaly=log(runif())+aux

c+++++++ Step b
         call dbdpstepbslice(0,gp,gw,gL,gR)
         do i=1,maxn
            do j=1,p
               call dbdpstepbslice(0,beta(i,j),betaw,
     &                             betal(i,j),betar(i,j))
            end do
         end do
         do i=1,(maxn-1)
            call dbdpstepbslice(1,v(i),vw,vL(i),vR(i))
         end do

c+++++++ Step c
         if(tau1.gt.0.d0)then
            call dbdpstepc1slice(gl,gr,g1)
           else 
            g1=gp
         end if

         do i=1,maxn
            do j=1,p
               call dbdpstepc1slice(betal(i,j),betar(i,j),beta1(i,j))
            end do
         end do

         do i=1,(maxn-1)
            call dbdpstepc1slice(vl(i),vr(i),v1(i))
         end do
         v1(maxn)=v(maxn)

         call dbdplogposteri(kk,maxn,nrec,p,
     &                       alpha,beta1,g1,v1,lambda,
     &                       tau1,tau2,mub,sbinv,x,y,
     &                       workdpw,weight,evalf)

         do while(evaly.gt.evalf)


c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            acrate(1)=acrate(1)+1.d0

            if(tau1.gt.0.d0)then
               call dbdpstepc2slice(0,gp,g1,gl,gr)
               call dbdpstepc1slice(gl,gr,g1)
              else
               g1=gp
            end if
            
            do i=1,maxn
               do j=1,p
                  call dbdpstepc2slice(0,beta(i,j),beta1(i,j),
     &                                 betal(i,j),betar(i,j))
                  call dbdpstepc1slice(betal(i,j),betar(i,j),
     &                                 beta1(i,j))
               end do
            end do
            do i=1,(maxn-1)
               call dbdpstepc2slice(1,v(i),v1(i),vl(i),vr(i))
               call dbdpstepc1slice(vl(i),vr(i),v1(i))
            end do
        
            call dbdplogposteri(kk,maxn,nrec,p,
     &                          alpha,beta1,g1,v1,lambda,
     &                          tau1,tau2,mub,sbinv,x,y,
     &                          workdpw,weight,evalf)

         end do
        
c+++++++ update parameters

         gp=g1
         do i=1,maxn
            do j=1,p
               beta(i,j)=beta1(i,j)
            end do
         end do
         do i=1,maxn
            v(i)=v1(i)
         end do
 
c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ sampling kk
c+++++++++++++++++++++++++++++++++++++++++++++++++
         if(lambda.gt.0.d0)then

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            kk1=kk+1
            if(kk.gt.1)then
               if(runif().gt.0.5)then
                 kk1=kk-1
               end if
            end if

            call dbdploglike(kk,maxn,nrec,p,beta,x,y,v,workdpw,
     &                       weight,evaly)

            call dbdploglike(kk1,maxn,nrec,p,beta,x,y,v,workdpw,
     &                       weight,evalf)

            aux=evalf-evaly+
     &          dpoiss(dble(kk1),dble(lambda),1)-
     &          dpoiss(dble(kk) ,dble(lambda),1)

            if(log(dble(runif())).lt.aux)then
               acrate(2)=acrate(2)+1.d0
               kk=kk1
            end if
         end if

c         call intpr("k",-1,kk,1)

c+++++++++++++++++++++++++++++++++++
c+++++++ baseline mean           +++
c+++++++++++++++++++++++++++++++++++

         do i=1,p
            workv1(i)=s0invm(i)
         end do

         do i=1,p
            do j=1,p
               workm1(i,j)=s0inv(i,j)+dble(maxn)*sbinv(i,j)
            end do
         end do
         call inverse(workm1,p,iflag)

         do ii=1,maxn      
            do i=1,p
               workv2(i)=beta(ii,i)
            end do
            workv1=workv1+matmul(sbinv,workv2)
         end do

         workv2=matmul(workm1,workv1)
         call rmvnorm(p,workv2,workm1,workmh1,workv1,mub)

c++++++++++++++++++++++++++++++++++++++
c+++++++ baseline covariance matrix +++
c++++++++++++++++++++++++++++++++++++++

         if(tau1.le.0.d0)then
            do i=1,p
               do j=1,p
                  sb(i,j)=psiinv(i,j)
                  workm1(i,j)=0.d0
               end do
               workv1(i)=0.d0
               iflag(i)=0
            end do
  
            do ii=1,maxn
               do i=1,p
                  do j=1,p
                     sb(i,j)=sb(i,j)+(beta(ii,i)-mub(i))*
     &                               (beta(ii,j)-mub(j))
                  end do
               end do
            end do

            call riwishart(p,nu+maxn,sb,sbinv,workm1,workv1,workmh1,
     &                     workmh2,iflag)

         end if

c         call dblepr("sb",-1,sb,p*p)
c         call dblepr("sbinv",-1,sbinv,p*p)

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

               count=0

c+++++++++++++ kk
               count=count+1
               thetasave(isave,count)=kk

c+++++++++++++ mub, sbinv
               do i=1,p
                  count=count+1
                  thetasave(isave,count)=mub(i)
               end do
               do i=1,p
                  do j=i,p
                     count=count+1
                     thetasave(isave,count)=sb(i,j)
                  end do 
               end do

c+++++++++++++ g-prior
               count=count+1
               thetasave(isave,count)=gp

c+++++++++++++ beta
               count=0
               do i=1,maxn
                  do j=1,p
                     count=count+1
                     randsave(isave,count)=beta(i,j)
                  end do
               end do

c+++++++++++++ stick-brealking weights
               do i=1,maxn
                  count=count+1
                  randsave(isave,count)=v(i)
               end do

c+++++++++++++ cpo
               call sickbreak(maxn,v,workdpw,weight)

               do i=1,nrec
                  fw(i)=0.d0
                  do j=1,maxn

                     aux=0.d0
                     do k=1,p
                        aux=aux+(x(i,k)*beta(j,k))
                     end do
                     tmp1=exp(aux)/(1.d0+exp(aux))

                     call jcomponentbd(tmp1,kk,k)

                     fw(i)=fw(i)+weight(j)*
     &                     dbet(y(i),dble(k),dble(kk-k+1),0)
                  end do
               end do

               do i=1,nrec
                   cpo(i,1)=cpo(i,1)+1.d0/fw(i)
                   cpo(i,2)=cpo(i,2)+fw(i)
               end do 
               
c+++++++++++++ predictions

               do i1=1,npred
                  fm(i1)=0.d0

                  do j1=1,ngrid
                      fw2(i1,j1)=0.d0

                      do j=1,maxn
                         aux=0.d0
                         do k=1,p
                            aux=aux+(xpred(i1,k)*beta(j,k))
                         end do
                         tmp1=exp(aux)/(1.d0+exp(aux))
                         call jcomponentbd(tmp1,kk,k)
                     
                         fw2(i1,j1)=fw2(i1,j1)+weight(j)*
     &                        dbet(grid(j1),dble(k),dble(kk-k+1),0)
                         
                         if(j1.eq.1)then 
                            fm(i1)=fm(i1)+weight(j)*dble(k)/dble(kk+1)
                         end if
                      end do
                  
                      fmean(i1,j1)=fmean(i1,j1)+fw2(i1,j1)
                  end do

                  meanfpm(i1)=meanfpm(i1)+fm(i1)

                  write(1) (fw2(i1,j1),j1=1,ngrid)
               end do

               write(2) (fm(i),i=1,npred)

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

      do i1=1,npred
         meanfpm(i1)=meanfpm(i1)/dble(nsave) 
         do j1=1,ngrid
            fmean(i1,j1)=fmean(i1,j1)/dble(nsave)
         end do
      end do 

      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do

      do i=1,2
         acrate(i)=acrate(i)/dble(nscan)
      end do

      close(unit=1)
      close(unit=2)

      if(cband.eq.1)then
         call hpddensreg(nsave,npred,ngrid,0.05d0,tband,worksam,fs,
     &                   flow,fupp)

         call hpddensregmf(nsave,npred,0.05d0,tband,worksam,
     &                     meanfpl,meanfph)

      end if

      return
      end



c=======================================================================
       subroutine dbdploglike(kk,maxn,nrec,p,beta,x,y,v,workv,weight,
     &                        eval)
c=======================================================================
c      log-Likelihood for the BDDP 
c=======================================================================
       implicit none

c+++++ Input
       integer kk,maxn,nrec,p
       double precision beta(maxn,p)
       double precision y(nrec),x(nrec,p),v(maxn)         

c+++++ External working space
       double precision workv(maxn+1),weight(maxn)

c+++++ Internal working space
       integer i,j,k
       double precision dbet,tmp1,tmp2

c+++++ Output
       double precision eval

c+++++ Algorithm

       call sickbreak(maxn,v,workv,weight)

       eval=0.d0
       do i=1,nrec
          tmp2=0.d0
          do j=1,maxn
             tmp1=0.d0
             do k=1,p
                tmp1=tmp1+(x(i,k)*beta(j,k))
             end do
             tmp1=exp(tmp1)/(1.d0+exp(tmp1))
             call jcomponentbd(tmp1,kk,k)
             tmp2=tmp2+weight(j)*dbet(y(i),dble(k),dble(kk-k+1),0)
          end do
          eval=eval+log(tmp2)
       end do

       return
       end
 

c=======================================================================
       subroutine dbdplogpriori(kk,maxn,p,
     &                          alpha,beta,gp,lambda,tau1,tau2,v,                 
     &                          mb,sbinv,eval)
c=======================================================================
c      log-priori distribution 
c=======================================================================
       implicit none

c+++++ Input
       integer kk,maxn,p
       double precision alpha,beta(maxn,p),gp,lambda,tau1,tau2,v(maxn)
       double precision mb(p),sbinv(p,p)

c+++++ Output 
       double precision eval

c+++++ Internal working space
       integer i,j,k
       double precision dpoiss,dgamm,dbet,tmp1

c+++++ Algorithm

c      kk
       eval=dpoiss(dble(kk),lambda,1)
       eval=eval-log(1.d0-dpoiss(0.d0,lambda,0))

c      g
       if(tau1.gt.0.d0)then
          eval=eval+dgamm(gp,tau1,tau2,1)       
       end if

c      v
       do i=1,maxn-1
          eval=eval+dbet(v(i),1.d0,alpha,1)
       end do

c      beta
       do i=1,maxn
          tmp1=0.d0
          do j=1,p
             do k=1,p
                tmp1=tmp1+(beta(i,j)-mb(j))*
     &                     sbinv(j,k)*
     &                    (beta(i,k)-mb(k))
             end do
          end do
          eval=eval-0.5d0*tmp1
       end do

       return
       end

c=======================================================================
       subroutine dbdplogposteri(kk,maxn,nrec,p,
     &                           alpha,beta,gp,v,lambda,
     &                           tau1,tau2,mb,sbinv,x,y,
     &                           workv,weight,eval)
c=======================================================================
c      log-posteriori distribution.
c=======================================================================
       implicit none

c+++++ Input
       integer kk,maxn,nrec,p
       double precision alpha
       double precision beta(maxn,p),gp,v(maxn)
       double precision lambda,tau1,tau2
       double precision mb(p)
       double precision sbinv(p,p)
       double precision x(nrec,p),y(nrec)

c+++++ External working space
       double precision workv(maxn+1),weight(maxn)

c+++++ Internal working space
       double precision tmp1

c+++++ Output
       double precision eval

c+++++ Algorithm

       call dbdploglike(kk,maxn,nrec,p,beta,x,y,v,workv,
     &                  weight,eval)

       call dbdplogpriori(kk,maxn,p,
     &                    alpha,beta,gp,lambda,tau1,tau2,v,                 
     &                    mb,sbinv,tmp1)

       eval=eval+tmp1

       return
       end


c=======================================================================
       subroutine dbdpstepbslice(limit,x0,w,l,r)
c=======================================================================
c      step b) slice sampler - Neal (2003) Fig 8
c=======================================================================
       implicit none

c+++++ Input
       integer limit
       double precision x0,w

c+++++ Output
       double precision l,r

c+++++ Internal working space
       real runif

c+++++ Algorithm

       l=x0-(w*dble(runif()))
       r=l+w

       if(limit.eq.1)then
          if(l.lt.0.d0)then
             l=0.d0   
          end if
          if(r.gt.1.d0)then
             r=1.d0
          end if
       end if

       return
       end

c=======================================================================
       subroutine dbdpstepc1slice(l,r,x1)
c=======================================================================
c      step b)-1 slice sampler - Neal (2003) Fig 8
c=======================================================================
       implicit none

c+++++ Input
       double precision l,r

c+++++ Output
       double precision x1

c+++++ Internal working space
       real runif

c+++++ Algorithm

       x1=l+(dble(runif())*(r-l))

       return
       end


c=======================================================================
       subroutine  dbdpstepc2slice(limit,x0,x1,l,r)
c=======================================================================
c      step b)-2 slice sampler - Neal (2003) Fig 8
c=======================================================================
       implicit none

c+++++ Input
       integer limit
       double precision x0,x1

c+++++ Output
       double precision l,r

c+++++ Algorithm

       if(x1.lt.x0)then
          l=x1
       end if

       if(x1.ge.x0)then
          r=x1
       end if

       if(limit.eq.1)then
          if(l.lt.0d0)then
             l=0.d0   
          end if
          if(r.gt.1.d0)then
             r=1.d0
          end if
       end if

       return
       end


