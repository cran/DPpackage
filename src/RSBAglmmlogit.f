c=======================================================================
      subroutine rsbaglmmlogit(maxni,nrec,p,x,y,
     &                         nsubject,subject,
     &                         datastr,
     &                         ngrid,grid,
     &                         prec,sb, 
     &                         maxn,maxr,maxc,sigma21,
     &                         tau1,taus1,taus2,tauh1,tauh2,
     &                         mcmc,nsave,
     &                         ss,beta,b,sba,sigma2,tau2,
     &                         acrate,cpo,densave,mc,thetasave,
     &                         randsave,sbasave,
     &                         seed,
     &                         ccluster,cstrt,possv,fw,prob,theta,w,
     &                         iflagp,betac,workmh1,workv1,
     &                         xtx,xty,
     &                         betasave,bsave)
c=======================================================================
c #  arguments
c=======================================================================
      implicit none

c+++++Data
      integer maxni,nrec,nsubject,p
      integer subject(nrec),y(nrec,2)
      integer datastr(nsubject,maxni+1)
      real*8 x(nrec,p)

c+++++Density
      integer ngrid
      real*8 grid(ngrid)
 
c+++++Prior
      integer maxn,maxr,maxc
      real*8 prec(p,p),sb(p,2)
      real*8 mu1,sigma21
      real*8 tau1
      real*8 taus1,taus2
      real*8 tauh1,tauh2

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Current value of the parameters
      integer ss(nsubject)
      real*8 beta(p)
      real*8 b(nsubject)
      real*8 sba(maxn,maxr)
      real*8 sigma2(maxc)
      real*8 tau2

c+++++Output
      real*8 acrate
      real*8 cpo(nrec,2)
      real*8 densave(nsave,ngrid) 
      real*8 mc(5)
      real*8 randsave(nsave,nsubject)
      real*8 sbasave(nsave,maxr+maxc)
      real*8 thetasave(nsave,p+2)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++External working space
      integer ccluster(maxc)
      integer cstrt(maxc,nsubject)
      integer possv(maxn)
      real*8 fw(maxn-1,maxc-1)
      real*8 prob(maxc)
      real*8 theta(maxc)
      real*8 w(maxc)
 
      integer iflagp(p)
      real*8 betac(p)
      real*8 workmh1(p*(p+1)/2)
      real*8 workv1(p)
      real*8 xtx(p,p)
      real*8 xty(p)

      real*8 betasave(p)
      real*8 bsave(nsubject)

c+++++Internal working space
      integer i,j,k,count
      integer dispcount
      integer iscan
      integer isave
      integer nscan
      integer ni,nij
      integer since
      integer skipcount
      integer sprint
      integer yij 

      real*8 dbin,dnrm
      real*8 lratio

      real*8 mu0,sigma20

      real*8 eta,gprime
      real*8 mean
      real*8 logcgkn,logcgko
      real*8 loglikn,logliko
      real*8 logpriorn,logprioro
      real*8 offset
      real*8 ratio
      real*8 tmp1,tmp2,tmp3
      real*8 ytilde
 
      real runif

c++++ Working space slice sampling
      integer evali
      real*8 rexpo,re,uwork
      real*8 logy,xx0,xx1,llim,rlim
      real*8 grlim,gllim,gxx0,gxx1

c++++ model's performance
      real*8 dbarc,dbar,dhat,pd,lpml

c+++++CPU time
      real*8 sec00,sec0,sec1,sec

c++++++++++++++++++++++++++
c     initialize variables
c++++++++++++++++++++++++++

c++++ rsba
      mu0=0.d0
      sigma20=0.d0
      mu1=0.d0

c++++ mcmc, priors and "zipped"

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
   
      seed1=seed(1)
      seed2=seed(2)

c++++ set random number generator

      call setall(seed1,seed2)
      
c++++ start the MCMC algorithm

      dbar=0.d0
      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*nsave

      call cpu_time(sec0)
      sec00=0.0

c++++ cluster structure

      do i=1,nsubject
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
      end do

c++++ compute weights and atoms from SBA

      call discrsba(maxn,maxr,maxc,sba,fw,theta,w)

c++++ mcmc

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++
c+++++++ fixed effects        +++
c++++++++++++++++++++++++++++++++

         do i=1,p
            do j=1,p
               xtx(i,j)=prec(i,j)
            end do
            xty(i)=sb(i,1)
         end do

         logliko=0.d0

         do i=1,nrec
            eta=0.d0
            offset=0.d0
            mean=0.d0
            gprime=0.d0
            tmp1=0.d0

            yij=y(i,1)
            nij=y(i,2)
            
            do j=1,p
               eta=eta+x(i,j)*beta(j)
            end do
               
            eta=eta+b(subject(i)) 
            offset=b(subject(i)) 

            mean=exp(eta)/(1.d0+exp(eta))
            logliko=logliko+dbin(dble(yij),dble(nij),mean,1)

            tmp1=mean*(1.0d0-mean)
            gprime=1.d0/(dble(nij)*tmp1)

            mean=dble(nij)*exp(eta)/(1.d0+exp(eta))
            ytilde=eta+(dble(yij)-mean)*gprime-offset

            do j=1,p
               do k=1,p
                  xtx(j,k)=xtx(j,k)+x(i,j)*x(i,k)/gprime
               end do
               xty(j)=xty(j)+x(i,j)*ytilde/gprime
            end do
         end do

         call inverse(xtx,p,iflagp)      
         do i=1,p
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+xtx(i,j)*xty(j) 
            end do
            workv1(i)=tmp1
         end do

         call rmvnorm(p,workv1,xtx,workmh1,xty,betac)

c         call dblepr("betac",-1,betac,p)

c+++++++ evaluating the candidate generating kernel

         call dmvnd(p,betac,workv1,xtx,logcgko,iflagp)

c+++++++ evaluating the likelihood

         do i=1,p
            do j=1,p
               xtx(i,j)=prec(i,j)
            end do
            xty(i)=sb(i,1)
         end do

         loglikn=0.d0

         do i=1,nrec
            eta=0.d0
            offset=0.d0
            mean=0.d0
            gprime=0.d0
            tmp1=0.d0

            yij=y(i,1)
            nij=y(i,2)
            
            do j=1,p
               eta=eta+x(i,j)*betac(j)
            end do
               
            eta=eta+b(subject(i)) 
            offset=b(subject(i)) 
               
            mean=exp(eta)/(1.d0+exp(eta))
            loglikn=loglikn+dbin(dble(yij),dble(nij),mean,1)

            tmp1=mean*(1.0d0-mean)
            gprime=1.d0/(dble(nij)*tmp1)

            mean=dble(nij)*exp(eta)/(1.d0+exp(eta))
            ytilde=eta+(dble(yij)-mean)*gprime-offset

            do j=1,p
               do k=1,p
                  xtx(j,k)=xtx(j,k)+x(i,j)*x(i,k)/gprime
               end do
               xty(j)=xty(j)+x(i,j)*ytilde/gprime
            end do
         end do

         call inverse(xtx,p,iflagp)      
         do i=1,p
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+xtx(i,j)*xty(j) 
            end do
            workv1(i)=tmp1
         end do

c+++++++ evaluating the candidate generating kernel

         call dmvnd(p,beta,workv1,xtx,logcgkn,iflagp)

c+++++++ prior ratio
         logprioro=0.d0
         logpriorn=0.d0
         
         do i=1,p
            do j=1,p
               logpriorn=logpriorn+(betac(i)-sb(i,2))* 
     &                      prec(i,j)       *
     &                      (betac(j)-sb(j,2))

               logprioro=logprioro+(beta(i)-sb(i,2))* 
     &                      prec(i,j)      *
     &                      (beta(j)-sb(j,2))

            end do
         end do
      
         logpriorn=-0.5d0*logpriorn
         logprioro=-0.5d0*logprioro

c+++++++ mh step

         ratio=loglikn-logliko+logcgkn-logcgko+
     &         logpriorn-logprioro

         if(log(dble(runif())).lt.ratio)then
            acrate=acrate+1.d0
            do i=1,p
               beta(i)=betac(i) 
            end do
         end if

c         call dblepr("beta",-1,beta,p)
c         go to 9999

c++++++++++++++++++++++++++++++++++
c+++++++ Update random effects  +++
c++++++++++++++++++++++++++++++++++

         do i=1,nsubject

            evali=1
            xx0=b(i)
            re=rexpo(1.d0)
            logy=-re

            call rsbalpostlogit(i,
     &                          maxc,maxni,nrec,nsubject,p,
     &                          datastr,y,
     &                          beta,b,x,
     &                          theta,sigma2,w,tmp1)

            logy=logy+tmp1
            uwork=dble(runif())*0.25d0  
            llim=xx0-uwork
            rlim=llim+0.25d0

            evali=evali+1
            b(i)=llim
            call rsbalpostlogit(i,
     &                          maxc,maxni,nrec,nsubject,p,
     &                          datastr,y,
     &                          beta,b,x,
     &                          theta,sigma2,w,gllim)


            evali=evali+1
            b(i)=rlim
            call rsbalpostlogit(i,
     &                          maxc,maxni,nrec,nsubject,p,
     &                          datastr,y,
     &                          beta,b,x,
     &                          theta,sigma2,w,grlim)

            do while(gllim.gt.logy)
               llim=llim-0.25d0

               b(i)=llim
               call rsbalpostlogit(i,
     &                             maxc,maxni,nrec,nsubject,p,
     &                             datastr,y,
     &                             beta,b,x,
     &                             theta,sigma2,w,gllim)

            end do 

            do while(grlim.gt.logy)
               rlim=rlim+0.25d0

               b(i)=rlim
               call rsbalpostlogit(i,
     &                             maxc,maxni,nrec,nsubject,p,
     &                             datastr,y,
     &                             beta,b,x,
     &                             theta,sigma2,w,grlim)

            end do 

            xx1=llim+(rlim-llim)*dble(runif())
            evali=evali+1
            b(i)=xx1
            call rsbalpostlogit(i,
     &                          maxc,maxni,nrec,nsubject,p,
     &                          datastr,y,
     &                          beta,b,x,
     &                          theta,sigma2,w,gxx1)

            do while(gxx1.lt.logy)
               xx1=llim+(rlim-llim)*dble(runif())
               b(i)=xx1

               call rsbalpostlogit(i,
     &                             maxc,maxni,nrec,nsubject,p,
     &                             datastr,y,
     &                             beta,b,x,
     &                             theta,sigma2,w,gxx1)

            end do
         end do

c         call dblepr("b",-1,b,nsubject)

c++++++++++++++++++++++++++++++++++
c+++++++ clustering structure   +++
c++++++++++++++++++++++++++++++++++

         call rsbauclus(nsubject,b,maxn,maxc,w,theta,sigma2,prob,
     &                  ccluster,cstrt,ss)

c         call intpr("ss",-1,ss,nsubject)

c++++++++++++++++++++++++++++++++++
c+++++++ barycenter parameters  +++
c++++++++++++++++++++++++++++++++++

         call rsbauba(nsubject,b,
     &                maxn,maxr,maxc,
     &                mu0,mu1,sigma20,sigma21,
     &                ccluster,cstrt,sigma2,
     &                fw,possv,
     &                sba,theta,w)

c         call dblepr("theta",-1,theta,maxc)
c         call dblepr("w",-1,w,maxc)

c++++++++++++++++++++++++++++++++++
c+++++++ kernel variances       +++
c++++++++++++++++++++++++++++++++++

         call rsbaukvar(nsubject,b,maxn,maxc,ccluster,cstrt,
     &                  ss,theta,
     &                  tau1,tau2,
     &                  sigma2)

c         call dblepr("sigma2",-1,sigma2,maxc)

c++++++++++++++++++++++++++++++++++
c+++++++ gamma parameter        +++
c++++++++++++++++++++++++++++++++++

         call rsbaugamma(maxc,sigma2,
     &                   tau1,taus1,taus2,
     &                   tau2)

c         call dblepr("tau2",-1,tau2,1)


c++++++++++++++++++++++++++++++++++
c+++++++ SS Update sd1          +++
c++++++++++++++++++++++++++++++++++

         if(tauh1.gt.0.d0)then

            evali=1
            xx0=sigma21
            re=rexpo(1.d0)
            logy=-re
            call lposth1ssba(maxn,maxr,mu1,
     &                       sigma21,sba,tauh1,
     &                       tauh2,tmp1)
            logy=logy+tmp1
            uwork=dble(runif())*0.25d0  
            llim=xx0-uwork
            rlim=llim+0.25d0

            if(llim.lt.0.01d0)llim=0.01d0 

            evali=evali+1
            tmp1=llim
            call lposth1ssba(maxn,maxr,mu1,
     &                       tmp1,sba,tauh1,
     &                       tauh2,gllim)

            evali=evali+1
            tmp1=rlim
            call lposth1ssba(maxn,maxr,mu1,
     &                       tmp1,sba,tauh1,
     &                       tauh2,grlim)

            do while(gllim.gt.logy)
               llim=llim-0.25d0
               if(llim.lt.0.01d0)then
                  llim=0.01d0
                  gllim=logy-1.d0
                 else
                  evali=evali+1
                  tmp1=llim
                  call lposth1ssba(maxn,maxr,
     &                           mu1,tmp1,sba,
     &                           tauh1,tauh2,
     &                           gllim)
               end if
            end do 

            do while(grlim.gt.logy)
               rlim=rlim+0.25d0
               tmp1=rlim
               call lposth1ssba(maxn,maxr,mu1,
     &                          tmp1,sba,tauh1,
     &                          tauh2,grlim)
            end do 

            xx1=llim+(rlim-llim)*dble(runif())
            evali=evali+1
            tmp1=xx1
            call lposth1ssba(maxn,maxr,mu1,
     &                       tmp1,sba,tauh1,
     &                       tauh2,gxx1)


            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1
               if(llim.lt.0.01d0)llim=0.01d0 

               evali=evali+1
               xx1=llim+(rlim-llim)*dble(runif())
               tmp1=xx1
               call lposth1ssba(maxn,maxr,mu1,
     &                          tmp1,sba,tauh1,
     &                          tauh2,gxx1)
            end do

            sigma21=xx1

c            call dblepr("sigma21",-1,sigma21,1)

         end if


c+++++++++++++++++++++++++++++++++++
c+++++++ saving samples          +++
c+++++++++++++++++++++++++++++++++++

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ fixed effects
               do i=1,p
                  thetasave(isave,i)=beta(i) 
                  betasave(i)=betasave(i)+beta(i)
               end do

c+++++++++++++ random effecs
               do i=1,nsubject
                  randsave(isave,i)=b(i)
                  bsave(i)=bsave(i)+b(i)
               end do

c+++++++++++++ sba parameters
               count=0
               do i=1,maxr
                  count=count+1
                  sbasave(isave,count)=sba(maxn,i)
               end do          

c+++++++++++++ kernel variances
               do i=1,maxc
                  count=count+1
                  sbasave(isave,count)=sigma2(i)
               end do          

c+++++++++++++ gamma parameter
               thetasave(isave,p+1)=tau2

c+++++++++++++ H1 parameter
               thetasave(isave,p+2)=sigma21

c+++++++++++++ density
               do i=1,ngrid
                  tmp1=0.d0
                  do j=1,maxc
                     tmp1=tmp1+w(j)*dnrm(grid(i),theta(j),
     &                                   sqrt(sigma2(j)),0)
                  end do
                  densave(isave,i)=tmp1
               end do
     
c+++++++++++++ cpo's

               dbarc=0.d0
               do i=1,nrec
                  yij=y(i,1)
                  nij=y(i,2)
                  eta=0.d0
                  do j=1,p
                     eta=eta+x(i,j)*beta(j)
                  end do
                  eta=eta+b(subject(i))
                  mean=exp(eta)/(1.d0+exp(eta))

                  tmp1=dbin(dble(yij),dble(nij),mean,0)
                  cpo(i,1)=cpo(i,1)+1.0d0/tmp1
                  cpo(i,2)=cpo(i,2)+tmp1
                  
                  tmp1=dbin(dble(yij),dble(nij),mean,1)
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

      acrate=acrate/dble(nscan)

      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)
      end do
            
      do i=1,p
         betasave(i)=betasave(i)/dble(nsave)
      end do

      do i=1,nsubject
         bsave(i)=bsave(i)/dble(nsave)
      end do   

      dhat=0.d0
      lpml=0.d0
      do i=1,nrec
         yij=y(i,1)
         nij=y(i,2)
         eta=0.d0
         do j=1,p
            eta=eta+x(i,j)*betasave(j)
         end do
         eta=eta+bsave(subject(i))
         mean=exp(eta)/(1.d0+exp(eta))
         dhat=dhat+dbin(dble(yij),dble(nij),mean,1)
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


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
      subroutine rsbalpostlogit(ind,
     &                          maxc,maxni,nrec,nsubject,p,
     &                          datastr,y,
     &                          beta,b,x,
     &                          theta,sigma2,w,lpost)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
      implicit none

c++++ Input
      integer ind
      integer maxc,maxni,nrec,nsubject,p
      integer datastr(nsubject,maxni+1)
      integer y(nrec,2)
      real*8 beta(p)
      real*8 b(nsubject)
      real*8 x(nrec,p)
      real*8 theta(maxc),sigma2(maxc),w(maxc)

c++++ Output
      real*8 lpost   

c++++ Internal working space 
      integer i,j,jj,k,nn,nij,yij
      real*8 dnrm,tmp1
      real*8 dbin,mean

      lpost=0.d0

c++++ likelihood contribution      
      nn=datastr(ind,1) 
      do j=1,nn

         jj=j+1

         yij=y(datastr(ind,jj),1)            
         nij=y(datastr(ind,jj),2)            

         tmp1=0
         do k=1,p
            tmp1=tmp1+x(datastr(ind,jj),k)*beta(k)
         end do
         tmp1=tmp1+b(ind)
         mean=exp(tmp1)/(1.d0+exp(tmp1))
         lpost=lpost+dbin(dble(yij),dble(nij),mean,1)
      end do

c++++ prior contribution      
      tmp1=0.d0
      do j=1,maxc
         tmp1=tmp1+w(j)*dnrm(b(ind),theta(j),
     &                       sqrt(sigma2(j)),0)
      end do
      lpost=lpost+log(tmp1)

      return
      end

