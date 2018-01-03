c=======================================================================
       subroutine ldtfpglmmpois(maxni,nrec,nsubject,p,
     &                          subject,datastr,x,y,
     &                          ptf,xtf,
     &                          ngrid,npred,grid,
     &                          xpred,xtfpred,
     &                          quans,
     &                          prec,sb,
     &                          maxm,ntprob,ntlr,a0b0,
     &                          gprior,taub,
     &                          as,b,
     &                          beta,betatf,
     &                          mcmc,nsave,seed,
     &                          cpo,densm,densl,densu,
     &                          qmm,qml,qmu,
     &                          thetasave,randsave,tfpsave,
     &                          betac,workmp1,workmp2,
     &                          iflag,iflagtf,
     &                          nobsbc,obsbc,c0,
     &                          workm1,workvh1,workv1,workvp1,
     &                          worksam,worksam2,worksam3,fs,
     &                          workm2,workvh2,
     &                          workv2,workv3,workv4,k,
     &                          prob,probc)
c=======================================================================
c     # 65 arguments
c
c     Subroutine `ldtfpglmmpos' to run a Markov chain for a
c     semiparametric GLMM when a linear dependent tailfree process 
c     prior is used for the random effects. A Poisson conditional
c     distribution is fit in this functions.
c
c     Copyright: Alejandro Jara, 2012.
c
c     Version 1.0: 
c
c     Last modification: 26-06-2012.
c     
c=======================================================================                  

      implicit none
c++++ input - glmm

      integer maxni,nrec,nsubject,p
      integer subject(nrec)
      integer datastr(nsubject,maxni+1)
      integer y(nrec)

      double precision x(nrec,p+1) 

c++++ input - ldtfp
      integer ptf
      double precision xtf(nsubject,ptf)

c++++ prediction
      integer ngrid,npred
      double precision grid(ngrid)
      double precision xpred(npred,p)
      double precision xtfpred(npred,ptf)
      double precision quans(3)

c++++ prior - glmm
      double precision prec(p,p)
      double precision sb(p,2)

c++++ prior - ldtfp
      integer maxm,ntprob,ntlr
      double precision a0b0(2)
      double precision gprior(ptf,ptf)
      double precision taub(2),taub1,taub2

c++++ current value
      double precision as(3)
      double precision alpha
      double precision b(nsubject)
      double precision beta(p)
      double precision betatf(ntlr,ptf)
      double precision sigma2b

c++++ mcmc
      integer mcmc(5),nburn,nskip,nsave,ndisplay

c++++ bands
      integer cband,tband

c++++ seeds
      integer seed1,seed2,seed(2)

c++++ output
      double precision cpo(nrec,2)
      double precision densm(npred,ngrid)
      double precision densl(npred,ngrid)
      double precision densu(npred,ngrid)

      double precision qmm(npred,3)
      double precision qml(npred,3)
      double precision qmu(npred,3)

      double precision thetasave(nsave,p+2)
      double precision randsave(nsave,nsubject)
      double precision tfpsave(nsave,(ntlr-1)*ptf)

c++++ external working space - vector and matrices
      integer iflagtf(ptf)
      integer nobsbc(ntprob)
      integer obsbc(ntprob,nsubject)
      double precision c0(ptf,ptf)

      integer iflag(p)
      double precision betac(p)
      double precision workm1(p,p)
      double precision workmp1(p,p)
      double precision workmp2(p,p)
      double precision workvh1(p*(p+1)/2)
      double precision workv1(p)
      double precision workvp1(p)

      double precision worksam(nsave)
      double precision worksam2(nsave,ngrid)
      double precision worksam3(nsave,npred)
      double precision fs(ngrid)

      double precision workm2(ptf,ptf)
      double precision workvh2(ptf*(ptf+1)/2)
      double precision workv2(ptf)
      double precision workv3(ptf)
      double precision workv4(ptf)

      integer k(maxm)

      double precision prob(2**maxm)
      double precision probc(2**maxm)
 
c++++ internal working space
      integer accept
      integer i,ii,i1,j,jj,j1,j2,k1,ll,mm
      integer d1,d2
      integer dispcount,iscan,isave
      integer ntot,n1,n2
      integer nscan,skipcount
      integer sprint
      double precision acrate
      double precision loglikn,logliko
      double precision logpriorn,logprioro
      double precision logcgkn,logcgko
      double precision ratio
      double precision rgamma
      real runif
      double precision tmp1,tmp2
      double precision uni
  
      integer yij
      double precision dpoiss
      double precision eta
      double precision offset
      double precision mean
      double precision gprime
      double precision ytilde

c++++ CPU time
      double precision sec00,sec0,sec1,sec

c++++ Working space slice sampling
      integer evali
      double precision re,uwork
      double precision logy,xx0,xx1,llim,rlim
      double precision grlim,gllim,gxx1

c++++ initializing variables

      open(unit=1,file='dppackage_dftp1.out',status='unknown',
     &     form='unformatted')
      open(unit=2,file='dppackage_dftp2.out',status='unknown',
     &     form='unformatted')
      open(unit=3,file='dppackage_dftp3.out',status='unknown',
     &     form='unformatted')
      open(unit=4,file='dppackage_dftp4.out',status='unknown',
     &     form='unformatted')
      open(unit=5,file='dppackage_dftp5.out',status='unknown',
     &     form='unformatted')
 

      alpha=as(1)
      sigma2b=as(2)

      taub1=taub(1)
      taub2=taub(2)  

      call loglikldtfpre(maxm,ntlr,ntprob,nsubject,ptf,betatf,b,
     &                   xtf,sigma2b,nobsbc,obsbc,logliko,k)

      i1=1
      do i=2,maxm
         j1=2**(i-1)

         do j=1,j1
            k1=i1+j

            c0=gprior/(alpha*dble(i**2))

            ii=(k1-1)*2+1
            jj=(k1-1)*2+2
            n1=nobsbc(ii)
            n2=nobsbc(jj)
            ntot=n1+n2

            if(n1.gt.0.and.n2.gt.0)then

               call startlrcoefldtfp(5,k1,ii,jj,n1,n2,maxm,ntlr,
     &                               ntprob,nsubject,ptf,obsbc,betatf,
     &                               xtf,c0,iflagtf,workv2,workm2,
     &                               workv3)
            end if
         end do
         i1=i1+j1
      end do

      call loglikldtfpre(maxm,ntlr,ntprob,nsubject,ptf,betatf,b,
     &                   xtf,sigma2b,nobsbc,obsbc,logliko,k)


      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
   
      cband=mcmc(4)
      tband=mcmc(5)

      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)

c++++ MCMC
      acrate=0.d0
      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*nsave

      call cpu_time(sec0)
      sec00=0.d0

      do iscan=1,nscan
 
c++++++++++++++++++++++++++++++++++++++
c+++++++ fixed effects  +++
c++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,p
            do j=1,p
               workm1(i,j)=prec(i,j) 
            end do
            workv1(i)=sb(i,1)
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
               
            eta=eta+b(subject(i)) 
            offset=b(subject(i)) 

            eta=eta+x(i,p+1)
            offset=offset+x(i,p+1)

            mean=exp(eta)
            gprime=exp(-eta)

            ytilde=eta+(dble(yij)-mean)*gprime-offset
               
            do j1=1,p
               do j2=1,p
                  workm1(j1,j2)=workm1(j1,j2)+
     &                          x(i,j1)*x(i,j2)*mean
               end do
               workv1(j1)=workv1(j1)+
     &                    x(i,j1)*ytilde*mean
            end do
               
            logliko=logliko+dpoiss(dble(yij),mean,1)
         end do

         call inverse(workm1,p,iflag)      

         do i=1,p
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+workm1(i,j)*workv1(j) 
            end do
            workvp1(i)=tmp1
         end do

         call rmvnorm(p,workvp1,workm1,workvh1,workv1,betac)

c+++++++ evaluating the candidate generating kernel

         call dmvn2(p,betac,workvp1,workm1,logcgko,
     &              workv1,workmp1,workmp2,iflag)                 


c+++++++ evaluating the likelihood

         do i=1,p
            do j=1,p
               workm1(i,j)=prec(i,j) 
            end do
            workv1(i)=sb(i,1)
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
               
            eta=eta+b(subject(i)) 
            offset=b(subject(i)) 

            eta=eta+x(i,p+1)
            offset=offset+x(i,p+1)

            mean=exp(eta)
            gprime=exp(-eta)

            ytilde=eta+(dble(yij)-mean)*gprime-offset
               
            do j1=1,p
               do j2=1,p
                  workm1(j1,j2)=workm1(j1,j2)+
     &                          x(i,j1)*x(i,j2)*mean
               end do
               workv1(j1)=workv1(j1)+
     &                    x(i,j1)*ytilde*mean
            end do
               
            loglikn=loglikn+dpoiss(dble(yij),mean,1)
         end do

         call inverse(workm1,p,iflag)      

         do i=1,p
            tmp1=0.d0
            do j=1,p
               tmp1=tmp1+workm1(i,j)*workv1(j) 
            end do
            workvp1(i)=tmp1
         end do

c+++++++ evaluating the candidate generating kernel

         call dmvn2(p,beta,workvp1,workm1,logcgkn,
     &              workv1,workmp1,workmp2,iflag)                 


c+++++++ prior ratio
         logprioro=0.d0
         logpriorn=0.d0
         
         do i=1,p
            do j=1,p
               logpriorn=logpriorn+(betac(i)-sb(i,2))* 
     &                      prec(i,j)       *
     &                   (betac(j)-sb(j,2))

               logprioro=logprioro+(beta(i)-sb(i,2))* 
     &                      prec(i,j)      *
     &                   (beta(j)-sb(j,2))

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

c+++++++++++++++++++++++++++++++++
c+++++++ random effects        +++ 
c+++++++++++++++++++++++++++++++++

         do i=1,nsubject


c++++++++++ check if the user has requested an interrupt
            call rchkusr()
       
            evali=1
            xx0=b(i)

            call lpostdtfprepoi(xx0,i,
     &                          maxm,maxni,nrec,nsubject,
     &                          ntlr,ntprob,p,
     &                          datastr,y,
     &                          beta,x,ptf,betatf,xtf,
     &                          sigma2b,tmp1,k)


            uni=runif()
            re=-log(uni)
            logy=tmp1-re

            uwork=dble(runif())*0.25d0  
            llim=xx0-uwork
            rlim=xx0+(0.25d0-uwork)

            evali=evali+1
            call lpostdtfprepoi(llim,i,
     &                          maxm,maxni,nrec,nsubject,
     &                          ntlr,ntprob,p,
     &                          datastr,y,
     &                          beta,x,ptf,betatf,xtf,
     &                          sigma2b,gllim,k)

            evali=evali+1
            call lpostdtfprepoi(rlim,i,
     &                          maxm,maxni,nrec,nsubject,
     &                          ntlr,ntprob,p,
     &                          datastr,y,
     &                          beta,x,ptf,betatf,xtf,
     &                          sigma2b,grlim,k)

            do while(gllim.gt.logy)
               llim=llim-0.25d0
               evali=evali+1
               call lpostdtfprepoi(llim,i,
     &                          maxm,maxni,nrec,nsubject,
     &                          ntlr,ntprob,p,
     &                          datastr,y,
     &                          beta,x,ptf,betatf,xtf,
     &                          sigma2b,gllim,k)

            end do 

            do while(grlim.gt.logy)
               rlim=rlim+0.25d0
               evali=evali+1
               call lpostdtfprepoi(rlim,i,
     &                          maxm,maxni,nrec,nsubject,
     &                          ntlr,ntprob,p,
     &                          datastr,y,
     &                          beta,x,ptf,betatf,xtf,
     &                          sigma2b,grlim,k)
            end do 

            xx1=llim+(rlim-llim)*dble(runif())

            evali=evali+1
            call lpostdtfprepoi(xx1,i,
     &                          maxm,maxni,nrec,nsubject,
     &                          ntlr,ntprob,p,
     &                          datastr,y,
     &                          beta,x,ptf,betatf,xtf,
     &                          sigma2b,gxx1,k)

            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1

               xx1=llim+(rlim-llim)*dble(runif())
               evali=evali+1
               call lpostdtfprepoi(xx1,i,
     &                          maxm,maxni,nrec,nsubject,
     &                          ntlr,ntprob,p,
     &                          datastr,y,
     &                          beta,x,ptf,betatf,xtf,
     &                          sigma2b,gxx1,k)

            end do
            b(i)=xx1

         end do

c         call dblepr("b",-1,b,nsubject)


c++++++++++++++++++++++++++++++++++++++
c+++++++ baseline variance          +++
c++++++++++++++++++++++++++++++++++++++


c+++++++ check if the user has requested an interrupt
         call rchkusr()

         evali=1
         xx0=sigma2b

         call loglikldtfpre(maxm,ntlr,ntprob,nsubject,
     &                      ptf,betatf,b,
     &                      xtf,xx0,nobsbc,obsbc,loglikn,k)

         if(taub1>0)then
            logpriorn=-(0.5d0*taub1+1.d0)*log(xx0)-0.5d0*taub2/xx0
          else
            logpriorn=-log(xx0)
         end if      
         tmp1=loglikn+logpriorn

         uni=runif()
         re=-log(uni)
         logy=tmp1-re

         uwork=dble(runif())*0.05d0  
         llim=xx0-uwork
         rlim=xx0+(0.05d0-uwork)
           
         if(llim.lt.0.00001d0)llim=0.00001d0 

         evali=evali+1
         call loglikldtfpre(maxm,ntlr,ntprob,nsubject,
     &                      ptf,betatf,b,
     &                      xtf,llim,nobsbc,obsbc,loglikn,k)
         if(taub1>0)then
            logpriorn=-(0.5d0*taub1+1.d0)*log(llim)-
     &                  0.5d0*taub2/llim
          else
            logpriorn=-log(llim)
         end if      
         gllim=loglikn+logpriorn

         evali=evali+1
         call loglikldtfpre(maxm,ntlr,ntprob,nsubject,
     &                      ptf,betatf,b,
     &                      xtf,rlim,nobsbc,obsbc,loglikn,k)
         if(taub1>0)then
            logpriorn=-(0.5d0*taub1+1.d0)*log(rlim)-
     &                  0.5d0*taub2/rlim
          else
            logpriorn=-log(rlim)
         end if      
         grlim=loglikn+logpriorn

         do while(gllim.gt.logy)
            llim=llim-0.05d0

            if(llim.lt.0.00001d0)then
               llim=0.00001d0 
               gllim=logy-1.d0
              else   
               evali=evali+1

               call loglikldtfpre(maxm,ntlr,ntprob,nsubject,
     &                      ptf,betatf,b,
     &                      xtf,llim,nobsbc,obsbc,loglikn,k)
               if(taub1>0)then
                  logpriorn=-(0.5d0*taub1+1.d0)*log(llim)-
     &                        0.5d0*taub2/llim
                 else
                  logpriorn=-log(llim)
               end if      
               gllim=loglikn+logpriorn
             end if
         end do 

         do while(grlim.gt.logy)
            rlim=rlim+0.05d0

            evali=evali+1
            call loglikldtfpre(maxm,ntlr,ntprob,nsubject,
     &                      ptf,betatf,b,
     &                      xtf,rlim,nobsbc,obsbc,loglikn,k)

            if(taub1>0)then
               logpriorn=-(0.5d0*taub1+1.d0)*log(rlim)-
     &                     0.5d0*taub2/rlim
             else
               logpriorn=-log(rlim)
            end if      
            grlim=loglikn+logpriorn
         end do 

         xx1=llim+(rlim-llim)*dble(runif())
 
         evali=evali+1
         call loglikldtfpre(maxm,ntlr,ntprob,nsubject,
     &                      ptf,betatf,b,
     &                      xtf,xx1,nobsbc,obsbc,loglikn,k)

         if(taub1>0)then
           logpriorn=-(0.5d0*taub1+1.d0)*log(xx1)-
     &                 0.5d0*taub2/xx1
         else
           logpriorn=-log(xx1)
        end if      
        gxx1=loglikn+logpriorn

        do while(gxx1.lt.logy)
           if(xx1.gt.xx0)rlim=xx1
           if(xx1.lt.xx0)llim=xx1

           if(llim.lt.0.00001d0)llim=0.00001d0 

           xx1=llim+(rlim-llim)*dble(runif())

           evali=evali+1

           call loglikldtfpre(maxm,ntlr,ntprob,nsubject,
     &                      ptf,betatf,b,
     &                      xtf,xx1,nobsbc,obsbc,loglikn,k)

           if(taub1>0)then
              logpriorn=-(0.5d0*taub1+1.d0)*log(xx1)-
     &                    0.5d0*taub2/xx1
            else
              logpriorn=-log(xx1)
           end if
           gxx1=loglikn+logpriorn
        end do

        sigma2b=xx1

c        call dblepr("sigma2b",-1,sigma2b,1)


c+++++++++++++++++++++++++++++++++++++
c++++++ tf logistic regressions    +++
c+++++++++++++++++++++++++++++++++++++

        i1=1
        do i=2,maxm
           j1=2**(i-1)

           do j=1,j1

c++++++++++++ check if the user has requested an interrupt
              call rchkusr()

              k1=i1+j

              do d1=1,ptf
                 do d2=1,ptf 
                    c0(d1,d2)=gprior(d1,d2)/
     &                        (alpha*real(i**2))
                 end do
              end do
              ii=(k1-1)*2+1
              jj=(k1-1)*2+2
              n1=nobsbc(ii)
              n2=nobsbc(jj)
              ntot=n1+n2

              if(ntot.gt.0)then

                 call updatelrcoefldtfpss(k1,ii,jj,n1,n2,maxm,
     &                            ntlr,ntprob,
     &                            nsubject,ptf,obsbc,betatf,xtf,c0,
     &                            accept,iflagtf,workv2,workm2)

               else
                 call updatelrcoefldtfp0(k1,maxm,ntlr,ntprob,
     &                                   ptf,betatf,c0,
     &                                   iflagtf,workv2,workv3,
     &                                   workm2,workvh2,workv4)
              end if

           end do
           i1=i1+j1
        end do


c+++++++++++++++++++++++++++++++++++++
c++++++ alpha                      +++
c+++++++++++++++++++++++++++++++++++++
     
        if(a0b0(1)>0)then

c+++++++++ check if the user has requested an interrupt
           call rchkusr()

           i1=1
           tmp1=0.d0

           do i=2,maxm
              j1=2**(i-1)
        
              do j=1,j1
                 k1=i1+j

                 do d1=1,ptf
                    do d2=1,ptf 
                       c0(d1,d2)=gprior(d1,d2)/
     &                        (dble(i**2))
                    end do
                 end do

                 call inverse(c0,ptf,iflagtf)
                 tmp2=0.d0 
                 do ll=1,ptf
                    do mm=1,ptf
                       tmp2=tmp2+betatf(k1,ll)*c0(ll,mm)*betatf(k1,mm)
                    end do
                 end do
                 tmp1=tmp1+tmp2
              end do
              i1=i1+j1
           end do
           tmp1=0.5d0*tmp1
           alpha=rgamma(a0b0(1)+dble(ptf)*dble(ntlr-1),a0b0(2)+tmp1)

c           call dblepr("alpha",-1,alpha,1)     
        end if

c        return

c++++++++++++++++++++++++++++++++         
c++++++ save samples
c++++++++++++++++++++++++++++++++         

        if(iscan.gt.nburn)then
           skipcount=skipcount+1
           if(skipcount.gt.nskip)then
              isave=isave+1
              dispcount=dispcount+1

c++++++++++++ fixed effects

              do i=1,p
                 thetasave(isave,i)=beta(i)
              end do

c++++++++++++ centering variance

              thetasave(isave,p+1)=sigma2b

c++++++++++++ precision parameter

              thetasave(isave,p+2)=alpha

c++++++++++++ random effects

              do i=1,nsubject
                  randsave(isave,i)=b(i)
              end do

c++++++++++++ tailfree regression coefficients

              mm=0
              do i=2,ntlr
                 do j=1,ptf
                    mm=mm+1 
                    tfpsave(isave,mm)=betatf(i,j)
                 end do
              end do

c++++++++++++ density

              call densldtfp(maxm,ngrid,ntlr,ntprob,npred,
     &                       p,ptf,beta,
     &                       betatf,grid,xpred,xtfpred,
     &                       sigma2b,densl,densm,k)

              do i=1,npred
                 write(1) (densl(i,j),j=1,ngrid)
              end do   

c++++++++++++ cpo

              do i=1,nrec

c+++++++++++++++ check if the user has requested an interrupt
                 call rchkusr()

                 yij=y(i)
                 eta=0.d0
                 do j=1,p
                    eta=eta+x(i,j)*beta(j)
                 end do
                 eta=eta+b(subject(i)) 
                 eta=eta+x(i,p+1)
                 mean=exp(eta)

                 tmp2=dpoiss(dble(yij),mean,1)

                 cpo(i,1)=cpo(i,1)+1.d0/exp(tmp2)
                 cpo(i,2)=cpo(i,2)+exp(tmp2)  
              end do

c++++++++++++ quantiles

              call quantileldtfp(maxm,ntlr,npred,p,ptf,
     &                           beta,betatf,
     &                           sigma2b,xtfpred,xpred,
     &                           qmm,qml,k,prob,probc,quans)

              write(3) (qml(j,1),j=1,npred)
              write(4) (qml(j,2),j=1,npred)
              write(5) (qml(j,3),j=1,npred)

c++++++++++++ print
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

      as(1)=alpha
      as(2)=sigma2b
      as(3)=acrate

      do i=1,npred
         do j=1,ngrid
            densm(i,j)=densm(i,j)/dble(nsave)
         end do
      end do

      do i=1,npred
         do j=1,3
            qmm(i,j)=qmm(i,j)/dble(nsave)
         end do
      end do

      do i=1,nrec
         cpo(i,1)=dble(nsave)/cpo(i,1)
         cpo(i,2)=cpo(i,2)/dble(nsave)                                    
      end do

      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=5)

      if(cband.eq.1)then

      call hpddensregdtfp(nsave,npred,ngrid,worksam,
     &                    worksam2,fs,densl,densu,tband)

      call hpdldtfpq1(nsave,npred,worksam,worksam3,
     &                qml,qmu,tband)

      call hpdldtfpq2(nsave,npred,worksam,worksam3,
     &                qml,qmu,tband)

      call hpdldtfpq3(nsave,npred,worksam,worksam3,
     &                qml,qmu,tband)

      end if

      return
      end

