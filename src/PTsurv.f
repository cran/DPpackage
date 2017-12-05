c=======================================================================                      
      subroutine ptsurv
     &                (nrec,nsubject,nvar,p,maxni,datastr,              
     &                 type,left,right,                                 
     &                 x,subject,options,                               
     &                 priorv,prec,sb,tinvb,tauv,                       
     &                 mcmc,nsave,                                      
     &                 acrate,randsave,thetasave,                       
     &                 cpar,beta,b,sigmab,sigmav,y,                     
     &                 iflagp,iflagr,parti,poss,typeint,whicho,whichn,  
     &                 betac,linf,lsup,linf2,lsup2,                     
     &                 sigmainvb,sigmainvv,sigmainvvc,                  
     &                 sigmabc,sigmavc,ssb,ssbc,theta,thetac,vz,vzc,    
     &                 workmr,workmr2,workmhp,                          
     &                 workmhr,workmhr2,workvp,workvr,                  
     &                 xtx,xty,yvec,yvecc,ztz,zty,                      
     &                 seed)                                            
     
c=======================================================================                      
c     # 65 arguments
c
c     Subroutine `ptsurv' to run a Markov chain in a semiparametric 
c     AFT effect model, using a Mixture of Multivariate Polya 
c     trees prior for the distribution of the errors and the frailty 
c     terms.
c
c     Copyright: Alejandro Jara and Tim Hanson, 2009-2010.
c
c     Version 1.0:
c
c     Last modification: 20-06-2007.
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
c=======================================================================                      
      implicit none

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Data
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer maxni,nrec,nsubject,nvar,p
      integer datastr(nsubject,maxni+1)
      integer type(nrec,nvar)
      integer subject(nrec)
      double precision left(nrec,nvar),right(nrec,nvar)
      double precision x(nrec,p)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Options
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer frailty,latent,options(2)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Prior 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer m1,m2
      double precision prec(p*nvar,p*nvar),sb(p*nvar,2)
      double precision nu0b,tinvb(nvar,nvar)
      double precision tauv(nvar,2)
      double precision priorv(3)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++MCMC
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer mcmc(5),nburn,nskip,nsave,ndisplay
      double precision tune1,tune2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Output
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision acrate(5)
      double precision randsave(nsave,nsubject*nvar+2*nvar)
      double precision thetasave(nsave,p*nvar+nvar+(nvar*(nvar+1))/2+2)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Current values of the parameters
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      double precision cpar(2),cparb,cparv,beta(p*nvar),b(nsubject,nvar)
      double precision sigmab(nvar,nvar),sigmav(nvar,nvar)
      double precision y(nrec,nvar)
      
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer iflagp(p*nvar)
      integer iflagr(nvar)
      integer parti(nvar)
      integer poss(nvar,p)
      integer typeint(nvar)
      integer whicho(nrec),whichn(nrec)      

      double precision betac(p*nvar)
      double precision linf(nvar),lsup(nvar)
      double precision linf2(nvar),lsup2(nvar)
      double precision sigmainvb(nvar,nvar)
      double precision sigmainvv(nvar,nvar)
      double precision sigmainvvc(nvar,nvar)
      double precision sigmabc(nvar,nvar)      
      double precision sigmavc(nvar,nvar)
      double precision ssb(nsubject)
      double precision ssbc(nsubject)      
      double precision theta(nvar)
      double precision thetac(nvar)
      double precision vz(nrec,nvar),vzc(nrec,nvar)
      double precision workmr(nvar,nvar),workmr2(nvar,nvar)
      double precision workmhp(p*nvar*(p*nvar+1)/2)
      double precision workmhr(nvar*(nvar+1)/2),
     1  workmhr2(nvar*(nvar+1)/2)
      double precision workvp(p*nvar)
      double precision workvr(nvar)
      double precision xtx(p*nvar,p*nvar)
      double precision xty(p*nvar)
      double precision yvec(nvar),yvecc(nvar)
      double precision ztz(nvar,nvar)
      double precision zty(nvar)

c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Working space - RNG
c+++++++++++++++++++++++++++++++++++++++++++++++++
      integer seed1,seed2,seed(2)

c+++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++
      integer count
      integer dispcount
      integer i,j,k,l
      integer i1,j1
      integer ii,jj,kk,ll,mm
      integer iscan,isave
      integer marea
      parameter (marea=2**20)
      integer massi(marea)
      integer narea
      integer ni
      integer nscan
      integer nu
      integer skipcount
      integer sprint

      double precision acrate2      
      double precision detlogb,detlogbc,detlogv,detlogvc
      double precision logcgkn,logcgko
      double precision loglikn,loglikn2,logliko,logliko2
      double precision logpriorn,logprioro
      double precision mass(marea)
      double precision ratio
      double precision sd,sdc
      double precision ssevaln,ssevalo
      double precision tmp1,tmp2
      
      integer maxp
      parameter(maxp=200)
      double precision mumh(maxp),sigmamh(maxp,maxp)
      
      double precision dlnrm,rtlnorm
      real runif

c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++CPU time
c+++++++++++++++++++++++++++++++++++++++++++++++++
      double precision sec00,sec0,sec1,sec

c+++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Initialize variables
c+++++++++++++++++++++++++++++++++++++++++++++++++

c+++++checking dimensions
      narea=2**nvar  
      if(narea.gt.marea)then
         call rexit("increase dimension marea in subroutine ptsurv")
      end if

      if(nvar*p.gt.maxp)then
         call rexit("increase dimension maxp in subroutine ptsurv")
      end if
      
      do i=1,marea
         massi(i)=0
         mass(i)=0.d0
      end do
      
      do i=1,maxp
         do j=1,maxp
            sigmamh(i,j)=0.d0
         end do
         mumh(i)=0.d0
      end do

c++++ options
      frailty=options(1)
      latent=options(2)

c++++ mcmc, priors and "zipped"

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      m1=mcmc(4)
      m2=mcmc(5)

      nu0b=priorv(1)
      tune1=priorv(2)
      tune2=priorv(3)
         
      seed1=seed(1)
      seed2=seed(2)
      
      cparv=cpar(1)
      cparb=cpar(2)
      
c++++ set random number generator

      call setall(seed1,seed2)
      
c++++ set possition of covariates
     
      count=0 
      do i=1,nvar
         do j=1,p
            count=count+1
            poss(i,j)=count    
         end do
      end do   

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*nsave

      do i=1,nvar*p
         mumh(i)=beta(i)
         do j=1,nvar*p
            sigmamh(i,j)=0.d0
         end do
         sigmamh(i,i)=1.d0
      end do
      
      call cpu_time(sec0)
      sec00=0.d0

c++++ first evaluation of the likelihood

      detlogv=0.d0
      do i=1,nvar
         sigmainvv(i,i)=1.d0/sigmav(i,i)
         detlogv=detlogv+log(sigmav(i,i))
      end do

      call loglikpt_betacan_mz(m1,nvar,nrec,nsubject,subject,p,parti,
     &                         whicho,whichn,y,vz,cparv,detlogv,
     &                         linf,lsup,theta,sigmainvv,
     &                         x,beta,poss,b,logliko)


c      call dblepr("logliko",-1,logliko,1) 
c      call dblepr("beta",-1,beta,p*nvar) 
c      call dblepr("sigmav",-1,sigmav,nvar*nvar) 
c      call intpr("frailty",-1,frailty,1) 
c      call intpr("m1",-1,m1,1) 
c      call intpr("m2",-1,m2,1) 

      if(frailty.eq.1)then
         do i=1,nvar
            do j=1,nvar
               sigmainvb(i,j)=sigmab(i,j)
            end do
         end do
         call inversedet(sigmainvb,nvar,iflagr,detlogb)
      end if


      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c         call intpr("Step1",-1,latent,1)

c++++++++++++++++++++++++++++++++
c+++++++ latent variables     +++
c++++++++++++++++++++++++++++++++
         if(latent.eq.1)then
          
         acrate2=0.d0
         
         do i=1,nrec
         
c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do j=1,nvar
               yvec(j)=y(i,j)
               linf(j)=left(i,j)
               lsup(j)=right(i,j)
               linf2(j)=left(i,j)
               lsup2(j)=right(i,j)
               typeint(j)=type(i,j)
            end do

            do j=1,nvar
               tmp1=0.d0
               do k=1,p
                  l=poss(j,k) 
                  tmp1=tmp1-x(i,k)*beta(l)
               end do
               theta(j)=tmp1+b(subject(i),j)
            end do

c++++++++++ generating the candidate value

            call rtmvn(nvar,theta,sigmav,linf,lsup,workmr,workmr2,
     &                 typeint,yvecc,logcgko)
     
            call dtmvn(nvar,theta,sigmav,linf,lsup,workmr,workmr2,
     &                 typeint,yvecc,yvec,logcgkn)

c++++++++++ evaluating likelihood

            call loglikpt_ycan_mz(i,m1,nvar,nrec,nsubject,subject,
     &                            p,parti,whicho,whichn,yvec,zty,
     &                            vz,cparv,detlogv,linf,lsup,theta,
     &                            sigmainvv,x,beta,poss,b,logliko)   

            call loglikpt_ycan_mz(i,m1,nvar,nrec,nsubject,subject,
     &                            p,parti,whicho,whichn,yvecc,zty,
     &                            vz,cparv,detlogv,linf,lsup,theta,
     &                            sigmainvv,x,beta,poss,b,loglikn)   

c++++++++++ mh step
  
            ratio=loglikn-logliko+logcgkn-logcgko

            if(log(dble(runif())).lt.ratio)then
               acrate2=acrate2+1.d0
               do j=1,nvar
                  y(i,j)=yvecc(j)
                  vz(i,j)=zty(j)
               end do
            end if
         end do
         acrate(1)=acrate(1)+acrate2/dble(nrec)
         
c+++++++ evaluating the current value of the likelihood

         call loglikpt_betacan_mz(m1,nvar,nrec,nsubject,subject,p,parti,
     &                         whicho,whichn,y,vz,cparv,detlogv,
     &                         linf,lsup,theta,sigmainvv,
     &                         x,beta,poss,b,logliko)

         end if

c         call dblepr("y",-1,y,nsubject*nvar)

c++++++++++++++++++++++++++++++++
c+++++++ fixed effects        +++
c++++++++++++++++++++++++++++++++

         do i=1,p*nvar
            do j=1,p*nvar
               xtx(i,j)=prec(i,j)
            end do
            xty(i)=sb(i,1)            
         end do   

         do ii=1,nrec
            do jj=1,nvar
               yvec(jj)=y(ii,jj)-b(subject(ii),jj)
            end do
         
            do jj=1,nvar
               do kk=1,p
                  i1=poss(jj,kk)
                  do ll=1,nvar
                     do mm=1,p
                        j1=poss(ll,mm)
                        xtx(i1,j1)=xtx(i1,j1)+sigmainvv(jj,ll)*
     &                             x(ii,kk)*x(ii,mm)
                     end do
                     xty(i1)=xty(i1)-sigmainvv(jj,ll)*yvec(ll)*x(ii,kk)
                  end do
               end do
            end do
         end do
         call inverse(xtx,p*nvar,iflagp)      

         if(iscan.lt.nburn/2)then
            do i=1,p*nvar
               do j=1,p*nvar 
                  xtx(i,j)=0.01d0*xtx(i,j)
               end do   
            end do   
            call rmvnorm(p*nvar,beta,xtx,workmhp,workvp,betac)

          else
            ratio=dble(runif())
            if(ratio.le.0.25d0)then
               do i=1,p*nvar
                  do j=1,p*nvar
                      xtx(i,j)=(5.4264d0/dble(p*nvar))*sigmamh(i,j)
                  end do
               end do
               call rmvnorm(p*nvar,beta,xtx,workmhp,workvp,betac)
             else if(ratio.le.0.5d0)then
               do i=1,p*nvar
                  do j=1,p*nvar
                      xtx(i,j)=(0.001d0)*xtx(i,j)
                  end do
               end do
               call rmvnorm(p*nvar,beta,xtx,workmhp,workvp,betac)
             else if(ratio.le.0.75d0)then
               do i=1,p*nvar
                  do j=1,p*nvar
                      xtx(i,j)=(0.001d0)*xtx(i,j)
                  end do
               end do
               call rmvnorm(p*nvar,beta,xtx,workmhp,workvp,betac)
             else 
               do i=1,p*nvar
                  do j=1,p*nvar 
                     xtx(i,j)=xtx(i,j)
                  end do   
               end do   
               call rmvnorm(p*nvar,beta,xtx,workmhp,workvp,betac)
            end if
         end if


c+++++++ prior ratio
         logprioro=0.d0
         logpriorn=0.d0
         
         do i=1,p*nvar
            do j=1,p*nvar
               logpriorn=logpriorn+(betac(i)-sb(i,2))* 
     &                   prec(i,j)       *
     &                  (betac(j)-sb(j,2))

               logprioro=logprioro+(beta(i)-sb(i,2))* 
     &                   prec(i,j)      *
     &                   (beta(j)-sb(j,2))

            end do
         end do
      
         logpriorn=-0.5d0*logpriorn
         logprioro=-0.5d0*logprioro

c+++++++ evaluating likelihood for betac

         call loglikpt_betacan_mz(m1,nvar,nrec,nsubject,subject,p,parti,
     &                         whicho,whichn,y,vzc,cparv,detlogv,
     &                         linf,lsup,theta,sigmainvv,
     &                         x,betac,poss,b,loglikn)

c         call dblepr("logliko",-1,logliko,1) 
c         call dblepr("loglikn",-1,loglikn,1) 

c+++++++ acceptance step

         ratio=loglikn-logliko+logpriorn-logprioro

         if(log(dble(runif())).lt.ratio)then
            acrate(2)=acrate(2)+1.d0
            do i=1,p*nvar
               beta(i)=betac(i) 
            end do
            logliko=loglikn
            do i=1,nrec
               do j=1,nvar
                  vz(i,j)=vzc(i,j)
               end do
            end do
         end if

c+++++++ adapting the parameters for the MH algorithm

         do i=1,p*nvar
            do j=1,p*nvar
               tmp1=(beta(i)-mumh(i))*(beta(j)-mumh(j)) 
               
               if(tmp1.eq.0.d0)then
                  tmp2=sigmamh(i,j)/dble(iscan)
                else
                  tmp2=
     &            (tmp1- 
     &             (dble(iscan+1)/dble(iscan))*sigmamh(i,j)
     &            )/dble(iscan+1)
               end if
               sigmamh(i,j)=sigmamh(i,j)+tmp2
            end do
         end do
            
         do i=1,p*nvar
            tmp1=(beta(i)-mumh(i))/dble(iscan+1)
            mumh(i)=mumh(i)+tmp1
         end do

         call dblepr("beta",-1,beta,p*nvar)

c         call dblepr("betac",-1,betac,p*nvar)
c         call dblepr("mumh",-1,mumh,p*nvar)
c         do i=1,p*nvar
c            do j=i,p*nvar
c                call dblepr("sigmamh",-1,sigmamh(i,j),1)
c            end do
c         end do   

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating sigmav using a MH step                    +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ generating the candidate value/
c+++++++ evaluating the candidate generating kernel/
c+++++++ evaluate log-prior

         logcgkn=0.d0
         logcgko=0.d0
         detlogvc=0.d0
         logpriorn=0.d0
         logprioro=0.d0

         do i=1,nvar
            sd=sigmainvv(i,i)
            sdc=rtlnorm(log(sd),tune1*0.01d0,0,0,.true.,.true.)
            sigmainvvc(i,i)=sdc

            sigmavc(i,i)=1.d0/sigmainvvc(i,i)
            detlogvc=detlogvc+log(sigmavc(i,i))

            logcgkn=logcgkn+dlnrm(sd,log(sdc),tune1*0.01d0,1) 
            logcgko=logcgko+dlnrm(sdc,log(sd ),tune1*0.01d0,1) 

            call dgamma2(sigmainvvc(i,i),0.5d0*tauv(i,1),
     &                   0.5d0*tauv(i,2),tmp1)  
            logpriorn=logpriorn+tmp1 

            call dgamma2(sigmainvv(i,i) ,0.5d0*tauv(i,1),
     &                   0.5d0*tauv(i,2),tmp1)  
            logprioro=logprioro+tmp1 
         end do

c+++++++ evaluating likelihood

         call loglikpt_betacan_mz(m1,nvar,nrec,nsubject,subject,p,parti,
     &                            whicho,whichn,y,vzc,cparv,detlogvc,
     &                            linf,lsup,yvec,sigmainvvc,
     &                            x,beta,poss,b,loglikn)

c+++++++ acceptance step
         
         ratio=loglikn-logliko+logcgkn-logcgko+
     &         logpriorn-logprioro

         if(log(dble(runif())).lt.ratio)then
            do i=1,nvar
               sigmav(i,i)=sigmavc(i,i)
               sigmainvv(i,i)=sigmainvvc(i,i)
            end do
            detlogv=detlogvc
            do i=1,nrec
               do j=1,nvar
                  vz(i,j)=vzc(i,j) 
               end do
            end do
            logliko=loglikn
            acrate(3)=acrate(3)+1.d0
         end if

         call dblepr("sigmav",-1,sigmav,nvar*nvar)
         
        
c++++++++++++++++++++++++++++++++
c+++++++ frailty terms        +++
c++++++++++++++++++++++++++++++++

         if(frailty.eq.1)then
            acrate2=0.d0
            do i=1,nsubject

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               do j=1,nvar
                  theta(j)=b(i,j)
               end do   

c+++++++++++++ generating a candidate
               do j=1,nvar
                  tmp1=0.d0
                  do k=1,nvar
                     ztz(j,k)=sigmainvb(j,k)
                  end do
                  zty(j)=0.d0
               end do

               ni=datastr(i,1) 

               do j=1,ni
                  do k=1,nvar
                     tmp1=0.d0
                     do l=1,p
                        jj=poss(j,k) 
                        tmp1=tmp1-x(datastr(i,j+1),l)*beta(jj)
                     end do
                     thetac(k)=tmp1
                  end do
   
                  do k=1,nvar
                     ztz(k,k)=ztz(k,k)+sigmainvv(k,k)
                     zty(k)=zty(k)+(y(datastr(i,j+1),k)-thetac(k))*
     &                      sigmainvv(k,k)              
                  end do
               end do   

               call inverse(ztz,nvar,iflagr)      

               call rmvnorm(nvar,theta,ztz,workmhr,workvr,thetac)

c+++++++++++++ evaluating the prior

               logprioro=0.d0
               logpriorn=0.d0

               call condptpriorell(i,m2,nrec,nsubject,nvar,theta,ssb,
     &                          sigmainvb,detlogb,cparb,
     &                          whicho,whichn,
     &                          logprioro,ssevalo)

               call condptpriorell(i,m2,nrec,nsubject,nvar,thetac,ssb,
     &                          sigmainvb,detlogb,cparb,
     &                          whicho,whichn,
     &                          logpriorn,ssevaln)

c+++++++++++++ evaluating likelihood 

               do j=1,nvar
                  b(i,j)=theta(j)
               end do

               call loglikpt_bcan_mz(m1,nvar,nrec,nsubject,subject,p,
     &                         parti,whicho,whichn,y,vzc,cparv,detlogv,
     &                         linf,lsup,yvec,sigmainvv,
     &                         x,beta,poss,b,logliko)

               do j=1,nvar
                  b(i,j)=thetac(j)
               end do

               call loglikpt_bcan_mz(m1,nvar,nrec,nsubject,subject,p,
     &                         parti,whicho,whichn,y,vzc,cparv,detlogv,
     &                         linf,lsup,yvec,sigmainvv,
     &                         x,beta,poss,b,loglikn)

c+++++++++++++ mh step
  
               ratio=loglikn-logliko+
     &               logcgkn-logcgko+
     &               logpriorn-logprioro

               if(log(dble(runif())).lt.ratio)then
                  acrate2=acrate2+1.d0
                  do j=1,nvar
                     b(i,j)=thetac(j)
                     vz(i,j)=vzc(i,j)
                  end do
                  ssb(i)=ssevaln
                  logliko=loglikn
                 else
                  do j=1,nvar
                     b(i,j)=theta(j)
                  end do
               end if
            end do
            acrate(4)=acrate(4)+acrate2/dble(nsubject)
         end if
         

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Updating sigmab using a MH step                    +++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(frailty.eq.1)then

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

c++++++++++ generating the candidate value

            nu=int((dble(nsubject))*tune2)
         
            do i=1,nvar
               do j=1,nvar
                  sigmabc(i,j)=dble(nu-nvar-1)*sigmab(i,j)
               end do
            end do

            call riwishart(nvar,nu,sigmabc,sigmainvvc,workmr,workvr,
     &                     workmhr,workmhr2,iflagr)


c++++++++++ evaluating the candidate generating kernel

            do i=1,nvar
               do j=1,nvar
                  ztz(i,j)=dble(nu-nvar-1)*sigmab(i,j)
               end do
            end do

            call diwishart(nvar,nu,sigmabc,ztz,workmr,workmr2,workvr,
     &                  iflagr,logcgko)        

            do i=1,nvar
               do j=1,nvar
                  ztz(i,j)=dble(nu-nvar-1)*sigmabc(i,j)
               end do
            end do

            call diwishart(nvar,nu,sigmab,ztz,workmr,workmr2,workvr,
     &                  iflagr,logcgkn)        


c++++++++++ evaluating the prior

            call diwishart(nvar,int(nu0b),sigmabc,tinvb,workmr,
     &                  workmr2,workvr,iflagr,logpriorn)        
 
            call diwishart(nvar,int(nu0b),sigmab,tinvb,workmr,
     &                  workmr2,workvr,iflagr,logprioro)        

c++++++++++ evaluating likelihood

            call loglik_ell(m2,nrec,nsubject,nvar,b,ssb,
     &                      sigmainvb,detlogb,cparb,
     &                      whicho,whichn,
     &                      logliko2)

            do i=1,nvar
               do j=1,nvar
                  sigmainvvc(i,j)=sigmabc(i,j)
               end do
            end do

            call inversedet(sigmainvvc,nvar,iflagr,detlogbc) 

            call loglik_ell(m2,nrec,nsubject,nvar,b,ssbc,
     &                      sigmainvvc,detlogbc,cparb,
     &                      whicho,whichn,
     &                      loglikn2)

c++++++++++ acceptance step
         
            ratio=loglikn2-logliko2+logcgkn-logcgko+
     &            logpriorn-logprioro

            if(log(dble(runif())).lt.ratio)then
               do i=1,nvar
                  do j=1,nvar
                     sigmab(i,j)=sigmabc(i,j)
                     sigmainvb(i,j)=sigmainvvc(i,j)
                  end do
               end do
               detlogb=detlogbc
               do i=1,nsubject
                  ssb(i)=ssbc(i)
               end do
               acrate(5)=acrate(5)+1.d0
            end if

         end if
         
c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         

         cpar(1)=cparv
         cpar(2)=cparb

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ regression coefficient information

               do i=1,p*nvar
                  thetasave(isave,i)=beta(i)
               end do

c+++++++++++++ baseline covariance for the errors

               do i=1,nvar
                  thetasave(isave,p*nvar+i)=sigmav(i,i)
               end do

c+++++++++++++ baseline covariance for the frailty

               if(frailty.eq.1)then
                  k=0
                  do i=1,nvar
                     do j=i,nvar
                        k=k+1
                        thetasave(isave,p*nvar+nvar+k)=sigmab(i,j)
                     end do
                  end do
                else
                  k=nvar*(nvar+1)/2
               end if   

c+++++++++++++ precision parameter information

               thetasave(isave,p*nvar+nvar+k+1)=cparv
               if(frailty.eq.1)then
                  thetasave(isave,p*nvar+nvar+k+2)=cparb
               end if   

c+++++++++++++ frailty

               if(frailty.eq.1)then
                  k=0
                  do i=1,nsubject
                     do j=1,nvar
                        k=k+1
                        randsave(isave,k)=b(i,j)
                     end do   
                  end do
                 else
                  k=nsubject*nvar
               end if   

c+++++++++++++ error predictive information

               call sampredptun(marea,nvar,nrec,parti,m1,mass,massi,
     &                          iflagr,typeint,whichn,whicho,vz,
     &                          cparv,thetac,linf,lsup,workmr,
     &                          workmhr,workvr,theta,1)   

               do i=1,nvar
                  tmp1=0.d0
                  theta(i)=sqrt(sigmav(i,i))*theta(i)
                  k=k+1
                  randsave(isave,k)=theta(i)
               end do


c+++++++++++++ frailty predictive information
               if(frailty.eq.1)then
                  do i=1,nvar
                     thetac(i)=0.d0
                  end do

                  call sampredellpt(nrec,nsubject,nvar,m2,ssb,cparb,
     &                              whicho,whichn,workmr,workmhr,workvr,
     &                              thetac,sigmab,theta)   
                  do i=1,nvar
                     k=k+1
                     randsave(isave,k)=theta(i)
                  end do
               end if   

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

      do i=1,5
         acrate(i)=acrate(i)/dble(nscan)
      end do

      
      return

      end


c=======================================================================                  
      subroutine loglikpt_ycan_mz(ind,m,nvar,nrec,nsubject,subject,
     &                            p,parti,
     &                            whicho,whichn,yc,tvec,vz,cpar,detlogl,
     &                            linf,lsup,theta,sigmainv,
     &                            x,beta,poss,b,loglikc)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     value of the survival time in a marginal Multivariate PT.
c
c     Alejandro Jara, 2006
c======================================================================= 
      implicit none 

c-----Input
      integer ind,m,nvar,nrec,nsubject,p
      integer parti(nvar),poss(nvar,p)
      integer subject(nrec)
      integer whicho(nrec),whichn(nrec)
      double precision beta(nvar*p),b(nsubject,nvar),x(nrec,p)
      double precision yc(nvar),tvec(nvar),vz(nrec,nvar),cpar,detlogl
      double precision linf(nvar),lsup(nvar)
      double precision theta(nvar),sigmainv(nvar,nvar)

c-----Output
      double precision loglikc

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan,tmp1

c-----Routine

      loglikc=0.d0

      do i=1,nvar
         tmp1=0.d0
         do j=1,p
            k=poss(i,j) 
            tmp1=tmp1-x(ind,j)*beta(k)
         end do
         theta(i)=tmp1+b(subject(ind),i)
      end do

      do i=1,nvar
         tmp1=0.d0
         do j=1,nvar
            tmp1=tmp1+sigmainv(i,j)*(yc(j)-theta(j))
         end do
         tvec(i)=tmp1
      end do

      nint=2
      prob=1.d0/dble(nint)
      quan=invcdfnorm(prob,0.d0,1.d0,1,0)

      countero=0
      
      do j=1,nvar
         if(tvec(j).le.quan)then
            linf(j)=-999999.d0
            lsup(j)=quan
            parti(j)=1
          else
            linf(j)=quan
            lsup(j)= 999999.d0
            parti(j)=2
         end if
      end do
      
      do l=1,nrec
         if(l.ne.ind)then 
            final=1
            do j=1,nvar
               if(vz(l,j).gt.lsup(j).or.vz(l,j).lt.linf(j))then
                  final=0
               end if
c               if(y(l,j).gt.lsup2(j).or.y(l,j).lt.linf2(j))then
c                  final=0
c               end if
            end do

            if(final.eq.1)then
               countero=countero+1
               whicho(countero)=l
            end if   
         end if
      end do

      if(countero.eq.0) go to 1

      ok=1
      j=2
      do while(ok.eq.1.and.j.le.m)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)

         do k=1,nvar
            k1=2*(parti(k)-1)+1
            k2=2*(parti(k)-1)+2
            quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
            if(tvec(k).le.quan)then
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
               if(vz(whicho(l),k).gt.lsup(k).or.
     &            vz(whicho(l),k).lt.linf(k)    )then
                  final=0 
               end if   
            end do
                  
            if(final.eq.1)then
              countern=countern+1
              whichn(countern)=whicho(l)
            end if
         end do

         loglikc=loglikc+
     &           log((2.d0**nvar)*cpar*dble(je2)+
     &               dble((2.d0**nvar)*countero))-
     &           log((2.d0**nvar)*cpar*dble(je2)+dble(i-1))

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

1     continue

      loglikc=loglikc-0.5d0*detlogl
      do j=1,nvar
         loglikc=loglikc+dnrm(tvec(j),0.d0, 1.d0, 1)
      end do   

      return
      end



c=======================================================================                  
      subroutine loglikpt_betacan_mz(m,nvar,nrec,nsubject,subject,
     &                               p,parti,whicho,whichn,y,vzc,cpar,
     &                               detlogl,linf,lsup,theta,sigmainv,
     &                               x,betac,poss,b,loglikc)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     value of the regression coefficients using a marginal 
c     Multivariate PT in a median zero AFT regression model.
c
c     Alejandro Jara, 2006
c======================================================================= 
      implicit none 

c-----Input
      integer m,nvar,nrec,nsubject,p
      integer parti(nvar),poss(nvar,p)
      integer subject(nrec)
      integer whicho(nrec),whichn(nrec)
      double precision betac(nvar*p),b(nsubject,nvar),x(nrec,p)
      double precision y(nrec,nvar),vzc(nrec,nvar),cpar,detlogl
      double precision linf(nvar),lsup(nvar)
      double precision theta(nvar),sigmainv(nvar,nvar)

c-----Output
      double precision loglikc

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan,tmp1

c-----Routine

      loglikc=0.d0

      do i=1,nrec
         do j=1,nvar
            tmp1=0.d0
            do k=1,p
               l=poss(j,k) 
               tmp1=tmp1-x(i,k)*betac(l)
            end do
            theta(j)=tmp1+b(subject(i),j)
         end do
      
         do j=1,nvar
            tmp1=0.d0
            do k=1,nvar
               tmp1=tmp1+sigmainv(j,k)*(y(i,k)-theta(k))
            end do
            vzc(i,j)=tmp1
         end do
         
c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first subject
         if(i.eq.1)then
            loglikc=-0.5d0*detlogl
            do j=1,nvar
               loglikc=loglikc+dnrm(vzc(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following subjects
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nvar
               if(vzc(i,j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nvar
                  if(vzc(l,j).gt.lsup(j).or.vzc(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do

            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nvar
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(vzc(i,k).le.quan)then
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
                     if(vzc(whicho(l),k).gt.lsup(k).or.
     &                  vzc(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               loglikc=loglikc+
     &            log((2.d0**nvar)*cpar*dble(je2)+
     &                 dble((2.d0**nvar)*countero))-
     &            log((2.d0**nvar)*cpar*dble(je2)+dble(i-1))

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

            loglikc=loglikc-0.5d0*detlogl
            do j=1,nvar
               loglikc=loglikc+dnrm(vzc(i,j),0.d0, 1.d0, 1)
            end do   
         end if
      end do   

      return
      end


c=======================================================================                  
      subroutine loglikpt_bcan_mz(m,nvar,nrec,nsubject,subject,p,parti,
     &                            whicho,whichn,y,vzc,cpar,detlogl,
     &                            linf,lsup,muc,sigmainv,
     &                            x,beta,poss,b,loglikc)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     value of the frailty term in a marginal Multivariate PT.
c
c     Alejandro Jara, 2006
c======================================================================= 
      implicit none 

c-----Input
      integer m,nvar,nrec,nsubject,p
      integer parti(nvar),poss(nvar,p)
      integer subject(nrec)
      integer whicho(nrec),whichn(nrec)
      double precision beta(nvar*p),b(nsubject,nvar),x(nrec,p)
      double precision y(nrec,nvar),vzc(nrec,nvar),cpar,detlogl
      double precision linf(nvar),lsup(nvar)
      double precision muc(nvar),sigmainv(nvar,nvar)

c-----Output
      double precision loglikc

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan,tmp1

c-----Routine

      loglikc=0.d0

      do i=1,nrec
         
         do j=1,nvar
            tmp1=0.d0
            do k=1,p
               l=poss(j,k) 
               tmp1=tmp1-x(i,k)*beta(l)
            end do
            muc(j)=tmp1+b(subject(i),j)
         end do
      
         do j=1,nvar
            tmp1=0.d0
            do k=1,nvar
               tmp1=tmp1+sigmainv(j,k)*(y(i,k)-muc(k))
            end do
            vzc(i,j)=tmp1
         end do
         
c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first subject
         if(i.eq.1)then
            loglikc=-0.5d0*detlogl
            do j=1,nvar
               loglikc=loglikc+dnrm(vzc(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following subjects
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nvar
               if(vzc(i,j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nvar
                  if(vzc(l,j).gt.lsup(j).or.vzc(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do

            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nvar
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(vzc(i,k).le.quan)then
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
                     if(vzc(whicho(l),k).gt.lsup(k).or.
     &                  vzc(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               loglikc=loglikc+
     &            log((2.d0**nvar)*cpar*dble(je2)+
     &                dble((2.d0**nvar)*countero))-
     &            log((2.d0**nvar)*cpar*dble(je2)+dble(i-1))

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

            loglikc=loglikc-0.5d0*detlogl
            do j=1,nvar
               loglikc=loglikc+dnrm(vzc(i,j),0.d0, 1.d0, 1)
            end do   
         end if
      end do   

      return
      end


