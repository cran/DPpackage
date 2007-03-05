c=======================================================================                      
      subroutine DDPrasch(ngrid,npred,nsubject,p,q,y,grid,z,zpred,      #9
     &                    murand,a0b0,b0,prec1,sb,mu0,prec2,smu,tinv,   #9
     &                    acrate,cpo,densr,randsave,thetasave,          #5
     &                    ncluster,ss,alpha,m,mu,beta,b,sigma,          #8
     &                    sigmainv,sigma2b,                             #2           
     &                    mcmc,nsave,                                   #2
     &                    iflagp,betac,xtx,xty,workmhp1,                #5 
     &                    workmp1,workmp2,workmp3,workvp1,workvp2,      #5 
     &                    ccluster,iflagq,alphac,prob,quadf,            #5
     &                    workmhq1,workmhq2,                            #2
     &                    workmq1,workmq2,workmq3,workvq1,workvq2,      #5
     &                    wtw,wty,                                      #2
     &                    seed)                                         #1 
     
c=======================================================================                  
c 60 elements
c=======================================================================

      implicit none 

c+++++Data
      integer ngrid,npred,nsubject,p,q
      integer y(nsubject,p)
      real*8 grid(ngrid),z(nsubject,q),zpred(npred,q)

c+++++Prior 
      integer murand 
      real*8 aa0,ab0,a0b0(5)
      real*8 b0(p-1),prec1(p-1,p-1),sb(p-1)
      real*8 mu0(q),prec2(q,q),smu(q)
      real*8 tau1,tau2
      real*8 nu0,tinv(q,q) 

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      real*8 acrate(2),cpo(nsubject,p),densr(npred,ngrid)
      real*8 randsave(nsave,nsubject+npred)
      real*8 thetasave(nsave,p+q+(q*(q+1)/2)+2)

c+++++Current values of the parameters
      integer ncluster,ss(nsubject) 
      real*8 alpha(nsubject,q),m,mu(q),beta(p-1),b(nsubject)
      real*8 sigma(q,q),sigmainv(q,q),sigma2b


c+++++External Working space - Difficulty parameters
      integer iflagp(p-1)
      real*8 betac(p-1)
      real*8 xtx(p-1,p-1),xty(p-1)
      real*8 workmhp1((p-1)*p/2)
      real*8 workmp1(p-1,p-1),workmp2(p-1,p-1),workmp3(p-1,p-1)
      real*8 workvp1(p-1),workvp2(p-1)

c+++++External Working space - DDP part
      integer ccluster(nsubject),iflagq(q)
      real*8 alphac(q)
      real*8 prob(nsubject+2)
      real*8 quadf(q,q)
      real*8 workmhq1(q*(q+1)/2),workmhq2(q*(q+1)/2)
      real*8 workmq1(q,q),workmq2(q,q),workmq3(q,q)
      real*8 workvq1(q),workvq2(q)
      real*8 wtw(q,q),wty(q)

c+++++External Working space - RNG
      integer seed(2),seed1,seed2

c+++++Internal Working space
      integer counter
      integer dispcount
      integer evali
      integer i,ii
      integer iscan,isave
      integer j,k
      integer ns,nscan
      integer since
      integer skipcount
      integer sprint
      integer yij
      real*8 acrate2
      real*8 dbin,detlog,dnrm
      real*8 eta,gprime
      real*8 logcgkn,logcgko
      real*8 loglikn,logliko
      real*8 logpriorn,logprioro
      real*8 mean
      real*8 offset
      real*8 ratio,rgamma,rnorm,runif
      real*8 ssb
      real*8 thetac
      real*8 tmp1,tmp2
      real*8 winv
      real*8 ytilde
      real*8 ztz,ztzinv,zty

c+++++CPU time
      real*8 sec00,sec0,sec1,sec

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

      aa0=a0b0(1)
      ab0=a0b0(2)
      tau1=a0b0(3)
      tau2=a0b0(4)
      nu0=a0b0(5)

c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      
      call setall(seed1,seed2)

c++++ cluster structure
      do i=1,nsubject
         ccluster(ss(i))=ccluster(ss(i))+1
      end do
     
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)
      
      call cpu_time(sec0)
      sec00=0.d0

      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ difficulty parameters
c++++++++++++++++++++++++++++++++++

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec1(i,j)
               workmp1(i,j)=0.d0
               workmp2(i,j)=0.d0
               workmp3(i,j)=0.d0
            end do
            xty(i)=sb(i)
            workvp1(i)=0.d0
            workvp2(i)=0.d0
            iflagp(i)=0
         end do
         
         logliko=0.d0
         
         do i=1,nsubject
            do j=1,p-1
               eta=b(i)-beta(j) 
               offset=b(i)
               mean=exp( eta - log(1.d0+exp(eta)))
               gprime=exp( -log(mean)-log(1.0-mean))
               winv=exp( -log(mean)-log(1.0-mean))
               
               ytilde=eta+(dble(y(i,j+1))-mean)*gprime-offset

               xtx(j,j)=xtx(j,j)+1.d0/winv
               xty(j)=xty(j)-ytilde/winv
               
               logliko=logliko+dbin(dble(y(i,j+1)),1.d0,mean,1)
               
            end do
         end do

c         call dblepr("xtx",-1,xtx,(p-1)*(p-1))
         
         call invdet(xtx,p-1,workmp1,detlog,iflagp,workvp1)

         do i=1,p-1
            tmp1=0.d0
            workvp1(i)=0.d0
            do j=1,p-1
               tmp1=tmp1+workmp1(i,j)*xty(j) 
            end do
            workvp1(i)=tmp1
         end do

            
         call rmvnorm(p-1,workvp1,workmp1,workmhp1,workvp2,betac)


         call dmvn2(p-1,betac,workvp1,workmp1,logcgko,
     &              workvp2,workmp2,workmp3,iflagp)                 

c+++++++ prior ratio

         tmp1=0.d0
         tmp2=0.d0
         
         do i=1,p-1
            do j=1,p-1
               tmp1=tmp1+(betac(i)-b0(i))* 
     &                    prec1(i,j)      *
     &                   (betac(j)-b0(j))

               tmp2=tmp2+(beta(i) -b0(i))* 
     &                    prec1(i,j)      *
     &                   (beta(j) -b0(j))
            end do
         end do

         logpriorn=-0.5d0*tmp1
         logprioro=-0.5d0*tmp2

            
c+++++++ candidate generating kernel contribution

         do i=1,p-1
            do j=1,p-1
               xtx(i,j)=prec1(i,j)
               workmp1(i,j)=0.d0
               workmp2(i,j)=0.d0
               workmp3(i,j)=0.d0
            end do
            xty(i)=sb(i)
            workvp1(i)=0.d0
            workvp2(i)=0.d0
            iflagp(i)=0
         end do

         loglikn=0.d0

         do i=1,nsubject
            do j=1,p-1
               eta=b(i)-betac(j) 
               offset=b(i)
               mean=exp( eta - log(1.d0+exp(eta)))
               gprime=exp( -log(mean)-log(1.0-mean))
               winv=exp( -log(mean)-log(1.0-mean))
               
               ytilde=eta+(dble(y(i,j+1))-mean)*gprime-offset
               
               loglikn=loglikn+dbin(dble(y(i,j+1)),1.d0,mean,1)
               
               xtx(j,j)=xtx(j,j)+1.d0/winv
               xty(j)=xty(j)-ytilde/winv
            end do
         end do
         
         call invdet(xtx,p-1,workmp1,detlog,iflagp,workvp1)

         do i=1,p-1
            tmp1=0.d0
            workvp1(i)=0.d0
            do j=1,p-1
               tmp1=tmp1+workmp1(i,j)*xty(j) 
            end do
            workvp1(i)=tmp1
         end do
            
         call dmvn2(p-1,beta,workvp1,workmp1,logcgkn,
     &              workvp2,workmp2,workmp3,iflagp)                 

 
c+++++++ mh step

         ratio=dexp(loglikn-logliko+logcgkn-logcgko+
     &              logpriorn-logprioro)

         if(dble(runif()).lt.ratio)then
            acrate(1)=acrate(1)+1.d0
            do i=1,p-1
               beta(i)=betac(i) 
            end do
         end if

c         call dblepr("ratio",-1,ratio,1)
c         call dblepr("loglikn",-1,loglikn,1)
c         call dblepr("logliko",-1,logliko,1)
c         call dblepr("logcgkn",-1,logcgkn,1)
c         call dblepr("logcgko",-1,logcgko,1)
c         call dblepr("logpriorn",-1,logpriorn,1)
c         call dblepr("logprioro",-1,logprioro,1)
c         call dblepr("beta",-1,beta,p-1)
c         call dblepr("betac",-1,betac,p-1)

c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         acrate2=0.d0

         do i=1,nsubject

            tmp1=0.d0
            do j=1,q
               tmp1=tmp1+z(i,j)*alpha(ss(i),j)
            end do

            logprioro=dnrm(b(i),tmp1,sqrt(sigma2b),1)

            ztz=1.d0/sigma2b
            zty=tmp1/sigma2b

            logliko=0.d0
            
            do j=1,p
               if(j.eq.1)then
                  eta=b(i)
                  offset=0.d0
                 else
                  eta=b(i)-beta(j-1)
                  offset=-beta(j-1)
               end if
               
               yij=y(i,j)
               mean=exp(eta)/(1.d0+exp(eta))
               gprime=1.d0/(mean*(dble(1.d0)-mean))
               winv=1.d0/(mean*(dble(1.d0)-mean))

               ytilde=eta+(dble(yij)-mean)*gprime-offset    

               logliko=logliko+dbin(dble(yij),1.d0,mean,1)
               
               ztz=ztz+1.d0/winv
               zty=zty+ytilde/winv
            end do

            ztzinv=1.d0/ztz
            tmp1=ztzinv*zty
 
            thetac=rnorm(tmp1,sqrt(ztzinv))

            logcgko=dnrm(thetac,tmp1,sqrt(ztzinv),1)

c++++++++++ candidate generating kernel contribution

            tmp1=0.d0
            do j=1,q
               tmp1=tmp1+z(i,j)*alpha(ss(i),j)
            end do

            logpriorn=dnrm(thetac,tmp1,sqrt(sigma2b),1)

            ztz=1.d0/sigma2b
            zty=tmp1/sigma2b

            loglikn=0.d0
            
            do j=1,p
               if(j.eq.1)then
                  eta=thetac
                  offset=0.d0
                 else
                  eta=thetac-beta(j-1)
                  offset=-beta(j-1)
               end if
               
               yij=y(i,j)
               mean=exp(eta)/(1.d0+exp(eta))
               gprime=1.d0/(mean*(dble(1.d0)-mean))
               winv=1.d0/(mean*(dble(1.d0)-mean))

               ytilde=eta+(dble(yij)-mean)*gprime-offset    

               loglikn=loglikn+dbin(dble(yij),1.d0,mean,1)
               
               ztz=ztz+1.d0/winv
               zty=zty+ytilde/winv
            end do

            ztzinv=1.d0/ztz
            tmp1=ztzinv*zty
 
            logcgkn=dnrm(b(i),tmp1,sqrt(ztzinv),1)

c++++++++++ mh step
            ratio=dexp(loglikn-logliko+logcgkn-logcgko+
     &                 logpriorn-logprioro)

            if(dble(runif()).lt.ratio)then
               acrate2=acrate2+1.d0
               b(i)=thetac
            end if
         end do

         acrate(2)=acrate(2)+acrate2/dble(nsubject)


c++++++++++++++++++++++++++++++++++         
c+++++++ DDP part
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn 
c++++++++++++++++++++++++++++++

         do i=1,nsubject
         
            ns=ccluster(ss(i))

c++++++++++ observation in cluster with more than 1 element
             
            if(ns.gt.1)then

               ccluster(ss(i))=ccluster(ss(i))-1 

               do j=1,ncluster
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+z(i,k)*alpha(j,k)
                  end do

                  tmp2=dnrm(b(i),tmp1,sqrt(sigma2b),1)

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        tmp2)
               end do

               call rmvnorm(q,mu,sigma,workmhq1,workvq1,alphac)
               
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+z(i,k)*alphac(k)
               end do
               tmp2=dnrm(b(i),tmp1,sqrt(sigma2b),1)
               prob(ncluster+1)=exp(log(m)+tmp1)
               
               call simdisc(prob,nsubject+2,ncluster+1,evali)

               if(evali.le.ncluster)then
                  ss(i)=evali
                  ccluster(evali)=ccluster(evali)+1
               end if   
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  ss(i)=ncluster
                  ccluster(ncluster)=1
	          do j=1,q
	             alpha(ncluster,j)=alphac(j)
	          end do
               end if               
            end if

c++++++++++ subject in cluster with only 1 observation
             
            if(ns.eq.1)then
                
               since=ss(i)
               
               if(since.lt.ncluster)then
                   call relabelDDP(i,since,nsubject,q,ncluster,
     &                             ccluster,ss,alpha,alphac)
	       end if

               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               do j=1,ncluster
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+z(i,k)*alpha(j,k)
                  end do

                  tmp2=dnrm(b(i),tmp1,sqrt(sigma2b),1)

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        tmp2)
               end do
               
               do j=1,q
                  alphac(j)=alpha(ncluster+1,j)
               end do

               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+z(i,k)*alphac(k)
               end do
               tmp2=dnrm(b(i),tmp1,sqrt(sigma2b),1)
               prob(ncluster+1)=exp(log(m)+tmp1)
               
               call simdisc(prob,nsubject+2,ncluster+1,evali)

               if(evali.le.ncluster)then
                  ss(i)=evali
                  ccluster(evali)=ccluster(evali)+1
               end if   
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  ss(i)=ncluster
                  ccluster(ncluster)=1
	          do j=1,q
	             alpha(ncluster,j)=alphac(j)
	          end do
               end if      
	    
	    end if

         end do

c++++++++++++++++++++++++++++++
c+++++++ b) Resampling step
c++++++++++++++++++++++++++++++

         acrate2=0.d0
         
         do ii=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do i=1,q
               tmp1=0.d0
               do j=1,q
                  wtw(i,j)=sigmainv(i,j)
                  tmp1=sigmainv(i,j)*mu(j)
               end do
               wty(i)=tmp1
            end do

            do i=1,nsubject
               if(ii.eq.ss(i))then
                  do j=1,q
                     do k=1,q
                        tmp1=z(i,j)*z(i,k)/sigma2b
                        wtw(j,k)=wtw(j,k)+tmp1
                     end do
                  end do
                  do j=1,q
                     tmp1=z(i,j)*b(i)/sigma2b
                     wty(j)=wty(j)+tmp1
                  end do                  
               end if
            end do
            
            call invdet(wtw,q,workmq1,detlog,iflagq,workvq1)

            do i=1,q
               tmp1=0.d0
               workvq1(i)=0.d0
               do j=1,q
                  tmp1=tmp1+workmq1(i,j)*wty(j) 
               end do
               workvq1(i)=tmp1
            end do
            
            call rmvnorm(q,workvq1,workmq1,workmhq1,workvq2,alphac)

            do i=1,q
               alpha(ii,i)=alphac(i)
            end do
         end do

c++++++++++++++++++++++++++++++++++         
c+++++++ sigma2b
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(tau1.gt.0.d0)then
            ssb=0.d0 
            do i=1,nsubject
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+z(i,j)*alpha(ss(i),j)
               end do
               tmp2=b(i)-tmp1
               ssb=ssb+tmp2*tmp2
            end do

            sigma2b=1.d0/
     &           rgamma(0.5d0*(dble(nsubject)+tau1),0.5d0*(ssb+tau2))
     
         end if


c++++++++++++++++++++++++++++++++++         
c+++++++ Baseline mean
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()
         
         if(murand.eq.1)then

            do i=1,q
               workvq1(i)=smu(i)
               do j=1,q
                  workmq1(i,j)=(sigmainv(i,j)*dble(ncluster))+
     &                          prec2(i,j)
               end do
            end do

            call invdet(workmq1,q,workmq2,detlog,iflagq,workvq2)

            do i=1,ncluster
               do j=1,q
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+sigmainv(j,k)*alpha(i,k)
                  end do
                  workvq1(j)=workvq1(j)+tmp1
               end do
            end do
     
            do i=1,q
               workvq2(i)=0.d0
            end do
     
            do i=1,q
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+workmq2(i,j)*workvq1(j)
               end do
               workvq2(i)=tmp1
            end do
          
            call rmvnorm(q,workvq2,workmq2,workmhq2,workvq1,mu)

         end if

         
c++++++++++++++++++++++++++++++++++         
c+++++++ Baseline variance
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(nu0.gt.0.d0)then

            do i=1,q
               do j=1,q
                  quadf(i,j)=0.d0
               end do
            end do
         
            do i=1,ncluster
               do j=1,q
                  do k=1,q
                     quadf(j,k)=quadf(j,k)+               
     &                         (alpha(i,j)-mu(j))*(alpha(i,k)-mu(k))                   
                  end do
               end do
            end do

            do i=1,q
               do j=1,q
                  quadf(i,j)=quadf(i,j)+tinv(i,j)
               end do
            end do


            call riwishart(q,int(nu0+ncluster),quadf,workmq1,workmq2,
     &                     workvq1,workmhq1,workmhq2,iflagq)

            do i=1,q
               do j=1,q
                  sigma(i,j)=quadf(i,j)
                  sigmainv(i,j)=workmq1(i,j)
               end do
            end do

         end if 


c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++
         if(aa0.gt.0.d0)then
            call samalph(m,aa0,ab0,ncluster,nsubject)
         end if 


9999     continue        

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ difficulty parameters

               do i=1,p-1
                  thetasave(isave,i)=beta(i)
               end do

c+++++++++++++ random effects variance
               thetasave(isave,p)=sigma2b


c+++++++++++++ baselihne information
               
               do i=1,q
                  thetasave(isave,p+i)=mu(i)
               end do   
              
               counter=0
               do i=1,q
                  do j=1,i
                     counter=counter+1
                     thetasave(isave,p+q+counter)=sigma(i,j)
                  end do   
               end do

c+++++++++++++ cluster information

               thetasave(isave,p+q+counter+1)=ncluster
               thetasave(isave,p+q+counter+2)=m

c+++++++++++++ random effects

               do i=1,nsubject
                  randsave(isave,i)=b(i)
               end do

c+++++++++++++ predictive information

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(m+dble(nsubject))
               end do
               prob(ncluster+1)=m/(m+dble(nsubject))
               
               call rmvnorm(q,mu,sigma,workmhq1,workvq1,alphac)

               do i=1,npred

                  call simdisc(prob,nsubject+2,ncluster+1,evali)

                  do j=1,ncluster
                     tmp1=0.d0 
                     do k=1,q
                        tmp1=tmp1+zpred(i,k)*alpha(j,k) 
                     end do
                     
                     do k=1,ngrid
                        tmp2=prob(j)*
     &                       dnrm(grid(k),tmp1,sqrt(sigma2b),0)
                       densr(i,k)=densr(i,k)+tmp2
                     end do   
                  end do

                  tmp1=0.d0 
                  do k=1,q
                     tmp1=tmp1+zpred(i,k)*alphac(k) 
                  end do                
                  
                  do k=1,ngrid
                     tmp2=prob(ncluster+1)*
     &                    dnrm(grid(k),tmp1,sqrt(sigma2b),0)
                     densr(i,k)=densr(i,k)+tmp2
                  end do   

                  
                  if(evali.le.ncluster)then
                     tmp1=0.d0 
                     do k=1,q
                        tmp1=tmp1+zpred(i,k)*alpha(evali,k) 
                     end do
                    else
                     tmp1=0.d0 
                     do k=1,q
                        tmp1=tmp1+zpred(i,k)*alphac(k) 
                     end do
                  end if                      
                  thetac=rnorm(tmp1,sqrt(sigma2b))
                  randsave(isave,nsubject+i)=thetac
               end do
               

c+++++++++++++ cpo

               do i=1,nsubject
                  do j=1,p
                     yij=y(i,j)
                     if(j.eq.1)then
                       eta=b(i)
                      else
                       eta=b(i)-beta(j-1)
                     end if  
                     tmp1=dbin(dble(yij),1.d0,mean,0)
                     cpo(i,j)=cpo(i,j)+1.0d0/tmp1
                  end do
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
      
      do i=1,2
         acrate(i)=acrate(i)/dble(nscan)    
      end do   
      
      do i=1,nsubject
         do j=1,p
            cpo(i,j)=dble(nsave)/cpo(i,j)
         end do   
      end do

      do i=1,npred
         do j=1,ngrid
            densr(i,j)=densr(i,j)/dble(nsave)
         end do   
      end do

      return
      end


c=======================================================================      
      subroutine relabelDDP(ind,since,nsubject,q,ncluster,ccluster,ss,
     &                      alpha,alphac)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2007
      implicit none
      integer i,j,ind,since,nsubject,q,ncluster,ccluster(nsubject)
      integer ss(nsubject)
      real*8 alpha(nsubject,q),alphac(q)

      do i=1,q
         alphac(i)=alpha(since,i)
      end do
      
      do i=since+1,ncluster
         do j=1,nsubject
            if(ss(j).eq.i)then
               ss(j)=i-1 
            end if
         end do
         do j=1,q
            alpha(i-1,j)=alpha(i,j)
         end do
         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      
      do i=1,q
         alpha(ncluster,i)=alphac(i)
      end do
      ccluster(ncluster)=1
      
      return
      end  
      