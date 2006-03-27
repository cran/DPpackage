
c=======================================================================                      
      subroutine sppoislme(datastr,maxni,nrec,nsubject,nfixed,p,q,
     &                     subject,x,y,z,roffset,
     &                     a0b0,b0,nu0,prec,psiinv,sb,smu,tinv,
     &                     mcmc,nsave,
     &                     acrate,cpo,randsave,thetasave,
     &                     alpha,b,bclus,beta,betar,mu,ncluster,sigma,
     &                     sigmainv,ss,
     &                     betac,ccluster,iflag,iflagb,prob,quadf,seed,
     &                     theta,thetac,work1,work2,work3,workb1,
     &                     workb2,workmh1,workmh2,workmh3,workv1,workv2,
     &                     workv3,workvb1,workvb2,workvb3,xtx,xty,
     &                     zty,ztz,
     &                     ztzinv)
c=======================================================================                      
c
c     Version 1.0: 
c     Last modification: 27-04-2006.
c
c     Subroutine `sppoislme' to run a Markov chain in the  
c     semiparametric poison mixed model. In this routine, inference 
c     is based on  the Polya urn representation of Dirichlet process.
c     The algorithm 8 with m=1 of Neal (2000) is used to sample the 
c     configurations.
c
c     Author: A.J.V.
c
c---- Data -------------------------------------------------------------
c 
c        datastr     :  integer matrix giving the number of measurements
c                       and the location in y of the observations for 
c                       each subject, datastr(nsubject,maxni+1)
c        maxni       :  integer giving the maximum number of 
c                       measurements for subject.
c        nrec        :  integer giving the number of observations.
c        nsubject    :  integer giving the number of subjects.
c        nfixed      :  integer giving the number of fixed effects,
c                       if nfixed is 0 then p=1.
c        p           :  integer giving the number of fixed coefficients.
c        q           :  integer giving the number of random effects.
c        roffset     :  real vector giving the real offset for each
c                       observation.
c        subject     :  integer vector giving the subject for each.
c                       observation, subject(nsubject).
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        y           :  integer matrix giving the response variable,
c                       y(nrec).
c        z           :  real matrix giving the design matrix for the 
c                       random effects, z(nrec,q). 
c-----------------------------------------------------------------------
c
c---- Prior information ------------------------------------------------
c 
c        aa0, ab0    :  real giving the hyperparameters of the prior
c                       distribution for the precision parameter,
c                       alpha ~ Gamma(aa0,ab0). If aa0<0 the precision 
c                       parameter is considered as a constant.
c        b0          :  real vector giving the prior mean of fixed
c                       effects, b0(p).
c        nu0         :  integer giving the degrees of freedom for the
c                       inverted-Wishart prior distribution for the
c                       covariance matrix of the random effects
c                       (This is for the base line).
c        prec        :  real matrix giving the prior precision matrix
c                       for the fixed effects, prec(p,p).
c        psiinv      :  real matrix giving the prior precision matrix
c                       for the baseline mean, psiinv(q,q).
c        sb          :  real vector giving the product of the prior 
c                       precision and prior mean for the fixed effects,
c                       sb(p).
c        smu         :  real vector giving the product of the prior 
c                       precision and prior mean for the baseline mean,
c                       smu(q).
c        tinv        :  real matrix giving the scale matrix for the
c                       inverted-Wishart prior distribution for the
c                       covariance matrix of the random effects, 
c                       sigma ~ Inv-Wishart(nu0,tinv^{-1}), such that 
c                       E(sigma)=(1/(nu0-q-1)) * tinv 
c                       (This is for the base line distribution)
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
c        acrate      :  real vector giving the MH acceptance rate. 
c        cpo         :  real giving the cpo, acrate(2).
c        randsave    :  real matrix containing the mcmc samples for
c                       the random effects and prediction,
c                       randsave(nsave,q*(nsubject+1))
c                       thetsave(nsave,q+nfixed+1+q+nuniq(Sigma)+2).
c        thetasave   :  real matrix containing the mcmc samples for
c                       the averaged random effects, fixed effects, 
c                       error variance, and mean and covariance of
c                       the baseline distribution, 
c                       thetsave(nsave,q+nfixed+q+nuniq(Sigma)+2).
c
c-----------------------------------------------------------------------
c
c---- Current value of the parameters ----------------------------------
c
c        alpha       :  real giving the current value of the precision
c                       parameter of the Dirichlet process.
c        b           :  real matrix giving the current value of the 
c                       random effects, b(nsubject,q).
c        bclus       :  real matrix giving the current value of the 
c                       different values of random effects, 
c                       bclus(nsubject,q).
c        beta        :  real vector giving the current value of the 
c                       fixed effects, beta(p).
c        betar       :  real vector giving the current value of the 
c                       averaged random effects, betar(q).
c        mu          :  real vector giving the mean of the normal 
c                       base line distribution for the random effects,
c                       mu(q).
c        ncluster    :  integer giving the number of clusters in the
c                       random effects.
c        sigma       :  real matrix giving the current value of the
c                       covariance matrix for normal base line 
c                       distribution for the random effects,
c                       sigma(q,q).
c        sigmainv    :  real matrix used to save the base line 
c                       covariance matrix for the random effects,
c                       sigmainv(q,q).
c        ss          :  integer vector giving the cluster label for 
c                       each subject, ss(nsubject).
c-----------------------------------------------------------------------
c
c---- Working space ----------------------------------------------------
c
c        acrate2     :  real used to calculate the acceptance rate. 
c        betac       :  real vector giving the current value of the 
c                       candidate for fixed effects, betac(p).
c        ccluster    :  integer vector indicating the number of
c                       subjects in each cluster, ccluster(nsubject).
c        detlog      :  real used to save the log-determinant in a
c                       matrix inversion process.
c        dispcount   :  index. 
c        evali       :  integer indicator used in updating the state.
c        i           :  index. 
c        ii          :  index. 
c        iflag       :  integer vector used to invert the of the lhs
c                       least square solution for the fixed effects,
c                       iflag(p).
c        iflagb      :  integer vector used to invert the of the lhs
c                       least square solution for the random effects,
c                       iflagb(q).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index. 
c        l           :  index.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nsubject+2).
c        quadf       :  real matrix used to save the bilinear product
c                       of random effects, quadf(q,q).
c        ni          :  integer indicator used in updating the state. 
c        ns          :  integer indicator used in updating the state. 
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        runif       :  uniform random number generator.
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        seed3       :  seed for random number generation.
c        since       :  index.
c        skipcount   :  index. 
c        theta       :  real vector used to save randomnly generated
c                       random effects, theta(q).
c        thetac      :  real vector used to save randomnly generated
c                       random effects, thetac(q).
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        work1       :  real matrix used to update the fixed effects,
c                       work1(p,p).
c        work2       :  real matrix used to update the fixed effects,
c                       work2(p,p).
c        work3       :  real matrix used to update the fixed effects,
c                       work3(p,p).
c        workb1      :  real matrix used to update the random effects,
c                       workb1(q,q).
c        workb2      :  real matrix used to update the random effects,
c                       workb2(q,q).
c        workmh1     :  real vector used to update the fixed effects,
c                       workmh1(p*(p+1)/2).
c        workmh2     :  real vector used to update the random effects,
c                       workmh2(q*(q+1)/2).
c        workmh3     :  real vector used to update the random effects,
c                       workmh3(q*(q+1)/2).
c        workv1      :  real vector used to update the fixed effects,
c                       workv1(p).
c        workv2      :  real vector used to update the fixed effects,
c                       workv2(p).
c        workv3      :  real vector used to update the fixed effects,
c                       workv3(p).
c        workvb1     :  real vector used to update the random effects,
c                       workvb1(q).
c        workvb2     :  real vector used to update the random effects,
c                       workvb2(q).
c        workvb3     :  real vector used to update the random effects,
c                       workvb3(q).
c        xtx         :  real matrix givind the product X^tX, xtx(p,p).
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p).
c        zty         :  real vector used to save the product 
c                       Zt(Y-Xbeta), zty(q).
c        ztz         :  real matrix used to save the product 
c                       ZtSigma^1Z, ztz(q,q).
c        ztzinv      :  real matrix used to save the inverted 
c                       ztz, ztzinv(q,q).
c=======================================================================                  
      implicit none 

c+++++Data
      integer maxni,nrec,nsubject,nfixed,p,q,subject(nrec)
      integer datastr(nsubject,maxni+1),y(nrec)
      real*8 roffset(nrec),x(nrec,p),z(nrec,q)	
      
c+++++Prior 
      integer nu0
      real*8 aa0,ab0,a0b0(2),b0(p),prec(p,p),psiinv(q,q)
      real*8 sb(p),smu(q)
      real*8 tinv(q,q)      

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      real*8 acrate(2)
      real*8 cpo(nrec)
      real*8 randsave(nsave,q*(nsubject+1))
      real*8 thetasave(nsave,q+nfixed+q+(q*(q+1)/2)+2)

c+++++Current values of the parameters
      integer ncluster,ss(nsubject)
      real*8 alpha,beta(p),b(nsubject,q)
      real*8 betar(q),bclus(nsubject,q)
      real*8 mu(q),sigma(q,q),sigmainv(q,q)

c+++++Working space - Loops
      integer ii,i,j,k,l

c+++++Working space - Random effects
      integer iflagb(q)
      real*8 quadf(q,q)
      real*8 thetac(q)
      real*8 zty(q),ztz(q,q),ztzinv(q,q)
      real*8 workb1(q,q),workb2(q,q)
      real*8 workmh2(q*(q+1)/2),workmh3(q*(q+1)/2)
      real*8 workvb1(q),workvb2(q),workvb3(q)

c+++++Working space - RNG
      integer seed(3),seed1,seed2,seed3
      real runif

c+++++Working space - MCMC
      integer iscan,isave,nscan
      integer sprint,skipcount,dispcount
      
c+++++Working space - Configurations
      integer ccluster(nsubject),evali 
      integer ni,ns
      integer since
      real*8 dgamlog,prob(nsubject+2)
      real*8 tmp1,tmp2
      real*8 theta(q)

c+++++Working space - Fixed effects
      integer iflag(p)
      real*8 betac(p)
      real*8 detlog
      real*8 xtx(p,p),xty(p)
      real*8 workmh1(p*(p+1)/2)
      real*8 work1(p,p),work2(p,p),work3(p,p)
      real*8 workv1(p),workv2(p),workv3(p)

c+++++Working space - GLM part
      integer yij
      real*8 acrate2
      real*8 eta,etac,gprime,gprimec,mean,meanc,offset,ytilde,ytildec
      real*8 logp

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
      aa0=a0b0(1)
      ab0=a0b0(2)
      
c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      seed3=seed(3)
      call setrand(seed1,seed2,seed3)
     
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
      
      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++
c+++++++ fixed effects
c++++++++++++++++++++++++++++++++++

         if(nfixed.eq.0)go to 1
         
            do i=1,p
               do j=1,p
                  xtx(i,j)=0.d0
                  work1(i,j)=0.d0
                  work2(i,j)=0.d0
                  work3(i,j)=0.d0
               end do
               xty(i)=sb(i)
               workv1(i)=0.d0
               workv2(i)=0.d0
               workv3(i)=0.d0
               iflag(i)=0
            end do
            
            do i=1,nrec
               eta=0.d0
               offset=0.d0
               mean=0.d0
               gprime=0.d0
               
               do j=1,p
                  eta=eta+x(i,j)*beta(j)
               end do
               
               do j=1,q
                  eta=eta+z(i,j)*b(subject(i),j) 
                  offset=offset+z(i,j)*b(subject(i),j) 
               end do
               
               eta=eta+roffset(i)
               
               offset=offset+roffset(i)
               
               mean=exp(eta)
               
               gprime=exp(-eta)
               
               ytilde=eta+(dble(y(i))-mean)*gprime-offset
               
               do j=1,p
                  do k=1,p
                     xtx(j,k)=xtx(j,k)+x(i,j)*x(i,k)/gprime
                  end do
                  xty(j)=xty(j)+x(i,j)*ytilde/gprime
               end do
            end do

            do i=1,p
               do j=1,p
                  work1(i,j)=xtx(i,j)+prec(i,j)          
               end do
            end do

            call invdet(work1,p,work2,detlog,iflag,workv1)
            
            do i=1,p
               tmp1=0.d0
               do j=1,p
c                  work2(i,j)=work2(i,j)*1.5
                  tmp1=tmp1+work2(i,j)*xty(j) 
               end do
               workv2(i)=tmp1
            end do

            
            call rmvnorm(p,workv2,work2,workmh1,workv3,betac)

            call dmvn(p,betac,workv2,work2,tmp1,
     &                workv1,work1,work3,workv3,iflag)                 

            logp=0.d0            
            logp=logp-tmp1
            

c++++++++++ likelihood ratio

            do i=1,nrec
               eta=0.d0
               etac=0.d0
               
               do j=1,p
                  eta=eta+x(i,j)*beta(j)
                  etac=etac+x(i,j)*betac(j)
               end do
               
               do j=1,q
                  eta=eta+z(i,j)*b(subject(i),j) 
                  etac=etac+z(i,j)*b(subject(i),j) 
               end do
               
               eta=eta+roffset(i)
               etac=etac+roffset(i)

               logp=logp+
     &              dble(y(i))*(etac-eta)-exp(etac)+exp(eta)
     
            end do

c++++++++++ prior ratio

            tmp1=0.d0
            tmp2=0.d0
            
            do i=1,p
               do j=1,p
                  tmp1=tmp1+(betac(i)-b0(i))* 
     &                       prec(i,j)      *
     &                      (betac(j)-b0(j))

                  tmp2=tmp2+(beta(i) -b0(i))* 
     &                       prec(i,j)      *
     &                      (beta(j) -b0(j))
               end do
            end do
            
            logp=logp-0.5d0*tmp1+0.5d0*tmp2
            
c++++++++++ candidate generating kernel contribution

            do i=1,p
               do j=1,p
                  xtx(i,j)=0.d0
                  work1(i,j)=0.d0
                  work2(i,j)=0.d0
                  work3(i,j)=0.d0
               end do
               xty(i)=sb(i)
               workv1(i)=0.d0
               workv2(i)=0.d0
               workv3(i)=0.d0
               iflag(i)=0
            end do
        
            do i=1,nrec
               etac=0.d0
               offset=0.d0
               meanc=0.d0
               gprimec=0.d0
               
               do j=1,p
                  etac=etac+x(i,j)*betac(j)
               end do
               
               do j=1,q
                  etac=etac+z(i,j)*b(subject(i),j) 
                  offset=offset+z(i,j)*b(subject(i),j) 
               end do

               etac=etac+roffset(i)

               offset=offset+roffset(i)
               
               meanc=exp(etac)
               
               gprimec=exp(-etac)
               
               ytildec=etac+(dble(y(i))-meanc)*gprimec-offset
               
               do j=1,p
                  do k=1,p
                     xtx(j,k)=xtx(j,k)+x(i,j)*x(i,k)/gprimec
                  end do
                  xty(j)=xty(j)+x(i,j)*ytildec/gprimec
               end do
            end do

            do i=1,p
               do j=1,p
                  work1(i,j)=xtx(i,j)+prec(i,j)          
               end do
            end do

            call invdet(work1,p,work2,detlog,iflag,workv1)

            do i=1,p
               tmp1=0.d0
               do j=1,p
c                  work2(i,j)=work2(i,j)*1.5
                  tmp1=tmp1+work2(i,j)*xty(j) 
               end do
               workv2(i)=tmp1
            end do
            
            call dmvn(p,beta,workv2,work2,tmp1,
     &                workv1,work1,work3,workv3,iflag)                 
 
            logp=logp+tmp1

c++++++++++ mh step

c            call dblepr("logp",-1,logp,1)  

            if(log(dble(runif())).lt.logp)then
               acrate(1)=acrate(1)+1.d0
               do i=1,p
                  beta(i)=betac(i) 
               end do
            end if
            
1        continue            

        
         
c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn 
c++++++++++++++++++++++++++++++


         do i=1,nsubject
         
            ns=ccluster(ss(i))
            ni=datastr(i,1) 

c++++++++++ observation in cluster with more than 1 element
             
            if(ns.gt.1)then
 
               ccluster(ss(i))=ccluster(ss(i))-1 
               
               do j=1,ncluster
                  tmp1=0.d0
                  do k=1,ni
                     yij=y(datastr(i,k+1))
                     
                     eta=0.d0
                     do l=1,p
                        eta=eta+x(datastr(i,k+1),l)*beta(l)
                     end do
		     do l=1,q
		        eta=eta+z(datastr(i,k+1),l)*bclus(j,l)
                     end do
                     
                     eta=eta+roffset(datastr(i,k+1))

                     tmp1=tmp1+dble(yij)*eta-exp(eta)-
     &                         dgamlog(dble(yij+1))

                  end do

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        tmp1)
               end do
               
               call rmvnorm(q,mu,sigma,workmh2,workvb1,theta)
               
               tmp1=0.d0
               do k=1,ni
                  yij=y(datastr(i,k+1))
                      
                  eta=0.d0
                  do l=1,p
                     eta=eta+x(datastr(i,k+1),l)*beta(l)
                  end do
 		  do l=1,q
 		     eta=eta+z(datastr(i,k+1),l)*theta(l)
                  end do
                  eta=eta+roffset(datastr(i,k+1))
                      
                  tmp1=tmp1+dble(yij)*eta-exp(eta)-
     &                      dgamlog(dble(yij+1))
               end do
 
               prob(ncluster+1)=exp(log(alpha)+tmp1)
               
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
	             bclus(ncluster,j)=theta(j)
	          end do
               end if               
            end if

c++++++++++ subject in cluster with only 1 observation
             
            if(ns.eq.1)then
                
               since=ss(i)
                
               if(since.lt.ncluster)then
                   call relabelg(i,since,nsubject,q,ncluster,
     &                          ccluster,ss,bclus,theta)                   
	       end if

               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               do j=1,ncluster
                  tmp1=0.d0
                  do k=1,ni
                     yij=y(datastr(i,k+1))
                     
                     eta=0.d0
                     do l=1,p
                        eta=eta+x(datastr(i,k+1),l)*beta(l)
                     end do
		     do l=1,q
		        eta=eta+z(datastr(i,k+1),l)*bclus(j,l)
                     end do
                     
                     eta=eta+roffset(datastr(i,k+1))

                     tmp1=tmp1+dble(yij)*eta-exp(eta)-
     &                         dgamlog(dble(yij+1))                     

                  end do

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        tmp1)
               end do

               call rmvnorm(q,mu,sigma,workmh2,workvb1,theta)
               
               tmp1=0.d0
               do k=1,ni
                  yij=y(datastr(i,k+1))
                      
                  eta=0.d0
                  do l=1,p
                     eta=eta+x(datastr(i,k+1),l)*beta(l)
                  end do
 		  do l=1,q
 		     eta=eta+z(datastr(i,k+1),l)*theta(l)
                  end do
                  
                  eta=eta+roffset(datastr(i,k+1))
                      
                  tmp1=tmp1+dble(yij)*eta-exp(eta)-
     &                      dgamlog(dble(yij+1))                  

               end do
 
               prob(ncluster+1)=exp(log(alpha)+tmp1)
               
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
	             bclus(ncluster,j)=theta(j)
	          end do
               end if      
	    
	    end if

         end do


c         call intpr("ncluster",-1,ncluster,1)  

c++++++++++++++++++++++++++++++
c+++++++ b) Resampling step
c++++++++++++++++++++++++++++++

         acrate2=0.d0
         
         do ii=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do i=1,q
               do j=1,q
                  ztz(i,j)=0.d0
                  ztzinv(i,j)=0.d0
                  workb1(i,j)=0.d0
                  workb2(i,j)=0.d0
               end do
               zty(i)=0.d0
               workvb1(i)=0.d0
               workvb2(i)=0.d0
               workvb3(i)=0.d0
               iflagb(i)=0
            end do
            
            do i=1,nsubject
               if(ii.eq.ss(i))then
                  ni=datastr(i,1) 
                  do j=1,ni
                     eta=0.d0
                     offset=0.d0
                     mean=0.d0
                     gprime=0.d0

                     yij=y(datastr(i,j+1))
                     
                     do k=1,p
                        eta=eta+x(datastr(i,j+1),k)*beta(k)
                        offset=offset+x(datastr(i,j+1),k)*beta(k)
                     end do
                     do k=1,q
                        eta=eta+z(datastr(i,j+1),k)*bclus(ii,k)
                     end do
                   
                     eta=eta+roffset(datastr(i,j+1))
                     
                     offset=offset+roffset(datastr(i,j+1))

                     mean=exp(eta)
               
                     gprime=exp(-eta)

                     ytilde=eta+(dble(yij)-mean)*gprime-offset
                     
                     do k=1,q
                        do l=1,q
                           ztz(k,l)=ztz(k,l)+z(datastr(i,j+1),k)*
     &                                       z(datastr(i,j+1),l)/gprime
                        end do
                        zty(k)=zty(k)+z(datastr(i,j+1),k)*
     &                         ytilde/gprime
                     end do
                  end do
               end if
            end do
    
            do i=1,q
               do j=1,q
                  ztz(i,j)=ztz(i,j)+sigmainv(i,j)
               end do
            end do

            call invdet(ztz,q,ztzinv,detlog,iflagb,workvb2)

            do i=1,q
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+sigmainv(i,j)*mu(j)   
               end do
               zty(i)=zty(i)+tmp1
            end do
          
            do i=1,q
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+ztzinv(i,j)*zty(j) 
               end do
               workvb1(i)=tmp1
            end do
  
            call rmvnorm(q,workvb1,ztzinv,workmh2,workvb2,thetac)

            call dmvn(q,thetac,workvb1,ztzinv,tmp1,
     &                workvb2,workb1,workb2,workvb3,iflagb)                 
 
            logp=0.d0
            logp=logp-tmp1


c++++++++++ likelihood ratio

            do i=1,nsubject
               if(ii.eq.ss(i))then
                  ni=datastr(i,1) 
                  do j=1,ni
                     eta=0.d0
                     etac=0.d0

                     yij=y(datastr(i,j+1))
                     
                     do k=1,p
                        eta=eta+x(datastr(i,j+1),k)*beta(k)
                        etac=etac+x(datastr(i,j+1),k)*beta(k)
                     end do
                     do k=1,q
                        eta=eta+z(datastr(i,j+1),k)*bclus(ii,k)
                        etac=etac+z(datastr(i,j+1),k)*thetac(k)
                     end do
                     
                     eta=eta+roffset(datastr(i,j+1))
                     etac=etac+roffset(datastr(i,j+1))

                     logp=logp+
     &                    dble(yij)*(etac-eta)-exp(etac)+exp(eta)                     
                  end do
               end if
            end do
 
c++++++++++ prior ratio
 
            tmp1=0.d0
            tmp2=0.d0
            
            do i=1,q
               do j=1,q
                  tmp1=tmp1+(thetac(i)-mu(i))* 
     &                       sigmainv(i,j)      *
     &                      (thetac(j)-mu(j))

                  tmp2=tmp2+(bclus(ii,i) -mu(i))* 
     &                       sigmainv(i,j)      *
     &                      (bclus(ii,j) -mu(j))
               end do
            end do
            
            logp=logp-0.5d0*tmp1+0.5d0*tmp2 

c++++++++++ candidate generating kernel contribution


            do i=1,q
               do j=1,q
                  ztz(i,j)=0.d0
                  ztzinv(i,j)=0.d0
                  workb1(i,j)=0.d0
                  workb2(i,j)=0.d0
               end do
               zty(i)=0.d0
               workvb1(i)=0.d0
               workvb2(i)=0.d0
               workvb3(i)=0.d0
               iflagb(i)=0
            end do
   
                 
            do i=1,nsubject
               if(ii.eq.ss(i))then
                  ni=datastr(i,1) 
                  do j=1,ni
                     etac=0.d0
                     offset=0.d0
                     meanc=0.d0
                     gprimec=0.d0

                     yij=y(datastr(i,j+1))
                     
                     do k=1,p
                        etac=etac+x(datastr(i,j+1),k)*beta(k)
                        offset=offset+x(datastr(i,j+1),k)*beta(k)
                     end do
                     do k=1,q
                        etac=etac+z(datastr(i,j+1),k)*thetac(k)
                     end do

                     etac=etac+roffset(datastr(i,j+1))
                     
                     offset=offset+roffset(datastr(i,j+1))

                     meanc=exp(etac)
 
                     gprimec=exp(-etac)

                     ytildec=etac+(dble(yij)-meanc)*gprimec-offset
                     
                     do k=1,q
                        do l=1,q
                           ztz(k,l)=ztz(k,l)+z(datastr(i,j+1),k)*
     &                          z(datastr(i,j+1),l)/gprimec
                        end do
                        zty(k)=zty(k)+z(datastr(i,j+1),k)*
     &                         ytildec/gprimec
                     end do
                  end do
               end if
            end do
    
            do i=1,q
               do j=1,q
                  ztz(i,j)=ztz(i,j)+sigmainv(i,j)
               end do
            end do

            call invdet(ztz,q,ztzinv,detlog,iflagb,workvb2)

            do i=1,q
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+sigmainv(i,j)*mu(j)   
               end do
               zty(i)=zty(i)+tmp1
            end do
          
            do i=1,q
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+ztzinv(i,j)*zty(j) 
               end do
               workvb1(i)=tmp1
               theta(i)=bclus(ii,i)
            end do
  
            call dmvn(q,theta,workvb1,ztzinv,tmp1,
     &                workvb2,workb1,workb2,workvb3,iflagb)                 
 
            logp=logp+tmp1


c++++++++++ mh step

            if(log(dble(runif())).lt.logp)then
               acrate2=acrate2+1.d0
               do i=1,q
                  bclus(ii,i)=thetac(i) 
               end do
            end if
         end do

         acrate(2)=acrate(2)+acrate2/dble(ncluster)

         do i=1,q
            betar(i)=0.d0
         end do

         do i=1,nsubject
            do j=1,q
               b(i,j)=bclus(ss(i),j)
               betar(j)=betar(j)+b(i,j)
            end do
         end do

         do i=1,q
            betar(i)=betar(i)/dble(nsubject)
         end do

c         call dblepr("betar",-1,betar,q)  


c++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,q
            workvb1(i)=smu(i)
            do j=1,q
               workb1(i,j)=(sigmainv(i,j)*dble(ncluster))+psiinv(i,j)
            end do
         end do

         call invdet(workb1,q,workb2,detlog,iflagb,workvb2)

         do i=1,ncluster
            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+sigmainv(j,k)*bclus(i,k)
               end do
               workvb1(j)=workvb1(j)+tmp1
            end do
         end do
     
         do i=1,q
            workvb2(i)=0.d0
         end do
     
         do i=1,q
            tmp1=0.d0
            do j=1,q
               tmp1=tmp1+workb2(i,j)*workvb1(j)
            end do
            workvb2(i)=tmp1
         end do
          
         call rmvnorm(q,workvb2,workb2,workmh2,workvb1,theta)


c+++++++ check if the user has requested an interrupt
         call rchkusr()
     
         do i=1,q
            mu(i)=theta(i)
            do j=1,q
               quadf(i,j)=0.d0
            end do
         end do

c         call dblepr("mu",-1,mu,q)  


         do i=1,ncluster
            do j=1,q
               do k=1,q
                  quadf(j,k)=quadf(j,k)+               
     &                       (bclus(i,j)-mu(j))*(bclus(i,k)-mu(k))                   
               end do
            end do
         end do

         do i=1,q
            do j=1,q
               quadf(i,j)=quadf(i,j)+tinv(i,j)
            end do
         end do

         call riwishart(q,nu0+ncluster,quadf,workb1,workb2,workvb1,
     &                  workmh2,workmh3,iflagb)

         do i=1,q
            do j=1,q
               sigma(i,j)=quadf(i,j)
               sigmainv(i,j)=workb1(i,j)
            end do
         end do

c         call dblepr("sigma",-1,sigma,q*q)  

c++++++++++++++++++++++++++++++++++         
c+++++++ Precision parameter
c++++++++++++++++++++++++++++++++++
         if(aa0.gt.0.d0)then
            call samalph(alpha,aa0,ab0,ncluster,nsubject)
         end if 

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ regression coefficient information

               do i=1,q
                  thetasave(isave,i)=betar(i)
               end do

               if(nfixed.gt.0)then
                  do i=1,p
                     thetasave(isave,q+i)=beta(i)
                  end do
               end if   

c+++++++++++++ baseline mean

               do i=1,q
                  thetasave(isave,q+nfixed+i)=mu(i)
               end do

c+++++++++++++ baseline covariance

               k=0
               do i=1,q
                  do j=i,q
                     k=k+1
                     thetasave(isave,q+nfixed+q+k)=sigma(i,j)
                  end do
               end do

c+++++++++++++ cluster information
               k=(q*(q+1)/2)  
               thetasave(isave,q+nfixed+q+k+1)=ncluster
               thetasave(isave,q+nfixed+q+k+2)=alpha


c+++++++++++++ random effects

               k=0
               do i=1,nsubject
                  do j=1,q
                     k=k+1
                     randsave(isave,k)=b(i,j)
                  end do   
               end do


c+++++++++++++ predictive information

               do i=1,ncluster
                  prob(i)=dble(ccluster(i))/(alpha+dble(nsubject))
               end do
               prob(ncluster+1)=alpha/(alpha+dble(nsubject))

               call simdisc(prob,nsubject+2,ncluster+1,evali)
               
               if(evali.le.ncluster)then
                  do j=1,q
                     theta(j)=bclus(evali,j)
                  end do
               end if
               if(evali.gt.ncluster)then
                  call rmvnorm(q,mu,sigma,workmh2,workvb2,theta)
               end if
               
               do i=1,q
                  k=k+1
                  randsave(isave,k)=theta(i) 
               end do

c+++++++++++++ cpo

               do i=1,nrec
                  yij=y(i)
                  
                  eta=0.d0

                  do j=1,p
                     eta=eta+x(i,j)*beta(j)
                  end do
		  do j=1,q
		     eta=eta+z(i,j)*b(subject(i),j)
                  end do
                  
                  eta=eta+roffset(i)

                  tmp1=dble(yij)*eta-exp(eta)-dgamlog(dble(yij+1))

                  cpo(i)=cpo(i)+1.0d0/exp(tmp1)  
               end do

c+++++++++++++ print
               skipcount = 0
               if(dispcount.ge.ndisplay)then
c                  call intpr("isave",5,isave,1)
                  tmp1=sprint(isave,nsave)
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
      subroutine relabelg(ind,since,nsubject,q,ncluster,ccluster,ss,
     &                   bclus,theta)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2005
      implicit none
      integer i,j,ind,since,nsubject,q,ncluster,ccluster(nsubject)
      integer ss(nsubject)
      real*8 bclus(nsubject,q),theta(q)

      do i=1,q
         theta(i)=bclus(since,i)
      end do
      
      do i=since+1,ncluster
         do j=1,nsubject
            if(ss(j).eq.i)then
               ss(j)=i-1 
            end if
         end do
         do j=1,q
            bclus(i-1,j)=bclus(i,j)
         end do
         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      
      do i=1,q
         bclus(ncluster,i)=theta(i)
      end do
      ccluster(ncluster)=1
      
      return
      end  
      