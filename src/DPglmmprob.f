
c=======================================================================                      
      subroutine spproblme(datastr,maxni,nrec,nsubject,nfixed,p,q,
     &             subject,x,xtx,y,yr,z,a0b0,nu0,prec,psiinv,sb,smu,
     &             tinv,mcmc,nsave,randsave,thetasave,cpo,alpha,b,bclus,
     &             beta,betar,mu,ncluster,sigma,ss,ccluster,
     &             iflag,iflag2,iflagb,prob,quadf,res,seed,
     &             sigmainv,theta,work1,work2,workb1,workb2,
     &             workmh1,workmh2,workmh3,workk1,workkv1,workkm1,
     &             workkm2,workv1,workv2,workvb1,workvb2,xty,ywork,zty,
     &             ztz,ztzinv)
c=======================================================================                      
c
c     Version 1.0: 
c     Last modification: 28-04-2006.
c
c     Subroutine `spproblme' to run a Markov chain in the semiparametric 
c     probit mixed model. In this routine, inference is based on the 
c     Polya urn representation of Dirichlet process.
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
c        subject     :  integer vector giving the subject for each.
c                       observation, subject(nsubject).
c        x           :  real matrix giving the design matrix for the 
c                       fixed effects, x(nrec,p). 
c        xtx         :  real matrix givind the product X^tX, xtx(p,p).
c        y           :  real vector giving the response variable,
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
c        cpo         :  real giving the cpo. 
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
c        iflag2      :  integer vector used to invert the of the 
c                       covariance matrix for each subject,
c                       iflag2(maxni).
c        iflagb      :  integer vector used to invert the of the lhs
c                       least square solution for the random effects,
c                       iflagb(q).
c        isave       :  index. 
c        iscan       :  index. 
c        j           :  index. 
c        k           :  index. 
c        l           :  index.
c        prob        :  real vector used to update the cluster 
c                       structure, prob(nsubject+1).
c        quadf       :  real matrix used to save the bilinear product
c                       of random effects, quadf(q,q).
c        ni          :  integer indicator used in updating the state. 
c        ns          :  integer indicator used in updating the state. 
c        nscan       :  integer indicating the total number of MCMC
c                       scans.
c        res         :  real vector used to save the residual effects,
c                       res(nrec).
c        seed1       :  seed for random number generation.
c        seed2       :  seed for random number generation.
c        seed3       :  seed for random number generation.
c        since       :  index.
c        skipcount   :  index. 
c        theta       :  real vector used to save randomnly generated
c                       random effects, theta(q).
c        tmp1        :  real used to accumulate quantities. 
c        tmp2        :  real used to accumulate quantities.
c        tmp3        :  real used to accumulate quantities.
c        tpi         :  real parameter used to evaluate the normal 
c                       density.
c        work1       :  real matrix used to update the fixed effects,
c                       work1(p,p).
c        work2       :  real matrix used to update the fixed effects,
c                       work2(p,p).
c        workb1      :  real matrix used to update the random effects,
c                       workb1(q,q).
c        workb2      :  real matrix used to update the random effects,
c                       workb2(q,q).
c        workmh1     :  real vector used to update the fixed effects,
c                       workmh1(p*(p+1)/2)
c        workmh2     :  real vector used to update the random effects,
c                       workmh2(q*(q+1)/2)
c        workmh3     :  real vector used to update the random effects,
c                       workmh3(q*(q+1)/2)
c        workk1      :  real matrix used to update the cluster 
c                       structure, workk1(maxni,q)
c        workkv1     :  real vector used to update the cluster 
c                       structure, workkv1(maxni)
c        workkm1     :  real matrix used to update the cluster 
c                       structure, workkm1(maxni,maxni).
c        workkm2     :  real matrix used to update the cluster 
c                       structure, workkm2(maxni,maxni).
c        workv1      :  real vector used to update the fixed effects,
c                       workv1(p)
c        workv2      :  real vector used to update the fixed effects,
c                       workv2(p)
c        workvb1     :  real vector used to update the random effects,
c                       workvb1(p)
c        workvb2     :  real vector used to update the random effects,
c                       workvb2(p)
c        xty         :  real vector used to save the product 
c                       Xt(Y-Zb), xty(p).
c        ywork       :  real vector used to save the measurement for
c                       a particular subject, ywork(maxni).
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
      integer datastr(nsubject,maxni+1),yr(nrec)
      real*8 y(nrec),x(nrec,p),z(nrec,q),xtx(p,p)	
      
c+++++Prior 
      integer nu0
      real*8 aa0,ab0,a0b0(2),prec(p,p),psiinv(q,q)
      real*8 sb(p),smu(q)
      real*8 tinv(q,q)      

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c+++++Output
      real*8 cpo(nrec)
      real*8 randsave(nsave,q*(nsubject+1))
      real*8 thetasave(nsave,q+nfixed+q+(q*(q+1)/2)+2)

c+++++Current values of the parameters
      integer ncluster,ss(nsubject)
      real*8 alpha,beta(p),b(nsubject,q)
      real*8 betar(q),bclus(nsubject,q)
      real*8 mu(q),sigma(q,q),sigmainv(q,q)

c+++++Working space
      integer ccluster(nsubject),evali,i,ii,iscan,isave,j,k,l,ni,ns
      integer iflag(p),iflag2(maxni),iflagb(q)
      integer nscan
      integer since,sprint
      integer seed(3),seed1,seed2,seed3,skipcount,dispcount
      real*8 cdfnorm,detlog
      real*8 prob(nsubject+1)
      real*8 quadf(q,q)
      real*8 res(nrec),rtnorm,sigma2e
      real*8 theta(q),tmp1,tmp2,tmp3,tpi
      real*8 work1(p,p),work2(p,p)
      real*8 workb1(q,q),workb2(q,q)
      real*8 workmh1(p*(p+1)/2),workmh2(q*(q+1)/2),workmh3(q*(q+1)/2)
      real*8 workk1(maxni,q)
      real*8 workkm1(maxni,maxni),workkm2(maxni,maxni)
      real*8 workkv1(maxni),workv1(p),workv2(p),workvb1(q),workvb2(q)
      real*8 xty(p)
      real*8 ywork(maxni)
      real*8 zty(q),ztz(q,q),ztzinv(q,q)
      logical ainf,asup
      parameter(tpi=6.283185307179586476925286766559d0)

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

      aa0=a0b0(1)
      ab0=a0b0(2)
      
      sigma2e=1.d0
      
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
               tmp1=tmp1+z(i,j)*b(subject(i),j) 
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


c++++++++++++++++++++++++++++++++++
c+++++++ fixed effects
c++++++++++++++++++++++++++++++++++

         if(nfixed.eq.0)go to 1
            do i=1,p
               xty(i)=sb(i)
               workv1(i)=0.d0
               workv2(i)=0.d0
            end do

            do i=1,nrec
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+z(i,j)*b(subject(i),j) 
               end do
               tmp1=y(i)-tmp1
             
               do j=1,p
                  xty(j)=xty(j)+x(i,j)*(tmp1/sigma2e)
               end do
            end do

            do i=1,p
               do j=1,p
                  work1(i,j)=xtx(i,j)/sigma2e+prec(i,j)          
               end do
            end do

            call invdet(work1,p,work2,detlog,iflag,workv1)

            do i=1,p
               tmp1=0.d0
               do j=1,p
                  tmp1=tmp1+work2(i,j)*xty(j) 
               end do
               workv2(i)=tmp1
            end do

            call rmvnorm(p,workv2,work2,workmh1,workv1,beta)
1        continue            

         
c++++++++++++++++++++++++++++++++++         
c+++++++ random effects 
c++++++++++++++++++++++++++++++++++

c++++++++++++++++++++++++++++++
c+++++++ a) Polya Urn 
c++++++++++++++++++++++++++++++

         if(nfixed.eq.0)then
             do i=1,nrec
                res(i)=y(i) 
             end do
           else
             do i=1,nrec
                tmp1=0.d0
                do j=1,p
                   tmp1=tmp1+x(i,j)*beta(j)    
                end do
                res(i)=y(i)-tmp1
             end do
         end if  


         do i=1,nsubject
         
            ns=ccluster(ss(i))
            ni=datastr(i,1) 


c++++++++++ subject in cluster with more than 1 observations
             
            if(ns.gt.1)then
          
               ccluster(ss(i))=ccluster(ss(i))-1 

               do j=1,ncluster
                  tmp3=0.d0
                  do k=1,ni
                     tmp1=0.d0
                     do l=1,q
                        tmp1=tmp1+z(datastr(i,k+1),l)*bclus(j,l)
                     end do
                     tmp1=res(datastr(i,k+1))-tmp1
                     
                     tmp3=tmp3+tmp1*tmp1
                  end do
                  tmp3=tmp3/sigma2e

                  tmp1=-(dble(ni)*log(tpi))

                  tmp2=dble(ni)*log(sigma2e)

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        (tmp1-tmp2-tmp3)/2.d0)

               end do
               
               do j=1,ni
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+z(datastr(i,j+1),k)*mu(k)
                  end do
                  ywork(j)=res(datastr(i,j+1))-tmp1
               end do
               
               do j=1,ni
                  do k=1,q
                     tmp1=0.d0
                     do l=1,q
                        tmp1=tmp1+z(datastr(i,j+1),l)*sigma(l,k)   
                     end do
                     workk1(j,k)=tmp1
                  end do
               end do
        
               do j=1,ni
                  do k=1,ni
                     tmp1=0.d0
                     do l=1,q
                       tmp1=tmp1+workk1(j,l)*z(datastr(i,k+1),l)
                     end do
                     workkm1(j,k)=tmp1
                  end do
               end do


               do j=1,ni
                  workkm1(j,j)=workkm1(j,j)+sigma2e
               end do

               call invdet2(workkm1,maxni,ni,workkm2,detlog,iflag2,
     &                      workkv1)

               tmp1=-(dble(ni)*log(tpi))
               tmp2=detlog

               tmp3=0.d0
               do j=1,ni
                  do k=1,ni
                     tmp3=tmp3+ywork(j)*workkm2(j,k)*ywork(k)                  
                  end do
               end do
               
               prob(ncluster+1)=exp(log(alpha)+(tmp1-tmp2-tmp3)/2.d0)

               call simdisc(prob,nsubject+1,ncluster+1,evali)

               ss(i)=evali
               
               ccluster(evali)=ccluster(evali)+1
               
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  
                  do j=1,q
                     do k=1,q
                        tmp1=0.d0
                        do l=1,ni
                           tmp1=tmp1+z(datastr(i,l+1),j)*
     &                          z(datastr(i,l+1),k)                           
                        end do
                        ztz(j,k)=tmp1
                     end do
                  end do
                  
                  do j=1,q
                     do k=1,q
                        ztz(j,k)=ztz(j,k)/sigma2e+sigmainv(j,k)
                     end do
                  end do   
                  
                  call invdet(ztz,q,ztzinv,detlog,iflagb,workvb2)
                  
                  do j=1,q
                     tmp1=0.d0
                     do k=1,ni
                        tmp1=tmp1+z(datastr(i,k+1),j)*
     &                        res(datastr(i,k+1))
                     end do
                     zty(j)=tmp1
                  end do
                 
                  do j=1,q
                     tmp1=0.d0
                     do k=1,q
                        tmp1=tmp1+sigmainv(j,k)*mu(k)   
                     end do
                     zty(j)=zty(j)/sigma2e+tmp1
                  end do
                 
                  do j=1,q
                     tmp1=0.d0
                     do k=1,q
                        tmp1=tmp1+ztzinv(j,k)*zty(k) 
                     end do
                     workvb1(j)=tmp1
                  end do

                  call rmvnorm(q,workvb1,ztzinv,workmh2,workvb2,theta)

                  do j=1,q
                     bclus(evali,j)=theta(j)
                  end do
                  
               end if
            end if


c++++++++++ subject in cluster with only 1 observation
             
            if(ns.eq.1)then
                
               since=ss(i)
                
               if(since.lt.ncluster)then
                   call relabel(i,since,nsubject,q,ncluster,
     &                          ccluster,ss,bclus,theta)                   
	       end if
	              
               ccluster(ncluster)=ccluster(ncluster)-1 
               ncluster=ncluster-1

               do j=1,ncluster
                  tmp3=0.d0
                  do k=1,ni
                     tmp1=0.d0
                     do l=1,q
                        tmp1=tmp1+z(datastr(i,k+1),l)*bclus(j,l)
                     end do
                     tmp1=res(datastr(i,k+1))-tmp1                     
                     tmp3=tmp3+tmp1*tmp1
                  end do
                  tmp3=tmp3/sigma2e

                  tmp1=-(dble(ni)*log(tpi))

                  tmp2=dble(ni)*log(sigma2e)

                  prob(j)=exp(log(dble(ccluster(j)))+
     &                        (tmp1-tmp2-tmp3)/2.d0)
     
               end do

               do j=1,ni
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+z(datastr(i,j+1),k)*mu(k)
                  end do
                  ywork(j)=res(datastr(i,j+1))-tmp1
               end do

               do j=1,ni
                  do k=1,q
                     tmp1=0.d0
                     do l=1,q
                        tmp1=tmp1+z(datastr(i,j+1),l)*sigma(l,k)   
                     end do
                     workk1(j,k)=tmp1
                  end do
               end do
        
               do j=1,ni
                  do k=1,ni
                     tmp1=0.d0
                     do l=1,q
                       tmp1=tmp1+workk1(j,l)*z(datastr(i,k+1),l)
                     end do
                     workkm1(j,k)=tmp1
                  end do
               end do

               do j=1,ni
                  workkm1(j,j)=workkm1(j,j)+sigma2e
               end do

               call invdet2(workkm1,maxni,ni,workkm2,detlog,iflag2,
     &                      workkv1)

               tmp1=-(dble(ni)*log(tpi))
               tmp2=detlog

               tmp3=0.d0
               do j=1,ni
                  do k=1,ni
                     tmp3=tmp3+ywork(j)*workkm2(j,k)*ywork(k)                  
                  end do
               end do
               
               prob(ncluster+1)=exp(log(alpha)+(tmp1-tmp2-tmp3)/2.d0)


               call simdisc(prob,nsubject+1,ncluster+1,evali)
               
               ss(i)=evali
               
               ccluster(evali)=ccluster(evali)+1
               
               if(evali.gt.ncluster)then
                  ncluster=ncluster+1
                  
                  do j=1,q
                     do k=1,q
                        tmp1=0.d0
                        do l=1,ni
                           tmp1=tmp1+z(datastr(i,l+1),j)*
     &                          z(datastr(i,l+1),k)                           
                        end do
                        ztz(j,k)=tmp1
                     end do
                  end do

                  do j=1,q
                     do k=1,q
                        ztz(j,k)=ztz(j,k)/sigma2e+sigmainv(j,k)
                     end do
                  end do   
 
                  call invdet(ztz,q,ztzinv,detlog,iflagb,workvb2)
                  
                  do j=1,q
                     tmp1=0.d0
                     do k=1,ni
                        tmp1=tmp1+z(datastr(i,k+1),j)*
     &                        res(datastr(i,k+1))
                     end do
                     zty(j)=tmp1
                  end do
                  
                  do j=1,q
                     tmp1=0.d0
                     do k=1,q
                        tmp1=tmp1+sigmainv(j,k)*mu(k)   
                     end do
                     zty(j)=zty(j)/sigma2e+tmp1
                  end do
                 
                  
                  do j=1,q
                     tmp1=0.d0
                     do k=1,q
                        tmp1=tmp1+ztzinv(j,k)*zty(k) 
                     end do
                     workvb1(j)=tmp1
                  end do

                  call rmvnorm(q,workvb1,ztzinv,workmh2,workvb2,theta)

                  do j=1,q
                     bclus(evali,j)=theta(j)
                  end do
                  
               end if
            end if


         end do

c         call intpr("ncluster",-1,ncluster,1)
c         do i=1,ncluster
c            call dblepr("blcus",-1,bclus(i,1),1)
c         end do
c         return 
c         call intpr("ss",-1,ss,nsubject)
c         call dblepr("mu",-1,mu,q)
c         call dblepr("betar",-1,betar,q)
         

c++++++++++++++++++++++++++++++
c+++++++ b) Resampling step
c++++++++++++++++++++++++++++++

         do ii=1,ncluster

c++++++++++ check if the user has requested an interrupt
            call rchkusr()
    
            do i=1,q
               do j=1,q
                  ztz(i,j)=0.d0
               end do
               zty(i)=0.d0
            end do

            do i=1,nsubject
               if(ii.eq.ss(i))then
                  ni=datastr(i,1) 
                  do j=1,q
                     do k=1,q
                        tmp1=0.d0
                        do l=1,ni
                           tmp1=tmp1+z(datastr(i,l+1),j)*
     &                          z(datastr(i,l+1),k)                           
                        end do
                        ztz(j,k)=ztz(j,k)+tmp1
                     end do
                  end do
                  do j=1,q
                     tmp1=0.d0
                     do k=1,ni
                        tmp1=tmp1+z(datastr(i,k+1),j)*
     &                        res(datastr(i,k+1))
                     end do
                     zty(j)=zty(j)+tmp1
                  end do                  
               end if
            end do
            
            do i=1,q
               do j=1,q
                  ztz(i,j)=ztz(i,j)/sigma2e+sigmainv(i,j)
               end do
            end do

            call invdet(ztz,q,ztzinv,detlog,iflagb,workvb2)
            
            do i=1,q
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+sigmainv(i,j)*mu(j)   
               end do
               zty(i)=zty(i)/sigma2e+tmp1
            end do
            
            do i=1,q
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+ztzinv(i,j)*zty(j) 
               end do
               workvb1(i)=tmp1
            end do
            
            call rmvnorm(q,workvb1,ztzinv,workmh2,workvb2,theta)

            do i=1,q
               bclus(ii,i)=theta(i)
            end do

         end do


c++++++++++++++++++++++++++++++++++         
c+++++++ update b
c++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

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

c         call dblepr("sigmainv",-1,tinv,q*q)


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

               call simdisc(prob,nsubject+1,ncluster+1,evali)
               
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
                 
                  tmp1=0.d0

                  do j=1,p
                     tmp1=tmp1+x(i,j)*beta(j)
                  end do
		  do j=1,q
		     tmp1=tmp1+z(i,j)*b(subject(i),j)
                  end do
                  
                  tmp1=cdfnorm(tmp1,0.d0,1.d0,1,0)

                  tmp1=dble(yr(i))  *log(tmp1)+
     &                 dble(1-yr(i))*log(1.d0-tmp1)
               
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
      
      do i=1,nrec
         cpo(i)=dble(nsave)/cpo(i)
      end do
            
      
      return
      end
         
