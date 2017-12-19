c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR RANDOM SEQUENTIAL BARYCENTER ARRAYS
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
c=======================================================================                  

c=======================================================================                  
c=======================================================================                  
c=======================================================================                  
c     AUXILIARY FUNCTIONS/COMMON FUNCTIONS
c=======================================================================                  
c=======================================================================                  
c=======================================================================                  



c=======================================================================
       subroutine gennsba(maxn,maxr,mu0,sd0,mu1,sd1,sba)
c=======================================================================
c      Generates a 'normal' sequential barycenter array.
c
c      Alejandro Jara, 2010. 
c=======================================================================
       implicit none

c+++++ Input
       integer maxn,maxr
       double precision mu0,sd0,mu1,sd1
 
c+++++ Output
       double precision sba(maxn,maxr)

c+++++ Internal working space
       integer i,j,nn
       double precision rnorm,rtnorm,ll,uu
       logical linf,lsup

c+++++ Algorithm

       if(sd0.eq.0.d0)then
          sba(1,1)=mu0
         else
          sba(1,1)=rnorm(mu0,sd0)   
       end if

       do i=2,maxn

          call rchkusr()

          nn=(2**(i-1))-1
          do j=1,nn
             sba(i,2*j)=sba(i-1,j)
          end do
   
          nn=2**(i-1)
          do j=1,nn
             if(j.eq.1)then
                linf=.true.
                lsup=.false.
                ll=0.d0
                uu=sba(i-1,j)
                sba(i,2*j-1)=rtnorm(mu1,sd1,ll,uu,linf,lsup)
              else if(j.eq.nn)then
                linf=.false.
                lsup=.true.
                ll=sba(i-1,j-1)
                uu=0.d0
                sba(i,2*j-1)=rtnorm(mu1,sd1,ll,uu,linf,lsup)
              else
                linf=.false.
                lsup=.false.
                ll=sba(i-1,j-1)
                uu=sba(i-1,j)
                sba(i,2*j-1)=rtnorm(mu1,sd1,ll,uu,linf,lsup)
             end if   
          end do
       end do 

       return
       end
       

c=======================================================================
       subroutine discrsba(maxn,maxr,maxc,sba,f,theta,w)
c=======================================================================
c      Computes the weights and point masses from a RSBA.
c
c      Alejandro Jara, 2010. 
c=======================================================================
       implicit none

c+++++ Input
       integer maxn,maxr,maxc
       double precision sba(maxn,maxr)
 
c+++++ External working space
       double precision f(maxn-1,maxc-1)

c+++++ Output
       double precision theta(maxc)
       double precision w(maxc)

c+++++ Internal working space
       integer count,i,j,nn
       double precision ll,uu

c+++++ Algorithm

       f(1,1)=(sba(2,3)-sba(2,2))/(sba(2,3)-sba(2,1))

       do i=2,maxn-1

          call rchkusr()

          nn=(2**(i-1))-1
          do j=1,nn
             f(i,2*j)=f(i-1,j)
          end do
   
          nn=2**(i-1)
          do j=1,nn
             if(j.eq.1)then
                ll=0.0
                uu=f(i-1,j)
              else if(j.eq.nn)then
                ll=f(i-1,j-1)
                uu=1.0
              else
                ll=f(i-1,j-1)
                uu=f(i-1,j)
             end if
             f(i,2*j-1)=ll+(uu-ll)*
     &                     exp(
     &                       log(sba(i+1,4*j-1)-sba(i+1,4*j-2))-
     &                       log(sba(i+1,4*j-1)-sba(i+1,4*j-3)))  
          end do
       end do 
 
       w(1)=f(maxn-1,1)
       w(maxc)=1.d0-f(maxn-1,maxc-1)

       do i=2,maxc-1
          w(i)=f(maxn-1,i)-f(maxn-1,i-1)
       end do

       count=1
       do i=1,maxc
          theta(i)=sba(maxn,count)
          count=count+2
       end do 

       return
       end

c=======================================================================      
      integer function ceil(x)
c=======================================================================      
c     Ceiling function.
c 
c     Alejandro Jara, 2010. 
c=======================================================================      
      implicit none
      double precision x
      ceil=x
      if(ceil.lt.x) ceil=ceil+1
      end  
       
c=======================================================================      
      subroutine isweight(x,out)
c=======================================================================      
c     Determines whether a barycenter just modifies the weights 
c     (out=1) or the atoms as well (out=0).
c
c     Alejandro Jara, 2010. 
c=======================================================================      
      implicit none
      integer x,y,ceil,out

      y=ceil(x/2.d0)*2
      out=0
      if(y.eq.x) out=1
      return
      end  


c=======================================================================      
      subroutine rsbauclus(nrec,y,maxn,maxc,w,theta,sigma2,prob,
     &                     ccluster,cstrt,ss)
c=======================================================================      
c      Update the cluster indicators.
c
c      Alejandro Jara, 2010. 
c=======================================================================      
      implicit none

c++++ Input
      integer maxn,maxc,nrec
      double precision y(nrec) 
      double precision w(maxc)
      double precision theta(maxc)
      double precision sigma2(maxc)

c++++ External working space
      double precision prob(maxc)

c++++ Output
      integer ccluster(maxc)
      integer cstrt(maxc,nrec)
      integer ss(nrec)

c++++ Internal working space
      integer i,ind,j
      double precision dnrm

c++++ Algorithm

      do i=1,maxc
         ccluster(i)=0
      end do

      do i=1,nrec
         call rchkusr()

         do j=1,maxc
            prob(j)=w(j)*dnrm(y(i),theta(j),sqrt(sigma2(j)),0)
         end do
         call simdisc(prob,maxc,maxc,ind)

         ss(i)=ind
         ccluster(ss(i))=ccluster(ss(i))+1     
         cstrt(ss(i),ccluster(ss(i)))=i
      end do

      return
      end  


c=======================================================================      
      subroutine rsbauclus2(nrec,y,maxn,maxc,w,theta,sigma2,prob,
     &                      ccluster,cstrt,ss)
c=======================================================================      
c      Update the cluster indicators.
c
c      Alejandro Jara, 2010. 
c=======================================================================      
      implicit none

c++++ Input
      integer maxn,maxc,nrec
      double precision y(nrec) 
      double precision w(maxc)
      double precision theta(maxc)
      double precision sigma2

c++++ External working space
      double precision prob(maxc)

c++++ Output
      integer ccluster(maxc)
      integer cstrt(maxc,nrec)
      integer ss(nrec)

c++++ Internal working space
      integer i,ind,j
      double precision dnrm

c++++ Algorithm

      do i=1,maxc
         ccluster(i)=0
      end do

      do i=1,nrec
         call rchkusr()

         do j=1,maxc
            prob(j)=w(j)*dnrm(y(i),theta(j),sqrt(sigma2),0)
         end do
         call simdisc(prob,maxc,maxc,ind)

         ss(i)=ind
         ccluster(ss(i))=ccluster(ss(i))+1     
         cstrt(ss(i),ccluster(ss(i)))=i
      end do

      return
      end  


c=======================================================================      
      subroutine rsbaukvar(nrec,y,maxn,maxc,ccluster,cstrt,ss,theta,
     &                     tau1,tau2,
     &                     sigma2)
c=======================================================================      
c     Updates the kernel variances.
c 
c     Alejandro Jara, 2010. 
c=======================================================================      
      implicit none

c++++ Input
      integer maxn,maxc,nrec
      integer ccluster(maxc)
      integer cstrt(maxc,nrec)
      integer ss(nrec)
      double precision y(nrec) 
      double precision theta(maxc)
      double precision tau1,tau2

c++++ Output
      double precision sigma2(maxc)

c++++ Internal working space
      integer i,j,ns
      double precision rgamma,sse,tmp1

c++++ Algorithm

      do i=1,maxc
         call rchkusr()

         ns=ccluster(i)
         sse=0.d0
 
         if(ns.gt.0)then
            do j=1,ns
               tmp1=y(cstrt(i,j))-theta(i)
               sse=sse+tmp1*tmp1
            end do
         end if
         sigma2(i)=1.d0/rgamma(0.5d0*(tau1+dble(ns)),
     &                         0.5d0*(tau2+sse))
      end do

      return
      end  


c=======================================================================      
      subroutine rsbaukvar2(nrec,y,maxn,maxc,ccluster,cstrt,ss,theta,
     &                      tau1,tau2,
     &                      sigma2)
c=======================================================================      
c     Updates the kernel variance.
c 
c     Alejandro Jara, 2010. 
c=======================================================================      
      implicit none

c++++ Input
      integer maxn,maxc,nrec
      integer ccluster(maxc)
      integer cstrt(maxc,nrec)
      integer ss(nrec)
      double precision y(nrec) 
      double precision theta(maxc)
      double precision tau1,tau2

c++++ Output
      double precision sigma2

c++++ Internal working space
      integer i,j,ns
      double precision rgamma,sse,tmp1

c++++ Algorithm

      do i=1,maxc
         call rchkusr()

         ns=ccluster(i)
         sse=0.d0
 
         if(ns.gt.0)then
            do j=1,ns
               tmp1=y(cstrt(i,j))-theta(i)
               sse=sse+tmp1*tmp1
            end do
         end if
      end do

      sigma2=1.d0/rgamma(0.5d0*(tau1+dble(nrec)),
     &                   0.5d0*(tau2+sse))

      return
      end  

c=======================================================================      
      subroutine rsbaugamma(maxc,sigma2,
     &                      tau1,taus1,taus2,
     &                      tau2)
c=======================================================================      
c      Update the gamma hyperparameter for the kernel variances.
c
c      Alejandro Jara, 2010. 
c=======================================================================      
      implicit none

c++++ Input
      integer maxc
      double precision sigma2(maxc)
      double precision tau1,taus1,taus2

c++++ Output
      double precision tau2

c++++ Internal working space
      integer i
      double precision tmp1,rgamma

c++++ Algorithm

      tmp1=0.d0
      do i=1,maxc
         call rchkusr()
         tmp1=tmp1+1.d0/sigma2(i)
      end do 

      tau2=rgamma(0.5d0*(dble(maxc)*tau1+taus1),
     &            0.5d0*(tmp1+taus2))   

      return
      end

c=======================================================================      
      subroutine rsbaposs(maxn,ii,ismean,nb,possv)
c=======================================================================      
c     Determines the elements that need to be changed when
c     updating the elements of the last level of the barycenter array.
c     The output includes the vector of positions (possv), 
c     how many levels need to be updated (nb) and whether
c     the barycenter corresponds to the mean (ismean=1) or not
c     (ismean=0).
c 
c     Alejandro Jara, 2010.
c=======================================================================      

      implicit none
c++++ Input
      integer maxn,ii

c++++ Output
      integer ismean,nb,possv(maxn)

c++++ Internal working space
      integer i,ind,jj,out

c++++ Algorithm
   
      out=0
      nb=0
      ismean=0
      do i=1,maxn
         possv(i)=0
      end do

      possv(maxn)=ii
      i=0
      do while(out.eq.0.and.i.lt.(maxn-1))
         i=i+1 
         jj=possv(maxn-i+1)
         call isweight(jj,ind)
         if(ind.eq.0)then
            out=1
           else
            nb=nb+1
            possv(maxn-i)=int(jj/2)
         end if
      end do

      if(nb.eq.(maxn-1))ismean=1
      return
      end


c=======================================================================      
      subroutine rsbawlogpost(ii,cind,xx,sigma2c,isw,ism,
     &                        nrec,y,
     &                        maxn,maxc,w,
     &                        ccluster,cstrt,
     &                        mu0,sigma20,mu1,sigma21,
     &                        out)
c=======================================================================      
      implicit none
c++++ Input
      integer cind,ii,maxn,maxc,nrec,isw,ism
      integer ccluster(maxc)
      integer cstrt(maxc,nrec)
      double precision xx,sigma2c,w(maxc)
      double precision mu0,sigma20,mu1,sigma21
      double precision y(nrec)

c++++ Output
      double precision out
 
c++++ Internal working space
      integer i,ns
      double precision dnrm

c++++ Algorithm

      out=0.d0
      do i=1,maxc
         out=out+dble(ccluster(i))*log(w(i))
      end do

      if(ism.eq.1)then
         out=out+dnrm(xx,mu0,sqrt(sigma20),1) 
       else
         out=out+dnrm(xx,mu1,sqrt(sigma21),1) 
      end if  

      if(isw.eq.0)then
         ns=ccluster(cind)

         if(ns.gt.0)then
            do i=1,ns
               out=out+
     &         dnrm(y(cstrt(cind,i)),xx,sqrt(sigma2c),1)
            end do
         end if
      end if

      return
      end


c=======================================================================      
      subroutine rsbauba(nrec,y,
     &                   maxn,maxr,maxc,
     &                   mu0,mu1,sigma20,sigma21,
     &                   ccluster,cstrt,sigma2,
     &                   f,possv,
     &                   sba,theta,w)
c=======================================================================      
c      Update the RSBA parameters.
c
c      Alejandro Jara, 2010. 
c=======================================================================      

      implicit none

c++++ Data
      integer nrec
      double precision y(nrec)

c++++ Prior
      integer maxn,maxr,maxc
      double precision mu0,sigma20
      double precision mu1,sigma21

c++++ Current state of the chain (cluster info)
      integer ccluster(maxc)
      integer cstrt(maxc,nrec)
      double precision sigma2(maxc)

c++++ External working space
      integer possv(maxn)
      double precision f(maxn-1,maxc-1)

c++++ Input/Output (RSBA elements)
      double precision sba(maxn,maxr)
      double precision theta(maxc)
      double precision w(maxc)

c++++ Internal working space
      integer cind
      integer i,isw,ism,j,k,nb
      double precision tmp1

c++++ RNG and distributions
      real runif

c++++ Working space slice sampling
      integer evali
      double precision rexpo,re,uwork
      double precision logy,xx0,xx1,llim,rlim
      double precision grlim,gllim,gxx1

c++++ Algorithm

      cind=0
            
      do i=1,maxr

         call rchkusr()

         call isweight(i,isw)
         call rsbaposs(maxn,i,ism,nb,possv)

         if(isw.eq.0)cind=cind+1

         if(ism.eq.0.or.(ism.eq.1.and.sigma20.gt.0.d0))then
            evali=1
            xx0=sba(maxn,i)
            re=rexpo(1.d0)
            logy=-re

            call rsbawlogpost(i,cind,xx0,sigma2(cind),isw,ism,
     &                        nrec,y,
     &                        maxn,maxc,w,
     &                        ccluster,cstrt,
     &                        mu0,sigma20,mu1,sigma21,
     &                        tmp1)

            logy=logy+tmp1
            uwork=dble(runif())*0.25d0  
            llim=xx0-uwork
            rlim=llim+0.25d0

            if(i.eq.1)then
               if(rlim.gt.sba(maxn,i+1))rlim=sba(maxn,i+1) 
              else if(i.eq.maxr)then
               if(llim.lt.sba(maxn,i-1))llim=sba(maxn,i-1) 
              else
               if(llim.lt.sba(maxn,i-1))llim=sba(maxn,i-1) 
               if(rlim.gt.sba(maxn,i+1))rlim=sba(maxn,i+1) 
            end if


            evali=evali+1
            sba(maxn,i)=llim
            if(nb.gt.0)then
               do j=1,nb
                  k=maxn-j
                  sba(k,possv(k))=llim                  
               end do
            end if
            call discrsba(maxn,maxr,maxc,sba,f,theta,w)

            call rsbawlogpost(i,cind,llim,sigma2(cind),isw,ism,
     &                        nrec,y,
     &                        maxn,maxc,w,
     &                        ccluster,cstrt,
     &                        mu0,sigma20,mu1,sigma21,
     &                        gllim)


            evali=evali+1
            sba(maxn,i)=rlim
            if(nb.gt.0)then
               do j=1,nb
                  k=maxn-j
                  sba(k,possv(k))=rlim                  
               end do
            end if
            call discrsba(maxn,maxr,maxc,sba,f,theta,w)

            call rsbawlogpost(i,cind,rlim,sigma2(cind),isw,ism,
     &                        nrec,y,
     &                        maxn,maxc,w,
     &                        ccluster,cstrt,
     &                        mu0,sigma20,mu1,sigma21,
     &                        grlim)


            do while(gllim.gt.logy)
               llim=llim-0.25d0

               if(i.gt.1.and.llim.lt.sba(maxn,i-1))then
                  llim=sba(maxn,i-1) 
                  gllim=logy-1.d0
                else
                  evali=evali+1
                  sba(maxn,i)=llim
                  if(nb.gt.0)then
                     do j=1,nb
                        k=maxn-j
                        sba(k,possv(k))=llim                  
                     end do
                  end if
                  call discrsba(maxn,maxr,maxc,sba,f,theta,w)

                  call rsbawlogpost(i,cind,llim,sigma2(cind),
     &                              isw,ism,
     &                              nrec,y,
     &                              maxn,maxc,w,
     &                              ccluster,cstrt,
     &                              mu0,sigma20,mu1,sigma21,
     &                              gllim)
               end if
            end do 


            do while(grlim.gt.logy)
               rlim=rlim+0.25d0

               if(i.lt.maxr.and.rlim.gt.sba(maxn,i+1))then
                  rlim=sba(maxn,i+1) 
                  grlim=logy-1.d0
                else
                  evali=evali+1
                  sba(maxn,i)=rlim
                  if(nb.gt.0)then
                     do j=1,nb
                        k=maxn-j
                        sba(k,possv(k))=rlim                  
                     end do
                  end if
                  call discrsba(maxn,maxr,maxc,sba,f,theta,w)

                  call rsbawlogpost(i,cind,rlim,sigma2(cind),
     &                              isw,ism,
     &                              nrec,y,
     &                              maxn,maxc,w,
     &                              ccluster,cstrt,
     &                              mu0,sigma20,mu1,sigma21,
     &                              grlim)
               end if
            end do 


            xx1=llim+(rlim-llim)*dble(runif())
            evali=evali+1

            sba(maxn,i)=xx1
            if(nb.gt.0)then
               do j=1,nb
                  k=maxn-j
                  sba(k,possv(k))=xx1                  
               end do
            end if
            call discrsba(maxn,maxr,maxc,sba,f,theta,w)

            call rsbawlogpost(i,cind,xx1,sigma2(cind),isw,ism,
     &                        nrec,y,
     &                        maxn,maxc,w,
     &                        ccluster,cstrt,
     &                        mu0,sigma20,mu1,sigma21,
     &                        gxx1)


            do while(gxx1.lt.logy)
               if(xx1.gt.xx0)rlim=xx1
               if(xx1.lt.xx0)llim=xx1

               if(i.eq.1)then
                  if(rlim.gt.sba(maxn,i+1))rlim=sba(maxn,i+1) 
                 else if(i.eq.maxr)then
                  if(llim.lt.sba(maxn,i-1))llim=sba(maxn,i-1) 
                 else
                  if(llim.lt.sba(maxn,i-1))llim=sba(maxn,i-1) 
                  if(rlim.gt.sba(maxn,i+1))rlim=sba(maxn,i+1) 
               end if

               xx1=llim+(rlim-llim)*dble(runif())

               sba(maxn,i)=xx1
               if(nb.gt.0)then
                  do j=1,nb
                     k=maxn-j
                     sba(k,possv(k))=xx1                  
                  end do
               end if
               call discrsba(maxn,maxr,maxc,sba,f,theta,w)

               call rsbawlogpost(i,cind,xx1,sigma2(cind),isw,ism,
     &                           nrec,y,
     &                           maxn,maxc,w,
     &                           ccluster,cstrt,
     &                           mu0,sigma20,mu1,sigma21,
     &                           gxx1)

c               call intpr("evaluation c#",-1,evali,1)
c               call dblepr("xx1",-1,xx1,1)
c               call dblepr("gxx1",-1,gxx1,1)
            end do

         end if

      end do

      return
      end



c=======================================================================
       subroutine lposth1msba(maxn,maxr,mu1,sigma21,sba,m1,s1,lpost)
c=======================================================================
c      Computes the log posterior for the mean of the normal
c      measure of the conditional means in a 'normal' sequential 
c      barycenter array.
c
c      Alejandro Jara, 2011. 
c=======================================================================
       implicit none

c+++++ Input
       integer maxn,maxr
       double precision mu1,sigma21
       double precision sba(maxn,maxr)
       double precision m1,s1
 
c+++++ Output
       double precision lpost

c+++++ Internal working space
       integer i,j,nn
       double precision ll,uu,xx
       double precision tmp1,tmp2,tmp3,tmp4
       double precision cdfnorm,dnrm
       double precision sd1

c+++++ Algorithm

       sd1=dsqrt(sigma21)
       lpost=0.d0

       do i=2,maxn

          call rchkusr()

          do j=1,nn
             xx=sba(i,2*j-1)
             tmp1=dnrm(xx,mu1,sd1,0)

             if(j.eq.1)then
                uu=sba(i-1,j)
                tmp2=cdfnorm(uu,mu1,sd1,1,0)
                lpost=lpost+dlog(tmp1/tmp2)

              else if(j.eq.nn)then
                ll=sba(i-1,j-1)
                tmp2=cdfnorm(ll,mu1,sd1,1,1)
                lpost=lpost+dlog(tmp1/tmp2)
              else
                ll=sba(i-1,j-1)
                uu=sba(i-1,j)
                tmp2=cdfnorm(ll,mu1,sd1,1,0)
                tmp3=cdfnorm(uu,mu1,sd1,1,0)
                tmp4=tmp3-tmp2
                lpost=lpost+dlog(tmp1/tmp4)
             end if   
          end do
       end do 

       lpost=lpost+dnrm(mu1,m1,s1,1)

       return
       end
       

c=======================================================================
       subroutine lposth1ssba(maxn,maxr,mu1,sigma21,sba,tauh1,
     &                        tauh2,lpost)
c=======================================================================
c      Computes the log posterior for the variance of the normal
c      measure of the conditional means in a 'normal' sequential 
c      barycenter array.
c
c      Alejandro Jara, 2011. 
c=======================================================================
       implicit none

c+++++ Input
       integer maxn,maxr
       double precision mu1,sigma21
       double precision sba(maxn,maxr)
       double precision tauh1,tauh2
 
c+++++ Output
       double precision lpost

c+++++ Internal working space
       integer i,j,nn
       double precision ll,uu,xx
       double precision tmp1,tmp2,tmp3,tmp4
       double precision cdfnorm,dnrm
       double precision sd1
c+++++ Algorithm

       lpost=0.d0
       sd1=dsqrt(sigma21)

       do i=2,maxn

          call rchkusr()

          do j=1,nn
             xx=sba(i,2*j-1)
             tmp1=dnrm(xx,mu1,sd1,0)

             if(j.eq.1)then
                uu=sba(i-1,j)
                tmp2=cdfnorm(uu,mu1,sd1,1,0)
                lpost=lpost+dlog(tmp1/tmp2)

              else if(j.eq.nn)then
                ll=sba(i-1,j-1)
                tmp2=cdfnorm(ll,mu1,sd1,1,1)
                lpost=lpost+dlog(tmp1/tmp2)
              else
                ll=sba(i-1,j-1)
                uu=sba(i-1,j)
                tmp2=cdfnorm(ll,mu1,sd1,1,0)
                tmp3=cdfnorm(uu,mu1,sd1,1,0)
                tmp4=tmp3-tmp2
                lpost=lpost+dlog(tmp1/tmp4)
             end if   
          end do
       end do 

       xx=sigma21
       tmp1=-(0.5d0*tauh1+1.d0)*dlog(xx)-0.5d0*tauh2/xx
       lpost=lpost+tmp1
       return
       end



c=======================================================================
       subroutine gennsbad(nrec,y,maxn,maxr,mu0,sd0,mu1,sd1,sba)
c=======================================================================
c      Generates a 'normal' sequential barycenter array from the data.
c
c      Alejandro Jara, 2011. 
c=======================================================================
       implicit none

c+++++ Input
       integer maxn,maxr,nrec
       double precision mu0,sd0,mu1,sd1
       double precision y(nrec)
 
c+++++ Output
       double precision sba(maxn,maxr)

c+++++ Internal working space
       integer i,ii,j,nn
       integer count
       double precision rtnorm,ll,uu
       double precision tmp1
       logical linf,lsup

c+++++ Algorithm

       if(sd0.eq.0.d0)then
          sba(1,1)=mu0
         else
          tmp1=0.d0
          do ii=1,nrec
             tmp1=tmp1+y(ii)
          end do
          sba(1,1)=tmp1/dble(nrec)   
       end if

       do i=2,maxn

          call rchkusr()

          nn=(2**(i-1))-1
          do j=1,nn
             sba(i,2*j)=sba(i-1,j)
          end do
   
          nn=2**(i-1)
          do j=1,nn
             if(j.eq.1)then
                linf=.true.
                lsup=.false.
                ll=0.d0
                uu=sba(i-1,j)

                tmp1=0.d0
                count=0
                do ii=1,nrec
                   if(y(ii).le.uu)then
                      tmp1=tmp1+y(ii)
                      count=count+1
                   end if
                end do

                if(count.gt.0)then
                   sba(i,2*j-1)=tmp1/dble(count)
                 else    
                   sba(i,2*j-1)=
     &             rtnorm(mu1,sd1,ll,uu,linf,lsup)
                end if
              else if(j.eq.nn)then
                linf=.false.
                lsup=.true.
                ll=sba(i-1,j-1)
                uu=0.d0

                tmp1=0.d0
                count=0
                do ii=1,nrec
                   if(y(ii).gt.ll)then
                      tmp1=tmp1+y(ii)
                      count=count+1
                   end if
                end do

                if(count.gt.0)then
                   sba(i,2*j-1)=tmp1/dble(count)
                 else    
                   sba(i,2*j-1)=
     &             rtnorm(mu1,sd1,ll,uu,linf,lsup)
                end if
              else
                linf=.false.
                lsup=.false.
                ll=sba(i-1,j-1)
                uu=sba(i-1,j)

                tmp1=0.d0
                count=0
                do ii=1,nrec
                   if(y(ii).gt.ll.and.y(ii).le.uu)then
                      tmp1=tmp1+y(ii)
                      count=count+1
                   end if
                end do

                if(count.gt.0)then
                   sba(i,2*j-1)=tmp1/dble(count)
                 else    
                   sba(i,2*j-1)=
     &             rtnorm(mu1,sd1,ll,uu,linf,lsup)
                end if
             end if   
          end do
       end do 

       return
       end


