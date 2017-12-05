
c=======================================================================
      subroutine lddpmregconfuni(nsubject,p,
     &                           ncluster,ccluster,cstrt,ss,
     &                           y,x,sigma2k,
     &                           mub,sigmab,sigmabinv,alpha,
     &                           iflagc,prob,theta,
     &                           ztz,zty)
c=======================================================================      
c     Sample the configuration for the regression coefficients 
c     using the collapsed-state proposed by MacEachern (1998).
c     This function is for univariate responses.
c
c     Alejandro Jara, 2008-2010.
c=======================================================================            
      implicit none
c+++++Parameters
      integer nsubject,p
      integer ncluster
      integer ccluster(nsubject)
      integer cstrt(nsubject,nsubject)
      integer ss(nsubject)
      
      double precision alpha
      double precision y(nsubject)
      double precision x(nsubject,p)
      double precision mub(p)
      double precision sigma2k
      double precision sigmab(p,p)
      double precision sigmabinv(p,p)

c+++++Working external
      integer iflagc(p)
      double precision prob(nsubject+1)
      double precision theta(p)
      double precision ztz(p,p)
      double precision zty(p)

c+++++Working internal
      integer evali,i,j,k,l,m,ok
      integer ns
      integer since
      double precision tmp1,tmp2,tmp3
      double precision dnrm

c+++++Algorithm

      do i=1,p
         do j=1,p
            sigmabinv(i,j)=sigmab(i,j)
         end do
      end do 
      call inverse(sigmabinv,p,iflagc)      

      do i=1,nsubject

c+++++++ check if the user has requested an interrupt
         call rchkusr()
        
         ns=ccluster(ss(i))

c+++++++ subject in cluster with only 1 observation

         if(ns.eq.1)then
            since=ss(i)
            if(since.lt.ncluster)then
               call relabel_dpm(i,since,nsubject,ncluster,
     &                          ccluster,ss,cstrt)                   
            end if

            ccluster(ncluster)=ccluster(ncluster)-1 
            ncluster=ncluster-1
         end if   

c+++++++ subject in cluster with more than 1 observations
             
         if(ns.gt.1)then
            j=1
            ok=0
            do while(ok.eq.0.and.j.le.ns)
               if(cstrt(ss(i),j).eq.i)ok=j
               j=j+1
            end do
   
            do j=ok,ns-1
               cstrt(ss(i),j)=cstrt(ss(i),j+1)
            end do
          
            ccluster(ss(i))=ccluster(ss(i))-1 
         end if

c+++++++ updating configurations

         do j=1,ncluster

            do k=1,p
               tmp1=0.d0
               do l=1,p
                  ztz(k,l)=sigmabinv(k,l)
                  tmp1=tmp1+sigmabinv(k,l)*mub(l)
               end do
               zty(k)=tmp1
            end do

            do k=1,ccluster(j)
               do l=1,p
                  do m=1,p   
                     ztz(l,m)=ztz(l,m)+
     &                        x(cstrt(j,k),l)*(1.d0/sigma2k)* 
     &                        x(cstrt(j,k),m)
                  end do
                  zty(l)=zty(l)+
     &                      x(cstrt(j,k),l)*(1.d0/sigma2k)* 
     &                      y(cstrt(j,k))
               end do
            end do 

            call inverse(ztz,p,iflagc) 
               
            do k=1,p
               tmp1=0.d0
               do l=1,p
                  tmp1=tmp1+ztz(k,l)*zty(l)                     
               end do
               theta(k)=tmp1
            end do

            tmp1=0.d0
            do l=1,p
               tmp1=tmp1+x(i,l)*theta(l)
            end do
            
            tmp2=0.d0
            do l=1,p
               do m=1,p   
                  tmp2=tmp2+
     &                 x(i,l)*ztz(l,m)*x(i,m)
               end do
            end do
            tmp2=tmp2+sigma2k
            tmp3=dnrm(y(i),tmp1,sqrt(tmp2),0)
            prob(j)=dble(ccluster(j))*tmp3
         end do

         tmp1=0.d0
         do l=1,p
            tmp1=tmp1+x(i,l)*mub(l)
         end do

         tmp2=0.d0
         do l=1,p
            do m=1,p   
               tmp2=tmp2+
     &              x(i,l)*ztz(l,m)*x(i,m)
            end do
         end do
         tmp2=tmp2+sigma2k
         tmp3=dnrm(y(i),tmp1,sqrt(tmp2),0)
         prob(ncluster+1)=alpha*tmp3

         call simdisc(prob,nsubject+1,ncluster+1,evali)

         if(evali.lt.1.or.evali.gt.(ncluster+1))then
            call intpr("ncluster",-1,ncluster,1)  
            call intpr("evali",-1,evali,1)  
            call dblepr("prob",-1,prob,ncluster+1)  
            call rexit("error in clustering location")         
         end if
         
         ss(i)=evali
         ccluster(evali)=ccluster(evali)+1
         cstrt(evali,ccluster(evali))=i
          
         if(evali.gt.ncluster)then
            ncluster=ncluster+1
         end if

      end do

      return
      end


c=======================================================================
      subroutine lddpmreglocuni(nsubject,p,
     &                          ncluster,ccluster,cstrt,
     &                          y,x,sigma2k,
     &                          mub,sigmabinv,b,
     &                          iflagc,theta,
     &                          workvc,workmhc,
     &                          ztz,zty)
c=======================================================================      
c     Resample the clustered regression coefficients given the 
c     cluster configurations.
c     This function if for univariate responses. 
c
c     Alejandro Jara, 2008
c=======================================================================            
      implicit none
c+++++Parameters
      integer nsubject,p
      integer ncluster
      integer ccluster(nsubject)
      integer cstrt(nsubject,nsubject)
      
      double precision y(nsubject)
      double precision x(nsubject,p)
      
      double precision b(nsubject,p)
      double precision sigma2k
      double precision mub(p)
      double precision sigmabinv(p,p)

c+++++Working external
      integer iflagc(p)
      double precision theta(p)
      double precision workvc(p)
      double precision workmhc(p*(p+1)/2)
      double precision ztz(p,p)
      double precision zty(p)

c+++++Working internal
      integer j,k,l,m
      integer ns
      double precision tmp1

c+++++Algorithm

      do j=1,ncluster

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         ns=ccluster(j)

         do k=1,p
            tmp1=0.d0
            do l=1,p
               ztz(k,l)=sigmabinv(k,l)
               tmp1=tmp1+sigmabinv(k,l)*mub(l)
            end do
            zty(k)=tmp1
         end do

         do k=1,ns
            do l=1,p
               do m=1,p   
                  ztz(l,m)=ztz(l,m)+
     &                     x(cstrt(j,k),l)*(1.d0/sigma2k)* 
     &                     x(cstrt(j,k),m)
               end do
               zty(l)=zty(l)+
     &                    x(cstrt(j,k),l)*(1.d0/sigma2k)* 
     &                    y(cstrt(j,k))
            end do
         end do 

         call inverse(ztz,p,iflagc) 
               
         do k=1,p
            tmp1=0.d0
            do l=1,p
               tmp1=tmp1+ztz(k,l)*zty(l)                     
            end do
            workvc(k)=tmp1
         end do

         call rmvnorm(p,workvc,ztz,workmhc,zty,theta)

         do k=1,p
            b(j,k)=theta(k)
         end do

      end do

      return
      end


c=======================================================================      
      subroutine relabel_dpm(ind,since,nsubject,ncluster,ccluster,ss,
     &                      cstrt)
c=======================================================================
c     relabel the clusters after the elimination of one of them
c     A.J.V., 2007
      implicit none
      integer i,j,ind,since,nsubject,ncluster,ccluster(nsubject)
      integer ss(nsubject),cstrt(nsubject,nsubject)
      
      do i=since+1,ncluster
         do j=1,nsubject
            if(ss(j).eq.i)then
               ss(j)=i-1 
            end if
         end do
         do j=1,ccluster(i)
            cstrt(i-1,j)=cstrt(i,j) 
         end do
         ccluster(i-1)=ccluster(i)
      end do
      
      ss(ind)=ncluster
      ccluster(ncluster)=1
      
      return
      end  


c=======================================================================
      subroutine lddpmregconfmult(nsubject,nvar,p,
     &                            ncluster,ccluster,cstrt,ss,
     &                            possi,y,x,sigmak,sigmakinv,
     &                            mub,sigmab,sigmabinv,alpha,
     &                            iflagn,iflagc,prob,theta,
     &                            workmn,workvn,workvn2,
     &                            ztz,zty)
c=======================================================================      
c     Sample the configuration for the regression coefficients 
c     using the collapsed-state proposed by MacEachern (1998).
c     This function is for multivariate responses.
c     The function assumes that the same vector of predictors is
c     used for each coordinate.
c
c     Alejandro Jara, 2008 
c=======================================================================            
      implicit none
c+++++Parameters
      integer nsubject,nvar,p
      integer ncluster
      integer ccluster(nsubject)
      integer cstrt(nsubject,nsubject)
      integer ss(nsubject)
      integer possi(nvar,p)
      
      double precision alpha
      double precision y(nsubject,nvar)
      double precision x(nsubject,p)
      double precision mub(nvar*p)
      double precision sigmak(nvar,nvar)      
      double precision sigmakinv(nvar,nvar)            
      double precision sigmab(nvar*p,nvar*p)
      double precision sigmabinv(nvar*p,nvar*p)

c+++++Working external
      integer iflagn(nvar)
      integer iflagc(nvar*p)
      double precision prob(nsubject+1)
      double precision theta(nvar*p)
      double precision workmn(nvar,nvar)      
      double precision workvn(nvar)
      double precision workvn2(nvar)
      double precision ztz(nvar*p,nvar*p)
      double precision zty(nvar*p)

c+++++Working internal
      integer evali,i,ii,j,jj,k,l,m,ok
      integer ns
      integer pos1,pos2
      integer since
      double precision tmp1,tmp2

c+++++Algorithm

      do i=1,nsubject

c+++++++ check if the user has requested an interrupt
         call rchkusr()
        
         ns=ccluster(ss(i))

c+++++++ subject in cluster with only 1 observation

         if(ns.eq.1)then
            since=ss(i)
            if(since.lt.ncluster)then
               call relabel_dpm(i,since,nsubject,ncluster,
     &                          ccluster,ss,cstrt)                   
            end if

            ccluster(ncluster)=ccluster(ncluster)-1 
            ncluster=ncluster-1
         end if   

c+++++++ subject in cluster with more than 1 observations
             
         if(ns.gt.1)then
            j=1
            ok=0
            do while(ok.eq.0.and.j.le.ns)
               if(cstrt(ss(i),j).eq.i)ok=j
               j=j+1
            end do
   
            do j=ok,ns-1
               cstrt(ss(i),j)=cstrt(ss(i),j+1)
            end do
          
            ccluster(ss(i))=ccluster(ss(i))-1 
         end if

c+++++++ updating configurations

         do j=1,ncluster

            do k=1,nvar*p
               tmp1=0.d0
               do l=1,nvar*p
                  ztz(k,l)=sigmabinv(k,l)
                  tmp1=tmp1+sigmabinv(k,l)*mub(l)
               end do
               zty(k)=tmp1
            end do

            do k=1,ccluster(j)
               do ii=1,nvar
                  do l=1,p
                     pos1=possi(ii,l)
                     do jj=1,nvar 
                        do m=1,p   
                           pos2=possi(jj,m)
                           ztz(pos1,pos2)=ztz(pos1,pos2)+
     &                        x(cstrt(j,k),l)*sigmakinv(ii,jj)* 
     &                        x(cstrt(j,k),m)
                        end do
                        zty(pos1)=zty(pos1)+
     &                       x(cstrt(j,k),l)*sigmakinv(ii,jj)* 
     &                       y(cstrt(j,k),jj)
                     end do
                  end do
               end do
            end do 

            call inverse(ztz,nvar*p,iflagc) 
               
            do k=1,nvar*p
               tmp1=0.d0
               do l=1,nvar*p
                  tmp1=tmp1+ztz(k,l)*zty(l)                     
               end do
               theta(k)=tmp1
            end do

            do k=1,nvar
               workvn(k)=y(i,k)
               tmp1=0.d0
               do l=1,p
                  tmp1=tmp1+x(i,l)*theta(possi(k,l))
               end do
               workvn2(k)=tmp1
               
               do l=1,nvar
                  workmn(k,l)=sigmak(k,l)
               end do
            end do

            do ii=1,nvar
               do l=1,p
                  pos1=possi(ii,l)
                  do jj=1,nvar 
                     do m=1,p
                        pos2=possi(jj,m)
                        workmn(ii,jj)=workmn(ii,jj)+
     &                          x(i,l)*ztz(pos1,pos2)* 
     &                          x(i,m)
                     end do
                  end do
               end do
            end do

            tmp2=0.d0
            call dmvnd(nvar,workvn,workvn2,workmn,tmp2,iflagn)  
            prob(j)=dble(ccluster(j))*exp(tmp2)
         end do

         do k=1,nvar
            workvn(k)=y(i,k)
            
            tmp1=0.d0
            do l=1,p
               tmp1=tmp1+x(i,l)*mub(possi(k,l))
            end do
            workvn2(k)=tmp1

            do l=1,nvar
               workmn(k,l)=sigmak(k,l)
            end do
         end do

         do ii=1,nvar
            do l=1,p
               pos1=possi(ii,l)
               do jj=1,nvar 
                  do m=1,p
                     pos2=possi(jj,m)
                     workmn(ii,jj)=workmn(ii,jj)+
     &                       x(i,l)*sigmab(pos1,pos2)* 
     &                       x(i,m)
                  end do
               end do
            end do
         end do
         
         tmp2=0.d0
         call dmvnd(nvar,workvn,workvn2,workmn,tmp2,iflagn)        
           
         prob(ncluster+1)=alpha*exp(tmp2)

         call simdisc(prob,nsubject+1,ncluster+1,evali)

         if(evali.lt.1.or.evali.gt.(ncluster+1))then
            call intpr("ncluster",-1,ncluster,1)  
            call intpr("evali",-1,evali,1)  
            call dblepr("prob",-1,prob,ncluster+1)  
            call rexit("error in clustering location")         
         end if
         
         ss(i)=evali
            
         ccluster(evali)=ccluster(evali)+1
            
         cstrt(evali,ccluster(evali))=i
          
         if(evali.gt.ncluster)then
            ncluster=ncluster+1
         end if

      end do

      return
      end


c=======================================================================
      subroutine lddpmreglocmult(nsubject,nvar,p,
     &                           ncluster,ccluster,cstrt,
     &                           possi,y,x,sigmakinv,
     &                           mub,sigmabinv,b,
     &                           iflagc,theta,
     &                           workvc,workmhc,
     &                           ztz,zty)
c=======================================================================      
c     Resample the clustered regression coefficients given the 
c     cluster configurations.
c     This function if for multivariate responses. 
c
c     Alejandro Jara, 2008
c=======================================================================            
      implicit none
c+++++Parameters
      integer nsubject,nvar,p
      integer ncluster
      integer ccluster(nsubject)
      integer cstrt(nsubject,nsubject)
      integer possi(nvar,p)
      
      double precision y(nsubject,nvar)
      double precision x(nsubject,p)
      
      double precision b(nsubject,nvar*p)
      double precision sigmakinv(nvar,nvar)            
      double precision mub(nvar*p)
      double precision sigmabinv(nvar*p,nvar*p)

c+++++Working external
      integer iflagc(nvar*p)
      double precision theta(nvar*p)
      double precision workvc(nvar*p)
      double precision workmhc(nvar*p*(nvar*p+1)/2)
      double precision ztz(nvar*p,nvar*p)
      double precision zty(nvar*p)

c+++++Working internal
      integer ii,j,jj,k,l,m
      integer ns
      integer pos1,pos2
      double precision tmp1

c+++++Algorithm

      do j=1,ncluster

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         ns=ccluster(j)

         do k=1,nvar*p
            tmp1=0.d0
            do l=1,nvar*p
               ztz(k,l)=sigmabinv(k,l)
               tmp1=tmp1+sigmabinv(k,l)*mub(l)
            end do
            zty(k)=tmp1
         end do

         do k=1,ns
            do ii=1,nvar
               do l=1,p
                  pos1=possi(ii,l)
                  do jj=1,nvar 
                     do m=1,p   
                        pos2=possi(jj,m)
                        ztz(pos1,pos2)=ztz(pos1,pos2)+
     &                          x(cstrt(j,k),l)*sigmakinv(ii,jj)* 
     &                          x(cstrt(j,k),m)
                     end do
                     zty(pos1)=zty(pos1)+
     &                    x(cstrt(j,k),l)*sigmakinv(ii,jj)* 
     &                    y(cstrt(j,k),jj)
                  end do
               end do
            end do
         end do 

         call inverse(ztz,nvar*p,iflagc) 
               
         do k=1,nvar*p
            tmp1=0.d0
            do l=1,nvar*p
               tmp1=tmp1+ztz(k,l)*zty(l)                     
            end do
            workvc(k)=tmp1
         end do

         call rmvnorm(nvar*p,workvc,ztz,workmhc,zty,theta)

         do k=1,nvar*p
            b(j,k)=theta(k)
         end do

      end do

      return
      end

