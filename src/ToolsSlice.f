c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR SLICE SAMPLING
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
      subroutine vertl_sigma2(x0,nsubject,kk,b,z,mu,tau1,tau2,zz)
c=======================================================================
c     picks a value for the auxiliary variable of the salice sampler
c     for sigma2 in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================
      implicit none  

c+++++Arguments function
      integer nsubject,kk
      double precision x0,b(nsubject),z(nsubject),mu,zz
      double precision tau1,tau2

c+++++Slice sampling working variables
      double precision eval,rexp

c+++++Algorithm      

      call logpsigma2(nsubject,kk,b,z,mu,x0,tau1,tau2,eval)      

      zz=eval-rexp(1.d0)
      
      return
      end     


c=======================================================================      
      subroutine dbling_sigma2(x0,zz,w,p,nsubject,kk,b,z,mu,
     &                         tau1,tau2,l,r)
c=======================================================================      
c     Radford Neal's (2000) doubling procedure for sigma2 in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none

c+++++Arguments function
      integer nsubject,kk
      double precision x0,b(nsubject),z(nsubject),mu,tau1,tau2

c+++++Slice sampling working variables
      integer p,k
      double precision zz,w
      double precision u,v,l,r,tl,tr
      double precision ml
      real runif

c+++++Algorithm
      u=dble(runif())
      l=x0-w*u
      r=l+w
      k=p
      
      ml=0
      
      if(l.lt.0.000001d0)then
        l=0.000001d0
        ml=1
      end if  
      
      call logpsigma2(nsubject,kk,b,z,mu,l,tau1,tau2,tl)      
      call logpsigma2(nsubject,kk,b,z,mu,r,tau1,tau2,tr)      

      do while(k.gt.0.and.(zz.lt.tl.or.zz.lt.tr))
         v=dble(runif())
         if(v.lt.0.5d0)then
            if(ml.ne.1)then 
               l=l-(r-l)
               if(l.lt.0.000001d0)then
                 l=0.000001d0
                 ml=1
               end if  
               call logpsigma2(nsubject,kk,b,z,mu,l,tau1,tau2,tl)      
            end if   
         end if   
         if(v.ge.0.5d0)then
            r=r+(r-l)
            call logpsigma2(nsubject,kk,b,z,mu,r,tau1,tau2,tr)
         end if   
         k=k-1
      end do   
 
      return
      end

c=======================================================================      
      subroutine shrink_sigma2(x0,zz,l,r,w,nsubject,kk,b,z,mu,tau1,
     &                         tau2,x1)
c=======================================================================      
c     Radford Neal's (2000) shrinkage procedure for sigma2 in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none

c+++++Arguments function
      integer nsubject,kk
      double precision x0,b(nsubject),z(nsubject),mu
      double precision tau1,tau2

c+++++Slice sampling working variables
      integer accept,out
      double precision zz,l,r,w,u
      double precision x1,lbar,rbar,tmp1
      real runif

c+++++Algorithm

      lbar=l
      rbar=r
      accept=0

      do while(accept.eq.0)
         u=dble(runif())
         x1=lbar+u*(rbar-lbar)

         call logpsigma2(nsubject,kk,b,z,mu,x1,tau1,tau2,tmp1)
         if(zz.lt.tmp1)accept=1

         if(accept.eq.1)then  
            call test_sigma2(x0,x1,zz,l,r,w,nsubject,kk,b,z,mu,
     &                       tau1,tau2,out)
            if(out.eq.1)then
               accept=1
              else
               accept=0
            end if   
         end if
         
         if(accept.ne.1)then
           if(x1.lt.x0)then
              lbar=x1  
            else  
              rbar=x1
           end if   
         end if
      end do   
 
      return
      end

c=======================================================================      
      subroutine test_sigma2(x0,x1,zz,l,r,w,nsubject,kk,b,z,mu,
     &                       tau1,tau2,out)
c=======================================================================      
c     Radford Neal's (2000) testing procedure for sigma2 in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none

c+++++Arguments function
      integer nsubject,kk,out
      double precision x0,x1,b(nsubject),z(nsubject),mu
      double precision tau1,tau2
      double precision zz,l,r,w

c+++++Slice sampling working variables
      integer out2 
      double precision m,lbar,rbar,tmp1,tmp2

c+++++Algorithm
      lbar=l
      rbar=r
      out=0
      out2=0
      
      call logpsigma2(nsubject,kk,b,z,mu,rbar,tau1,tau2,tmp1)
      call logpsigma2(nsubject,kk,b,z,mu,lbar,tau1,tau2,tmp2)

      do while((rbar-lbar).gt.1.1d0*w.and.out2.eq.0)
         m=(rbar+lbar)/2.d0
         
         if((x0.lt.m.and.x1.ge.m).or.(x0.ge.m.and.x1.lt.m))out=1
         
         if(x1.lt.m)then
            rbar=m
            call logpsigma2(nsubject,kk,b,z,mu,rbar,tau1,tau2,tmp1)
          else
            lbar=m
            call logpsigma2(nsubject,kk,b,z,mu,lbar,tau1,tau2,tmp2)
         end if
         
         if(out.eq.1)then
            if(zz.ge.tmp2.and.zz.ge.tmp1)out2=1
         end if
      end do
      
      if(out2.eq.1)out=0
      if(out2.eq.0)out=1
 
      return
      end


c=======================================================================      
      subroutine slice_sigma2(sigma2,w,p,nsubject,kk,b,z,mu,
     &                        tau1,tau2,l,r)
c=======================================================================      
c     Radford Neal's (2000) slice sampling for sigma2 in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none

c+++++Arguments function
      integer nsubject,kk,p
      double precision sigma2,b(nsubject),z(nsubject),mu,tau1,tau2,w
      double precision l,r

c+++++Slice sampling working variables
      double precision x0,zz
 
c+++++Algorithm
      x0=sigma2

      call vertl_sigma2(x0,nsubject,kk,b,z,mu,tau1,tau2,zz)        

      call dbling_sigma2(x0,zz,w,p,nsubject,kk,b,z,mu,tau1,tau2,l,r)
      
      call shrink_sigma2(x0,zz,l,r,w,nsubject,kk,b,z,mu,
     &                   tau1,tau2,sigma2)      
 
      return
      end


c=======================================================================
      subroutine vertl_sigma(x0,nsubject,kk,b,z,mu,zz)        
c=======================================================================
c     picks a value for the auxiliary variable of the salice sampler
c     for sigma in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================
      implicit none  

c+++++Arguments function
      integer nsubject,kk
      double precision x0,b(nsubject),z(nsubject),mu,zz

c+++++Slice sampling working variables
      double precision eval,rexp

c+++++Algorithm      

      call logpsigma(nsubject,kk,b,z,mu,x0,eval)      

      zz=eval-rexp(1.d0)
      
      return
      end     

c=======================================================================      
      subroutine dbling_sigma(x0,zz,w,p,nsubject,kk,b,z,mu,
     &                        maxsigma,l,r)
c=======================================================================      
c     Radford Neal's (2000) doubling procedure for sigma in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none

c+++++Arguments function
      integer nsubject,kk
      double precision x0,b(nsubject),z(nsubject),mu,maxsigma

c+++++Slice sampling working variables
      integer p,k
      double precision zz,w
      double precision u,v,l,r,tl,tr
      double precision ml,mr
      real runif

c+++++Algorithm
      u=dble(runif())
      l=x0-w*u
      r=l+w
      k=p
      
      ml=0
      mr=0
      
      if(l.lt.0.000001d0)then
        l=0.000001d0
        ml=1
      end if  
      if(r.ge.maxsigma)then
        r=maxsigma
        mr=1
      end if  
      
      call logpsigma(nsubject,kk,b,z,mu,l,tl)      
      call logpsigma(nsubject,kk,b,z,mu,r,tr)      
      
      do while(k.gt.0.and.(zz.lt.tl.or.zz.lt.tr).and.(ml+mr).lt.2)
         v=dble(runif())
         if(v.lt.0.5d0)then
            if(ml.ne.1)then
               l=l-(r-l)
               if(l.lt.0.000001d0)then
                 l=0.000001d0
                 ml=1
               end if  
               call logpsigma(nsubject,kk,b,z,mu,l,tl)      
            end if   
         end if   
         if(v.ge.0.5d0)then
            if(mr.ne.1)then
               r=r+(r-l)
               if(r.gt.maxsigma)then
                  r=maxsigma
                  mr=1
               end if   
               call logpsigma(nsubject,kk,b,z,mu,r,tr)
            end if   
         end if   
         k=k-1
      end do   
 
      return
      end


c=======================================================================      
      subroutine shrink_sigma(x0,zz,l,r,w,nsubject,kk,b,z,mu,x1)
c=======================================================================      
c     Radford Neal's (2000) shrinkage procedure for sigma in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none

c+++++Arguments function
      integer nsubject,kk
      double precision x0,b(nsubject),z(nsubject),mu

c+++++Slice sampling working variables
      integer accept,out
      double precision zz,l,r,w,u
      double precision x1,lbar,rbar,tmp1
      real runif

c+++++Algorithm

      lbar=l
      rbar=r
      accept=0

      do while(accept.eq.0)
         u=dble(runif())
         x1=lbar+u*(rbar-lbar)

         call logpsigma(nsubject,kk,b,z,mu,x1,tmp1)
         if(zz.lt.tmp1)accept=1

         call test_sigma(x0,x1,zz,l,r,w,nsubject,kk,b,z,mu,out)
         if(accept.eq.1.and.out.eq.1)accept=1 
         
         if(accept.ne.1)then
           if(x1.lt.x0)lbar=x1  
           if(x1.ge.x0)rbar=x1
         end if
      end do   
 
      return
      end

c=======================================================================      
      subroutine test_sigma(x0,x1,zz,l,r,w,nsubject,kk,b,z,mu,out)
c=======================================================================      
c     Radford Neal's (2000) testing procedure for sigma in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none

c+++++Arguments function
      integer nsubject,kk,out
      double precision x0,x1,b(nsubject),z(nsubject),mu
      double precision zz,l,r,w

c+++++Slice sampling working variables
      integer out2 
      double precision m,lbar,rbar,tmp1,tmp2


c+++++Algorithm
      lbar=l
      rbar=r
      out=0
      out2=0

      do while((rbar-lbar).gt.1.1d0*w.and.out2.eq.0)
         m=(rbar+lbar)/2.d0
         
         if((x0.lt.m.and.x1.ge.m).or.(x0.ge.m.and.x1.lt.m))out=1
         
         if(x1.lt.m)then
            rbar=m
            call logpsigma(nsubject,kk,b,z,mu,rbar,tmp1)
          else
            lbar=m
            call logpsigma(nsubject,kk,b,z,mu,lbar,tmp2)
         end if
         
         if(out.eq.1)then
            if(zz.ge.tmp2.and.zz.ge.tmp1)out2=1
         end if
      end do
      
      if(out2.eq.1)out=0
      if(out2.eq.0)out=1
 
      return
      end


c=======================================================================      
      subroutine slice_sigma(sigma,w,p,nsubject,kk,b,z,mu,
     &                       maxsigma,l,r)
c=======================================================================      
c     Radford Neal's (2000) slice sampling for sigma in a BDP.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none

c+++++Arguments function
      integer nsubject,kk,p
      double precision sigma,b(nsubject),z(nsubject),mu,maxsigma,w
      double precision l,r

c+++++Slice sampling working variables
      double precision x0,zz
 
c+++++Algorithm
      x0=sigma

      call vertl_sigma(x0,nsubject,kk,b,z,mu,zz)        

      call dbling_sigma(x0,zz,w,p,nsubject,kk,b,z,mu,maxsigma,l,r)
      
      call shrink_sigma(x0,zz,l,r,w,nsubject,kk,b,z,mu,sigma)      
 
      return
      end

