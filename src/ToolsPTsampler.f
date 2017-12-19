c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR PTsampler
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
      subroutine ptsamplerwe3(np,nvar,nlevel,ntint,cpar,detlogl,z,
     &                        kphi,kphi2,kcount,kmat)
c=======================================================================
      implicit none
 
c++++ input
      integer np,nvar,nlevel,ntint
      double precision cpar,detlogl,z(np,nvar)

c++++ external working
      integer kphi(nvar)
      integer kphi2(nlevel)
      integer kcount(ntint,nlevel)

c++++ output
      integer kmat(np,nlevel)

c++++ internal workspace
      integer i,ind,j,j1,j2,k1
      double precision cdfnorm
      double precision tmp1,tmp2

c++++ algorithm

      do i=1,np

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do j=1,nvar
            tmp1=z(i,j)

            if(tmp1.gt.4.0d0)then
               tmp2=0.999968d0
             else if(tmp1.lt.-4.0d0)then
               tmp2=0.000032d0
             else 
               tmp2=cdfnorm(tmp1,0.d0,1.d0,1,0)
            end if
            kphi(j)=int((2**nlevel)*tmp2)
         end do

         do j1=1,nlevel
            j2=nlevel-j1+1

            ind=0
            do k1=1,nvar
               ind=ind+kphi(k1)*(2**((k1-1)*j2))
            end do
            ind=ind+1 

            if(ind.lt.1.or.ind.gt.ntint)then
               call rexit("error in set")
            end if

            kcount(ind,j2)=kcount(ind,j2)+1
            kmat(i,j2)=ind
 
            do k1=1,nvar
               kphi(k1)=int(kphi(k1)/2.0)
            end do
         end do
      end do


      return
      end
       

c=======================================================================      
      subroutine ptsamplerwe4(np,nvar,nlevel,ntint,cpar,z,detlogl,
     &                        kcount,kmat,lgw)
c=======================================================================
      implicit none
 
c++++ input
      integer np,nvar,nlevel,ntint
      double precision cpar,z(np,nvar),detlogl
      integer kcount(ntint,nlevel)
      integer kmat(np,nlevel)

c++++ output
      double precision lgw(np)

c++++ internal workspace
      integer n1,n2
      integer i,j
      double precision loglikn

c++++ algorithm

      do i=1,np

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         loglikn=0.d0
 
         do j=2,nlevel
            n1=kcount(kmat(i,j),j)   
            n2=kcount(kmat(i,j-1),j-1)
            loglikn=loglikn+
     &              log((2.0**nvar)*cpar*(j**2)+(2.0**nvar)*n1)-
     &              log((2.0**nvar)*cpar*(j**2)+n2)
         end do

         n1=kcount(kmat(i,1),1)
         loglikn=loglikn+
     &           log((2.0**nvar)*cpar+(2.0**nvar)*n1)-
     &           log((2.0**nvar)*cpar+np)

         loglikn=loglikn-0.5d0*detlogl

         do j=1,nvar
            loglikn=loglikn-0.5d0*z(i,j)*z(i,j)-0.9189385d0
         end do
         lgw(i)=loglikn
      end do

      return
      end
       

c=======================================================================      
      subroutine ptsamplerqe(np,nvar,nlevel,cpar,kmat,eval)
c=======================================================================
      implicit none
 
c++++ input
      integer np,nvar,nlevel
      double precision cpar
      integer kmat(np,nlevel)

c++++ output
      double precision eval

c++++ internal workspace
      integer n1,n2
      integer i,j

c++++ algorithm

      eval=0.d0

      do i=2,np
c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do j=2,nlevel
            n1=count(kmat(1:i-1,j)==kmat(i,j))   
            n2=count(kmat(1:i-1,j-1)==kmat(i,j-1))

            eval=eval+
     &              log((2.0**nvar)*cpar*(j**2)+(2.0**nvar)*n1)-
     &              log((2.0**nvar)*cpar*(j**2)+n2)
         end do

         n1=count(kmat(1:i-1,1)==kmat(i,1))

         eval=eval+
     &           log((2.0**nvar)*cpar+(2.0**nvar)*n1)-
     &           log((2.d0**nvar)*cpar+i-1)
      end do

      return
      end


c=======================================================================      
      subroutine ptsamplersam(narea,nvar,np,nlevel,cpar,m,u,z,
     &                        massi,mass,parti,pattern,patterns,whicho,
     &                        whichn,linf,lsup,limw,zwork,xwork)
c=======================================================================
      implicit none

c++++ input
      integer narea,nvar,np,nlevel
      double precision cpar,m(nvar),u(nvar,nvar),z(np,nvar)

c++++ external working space
      integer massi(narea)
      integer parti(nvar)
      integer pattern(nvar)
      integer patterns(nvar)
      integer whicho(np),whichn(np)

      double precision mass(narea)
      double precision linf(nvar),lsup(nvar),limw(nvar)
      double precision zwork(nvar),xwork(nvar)

c++++ internal working space
      integer binaryrep
      integer countero,countern
      integer evali,evali2
      integer nint
      integer i,j,je2,k,k1,k2,l,final

      double precision invcdfnorm
      double precision prob,quan,rtnorm
      double precision tmp1

c++++ algorithm
      nint=2
      prob=1.d0/dble(nint)
      quan=invcdfnorm(prob,0.d0,1.d0,1,0)

      countero=0
               
      do i=1,narea
         massi(i)=0
         mass(i)=0.d0
      end do
               
      do l=1,np
         do j=1,nvar
            evali=1 
            if(z(l,j).le.quan)evali=0  
            pattern(j)=evali
            zwork(j)=z(l,j)
         end do
         evali=binaryrep(nvar,pattern)
         massi(evali)=massi(evali)+1
      end do   

      do i=1,narea
         mass(i)=(cpar+dble(massi(i)))/
     &           ((2**nvar)*cpar+dble(np))
      end do

      call simdisc(mass,narea,narea,evali)  
      evali2=evali
      call binaryrepinv(nvar,evali2,patterns)

      do l=1,np
         final=1
         do j=1,nvar
            evali=1 
            if(z(l,j).le.quan)evali=0  
            pattern(j)=evali
            if(pattern(j).ne.patterns(j))final=0
         end do
               
         if(final.eq.1)then
            countero=countero+1
            whicho(countero)=l
         end if   
      end do 
               
      do i=1,nvar
         if(patterns(i).eq.0)then
            linf(i)=-999.d0
            lsup(i)=quan
            parti(i)=1
          else
            linf(i)=quan
            lsup(i)=999.d0
            parti(i)=2
         end if 
      end do

      countern=0
      do j=2,nlevel
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
                  
         do k=1,nvar
                     
            k1=2*(parti(k)-1)+1
            k2=2*(parti(k)-1)+2
            quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)

            limw(k)=quan
                     
            if(quan.gt.lsup(k).or.quan.lt.linf(k))then
               call intpr("j",-1,j,1)
               call dblepr("linf",-1,linf,nvar)
               call dblepr("lsup",-1,lsup,nvar)
               call dblepr("limw",-1,limw,nvar)
               call rexit("Errors in limits")
            end if
         end do   
                     
         do k=1,narea
            massi(k)=0
            mass(k)=0.d0
         end do
         
         if(countero.gt.0)then
                  do l=1,countero
               do k=1,nvar
                  evali=1 
                  if(z(whicho(l),k).le.limw(k))evali=0  
                  pattern(k)=evali
               end do
               evali=binaryrep(nvar,pattern)
               massi(evali)=massi(evali)+1
            end do                      
         end if 

         do l=1,narea
            mass(l)=(cpar*dble(je2)+dble(massi(l)))/
     &               ((2**nvar)*cpar*dble(je2)+dble(countero))
         end do

         call simdisc(mass,narea,narea,evali)  
         evali2=evali
         call binaryrepinv(nvar,evali2,patterns)

         countern=0
         if(countero.gt.0)then
             do l=1,countero
               final=1
               do k=1,nvar
                  evali=1 
                  if(z(whicho(l),k).le.limw(k))evali=0  
                  pattern(k)=evali
                  if(pattern(k).ne.patterns(k))final=0
               end do
               if(final.eq.1)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do                      
         end if
 
         do k=1,nvar
            if(patterns(k).eq.0)then
               parti(k)=2*(parti(k)-1)+1
               lsup(k)=limw(k)
             else
               parti(k)=2*(parti(k)-1)+2
               linf(k)=limw(k)
            end if 
         end do

         countero=countern
         if(countern.gt.0)then
            do l=1,countern
               whicho(l)=whichn(l)
            end do
         end if   
      end do


      do k=1,nvar
         if(linf(k).ge.lsup(k))then
            call intpr("j",-1,j,1)
            call dblepr("linf",-1,linf,nvar)
            call dblepr("lsup",-1,lsup,nvar)
            call rexit("Errors in limits")
         end if
      end do   

      do i=1,nvar
         zwork(i)=rtnorm(0.d0,1.d0,linf(i),
     &                   lsup(i),.false.,.false.)
      end do

      do j=1,nvar
         tmp1=0.d0
         do k=1,nvar
            tmp1=tmp1+u(j,k)*zwork(k)   
         end do
         xwork(j)=m(j)+tmp1
      end do

      return
      end


c=======================================================================      
      subroutine ptsamplermr(np,nvar,nlevel,ntint,cpar,detlogl,x,m,uinv,
     &                       kphi,kphi2,kcount,theta,eval)
c=======================================================================
      implicit none
 
c++++ input
      integer np,nvar,nlevel,ntint
      double precision cpar,detlogl,x(np,nvar),m(nvar),uinv(nvar,nvar)
      double precision theta(nvar)

c++++ external working
      integer kcount(ntint,nlevel)
      integer kphi(nvar)
      integer kphi2(nlevel)

c++++ output
      integer eval

c++++ internal workspace
      integer n1,n2
      integer i,ind,j,j1,j2,k,k1
      double precision cdfnorm
      double precision tmp1,tmp2

c++++ algorithm

      do i=1,np

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do j=1,nvar
            tmp1=0.d0
            do k=1,2
               tmp1=tmp1+uinv(j,k)*(x(i,k)-m(k))
            end do

            if(tmp1.gt.4.0d0)then
               tmp2=0.999968d0
             else if(tmp1.lt.-4.0d0)then
               tmp2=0.000032d0
             else 
               tmp2=cdfnorm(tmp1,0.d0,1.d0,1,0)
            end if
            kphi(j)=int(dble(2**nlevel)*tmp2)
         end do

         do j1=1,nlevel
            j2=nlevel-j1+1

            ind=0
            do k1=1,nvar
               ind=ind+kphi(k1)*(2**((k1-1)*j2))
            end do
            ind=ind+1 

            if(ind.lt.1.or.ind.gt.ntint)then
               call rexit("error in set")
            end if

            kcount(ind,j2)=kcount(ind,j2)+1
 
            do k1=1,nvar
               kphi(k1)=int(dble(kphi(k1))/2.0)
            end do
         end do
      end do

c++++ theta
      eval=0.d0

      do j=1,nvar
         tmp1=0.d0
         do k=1,2
            tmp1=tmp1+uinv(j,k)*(theta(k)-m(k))
         end do

         eval=eval-0.9189385d0+0.5d0*tmp1*tmp1

         if(tmp1.gt.4.0d0)then
            tmp2=0.999968d0
          else if(tmp1.lt.-4.0d0)then
            tmp2=0.000032d0
          else 
            tmp2=cdfnorm(tmp1,0.d0,1.d0,1,0)
         end if
         kphi(j)=int(dble(2**nlevel)*tmp2)
      end do

      do j1=1,nlevel
         j2=nlevel-j1+1

         ind=0
         do k1=1,nvar
            ind=ind+kphi(k1)*(2**((k1-1)*j2))
         end do
         ind=ind+1 

         if(ind.lt.1.or.ind.gt.ntint)then
            call rexit("error in set")
         end if

         kphi2(j2)=ind
 
         do k1=1,nvar
            kphi(k1)=int(dble(kphi(k1))/2.0)
         end do
      end do


      do j=2,nlevel
         n1=kcount(kphi2(j),j)   
         n2=kcount(kphi2(j-1),j-1)

         eval=eval+
     &        log((2.0**nvar)*cpar*(j**2)+(2.0**nvar)*n1)-
     &        log((2.0**nvar)*cpar*(j**2)+n2)
       end do

      n1=kcount(kphi2(1),1)
      eval=eval+
     &     log((2.0**nvar)*cpar+(2.0**nvar)*n1)-
     &     log((2.d0**nvar)*cpar+np)


      eval=eval-0.5d0*detlogl

      return
      end



