c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR ARS
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
      double precision function ssexp(x)
c=======================================================================
      real*8 x,emax
      parameter(emax=50.d0)

      if(x.lt.-emax)then
        ssexp=0.d0
       else 
        ssexp=dexp(min(x,emax))
      end if
      return
      end


c=======================================================================                  
c=======================================================================                  
c=======================================================================                  
c     DERIVATIVE-DEPENDENT (WORKING OKAY)
c=======================================================================                  
c=======================================================================                  
c=======================================================================                  

c=======================================================================
      subroutine arsCompuHunD(n,maxn,maxl,maxu,lHb,lHm,lHl,lHr,
     &                        uHb,uHm,uHl,uHr,uHpr,grid,fs,fps,fsmax)
c=======================================================================
c     computes the picewise-linear hull in an ARS scheme for a
c     distribution with unbounded support based on derivatives.
c     A.J.V., 2007

      implicit none

c++++ dimensions
      integer n
      integer maxn,maxl,maxu
      
c++++ lower Hull
      real*8 lHb(maxl),lHm(maxl)
      real*8 lHl(maxl),lHr(maxl)

c++++ upper Hull
      real*8 uHb(maxu),uHm(maxu)
      real*8 uHl(maxu),uHr(maxu)
      real*8 uHpr(maxu)

c++++ evaluations
      real*8 grid(maxn),fs(maxn),fps(maxn)

c++++ working space
      integer i,j,totall,totalu,maxcheck
      parameter(maxcheck=2000)
      integer checki(maxcheck),totali,counter
      real*8 accum,dh,b,m,pr
      real*8 b1,b2,m1,m2
      real*8 ix,hix,pr1,pr2
      real*8 tmp1,tmp2,ssexp
      real*8 zero,fsmax,emax
      parameter(emax=50.d0)
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zero=ssexp(-emax)

      if(maxn.gt.maxcheck)then
        call rexit("Error in 'arsCompuHnu' (maxcheck)")
      end if

c++++ checking -Inf before start
      totali=0
      do i=1,n
         checki(i)=0
         if(fs(i).lt.-9999999999999999999999.d0)then
            checki(i)=1
            totali=totali+1
         end if
      end do

      counter=0
      if(totali.gt.0)then
        do i=1,n
           if(checki(i).eq.0)then
             counter=counter+1
             grid(counter)=grid(i)
             fs(counter)=fs(i)
           end if
        end do
        n=n-totali
      end if

      if(n.lt.3)then
        call rexit("Error in 'arsCompuHnu': number of points < 3")
      end if

c++++ end checking
      
      totall=n-1
      totalu=2+2*(n-1) 

      if(totalu.gt.maxu)then
        call rexit("Error in 'arsCompuHnu' (maxu)")
      end if

      if(totall.gt.maxl)then
        call rexit("Error in 'arsCompuHnu' (maxl)")
      end if

      do i=1,maxu
         uHm(i)=0.d0
         uHb(i)=0.d0
         uHl(i)=0.d0
         uHr(i)=0.d0
         uHpr(i)=0.d0
      end do

      do i=1,maxl
         lHb(i)=0.d0
         lHm(i)=0.d0
         lHl(i)=0.d0
         lHr(i)=0.d0
      end do

c++++ computing the lines of the Lower Hull

      do i=1,totall
         lHm(i)=(fs(i+1)-fs(i))/(grid(i+1)-grid(i))
         lHb(i)=fs(i)-lHm(i)*grid(i)
         lHl(i)=grid(i)
         lHr(i)=grid(i+1)
      end do

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ computing the lines of the Upper Hull
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      fsmax=fs(1)
c++++ first line
      m=fps(1)
      b=fs(1)-m*grid(1)

      uHm(1)=m
      uHb(1)=b
      uHl(1)=-999999.d0
      uHr(1)=grid(1)

c++++ interior lines 
      j=1
      if(n.gt.2)then
        do i=1,(n-1)
c+++++++++ first interior line
           m1= fps(i)
           b1= fs(i)-m1*grid(i)
           fsmax=max(fs(i),fsmax)

c+++++++++ second interior line
           m2= fps(i+1)
           b2= fs(i+1)-m2*grid(i+1)
           fsmax=max(fs(i+1),fsmax)          

c+++++++++ intersection
           if(abs(m2-m1).lt.zero)then
             ix=0.5*(grid(i)+grid(i+1))
             hix=0.5*(fs(i)+fs(i+1))
            else if(abs(m1).lt.abs(m2))then
             ix=grid(i+1)+
     &          (fs(i)-fs(i+1)+m1*(grid(i+1)-grid(i)))/(m2-m1)
             hix=m1*(ix-grid(i))+fs(i)
            else            
             ix=grid(i)+
     &          (fs(i)-fs(i+1)+m2*(grid(i+1)-grid(i)))/(m2-m1)
             hix=m2*(ix-grid(i+1))+fs(i+1)
           end if  

           if(ix.lt.grid(i).or.ix.gt.grid(i+1))then
              call dblepr("ix",-1,ix,1)
              call dblepr("xi",-1,grid(i),1)
              call dblepr("xi+1",-1,grid(i+1),1)
              call rexit("Error in intersection")
           end if
           fsmax=max(hix,fsmax)

c+++++++++ update first interior line           
           j=j+1
           uHm(j)=m1 
           uHb(j)=b1
           uHl(j)=grid(i)
           uHr(j)=ix

c+++++++++ update second interior line           
           j=j+1
           uHm(j)=m2 
           uHb(j)=b2
           uHl(j)=ix
           uHr(j)=grid(i+1)
        end do
      end if

c++++ last line
      m=fps(n)
      b=fs(n)-m*grid(n)
      fsmax=max(fs(n),fsmax)
      j=j+1
      uHm(j)=m
      uHb(j)=b
      uHl(j)=grid(n)
      uHr(j)=+999999.d0

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ computing the areas of the Upper Hull
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      accum=0.d0
c++++ first line

      m=fps(1)
      if(abs(m).lt.zero)then
         pr=0.d0
        else
         pr=ssexp(fs(1)-fsmax)/m
      end if
      
      uHpr(1)=pr
      accum=accum+pr

c++++ interior lines 
      j=1
      if(n.gt.2)then
        do i=1,(n-1)

c+++++++++ intersection

           m1= fps(i)
           m2= fps(i+1)

           if(abs(m2-m1).lt.zero)then
             ix=0.5*(grid(i)+grid(i+1))
             hix=0.5*(fs(i)+fs(i+1))
            else if(abs(m1).lt.abs(m2))then
             ix=grid(i+1)+
     &          (fs(i)-fs(i+1)+m1*(grid(i+1)-grid(i)))/(m2-m1)
             hix=m1*(ix-grid(i))+fs(i)
            else            
             ix=grid(i)+
     &          (fs(i)-fs(i+1)+m2*(grid(i+1)-grid(i)))/(m2-m1)
             hix=m2*(ix-grid(i+1))+fs(i+1)
           end if  

           if(ix.lt.grid(i).or.ix.gt.grid(i+1))then
              call dblepr("ix",-1,ix,1)
              call dblepr("xi",-1,grid(i),1)
              call dblepr("xi+1",-1,grid(i+1),1)
              call rexit("Error in intersection")
           end if

c+++++++++ areas

           tmp1=fs(i)
           tmp2=hix
           dh=tmp1-tmp2

           if(abs(m1).lt.zero)then
              pr1=(ix-grid(i))*ssexp((tmp2+tmp1)*0.5-fsmax)
             else 
              if(dh.lt.emax)then 
                pr1=ssexp(tmp2-fsmax)*(1-ssexp(dh))/m1
               else
                pr1=-ssexp(tmp1-fsmax)/m1
              end if 
           end if             

           tmp1=hix
           tmp2=fs(i+1)
           dh=tmp1-tmp2

           if(abs(m2).lt.zero)then
              pr2=(grid(i+1)-ix)*ssexp((tmp2+tmp1)*0.5-fsmax)
             else 
              if(dh.lt.emax)then 
                pr2=ssexp(tmp2-fsmax)*(1-ssexp(dh))/m2
               else
                pr2=-ssexp(tmp1-fsmax)/m2
              end if 
           end if             

c+++++++++ update first interior line           
           j=j+1
           uHpr(j)=pr1
           accum=accum+pr1

c+++++++++ update second interior line           
           j=j+1
           uHpr(j)=pr2
           accum=accum+pr2

        end do
      end if

c++++ last line

      m=fps(n)
      if(abs(m).lt.zero)then
         pr=0.d0
        else
         pr=-ssexp(fs(n)-fsmax)/m
      end if

      j=j+1
      uHpr(j)=pr
      accum=accum+pr

      if(accum.lt.0.d0)then
        call dblepr("s",-1,grid,n)
        call dblepr("fs",-1,fs,n)
        call dblepr("upperHpr",-1,uHpr,totalu)
        call dblepr("upperHm",-1,uHm,totalu)
        call dblepr("upperHb",-1,uHb,totalu)
        call dblepr("fsmax",-1,fsmax,1)
        call rexit("Negative total area in ARS")
      end if

c++++ standardizing the areas

      do i=1,totalu
         if(uHpr(i).lt.0.d0)then
           call dblepr("s",-1,grid,n)
           call dblepr("fs",-1,fs,n)
           call dblepr("upperHpr",-1,uHpr,totalu)
           call dblepr("upperHm",-1,uHm,totalu)
           call dblepr("upperHb",-1,uHb,totalu)
           call dblepr("fsmax",-1,fsmax,1)
           call rexit("Negative area in ARS")
         end if
         
         uHpr(i)=uHpr(i)/accum
      end do

      if(accum.eq.0.d0)then
         call rexit("Total area equal to 0 in ARS")      
         do i=1,totalu
            uHpr(i)=1.d0/dble(totalu)
         end do
      end if

      return
      end

c=======================================================================
      subroutine arsUpdateSD(n,maxn,grid,fs,fps,sx,fsx,fpsx,fsmax,err)
c=======================================================================
c     update the support points in a ARS using derivatives
c     A.J.V., 2007
      implicit none

c++++ Input      
      integer maxn,n
      real*8 sx,fsx,fpsx,fsmax
      real*8 grid(maxn),fs(maxn),fps(maxn)

c++++ working space
      integer maxn2
      parameter(maxn2=10000)
      real*8 sxfsx(maxn2,3),ssexp,zero
      integer err,fin,i,j,errord,ind,step
      real*8 xeps,yeps
      parameter(xeps=0.00001,yeps=0.1)

c++++ algorithm
      zero=ssexp(-50.d0)
      
      if(maxn.gt.maxn2)then
        call rexit("Error in 'arsUpdateS': Increase maxn2")
      end if

      if((n+1).gt.maxn)then
        call rexit("Error in 'arsUpdateS' (maxn)")
      end if

      err=0
      do i=1,n
         if(abs(sx-grid(i)).lt.xeps)then
            err=1
         end if
      end do

      if(fsx-fsmax.lt.-9999999999999.d0)then
         err=1
      end if

      if(err.eq.0)then

         if(sx.lt.grid(1))then
            ind=1
            step=1
            go to 100
         end if

         if(sx.gt.grid(n))then
           ind=n+1
           step=2           
           go to 100
         end if
         
         fin=0
         ind=1
         step=2
         do while(fin.eq.0.and.ind.lt.n)
            ind=ind+1
            if(sx.lt.grid(ind))fin=1
         end do

         if(sx.gt.grid(ind).or.sx.lt.grid(ind-1))then
           call rexit("Error in the ordering 1")                  
         end if

100      continue

         if(ind.eq.1)then
            if(abs(fs(1)-fsx).lt.yeps)then
               err=1
               go to 200
            end if
            sxfsx(1,1)=sx
            sxfsx(1,2)=fsx
            sxfsx(1,3)=fpsx
            do j=1,n
               sxfsx(j+1,1)=grid(j)
               sxfsx(j+1,2)=fs(j)
               sxfsx(j+1,3)=fps(j)
            end do
            n=n+1
          else if(ind.eq.(n+1))then
            if(abs(fsx-fs(n)).lt.yeps)then
               err=1
               go to 200
            end if
            do j=1,n
               sxfsx(j,1)=grid(j)
               sxfsx(j,2)=fs(j)
               sxfsx(j,3)=fps(j)
            end do
            sxfsx(n+1,1)=sx
            sxfsx(n+1,2)=fsx
            sxfsx(n+1,3)=fpsx
            n=n+1
          else
            if(abs(fsx-fs(ind-1)).lt.yeps)then
               err=1
               go to 200
            end if
            if(abs(fsx-fs(ind)).lt.yeps)then
               err=1
               go to 200
            end if
            do j=1,ind-1
               sxfsx(j,1)=grid(j)
               sxfsx(j,2)=fs(j)
               sxfsx(j,3)=fps(j)
            end do
            sxfsx(ind,1)=sx
            sxfsx(ind,2)=fsx
            sxfsx(ind,3)=fpsx
            do j=ind,n
               sxfsx(j+1,1)=grid(j)
               sxfsx(j+1,2)=fs(j)
               sxfsx(j+1,3)=fps(j)
            end do
            n=n+1
         end if

         do j=1,n
            grid(j)=sxfsx(j,1)
            fs(j)=sxfsx(j,2)
            fps(j)=sxfsx(j,3)
         end do   

         errord=0
         do i=2,n
            if(grid(i).le.grid(i-1))then
               call intpr("ind",-1,ind,1)
               call intpr("step",-1,step,1)               
               call dblepr("grid",-1,grid,n)
               call dblepr("fs",-1,fs,n)
               call rexit("Error in the ordering 2")         
            end if
         end do 
      end if   

200   continue

      return
      end

c=======================================================================
      subroutine arsSampleHD(x,n,maxu,uHm,uHl,uHr,uHpr)
c=======================================================================
c     generates a sample x from the Upper Hull in an ARS scheme
c     for a log-concave distrtribution using derivatives
c     A.J.V., 2007
      implicit none

c++++ dimensions
      integer n 
      integer maxu
      
c++++ upper Hull
      real*8 uHm(maxu)
      real*8 uHl(maxu),uHr(maxu)
      real*8 uHpr(maxu)

c++++ output
      real*8 x

c++++ working space
      integer i,ind,totalu,ok 
      real*8 dh,emax,u,left,m,right,ssexp
      real*8 total,tmp1
      real*8 zero
      real runif
      parameter(emax=50.d0)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zero=ssexp(-emax)
      
      totalu=2+2*(n-1) 

      if(totalu.gt.maxu)then
        call rexit("Error in 'arsSampleH' (maxu)")
      end if

c++++ randomly choose a line segment
      u=dble(runif())

      total=0.d0
      do i=1,totalu
         total=total+uHpr(i)
      end do

      if(total.eq.0.d0)then
        call rexit("Probabilities add to 0 in 'arsSampleH'") 
      end if  

      i=1
      ok=1
      ind=0
      tmp1=0.0
      do while(ok.eq.1.and.i.le.totalu)
         tmp1=tmp1+uHpr(i)
         if(u.lt.tmp1)then
            ok=0 
            ind=i
         end if
         i=i+1
      end do
      
      if(ind.lt.1.or.ind.gt.totalu)then
          call intpr("selected line",-1,ind,1)
          call dblepr("Prob",-1,uHpr,totalu)
          call dblepr("Slopes",-1,uHm,totalu)
          call dblepr("left",-1,uHl,totalu)
          call dblepr("right",-1,uHr,totalu)
          call rexit("Error in ARS generation: line outside range")
      end if

c++++ sample along that (ind) line segment

      u=dble(runif())

      m=uHm(ind)
      left=uHl(ind)
      right=uHr(ind)

      if(right.eq.left)then
        x=right
        go to 100
      end if
   

      if(ind.eq.1)then
         x=right+log(u)/m
        else if(ind.eq.totalu)then
         x=left+log(u)/m
        else
         dh=m*(left-right)
         if(abs(m).lt.zero)then
            x=left+u*(right-left)   
          else if(dh.lt.emax)then
              x=right+log(u*(1-ssexp(dh))+ssexp(dh))/m
            else
              x=left+log(u)/m
         end if
      end if   
   
100   continue
      
      if(ind.eq.1)then
         if(x.gt.uHr(ind))then
            call dblepr("x",-1,x,1)
            call dblepr("lsup",-1,uHr(ind),1)
            call intpr("selected line",-1,ind,1)
            call dblepr("Prob",-1,uHpr,totalu)
            call dblepr("Slopes",-1,uHm,totalu)
            call dblepr("left",-1,uHl,totalu)
            call dblepr("right",-1,uHr,totalu)
            call rexit("Error in ARS generation: outside limits")      
         end if
      else if(ind.eq.totalu)then   
         if(x.lt.uHl(ind))then
            call dblepr("x",-1,x,1)
            call dblepr("linf",-1,uHl(ind),1)
            call intpr("selected line",-1,ind,1)
            call dblepr("Prob",-1,uHpr,totalu)
            call dblepr("Slopes",-1,uHm,totalu)
            call dblepr("left",-1,uHl,totalu)
            call dblepr("right",-1,uHr,totalu)
            call rexit("Error in ARS generation: outside limits")      
         end if
      else   
         if(x.lt.uHl(ind).or.x.gt.uHr(ind))then
            call dblepr("x",-1,x,1)
            call dblepr("linf",-1,uHl(ind),1)
            call dblepr("lsup",-1,uHr(ind),1)
            call intpr("selected line",-1,ind,1)
            call dblepr("Prob",-1,uHpr,totalu)
            call dblepr("Slopes",-1,uHm,totalu)
            call dblepr("left",-1,uHl,totalu)
            call dblepr("right",-1,uHr,totalu)
            call rexit("Error in ARS generation: outside limits")      
         end if
      end if   
      return
      end

c=======================================================================
      subroutine arsEvalHD(x,n,hlower,hupper,maxn,maxl,maxu,
     &                     lHb,lHm,lHl,lHr,
     &                     uHb,uHm,uHl,uHr)
c=======================================================================
c     Eval the Hull at x in an ARS scheme.
c     A.J.V., 2007

      implicit none

c++++ dimensions
      integer n
      integer maxn,maxl,maxu
      
c++++ lower Hull
      real*8 lHb(maxl),lHm(maxl)
      real*8 lHl(maxl),lHr(maxl)

c++++ upper Hull
      real*8 uHb(maxu),uHm(maxu)
      real*8 uHl(maxu),uHr(maxu)

c++++ input/output
      real*8 x,hlower,hupper

c++++ working space
      integer i,fin,totall,totalu
      real*8 left,right

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      totall=n-1
      totalu=2+2*(n-1) 

      if(totalu.gt.maxu)then
        call rexit("Error in 'arsEvalH' (maxu)")
      end if

      if(totall.gt.maxl)then
        call rexit("Error in 'arsEvalH' (maxl)")
      end if

c++++ lower bound
      
      if(x<lHl(1))then
         hlower=-999999999999999999999999999999999.d0
       else if(x>lHr(totall))then   
         hlower=-999999999999999999999999999999999.d0
       else
         i=1
         fin=0
         do while(fin.eq.0.and.i.le.totall)
            left=lHl(i)
            right=lHr(i)
            
            if(x.ge.left.and.x.le.right)then
               fin=1
               hlower=lHm(i)*x+lHb(i)
            end if
            i=i+1
         end do
      end if

c++++ upper bound

      fin=0
      i=1
      do while(fin.eq.0.and.i.le.totalu)
         left=uHl(i)
         right=uHr(i)
         if(x.ge.left.and.x.le.right)then
            fin=1
            hupper=uHm(i)*x+uHb(i)
         end if
         i=i+1
      end do

      return
      end      


c=======================================================================
c=======================================================================
      subroutine testarsD(nsamp,nstart,start,x,evali,seed1,seed2,
     $                    mu,sd)
c=======================================================================
c=======================================================================
c     Test for the ARS sampling algorithm using derivatives
c     A.J.V., 2007
      implicit none
      integer nsamp,nstart,seed1,seed2,evali(nsamp)
      real*8 x(nsamp),start(nstart),mu,sd
      
      integer i
      real*8 tmp1

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Working space - ARS
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer maxn,maxl,maxu,maxeval
      parameter(maxn=2000,maxl=maxn-1,maxu=2+2*(maxn-1),maxeval=10000)
      real*8 lHb(maxl),lHm(maxl)
      real*8 lHl(maxl),lHr(maxl)
      real*8 uHb(maxu),uHm(maxu)
      real*8 uHl(maxu),uHr(maxu)
      real*8 uHpr(maxu)
      real*8 grid(maxn),fs(maxn),fps(maxn)

      integer accept,neval,counter,err
      real*8 hlower,hupper
      real*8 fsx,fpsx,fsmax,luinf
      real*8 quan(3)
      real runif
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++End Working space - ARS
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c++++ set random number generator
      call setall(seed1,seed2)

c++++ algorithm

      neval=nstart
      do i=1,nstart
         grid(i)=start(i)
      end do
      
      do i=1,neval
         fs(i)=-0.5*( (grid(i)-mu)/sd )**2-100000.d0
         fps(i)=-((1/sd)**2)*(grid(i)-mu)
      end do
      
      call arsCompuHunD(neval,maxn,maxl,maxu,
     &                 lHb,lHm,
     &                 lHl,lHr,
     &                 uHb,uHm,
     &                 uHl,uHr,
     &                 uHpr,grid,fs,fps,fsmax)

      do i=1,nsamp
         counter=0 
         accept=0
         
         do while(accept.eq.0.and.counter.lt.maxeval)
            call rchkusr()
            counter=counter+1

            call arsSampleHD(tmp1,neval,maxn,maxu,uHm,uHl,uHr,uHpr)

            call arsEvalHD(tmp1,neval,hlower,hupper,maxn,maxl,maxu,
     &                     lHb,lHm,lHl,lHr,
     &                     uHb,uHm,uHl,uHr,uHpr,grid,fs)
     

            luinf=dlog(dble(runif()))
            hlower=hlower-fsmax
            hupper=hupper-fsmax
            
            if(tmp1.gt.grid(1).and.tmp1.lt.grid(neval))then
                if(luinf.le.(hlower-hupper))accept=1
            end if
            if(accept.eq.0)then
               fsx=-0.5*( (tmp1-mu)/sd )**2-100000.d0
               fpsx=-((1/sd)**2)*(tmp1-mu)
               
               if(luinf.le.(fsx-fsmax-hupper))accept=1
               
               if(accept.eq.0)then
                  if((neval+1).le.maxn)then 
                    call arsUpdateSD(neval,maxn,grid,fs,fps,tmp1,
     &                               fsx,fpsx,fsmax,err)

                    if(err.eq.0)then

                       call arsCompuHunD(neval,maxn,maxl,maxu,
     &                           lHb,lHm,lHl,lHr,uHb,uHm,uHl,uHr,
     &                           uHpr,grid,fs,fps,fsmax)
                    end if 
                  end if
               end if
            end if   

            if(accept.eq.0.and.counter.eq.maxeval)then
              call rexit("Maximum # of evaluations reached")
            end if
            
         end do 
         evali(i)=counter
         x(i)=tmp1
         call arsQuanH(neval,quan,maxn,maxu,uHm,uHl,uHr,uHpr)
         
      end do
      
      return
      end



c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c    DERIVARIVE-FREE FUNCTIONS (STILL UNDER DEVELOPMENT)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c=======================================================================
      subroutine arsCompuHunDF(n,maxn,maxl,maxu,lHb,lHm,lHl,lHr,
     &                       uHb,uHm,uHl,uHr,uHpr,grid,fs,fsmax)
c=======================================================================
c     computes the picewise-linear hull in a derivative free
c     ARS scheme for a distribution with unbounded support 
c     (DON'T USE IT)
c     A.J.V., 2007

      implicit none

c++++ dimensions
      integer n
      integer maxn,maxl,maxu
      
c++++ lower Hull
      real*8 lHb(maxl),lHm(maxl)
      real*8 lHl(maxl),lHr(maxl)

c++++ upper Hull
      real*8 uHb(maxu),uHm(maxu)
      real*8 uHl(maxu),uHr(maxu)
      real*8 uHpr(maxu)

c++++ evaluations
      real*8 grid(maxn),fs(maxn)

c++++ working space
      integer i,j,totall,totalu,maxcheck
      parameter(maxcheck=2000)
      integer checki(maxcheck),totali,counter
      real*8 b,m,pr
      real*8 b1,b2,m1,m2
      real*8 ix,hix,pr1,pr2
      real*8 tmp1,tmp2,tmp3,ssexp
      real*8 part1,part2,part3
      real*8 zero,fsmax,emax
      parameter(emax=50.d0)
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zero=ssexp(-emax)

      if(maxn.gt.maxcheck)then
        call rexit("Error in 'arsCompuHnu' (maxcheck)")
      end if

c++++ checking -Inf before start
      totali=0
      do i=1,n
         checki(i)=0
         if(fs(i).lt.-9999999999999999999999.d0)then
            checki(i)=1
            totali=totali+1
         end if
      end do

      counter=0
      if(totali.gt.0)then
        do i=1,n
           if(checki(i).eq.0)then
             counter=counter+1
             grid(counter)=grid(i)
             fs(counter)=fs(i)
           end if
        end do
        n=n-totali
      end if

      if(n.lt.3)then
        call rexit("Error in 'arsCompuHnu': number of points < 3")
      end if

c++++ end checking
      
      totall=n-1
      if(n.gt.3)then
        totalu=4+2*(n-3) 
       else
        totalu=4
      end if

      if(totalu.gt.maxu)then
        call rexit("Error in 'arsCompuHnu' (maxu)")
      end if

      if(totall.gt.maxl)then
        call rexit("Error in 'arsCompuHnu' (maxl)")
      end if

      do i=1,maxu
         uHm(i)=0.d0
         uHb(i)=0.d0
         uHl(i)=0.d0
         uHr(i)=0.d0
         uHpr(i)=0.d0
      end do

      do i=1,maxl
         lHb(i)=0.d0
         lHm(i)=0.d0
         lHl(i)=0.d0
         lHr(i)=0.d0
      end do

      fsmax=fs(1)
      do i=2,n
         fsmax=max(fs(i),fsmax)
      end do   

c++++ computing the lines of the Lower Hull

      do i=1,totall
         lHm(i)=(fs(i+1)-fs(i))/(grid(i+1)-grid(i))
         lHb(i)=fs(i)-lHm(i)*grid(i)
         lHl(i)=grid(i)
         lHr(i)=grid(i+1)
      end do

c++++ computing the lines of the Upper Hull

c++++ first line
      m=(fs(2)-fs(1))/(grid(2)-grid(1))
      b=fs(1)-m*grid(1)

      tmp1=b-fsmax
      tmp2=m*grid(1)
      tmp3=m*(-999999.d0)
      part1=ssexp(tmp1)
      part2=ssexp(tmp2)
      part3=ssexp(tmp3)

      if(abs(m).lt.zero)then
         pr=(part1)*(part2+part3)*(grid(1)+999999.d0)/2.d0
        else
         if(m.lt.0.d0)then
            pr=(part1/m)*(part2-part3)
          else
            pr=(part1/m)*part2
         end if
      end if
      
      uHm(1)=m
      uHb(1)=b
      uHpr(1)=pr
      uHl(1)=-999999.d0
      uHr(1)=grid(1)

c++++ second line      
      m=(fs(3)-fs(2))/(grid(3)-grid(2))
      b=fs(2)-m*grid(2)
      tmp1=b-fsmax
      tmp2=m*grid(2)
      tmp3=m*grid(1)

      part1=ssexp(tmp1)
      part2=ssexp(tmp2)
      part3=ssexp(tmp3)

      if(abs(m).lt.zero)then
         pr=(part1/2.d0)*(part2+part3)*(grid(2)-grid(1))
       else
        pr=(part1/m)*(part2-part3)
      end if         

      uHm(2)=m
      uHb(2)=b
      uHpr(2)=pr
      uHl(2)=grid(1)
      uHr(2)=grid(2)

c++++ interior lines 
      j=2
      if(n.gt.3)then
        do i=2,(n-2)
c+++++++++ first interior line
           m1= (fs(i)-fs(i-1))/(grid(i)-grid(i-1))
           b1= fs(i)-m1*grid(i)

c+++++++++ second interior line
           m2= (fs(i+2)-fs(i+1))/(grid(i+2)-grid(i+1))
           b2= fs(i+1)-m2*grid(i+1)

c+++++++++ intersection
           if(abs(m2-m1).lt.zero)then
             ix=0.5*(grid(i)+grid(i+1))
             hix=0.5*(fs(i)+fs(i+1))
            else if(abs(m1).lt.abs(m2))then
             ix=grid(i+1)+
     &          (fs(i)-fs(i+1)+m1*(grid(i+1)-grid(i)))/(m2-m1)
             hix=m1*(ix-grid(i))+fs(i)
            else            
             ix=grid(i)+
     &          (fs(i)-fs(i+1)+m2*(grid(i+1)-grid(i)))/(m2-m1)
             hix=m2*(ix-grid(i+1))+fs(i+1)
           end if  

           if(ix.lt.grid(i).or.ix.gt.grid(i+1))then
              call dblepr("ix",-1,ix,1)
              call dblepr("xi",-1,grid(i),1)
              call dblepr("xi+1",-1,grid(i+1),1)
              call rexit("Error in intersection")
c              ix=0.5*(grid(i)+grid(i+1))
           end if

c+++++++++ areas

           tmp1=b1-fsmax
           tmp2=m1*ix
           tmp3=m1*grid(i)
           part1=ssexp(tmp1)
           part2=ssexp(tmp2)
           part3=ssexp(tmp3)

           if(abs(m1).lt.zero)then
              pr1=(part1)*(part2+part3)*(ix-grid(i))/2.d0
             else 
              pr1=(part1/m1)*(part2-part3)
           end if             

           tmp1=b2-fsmax
           tmp2=m2*grid(i+1)
           tmp3=m2*ix
           part1=ssexp(tmp1)
           part2=ssexp(tmp2)
           part3=ssexp(tmp3)

           if(abs(m2).lt.zero)then
              pr2=(part1)*(part2+part3)*(grid(i+1)-ix)/2.d0
             else 
              pr2=(part1/m2)*(part2-part3)
           end if    

c+++++++++ update first interior line           
           j=j+1
           uHm(j)=m1 
           uHb(j)=b1
           uHpr(j)=pr1
           uHl(j)=grid(i)
           uHr(j)=ix

c+++++++++ update second interior line           
           j=j+1
           uHm(j)=m2 
           uHb(j)=b2
           uHpr(j)=pr2
           uHl(j)=ix
           uHr(j)=grid(i+1)
        end do
      end if

c++++ second last line
      m=(fs(n-1)-fs(n-2))/(grid(n-1)-grid(n-2))
      b=fs(n-1)-m*grid(n-1)
      tmp1=b-fsmax
      tmp2=m*grid(n)
      tmp3=m*grid(n-1)
      part1=ssexp(tmp1)
      part2=ssexp(tmp2)
      part3=ssexp(tmp3)

      if(abs(m).lt.zero)then
         pr=(part1)*(part2+part3)*(grid(n)-grid(n-1))/2.d0
        else 
         pr=(part1/m)*(part2-part3)
      end if
      
      j=j+1
      uHm(j)=m
      uHb(j)=b
      uHpr(j)=pr
      uHl(j)=grid(n-1)
      uHr(j)=grid(n)

c++++ last line
      m=(fs(n)-fs(n-1))/(grid(n)-grid(n-1))
      b=fs(n)-m*grid(n)
      tmp1=b-fsmax
      tmp2=m*(999999.d0)
      tmp3=m*grid(n)
      part1=ssexp(tmp1)
      part2=ssexp(tmp2)
      part3=ssexp(tmp3)

      if(abs(m).lt.zero)then
         pr=(part1)*(part2+part3)*(999999.d0-grid(n))/2.d0
        else
         if(m.gt.0.d0)then
            pr=(part1/m)*(part2-part3)
          else
            pr=(part1/m)*(-part3)
         end if
      end if
      
      j=j+1
      uHm(j)=m
      uHb(j)=b
      uHpr(j)=pr
      uHl(j)=grid(n)
      uHr(j)=+999999.d0

c++++ standardizing the areas
      tmp1=0.d0
      do i=1,totalu
         if(uHpr(i).lt.0.d0)then
           call dblepr("s",-1,grid,n)
           call dblepr("fs",-1,fs,n)
           call dblepr("upperHpr",-1,uHpr,totalu)
           call dblepr("upperHm",-1,uHm,totalu)
           call rexit("Negative area in ARS")
         end if
         tmp1=tmp1+uHpr(i)
      end do

      do i=1,totalu
         uHpr(i)=uHpr(i)/tmp1
      end do

      if(tmp1.eq.0.d0)then
         call rexit("Total area equal to 0 in ARS")      
         do i=1,totalu
            uHpr(i)=1.d0/dble(totalu)
         end do
      end if

      return
      end

c=======================================================================
      subroutine arsSampleHDF(x,n,maxu,uHm,uHl,uHr,uHpr)
c=======================================================================
c     generates a sample x from the Upper Hull in a derivative free 
c     ARS scheme for a log-concave distrtribution (DON'T USE IT)
c     A.J.V., 2007
      implicit none

c++++ dimensions
      integer n 
      integer maxu
      
c++++ upper Hull
      real*8 uHm(maxu)
      real*8 uHl(maxu),uHr(maxu)
      real*8 uHpr(maxu)

c++++ output
      real*8 x

c++++ working space
      integer i,ind,totalu,ok 
      real*8 e1,e2,emax,u,left,m,right,ssexp
      real*8 part1,part2,total,tmp1
      real*8 zero
      real runif
      parameter(emax=50.d0)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      zero=ssexp(-emax)
      
      if(n.ge.3)then
        totalu=4+2*(n-3) 
       else
        totalu=4
      end if

      if(totalu.gt.maxu)then
        call rexit("Error in 'arsSampleH' (maxu)")
      end if

c++++ randomly choose a line segment
      u=dble(runif())

      total=0.d0
      do i=1,totalu
         total=total+uHpr(i)
      end do

      if(total.eq.0.d0)then
        call rexit("Probabilities add to 0 in 'arsSampleH'") 
      end if  

      i=1
      ok=1
      ind=0
      tmp1=0.0
      do while(ok.eq.1.and.i.le.totalu)
         tmp1=tmp1+uHpr(i)
         if(u.lt.tmp1)then
            ok=0 
            ind=i
         end if
         i=i+1
      end do
      
      if(ind.lt.1.or.ind.gt.totalu)then
          call intpr("selected line",-1,ind,1)
          call dblepr("Prob",-1,uHpr,totalu)
          call dblepr("Slopes",-1,uHm,totalu)
          call dblepr("left",-1,uHl,totalu)
          call dblepr("right",-1,uHr,totalu)
          call rexit("Error in ARS generation: line outside range")
      end if

c++++ sample along that (ind) line segment

      u=dble(runif())

      m=uHm(ind)
      left=uHl(ind)
      right=uHr(ind)

      if(right.eq.left)then
        x=right
        go to 100
      end if
   
      e1=m*right
      e2=m*left
      part1=ssexp(e1)
      part2=ssexp(e2)

      if(e1.gt.emax.and.e2.gt.emax)then
         x=left 
        else if(abs(m).lt.zero)then
         x=left+u*(right-left)   
        else if(part1.eq.part2)then
         x=left 
        else
         x=log(part2+u*(part1-part2))/m
      end if   
   
100   continue

      if(x.lt.uHl(ind).or.x.gt.uHr(ind))then
         call dblepr("x",-1,x,1)
         call dblepr("linf",-1,uHl(ind),1)
         call dblepr("lsup",-1,uHr(ind),1)
         call intpr("selected line",-1,ind,1)
         call dblepr("Prob",-1,uHpr,totalu)
         call dblepr("Slopes",-1,uHm,totalu)
         call dblepr("left",-1,uHl,totalu)
         call dblepr("right",-1,uHr,totalu)
         call rexit("Error in ARS generation: outside limits")
      end if
      return
      end


c=======================================================================
      subroutine arsUpdateSDF(n,maxn,grid,fs,sx,fsx,err)
c=======================================================================
c     update the support points in a derivative free ARS scheme
c     A.J.V., 2007
      implicit none

c++++ Input      
      integer maxn,n
      real*8 sx,fsx
      real*8 grid(maxn),fs(maxn)

c++++ working space
      integer maxn2
      parameter(maxn2=1000)
      real*8 sxfsx(maxn2,2),ssexp,zero
      integer err,fin,i,j,errord,ind,step
      real*8 xeps,yeps
      parameter(xeps=0.00001,yeps=0.1)

c++++ algorithm
      zero=ssexp(-50.d0)
      
      if(maxn.gt.maxn2)then
        call rexit("Error in 'arsUpdateS': Increase maxn2")
      end if

      if((n+1).gt.maxn)then
        call rexit("Error in 'arsUpdateS' (maxn)")
      end if

      err=0
      do i=1,n
         if(abs(sx-grid(i)).lt.xeps)then
            err=1
         end if
      end do

      if(fsx.lt.-9999999999999999999999.d0)then
         err=1
      end if

      if(err.eq.0)then

         if(sx.lt.grid(1))then
            ind=1
            step=1
            go to 100
         end if

         if(sx.gt.grid(n))then
           ind=n+1
           step=2           
           go to 100
         end if
         
         fin=0
         ind=1
         step=2
         do while(fin.eq.0.and.ind.lt.n)
            ind=ind+1
            if(sx.lt.grid(ind))fin=1
         end do

         if(sx.gt.grid(ind).or.sx.lt.grid(ind-1))then
           call rexit("Error in the ordering 1")                  
         end if

100      continue

         if(ind.eq.1)then
            if(abs(fs(1)-fsx).lt.yeps)then
               err=1
               go to 200
            end if
            sxfsx(1,1)=sx
            sxfsx(1,2)=fsx
            do j=1,n
               sxfsx(j+1,1)=grid(j)
               sxfsx(j+1,2)=fs(j)
            end do
            n=n+1
          else if(ind.eq.(n+1))then
            if(abs(fsx-fs(n)).lt.yeps)then
               err=1
               go to 200
            end if
            do j=1,n
               sxfsx(j,1)=grid(j)
               sxfsx(j,2)=fs(j)
            end do
            sxfsx(n+1,1)=sx
            sxfsx(n+1,2)=fsx
            n=n+1
          else
            if(abs(fsx-fs(ind-1)).lt.yeps)then
               err=1
               go to 200
            end if
            if(abs(fsx-fs(ind)).lt.yeps)then
               err=1
               go to 200
            end if
            do j=1,ind-1
               sxfsx(j,1)=grid(j)
               sxfsx(j,2)=fs(j)
            end do
            sxfsx(ind,1)=sx
            sxfsx(ind,2)=fsx
            do j=ind,n
               sxfsx(j+1,1)=grid(j)
               sxfsx(j+1,2)=fs(j)
            end do
            n=n+1
         end if

         do j=1,n
            grid(j)=sxfsx(j,1)
            fs(j)=sxfsx(j,2)
         end do   

         errord=0
         do i=2,n
            if(grid(i).le.grid(i-1))then
               call intpr("ind",-1,ind,1)
               call intpr("step",-1,step,1)               
               call dblepr("grid",-1,grid,n)
               call dblepr("fs",-1,fs,n)
               call rexit("Error in the ordering 2")         
            end if
         end do 
      end if   

200   continue

      return
      end

c=======================================================================
      subroutine arsQuanH(n,quan,maxu,uHm,uHl,uHr,uHpr)
c=======================================================================
c     compute quantiles from the Upper Hull in a derivative free 
c     ARS scheme for a log-concave distrtribution
c     A.J.V., 2007
      implicit none

c++++ dimensions
      integer n 
      integer maxu
      
c++++ upper Hull
      real*8 uHm(maxu)
      real*8 uHl(maxu),uHr(maxu)
      real*8 uHpr(maxu)

c++++ output
      real*8 quan(3)

c++++ working space
      integer i,i1,i2,i3,ind,totalu 
      real*8 e1,e2,u,left,m,right,ssexp
      real*8 ok,part1,part2
      real*8 tmp1
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(n.ge.3)then
        totalu=4+2*(n-3) 
       else
        totalu=4
      end if

      if(totalu.gt.maxu)then
        call rexit("Error in 'arsQuanHun' (maxu)")
      end if

c++++ look for the line segment

      i=1
      ok=1
      ind=0
      tmp1=0.d0
      do while(ok.eq.1.and.i.le.totalu)
         tmp1=tmp1+uHpr(i)
         if(0.05.lt.tmp1)then
            ok=0 
            ind=i
         end if
         i=i+1
      end do
      i1=ind

      i=1
      ok=1
      ind=0
      tmp1=0.d0
      do while(ok.eq.1.and.i.le.totalu)
         tmp1=tmp1+uHpr(i)
         if(0.5.lt.tmp1)then
            ok=0 
            ind=i
         end if
         i=i+1
      end do
      i2=ind

      i=1
      ok=1
      ind=0
      tmp1=0.d0
      do while(ok.eq.1.and.i.le.totalu)
         tmp1=tmp1+uHpr(i)
         if(0.95.lt.tmp1)then
            ok=0 
            ind=i
         end if
         i=i+1
      end do
      i3=ind
      
c++++ look for q1 along that (i1) line segment

      m=uHm(i1)
      left=uHl(i1)
      right=uHr(i1)

      e1=m*right
      e2=m*left
      part1=ssexp(e1)
      part2=ssexp(e2)

      tmp1=0.d0
      if(i1.gt.1)then
         do i=1,i1-1
            tmp1=tmp1+uHpr(i)
         end do
      end if
      
      u=0.05-tmp1
      quan(1)=log(u*(part1-part2)+part2)/m
 
c++++ look for q2 along that (i2) line segment

      m=uHm(i2)
      left=uHl(i2)
      right=uHr(i2)

      e1=m*right
      e2=m*left
      part1=ssexp(e1)
      part2=ssexp(e2)

      tmp1=0.d0
      if(i2.gt.1)then
         do i=1,i2-1
            tmp1=tmp1+uHpr(i)
         end do
      end if
      
      u=0.5-tmp1
      quan(2)=log(u*(part1-part2)+part2)/m

c++++ look for q3 along that (i3) line segment

      m=uHm(i3)
      left=uHl(i3)
      right=uHr(i3)

      e1=m*right
      e2=m*left
      part1=ssexp(e1)
      part2=ssexp(e2)

      tmp1=0.d0
      if(i3.gt.1)then
         do i=1,i3-1
            tmp1=tmp1+uHpr(i)
         end do
      end if
      
      u=0.95-tmp1
      quan(3)=log(u*(part1-part2)+part2)/m

      if(quan(1).lt.uHl(i1).and.quan(1).gt.uHr(i1))then
          call rexit("Error in ARS quantile computation")      
      end if
 
      if(quan(2).lt.uHl(i2).and.quan(2).gt.uHr(i2))then
          call rexit("Error in ARS quantile computation")      
      end if

      if(quan(3).lt.uHl(i3).and.quan(3).gt.uHr(i3))then
          call rexit("Error in ARS quantile computation")      
      end if
 
      return
      end

c=======================================================================
      subroutine arsEvalHDF(x,n,hlower,hupper,maxn,maxl,maxu,
     &                      lHb,lHm,lHl,lHr,
     &                      uHb,uHm,uHl,uHr)
c=======================================================================
c     Eval the Hull at x in an ARS scheme.
c     A.J.V., 2007

      implicit none

c++++ dimensions
      integer n
      integer maxn,maxl,maxu
      
c++++ lower Hull
      real*8 lHb(maxl),lHm(maxl)
      real*8 lHl(maxl),lHr(maxl)

c++++ upper Hull
      real*8 uHb(maxu),uHm(maxu)
      real*8 uHl(maxu),uHr(maxu)

c++++ input/output
      real*8 x,hlower,hupper

c++++ working space
      integer i,fin,totall,totalu
      real*8 left,right

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ algorithm
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      totall=n-1
      if(n.ge.3)then
        totalu=4+2*(n-3) 
       else
        totalu=4
      end if

      if(totalu.gt.maxu)then
        call rexit("Error in 'arsEvalH' (maxu)")
      end if

      if(totall.gt.maxl)then
        call rexit("Error in 'arsEvalH' (maxl)")
      end if

c++++ lower bound
      
      if(x<lHl(1))then
         hlower=-999999999999999999999999999999999.d0
       else if(x>lHr(totall))then   
         hlower=-999999999999999999999999999999999.d0
       else
         i=1
         fin=0
         do while(fin.eq.0.and.i.le.totall)
            left=lHl(i)
            right=lHr(i)
            
            if(x.ge.left.and.x.le.right)then
               fin=1
               hlower=lHm(i)*x+lHb(i)
            end if
            i=i+1
         end do
      end if

c++++ upper bound

      fin=0
      i=1
      do while(fin.eq.0.and.i.le.totalu)
         left=uHl(i)
         right=uHr(i)
         if(x.ge.left.and.x.le.right)then
            fin=1
            hupper=uHm(i)*x+uHb(i)
         end if
         i=i+1
      end do

      return
      end      


c=======================================================================
c=======================================================================
      subroutine testarsDF(nsamp,nstart,start,x,evali,seed1,seed2,
     $                     mu,sd)
c=======================================================================
c=======================================================================
c     Test for the derivative free ARS sampling algorithm
c     A.J.V., 2007
      implicit none
      integer nsamp,nstart,seed1,seed2,evali(nsamp)
      real*8 x(nsamp),start(nstart),mu,sd
      
      integer i
      real*8 tmp1

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Working space - ARS
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer maxn,maxl,maxu,maxeval
      parameter(maxn=500,maxl=maxn-1,maxu=4+2*(maxn-3),maxeval=10000)
      real*8 lHb(maxl),lHm(maxl)
      real*8 lHl(maxl),lHr(maxl)
      real*8 uHb(maxu),uHm(maxu)
      real*8 uHl(maxu),uHr(maxu)
      real*8 uHpr(maxu)
      real*8 grid(maxn),fs(maxn)

      integer accept,neval,counter,err
      real*8 hlower,hupper
      real*8 fsx,fpsx,fsmax,luinf
      real*8 quan(3)
      real runif
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++End Working space - ARS
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c++++ set random number generator
      call setall(seed1,seed2)

c++++ algorithm

      neval=nstart
      do i=1,nstart
         grid(i)=start(i)
      end do
      
      do i=1,neval
         fs(i)=-0.5*( (grid(i)-mu)/sd )**2-100000.d0
      end do
      
      call arsCompuHunDF(neval,maxn,maxl,maxu,
     &                 lHb,lHm,
     &                 lHl,lHr,
     &                 uHb,uHm,
     &                 uHl,uHr,
     &                 uHpr,grid,fs,fsmax)

      do i=1,nsamp
         counter=0 
         accept=0
         
         do while(accept.eq.0.and.counter.lt.maxeval)
            call rchkusr()
            counter=counter+1

            call arsSampleHDF(tmp1,neval,maxn,maxu,uHm,uHl,uHr,uHpr)

            call arsEvalHDF(tmp1,neval,hlower,hupper,maxn,maxl,maxu,
     &                      lHb,lHm,lHl,lHr,
     &                      uHb,uHm,uHl,uHr,uHpr,grid,fs)
     

            luinf=dlog(dble(runif()))
            hlower=hlower-fsmax
            hupper=hupper-fsmax
            
            if(tmp1.gt.grid(1).and.tmp1.lt.grid(neval))then
                if(luinf.le.(hlower-hupper))accept=1
            end if
            if(accept.eq.0)then
               fsx=-0.5*( (tmp1-mu)/sd )**2-100000.d0
               fpsx=-((1/sd)**2)*(tmp1-mu)
               
               if(luinf.le.(fsx-fsmax-hupper))accept=1
               
               if(accept.eq.0)then
                  if((neval+1).le.maxn)then 
                    call arsUpdateSDF(neval,maxn,grid,fs,tmp1,fsx,err)

                    if(err.eq.0)then
                       call arsCompuHunDF(neval,maxn,maxl,maxu,
     &                           lHb,lHm,lHl,lHr,uHb,uHm,uHl,uHr,
     &                           uHpr,grid,fs,fsmax)

                    end if 
                  end if
               end if
            end if   

            if(accept.eq.0.and.counter.eq.maxeval)then
              call rexit("Maximum # of evaluations reached")
            end if
            
         end do 
         evali(i)=counter
         x(i)=tmp1
         call arsQuanH(neval,quan,maxn,maxu,uHm,uHl,uHr,uHpr)
         
      end do
      
      return
      end


