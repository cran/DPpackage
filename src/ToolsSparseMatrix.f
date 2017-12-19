c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR SPARSE MATRIX OPERATIONS
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
      subroutine hashm(y,j1,k1,ind,m,nr)
c=======================================================================      
c     stores or acculuates a spare matrix element y(j1,k1) in a hash
c     table ind, which rank is m x 3. nr returns the number of stored 
c     elements
c
c     A.J.V. 2007
c=======================================================================
      integer j1,k1,m,nr
      double precision y,ind(m,3)
      integer iaddr
      integer ie,ieq,izer,iaddress,k,j
      data ie/641/

      iaddr=433*j1+53*k1
      iaddress=mod(iabs(iaddr),m)+1
 
      do 10 k=1,200
         j=iaddress
         if (ind(iaddress,1).ne.j1.or.ind(iaddress,2).ne.k1)then
              ieq=1
            else
              ieq=0
         endif
         if (ind(iaddress,1).ne.0) then
             izer=1
           else
             izer=0
         end if
 
         if(izer.eq.0.or.ieq.eq.0)then
            if(izer.eq.0)then
               ind(iaddress,1)=j1
               ind(iaddress,2)=k1
               ind(iaddress,3)=y
               nr=nr+1
              else
               ind(iaddress,3)=ind(iaddress,3)+y
            end if
            return
         end if
         iaddress=mod(iaddress+ie-1,m)+1
10    continue
 
      call intpr("nr",-1,nr,1)  
      call intpr("m",-1,m,1)  
      call intpr("j1",-1,j1,1)  
      call intpr("k1",-1,k1,1)  
      call rexit("hash matrix too small,increase m")
      return
      end

c=======================================================================
      subroutine hashiajaa(x,nhash,n,ia,ja,a,m,tmp)
c=======================================================================
c     copies data form hash to the linked list ia-ja-a form. tmp is a 
c     temporary array of size n. Entries are sorted within rows
c
c     A.J.V. 2007
c=======================================================================
      integer nhash,n,m,ia(n+1),ja(m),tmp(n)
      double precision x(nhash,3),a(m)
      integer i,imin,imax,j,k,maxcol,maxrow,size

c++++ count the number of entries in each column of a
      call zero(n,tmp)
      maxrow=0
      maxcol=0
      do i=1,nhash
         j=x(i,1)
         if (j.ne.0) then
            if (j.gt.n) then
                 maxrow=max(maxrow,int(x(i,1)))
                 x(i,1)=0
            elseif (x(i,2).gt.n) then
                 maxcol=max(maxcol,int(x(i,2)))
                 x(i,1)=0
            else
                 tmp(j)=tmp(j)+1
            endif
         endif
      enddo

c++++ create the row count for b
      ia(1)=1
      do i=1,n
         ia(i+1)=ia(i)+tmp(i)
      enddo

      if(ia(n+1)-1.gt.m)then
         call intpr("ia(n+1)",-1,ia(n+1),1)  
         call rexit("Too small parameter m in hashiajaa")
      endif

c++++ load a into b
      call zero(n,tmp)
      do i=1,nhash
            j=x(i,1)
            if (j.ne.0) then
            k=ia(j)+tmp(j)
            ja(k)=x(i,2)
            a(k)=x(i,3)
            tmp(j)=tmp(j)+1
         endif
      enddo

c++++ sort columns
      do i=1,n
         imin=ia(i)
         imax=ia(i+1)-1
         size=imax-imin+1
         if (size.gt.1) then
            call sortjsp(ja(imin),a(imin),size)
         endif
      enddo
      end

c=======================================================================
      subroutine sortjsp(ja,xa,n)
c=======================================================================
c     sorts intger vector ja in ascending order, and orders the double precision 
c     vector a accordingly
c
c     A.J.V. 2007
c=======================================================================
      integer i,j,s,l,r,j1,a1,n,x
      integer ja(n),stack(50,2)
      double precision xa(n),xb

      s=1
      stack(1,1)=1
      stack(1,2)=n
10    l=stack(s,1)
      r=stack(s,2)
      s=s-1
20    i=l
      j=r
      j1=(l+r)/2
      x=ja(j1)

30    if (ja(i).lt.x) then
      i=i+1
      goto 30
      endif

40     if (ja(j).gt.x) then
      j=j-1
      goto 40
      endif
      if (i.le.j) then
      a1=ja(i)
      ja(i)=ja(j)
      ja(j)=a1
      xb=xa(i)
      xa(i)=xa(j)
      xa(j)=xb
      i=i+1
      j=j-1
      endif
      if (i.le.j) goto 30
      if (i.lt.r) then
      s=s+1
      stack(s,1)=i
      stack(s,2)=r
      endif
      r=j
      if (l.lt.r) goto 20
      if (s.ne.0) goto 10
      end

c=======================================================================
      subroutine quadrf(qaq,a,ia,ja,q,ibeg,iend,ivbeg)
c=======================================================================      
c     calculates the quadratic form qAq of part of vector q and part of 
c     sparse matrix  a (a,ia,ja). 
c     Quadratic form is taken from row-column
c     ibeg to iend.  Vector q is used from ivbeg+1 location
c
c     A.J.V. 2007
c=======================================================================
      integer ia(1),ja(1),ibeg,iend,ivbeg,i,j
      double precision qaq
      double precision a(1),q(1)
      
      qaq=0
      do 10 i=ibeg,iend
        do 10 j=ia(i),ia(i+1)-1
10       if (ja(j).ge.ibeg .and. ja(j).le.iend) 
     $       qaq=qaq+a(j)*q(i+ivbeg)*q(ja(j)+ivbeg)
      return
      end

c=======================================================================
      subroutine adjadd(a,ia,ja,i1,i2,b,ib,jb,j1,mult,stora,storb)
c=======================================================================
c     to sparse matrix a (a,ia,ja), rows j1:j2 add sparse matrix 
c     b (b,ib,jb) beginning from row j1, multiplied by mult.It is 
c     assumed that elements of a that are being added already exist. 
c
c     A.J.V. 2007
c=======================================================================
      integer storhalf,storfull
      parameter(storfull=1,storhalf=2)
      integer ia(1),ja(1),i1,i2,ib(1),jb(1),j1,stora,storb,i,j,rowb,
     $ rowb1 
      double precision a(1),b(1),mult
      
      if (stora.ne.storhalf.or.storb.ne.storfull) then
         call rexit("adjadd: option not implemented")
      end if
            
      do 10 i=i1,i2
        rowb=ib(i+j1-i1)
        rowb1=ib(i+j1-i1+1)
        do 10 j=ia(i),ia(i+1)-1
20        if (rowb.lt.rowb1) then
              if (jb(rowb)+i1-j1.eq.ja(j))  then
                  a(j)=a(j)+mult*b(rowb)
                  rowb=rowb+1
                else
                  rowb=rowb+1
                  goto 20
              endif      
          endif    
10    continue            
      return
      end

c=======================================================================      
      subroutine addspar(a,ia,ja,k,l,x)
c=======================================================================
c     adds x to (k,l) element of sparse matrix (a,ia,ja). The sparse 
c     matrix must already contain (k,l) location      
c
c     A.J.V. 2007
c=======================================================================
      double precision a(1),x
      integer ia(1),ja(1),k,l,i
      
      do 10 i=ia(k),ia(k+1)-1
        if (ja(i).eq.l) then
            a(i)=a(i)+x
            goto 20
        endif
10    continue

      call intpr("k",-1,k,1)  
      call intpr("l",-1,l,1)  
      call rexit("addspar:indices not found")
20    return
      end      

c=======================================================================      
      subroutine prhash(h,id,nhash)
c=======================================================================      
c     prints matrix in the hash table, must be <13
c
c     A.J.V. 2007
c=======================================================================
      integer nhash,i,j,isize,id
      double precision h(id,3),x(12,12)

      isize=0
      do 10 i=1,12
        do 10 j=1,12
10        x(i,j)=0
      do 15 i=1,nhash
      if (h(i,1).ne.0) then
            if(h(i,1).le.12 .and. h(i,2).le.12) then
               isize=max(int(h(i,1)),int(h(i,2)),isize)
               x(int(h(i,1)),int(h(i,2)))=h(i,3)
             else
               do j=1,3
                  call dblepr("prhash:elements too large>",-1,x(i,j),1) 
               end do  
            end if
      end if
15    continue
      do i=1,isize
         do j=1,isize
            call dblepr("x(i,j)",-1,x(i,j),1) 
         end do  
      end do
 
      return
      end      
            
c=======================================================================      
      subroutine prset(h,nh,isize)
c=======================================================================      
c     prints a small sparse matrix of rank isize<=12 stored in hash 
c     table h(nh,3)
c
c     A.J.V. 2007
c=======================================================================
      integer nh,isize
      double precision x(12,12),h(nh,3)
      integer i,j,ii,jj
      
      if(isize.gt.12)then
         call rexit("matrix too large to print")
      end if
      do 10 i=1,isize
        do 10 j=1,isize
10        x(i,j)=0
      do 20 i=1,isize
         ii=int(h(i,1))
         jj=int(h(i+1,1))   
         do 20 j=ii,jj-1
          if(h(j,2).gt.isize) then
             call intpr("element outside bound",-1,i,1) 
             call dblepr("element outside bound",-1,h(j,2),1) 
             call dblepr("element outside bound",-1,h(j,3),1) 
            else
            x(i,int(h(j,2)))=h(j,3)
          end if
20    continue              

      do i=1,isize
         do j=1,isize
            call dblepr("x(i,j)",-1,x(i,j),1) 
         end do  
      end do
      return
      end      


c=======================================================================      
      subroutine prspar(a,ia,ja,isize)
c=======================================================================      
c     prints a small sparse  matrix (a,ia,ja) of rank isize<=12
c
c     A.J.V. 2007
c=======================================================================
      double precision a(*),x(12,12)
      integer ia(*),ja(*),isize,i,j

      if(isize.gt.12) then
         call rexit("matrix too large to print")
      end if
      
      do 10 i=1,isize
        do 10 j=1,isize
10        x(i,j)=0
      do 20 i=1,isize
        do 20 j=ia(i),ia(i+1)-1
           if(ja(j).gt.isize) then
              call intpr("element outside bound",-1,i,1) 
              call intpr("element outside bound",-1,ja(j),1) 
              call dblepr("element outside bound",-1,a(j),1) 
             else
             if (ja(j).ne.0) x(i,ja(j))=a(j)
           end if
20    continue

      do i=1,isize
         do j=1,isize
            call dblepr("x(i,j)",-1,x(i,j),1) 
         end do  
      end do
30    continue

      return
      end

c=======================================================================      
      subroutine prsparu(a,ia,ju,iju,isize)
c=======================================================================      
c     prints a small sparse  matrix (a,ia,ja) of rank isize<=12
c
c     A.J.V. 2007
c=======================================================================
      double precision a(*),x(12,12)
      integer ia(*),ju(*),iju(*),isize,i,j

      if(isize.gt.12) then
         call rexit("matrix too large to print")      
      end if
      
      do 10 i=1,isize
        do 10 j=1,isize
10        x(i,j)=0
      do 20 i=1,isize
        do 20 j=ia(i),ia(i+1)-1
           if(ju(iju(i)-ia(i)+j).gt.isize) then
             call intpr("element outside bound",-1,i,1)
             call intpr("element outside bound",-1,ju(iju(i)-ia(i)+j),1)
             call dblepr("element outside bound",-1,a(j),1)
            else
              if (ju(iju(i)-ia(i)+j).ne.0) x(i,ju(iju(i)-ia(i)+j))=a(j)
          endif
20    continue

      do i=1,isize
         do j=1,isize
            call dblepr("x(i,j)",-1,x(i,j),1) 
         end do  
      end do

      return
      end

c=======================================================================      
      subroutine prindspar(ia,ja,isize)
c=======================================================================      
c     prints indices of a small sparse  matrix (a,ia,ja) of 
c     rank isize<=12
c
c     A.J.V. 2007
c=======================================================================
      integer ia(*),ja(*),isize,i,j,x(12,12)

      if(isize.gt.12) then
         call rexit("matrix too large to print")
      end if
      
      do 10 i=1,isize
        do 10 j=1,isize
10        x(i,j)=0
      do 20 i=1,isize
        do 20 j=ia(i),ia(i+1)-1
           if(ja(j).gt.isize) then
              call intpr("element outside bound",-1,i,1) 
              call intpr("element outside bound",-1,ja(j),1) 
            else
              if(ja(j).ne.0) x(i,ja(j))=x(i,ja(j))+1
          endif
20    continue

      do i=1,isize
         do j=1,isize
            call intpr("x(i,j)",-1,x(i,j),1) 
         end do  
      end do

      return
      end

c=======================================================================      
      subroutine checkdiag(n,ia,ja,a)
c=======================================================================      
c     quits if sparse matrix ia-ja-a of rank n has missing or <=0 
c     diagonals
c
c     A.J.V. 2007
c=======================================================================
      integer n,ia(*),ja(*),i,j,inddiag
      double precision a(*)

      do i=1,n
         inddiag=0
         if(ia(i+1)-ia(i).le.0) then
            call rexit("equation empty")         
         end if
         do j=ia(i), ia(i+1)-1
            if (ja(j).eq.i) inddiag=j
         end do
         if(inddiag.eq.0) then
            call rexit("missing diagonal element")                  
         else if(a(inddiag).le.0) then
            call rexit("<=0 diagonal element")
         end if
      end do
      end

c=======================================================================      
      function tracediag(n,ia,ja,a,l1,l2)
c=======================================================================      
c     computes the sum of diagonal elements l1-l2 in sparse matrix 
c     ia-ja-a.
c
c     A.J.V. 2007
c=======================================================================
      integer n,ia(*),ja(*),l1,l2,i,j,msglev
      double precision tracediag,a(1)
      data msglev/0/

      tracediag=0
      if(l1.le.0 .or. l2.gt.n) then
         call rexit("l1 or l2 out of range")
      end if
      do i=l1,l2
         do j=ia(i),ia(i+1)-1
            if(ja(j).eq.i) then
               tracediag=tracediag+a(j)
               goto 10
            end if
         end do
         call rexit("tracediag: missing diagonal element")

10       continue
      end do
      if (msglev.ge.1) then
          call intpr("TRACEDIAG: input matrix",-1,1,1) 
          call prspar(a,ia,ja,n)
          call intpr("racediag,l1",-1,l1,1) 
          call intpr("racediag,l2",-1,l2,1) 
      end if
      end

c=======================================================================      
      function tracematrix(n,ia,ja,a,n1,ib,jb,b,ix,iy,tmp)
c=======================================================================      
c     calculates trace of sparse matrices a-ia-ja of order
c     n and ib-jb-b of order n1. Columns of ib-jb-b over n1 are ignored.
c     ix and iy specify offsets for the second matrix.
c     tmp is a real(*8) matrix of order n.
c
c     The matrices may be full or half stored. In the later case,
c     the trace is computed as a sum of trace due to diagonal elements
c     plus 2 * sum of trace due to off-diagonals
c
c
c     A.J.V. 2007
c=======================================================================
      integer n,ia(*),ja(*),n1,ib(*),jb(*),ix,iy,i,j,msglev
      double precision tracematrix,a(*),b(*),tmp(*),tl,td,tu
      data msglev/0/

      do i=1,n
         tmp(i)=0
      enddo

      tl=0
      td=0
      tu=0
      tracematrix=0.d0

c scatter row i of b
      do i=1,n1
         do j=ib(i),ib(i+1)-1
            if (jb(j).le.n1) then
               tmp(jb(j)+iy)=b(j)
            endif
         enddo
c actual trace computation
         do j=ia(i+ix),ia(i+ix+1)-1
            if (ja(j).gt.i+ix) then
               tl=tl+a(j)*tmp(ja(j))
              elseif (ja(j).eq.i+ix) then
               td=td+a(j)*tmp(ja(j))
              else
               tu=tu+a(j)*tmp(ja(j))
             endif
         enddo

c zero row in tmp
         do j=ib(i+ix),ib(i+ix+1)-1
            if (jb(j).le.n1) then
               tmp(jb(j)+iy)=0
            endif
         enddo
      enddo

c Matrices full stored if upper and lower traces equal
      if (tu.eq.0) then
          tracematrix=td+2*tl
         elseif (tl.eq.0) then
          tracematrix=td+2*tu
         else if (abs((tl-tu)/(tl+tu)).lt.1e-6) then
          tracematrix=td+2*tu
         else
          call dblepr("tl",-1,tl,1)          
          call dblepr("td",-1,td,1)          
          call dblepr("tu",-1,tu,1)                    
          call rexit("tracematrix")                                    
       endif
      if (msglev.ge.1) then
          call intpr("TRACEDIAG: input matrix A",-1,1,1) 
          call dblepr("tl",-1,tl,1)          
          call dblepr("td",-1,td,1)          
          call dblepr("tu",-1,tu,1)                    
          call prspar(a,ia,ja,n)
          call intpr("TRACEDIAG: input matrix B",-1,1,1) 
          call prspar(b,ib,jb,n1)
          call dblepr("tracematrix",-1,tracematrix,1)          
          call intpr("ix",-1,ix,1)          
          call intpr("iy",-1,iy,1)                    
      end if
      end

c=======================================================================
      subroutine zero(n,x)
c=======================================================================
c     subroutine that setting to zero the elements in the vector x
c     A.J.V. 2007
c=======================================================================
      integer n,x(n),i
      do i=1,n
         x(i)=0
      end do
      end

