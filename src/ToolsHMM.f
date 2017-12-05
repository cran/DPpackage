c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR HIDDEN MARKOV MODELS
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
      subroutine hmm_sim1(maxm,n,maxnsi,alpha,linf,lsup,u,v,w,workv,
     &                    theta,y)         
c=======================================================================
c     Subroutine to gererate n samples from a DP(alpha U(linf,lsup) 
c     prior in a Semiparametric HMM.
c
c     A.J.V. 2008 
c=======================================================================
      implicit none

c+++++Input
      integer maxm,n,maxnsi
      double precision alpha,linf,lsup

c+++++External working space
      double precision theta(maxm) 
      double precision v(maxm),w(maxm),workv(maxm+1)
      double precision u(maxnsi)

c+++++Output
      double precision y(maxnsi)

c+++++Internal working space
      integer i,i1,mrand,ok
      double precision mass,maxu
      double precision tmp1,tmp2
      real runif 

      if(lsup.le.linf)then
         call rexit("Errors in limits in 'hmm_sim1'")
      end if

      mass=alpha*(lsup-linf)

      maxu=0.d0
      do i=1,n
         u(i)=dble(runif())
         if(u(i).gt.maxu)maxu=u(i)
      end do

      tmp1=0.d0
      tmp2=0.d0

      mrand=10
      call hmm_dpsim(0,maxm,mrand,mass,workv,w,v,tmp2,tmp1)

      do while(maxu.ge.tmp2.and.mrand.lt.maxm)
         mrand=mrand+1
         if(mrand.gt.maxm)then
            call rexit("Increase maxn in 'DPhmm.f'")
         end if  
         call hmm_dpsim(mrand-1,maxm,mrand,mass,workv,w,v,tmp2,tmp1)
      end do

      do i=1,mrand
         theta(i)=linf+(lsup-linf)*dble(runif())
      end do
      
      do i=1,n 
         i1=0
         tmp1=0.d0
         ok=1
         do while(ok.eq.1.and.i1.lt.mrand)
            i1=i1+1 
            tmp1=tmp1+w(i1)
            if(u(i).lt.tmp1)then
               ok=0 
            end if
         end do
         y(i)=theta(i1)
      end do

      return
      end


c=======================================================================
      subroutine hmm_sim2(maxm,n,nso,maxnsi,alpha,linf,lsup,u,v,w,workv,
     &                    theta,yo,y)         
c=======================================================================
c     Subroutine to gererate n samples from a DP(alpha U(linf,lsup) 
c     prior in a Semiparametric HMM.
c
c     A.J.V. 2008 
c=======================================================================
      implicit none

c+++++Input
      integer maxm,n,nso,maxnsi
      double precision alpha,linf,lsup

c+++++External working space
      double precision theta(maxm) 
      double precision v(maxm),w(maxm),workv(maxm+1)
      double precision u(maxnsi)
      double precision yo(maxnsi)

c+++++Output
      double precision y(maxnsi)

c+++++Internal working space
      integer i,i1,mrand,ok
      double precision mass,mass2,maxu
      double precision tmp1,tmp2
      real runif 

      if(lsup.le.linf)then
         call rexit("Errors in limits in 'hmm_sim1'")
      end if

      mass=(alpha+dble(nso))*(lsup-linf)
      mass2=(alpha+dble(nso))

      maxu=0.d0
      do i=1,n
         u(i)=dble(runif())
         if(u(i).gt.maxu)maxu=u(i)
      end do

      tmp1=0.d0
      tmp2=0.d0

      mrand=10
      call hmm_dpsim(0,maxm,mrand,mass,workv,w,v,tmp2,tmp1)

      do while(maxu.ge.tmp2)
         mrand=mrand+1
         if(mrand.gt.maxm)then
            call rexit("Increase maxn in 'DPhmm.f'")
         end if  
         call hmm_dpsim(mrand-1,maxm,mrand,mass,workv,w,v,tmp2,tmp1)
      end do

      do i=1,mrand
         if(dble(runif()).le.(1.d0-(alpha/mass2)))then
            call rdisc(1,nso,i1)
            theta(i)=yo(i1) 
           else
            theta(i)=linf+(lsup-linf)*dble(runif())
         end if
      end do
      
      do i=1,n 
         i1=0
         tmp1=0.d0
         ok=1
         do while(ok.eq.1.and.i1.lt.mrand)
            i1=i1+1 
            tmp1=tmp1+w(i1)
            if(u(i).lt.tmp1)then
               ok=0 
            end if
         end do
         y(i)=theta(i1)
      end do

      return
      end


c=======================================================================
      subroutine hmm_dpsim(indi,nr,n,alpha,workv,w,v,accsum,accsumw)
c=======================================================================
c     Subroutine to compute the stick breaking weights up to
c     finite level n.
c
c     A.J.V. 2007 
c=======================================================================
      implicit none

c+++++Input
      integer indi,nr,n
      double precision alpha

c+++++Output
      double precision accsum,accsumw
      double precision v(nr),w(nr),workv(nr+1)

c+++++Internal working space
      integer i 
      double precision rbeta
      double precision tmp1

c+++++Algorihtm

      workv(1)=0.d0
      do i=(indi+1),n
         v(i)=rbeta(1.0d0,alpha) 
         tmp1=1.d0-v(i)  
         accsumw=accsumw-dlog(tmp1) 
         workv(i+1)=accsumw 
      end do 
      workv(n+1)=1.0D300

      do i=(indi+1),n
         if(i.lt.n)then
            w(i)=dexp(-workv(i))-dexp(-workv(i+1))
          else
            w(i)=dexp(-workv(i)+dlog(v(n)))-dexp(-workv(i+1))
         end if
         accsum=accsum+w(i)
      end do  
      workv(n+1)=accsumw 

      return
      end


c=======================================================================
      subroutine hmm_zassign(nrec,z,maxint,npoints,endp,countv1,countv2,
     &                       countm1,countm2,kk,kkc,zindi,maxkk,probbig,
     &                       maxnsi,countiv,countim)         
c=======================================================================
c     Subroutine to the z's to the corresponding intervals and does the
c     counting in a Semiparametric HMM.
c
c     A.J.V. 2008 
c=======================================================================
      implicit none

c+++++Input
      integer nrec,maxint,maxkk,maxnsi,npoints
      integer kk,kkc
      integer countv1(maxkk),countv2(maxkk)
      integer countm1(maxkk,maxint),countm2(maxkk,maxint) 
      integer countiv(maxkk*maxint),countim(maxkk*maxint,maxnsi)  
      double precision endp(maxint),z(nrec)

c+++++Output
      integer zindi(nrec)
      double precision probbig(maxkk,maxint) 

c+++++Internal working space
      integer i,iceil,j,k,l,ll,ok
      integer k1,k2

c+++++Body

      l=0 
      do i=1,nrec
         zindi(i)=0
         if(kk.gt.kkc)then
            j=iceil(z(i)*kk) 
            k=1
            ok=0
            do while(k.le.countv1(j).and.ok.eq.0) 
               ll=countm1(j,k) 
               if(z(i).le.endp(ll))then
                  ok=1
                  l=ll 
               end if
               k=k+1
            end do
            zindi(i)=l 

          else

            j=iceil(z(i)*kkc) 
            k=1
            ok=0
            do while(k.le.countv2(j).and.ok.eq.0) 
               ll=countm2(j,k) 
               if(z(i).le.endp(ll))then
                  ok=1
                  l=ll 
               end if
               k=k+1
            end do
            zindi(i)=l 

         end if

         if(zindi(i).gt.(npoints+1).or.zindi(i).lt.1)then
            call intpr("zindi(i)",-1,zindi(i),1)
            call rexit("Error in interval assignment")
         end if
 
         if(i.gt.1)then
            k1=iceil(z(i-1)*kk)
            k2=zindi(i)  

            probbig(k1,k2)=probbig(k1,k2)+1.d0

            k=(k1-1)*maxint+k2

            if(k.gt.(maxint*maxkk).or.k.lt.1)then
               call intpr("k",-1,k,1)
               call rexit("Error in interval assignment")
            end if

            countiv(k)=countiv(k)+1

            if(countiv(k).gt.maxnsi)then
               call intpr("maxnsi",-1,maxnsi,1)
               call intpr("countiv(k)",-1,countiv(k),1)      
               call rexit("Increase maxnsi in 'dphmm.f'")
            end if
            countim(k,countiv(k))=i
         end if 
         
      end do
      return
      end
  


c=======================================================================
      subroutine hmm_part(kk,kkc,maxkk,maxint,countv1,countv2,
     &                    countm1,countm2,npoints,endp,intprob)         
c=======================================================================
c     Subroutine to the partition intervals in a Semiparametric HMM.
c
c     A.J.V. 2008 
c=======================================================================
      implicit none

c+++++Input
      integer kk,kkc

c+++++External working space
      integer maxkk,maxint  

c+++++Output
      integer npoints
      double precision endp(maxint)
      integer countv1(maxkk),countv2(maxkk)
      integer countm1(maxkk,maxint),countm2(maxkk,maxint) 
      integer intprob(maxint,2)

c+++++Internal working space
      integer count,count1,count2,countc
      integer i,j,ok,okc
      integer iceil
      double precision tmp1,tmp2

c+++++Body

      do i=1,maxkk
         countv1(i)=0
         countv2(i)=0
      end do

      do i=1,maxint
         endp(i)=0.d0
      end do

      count=0
      count1=1
      count2=1
      tmp1=0.d0
      tmp2=0.d0
      ok=0
 
      if(mod(kk,kkc).eq.0)then

         if(maxint.lt.kk)then            
            call rexit("Error in the maximum number of intervals")
         end if   

         do i=1,kk
            endp(i)=i/dble(kk)

            j=iceil(endp(i)*kk)
            intprob(i,1)=j
            countv1(j)=countv1(j)+1
            countm1(j,countv1(j))=i

            j=iceil(endp(i)*kkc)
            intprob(i,2)=j
            countv2(j)=countv2(j)+1
            countm2(j,countv2(j))=i

         end do
         npoints=kk-1

       else if(mod(kkc,kk).eq.0)then

         if(maxint.lt.kkc)then
            call rexit("Error in the maximum number of intervals")
         end if   

         do i=1,kkc
            endp(i)=i/dble(kkc)

            j=iceil(endp(i)*kk)
            intprob(i,1)=j
            countv1(j)=countv1(j)+1
            countm1(j,countv1(j))=i

            j=iceil(endp(i)*kkc)
            intprob(i,2)=j
            countv2(j)=countv2(j)+1
            countm2(j,countv2(j))=i

         end do
         npoints=kkc-1
       else
         do while(ok.eq.0) 
            tmp1=count1/dble(kk)   
            tmp2=count2/dble(kkc)   

            if(tmp2.lt.tmp1)then
               count=count+1
               if(maxint.lt.(count+1))then
                 call rexit("Error in the maximum number of intervals")
               end if   

               endp(count)=tmp2

               j=iceil(endp(count)*kk)
               intprob(count,1)=j
               countv1(j)=countv1(j)+1
               countm1(j,countv1(j))=count
   
               j=iceil(endp(count)*kkc)
               intprob(count,2)=j
               countv2(j)=countv2(j)+1
               countm2(j,countv2(j))=count

               okc=0
               countc=count2
               do while(okc.eq.0)
                  countc=countc+1
                  tmp2=countc/dble(kkc)
                  if(tmp2.lt.tmp1)then
                     count2=count2+1
                     count=count+1
                     if(maxint.lt.(count+1))then
                       call rexit("Error maximum number of intervals")
                     end if   

                     endp(count)=tmp2 

                     j=iceil(endp(count)*kk)
                     intprob(count,1)=j
                     countv1(j)=countv1(j)+1
                     countm1(j,countv1(j))=count

                     j=iceil(endp(count)*kkc)
                     intprob(count,2)=j
                     countv2(j)=countv2(j)+1
                     countm2(j,countv2(j))=count

                   else
                     okc=1
                     count=count+1
                     if(maxint.lt.(count+1))then
                       call rexit("Error maximum number of intervals")
                     end if   

                     endp(count)=tmp1

                     j=iceil(endp(count)*kk)
                     intprob(count,1)=j
                     countv1(j)=countv1(j)+1
                     countm1(j,countv1(j))=count

                     j=iceil(endp(count)*kkc)
                     intprob(count,2)=j
                     countv2(j)=countv2(j)+1
                     countm2(j,countv2(j))=count

                     count1=count1+1 
                     count2=countc 
                  end if
               end do
             else if(tmp1.lt.tmp2)then

               count=count+1
               if(maxint.lt.(count+1))then
                  call rexit("Error maximum number of intervals")
               end if   

               endp(count)=tmp1

               j=iceil(endp(count)*kk)
               intprob(count,1)=j
               countv1(j)=countv1(j)+1
               countm1(j,countv1(j))=count

               j=iceil(endp(count)*kkc)
               intprob(count,2)=j
               countv2(j)=countv2(j)+1
               countm2(j,countv2(j))=count

               okc=0
               countc=count1
               do while(okc.eq.0)
                  countc=countc+1
                  tmp1=countc/dble(kk)
                  if(tmp1.lt.tmp2)then
                     count1=count1+1
                     count=count+1
                     if(maxint.lt.(count+1))then
                       call rexit("Error maximum number of intervals")
                     end if   

                     endp(count)=tmp1 

                     j=iceil(endp(count)*kk)
                     intprob(count,1)=j
                     countv1(j)=countv1(j)+1
                     countm1(j,countv1(j))=count

                     j=iceil(endp(count)*kkc)
                     intprob(count,2)=j
                     countv2(j)=countv2(j)+1
                     countm2(j,countv2(j))=count

                   else

                     okc=1
                     count=count+1
                     if(maxint.lt.(count+1))then
                       call rexit("Error maximum number of intervals")
                     end if   

                     endp(count)=tmp2

                     j=iceil(endp(count)*kk)
                     intprob(count,1)=j
                     countv1(j)=countv1(j)+1
                     countm1(j,countv1(j))=count

                     j=iceil(endp(count)*kkc)
                     intprob(count,2)=j
                     countv2(j)=countv2(j)+1
                     countm2(j,countv2(j))=count

                     count2=count2+1 
                     count1=countc 
                  end if
               end do   
            end if

            if(count1.ge.kk.and.count2.ge.kkc)ok=1
         end do
         npoints=count-1

      end if
      
      return
      end 



