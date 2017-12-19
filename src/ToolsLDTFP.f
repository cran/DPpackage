c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR LDTFP
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
      subroutine hpddensregdtfp(nsave,npred,ngrid,worksam,worksam2,fs,
     &                          llower,lupper,tband)
c=======================================================================
c     Compute CI for the conditional densities.
c     Alejandro Jara, 2007-2011
c=======================================================================
      implicit none 

c++++ External parameters
      integer nsave,npred,ngrid,tband

c++++ External working
      double precision fs(ngrid)
      double precision worksam(nsave)
      double precision worksam2(nsave,ngrid)

c++++ Output      
      double precision llower(npred,ngrid)
      double precision lupper(npred,ngrid)

c++++ Internal parameters
      double precision alpha
      double precision aupp(2),alow(2)

c++++ Internal working
      integer i,ii,j,l   

c++++ Algorithm

      alpha=0.05d0
  
      open(unit=1,file='dppackage_dftp1.out',status='old',
     &     form='unformatted')

       do ii=1,npred
          do i=1,nsave 
             call rchkusr()
             do j=1,npred
                read(1) (fs(l),l=1,ngrid)
                if(ii.eq.j)then
                   do l=1,ngrid
                      worksam2(i,l)=fs(l) 
                   end do
                end if
             end do
           end do  
           rewind(unit=1)
          
           do i=1,ngrid
              do j=1,nsave
                 worksam(j)=worksam2(j,i) 
              end do
          
              call hpd(nsave,alpha,worksam,alow,aupp)
           
              if(tband.eq.1)then
                 llower(ii,i)=alow(2)
                 lupper(ii,i)=aupp(2)
                else
                 llower(ii,i)=alow(1)
                 lupper(ii,i)=aupp(1)
              end if

           end do
       end do      

       close(unit=1)      

       return
       end


c=======================================================================      
      subroutine hpddensregdtfps(nsave,npred,ngrid,worksam,worksam2,fs,
     &                           llower,lupper,tband)
c=======================================================================
c     Compute CI for the conditional survival functions.
c     Alejandro Jara, 2007-2011
!=======================================================================
      implicit none 

c++++ External parameters
      integer nsave,npred,ngrid,tband

c++++ External working
      double precision fs(ngrid)
      double precision worksam(nsave)
      double precision worksam2(nsave,ngrid)

c++++ Output      
      double precision llower(npred,ngrid)
      double precision lupper(npred,ngrid)

c++++ Internal parameters
      double precision alpha
      double precision aupp(2),alow(2)

c++++ Internal working
      integer i,ii,j,l   

c++++ Algorithm

      alpha=0.05d0
  
      open(unit=2,file='dppackage_dftp2.out',status='old',
     &     form='unformatted')

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            do j=1,npred
               read(2) (fs(l),l=1,ngrid)
               if(ii.eq.j)then
                  do l=1,ngrid
                     worksam2(i,l)=fs(l) 
                  end do
               end if
            end do
          end do  
          rewind(unit=2)
         
          do i=1,ngrid
             do j=1,nsave
                worksam(j)=worksam2(j,i) 
             end do
         
             call hpd(nsave,alpha,worksam,alow,aupp)
          
             if(tband.eq.1)then
                llower(ii,i)=alow(2)
                lupper(ii,i)=aupp(2)
               else
                llower(ii,i)=alow(1)
                lupper(ii,i)=aupp(1)
             end if
          end do
      end do      

      close(unit=2)      
      return
      end

c=======================================================================
      subroutine quantileldtfp(maxm,ntlr,npred,pce,ptf,betace,betatf,
     &                         sigma2,xtfpred,xcepred,
     &                         qquanm,qquanw,k,prob,probc,qqnum)
c=======================================================================
c     Alejandro Jara, 2008-2011
c=======================================================================
      implicit none

c++++ Input
      integer maxm,ntlr,npred,pce,ptf
      double precision betace(pce)
      double precision betatf(ntlr,ptf)
      double precision sigma2
      double precision xtfpred(npred,ptf)
      double precision xcepred(npred,pce)
      double precision qqnum(3)

c++++ Output
      double precision qquanm(npred,3)
      double precision qquanw(npred,3)

c++++ External working space
      integer k(maxm)
      double precision prob(2**maxm)
      double precision probc(2**maxm)
    
c++++ working space - scalar
      integer i,ii,j,j1,j2,k1,k2,ll,m
      integer kphi
      integer nsets 
      integer poss1,poss2,poss3
      double precision invcdfnorm
      double precision loglik
      double precision qworkr
      double precision tmp1,tmp2,tmp3

c++++ Algorithm

      nsets=2**maxm
 
      do i=1,nsets
         prob(i)=0.0
         probc(i)=0.0         
      end do

      do i=1,maxm
         k(i)=0
      end do
 
      do ii=1,npred

         poss1=0
         poss2=0
         poss3=0
         do i=1,nsets
     
            if(i.eq.1)then
               do j1=1,maxm
                  k(j1)=1
               end do
             else
              tmp2=real(i-1)/real(nsets)
              do j1=1,maxm
                 kphi=int(real(2**j1)*tmp2+1)
                 k(j1)=kphi
              end do
           end if 

           loglik=0.d0  

c+++++++++ check if the user has requested an interrupt
           call rchkusr()
     
           j=0
           m=0
           tmp1=0.0
           do j1=1,maxm
              if(j1.eq.1)then
                 j2=k(j1)
                 k2=1
               else 
                 j2=j+k(j1)
                 k2=m+k(j1-1)
              end if   
              j=j+2**j1
              m=m+2**(j1-1)
  
              if(j1.eq.1)then
                 ll=1
                else
                 ll=(k(j1-1)-1)*2+1
              end if

              tmp2=0.0
              do k1=1,ptf
                 tmp2=tmp2+xtfpred(ii,k1)*betatf(k2,k1)
              end do
              tmp3=exp(tmp2)/(1.0+exp(tmp2))

              if(k(j1).eq.ll)then
                 loglik=loglik+log(tmp3) 
               else
                 tmp3=1.0-tmp3
                 loglik=loglik+log(tmp3) 
              end if   
           end do
        
           prob(i)=exp(loglik)

           if(i.eq.1)then
              probc(i)=prob(i)
             else 
              probc(i)=prob(i)+probc(i-1)
           end if   

           if(poss1.eq.0)then
              if(qqnum(1).lt.probc(i))poss1=i
           end if
        
           if(poss2.eq.0)then
              if(qqnum(2).lt.probc(i))poss2=i
           end if

           if(poss3.eq.0)then
              if(qqnum(3).lt.probc(i))poss3=i
           end if
        
        end do

        tmp1=qqnum(1)-probc(poss1)+real(poss1)*prob(poss1)
        tmp2=real(nsets)*prob(poss1)
        tmp3=tmp1/tmp2

c        qworks=dinvnorm(tmp3)
c        call dblepr("tmp3",-1,tmp3,1)
c        call dblepr("qworks",-1,qworks,1)

        tmp1=0.d0
        do i=1,pce
           tmp1=tmp1+xcepred(ii,i)*betace(i)      
        end do
        qworkr=invcdfnorm(tmp3,tmp1,sqrt(sigma2),1,0) 
  
c        qworkr2=tmp1+sqrt(sigma2)*qworks
c        call dblepr("q1",-1,qworkr,1)
c        call dblepr("q1",-1,qworkr2,1)

        qquanm(ii,1)=qquanm(ii,1)+qworkr
        qquanw(ii,1)=qworkr

        tmp1=qqnum(2)-probc(poss2)+real(poss2)*prob(poss2)
        tmp2=real(nsets)*prob(poss2)
        tmp3=tmp1/tmp2
        tmp1=0.d0
        do i=1,pce
           tmp1=tmp1+xcepred(ii,i)*betace(i)      
        end do
        qworkr=invcdfnorm(tmp3,tmp1,sqrt(sigma2),1,0) 

        qquanm(ii,2)=qquanm(ii,2)+qworkr
        qquanw(ii,2)=qworkr

        tmp1=qqnum(3)-probc(poss3)+real(poss3)*prob(poss3)
        tmp2=real(nsets)*prob(poss3)
        tmp3=tmp1/tmp2 
        tmp1=0.d0
        do i=1,pce
           tmp1=tmp1+xcepred(ii,i)*betace(i)      
        end do
        qworkr=invcdfnorm(tmp3,tmp1,sqrt(sigma2),1,0) 

c        call dblepr("tmp3",-1,"tmp3",1)
c        call dblepr("q3",-1,"qworkr",1)

        qquanm(ii,3)=qquanm(ii,3)+qworkr
        qquanw(ii,3)=qworkr

      end do

      return
      end
 
c=======================================================================
      subroutine cdfldtfp(ii,val,maxm,ntlr,nrec,pce,ptf,betace,
     &                    betatf,sigma2,xtf,xce,cdfval,k,prob,probc)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none

c++++ Input
      integer ii  
      integer maxm,ntlr,nrec,pce,ptf
      double precision betace(pce)
      double precision betatf(ntlr,ptf)
      double precision sigma2
      double precision xtf(nrec,ptf)
      double precision xce(nrec,pce)
      double precision val

c++++ Output
      double precision cdfval

c++++ External working space
      integer k(maxm)
      double precision prob(2**maxm)
      double precision probc(2**maxm)

c++++ Internal working space - scalar
      integer i,j,j1,j2,k1,k2,ll,m
      integer kphi
      integer nsets 
      integer possi
      double precision cdfnorm
      double precision loglik
      double precision accsum
      double precision cdfw
      double precision tmp1,tmp2,tmp3
      double precision ee  

c++++ Algorithm

      nsets=2**maxm
      do i=1,nsets
         prob(i)=0.d0
         probc(i)=0.d0
      end do
      do i=1,maxm
         k(i)=0
      end do 
 
      tmp2=0.d0
      do i=1,pce
         tmp2=tmp2+xce(ii,i)*betace(i)
      end do
      tmp1=(val-tmp2)/dsqrt(sigma2)

      if(tmp1.gt.4.0)then
         tmp2=0.999968
       else if(tmp1.lt.-4.0)then
         tmp2=0.000032
       else 
         tmp2=cdfnorm(val,tmp2,sqrt(sigma2),1,0)
      end if
      cdfw=tmp2
 
      kphi=int(real(2**maxm)*tmp2+1)
      possi=kphi

      accsum=0.d0
      do i=1,nsets
     
         if(i.eq.1)then
            do j1=1,maxm
               k(j1)=1
            end do
          else
            tmp2=real(i-1)/real(nsets)
            do j1=1,maxm
               kphi=int(real(2**j1)*tmp2+1)
               k(j1)=kphi
            end do
         end if 

         loglik=0.d0  

c+++++++ check if the user has requested an interrupt
         call rchkusr()
     
         j=0
         m=0
         tmp1=0.d0
         do j1=1,maxm
            if(j1.eq.1)then
               j2=k(j1)
               k2=1
             else 
               j2=j+k(j1)
               k2=m+k(j1-1)
            end if   
            j=j+2**j1
            m=m+2**(j1-1)
  
            if(j1.eq.1)then
               ll=1
              else
               ll=(k(j1-1)-1)*2+1
            end if

            tmp2=0.d0
            do k1=1,ptf
               tmp2=tmp2+xtf(ii,k1)*betatf(k2,k1)
            end do
            tmp3=exp(tmp2)/(1.0+exp(tmp2))

            if(k(j1).eq.ll)then
               loglik=loglik+log(tmp3) 
             else
               tmp3=1.0-tmp3
               loglik=loglik+log(tmp3) 
            end if   
         end do
        
         prob(i)=dexp(loglik)
         accsum=accsum+prob(i)
      end do

      do i=1,nsets
         prob(i)=prob(i)/accsum
         if(i.eq.1)then
            probc(i)=prob(i)
           else 
            probc(i)=prob(i)+probc(i-1)
         end if   
      end do

      if(possi.eq.1)then
         cdfval=prob(possi)*(real(nsets)*cdfw-real(possi-1))
       else
         cdfval=probc(possi-1)+prob(possi)*
     &          (real(nsets)*cdfw-real(possi-1)) 
       end if
  
      if(cdfval.gt.1.0)then
         cdfval=1.d0

c         call intpr("possi",-1,possi,1)
c         call dblepr("betace",-1,betace,pce)
c         call dblepr("val",-1,val,1)
c         call dblepr("ress",-1,ee,1)
c         call dblepr("cdfw",-1,cdfw,1)
c         call dblepr("cdfval",-1,cdfval,1)
c         call rexit("error in cdf: upper")
      end if  

      if(cdfval.lt.0.d0)then
         call intpr("possi",-1,possi,1)
         call dblepr("betace",-1,betace,pce)
         call dblepr("val",-1,val,1)
         call dblepr("ress",-1,ee,1)
         call dblepr("cdfw",-1,cdfw,1)
         call dblepr("cdfval",-1,cdfval,1)
         call rexit("error in cdf: lower")
      end if    
  
      return
      end


c=======================================================================
      subroutine ldensldtfp(yc,ii,maxm,nrec,ntlr,ntprob,pce,
     &                      betace,xce,ptf,betatf,xtf,
     &                      sigma2,loglik,k)
c=======================================================================
      implicit none
c++++ Input
      integer ii
      integer maxm,nrec,ntlr,ntprob,pce,ptf
      double precision yc
      double precision betace(pce)
      double precision betatf(ntlr,ptf)
      double precision sigma2
      double precision xce(nrec,pce)  
      double precision xtf(nrec,ptf)

c++++ Output
      double precision loglik

c++++ External working space
      integer k(maxm)

c++++ Internal working space
      integer i,j,j1,j2,k1,k2,m,ll
      integer kphi
      double precision cdfnorm
      double precision dnrm
      double precision tmp1,tmp2,tmp3

c++++ Algorithm

      do i=1,maxm
         k(i)=0
      end do

c++++ check if the user has requested an interrupt
      call rchkusr()
  
      tmp2=0.d0
      do i=1,pce
         tmp2=tmp2+xce(ii,i)*betace(i)
      end do

      loglik=dnrm(yc,tmp2,sqrt(sigma2),1)

      tmp1=(yc-tmp2)/sqrt(sigma2)

      if(tmp1.gt.4.0)then
         tmp2=0.999968
        else if(tmp1.lt.-4.0)then
         tmp2=0.000032
        else 
         tmp2=cdfnorm(yc,tmp2,sqrt(sigma2),1,0)
       end if

       do j1=1,maxm
          kphi=int(real(2**j1)*tmp2+1)
          k(j1)=kphi
       end do

       j=0
       m=0
       do j1=1,maxm
          if(j1.eq.1)then
             j2=k(j1)
             k2=1
           else 
             j2=j+k(j1)
             k2=m+k(j1-1)
          end if   
          j=j+2**j1
          m=m+2**(j1-1)
          if(j1.eq.1)then
             ll=1
            else
             ll=(k(j1-1)-1)*2+1
          end if

          tmp2=0.d0
          do k1=1,ptf
             tmp2=tmp2+betatf(k2,k1)*xtf(ii,k1)
          end do
          tmp3=exp(tmp2)/(1.0+exp(tmp2))

          if(k(j1).eq.ll)then
             loglik=loglik+log(tmp3) 
           else
             tmp3=1.0-tmp3
             loglik=loglik+log(tmp3) 
          end if   
       end do
       loglik=loglik+real(maxm)*log(2.0) 
   
       return
       end


c=======================================================================
      subroutine densldtfpaft(maxm,ngrid,ntlr,ntprob,npred,pce,ptf,
     &                        betace,betatf,grid,xcepred,xtfpred,
     &                        sigma2,densw,densm,k)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none

c++++ Input
      integer maxm,ngrid,ntlr,ntprob,npred,pce,ptf
      double precision betace(pce)
      double precision betatf(ntlr,ptf)
      double precision grid(ngrid)
      double precision sigma2
      double precision xtfpred(npred,ptf)
      double precision xcepred(npred,pce)

c++++ Output
      double precision densw(npred,ngrid)
      double precision densm(npred,ngrid)

c++++ External working space
      integer k(maxm)

c++++ Internal working space
      integer ii,jj,i,j,j1,k1,k2,m,ll
      integer kphi
      double precision dnrm
      double precision cdfnorm
      double precision tmp
      double precision loglik
      double precision tmp1,tmp2,tmp3

c++++ Algorithm

      do i=1,maxm
         k(i)=0
      end do   

c++++ check if the user has requested an interrupt
      call rchkusr()
  
      do ii=1,npred

        do jj=1,ngrid
     
           tmp1=0.d0
           do j=1,pce
              tmp1=tmp1+xcepred(ii,j)*betace(j)      
           end do
           tmp=(log(grid(jj))-tmp1)/sqrt(sigma2)

           loglik=dnrm(log(grid(jj)),tmp1,sqrt(sigma2),1)

           if(tmp.gt.4.0)then
              tmp2=0.999968
             else if(tmp.lt.-4.0)then
              tmp2=0.000032
             else 
              tmp2=cdfnorm(log(grid(jj)),tmp1,sqrt(sigma2),1,0)
           end if

           do j1=1,maxm
              kphi=int(real(2**j1)*tmp2+1)
              k(j1)=kphi
           end do

c+++++++++ check if the user has requested an interrupt
           call rchkusr()
     
           j=0
           m=0
           tmp1=0.d0
           do j1=1,maxm
              if(j1.eq.1)then
                 k2=1
               else 
                 k2=m+k(j1-1)
              end if   
              m=m+2**(j1-1)
  
              if(j1.eq.1)then
                 ll=1
                else
                 ll=(k(j1-1)-1)*2+1
              end if

              tmp2=0.d0
              do k1=1,ptf
                 tmp2=tmp2+xtfpred(ii,k1)*betatf(k2,k1)
              end do
              tmp3=exp(tmp2)/(1.0+exp(tmp2))

              if(k(j1).eq.ll)then
                 loglik=loglik+log(tmp3) 
               else
                 tmp3=1.0-tmp3
                 loglik=loglik+log(tmp3) 
              end if   
           
           end do

           loglik=loglik + real(maxm)*log(2.0) 
           loglik=loglik - log(grid(jj))

           densw(ii,jj)=exp(loglik)
           densm(ii,jj)=densm(ii,jj)+exp(loglik)
       end do
      end do
      return
      end

c=======================================================================
      subroutine loglikldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,betace,
     &                       betatf,y,xce,xtf,sigma2,
     &                       nobsbc,obsbc,loglik,k)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none

c++++ Input
      integer maxm,ntlr,ntprob,nrec,pce,ptf
      double precision betace(pce)
      double precision betatf(ntlr,ptf)
      double precision y(nrec)
      double precision sigma2
      double precision xtf(nrec,ptf)
      double precision xce(nrec,pce)

c++++ Output
      integer nobsbc(ntprob)
      integer obsbc(ntprob,nrec)
      double precision loglik

c++++ External working space
      integer k(maxm)

c++++ Internal working space
      integer i,j,j1,j2,k1,k2,m,ll
      integer kphi
      double precision cdfnorm
      double precision dnrm
      double precision tmp1,tmp2,tmp3

c++++ Algorithm
 
      do i=1,maxm
         k(i)=0
      end do

      do i=1,ntprob  
         nobsbc(i)=0
      end do
      loglik=0.d0
   
      do i=1,nrec
 
c+++++++ check if the user has requested an interrupt
         call rchkusr()

         tmp1=0.d0
         do j=1,pce
            tmp1=tmp1+xce(i,j)*betace(j)      
         end do
         tmp3=(y(i)-tmp1)/dsqrt(sigma2)

         loglik=loglik+dnrm(y(i),tmp1,sqrt(sigma2),1)

         if(tmp3.gt.4.0)then
            tmp2=0.999968
           else if(tmp3.lt.-4.0)then
            tmp2=0.000032
           else 
            tmp2=cdfnorm(y(i),tmp1,sqrt(sigma2),1,0)
         end if

         do j1=1,maxm
            kphi=int(real(2**j1)*tmp2+1)
            k(j1)=kphi
         end do

c         call intpr("k",-1,k,maxm)
     
         j=0
         m=0
         tmp1=0.d0
         do j1=1,maxm
            if(j1.eq.1)then
               j2=k(j1)
               k2=1
             else 
               j2=j+k(j1)
               k2=m+k(j1-1)
            end if   
            j=j+2**j1
            m=m+2**(j1-1)

            nobsbc(j2)=nobsbc(j2)+1
            obsbc(j2,nobsbc(j2))=i

            if(j1.eq.1)then
               ll=1
              else
               ll=(k(j1-1)-1)*2+1
            end if

            tmp2=0.d0
            do k1=1,ptf
               tmp2=tmp2+xtf(i,k1)*betatf(k2,k1)
            end do
            tmp3=exp(tmp2)/(1.0+exp(tmp2))

            if(k(j1).eq.ll)then
               loglik=loglik+log(tmp3) 
             else
               tmp3=1.0-tmp3
               loglik=loglik+log(tmp3) 
            end if   
         end do
         loglik=loglik+real(maxm)*log(2.0) 
    
      end do

c      call intpr("nobsbc",-1,nobsbc,ntprob)
      return
      end


c=======================================================================
      subroutine startlrcoefldtfp(nri,kk,ii,jj,n1,n2,maxm,ntlr,
     &                            ntprob,nrec,ptf,obsbc,betatf,xtf,c0,
     &                            iflagp,beta,xtx,xty)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none

c++++ Input
      integer nri
      integer kk,ii,jj,n1,n2
      integer maxm,ntlr,ntprob,nrec,ptf
      integer obsbc(ntprob,nrec)
      double precision betatf(ntlr,ptf)
      double precision c0(ptf,ptf)
      double precision xtf(nrec,ptf)

c++++ External working space
      integer iflagp(ptf)
      double precision beta(ptf)
      double precision xtx(ptf,ptf)
      double precision xty(ptf)

c++++ Working space
      integer i,j,k,ll,nhr
      double precision eta,gprime,mu,tmp1,ytilde

c++++ Algorithm

      do i=1,ptf
         beta(i)=betatf(kk,i)
      end do
  
c++++ MLE with nri stept of N-R

      do nhr=1,nri

         do i=1,ptf
            do j=1,ptf 
               xtx(i,j)=c0(i,j)
            end do
            xty(i)=0.d0
         end do
         call inverse(xtx,ptf,iflagp)

         if(n1.gt.0)then
            do i=1,n1
               ll=obsbc(ii,i)
     
c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               eta=0.d0
               do j=1,ptf
                  eta=eta+xtf(ll,j)*beta(j)
               end do
               mu=dexp(eta)/(1.d0+dexp(eta))

               tmp1=dexp(eta)/((1.d0+dexp(eta))*
     &              (1.d0+dexp(eta)))
               gprime=1.d0/tmp1
               ytilde=eta+(1.d0-mu)*gprime 

               do j=1,ptf
                  do k=1,ptf
                     xtx(j,k)=xtx(j,k)+xtf(ll,j)*xtf(ll,k)*tmp1
                  end do
                  xty(j)=xty(j)+xtf(ll,j)*ytilde*tmp1
               end do
            end do
         end if

         if(n2.gt.0)then
            do i=1,n1
               ll=obsbc(jj,i)
     
c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               eta=0.d0
               do j=1,ptf
                  eta=eta+xtf(ll,j)*beta(j)
               end do
               mu=dexp(eta)/(1.d0+dexp(eta))

               tmp1=dexp(eta)/((1.d0+dexp(eta))*
     &              (1.d0+dexp(eta)))
               gprime=1.d0/tmp1
               ytilde=eta+(0.d0-mu)*gprime 

               do j=1,ptf
                  do k=1,ptf
                     xtx(j,k)=xtx(j,k)+
     &                        xtf(ll,j)*xtf(ll,k)*tmp1
                  end do
                  xty(j)=xty(j)+xtf(ll,j)*ytilde*tmp1
               end do
            end do
         end if

         call inverse(xtx,ptf,iflagp)

         do i=1,ptf
            tmp1=0.d0
            do j=1,ptf
               tmp1=tmp1+xtx(i,j)*xty(j) 
            end do
            beta(i)=tmp1
         end do
      end do

      do i=1,ptf
         betatf(kk,i)=beta(i)
      end do

      return
      end


c=======================================================================
      subroutine updatelrcoefldtfp0(kk,maxm,ntlr,ntprob,
     &                              ptf,betatf,c0,
     &                              iflagp,beta,betam,xtx,
     &                              workmhp,workvp)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none

c++++ Input
      integer kk
      integer maxm,ntlr,ntprob,ptf
      double precision betatf(ntlr,ptf)
      double precision c0(ptf,ptf)

c++++ Exetrnal working space
      integer iflagp(ptf)
      double precision beta(ptf)
      double precision betam(ptf)
      double precision xtx(ptf,ptf)
      double precision workmhp(ptf*(ptf+1)/2)
      double precision workvp(ptf)

c++++ Working space
      integer i,j

c++++ Algorithm
 
      do i=1,ptf
         do j=1,ptf
            xtx(i,j)=c0(i,j)
         end do
         betam(i)=0.d0
      end do

      call rmvnorm(ptf,betam,xtx,workmhp,workvp,beta)

      do i=1,ptf
         betatf(kk,i)=beta(i)
      end do
      return
      end
 
c=======================================================================
      subroutine updatelrcoefldtfp(kk,ii,jj,n1,n2,maxm,ntlr,ntprob,
     &                             nrec,ptf,obsbc,betatf,xtf,c0,
     &                             accept,iflagp,beta,betam,betac,
     &                             xtx,xty,workmhp)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none

c++++ Input
      integer kk,ii,jj,n1,n2
      integer maxm,ntlr,ntprob,nrec,ptf
      integer obsbc(ntprob,nrec)
      double precision betatf(ntlr,ptf)
      double precision xtf(nrec,ptf)
      double precision c0(ptf,ptf)

c++++ Output
      integer accept

c++++ External working space
      integer iflagp(ptf)
      double precision beta(ptf)
      double precision betam(ptf)
      double precision betac(ptf)
      double precision xtx(ptf,ptf)
      double precision xty(ptf)
      double precision workmhp(ptf*(ptf+1)/2)

c++++ Internal working space
      integer i,j,k,ll
      double precision eta,mu,tmp1,ytilde
      double precision lratio,luni
      double precision logliko,loglikn
      double precision logprioro,logpriorn
      double precision logcgko,logcgkn
      real runif

c++++ Algorithm

      logliko=0.d0
      loglikn=0.d0
      logprioro=0.d0
      logpriorn=0.d0
      logcgko=0.d0
      logcgkn=0.d0
  
      accept=0
  
      do i=1,ptf
         beta(i)=betatf(kk,i)
      end do
  
c++++ generating the candidate

      do i=1,ptf
         do j=1,ptf 
            xtx(i,j)=c0(i,j)
         end do
         xty(i)=0.d0
      end do

      call inverse(xtx,ptf,iflagp)

      if(n1.gt.0)then
         do i=1,n1
            ll=obsbc(ii,i)
     
c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            eta=0.d0
            do j=1,ptf
               eta=eta+xtf(ll,j)*beta(j)
            end do

            logliko=logliko+eta-log(1.0+exp(eta))
            mu=exp(eta)/(1.0+exp(eta))

            tmp1=mu*(1.0-mu)
            ytilde=eta+(1.0+exp(eta))/exp(eta)
 
            do j=1,ptf
               do k=1,ptf
                  xtx(j,k)=xtx(j,k)+
     &                     xtf(ll,j)*xtf(ll,k)*tmp1
               end do
               xty(j)=xty(j)+xtf(ll,j)*ytilde*tmp1
            end do
         end do
      end if

      if(n2.gt.0)then
         do i=1,n2
            ll=obsbc(jj,i)
     
c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            eta=0.d0
            do j=1,ptf
               eta=eta+xtf(ll,j)*beta(j)
            end do
            mu=exp(eta)/(1.0+exp(eta))

            logliko=logliko-log(1.0+exp(eta))

            tmp1=mu*(1.0-mu)
            ytilde=eta-(1.0+exp(eta))

            do j=1,ptf
               do k=1,ptf
                  xtx(j,k)=xtx(j,k)+
     &                     xtf(ll,j)*xtf(ll,k)*tmp1
               end do
               xty(j)=xty(j)+xtf(ll,j)*ytilde*tmp1
            end do
         end do
       end if

       call inverse(xtx,ptf,iflagp)

       do i=1,ptf
         tmp1=0.d0
         do j=1,ptf
             tmp1=tmp1+xtx(i,j)*xty(j) 
          end do
          betam(i)=tmp1
       end do

       call rmvnorm(ptf,betam,xtx,workmhp,xty,betac)

c+++++ evaluating the candidate generating kernel

       call dmvnd(ptf,betac,betam,xtx,logcgko,iflagp)

c+++++ evaluating the likelihood for the candidate

       do i=1,ptf
          do j=1,ptf 
             xtx(i,j)=c0(i,j)
          end do
          xty(i)=0.d0
       end do

       call inverse(xtx,ptf,iflagp)

       if(n1.gt.0)then
          do i=1,n1
             ll=obsbc(ii,i)
     
c+++++++++++ check if the user has requested an interrupt
             call rchkusr()

             eta=0.d0
             do j=1,ptf
                eta=eta+xtf(ll,j)*betac(j)
             end do
             mu=exp(eta)/(1.0+exp(eta))

             loglikn=loglikn+eta-log(1.0+exp(eta))

             tmp1=mu*(1.0-mu)
             ytilde=eta+(1.0+exp(eta))/exp(eta)

             do j=1,ptf
                do k=1,ptf
                   xtx(j,k)=xtx(j,k)+
     &                      xtf(ll,j)*xtf(ll,k)*tmp1
                end do
                xty(j)=xty(j)+xtf(ll,j)*ytilde*tmp1
             end do
          end do
       end if

       if(n2.gt.0)then
          do i=1,n2
             ll=obsbc(jj,i)
     
c+++++++++++ check if the user has requested an interrupt
             call rchkusr()

             eta=0.d0
             do j=1,ptf
                eta=eta+xtf(ll,j)*betac(j)
             end do
             mu=exp(eta)/(1.0+exp(eta))

             loglikn=loglikn-log(1.0+exp(eta))

             tmp1=mu*(1.0-mu)
             ytilde=eta-(1.0+exp(eta))
 
             do j=1,ptf
                do k=1,ptf
                   xtx(j,k)=xtx(j,k)+
     &                      xtf(ll,j)*xtf(ll,k)*tmp1
                end do
                xty(j)=xty(j)+xtf(ll,j)*ytilde*tmp1
             end do
          end do
       end if

       call inverse(xtx,ptf,iflagp)

       do i=1,ptf
          tmp1=0.d0
          do j=1,ptf
             tmp1=tmp1+xtx(i,j)*xty(j) 
          end do
          betam(i)=tmp1
       end do

       call dmvnd(ptf,beta,betam,xtx,logcgkn,iflagp)

c+++++ evaluating the prior

       do i=1,ptf
          do j=1,ptf 
             xtx(i,j)=c0(i,j)
          end do
          xty(i)=0.d0
       end do

       call inverse(xtx,ptf,iflagp)
   
       do i=1,ptf
          do j=1,ptf
             logpriorn=logpriorn+
     &                 betac(i)*xtx(i,j)*betac(j)
             logprioro=logprioro+
     &                 beta(i)*xtx(i,j)*beta(j)
          end do
       end do
       logpriorn=-0.5d0*logpriorn
       logprioro=-0.5d0*logprioro

c+++++ mh step
       lratio=loglikn-logliko+logcgkn-logcgko+
     &        logpriorn-logprioro
       luni=log(real(runif()))

       if(luni.lt.lratio)then
          accept=1
          do i=1,ptf
             betatf(kk,i)=betac(i)
          end do
       end if
      
       return
      end


c=======================================================================
      subroutine updatelrcoefldtfpss(kk,ii,jj,n1,n2,maxm,ntlr,ntprob,
     &                               nrec,ptf,obsbc,betatf,xtf,c0,
     &                               accept,iflagp,beta,xtx)
c=======================================================================
c     Alejandro Jara, 2011
c     Uses Slice Sampling instead of WLS MH.
c=======================================================================
      implicit none

c++++ Input
      integer kk,ii,jj,n1,n2
      integer maxm,ntlr,ntprob,nrec,ptf
      integer obsbc(ntprob,nrec)
      double precision betatf(ntlr,ptf)
      double precision xtf(nrec,ptf)
      double precision c0(ptf,ptf)

c++++ Output
      integer accept

c++++ External working space
      integer iflagp(ptf)
      double precision beta(ptf)
      double precision xtx(ptf,ptf)

c++++ Internal working space
      integer i,j
      double precision tmp1
      double precision uni
      real runif

c++++ Working space slice sampling
      integer evali
      double precision re,uwork
      double precision logy,xx0,xx1,llim,rlim
      double precision grlim,gllim,gxx1

c++++ Algorithm

      accept=1
  
      do i=1,ptf
         beta(i)=betatf(kk,i)
      end do
  
      do i=1,ptf
         do j=1,ptf 
            xtx(i,j)=c0(i,j)
         end do
      end do
      call inverse(xtx,ptf,iflagp)

      do i=1,ptf

c+++++++ check if the user has requested an interrupt
         call rchkusr()
       
         evali=1
         xx0=beta(i)

         call compullldtfp(ii,jj,ptf,n1,n2,nrec,ntprob,obsbc,
     &                     beta,xtx,xtf,tmp1)

c         call dblepr("tmp1",-1,tmp1,1)

         uni=runif()
         re=-log(uni)
         logy=tmp1-re

         uwork=dble(runif())*0.25d0  
         llim=xx0-uwork
         rlim=xx0+(0.25d0-uwork)
           
         evali=evali+1
         beta(i)=llim
         call compullldtfp(ii,jj,ptf,n1,n2,nrec,ntprob,obsbc,
     &                     beta,xtx,xtf,gllim)

c         call dblepr("llim",-1,llim,1)

         evali=evali+1
         beta(i)=rlim
         call compullldtfp(ii,jj,ptf,n1,n2,nrec,ntprob,obsbc,
     &                     beta,xtx,xtf,grlim)


c         call dblepr("rlim",-1,rlim,1)

         do while(gllim.gt.logy)
            llim=llim-0.25d0
            evali=evali+1
            beta(i)=llim
            call compullldtfp(ii,jj,ptf,n1,n2,nrec,
     &                        ntprob,obsbc,
     &                        beta,xtx,xtf,gllim)

c            call dblepr("llim",-1,llim,1)

         end do 

         do while(grlim.gt.logy)
            rlim=rlim+0.25d0
            beta(i)=rlim
            call compullldtfp(ii,jj,ptf,n1,n2,nrec,
     &                        ntprob,obsbc,
     &                        beta,xtx,xtf,grlim)

c            call dblepr("rlim",-1,rlim,1)

         end do 

         xx1=llim+(rlim-llim)*dble(runif())
         beta(i)=xx1
         call compullldtfp(ii,jj,ptf,n1,n2,nrec,ntprob,obsbc,
     &                     beta,xtx,xtf,gxx1)

c         call dblepr("xx1",-1,xx1,1)

         do while(gxx1.lt.logy)
            if(xx1.gt.xx0)rlim=xx1
            if(xx1.lt.xx0)llim=xx1

            xx1=llim+(rlim-llim)*dble(runif())
            beta(i)=xx1
            call compullldtfp(ii,jj,ptf,n1,n2,nrec,
     &                        ntprob,obsbc,
     &                        beta,xtx,xtf,gxx1)

c            call dblepr("xx1",-1,xx1,1)
c            call dblepr("gxx1",-1,gxx1,1)
c            call dblepr("logy",-1,logy,1)

         end do
         betatf(kk,i)=beta(i)

c         call intpr("evali",-1,evali,1)

      end do

      return
      end


c=======================================================================
      subroutine compullldtfp(ii,jj,ptf,n1,n2,nrec,ntprob,obsbc,
     &                        beta,xtx,xtf,logp)
c=======================================================================
c     Alejandro Jara, 2011
c=======================================================================
      implicit none

c++++ Input
      integer ii,jj,ptf,n1,n2,nrec,ntprob 
      integer obsbc(ntprob,nrec)
      double precision beta(ptf)
      double precision xtx(ptf,ptf)
      double precision xtf(nrec,ptf)

c++++ Output
      double precision logp
 
c++++ Working space
      integer i,j,ll
      double precision eta
      double precision loglikn
      double precision logpriorn

c++++ Algorithm

      loglikn=0.d0
      logpriorn=0.d0
 
      do i=1,ptf
         do j=1,ptf
            logpriorn=logpriorn+beta(i)*xtx(i,j)*beta(j)
         end do
      end do
      logpriorn=-0.5d0*logpriorn

      if(n1.gt.0)then
         do i=1,n1
            ll=obsbc(ii,i)
c++++++++++ check if the user has requested an interrupt
            call rchkusr()
            eta=0.d0
            do j=1,ptf
               eta=eta+xtf(ll,j)*beta(j)
            end do
            loglikn=loglikn+eta-log(1.0+exp(eta))
         end do
      end if
      if(n2.gt.0)then
         do i=1,n2
            ll=obsbc(jj,i)
c++++++++++ check if the user has requested an interrupt
            call rchkusr()
            eta=0.d0
            do j=1,ptf
               eta=eta+xtf(ll,j)*beta(j)
            end do
            loglikn=loglikn-log(1.0+exp(eta))
         end do
      end if

      logp=loglikn+logpriorn
          
      return
      end


c=======================================================================
      subroutine logposldtfp(maxm,ntlr,ntprob,nrec,pce,ptf,betace,
     &                       betatf,y,xce,xtf,sigma2,betacepm,precce,
     &                       nobsbc,obsbc,logpost,k)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none

c++++ Input
      integer maxm,ntlr,ntprob,nrec,pce,ptf
      double precision betace(pce)
      double precision betatf(ntlr,ptf)
      double precision y(nrec)
      double precision sigma2
      double precision xtf(nrec,ptf)
      double precision xce(nrec,pce)

      double precision betacepm(pce)
      double precision precce(pce,pce)

c++++ Output
      integer nobsbc(ntprob)
      integer obsbc(ntprob,nrec)
      double precision logpost

c++++ External working space
      integer k(maxm)

c++++ Internal working space
      integer i,j,j1,j2,k1,k2,m,ll
      integer kphi
      double precision cdfnorm
      double precision dnrm
      double precision tmp1,tmp2,tmp3
      double precision loglikn
      double precision logpriorn 

c++++ Algorithm
 
      do i=1,maxm
         k(i)=0
      end do

      do i=1,ntprob  
         nobsbc(i)=0
      end do
      loglikn=0.d0
   
      do i=1,nrec
 
c+++++++ check if the user has requested an interrupt
         call rchkusr()

         tmp1=0.d0
         do j=1,pce
            tmp1=tmp1+xce(i,j)*betace(j)      
         end do
         tmp3=(y(i)-tmp1)/dsqrt(sigma2)

         loglikn=loglikn+dnrm(y(i),tmp1,sqrt(sigma2),1)

         if(tmp3.gt.4.0)then
            tmp2=0.999968
           else if(tmp3.lt.-4.0)then
            tmp2=0.000032
           else 
            tmp2=cdfnorm(y(i),tmp1,sqrt(sigma2),1,0)
         end if

         do j1=1,maxm
            kphi=int(real(2**j1)*tmp2+1)
            k(j1)=kphi
         end do

c         call intpr("k",-1,k,maxm)
     
         j=0
         m=0
         tmp1=0.d0
         do j1=1,maxm
            if(j1.eq.1)then
               j2=k(j1)
               k2=1
             else 
               j2=j+k(j1)
               k2=m+k(j1-1)
            end if   
            j=j+2**j1
            m=m+2**(j1-1)

            nobsbc(j2)=nobsbc(j2)+1
            obsbc(j2,nobsbc(j2))=i

            if(j1.eq.1)then
               ll=1
              else
               ll=(k(j1-1)-1)*2+1
            end if

            tmp2=0.d0
            do k1=1,ptf
               tmp2=tmp2+xtf(i,k1)*betatf(k2,k1)
            end do
            tmp3=exp(tmp2)/(1.0+exp(tmp2))

            if(k(j1).eq.ll)then
               loglikn=loglikn+log(tmp3) 
             else
               tmp3=1.0-tmp3
               loglikn=loglikn+log(tmp3) 
            end if   
         end do
         loglikn=loglikn+real(maxm)*log(2.0) 
    
      end do

      logpriorn=0.d0
      do i=1,pce
         do j=1,pce
            logpriorn=logpriorn+(betace(i)-betacepm(i))*
     &                           precce(i,j)*
     &                          (betace(j)-betacepm(j))

         end do
      end do
      logpriorn=-0.5d0*logpriorn

      logpost=loglikn+logpriorn

      return
      end



c=======================================================================
      subroutine densldtfp(maxm,ngrid,ntlr,ntprob,npred,pce,ptf,
     &                     betace,betatf,grid,xcepred,xtfpred,
     &                     sigma2,densw,densm,k)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none

c++++ Input
      integer maxm,ngrid,ntlr,ntprob,npred,pce,ptf
      double precision betace(pce)
      double precision betatf(ntlr,ptf)
      double precision grid(ngrid)
      double precision sigma2
      double precision xtfpred(npred,ptf)
      double precision xcepred(npred,pce)

c++++ Output
      double precision densw(npred,ngrid)
      double precision densm(npred,ngrid)

c++++ External working space
      integer k(maxm)

c++++ Internal working space
      integer ii,jj,i,j,j1,k1,k2,m,ll
      integer kphi
      double precision dnrm
      double precision cdfnorm
      double precision tmp
      double precision loglik
      double precision tmp1,tmp2,tmp3

c++++ Algorithm

      do i=1,maxm
         k(i)=0
      end do   

c++++ check if the user has requested an interrupt
      call rchkusr()
  
      do ii=1,npred

        do jj=1,ngrid
     
           tmp1=0.d0
           do j=1,pce
              tmp1=tmp1+xcepred(ii,j)*betace(j)      
           end do
           tmp=(grid(jj)-tmp1)/sqrt(sigma2)

           loglik=dnrm(grid(jj),tmp1,sqrt(sigma2),1)

           if(tmp.gt.4.0)then
              tmp2=0.999968
             else if(tmp.lt.-4.0)then
              tmp2=0.000032
             else 
              tmp2=cdfnorm(grid(jj),tmp1,sqrt(sigma2),1,0)
           end if

           do j1=1,maxm
              kphi=int(real(2**j1)*tmp2+1)
              k(j1)=kphi
           end do

c+++++++++ check if the user has requested an interrupt
           call rchkusr()
     
           j=0
           m=0
           tmp1=0.d0
           do j1=1,maxm
              if(j1.eq.1)then
                 k2=1
               else 
                 k2=m+k(j1-1)
              end if   
              m=m+2**(j1-1)
  
              if(j1.eq.1)then
                 ll=1
                else
                 ll=(k(j1-1)-1)*2+1
              end if

              tmp2=0.d0
              do k1=1,ptf
                 tmp2=tmp2+xtfpred(ii,k1)*betatf(k2,k1)
              end do
              tmp3=exp(tmp2)/(1.0+exp(tmp2))

              if(k(j1).eq.ll)then
                 loglik=loglik+log(tmp3) 
               else
                 tmp3=1.0-tmp3
                 loglik=loglik+log(tmp3) 
              end if   
           
           end do

           loglik=loglik + real(maxm)*log(2.0) 
           densw(ii,jj)=exp(loglik)
           densm(ii,jj)=densm(ii,jj)+exp(loglik)
       end do
      end do
      return
      end


c=======================================================================      
      subroutine hpdldtfpq2(nsave,npred,worksam,worksam2,llower,
     &                      lupper,tband)
c=======================================================================
c     Compute CI for the median.
c     Alejandro Jara, 2011
c=======================================================================
      implicit none 
c++++ External parameters
      integer nsave,npred,tband

c++++ External working
      double precision worksam(nsave)
      double precision worksam2(nsave,npred)

c++++ Output      
      double precision llower(npred,3)
      double precision lupper(npred,3)

c++++ Internal parameters
      double precision alpha
      double precision aupp(2),alow(2)

c++++ Internal working
      integer i,ii,l   

c++++ Algorithm

      alpha=0.05

      open(unit=4,file='dppackage_dftp4.out',status='old',
     &     form='unformatted')

      do i=1,nsave
         call rchkusr()
         read(4) (worksam2(i,l),l=1,npred)
      end do

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            worksam(i)=worksam2(i,ii)
         end do  
         
         call hpd(nsave,alpha,worksam,alow,aupp)
          
         if(tband.eq.1)then  
            llower(ii,2)=alow(2)
            lupper(ii,2)=aupp(2)
           else
            llower(ii,2)=alow(1)
            lupper(ii,2)=aupp(1)
         end if
      end do      

      close(unit=4)      

      return
      end

c=======================================================================      
      subroutine hpdldtfpq1(nsave,npred,worksam,worksam2,llower,
     &                      lupper,tband)
c=======================================================================
c     Compute CI for Q1.
c     Alejandro Jara, 2011
c=======================================================================
      implicit none 
c++++ External parameters
      integer nsave,npred,tband

c++++ External working
      double precision worksam(nsave)
      double precision worksam2(nsave,npred)

c++++ Output      
      double precision llower(npred,3)
      double precision lupper(npred,3)

c++++ Internal parameters
      double precision alpha
      double precision aupp(2),alow(2)

c++++ Internal working
      integer i,ii,l   

c++++ Algorithm

      alpha=0.05
  
      open(unit=3,file='dppackage_dftp3.out',status='old',
     &     form='unformatted')

      do i=1,nsave
         call rchkusr()
         read(3) (worksam2(i,l),l=1,npred)
      end do

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            worksam(i)=worksam2(i,ii)
         end do  
         
         call hpd(nsave,alpha,worksam,alow,aupp)

         if(tband.eq.1)then  
            llower(ii,1)=alow(2)
            lupper(ii,1)=aupp(2)
           else
            llower(ii,1)=alow(1)
            lupper(ii,1)=aupp(1)
         end if

      end do      

      close(unit=3)      

      return
      end

c=======================================================================      
      subroutine hpdldtfpq3(nsave,npred,worksam,worksam2,llower,
     &                      lupper,tband)
c=======================================================================
c     Compute CI for Q3.
c     Alejandro Jara, 2011
c=======================================================================
      implicit none 
c++++ External parameters
      integer nsave,npred,tband

c++++ External working
      double precision worksam(nsave)
      double precision worksam2(nsave,npred)

c++++ Output      
      double precision llower(npred,3)
      double precision lupper(npred,3)

c++++ Internal parameters
      double precision alpha
      double precision aupp(2),alow(2)

c++++ Internal working
      integer i,ii,l   

c++++ Algorithm

      alpha=0.05
  
      open(unit=5,file='dppackage_dftp5.out',status='old',
     &     form='unformatted')

      do i=1,nsave
         call rchkusr()
         read(5) (worksam2(i,l),l=1,npred)
      end do

      do ii=1,npred
         do i=1,nsave 
            call rchkusr()
            worksam(i)=worksam2(i,ii)
         end do  
        
         call hpd(nsave,alpha,worksam,alow,aupp)

         if(tband.eq.1)then  
            llower(ii,3)=alow(2)
            lupper(ii,3)=aupp(2)
           else
            llower(ii,3)=alow(1)
            lupper(ii,3)=aupp(1)
         end if
      end do      

      close(unit=5)      

      return
      end

c+++++++++++++++++++++++++++++++++++++++
c a routine written by john herrero
c+++++++++++++++++++++++++++++++++++++++
      double precision function dinvnorm(p)
      double precision p,p_low,p_high
      double precision a1,a2,a3,a4,a5,a6
      double precision b1,b2,b3,b4,b5
      double precision c1,c2,c3,c4,c5,c6
      double precision d1,d2,d3,d4
      double precision z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/
     &((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/
     &(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/
     &((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   dinvnorm=z
      return
      end


c=======================================================================
      subroutine loglikldtfpre(maxm,ntlr,ntprob,nsubject,ptf,
     &                         betatf,b,xtf,sigma2b,
     &                         nobsbc,obsbc,loglik,k)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none

c++++ Input
      integer maxm,ntlr,ntprob,nsubject,ptf
      double precision betatf(ntlr,ptf)
      double precision b(nsubject)
      double precision sigma2b
      double precision xtf(nsubject,ptf)

c++++ Output
      integer nobsbc(ntprob)
      integer obsbc(ntprob,nsubject)
      double precision loglik

c++++ External working space
      integer k(maxm)

c++++ Internal working space
      integer i,j,j1,j2,k1,k2,m,ll
      integer kphi
      double precision cdfnorm
      double precision dnrm
      double precision tmp1,tmp2,tmp3

c++++ Algorithm
 
      do i=1,maxm
         k(i)=0
      end do

      do i=1,ntprob  
         nobsbc(i)=0
      end do
      loglik=0.d0
   
      do i=1,nsubject
 
c+++++++ check if the user has requested an interrupt
         call rchkusr()

         tmp3=b(i)/dsqrt(sigma2b)

         loglik=loglik+dnrm(b(i),0.d0,sqrt(sigma2b),1)

         if(tmp3.gt.4.0)then
            tmp2=0.999968
           else if(tmp3.lt.-4.0)then
            tmp2=0.000032
           else 
            tmp2=cdfnorm(b(i),0.d0,sqrt(sigma2b),1,0)
         end if

         do j1=1,maxm
            kphi=int(real(2**j1)*tmp2+1)
            k(j1)=kphi
         end do

         j=0
         m=0
         tmp1=0.d0
         do j1=1,maxm
            if(j1.eq.1)then
               j2=k(j1)
               k2=1
             else 
               j2=j+k(j1)
               k2=m+k(j1-1)
            end if   
            j=j+2**j1
            m=m+2**(j1-1)

            nobsbc(j2)=nobsbc(j2)+1
            obsbc(j2,nobsbc(j2))=i

            if(j1.eq.1)then
               ll=1
              else
               ll=(k(j1-1)-1)*2+1
            end if

            tmp2=0.d0
            do k1=1,ptf
               tmp2=tmp2+xtf(i,k1)*betatf(k2,k1)
            end do
            tmp3=exp(tmp2)/(1.0+exp(tmp2))

            if(k(j1).eq.ll)then
               loglik=loglik+log(tmp3) 
             else
               tmp3=1.0-tmp3
               loglik=loglik+log(tmp3) 
            end if   
         end do
         loglik=loglik+real(maxm)*log(2.0) 
    
      end do

      return
      end


c=======================================================================
      subroutine lpostdtfprepoi(bc,ii,
     &                          maxm,maxni,nrec,nsubject,
     &                          ntlr,ntprob,p,
     &                          datastr,y,
     &                          beta,x,ptf,betatf,xtf,
     &                          sigma2b,logpost,k)
c=======================================================================
c     Alejandro Jara, 2008
c=======================================================================
      implicit none
c++++ Input
      integer ii
      integer maxm,maxni,nrec,nsubject,ntlr,ntprob,p,ptf
      integer datastr(nsubject,maxni+1)
      integer y(nrec)
      double precision bc
      double precision beta(p)
      double precision betatf(ntlr,ptf)
      double precision sigma2b
      double precision x(nrec,p+1)  
      double precision xtf(nsubject,ptf)

c++++ Output
      double precision logpost

c++++ External working space
      integer k(maxm)

c++++ Internal working space
      integer i,j,j1,j2,k1,k2,m,ni,ll
      integer kphi
      integer yij
      double precision cdfnorm
      double precision dnrm
      double precision dpoiss
      double precision eta
      double precision tmp1,tmp2,tmp3
      double precision logprior,loglik

c++++ Algorithm

      logprior=0.d0
      loglik=0.d0

      do i=1,maxm
         k(i)=0
      end do

c++++ check if the user has requested an interrupt
      call rchkusr()
  
c++++ Prior contribution

      logprior=dnrm(bc,0.d0,sqrt(sigma2b),1)

      tmp1=bc/sqrt(sigma2b)

      if(tmp1.gt.4.0)then
         tmp2=0.999968
        else if(tmp1.lt.-4.0)then
         tmp2=0.000032
        else 
         tmp2=cdfnorm(bc,0.d0,sqrt(sigma2b),1,0)
      end if

      do j1=1,maxm
         kphi=int(real(2**j1)*tmp2+1)
         k(j1)=kphi
      end do

      j=0
      m=0
      do j1=1,maxm
          if(j1.eq.1)then
             j2=k(j1)
             k2=1
           else 
             j2=j+k(j1)
             k2=m+k(j1-1)
          end if   
          j=j+2**j1
          m=m+2**(j1-1)
          if(j1.eq.1)then
             ll=1
            else
             ll=(k(j1-1)-1)*2+1
          end if

          tmp2=0.d0
          do k1=1,ptf
             tmp2=tmp2+betatf(k2,k1)*xtf(ii,k1)
          end do
          tmp3=exp(tmp2)/(1.0+exp(tmp2))

          if(k(j1).eq.ll)then
             logprior=logprior+log(tmp3) 
           else
             tmp3=1.0-tmp3
             logprior=logprior+log(tmp3) 
          end if   
      end do
      logprior=logprior+real(maxm)*log(2.0) 
   
c++++ check if the user has requested an interrupt
      call rchkusr()
  
c++++ Likelihood contribution

      ni=datastr(ii,1) 

      loglik=0.d0
            
      do j=1,ni
         yij=y(datastr(ii,j+1))            
            
         eta=0.d0 
         do ll=1,p
            eta=eta+x(datastr(ii,j+1),ll)*beta(ll)
         end do
         eta=eta+bc+x(datastr(ii,j+1),p+1)            

         tmp1=dexp(eta)

         loglik=loglik+dpoiss(dble(yij),tmp1,1)
      end do

c++++ Output

      logpost=loglik+logprior

      return
      end




