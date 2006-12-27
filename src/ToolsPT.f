c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR MULTIVARIATES POLYA TREES
c=======================================================================                  
c=======================================================================                  

c=======================================================================                  
      subroutine loglikpt_mucan(m,nrand,nsubject,parti,
     &                          whicho,whichn,b,bzc,cpar,detlogl,
     &                          linf,lsup,muc,sigmainv,
     &                          vec,workmhr,workmr,workvr,
     &                          fixed,loglikc)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     value of the baseline mean in a marginal Multivariate PT.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      real*8 b(nsubject,nrand),bzc(nsubject,nrand),cpar,detlogl
      real*8 linf(nrand),lsup(nrand)
      real*8 muc(nrand),sigmainv(nrand,nrand)
      real*8 vec(nrand)
      real*8 workmhr(nrand*(nrand+1)/2),workmr(nrand,nrand)
      real*8 workvr(nrand)

c-----Output
      real*8 loglikc

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      real*8 dnrm,invcdfnorm,prob,quan,tmp1

c-----Routine

      loglikc=0.d0

      do i=1,nsubject
         do j=1,nrand
            tmp1=0.d0
            do k=1,nrand
               tmp1=tmp1+sigmainv(j,k)*(b(i,k)-muc(k))
            end do
            vec(j)=tmp1
         end do
         
         do j=1,nrand
            bzc(i,j)=vec(j)
         end do

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first subject
         if(i.eq.1)then
            loglikc=-0.5d0*detlogl
            do j=1,nrand
               loglikc=loglikc+dnrm(bzc(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following subjects
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nrand
               if(bzc(i,j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nrand
                  if(bzc(l,j).gt.lsup(j).or.bzc(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do

            if(fixed.ne.1)then 
              loglikc=loglikc+
     &         log((2.d0**nrand)*cpar+dble(2.d0**nrand)*dble(countero))-
     &         log((2.d0**nrand)*cpar+dble(i-1))
            end if 

            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nrand
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(bzc(i,k).le.quan)then
                    parti(k)=k1 
                    lsup(k)=quan
                   else 
                    parti(k)=k2
                    linf(k)=quan
                  end if
               end do                 
               
               countern=0
               do l=1,countero
                  final=1
                  do k=1,nrand
                     if(bzc(whicho(l),k).gt.lsup(k).or.
     &                  bzc(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               loglikc=loglikc+
     &           log((2.d0**nrand)*cpar*dble(je2)+
     &               dble(2.d0**nrand)*dble(countern))-
     &           log((2.d0**nrand)*cpar*dble(je2)+dble(countero))

               if(countern.eq.0)then
                  ok=0
                else  
                  countero=countern
                  do l=1,countern
                     whicho(l)=whichn(l)
                  end do
                  j=j+1
               end if   
            end do

1           continue

            loglikc=loglikc-0.5d0*detlogl
            do j=1,nrand
               loglikc=loglikc+dnrm(bzc(i,j),0.d0, 1.d0, 1)
            end do   
         end if
      end do   

      return
      end


c=======================================================================                  
      subroutine condptprior(vec,ind,nsubject,nrand,bz,cpar,m,detlogl,
     &                       linf,lsup,parti,whicho,whichn,
     &                       fixed,logprior)
c======================================================================= 
c     This subroutine evaluate the log-contional prior distribution,
c     arising in a marginal Multivariate PT, for subject 'ind' with 
c     values 'vec'. The values of random effects 'bz' and 'vec' must 
c     be in a standarized form.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,ind,nrand,nsubject,parti(nrand),m
      integer whicho(nsubject),whichn(nsubject)
      real*8 bz(nsubject,nrand),cpar,detlogl
      real*8 linf(nrand),lsup(nrand)
      real*8 vec(nrand)

c-----Output
      real*8 logprior

c-----Working
      integer countero,countern,final,i,j,je2,k,k1,k2,l,nint,ok
      real*8 dnrm,invcdfnorm,prob,quan

c-----Routine
      
      logprior=0.d0
      
c++++ check if the user has requested an interrupt
      call rchkusr()

      nint=2
      prob=1.d0/dble(nint)
      quan=invcdfnorm(prob,0.d0,1.d0,1,0)

      countero=0
             
      do i=1,nrand
         if(vec(i).le.quan)then
            linf(i)=-999999.d0
            lsup(i)=quan
            parti(i)=1
          else
            linf(i)=quan
            lsup(i)= 999999.d0
            parti(i)=2
         end if
      end do      

      do i=1,nsubject
         final=1
         if(i.ne.ind)then
         do j=1,nrand
            if(bz(i,j).gt.lsup(j).or.bz(i,j).lt.linf(j))then
              final=0
            end if
         end do
      
         if(final.eq.1)then
            countero=countero+1
            whicho(countero)=i
         end if   
         end if
      end do

      if(fixed.ne.1)then
         logprior=logprior+
     &     log((2.d0**nrand)*cpar+dble((2.d0**nrand)*countero))-
     &     log((2.d0**nrand)*cpar+dble(nsubject-1))
      end if

      if(countero.eq.0) go to 1

      ok=1
      j=2
      do while(ok.eq.1.and.j.le.m)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)

         do k=1,nrand
            k1=2*(parti(k)-1)+1
            k2=2*(parti(k)-1)+2
            quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
         
            if(vec(k).le.quan)then
              parti(k)=k1 
              lsup(k)=quan
             else 
              parti(k)=k2
              linf(k)=quan
            end if
         end do                 
               
         countern=0
         do l=1,countero
            final=1
            do k=1,nrand
               if(bz(whicho(l),k).gt.lsup(k).or.
     &            bz(whicho(l),k).lt.linf(k)    )then
                  final=0 
               end if   
            end do
         
            if(final.eq.1)then
              countern=countern+1
              whichn(countern)=whicho(l)
            end if
         end do

         logprior=logprior+
     &           log((2.d0**nrand)*cpar*dble(je2)+
     &               dble(2.d0**nrand)*dble(countern))-
     &           log((2.d0**nrand)*cpar*dble(je2)+dble(countero))

         if(countern.eq.0)then
            ok=0
          else  
            countero=countern
            do l=1,countern
               whicho(l)=whichn(l)
            end do
            j=j+1
         end if   
      end do

1     continue

      logprior=logprior-0.5d0*detlogl
      do j=1,nrand
         logprior=logprior+dnrm(vec(j),0.d0, 1.d0, 1)
      end do   

      return
      end


c=======================================================================                  
      subroutine loglikpt_update(ind,m,nrand,nsubject,parti,
     &                           whicho,whichn,bz,cpar,detlogl,
     &                           linf,lsup,
     &                           fixed,logliko)
c======================================================================= 
c     This subroutine evaluate sequentially the log-likelihood for 
c     the current value of the baseline parameters in a marginal 
c     Multivariate PT.
c     This sequential update is performed as the random effects are
c     updated.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,ind,m,nrand,nsubject
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      real*8 bz(nsubject,nrand),cpar,detlogl
      real*8 linf(nrand),lsup(nrand)

c-----Input/Output
      real*8 logliko

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      real*8 dnrm,invcdfnorm,prob,quan

c-----Routine

      do i=1,ind

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first subject
         if(i.eq.1)then
            logliko=-0.5d0*detlogl
            do j=1,nrand
               logliko=logliko+dnrm(bz(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following subjects
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nrand
               if(bz(i,j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nrand
                  if(bz(l,j).gt.lsup(j).or.bz(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do
            
            if(fixed.ne.1)then
              logliko=logliko+
     &         log((2.d0**nrand)*cpar+dble(2.d0**nrand)*dble(countero))-
     &         log((2.d0**nrand)*cpar+dble(i-1))
            end if 
            
            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nrand
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(bz(i,k).le.quan)then
                    parti(k)=k1 
                    lsup(k)=quan
                   else 
                    parti(k)=k2
                    linf(k)=quan
                  end if
               end do                 
               
               countern=0
               do l=1,countero
                  final=1
                  do k=1,nrand
                     if(bz(whicho(l),k).gt.lsup(k).or.
     &                  bz(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               logliko=logliko+
     &           log((2.d0**nrand)*cpar*dble(je2)+
     &               dble(2.d0**nrand)*dble(countern))-
     &           log((2.d0**nrand)*cpar*dble(je2)+dble(countero))


               if(countern.eq.0)then
                  ok=0
                else  
                  countero=countern
                  do l=1,countern
                     whicho(l)=whichn(l)
                  end do
                  j=j+1
               end if   
            end do

1           continue

            logliko=logliko-0.5d0*detlogl
            do j=1,nrand
               logliko=logliko+dnrm(bz(i,j),0.d0, 1.d0, 1)
            end do   
         end if
      end do   

      return
      end

c=======================================================================                  
      subroutine loglikpt_updatet(m,nrand,nsubject,parti,
     &                            whicho,whichn,bz,cpar,detlogl,
     &                            linf,lsup,
     &                            fixed,logliko)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for 
c     the current value of the baseline parameters in a marginal 
c     Multivariate PT.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      real*8 bz(nsubject,nrand),cpar,detlogl
      real*8 linf(nrand),lsup(nrand)

c-----Output
      real*8 logliko

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      real*8 dnrm,invcdfnorm,prob,quan

c-----Routine

      logliko=0.d0

      do i=1,nsubject

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first subject
         if(i.eq.1)then
            logliko=-0.5d0*detlogl
            do j=1,nrand
               logliko=logliko+dnrm(bz(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following subjects
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nrand
               if(bz(i,j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nrand
                  if(bz(l,j).gt.lsup(j).or.bz(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do
            
            if(fixed.ne.1)then
              logliko=logliko+
     &         log((2.d0**nrand)*cpar+dble(2.d0**nrand)*dble(countero))-
     &         log((2.d0**nrand)*cpar+dble(i-1))
            end if
            
            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nrand
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(bz(i,k).le.quan)then
                    parti(k)=k1 
                    lsup(k)=quan
                   else 
                    parti(k)=k2
                    linf(k)=quan
                  end if
               end do                 
               
               countern=0
               do l=1,countero
                  final=1
                  do k=1,nrand
                     if(bz(whicho(l),k).gt.lsup(k).or.
     &                  bz(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               logliko=logliko+
     &           log((2.d0**nrand)*cpar*dble(je2)+
     &               dble(2.d0**nrand)*dble(countern))-
     &           log((2.d0**nrand)*cpar*dble(je2)+dble(countero))

               if(countern.eq.0)then
                  ok=0
                else  
                  countero=countern
                  do l=1,countern
                     whicho(l)=whichn(l)
                  end do
                  j=j+1
               end if   
            end do

1           continue

            logliko=logliko-0.5d0*detlogl
            do j=1,nrand
               logliko=logliko+dnrm(bz(i,j),0.d0, 1.d0, 1)
            end do   
         end if
      end do   

      return
      end

c=======================================================================                  
      subroutine loglikpt_covarcan(m,nrand,nsubject,iflagr,parti,
     &                             whicho,whichn,b,bzc,cpar,detloglc,
     &                             linf,lsup,mu,sigmac,sigmainvc,
     &                             vec,workmhr,workmr,workvr,
     &                             loglikc,typep,workmr1,workmr2,fixed)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     of the baseline covariance matrix in a marginal Multivariate PT.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject,typep
      integer iflagr(nrand)
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      real*8 b(nsubject,nrand),bzc(nsubject,nrand),cpar,detloglc
      real*8 linf(nrand),lsup(nrand)
      real*8 mu(nrand),sigmac(nrand,nrand),sigmainvc(nrand,nrand)
      real*8 vec(nrand)
      real*8 workmhr(nrand*(nrand+1)/2),workmr(nrand,nrand)
      real*8 workmr1(nrand,nrand),workmr2(nrand,nrand)
      real*8 workvr(nrand)

c-----Output
      real*8 loglikc

c-----Working
      integer countero,countern,final
      integer i,ihmssf,j,je2,k,k1,k2,l,nint,ok
      real*8 dnrm,invcdfnorm,prob,quan,tmp1

c-----Routine

      loglikc=0.d0

      do i=1,nrand
         do j=1,nrand
            workmr(i,j)=sigmac(i,j)
         end do
      end do

      call invdet(workmr,nrand,sigmainvc,detloglc,iflagr,workvr)

      if(typep.eq.1)then
         do i=1,nrand
            do j=1,nrand
               workmr(i,j)=0.d0
               sigmainvc(i,j)=0.d0
            end do
         end do
         call cholesky(nrand,sigmac,workmhr)
         do i=1,nrand
            do j=1,i
               sigmainvc(i,j)=workmhr(ihmssf(i,j,nrand))
            end do
         end do
         call inverse(sigmainvc,nrand,iflagr)      
         
       else if(typep.eq.2)then

         do i=1,nrand
            do j=1,nrand
               workmr1(i,j)=0.d0
               sigmainvc(i,j)=0.d0
            end do
         end do
         call eigenv(nrand,nrand,sigmac,workvr,linf,workmr)

         do i=1,nrand
            workmr1(i,i)=sqrt(workvr(i))
         end do
         
         do i=1,nrand
            do j=1,nrand
               tmp1=0.d0
               do k=1,nrand
                  tmp1=tmp1+workmr(i,k)*workmr1(k,j)
               end do
               sigmainvc(i,j)=tmp1
            end do
         end do
         call inverse(sigmainvc,nrand,iflagr)      
      
       else
         do i=1,nrand
            do j=1,nrand
               workmr1(i,j)=0.d0
               workmr2(i,j)=0.d0
               sigmainvc(i,j)=0.d0
            end do
         end do
         call eigenv(nrand,nrand,sigmac,workvr,linf,workmr)

         do i=1,nrand
            workmr1(i,i)=sqrt(workvr(i))
         end do

         do i=1,nrand
            do j=1,nrand
               tmp1=0.d0
               do k=1,nrand
                  tmp1=tmp1+workmr(i,k)*workmr1(k,j)
               end do
               workmr2(i,j)=tmp1
            end do
         end do

         do i=1,nrand
            do j=1,nrand
               tmp1=0.d0
               do k=1,nrand
                  tmp1=tmp1+workmr2(i,k)*workmr(j,k)
               end do
               sigmainvc(i,j)=tmp1
            end do
         end do
         call inverse(sigmainvc,nrand,iflagr)      
      end if 
      
      do i=1,nsubject
         do j=1,nrand
            tmp1=0.d0
            do k=1,nrand
               tmp1=tmp1+sigmainvc(j,k)*(b(i,k)-mu(k))
            end do
            vec(j)=tmp1
         end do
         
         do j=1,nrand
            bzc(i,j)=vec(j)
         end do

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first subject
         if(i.eq.1)then
            loglikc=-0.5d0*detloglc
            do j=1,nrand
               loglikc=loglikc+dnrm(bzc(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following subjects
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nrand
               if(bzc(i,j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nrand
                  if(bzc(l,j).gt.lsup(j).or.bzc(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do
            
            if(fixed.ne.1)then
              loglikc=loglikc+
     &         log((2.d0**nrand)*cpar+(2.d0**nrand)*dble(countero))-
     &         log((2.d0**nrand)*cpar+dble(i-1))
            end if 
            
            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nrand
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(bzc(i,k).le.quan)then
                    parti(k)=k1 
                    lsup(k)=quan
                   else 
                    parti(k)=k2
                    linf(k)=quan
                  end if
               end do                 
               
               countern=0
               do l=1,countero
                  final=1
                  do k=1,nrand
                     if(bzc(whicho(l),k).gt.lsup(k).or.
     &                  bzc(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               loglikc=loglikc+
     &         log((2.d0**nrand)*cpar*dble(je2)+
     &             (2.d0**nrand)*dble(countern))-
     &         log((2.d0**nrand)*cpar*dble(je2)+dble(countero))

               if(countern.eq.0)then
                  ok=0
                else  
                  countero=countern
                  do l=1,countern
                     whicho(l)=whichn(l)
                  end do
                  j=j+1
               end if   
            end do

1           continue

            loglikc=loglikc-0.5d0*detloglc
            do j=1,nrand
               loglikc=loglikc+dnrm(bzc(i,j),0.d0, 1.d0, 1)
            end do   
         end if
      end do   

      return
      end

c=======================================================================                  
      subroutine loglikpt_cparcan(m,nrand,nsubject,iflagr,parti,
     &                            whicho,whichn,bz,cparc,detlogl,
     &                            linf,lsup,
     &                            vec,workmhr,workmr,workvr,
     &                            fixed,loglikn)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     value of the precision parameter in a marginal Multivariate PT.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer iflagr(nrand)
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      real*8 bz(nsubject,nrand),cparc,detlogl
      real*8 linf(nrand),lsup(nrand)
      real*8 vec(nrand)
      real*8 workmhr(nrand*(nrand+1)/2),workmr(nrand,nrand)
      real*8 workvr(nrand)

c-----Output
      real*8 loglikn

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      real*8 dnrm,invcdfnorm,prob,quan

c-----Routine

      loglikn=0.d0

      do i=1,nsubject

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first subject
         if(i.eq.1)then
            loglikn=-0.5d0*detlogl
            do j=1,nrand
               loglikn=loglikn+dnrm(bz(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following subjects
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nrand
               if(bz(i,j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nrand
                  if(bz(l,j).gt.lsup(j).or.bz(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do
            
            if(fixed.ne.1)then
              loglikn=loglikn+
     &         log((2.d0**nrand)*cparc+dble((2.d0**nrand)*countero))-
     &         log((2.d0**nrand)*cparc+dble(i-1))
            end if
            
            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nrand
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(bz(i,k).le.quan)then
                    parti(k)=k1 
                    lsup(k)=quan
                   else 
                    parti(k)=k2
                    linf(k)=quan
                  end if
               end do                 
               
               countern=0
               do l=1,countero
                  final=1
                  do k=1,nrand
                     if(bz(whicho(l),k).gt.lsup(k).or.
     &                  bz(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               loglikn=loglikn+
     &           log((2.d0**nrand)*cparc*dble(je2)+
     &               dble(2.d0**nrand)*dble(countern))-
     &           log((2.d0**nrand)*cparc*dble(je2)+dble(countero))


               if(countern.eq.0)then
                  ok=0
                else  
                  countero=countern
                  do l=1,countern
                     whicho(l)=whichn(l)
                  end do
                  j=j+1
               end if   
            end do

1           continue

            loglikn=loglikn-0.5d0*detlogl
            do j=1,nrand
               loglikn=loglikn+dnrm(bz(i,j),0.d0, 1.d0, 1)
            end do   
         end if
      end do   

      return
      end

c=======================================================================                  
      subroutine sampredpt(marea,nrand,nsubject,parti,m,mass,massi,
     &                     pattern,patterns,whichn,whicho,bz,mu,sigma,
     &                     cpar,limw,linf,lsup,workmr,workmhr,
     &                     workvr,vec,fixed)
c======================================================================= 
c     This subroutine generates a sample 'vec' from the predictive
c     distribution arising in a marginal Multivariate PT.
c     The values of random effects 'bz' must be in a standarized form.
c     The output 'vec' is in a normal form. It uses the Cholesky
c     transformation.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,marea,narea,nrand,nsubject,parti(nrand)
      integer m,massi(marea),pattern(nrand),patterns(nrand)
      integer whicho(nsubject),whichn(nsubject)      
      real*8 bz(nsubject,nrand),cpar
      real*8 limw(nrand),linf(nrand),lsup(nrand) 
      real*8 mu(nrand),sigma(nrand,nrand)
      real*8 mass(marea),rtnorm

      real*8 workmhr(nrand*(nrand+1)/2),workmr(nrand,nrand)
      real*8 workvr(nrand)
      
c-----Output
      real*8 vec(nrand)

c-----Working
      integer binaryrep,countero,countern,evali,evali2,final
      integer i,ihmssf,j,je2,k,k1,k2,l,nint,ok 
      real*8 invcdfnorm,prob,quan,tmp1

c-----Routine
        
      narea=2**nrand
  
      nint=2
      prob=1.d0/dble(nint)
      quan=invcdfnorm(prob,0.d0,1.d0,1,0)

      countero=0
      
      do i=1,narea
         massi(i)=0
         mass(i)=0.d0
      end do
      
      do i=1,nsubject
         do j=1,nrand
            evali=1 
            if(bz(i,j).le.quan)evali=0  
            pattern(j)=evali
         end do
         evali=binaryrep(nrand,pattern)
         massi(evali)=massi(evali)+1
      end do   

      if(fixed.ne.1)then 
         do l=1,narea
            mass(l)=(cpar+dble(massi(l)))/
     &              ((2**nrand)*cpar+dble(nsubject))
         end do
       else
         do l=1,narea
            mass(l)=1.d0/dble(narea)
         end do
      end if 

      call simdisc(mass,marea,narea,evali)  
      
      evali2=evali
      call binaryrepinv(nrand,evali2,patterns)
       
      do i=1,nsubject
         final=1
         do j=1,nrand
            evali=1 
            if(bz(i,j).le.quan)evali=0  
            pattern(j)=evali
            if(pattern(j).ne.patterns(j))final=0
         end do
      
         if(final.eq.1)then
            countero=countero+1
            whicho(countero)=i
         end if   
      end do 
      
      do i=1,nrand
         if(patterns(i).eq.0)then
           linf(i)=-999999.d0
           lsup(i)=quan
           parti(i)=1
          else
           linf(i)=quan
           lsup(i)= 999999.d0
           parti(i)=2
         end if 
      end do

      if(countero.eq.0) go to 1

      ok=1
      j=2
      countern=0
      
      do while(ok.eq.1.and.j.le.m)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
         
         do k=1,nrand
            
            k1=2*(parti(k)-1)+1
            k2=2*(parti(k)-1)+2
            quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)

            limw(k)=quan
            
            if(quan.gt.lsup(k).or.quan.lt.linf(k))then
               call intpr("j",-1,j,1)
               call dblepr("linf",-1,linf,nrand)
               call dblepr("lsup",-1,lsup,nrand)
               call dblepr("limw",-1,limw,nrand)
               call rexit("Errors in limits")
            end if
         end do   

         do k=1,narea
            massi(k)=0
            mass(k)=0.d0
         end do
         
         do l=1,countero
            do k=1,nrand
               evali=1 
               if(bz(whicho(l),k).le.limw(k))evali=0  
               pattern(k)=evali
            end do
            evali=binaryrep(nrand,pattern)
            massi(evali)=massi(evali)+1
         end do                      

         do l=1,narea
            mass(l)=(cpar*dble(je2)+dble(massi(l)))/
     &          ((2**nrand)*cpar*dble(je2)+dble(countero))
         end do

         call simdisc(mass,marea,narea,evali)  
         
         evali2=evali
         call binaryrepinv(nrand,evali2,patterns)

         countern=0
         do l=1,countero
            final=1
            do k=1,nrand
               evali=1 
               if(bz(whicho(l),k).le.limw(k))evali=0  
               pattern(k)=evali
               if(pattern(k).ne.patterns(k))final=0
            end do
            if(final.eq.1)then
               countern=countern+1
               whichn(countern)=whicho(l)
            end if   
         end do  

         do k=1,nrand
            if(patterns(k).eq.0)then
              parti(k)=2*(parti(k)-1)+1
              lsup(k)=limw(k)
             else
              parti(k)=2*(parti(k)-1)+2
              linf(k)=limw(k)
            end if 
         end do

         if(countern.eq.0)then
            ok=0
          else  
            countero=countern
            do l=1,countern
               whicho(l)=whichn(l)
            end do
            j=j+1
         end if   
      end do

1     continue 

      do k=1,nrand
         if(linf(k).ge.lsup(k))then
            call intpr("k",-1,k,1)
            call dblepr("linf",-1,linf,nrand)
            call dblepr("lsup",-1,lsup,nrand)
            call rexit("Errors in limits")
         end if
      end do   

      do i=1,nrand
         workvr(i)=rtnorm(0.d0,1.d0,linf(i),lsup(i),
     &                   .false.,.false.)
      end do


      do i=1,nrand
         do j=1,nrand
            workmr(i,j)=sigma(i,j)
         end do
      end do
        
      call cholesky(nrand,workmr,workmhr)

      do i=1,nrand
         do j=1,nrand
            workmr(i,j)=0.d0
         end do
      end do
       
      do i=1,nrand
         do j=1,i
            workmr(i,j)=workmhr(ihmssf(i,j,nrand))
         end do
      end do

      do i=1,nrand
         tmp1=0.d0
         do j=1,nrand
            tmp1=tmp1+workmr(i,j)*workvr(j)   
         end do
         vec(i)=tmp1+mu(i)
      end do

      return
      end

c=======================================================================                  
      subroutine sampredptun(marea,nrand,nsubject,parti,m,mass,massi,
     &                       pattern,patterns,whichn,whicho,bz,
     &                       cpar,limw,linf,lsup,workmr,workmhr,
     &                       workvr,vec,fixed)
c======================================================================= 
c     This subroutine generates a sample 'vec' from the predictive
c     distribution arising in a marginal Multivariate PT.
c     The values of random effects 'bz' must be in a standarized form.
c     The output 'vec' is also standarized.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,marea,narea,nrand,nsubject,parti(nrand)
      integer m,massi(marea),pattern(nrand),patterns(nrand)
      integer whicho(nsubject),whichn(nsubject)      
      real*8 bz(nsubject,nrand),cpar
      real*8 limw(nrand),linf(nrand),lsup(nrand) 
      real*8 mass(marea),rtnorm

      real*8 workmhr(nrand*(nrand+1)/2),workmr(nrand,nrand)
      real*8 workvr(nrand)
      
c-----Output
      real*8 vec(nrand)

c-----Working
      integer binaryrep,countero,countern,evali,evali2,final
      integer i,j,je2,k,k1,k2,l,nint,ok 
      real*8 invcdfnorm,prob,quan

c-----Routine
        
      narea=2**nrand
  
      nint=2
      prob=1.d0/dble(nint)
      quan=invcdfnorm(prob,0.d0,1.d0,1,0)

      countero=0
      
      do i=1,narea
         massi(i)=0
         mass(i)=0.d0
      end do
      
      do i=1,nsubject
         do j=1,nrand
            evali=1 
            if(bz(i,j).le.quan)evali=0  
            pattern(j)=evali
         end do
         evali=binaryrep(nrand,pattern)
         massi(evali)=massi(evali)+1
      end do   

      if(fixed.ne.1)then 
         do l=1,narea
            mass(l)=(cpar+dble(massi(l)))/
     &              ((2**nrand)*cpar+dble(nsubject))
         end do
       else
         do l=1,narea
            mass(l)=1.d0/dble(narea)
         end do
      end if 

      call simdisc(mass,marea,narea,evali)  
      
      evali2=evali
      call binaryrepinv(nrand,evali2,patterns)
       
      do i=1,nsubject
         final=1
         do j=1,nrand
            evali=1 
            if(bz(i,j).le.quan)evali=0  
            pattern(j)=evali
            if(pattern(j).ne.patterns(j))final=0
         end do
      
         if(final.eq.1)then
            countero=countero+1
            whicho(countero)=i
         end if   
      end do 
      
      do i=1,nrand
         if(patterns(i).eq.0)then
           linf(i)=-999999.d0
           lsup(i)=quan
           parti(i)=1
          else
           linf(i)=quan
           lsup(i)= 999999.d0
           parti(i)=2
         end if 
      end do

      if(countero.eq.0) go to 1

      ok=1
      j=2
      countern=0
      
      do while(ok.eq.1.and.j.le.m)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
         
         do k=1,nrand
            
            k1=2*(parti(k)-1)+1
            k2=2*(parti(k)-1)+2
            quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)

            limw(k)=quan
            
            if(quan.gt.lsup(k).or.quan.lt.linf(k))then
               call intpr("j",-1,j,1)
               call dblepr("linf",-1,linf,nrand)
               call dblepr("lsup",-1,lsup,nrand)
               call dblepr("limw",-1,limw,nrand)
               call rexit("Errors in limits")
            end if
         end do   

         do k=1,narea
            massi(k)=0
            mass(k)=0.d0
         end do
         
         do l=1,countero
            do k=1,nrand
               evali=1 
               if(bz(whicho(l),k).le.limw(k))evali=0  
               pattern(k)=evali
            end do
            evali=binaryrep(nrand,pattern)
            massi(evali)=massi(evali)+1
         end do                      

         do l=1,narea
            mass(l)=(cpar*dble(je2)+dble(massi(l)))/
     &          ((2**nrand)*cpar*dble(je2)+dble(countero))
         end do

         call simdisc(mass,marea,narea,evali)  
         
         evali2=evali
         call binaryrepinv(nrand,evali2,patterns)

         countern=0
         do l=1,countero
            final=1
            do k=1,nrand
               evali=1 
               if(bz(whicho(l),k).le.limw(k))evali=0  
               pattern(k)=evali
               if(pattern(k).ne.patterns(k))final=0
            end do
            if(final.eq.1)then
               countern=countern+1
               whichn(countern)=whicho(l)
            end if   
         end do  

         do k=1,nrand
            if(patterns(k).eq.0)then
              parti(k)=2*(parti(k)-1)+1
              lsup(k)=limw(k)
             else
              parti(k)=2*(parti(k)-1)+2
              linf(k)=limw(k)
            end if 
         end do

         if(countern.eq.0)then
            ok=0
          else  
            countero=countern
            do l=1,countern
               whicho(l)=whichn(l)
            end do
            j=j+1
         end if   
      end do

1     continue 

      do k=1,nrand
         if(linf(k).ge.lsup(k))then
            call intpr("k",-1,k,1)
            call dblepr("linf",-1,linf,nrand)
            call dblepr("lsup",-1,lsup,nrand)
            call rexit("Errors in limits")
         end if
      end do   

      do i=1,nrand
         vec(i)=rtnorm(0.d0,1.d0,linf(i),lsup(i),
     &                   .false.,.false.)
      end do

      return
      end


c=======================================================================                  
      subroutine loglikpt_cur(m,nrand,nsubject,parti,
     &                        whicho,whichn,bz,cpar,detlogl,
     &                        linf,lsup,
     &                        vec,workmhr,workmr,workvr,
     &                        fixed,logliko)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the current value
c     of the baseline parameters in a marginal Multivariate PT.
c     This function does not include the standarization of random 
c     effects
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      real*8 bz(nsubject,nrand),cpar,detlogl
      real*8 linf(nrand),lsup(nrand)
      real*8 vec(nrand)
      real*8 workmhr(nrand*(nrand+1)/2),workmr(nrand,nrand)
      real*8 workvr(nrand)

c-----Output
      real*8 logliko

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      real*8 dnrm,invcdfnorm,prob,quan

c-----Routine

      logliko=0.d0

      do i=1,nsubject

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++ first subject
         if(i.eq.1)then
            logliko=-0.5d0*detlogl
            do j=1,nrand
               logliko=logliko+dnrm(bz(i,j),0.d0, 1.d0, 1)
            end do   

c+++++++ following subjects
          else

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,nrand
               if(bz(i,j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do
           
            do l=1,i-1
               final=1
               do j=1,nrand
                  if(bz(l,j).gt.lsup(j).or.bz(l,j).lt.linf(j))then
                    final=0
                  end if
               end do
               
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do
            
            if(fixed.ne.1)then
              logliko=logliko+
     &         log((2.d0**nrand)*cpar+dble(2.d0**nrand)*dble(countero))-
     &         log((2.d0**nrand)*cpar+dble(i-1))
            end if
            
            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,nrand
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                  if(bz(i,k).le.quan)then
                    parti(k)=k1 
                    lsup(k)=quan
                   else 
                    parti(k)=k2
                    linf(k)=quan
                  end if
               end do                 
               
               countern=0
               do l=1,countero
                  final=1
                  do k=1,nrand
                     if(bz(whicho(l),k).gt.lsup(k).or.
     &                  bz(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
                  
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               logliko=logliko+
     &           log((2.d0**nrand)*cpar*dble(je2)+
     &               dble(2.d0**nrand)*dble(countern))-
     &           log((2.d0**nrand)*cpar*dble(je2)+dble(countero))


               if(countern.eq.0)then
                  ok=0
                else  
                  countero=countern
                  do l=1,countern
                     whicho(l)=whichn(l)
                  end do
                  j=j+1
               end if   
            end do

1           continue

            logliko=logliko-0.5d0*detlogl
            do j=1,nrand
               logliko=logliko+dnrm(bz(i,j),0.d0, 1.d0, 1)
            end do   
         end if
      end do   

      return
      end


c=======================================================================
      subroutine predictiveptb(m,nsubject,q,nsave,randsave,mumat,
     &                         sigmamat,cparvec,typep,
     &                         ngrid1,ngrid2,grid1,grid2,fs,
     &                         iflagr,parti,whicho,whichn,
     &                         b,bz,linf,lsup,mu,sigma,sigmainv,
     &                         theta,thetaz,workmr,workmr1,workmr2,
     &                         workmhr,workvr,workvr1,fixed)
c=======================================================================
c     computes the bivariate posterior predictive density from the
c     output of a PTfunction. This is used for random effects models.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c=======================================================================
      implicit none

c++++ input
      integer fixed,m,ngrid1,ngrid2,nsubject,nsave,q,typep
      real*8 cparvec(nsave),randsave(nsave,q*(nsubject+1))
      real*8 mumat(nsave,q)
      real*8 sigmamat(nsave,q*(q+1)/2)
      real*8 grid1(ngrid1),grid2(ngrid2)

c++++ output
      real*8 fs(ngrid1,ngrid2)

c++++ external working space
      integer iflagr(q),parti(q)
      integer whicho(nsubject),whichn(nsubject)
      real*8 b(nsubject,q),bz(nsubject,q)
      real*8 linf(q),lsup(q),mu(q),sigma(q,q),sigmainv(q,q)
      real*8 theta(q),thetaz(q)
      real*8 workmr(q,q),workmr1(q,q),workmr2(q,q)
      real*8 workmhr(q*(q+1)/2)
      real*8 workvr(q),workvr1(q)

c++++ internal working space
      integer countero,countern
      integer final,i,ii,ihmssf,j,jj,je2,k,kk,k1,k2,l
      integer nint,ok
      real*8 cpar,detlogl,dnrm
      real*8 invcdfnorm
      real*8 loglik,prob,quan,tmp1

c++++ algorithm      
      if(q.gt.2)then
        call rexit("Only bivariate evaluation supported")
      end if

      do ii=1,nsave

c+++++++ save elements

c+++++++ check if the user has requested an interrupt
         call rchkusr()
          
         cpar=cparvec(ii)
         
         do i=1,q
            mu(i)=mumat(ii,i)
            do j=1,q
               sigma(i,j)=sigmamat(ii,ihmssf(i,j,q))
            end do
         end do
         
         k=0
         do i=1,nsubject
            do j=1,q
               k=k+1
               b(i,j)=randsave(ii,k)
            end do
         end do
      
c+++++++ covariance matrix decomposition

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,q
            do j=1,q
               workmr(i,j)=sigma(i,j)
            end do
         end do
         call inversedet(workmr,q,iflagr,detlogl)

         if(typep.eq.1)then
            do i=1,q
               do j=1,q
                  workmr(i,j)=0.d0
                  sigmainv(i,j)=0.d0
               end do
            end do
            call cholesky(q,sigma,workmhr)
            do i=1,q
               do j=1,i
                  sigmainv(i,j)=workmhr(ihmssf(i,j,q))
               end do
            end do
            call inverse(sigmainv,q,iflagr)      
         
          else if(typep.eq.2)then
            do i=1,q
               do j=1,q
                  workmr(i,j)=0.d0
                  workmr1(i,j)=0.d0
                  sigmainv(i,j)=0.d0
               end do
            end do
            call eigenv(q,q,sigma,workvr,workvr1,workmr)
            do i=1,q
               workmr1(i,i)=sqrt(workvr(i))
            end do
         
            do i=1,q
               do j=1,q
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+workmr(i,k)*workmr1(k,j)
                  end do
                  sigmainv(i,j)=tmp1
               end do
            end do
            call inverse(sigmainv,q,iflagr)      
      
          else
            do i=1,q
               do j=1,q
                  workmr(i,j)=0.d0
                  workmr1(i,j)=0.d0
                  workmr2(i,j)=0.d0
                  sigmainv(i,j)=0.d0
               end do
            end do
            call eigenv(q,q,sigma,workvr,workvr1,workmr)
            do i=1,q
               workmr1(i,i)=sqrt(workvr(i))
            end do
 
            do i=1,q
               do j=1,q
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+workmr(i,k)*workmr1(k,j)
                  end do
                  workmr2(i,j)=tmp1
               end do
            end do

            do i=1,q
               do j=1,q
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+workmr2(i,k)*workmr(j,k)
                  end do
                  sigmainv(i,j)=tmp1
               end do
            end do
            call inverse(sigmainv,q,iflagr)      
         end if 

c+++++++ transformation of the random effects

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,nsubject
            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+sigmainv(j,k)*(b(i,k)-mu(k))
               end do
               bz(i,j)=tmp1
            end do
         end do  

c+++++++ evaluate the grid
         
         do jj=1,ngrid1
            theta(1)=grid1(jj)
            do kk=1,ngrid2
               theta(2)=grid2(kk)
            
c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()
               
               loglik=0.d0
               
               do i=1,q
                  tmp1=0.d0
                  do j=1,q
                     tmp1=tmp1+sigmainv(i,j)*(theta(j)-mu(j))
                  end do
                  thetaz(i)=tmp1
               end do

               nint=2
               prob=1.d0/dble(nint)
               quan=invcdfnorm(prob,0.d0,1.d0,1,0)

               countero=0
            
               do j=1,q
                  if(thetaz(j).le.quan)then
                     linf(j)=-999999.d0
                     lsup(j)=quan
                     parti(j)=1
                   else
                     linf(j)=quan
                     lsup(j)= 999999.d0
                     parti(j)=2
                  end if
               end do

               do l=1,nsubject
                  final=1
                  do j=1,q
                     if(bz(l,j).gt.lsup(j).or.bz(l,j).lt.linf(j))then
                        final=0
                     end if
                  end do
               
                  if(final.eq.1)then
                     countero=countero+1
                     whicho(countero)=l
                  end if   
               end do
               
               if(fixed.ne.1)then
                 loglik=loglik+
     &           log((2.d0**q)*cpar+dble(2.d0**q)*dble(countero))-
     &           log((2.d0**q)*cpar+dble(nsubject))
               end if 

               if(countero.eq.0) go to 1

               ok=1
               j=2
               do while(ok.eq.1.and.j.le.m)
                  nint=2**j
                  je2=j**2
                  prob=1.d0/dble(nint)

                  do k=1,q
                     k1=2*(parti(k)-1)+1
                     k2=2*(parti(k)-1)+2
                     quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
               
                     if(thetaz(k).le.quan)then
                       parti(k)=k1 
                       lsup(k)=quan
                      else 
                       parti(k)=k2
                       linf(k)=quan
                     end if
                  end do                 
               
                  countern=0
                  do l=1,countero
                     final=1
                     do k=1,q
                        if(bz(whicho(l),k).gt.lsup(k).or.
     &                     bz(whicho(l),k).lt.linf(k)    )then
                           final=0 
                        end if   
                     end do
                  
                     if(final.eq.1)then
                       countern=countern+1
                       whichn(countern)=whicho(l)
                     end if
                  end do

                  loglik=loglik+
     &               log((2.d0**q)*cpar*dble(je2)+
     &                   dble(2.d0**q)*dble(countern))-
     &               log((2.d0**q)*cpar*dble(je2)+dble(countero))

                  if(countern.eq.0)then
                     ok=0
                   else  
                     countero=countern
                     do l=1,countern
                        whicho(l)=whichn(l)
                     end do
                     j=j+1
                  end if   
               end do

1              continue

               loglik=loglik-0.5d0*detlogl
               do j=1,q
                  loglik=loglik+dnrm(thetaz(j),0.d0, 1.d0, 1)
               end do   
               fs(jj,kk)=fs(jj,kk)+exp(loglik)            
            end do
         end do
      end do

      do i=1,ngrid1
         do j=1,ngrid2
            fs(i,j)=fs(i,j)/dble(nsave)
         end do
      end do

      return
      end
      
c=======================================================================
      subroutine predictiveptu(m,nsubject,q,nsave,randsave,mumat,
     &                         sigmamat,cparvec,typep,
     &                         ngrid,grid,fs,
     &                         iflagr,parti,whicho,whichn,
     &                         b,bz,linf,lsup,mu,sigma,sigmainv,
     &                         theta,thetaz,workmr,workmr1,workmr2,
     &                         workmhr,workvr,workvr1,fixed)
c=======================================================================
c     computes the univariate posterior predictive density from the
c     output of a PTfunction. This is used for random effects models.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006
c     Last modification: 24-04-2007.
c=======================================================================
      implicit none

c++++ input
      integer fixed,m,ngrid,nsubject,nsave,q,typep
      real*8 cparvec(nsave),randsave(nsave,q*(nsubject+1))
      real*8 mumat(nsave,q)
      real*8 sigmamat(nsave,q*(q+1)/2)
      real*8 grid(ngrid)

c++++ output
      real*8 fs(ngrid)

c++++ external working space
      integer iflagr(q),parti(q)
      integer whicho(nsubject),whichn(nsubject)
      real*8 b(nsubject,q),bz(nsubject,q)
      real*8 linf(q),lsup(q),mu(q),sigma(q,q),sigmainv(q,q)
      real*8 theta(q),thetaz(q)
      real*8 workmr(q,q),workmr1(q,q),workmr2(q,q)
      real*8 workmhr(q*(q+1)/2)
      real*8 workvr(q),workvr1(q)

c++++ internal working space
      integer countero,countern
      integer final,i,ii,ihmssf,j,jj,je2,k,k1,k2,l
      integer nint,ok
      real*8 cpar,detlogl,dnrm
      real*8 invcdfnorm
      real*8 loglik,prob,quan,tmp1

c++++ algorithm      
      if(q.gt.1)then
        call rexit("Only univariate evaluation supported")
      end if

      do ii=1,nsave

c+++++++ save elements

c+++++++ check if the user has requested an interrupt
         call rchkusr()
          
         cpar=cparvec(ii)
         
         do i=1,q
            mu(i)=mumat(ii,i)
            do j=1,q
               sigma(i,j)=sigmamat(ii,ihmssf(i,j,q))
            end do
         end do
         
         k=0
         do i=1,nsubject
            do j=1,q
               k=k+1
               b(i,j)=randsave(ii,k)
            end do
         end do
      
c+++++++ covariance matrix decomposition

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,q
            do j=1,q
               workmr(i,j)=sigma(i,j)
            end do
         end do
         call inversedet(workmr,q,iflagr,detlogl)

         if(typep.eq.1)then
            do i=1,q
               do j=1,q
                  workmr(i,j)=0.d0
                  sigmainv(i,j)=0.d0
               end do
            end do
            call cholesky(q,sigma,workmhr)
            do i=1,q
               do j=1,i
                  sigmainv(i,j)=workmhr(ihmssf(i,j,q))
               end do
            end do
            call inverse(sigmainv,q,iflagr)      
         
          else if(typep.eq.2)then
            do i=1,q
               do j=1,q
                  workmr(i,j)=0.d0
                  workmr1(i,j)=0.d0
                  sigmainv(i,j)=0.d0
               end do
            end do
            call eigenv(q,q,sigma,workvr,workvr1,workmr)
            do i=1,q
               workmr1(i,i)=sqrt(workvr(i))
            end do
         
            do i=1,q
               do j=1,q
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+workmr(i,k)*workmr1(k,j)
                  end do
                  sigmainv(i,j)=tmp1
               end do
            end do
            call inverse(sigmainv,q,iflagr)      
      
          else
            do i=1,q
               do j=1,q
                  workmr(i,j)=0.d0
                  workmr1(i,j)=0.d0
                  workmr2(i,j)=0.d0
                  sigmainv(i,j)=0.d0
               end do
            end do
            call eigenv(q,q,sigma,workvr,workvr1,workmr)
            do i=1,q
               workmr1(i,i)=sqrt(workvr(i))
            end do
 
            do i=1,q
               do j=1,q
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+workmr(i,k)*workmr1(k,j)
                  end do
                  workmr2(i,j)=tmp1
               end do
            end do

            do i=1,q
               do j=1,q
                  tmp1=0.d0
                  do k=1,q
                     tmp1=tmp1+workmr2(i,k)*workmr(j,k)
                  end do
                  sigmainv(i,j)=tmp1
               end do
            end do
            call inverse(sigmainv,q,iflagr)      
         end if 

c+++++++ transformation of the random effects

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,nsubject
            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+sigmainv(j,k)*(b(i,k)-mu(k))
               end do
               bz(i,j)=tmp1
            end do
         end do  

c+++++++ evaluate the grid
         
         do jj=1,ngrid
            theta(1)=grid(jj)
c+++++++++++check if the user has requested an interrupt
            call rchkusr()
            
            loglik=0.d0
            
            do i=1,q
               tmp1=0.d0
               do j=1,q
                  tmp1=tmp1+sigmainv(i,j)*(theta(j)-mu(j))
               end do
               thetaz(i)=tmp1
            end do

            nint=2
            prob=1.d0/dble(nint)
            quan=invcdfnorm(prob,0.d0,1.d0,1,0)

            countero=0
            
            do j=1,q
               if(thetaz(j).le.quan)then
                  linf(j)=-999999.d0
                  lsup(j)=quan
                  parti(j)=1
                else
                  linf(j)=quan
                  lsup(j)= 999999.d0
                  parti(j)=2
               end if
            end do

            do l=1,nsubject
               final=1
               do j=1,q
                  if(bz(l,j).gt.lsup(j).or.bz(l,j).lt.linf(j))then
                     final=0
                  end if
               end do
            
               if(final.eq.1)then
                  countero=countero+1
                  whicho(countero)=l
               end if   
            end do
            
            if(fixed.ne.1)then
               loglik=loglik+
     &           log((2.d0**q)*cpar+dble(2.d0**q)*dble(countero))-
     &           log((2.d0**q)*cpar+dble(nsubject))
            end if
            
            if(countero.eq.0) go to 1

            ok=1
            j=2
            do while(ok.eq.1.and.j.le.m)
               nint=2**j
               je2=j**2
               prob=1.d0/dble(nint)

               do k=1,q
                  k1=2*(parti(k)-1)+1
                  k2=2*(parti(k)-1)+2
                  quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
            
                  if(thetaz(k).le.quan)then
                    parti(k)=k1 
                    lsup(k)=quan
                   else 
                    parti(k)=k2
                    linf(k)=quan
                  end if
               end do                 
            
               countern=0
               do l=1,countero
                  final=1
                  do k=1,q
                     if(bz(whicho(l),k).gt.lsup(k).or.
     &                  bz(whicho(l),k).lt.linf(k)    )then
                        final=0 
                     end if   
                  end do
               
                  if(final.eq.1)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                  end if
               end do

               loglik=loglik+
     &           log((2.d0**q)*cpar*dble(je2)+
     &               dble(2.d0**q)*dble(countern))-
     &           log((2.d0**q)*cpar*dble(je2)+dble(countero))

               if(countern.eq.0)then
                  ok=0
                else  
                  countero=countern
                  do l=1,countern
                     whicho(l)=whichn(l)
                  end do
                  j=j+1
                end if   
            end do

1           continue

            loglik=loglik-0.5d0*detlogl
            do j=1,q
               loglik=loglik+dnrm(thetaz(j),0.d0, 1.d0, 1)
            end do   
             
            fs(jj)=fs(jj)+exp(loglik)            

         end do
      end do

      do i=1,ngrid
         fs(i)=fs(i)/dble(nsave)
      end do

      return
      end
      

c=======================================================================                  
      integer function binaryrep(nvar,pattern)
c=======================================================================                  
c     function that return the integer number+1, based on its
c     binary representation
c     Alejandro Jara, 2006

      implicit none
      integer i,nvar
      integer pattern(nvar)
      
      binaryrep=0
      do i=1,nvar
         binaryrep=binaryrep+(2**(i-1))*pattern(i)
      end do
      binaryrep=binaryrep+1
      return
      end


c=======================================================================                  
      subroutine binaryrepinv(nvar,evali,pattern)
c=======================================================================                  
c     function that return the binary representation of number 
c     given evali=number+1
c     Alejandro Jara, 2006

      implicit none
      integer i,nvar,evali,tmp1
      integer pattern(nvar)
      integer evali2
      real*8 tmp2
      
      evali2=evali-1
      do i=1,nvar
         tmp1=int(evali2/2.d0)
         tmp2=dble(evali2)/2.d0
         if(abs(tmp2-tmp1).gt.0.d0)then
           pattern(i)=1
          else
           pattern(i)=0
         end if  
         evali2=tmp1
      end do
      return
      end


c=======================================================================
      subroutine locationptu(x,m,n,loca)
c=======================================================================
c     function that return the succesive location of x across the 
c     finite tree of length m. This is based on a standard normal 
c     distribution.
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer m,n,loca(n)
      integer i,k,k1,k2,nint
      real*8 invcdfnorm
      real*8 prob,quan
      real*8 x
      
      if(m.gt.n)then
        call rexit("Error in ´locationpt´")
      end if

      quan=0.d0
      if(x.le.quan)then
        k=1
       else
        k=2
      end if  

      loca(1)=k
      
      do i=2,m
         nint=2**i
         prob=1.d0/dble(nint)
         k1=2*(k-1)+1
         k2=2*(k-1)+2
         quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
         if(x.le.quan)then
           k=k1
          else
           k=k2
         end if  
         loca(i)=k
      end do
      return
      end

c=======================================================================
      subroutine locationptm(nvar,x,m,maxm,parti,pattern,loca)
c=======================================================================
c     function that return the succesive location of the vector x(nvar) 
c     across the finite tree of length m. This is based on a standard 
c     normal distribution.
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer m,maxm,nvar,loca(maxm),parti(nvar),pattern(nvar)
      integer evali,i,j,k1,k2,nint,binaryrep
      real*8 invcdfnorm
      real*8 prob,quan
      real*8 x(nvar)
      
      if(m.gt.maxm)then
        call rexit("Error in ´locationptm´")
      end if

      quan=0.d0
      do i=1,nvar
         if(x(i).le.quan)then
            pattern(i)=0  
            parti(i)=1
           else
            pattern(i)=1
            parti(i)=2
         end if
      end do
      evali=binaryrep(nvar,pattern)

      loca(1)=evali
      
      do i=2,m
         nint=2**i
         prob=1.d0/dble(nint)
         do j=1,nvar
            k1=2*(parti(j)-1)+1
            k2=2*(parti(j)-1)+2
            quan=invcdfnorm(dble(k1)*prob,0.d0,1.d0,1,0)
            if(x(j).le.quan)then
               pattern(j)=0  
               parti(j)=k1
              else
               pattern(j)=1
               parti(j)=k2
            end if
         end do            
         evali=binaryrep(nvar,pattern)
         loca(i)=(2**nvar)*(loca(i-1)-1)+evali
      end do
      return
      end


c=======================================================================                  
      integer function binaryrep2(maxvar,nvar,pattern)
c=======================================================================                  
c     function that return the integer number+1, based on its
c     binary representation
c     Alejandro Jara, 2006

      implicit none
      integer i,maxvar,nvar
      integer pattern(maxvar)
      
      binaryrep2=0
      do i=1,nvar
         binaryrep2=binaryrep2+(2**(i-1))*pattern(i)
      end do
      binaryrep2=binaryrep2+1
      return
      end


c=======================================================================
      subroutine accumcountpt(nsubject,nvar,m,z,maxm,loca,parti,pattern,
     &                        theta,ntotals,nhash,maxnzr,tmp,counts,nr,
     &                        ia,ja,a)
c=======================================================================
c     accumulate the number of subjects.  
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer i,j
      integer m,nsubject,nvar,maxm
      integer loca(maxm),parti(nvar),pattern(nvar)
      real*8 z(nsubject,nvar),theta(nvar),one
      parameter(one=1.d0)
      
      integer nhash,maxnzr,nr
      real*8 counts(nhash,3)
      integer ntotals,ia(ntotals+1),ja(maxnzr),tmp(ntotals)
      real*8 a(maxnzr)

      nr=0
      do i=1,nsubject
         do j=1,nvar
           theta(j)=z(i,j)
         end do
         call locationptm(nvar,theta,m,maxm,parti,pattern,loca)
         
         do j=1,m
            call hashm(one,loca(j),j,counts,nhash,nr)
         end do
      end do
      
      call hashiajaa(counts,nhash,ntotals,ia,ja,a,maxnzr,tmp)
      
      return
      end


c=======================================================================
      subroutine sprobpt(fixed,cpar,nvar,m,ntotals,nhash,maxnzr,
     &                   ia,ja,a,maxarea,mass,rvecs,tmp1,
     &                   tmp2)
c=======================================================================
c     generate probabilities for a finite PT.  
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer fixed,i,j,k,l,ista,iend
      integer m,maxarea,narea,nvar,ntotals,nhash,maxnzr
      integer ia(ntotals+1),ja(maxnzr)
      real*8 a(maxnzr)
      real*8 cpar
      real*8 mass(maxarea)
      real*8 rvecs(maxarea)
      real*8 tmp1(ntotals),tmp2(ntotals)

      integer ngroup
      real*8 je2
      
      narea=2**nvar 
      if(narea.gt.maxarea)then
        call rexit("Error in ´sprobpt´")      
      end if
      
c++++ First level probabilities
      ngroup=1
      do k=1,narea
         mass(k)=cpar
      end do

      ista=1
      iend=(2**nvar)
      
      if(fixed.ne.1)then
         do k=ista,iend
            do l=ia(k),ia(k+1)-1
               if(ja(l).eq.1)mass(k-ista+1)=mass(k-ista+1)+a(l)
            end do
         end do
         call dirichlet(mass,maxarea,narea,rvecs)
         do k=ista,iend
            tmp1(k)=rvecs(k)
            tmp2(k)=rvecs(k)
         end do
       else
         do k=ista,iend
            rvecs(k)=1.d0/dble(narea)
            tmp1(k)=rvecs(k)
            tmp2(k)=rvecs(k)
         end do
      end if
      
c      call dblepr("rvecs",-1,rvecs,narea)

c++++ Probabilities for the other levels
      do i=2,m
         je2=i**2
         ngroup=2**(nvar*(i-1))
         do j=1,ngroup

            do k=1,narea
               mass(k)=cpar*je2
            end do
            
            ista=narea*(j-1)+1
            iend=ista+narea-1
            
            do k=ista,iend
               do l=ia(k),ia(k+1)-1
                  if(ja(l).eq.i)mass(k-ista+1)=mass(k-ista+1)+a(l)
               end do
            end do
            call dirichlet(mass,maxarea,narea,rvecs)
            
c            call dblepr("rvecs",-1,rvecs,narea)
            
            do k=ista,iend
               tmp2(k)=rvecs(k-ista+1)*tmp1(j)
            end do
         end do
         
         if(i.lt.m)then
            do j=1,2**(nvar*i)
               tmp1(j)=tmp2(j)
            end do
         end if   
      end do
      
      return
      end

c=======================================================================      
      subroutine quandetpt(linf,lsup,x)
c=======================================================================
c     determine the quantile corresponding to the next partition of the
c     interval (linf,lsup)
c     Alejandro Jara, 2007 
c=======================================================================
      implicit none
      real*8 linf,lsup
      real*8 invcdfnorm,cdfnorm,x
      real*8 tmp1,tmp2,tmp3
      
      tmp1=cdfnorm(linf,0.d0,1.d0,1,0)
      tmp2=cdfnorm(lsup,0.d0,1.d0,1,0)
      tmp3=tmp1+(tmp2-tmp1)/2.d0

      x=invcdfnorm(tmp3,0.d0,1.d0,1,0)
      return
      end
      

c=======================================================================
      subroutine setspt(nvar,m,ntotals,pattern,maxvar,
     &                  liminf1,liminf2,limsup1,limsup2)
c=======================================================================
c     computes the limtis of the sets in a finite multivariate PT
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer maxvar,nvar,m,ntotals
      integer pattern(nvar)
      real*8 liminf1(ntotals,maxvar),liminf2(ntotals,maxvar)
      real*8 limsup1(ntotals,maxvar),limsup2(ntotals,maxvar)
      integer i,j,k,l,ngroup,ista,narea
      real*8 linf,lsup,quan

c++++ Algorithm

      narea=2**nvar  
      
c++++ First level sets

      do i=1,narea
         call binaryrepinv(nvar,i,pattern)
         do j=1,nvar
            if(pattern(j).eq.0)then
               liminf1(i,j)=-999999.d0
               limsup1(i,j)=0.d0
               liminf2(i,j)=-999999.d0
               limsup2(i,j)=0.d0
             else
               liminf1(i,j)=0.d0
               limsup1(i,j)=+999999.d0
               liminf2(i,j)=0.d0
               limsup2(i,j)=+999999.d0
            end if 
         end do
      end do

c++++ Sets for other levels

      do i=2,m

         ngroup=2**(nvar*(i-1))
         do j=1,ngroup

            ista=(2**nvar)*(j-1)+1

            do k=1,narea
               call binaryrepinv(nvar,k,pattern)

               do l=1,nvar
                  linf=liminf1(j,l)
                  lsup=limsup1(j,l)
                  call quandetpt(linf,lsup,quan)
                  
                  if(pattern(l).eq.0)then
                      liminf2(ista+k-1,l)=linf
                      limsup2(ista+k-1,l)=quan
                   else
                      liminf2(ista+k-1,l)=quan
                      limsup2(ista+k-1,l)=lsup
                  end if 
               end do
            end do
         end do
         
         if(i.lt.m)then
            do j=1,2**(nvar*i)
               do k=1,nvar   
                  liminf1(j,k)=liminf2(j,k)
                  limsup1(j,k)=limsup2(j,k)
               end do  
            end do
         end if   
      end do
      return
      end


c=======================================================================
      subroutine mompt(maxvar,nvar,m,ntotals,liminf,limsup,probs,
     &                 means,covs)
c=======================================================================
c     computes the moments of a finite multivariate PT
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer maxvar,nvar,m,ntotals
      real*8 liminf(ntotals,maxvar)
      real*8 limsup(ntotals,maxvar)
      real*8 probs(ntotals)
      real*8 means(nvar),covs(nvar,nvar)
      
      integer i,j,k,nsets
      real*8 linf,lsup,muwork1,muwork2
      real*8 cdfnorm,dnrm
      real*8 tmp1,tmp2,tmp3,tmp4
      real*8 tmass

c++++ Algorithm

      nsets=2**(nvar*m)  
      do i=1,nvar
         means(i)=0.d0 
         do j=1,nvar
            covs(i,j)=0.d0 
         end do
      end do
      
      if(nsets.gt.ntotals)then
         call rexit("Error in ´mompt´: nsets > ntotals")            
      end if
      
      tmass=0.d0
      do i=1,nsets
         do j=1,nvar
            linf=liminf(i,j)
            lsup=limsup(i,j)
            tmp1=dnrm(linf,0.d0,1.d0,0)
            tmp2=dnrm(lsup,0.d0,1.d0,0)
            tmp3=cdfnorm(linf,0.d0,1.d0,1,0)
            tmp4=cdfnorm(lsup,0.d0,1.d0,1,0)
            muwork1= - (tmp2-tmp1)/(tmp4-tmp3)
            means(j)=means(j)+probs(i)*muwork1
            
            do k=j+1,nvar
               linf=liminf(i,k)
               lsup=limsup(i,k)
               tmp1=dnrm(linf,0.d0,1.d0,0)
               tmp2=dnrm(lsup,0.d0,1.d0,0)
               tmp3=cdfnorm(linf,0.d0,1.d0,1,0)
               tmp4=cdfnorm(lsup,0.d0,1.d0,1,0)
               muwork2= - (tmp2-tmp1)/(tmp4-tmp3)
               covs(j,k)=covs(j,k)+ probs(i)*muwork1*muwork2
               covs(k,j)=covs(k,j)+ probs(i)*muwork1*muwork2
            end do

            muwork2=1.d0 - (lsup*tmp2-linf*tmp1)/(tmp4-tmp3) 
            covs(j,j)=covs(j,j)+probs(i)*muwork2
            
         end do
         tmass=tmass+probs(i)
      end do

      do i=1,nvar
         do j=1,nvar
            covs(i,j)=covs(i,j)-means(i)*means(j)   
         end do   
      end do  

      return
      end

c=======================================================================
      subroutine samplefuncpt(fixed,m,nrand,nsubject,cpar,bz,theta,
     &                        parti,pattern,means,covs)
c=======================================================================
      implicit none
      integer i,j
      integer fixed,m,nrand,nsubject
      integer parti(nrand),pattern(nrand)
      real*8 cpar,bz(nsubject,nrand),theta(nrand)
      real*8 means(nrand),covs(nrand,nrand)

c++++ parameters
      integer maxrand,maxm,maxarea,ntotals,ntotalp,nhash,maxnzr
      parameter(maxrand=10,maxm=20)
      parameter(maxarea=2**(maxrand),ntotals=2**(2*6),ntotalp=5460)
      parameter(nhash=ntotalp,maxnzr=ntotalp)
      integer loca(maxm),mwork,nr
      integer ia(ntotals+1),ja(maxnzr),tmp(ntotals)

      real*8 a(maxnzr)    
      real*8 counts(nhash,3)
      real*8 mass(maxarea)
      real*8 rvecs(maxarea)
      real*8 probw1(ntotals),probw2(ntotals)
      real*8 liminf1(ntotals,maxrand),liminf2(ntotals,maxrand)
      real*8 limsup1(ntotals,maxrand),limsup2(ntotals,maxrand)
      
      if(nrand.gt.maxrand)then
        call rexit("Error in ´samplefuncpt´: increase maxrand")      
      end if
      
      mwork=m
      do while(ntotals.lt.(2**(nrand*mwork)))
         mwork=mwork-1
      end do

      if(mwork.lt.1)then
        call rexit("Error in ´samplefuncpt´: maxm < 1")      
      end if
      
      if(mwork.gt.maxm)then
        call rexit("Error in ´samplefuncpt´: increase maxm")      
      end if

      nr=0
      do i=1,nhash
         do j=1,3
            counts(i,j)=0.d0
         end do   
      end do

c++++ set the MPT parameters      
      
      call accumcountpt(nsubject,nrand,mwork,bz,maxm,loca,parti,pattern,
     &                  theta,ntotals,nhash,maxnzr,tmp,counts,nr,
     &                  ia,ja,a)

c++++ sampling probabilities

      call sprobpt(fixed,cpar,nrand,mwork,ntotals,nhash,maxnzr,
     &             ia,ja,a,maxarea,mass,rvecs,probw1,
     &             probw2)

c++++ finding sets

      call setspt(nrand,mwork,ntotals,pattern,maxrand,
     &            liminf1,liminf2,limsup1,limsup2)

c++++ computing the means for each set

      call mompt(maxrand,nrand,mwork,ntotals,liminf2,limsup2,
     &           probw2,means,covs)

      return
      end




