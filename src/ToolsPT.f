c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR MULTIVARIATES POLYA TREES
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
      subroutine loglikpt_mucan(m,nrand,nsubject,parti,
     &                          whicho,whichn,b,bzc,cpar,detlogl,
     &                          linf,lsup,muc,sigmainv,
     &                          vec,fixed,loglikc)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     value of the baseline mean in a marginal Multivariate PT.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      double precision b(nsubject,nrand),bzc(nsubject,nrand),
     1  cpar,detlogl
      double precision linf(nrand),lsup(nrand)
      double precision muc(nrand),sigmainv(nrand,nrand)
      double precision vec(nrand)

c-----Output
      double precision loglikc

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan,tmp1

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
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,ind,nrand,nsubject,parti(nrand),m
      integer whicho(nsubject),whichn(nsubject)
      double precision bz(nsubject,nrand),cpar,detlogl
      double precision linf(nrand),lsup(nrand)
      double precision vec(nrand)

c-----Output
      double precision logprior

c-----Working
      integer countero,countern,final,i,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan

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
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,ind,m,nrand,nsubject
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      double precision bz(nsubject,nrand),cpar,detlogl
      double precision linf(nrand),lsup(nrand)

c-----Input/Output
      double precision logliko

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan

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
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      double precision bz(nsubject,nrand),cpar,detlogl
      double precision linf(nrand),lsup(nrand)

c-----Output
      double precision logliko

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan

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
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject,typep
      integer iflagr(nrand)
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      double precision b(nsubject,nrand),bzc(nsubject,nrand),
     1  cpar,detloglc
      double precision linf(nrand),lsup(nrand)
      double precision mu(nrand),sigmac(nrand,nrand),
     1  sigmainvc(nrand,nrand)
      double precision vec(nrand)
      double precision workmhr(nrand*(nrand+1)/2),workmr(nrand,nrand)
      double precision workmr1(nrand,nrand),workmr2(nrand,nrand)
      double precision workvr(nrand)

c-----Output
      double precision loglikc

c-----Working
      integer countero,countern,final
      integer i,ihmssf,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan,tmp1

c-----Routine

      loglikc=0.d0

      do i=1,nrand
         do j=1,nrand
            workmr(i,j)=sigmac(i,j)
         end do
      end do
      call inversedet(workmr,nrand,iflagr,detloglc)

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
     &                            vec,fixed,loglikn)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     value of the precision parameter in a marginal Multivariate PT.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer iflagr(nrand)
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      double precision bz(nsubject,nrand),cparc,detlogl
      double precision linf(nrand),lsup(nrand)
      double precision vec(nrand)

c-----Output
      double precision loglikn

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan

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
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,marea,narea,nrand,nsubject,parti(nrand)
      integer m,massi(marea),pattern(nrand),patterns(nrand)
      integer whicho(nsubject),whichn(nsubject)      
      double precision bz(nsubject,nrand),cpar
      double precision limw(nrand),linf(nrand),lsup(nrand) 
      double precision mu(nrand),sigma(nrand,nrand)
      double precision mass(marea),rtnorm

      double precision workmhr(nrand*(nrand+1)/2),workmr(nrand,nrand)
      double precision workvr(nrand)
      
c-----Output
      double precision vec(nrand)

c-----Working
      integer binaryrep,countero,countern,evali,evali2,final
      integer i,ihmssf,j,je2,k,k1,k2,l,nint,ok 
      double precision invcdfnorm,prob,quan,tmp1

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
     &                       cpar,limw,linf,lsup,vec,fixed)
c======================================================================= 
c     This subroutine generates a sample 'vec' from the predictive
c     distribution arising in a marginal Multivariate PT.
c     The values of random effects 'bz' must be in a standarized form.
c     The output 'vec' is also standarized.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,marea,narea,nrand,nsubject,parti(nrand)
      integer m,massi(marea),pattern(nrand),patterns(nrand)
      integer whicho(nsubject),whichn(nsubject)      
      double precision bz(nsubject,nrand),cpar
      double precision limw(nrand),linf(nrand),lsup(nrand) 
      double precision mass(marea),rtnorm

c-----Output
      double precision vec(nrand)

c-----Working
      integer binaryrep,countero,countern,evali,evali2,final
      integer i,j,je2,k,k1,k2,l,nint,ok 
      double precision invcdfnorm,prob,quan

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
     &                        vec,fixed,logliko)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the current value
c     of the baseline parameters in a marginal Multivariate PT.
c     This function does not include the standarization of random 
c     effects
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      double precision bz(nsubject,nrand),cpar,detlogl
      double precision linf(nrand),lsup(nrand)
      double precision vec(nrand)

c-----Output
      double precision logliko

c-----Working
      integer countero,countern,final
      integer i,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan

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
     &                         sigmamat,cparvec,typepvec,
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
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 24-04-2007.
c=======================================================================
      implicit none

c++++ input
      integer fixed,m,ngrid1,ngrid2,nsubject,nsave,q,typepvec(nsave)
      double precision cparvec(nsave),randsave(nsave,q*(nsubject+1))
      double precision mumat(nsave,q)
      double precision sigmamat(nsave,q*(q+1)/2)
      double precision grid1(ngrid1),grid2(ngrid2)

c++++ output
      double precision fs(ngrid1,ngrid2)

c++++ external working space
      integer iflagr(q),parti(q)
      integer whicho(nsubject),whichn(nsubject)
      double precision b(nsubject,q),bz(nsubject,q)
      double precision linf(q),lsup(q),mu(q),sigma(q,q),sigmainv(q,q)
      double precision theta(q),thetaz(q)
      double precision workmr(q,q),workmr1(q,q),workmr2(q,q)
      double precision workmhr(q*(q+1)/2)
      double precision workvr(q),workvr1(q)

c++++ internal working space
      integer countero,countern
      integer final,i,ii,ihmssf,j,jj,je2,k,kk,k1,k2,l
      integer nint,ok
      double precision cpar,detlogl,dnrm
      double precision invcdfnorm
      double precision loglik,prob,quan,tmp1

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

         if(typepvec(ii).eq.1)then
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
         
          else if(typepvec(ii).eq.2)then
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
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 24-04-2007.
c=======================================================================
      implicit none

c++++ input
      integer fixed,m,ngrid,nsubject,nsave,q,typep
      double precision cparvec(nsave),randsave(nsave,q*(nsubject+1))
      double precision mumat(nsave,q)
      double precision sigmamat(nsave,q*(q+1)/2)
      double precision grid(ngrid)

c++++ output
      double precision fs(ngrid)

c++++ external working space
      integer iflagr(q),parti(q)
      integer whicho(nsubject),whichn(nsubject)
      double precision b(nsubject,q),bz(nsubject,q)
      double precision linf(q),lsup(q),mu(q),sigma(q,q),sigmainv(q,q)
      double precision theta(q),thetaz(q)
      double precision workmr(q,q),workmr1(q,q),workmr2(q,q)
      double precision workmhr(q*(q+1)/2)
      double precision workvr(q),workvr1(q)

c++++ internal working space
      integer countero,countern
      integer final,i,ii,ihmssf,j,jj,je2,k,k1,k2,l
      integer nint,ok
      double precision cpar,detlogl,dnrm
      double precision invcdfnorm
      double precision loglik,prob,quan,tmp1

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
c     Alejandro Jara, 2006-2007-2008

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
c     Alejandro Jara, 2006-2007-2008

      implicit none
      integer i,nvar,evali,tmp1
      integer pattern(nvar)
      integer evali2
      double precision tmp2
      
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
c     Alejandro Jara, 2007-2008
c=======================================================================
      implicit none
      integer m,n,loca(n)
      integer i,k,k1,k2,nint
      double precision invcdfnorm
      double precision prob,quan
      double precision x
      
      if(m.gt.n)then
        call rexit("Error in 'locationpt'")
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
c     Alejandro Jara, 2007-2008
c=======================================================================
      implicit none
      integer m,maxm,nvar,loca(maxm),parti(nvar),pattern(nvar)
      integer evali,i,j,k1,k2,nint,binaryrep
      double precision invcdfnorm
      double precision prob,quan
      double precision x(nvar)
      
      if(m.gt.maxm)then
        call rexit("Error in 'locationptm'")
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
c     Alejandro Jara, 2006-2007-2008

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
c     Alejandro Jara, 2007-2008
c=======================================================================
      implicit none
      integer i,j
      integer m,nsubject,nvar,maxm
      integer loca(maxm),parti(nvar),pattern(nvar)
      double precision z(nsubject,nvar),theta(nvar),one
      parameter(one=1.d0)
      
      integer nhash,maxnzr,nr
      double precision counts(nhash,3)
      integer ntotals,ia(ntotals+1),ja(maxnzr),tmp(ntotals)
      double precision a(maxnzr)

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
c     Alejandro Jara, 2007-2008
c=======================================================================
      implicit none
      integer fixed,i,j,k,l,ista,iend
      integer m,maxarea,narea,nvar,ntotals,nhash,maxnzr
      integer ia(ntotals+1),ja(maxnzr)
      double precision a(maxnzr)
      double precision cpar
      double precision mass(maxarea)
      double precision rvecs(maxarea)
      double precision tmp1(ntotals),tmp2(ntotals)

      integer ngroup
      double precision je2
      
      narea=2**nvar 
      if(narea.gt.maxarea)then
        call rexit("Error in 'sprobpt'")      
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
c     Alejandro Jara, 2007-2008 
c=======================================================================
      implicit none
      double precision linf,lsup
      double precision invcdfnorm,cdfnorm,x
      double precision tmp1,tmp2,tmp3
      
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
c     Alejandro Jara, 2007-2008
c=======================================================================
      implicit none
      integer maxvar,nvar,m,ntotals
      integer pattern(nvar)
      double precision liminf1(ntotals,maxvar),liminf2(ntotals,maxvar)
      double precision limsup1(ntotals,maxvar),limsup2(ntotals,maxvar)
      integer i,j,k,l,ngroup,ista,narea
      double precision linf,lsup,quan

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
c     Alejandro Jara, 2007-2008
c=======================================================================
      implicit none
      integer maxvar,nvar,m,ntotals
      double precision liminf(ntotals,maxvar)
      double precision limsup(ntotals,maxvar)
      double precision probs(ntotals)
      double precision means(nvar),covs(nvar,nvar)
      
      integer i,j,k,nsets
      double precision linf,lsup,muwork1,muwork2
      double precision cdfnorm,dnrm
      double precision tmp1,tmp2,tmp3,tmp4
      double precision tmass

c++++ Algorithm

      nsets=2**(nvar*m)  
      do i=1,nvar
         means(i)=0.d0 
         do j=1,nvar
            covs(i,j)=0.d0 
         end do
      end do
      
      if(nsets.gt.ntotals)then
         call rexit("Error in 'mompt': nsets > ntotals")            
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
      double precision cpar,bz(nsubject,nrand),theta(nrand)
      double precision means(nrand),covs(nrand,nrand)

c++++ parameters
      integer maxrand,maxm,maxarea,ntotals,ntotalp,nhash,maxnzr
      parameter(maxrand=10,maxm=20)
      parameter(maxarea=2**(maxrand),ntotals=2**(2*6),ntotalp=5460)
      parameter(nhash=ntotalp,maxnzr=ntotalp)
      integer loca(maxm),mwork,nr
      integer ia(ntotals+1),ja(maxnzr),tmp(ntotals)

      double precision a(maxnzr)    
      double precision counts(nhash,3)
      double precision mass(maxarea)
      double precision rvecs(maxarea)
      double precision probw1(ntotals),probw2(ntotals)
      double precision liminf1(ntotals,maxrand),liminf2(ntotals,maxrand)
      double precision limsup1(ntotals,maxrand),limsup2(ntotals,maxrand)
      
      if(nrand.gt.maxrand)then
        call rexit("Error in 'samplefuncpt': increase maxrand")      
      end if
      
      mwork=m
      do while(ntotals.lt.(2**(nrand*mwork)))
         mwork=mwork-1
      end do

      if(mwork.lt.1)then
        call rexit("Error in 'samplefuncpt': maxm < 1")      
      end if
      
      if(mwork.gt.maxm)then
        call rexit("Error in 'samplefuncpt': increase maxm")      
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


c=======================================================================                  
      subroutine loglikpt_covarcan2(m,nrand,nsubject,iflagr,parti,
     &                              whicho,whichn,b,bzc,cpar,detloglc,
     &                              linf,lsup,mu,sigmac,sigmainvc,
     &                              ortho,vec,workmhr,workmr,
     &                              loglikc,fixed)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for the candidate 
c     of the baseline covariance matrix in a marginal Multivariate PT.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 22-06-2008.
c======================================================================= 
      implicit none 

c-----Input
      integer fixed,m,nrand,nsubject
      integer iflagr(nrand)
      integer parti(nrand)
      integer whicho(nsubject),whichn(nsubject)
      double precision b(nsubject,nrand),bzc(nsubject,nrand),
     1  cpar,detloglc
      double precision linf(nrand),lsup(nrand)
      double precision mu(nrand),sigmac(nrand,nrand),
     1  sigmainvc(nrand,nrand)
      double precision vec(nrand)
      double precision workmhr(nrand*(nrand+1)/2),workmr(nrand,nrand)
      double precision ortho(nrand,nrand)

c-----Output
      double precision loglikc

c-----Working
      integer countero,countern,final
      integer i,ihmssf,j,je2,k,k1,k2,l,nint,ok
      double precision dnrm,invcdfnorm,prob,quan,tmp1

c-----Routine

      loglikc=0.d0

      do i=1,nrand
         do j=1,nrand
            workmr(i,j)=sigmac(i,j)
         end do
      end do
      call inversedet(workmr,nrand,iflagr,detloglc)

      do i=1,nrand
         do j=1,nrand
            workmr(i,j)=0.d0
            sigmainvc(i,j)=0.d0
         end do
      end do
      call cholesky(nrand,sigmac,workmhr)
      do i=1,nrand
         do j=1,i
            workmr(i,j)=workmhr(ihmssf(i,j,nrand))
         end do
      end do

      do i=1,nrand
         do j=1,nrand
            tmp1=0.d0
            do k=1,nrand 
               tmp1=tmp1+workmr(i,k)*ortho(k,j)
            end do
            sigmainvc(i,j)=tmp1
         end do
      end do

      call inverse(sigmainvc,nrand,iflagr)      
         
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
      subroutine predictiveptb2(m,nsubject,q,nsave,randsave,mumat,
     &                          sigmamat,cparvec,typepmat,
     &                          ngrid1,ngrid2,grid1,grid2,fs,
     &                          iflagr,parti,whicho,whichn,
     &                          b,bz,linf,lsup,mu,sigma,sigmainv,
     &                          theta,thetaz,workmr,ortho,
     &                          workmhr,fixed)
c=======================================================================
c     computes the bivariate posterior predictive density from the
c     output of a PTfunction. This is used for random effects models.
c
c     Note that if fixed=1, the first level probabilities are fixed at 
c     (1/2)**nrand
c
c     Alejandro Jara, 2006-2007-2008
c     Last modification: 24-04-2007.
c=======================================================================
      implicit none

c++++ input
      integer fixed,m,ngrid1,ngrid2,nsubject,nsave,q
      double precision cparvec(nsave),randsave(nsave,q*(nsubject+1))
      double precision mumat(nsave,q)
      double precision sigmamat(nsave,q*(q+1)/2)
      double precision typepmat(nsave,q*q)
      double precision grid1(ngrid1),grid2(ngrid2)

c++++ output
      double precision fs(ngrid1,ngrid2)

c++++ external working space
      integer iflagr(q),parti(q)
      integer whicho(nsubject),whichn(nsubject)
      double precision b(nsubject,q),bz(nsubject,q)
      double precision linf(q),lsup(q),mu(q),sigma(q,q),sigmainv(q,q)
      double precision theta(q),thetaz(q)
      double precision workmr(q,q),ortho(q,q)
      double precision workmhr(q*(q+1)/2)

c++++ internal working space
      integer count 
      integer countero,countern
      integer final,i,ii,ihmssf,j,jj,je2,k,kk,k1,k2,l
      integer nint,ok
      double precision cpar,detlogl,dnrm
      double precision invcdfnorm
      double precision loglik,prob,quan,tmp1

c++++ algorithm      
      if(q.gt.2)then
        call rexit("Only bivariate evaluation supported")
      end if

      do ii=1,nsave

c+++++++ save elements

c+++++++ check if the user has requested an interrupt
         call rchkusr()
          
         cpar=cparvec(ii)
         
         count=0
         do i=1,q
            mu(i)=mumat(ii,i)
            do j=1,q
               count=count+1
               sigma(i,j)=sigmamat(ii,ihmssf(i,j,q))
               ortho(i,j)=typepmat(ii,count) 
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

         do i=1,q
           do j=1,q
              workmr(i,j)=0.d0
              sigmainv(i,j)=0.d0
           end do
         end do
         call cholesky(q,sigma,workmhr)
         do i=1,q
            do j=1,i
               workmr(i,j)=workmhr(ihmssf(i,j,q))
            end do
         end do

         do i=1,q
            do j=1,q
               tmp1=0.d0
               do k=1,q
                  tmp1=tmp1+workmr(i,k)*ortho(k,j) 
               end do 
               sigmainv(i,j)=tmp1
            end do
         end do

         call inverse(sigmainv,q,iflagr)      
         

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
      subroutine sampredellpt(nrec,nsubject,nrand,m,ssb,cparb,whicho,
     &                        whichn,workm,workmh,workv,mub,sigmab,b)   
c=======================================================================
c     This subroutine generates a sample 'b' from the predictive
c     distribution arising in a marginal Multivariate PT using a 
c     elliptical binary partition.
c     The output 'b' is in a normal form.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c+++++Input 
      integer nrec,nsubject,nrand,m
      integer whicho(nrec),whichn(nrec)
      double precision cparb,ssb(nsubject),b(nrand)
      double precision workm(nrand,nrand),workmh(nrand*(nrand+1)/2),
     1  workv(nrand)
      double precision mub(nrand),sigmab(nrand,nrand)

c+++++internal working variables
      integer countero,countern 
      integer i,ihmssf,j,je2,k1,k2  
      integer n1,n2
      integer nint
      integer ok
      integer parti
      
      double precision ellip
      double precision invcdfchisq
      double precision linf,lsup
      double precision prob
      double precision quan
      double precision rnorm
      double precision rtchisq
      double precision tmp1,tmp2
      
      real runif
      
      logical ainf,binf

c+++++algorithm

      j=1  
      nint=2
      prob=1.d0/dble(nint)
      quan=invcdfchisq(prob,dble(nrand),1,0)
      
      n1=0
      n2=0

      do i=1,nsubject
         if(ssb(i).le.quan)then
            n1=n1+1
           else
            n2=n2+1
         end if  
      end do

      tmp1=log(     cparb+dble(n1))-
     &     log(2.d0*cparb+dble(nsubject))
      
      tmp1=exp(tmp1)
      countero=0

      if(dble(runif()).le.tmp1)then
         parti=1
         do i=1,nsubject
            if(ssb(i).le.quan)then
               countero=countero+1
               whicho(countero)=i
            end if  
         end do
        else
         parti=2
         do i=1,nsubject
            if(ssb(i).gt.quan)then
               countero=countero+1
               whicho(countero)=i
             end if  
         end do
      end if

      if(countero.eq.0) go to 1  

      ok=1
      j=2
      do while(ok.eq.1.and.j.le.m)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
   
         k1=2*(parti-1)+1
         k2=2*(parti-1)+2

         quan=invcdfchisq(dble(k1)*prob,dble(nrand),1,0)

         countern=0

         do i=1,countero
            if(ssb(whicho(i)).le.quan)then
               countern=countern+1
            end if   
         end do

         tmp1=log(     cparb*dble(je2)+dble(countern))-
     &        log(2.d0*cparb*dble(je2)+dble(countero))
      
         tmp1=exp(tmp1)


         if(dble(runif()).le.tmp1)then
            parti=k1
            do i=1,countero
               if(ssb(whicho(i)).le.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(i)
               end if   
            end do
           else
            parti=k2
            do i=1,countero
               if(ssb(whicho(i)).le.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(i)
               end if   
            end do
         end if

         if(countern.eq.0)then
            ok=0
           else 
            countero=countern
            do i=1,countern
               whicho(i)=whichn(i)
            end do
            j=j+1
         end if   
      end do

1     continue

c++++ Now j indicates the level of the partition and parti the interval

      nint=2**j
      prob=1.d0/dble(nint)

      if(parti.eq.1)then
         quan=invcdfchisq(dble(parti)*prob,dble(nrand),1,0)

         ainf=.true.
         binf=.false.
         linf=0.d0
         lsup=quan
         ellip=rtchisq(dble(nrand),linf,lsup,ainf,binf)
         
       else if(parti.eq.nint)then
         quan=invcdfchisq(dble(parti-1)*prob,dble(nrand),1,0)

         ainf=.false.
         binf=.true.
         linf=quan
         lsup=0.d0
         ellip=rtchisq(dble(nrand),linf,lsup,ainf,binf)

       else
         tmp1=invcdfchisq(dble(parti-1)*prob,dble(nrand),1,0)
         tmp2=invcdfchisq(dble(parti  )*prob,dble(nrand),1,0)
         
         if(tmp1.ge.tmp2)then
            call rexit("Error in the limits")
         end if  

         ainf=.false.
         binf=.false.
         linf=tmp1
         lsup=tmp2
         ellip=rtchisq(dble(nrand),linf,lsup,ainf,binf)
      end if
      
      tmp1=0.d0
      do i=1,nrand
         workv(i)=rnorm(0.d0,1.d0)
         tmp1=tmp1+workv(i)**2
      end do
      
      do i=1,nrand
         workv(i)=sqrt(ellip)*workv(i)/sqrt(tmp1)
      end do


      call cholesky(nrand,sigmab,workmh)

      do i=1,nrand
         do j=1,nrand
            workm(i,j)=0.d0
         end do
      end do
        
      do i=1,nrand
         do j=1,i
            workm(i,j)=workmh(ihmssf(i,j,nrand))
         end do
      end do

      do i=1,nrand
         tmp1=0.d0
         do j=1,nrand
            tmp1=tmp1+workm(i,j)*workv(j)   
         end do
         b(i)=tmp1+mub(i)
      end do

      return
      end



c=======================================================================                  
      subroutine loglik_ell(m,nrec,nsubject,nrand,b,ssb,
     &                      sigmainvbc,detlogsbc,cparb,
     &                      whicho,whichn,
     &                      loglikc)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for a random effect
c     in a marginal Multivariate PT, using an elliptical binary 
c     partition.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none 

c-----Input
      integer m,nrand,nrec,nsubject
      integer whicho(nrec),whichn(nrec)
      double precision cparb,detlogsbc,b(nsubject,nrand)
      double precision ssb(nsubject)
      double precision sigmainvbc(nrand,nrand)
       
c-----Output
      double precision loglikc,sseval

c-----Working
      integer i,j,je2,k,k1,k2,l
      integer countero,countern,parti,nint,ok
      double precision invcdfchisq,prob,quan,tmp1,tmp2,tmp3,tpi

c-----Routine

      loglikc=0.d0
      tpi=6.283185307179586476925286766559d0      
      
      do i=1,nsubject

         
c+++++++ check if the user has requested an interrupt
         call rchkusr()

         tmp1=0.d0
         do k=1,nrand
            do l=1,nrand
               tmp1=tmp1+b(i,k)*sigmainvbc(k,l)*b(i,l)
            end do
         end do
         ssb(i)=tmp1
         sseval=tmp1

c+++++++ first subject
         if(i.eq.1)then

           tmp1=-(dble(nrand)*log(tpi))
           tmp2=detlogsbc
           tmp3=sseval
      
           loglikc=loglikc+(tmp1-tmp2-tmp3)/2.d0

c+++++++ following subjects
          else
          
           nint=2
           prob=1.d0/dble(nint)
           quan=invcdfchisq(prob,dble(nrand),1,0)
      
           countero=0
      
           if(sseval.le.quan)then
              parti=1
              do l=1,i-1
                if(ssb(l).le.quan)then
                   countero=countero+1
                   whicho(countero)=l
                 end if  
              end do
            else
              parti=2
              do l=1,i-1
                 if(ssb(l).gt.quan)then
                   countero=countero+1
                   whicho(countero)=l
                 end if  
              end do
           end if  
   
           loglikc=loglikc+
     &         log(2.d0*cparb+dble(2*countero))-
     &         log(2.d0*cparb+dble(i-1))

           if(countero.eq.0) go to 1  

           ok=1
           j=2
           do while(ok.eq.1.and.j.le.m)
              nint=2**j
              je2=j**2
              prob=1.d0/dble(nint)
   
              k1=2*(parti-1)+1
              k2=2*(parti-1)+2
              quan=invcdfchisq(dble(k1)*prob,dble(nrand),1,0)

              if(sseval.le.quan)then
                 parti=k1
                 k=k1
                else
                 parti=k2
                 k=k2
              end if  

              countern=0

              if(k.eq.1)then
                 do l=1,countero
                    if(ssb(whicho(l)).le.quan.and.
     &                 whicho(l).ne.i)then
                       countern=countern+1
                       whichn(countern)=whicho(l)
                    end if   
                 end do
               else if(k.eq.nint)then
                 quan=invcdfchisq(dble(k-1)*prob,dble(nrand),1,0)
                 do l=1,countero
                    if(ssb(whicho(l)).gt.quan.and.
     &                 whicho(l).ne.i)then
                       countern=countern+1
                       whichn(countern)=whicho(l)
                    end if   
                 end do
               else
                 tmp1=invcdfchisq(dble(k-1)*prob,dble(nrand),1,0)
                 tmp2=invcdfchisq(dble(k  )*prob,dble(nrand),1,0)

                 if(tmp1.ge.tmp2)then
                   call rexit("Error in the limits")
                 end if  
         
                 do l=1,countero
                    if(whicho(l).ne.i)then
                    if(ssb(whicho(l)).gt.tmp1.and.
     &                 ssb(whicho(l)).le.tmp2)then
                       countern=countern+1
                       whichn(countern)=whicho(l)
                    end if
                    end if
                 end do
              end if

              loglikc=loglikc+
     &               log(2.d0*cparb*dble(je2)+dble(2*countern))-
     &               log(2.d0*cparb*dble(je2)+dble(  countero))

              if(countern.eq.0)then
                 ok=0
               else  
                 countero=countern
                 j=j+1
              end if   
           end do

1          continue

           tmp1=-(dble(nrand)*log(tpi))
           tmp2=detlogsbc
           tmp3=sseval

           loglikc=loglikc+(tmp1-tmp2-tmp3)/2.d0
         end if
      end do

      return
      end
      

      
c=======================================================================                  
      subroutine condptpriorell(ind,m,nrec,nsubject,nrand,theta,ssb,
     &                          sigmainvb,detlogsb,cparb,
     &                          whicho,whichn,
     &                          logprior,sseval)
c======================================================================= 
c     This subroutine evaluate the log-prior for a random effect
c     in a marginal Multivariate PT, using an elliptical binary 
c     partition.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none 

c-----Input
      integer ind,m,nrand,nrec,nsubject
      integer whicho(nrec),whichn(nrec)
      double precision cparb,detlogsb,theta(nrand)
      double precision ssb(nsubject)
      double precision sigmainvb(nrand,nrand)
       
c-----Output
      double precision logprior,sseval

c-----Working
      integer i,j,je2,k,k1,k2,l
      integer countero,countern,parti,nint,ok
      double precision invcdfchisq,prob,quan,tmp1,tmp2,tmp3,tpi

c-----Routine
      logprior=0.d0
      tpi=6.283185307179586476925286766559d0      
      
c++++ check if the user has requested an interrupt
      call rchkusr()
      
      tmp1=0.d0
      do i=1,nrand
         do j=1,nrand
            tmp1=tmp1+theta(i)*sigmainvb(i,j)*theta(j)
         end do
      end do
      sseval=tmp1

      nint=2
      prob=1.d0/dble(nint)
      quan=invcdfchisq(prob,dble(nrand),1,0)
      
      countero=0
      
      if(sseval.le.quan)then
          parti=1
          do i=1,nsubject
             if(ssb(i).le.quan.and.i.ne.ind)then
               countero=countero+1
               whicho(countero)=i
             end if  
          end do
        else
          parti=2
          do i=1,nsubject
             if(ssb(i).gt.quan.and.i.ne.ind)then
               countero=countero+1
               whicho(countero)=i
             end if  
          end do
      end if  
   
      logprior=logprior+
     &     log(2.d0*cparb+dble(2*countero))-
     &     log(2.d0*cparb+dble(nsubject-1))

      if(countero.eq.0) go to 1

      ok=1
      j=2
      do while(ok.eq.1.and.j.le.m)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
   
         k1=2*(parti-1)+1
         k2=2*(parti-1)+2
         quan=invcdfchisq(dble(k1)*prob,dble(nrand),1,0)

         if(sseval.le.quan)then
            parti=k1
            k=k1
           else
            parti=k2
            k=k2
         end if  

         countern=0

         if(k.eq.1)then
            do l=1,countero
               if(ssb(whicho(l)).le.quan.and.
     &            whicho(l).ne.ind)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else if(k.eq.nint)then
            quan=invcdfchisq(dble(k-1)*prob,dble(nrand),1,0)
            do l=1,countero
               if(ssb(whicho(l)).gt.quan.and.
     &            whicho(l).ne.ind)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else
            tmp1=invcdfchisq(dble(k-1)*prob,dble(nrand),1,0)
            tmp2=invcdfchisq(dble(k  )*prob,dble(nrand),1,0)

            if(tmp1.ge.tmp2)then
              call rexit("Error in the limits")
            end if  
         
            do l=1,countero
               if(whicho(l).ne.ind)then
               if(ssb(whicho(l)).gt.tmp1.and.
     &            ssb(whicho(l)).le.tmp2)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if
               end if
            end do
         end if
        
         logprior=logprior+
     &          log(2.d0*cparb*dble(je2)+dble(2*countern))-
     &          log(2.d0*cparb*dble(je2)+dble(  countero))

         if(countern.eq.0)then
            ok=0
          else  
            countero=countern
            j=j+1
         end if   
      end do

1     continue

      tmp1=-(dble(nrand)*log(tpi))
      tmp2=detlogsb
      tmp3=sseval
      
      logprior=logprior+(tmp1-tmp2-tmp3)/2.d0

      return
      end



c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR UNIVARIATE POLYA TREES
c=======================================================================                  
c=======================================================================                  

c======================================================================= 
      subroutine loglik_unippt(nsubject,mdzero,maxm,alpha,mu,sigma,b,
     &                        whicho,whichn,logliko)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for 
c     the baseline parameters in a random effect model using
c     in a marginal univariate partially specified PT.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer mdzero,maxm,nsubject
      double precision alpha
      double precision mu,sigma 
      double precision b(nsubject)

c++++ Working External
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internal
      integer countero,countern
      integer i,j,je2,k,k1,k2,l,nint,parti
      integer ok
      double precision dnrm,invcdfnorm
      double precision prob
      double precision quan
      double precision tmp1,tmp2

c++++ Output
      double precision logliko

      logliko=0.d0
      do i=1,nsubject
         call rchkusr()   
c+++++++ first observation
         if(i.eq.1)then
              logliko=dnrm(b(1),mu,sqrt(sigma),1)
c+++++++ following observations
           else
              quan=mu
              countero=0
              if(b(i).le.quan) then
                  parti=1
                  do l=1,i-1
                     if(b(l).le.quan)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                  end do
               else
                  parti=2
                  do l=1,i-1
                     if(b(l).gt.quan)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                  end do
              end if  
              if(mdzero.ne.0)then 
                 logliko=logliko+
     &            log(2.d0*alpha+dble(2*countero))-
     &            log(2.d0*alpha+dble(i-1))
              end if 

              if(countero.eq.0)go to 1
              ok=1
              j=2
              do while(ok.eq.1.and.j.le.maxm)
                 nint=2**j
                 je2=j**2
                 prob=1.d0/dble(nint)
              
                 k1=2*(parti-1)+1
                 k2=2*(parti-1)+2
                 quan=invcdfnorm(dble(k1)*prob,mu,
     &                           sqrt(sigma),1,0)
                 if(b(i).le.quan)then
                   parti=k1
                   k=k1
                  else
                   parti=k2
                    k=k2
                 end if  
                 countern=0
                 if(k.eq.1)then
                    do l=1,countero
                       if(b(whicho(l)).le.quan)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                       end if   
                    end do
                  else if(k.eq.nint)then
                    quan=invcdfnorm(dble(k-1)*prob,mu,
     &                              sqrt(sigma),1,0) 
                    do l=1,countero
                       if(b(whicho(l)).gt.quan)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                       end if   
                    end do
                  else
                    tmp1=invcdfnorm(dble(k-1)*prob,mu,
     &                              sqrt(sigma),1,0)
                    tmp2=invcdfnorm(dble(k  )*prob,mu,
     &                              sqrt(sigma),1,0)

                    if(tmp1.ge.tmp2)then
                       call rexit("Error in the limits")
                    end if  
                 
                    do l=1,countero
                       if(b(whicho(l)).gt.tmp1.and.
     &                    b(whicho(l)).le.tmp2)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                       end if   
                    end do
                 end if
              
                 logliko=logliko+
     &                 log(2.d0*alpha*dble(je2)+dble(2*countern))-
     &                 log(2.d0*alpha*dble(je2)+dble(  countero))

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
1             continue
              logliko=logliko+dnrm(b(i),mu,sqrt(sigma),1)
         end if
      end do
      return
      end

c======================================================================= 
      subroutine loglik_unifpt(nsubject,mdzero,alpha,mu,sigma,b,
     &                         whicho,whichn,logliko)
c======================================================================= 
c     This subroutine evaluate the log-likelihood for 
c     the baseline parameters in a random effect model using
c     in a marginal univariate fully specified PT.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer mdzero,nsubject
      double precision alpha
      double precision mu,sigma 
      double precision b(nsubject)

c++++ Working External
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internal
      integer countero,countern
      integer i,j,je2,k,k1,k2,l,nint,parti
      integer ok
      double precision dnrm,invcdfnorm
      double precision prob
      double precision quan
      double precision tmp1,tmp2

c++++ Output
      double precision logliko

      logliko=0.d0
      do i=1,nsubject
         call rchkusr()   
c+++++++ first observation
         if(i.eq.1)then
              logliko=dnrm(b(1),mu,sqrt(sigma),1)
c+++++++ following observations
           else
              quan=mu
              countero=0
              if(b(i).le.quan) then
                  parti=1
                  do l=1,i-1
                     if(b(l).le.quan)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                  end do
               else
                  parti=2
                  do l=1,i-1
                     if(b(l).gt.quan)then
                        countero=countero+1
                        whicho(countero)=l
                     end if   
                  end do
              end if  
              if(mdzero.ne.0)then 
                 logliko=logliko+
     &            log(2.d0*alpha+dble(2*countero))-
     &            log(2.d0*alpha+dble(i-1))
              end if 

              if(countero.eq.0)go to 1
              ok=1
              j=2
              do while(ok.eq.1)
                 nint=2**j
                 je2=j**2
                 prob=1.d0/dble(nint)
              
                 k1=2*(parti-1)+1
                 k2=2*(parti-1)+2
                 quan=invcdfnorm(dble(k1)*prob,mu,
     &                           sqrt(sigma),1,0)
                 if(b(i).le.quan)then
                   parti=k1
                   k=k1
                  else
                   parti=k2
                    k=k2
                 end if  
                 countern=0
                 if(k.eq.1)then
                    do l=1,countero
                       if(b(whicho(l)).le.quan)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                       end if   
                    end do
                  else if(k.eq.nint)then
                    quan=invcdfnorm(dble(k-1)*prob,mu,
     &                              sqrt(sigma),1,0) 
                    do l=1,countero
                       if(b(whicho(l)).gt.quan)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                       end if   
                    end do
                  else
                    tmp1=invcdfnorm(dble(k-1)*prob,mu,
     &                              sqrt(sigma),1,0)
                    tmp2=invcdfnorm(dble(k  )*prob,mu,
     &                              sqrt(sigma),1,0)

                    if(tmp1.ge.tmp2)then
                       call rexit("Error in the limits")
                    end if  
                 
                    do l=1,countero
                       if(b(whicho(l)).gt.tmp1.and.
     &                    b(whicho(l)).le.tmp2)then
                          countern=countern+1
                          whichn(countern)=whicho(l)
                       end if   
                    end do
                 end if
              
                 logliko=logliko+
     &                 log(2.d0*alpha*dble(je2)+dble(2*countern))-
     &                 log(2.d0*alpha*dble(je2)+dble(  countero))

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
1             continue
              logliko=logliko+dnrm(b(i),mu,sqrt(sigma),1)
         end if
      end do
      return
      end


c=======================================================================                  
      subroutine condupptprior(theta,ii,maxm,mdzero,nsubject,alpha,mu,
     &                         sigma,b,
     &                         whicho,whichn,logprioro)
c======================================================================= 
c     This subroutine evaluate the log-contional prior distribution,
c     arising in a marginal univariate partially specified PT, 
c     for subject 'ii' with value 'theta'.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer maxm,mdzero,nsubject,ii
      double precision alpha,mu,sigma,b(nsubject),theta

c++++ Working Externals
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internals
      integer countero,countern
      integer j,je2,k,k1,k2,l,nint,parti
      integer ok
      double precision dnrm,invcdfnorm
      double precision prob
      double precision quan
      double precision tmp1,tmp2

c++++ Output
      double precision logprioro
      
      logprioro=0.d0
      
      quan=mu
      countero=0
      if(theta.le.quan) then
          parti=1
          do l=1,nsubject
             if(b(l).le.quan.and.l.ne.ii)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
        else
          parti=2
          do l=1,nsubject
             if(b(l).gt.quan.and.l.ne.ii)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
      end if  

      if(mdzero.ne.0)then 
         logprioro=logprioro+
     &    log(2.d0*alpha+dble(2*countero))-
     &    log(2.d0*alpha+dble(nsubject-1))
      end if 

      if(countero.eq.0)go to 1

      ok=1
      j=2
      do while(ok.eq.1.and.j.le.maxm)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
        
         k1=2*(parti-1)+1
         k2=2*(parti-1)+2
         quan=invcdfnorm(dble(k1)*prob,mu,sqrt(sigma),1,0)
      
         if(theta.le.quan)then
           parti=k1
           k=k1
          else
           parti=k2
           k=k2
         end if  
         
         countern=0

         if(k.eq.1)then
            do l=1,countero
               if(b(whicho(l)).le.quan.and.
     &            whicho(l).ne.ii)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else if(k.eq.nint)then
            quan=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0) 
            do l=1,countero
               if(b(whicho(l)).gt.quan.and.
     &            whicho(l).ne.ii)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else
            tmp1=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0)
            tmp2=invcdfnorm(dble(k  )/dble(nint),mu,
     &                      sqrt(sigma),1,0)

            if(tmp1.ge.tmp2)then
              call rexit("Error in the limits")
            end if  
         
            do l=1,countero
               if(whicho(l).ne.ii)then
               if(b(whicho(l)).gt.tmp1.and.
     &            b(whicho(l)).le.tmp2)then
                 countern=countern+1
                 whichn(countern)=whicho(l)
               end if
               end if
            end do
         end if
        
         logprioro=logprioro+
     &       log(2.d0)+                
     &       log(     alpha*dble(je2)+dble(  countern))-
     &       log(2.d0*alpha*dble(je2)+dble(  countero))

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

      logprioro=logprioro+dnrm(theta,mu,sqrt(sigma),1)
      
      return
      end

c=======================================================================                  
      subroutine condufptprior(theta,ii,mdzero,nsubject,alpha,mu,
     &                         sigma,b,
     &                         whicho,whichn,logprioro)
c======================================================================= 
c     This subroutine evaluate the log-contional prior distribution,
c     arising in a marginal univariate fully specified PT, 
c     for subject 'ii' with value 'theta'.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer mdzero,nsubject,ii
      double precision alpha,mu,sigma,b(nsubject),theta

c++++ Working Externals
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internals
      integer countero,countern
      integer j,je2,k,k1,k2,l,nint,parti
      integer ok
      double precision dnrm,invcdfnorm
      double precision prob
      double precision quan
      double precision tmp1,tmp2

c++++ Output
      double precision logprioro
      
      logprioro=0.d0
      
      quan=mu
      countero=0
      if(theta.le.quan) then
          parti=1
          do l=1,nsubject
             if(b(l).le.quan.and.l.ne.ii)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
        else
          parti=2
          do l=1,nsubject
             if(b(l).gt.quan.and.l.ne.ii)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
      end if  

      if(mdzero.ne.0)then 
         logprioro=logprioro+
     &    log(2.d0*alpha+dble(2*countero))-
     &    log(2.d0*alpha+dble(nsubject-1))
      end if 

      if(countero.eq.0)go to 1

      ok=1
      j=2
      do while(ok.eq.1)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
        
         k1=2*(parti-1)+1
         k2=2*(parti-1)+2
         quan=invcdfnorm(dble(k1)*prob,mu,sqrt(sigma),1,0)
      
         if(theta.le.quan)then
           parti=k1
           k=k1
          else
           parti=k2
           k=k2
         end if  
         
         countern=0

         if(k.eq.1)then
            do l=1,countero
               if(b(whicho(l)).le.quan.and.
     &            whicho(l).ne.ii)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else if(k.eq.nint)then
            quan=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0) 
            do l=1,countero
               if(b(whicho(l)).gt.quan.and.
     &            whicho(l).ne.ii)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else
            tmp1=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0)
            tmp2=invcdfnorm(dble(k  )/dble(nint),mu,
     &                      sqrt(sigma),1,0)

            if(tmp1.ge.tmp2)then
              call rexit("Error in the limits")
            end if  
         
            do l=1,countero
               if(whicho(l).ne.ii)then
               if(b(whicho(l)).gt.tmp1.and.
     &            b(whicho(l)).le.tmp2)then
                 countern=countern+1
                 whichn(countern)=whicho(l)
               end if
               end if
            end do
         end if
        
         logprioro=logprioro+
     &       log(2.d0)+                
     &       log(     alpha*dble(je2)+dble(  countern))-
     &       log(2.d0*alpha*dble(je2)+dble(  countero))

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

      logprioro=logprioro+dnrm(theta,mu,sqrt(sigma),1)
      
      return
      end


c=======================================================================                  
      subroutine sampupptpred(maxm,mdzero,nsubject,alpha,mu,
     &                       sigma,b,
     &                       whicho,whichn,theta)
c======================================================================= 
c     This subroutine get a sample from the predictive distribution
c     from a univariate partially specified PT.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer maxm,mdzero,nsubject
      double precision alpha,mu,sigma,b(nsubject)

c++++ Working Externals
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internals
      integer countero,countern
      integer j,je2,k,k1,k2,l,nint
      integer ok
      double precision invcdfnorm,rtnorm
      double precision linf,lsup
      double precision prob
      double precision quan
      double precision tmp1,tmp2,tmp3
      
      logical ainf,binf

      real runif
      
c++++ Output
      double precision theta

      quan=mu
      countern=0
      do l=1,nsubject
         if(b(l).le.quan)then
            countern=countern+1
         end if
      end do

      if(mdzero.ne.0)then 
         tmp3=exp(
     &     log(     alpha+dble(countern))-
     &     log(2.d0*alpha+dble(nsubject)))                  
        else
         tmp3=0.5d0
      end if   
      
      countero=0
      
      if(dble(runif()).le.tmp3)then
         k=1
         do l=1,nsubject
            if(b(l).le.quan)then
               countero=countero+1
               whicho(countero)=l
            end if
         end do
       else
         k=2
         do l=1,nsubject
            if(b(l).gt.quan)then
               countero=countero+1
               whicho(countero)=l
            end if   
         end do
      end if
      
      ok=1
      j=2
      do while(ok.eq.1.and.j.le.maxm)
         
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
         
         k1=2*(k-1)+1
         k2=2*(k-1)+2
         
         countern=0
           
         if(k1.eq.1)then
            quan=invcdfnorm(dble(k1)*prob,mu,sqrt(sigma),1,0) 
            do l=1,countero
               if(b(whicho(l)).le.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else if(k1.eq.nint)then
            quan=invcdfnorm(dble(k1-1)*prob,mu,
     &                      sqrt(sigma),1,0) 
            do l=1,countero
               if(b(whicho(l)).gt.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else
            tmp1=invcdfnorm(dble(k1-1)*prob,mu,
     &                      sqrt(sigma),1,0)
            tmp2=invcdfnorm(dble(k1  )*prob,mu,
     &                      sqrt(sigma),1,0)

            if(tmp1.ge.tmp2)then
               call rexit("Error in the limits")
            end if  
              
            do l=1,countero
                 if(b(whicho(l)).gt.tmp1.and.
     &              b(whicho(l)).le.tmp2)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                 end if   
            end do
         end if
         
         tmp3=exp(
     &        log(     alpha*dble(je2)+dble(countern))-
     &        log(2.d0*alpha*dble(je2)+dble(countero)))

         if(dble(runif()).le.tmp3)then
             k=k1
           else
             k=k2
             countern=countero-countern
         end if

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
      
      if(j.gt.maxm)j=maxm

c+++++Now j indicates the partition and k the interval

      nint=2**j
      prob=1.d0/dble(nint)

      if(k.eq.1)then
         quan=invcdfnorm(dble(k)*prob,mu,sqrt(sigma),1,0) 

         ainf=.true.
         binf=.false.
         linf=0.d0
         lsup=quan
         theta=rtnorm(mu,sqrt(sigma),linf,lsup,ainf,binf)
         
       else if(k.eq.nint)then
         quan=invcdfnorm(dble(k-1)*prob,mu,sqrt(sigma),1,0) 

         ainf=.false.
         binf=.true.
         linf=quan
         lsup=0.d0
         theta=rtnorm(mu,sqrt(sigma),linf,lsup,ainf,binf)

       else
         tmp1=invcdfnorm(dble(k-1)*prob,mu,sqrt(sigma),1,0)
         tmp2=invcdfnorm(dble(k  )*prob,mu,sqrt(sigma),1,0)

         if(tmp1.ge.tmp2)then
            call rexit("Error in the limits")
         end if  

         ainf=.false.
         binf=.false.
         linf=tmp1
         lsup=tmp2
         theta=rtnorm(mu,sqrt(sigma),linf,lsup,ainf,binf)
      end if
      
      return
      end

c=======================================================================                  
      subroutine sampufptpred(mdzero,nsubject,alpha,mu,
     &                        sigma,b,
     &                        whicho,whichn,theta)
c======================================================================= 
c     This subroutine get a sample from the predictive distribution
c     from a univariate fully specified PT.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer mdzero,nsubject
      double precision alpha,mu,sigma,b(nsubject)

c++++ Working Externals
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internals
      integer countero,countern
      integer j,je2,k,k1,k2,l,nint
      integer ok
      double precision invcdfnorm,rtnorm
      double precision linf,lsup
      double precision prob
      double precision quan
      double precision tmp1,tmp2,tmp3
      
      logical ainf,binf

      real runif
      
c++++ Output
      double precision theta

      quan=mu
      countern=0
      do l=1,nsubject
         if(b(l).le.quan)then
            countern=countern+1
         end if
      end do

      if(mdzero.ne.0)then 
         tmp3=exp(
     &     log(     alpha+dble(countern))-
     &     log(2.d0*alpha+dble(nsubject)))                  
        else
         tmp3=0.5d0
      end if   
      
      countero=0
      
      if(dble(runif()).le.tmp3)then
         k=1
         do l=1,nsubject
            if(b(l).le.quan)then
               countero=countero+1
               whicho(countero)=l
            end if
         end do
       else
         k=2
         do l=1,nsubject
            if(b(l).gt.quan)then
               countero=countero+1
               whicho(countero)=l
            end if   
         end do
      end if
      
      ok=1
      j=2
      do while(ok.eq.1)
         
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
         
         k1=2*(k-1)+1
         k2=2*(k-1)+2
         
         countern=0
           
         if(k1.eq.1)then
            quan=invcdfnorm(dble(k1)*prob,mu,sqrt(sigma),1,0) 
            do l=1,countero
               if(b(whicho(l)).le.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else if(k1.eq.nint)then
            quan=invcdfnorm(dble(k1-1)*prob,mu,
     &                      sqrt(sigma),1,0) 
            do l=1,countero
               if(b(whicho(l)).gt.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else
            tmp1=invcdfnorm(dble(k1-1)*prob,mu,
     &                      sqrt(sigma),1,0)
            tmp2=invcdfnorm(dble(k1  )*prob,mu,
     &                      sqrt(sigma),1,0)

            if(tmp1.ge.tmp2)then
               call rexit("Error in the limits")
            end if  
              
            do l=1,countero
                 if(b(whicho(l)).gt.tmp1.and.
     &              b(whicho(l)).le.tmp2)then
                    countern=countern+1
                    whichn(countern)=whicho(l)
                 end if   
            end do
         end if
         
         tmp3=exp(
     &        log(     alpha*dble(je2)+dble(countern))-
     &        log(2.d0*alpha*dble(je2)+dble(countero)))

         if(dble(runif()).le.tmp3)then
             k=k1
           else
             k=k2
             countern=countero-countern
         end if

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
      
c+++++Now j indicates the partition and k the interval

      nint=2**j
      prob=1.d0/dble(nint)

      if(k.eq.1)then
         quan=invcdfnorm(dble(k)*prob,mu,sqrt(sigma),1,0) 

         ainf=.true.
         binf=.false.
         linf=0.d0
         lsup=quan
         theta=rtnorm(mu,sqrt(sigma),linf,lsup,ainf,binf)
         
       else if(k.eq.nint)then
         quan=invcdfnorm(dble(k-1)*prob,mu,sqrt(sigma),1,0) 

         ainf=.false.
         binf=.true.
         linf=quan
         lsup=0.d0
         theta=rtnorm(mu,sqrt(sigma),linf,lsup,ainf,binf)

       else
         tmp1=invcdfnorm(dble(k-1)*prob,mu,sqrt(sigma),1,0)
         tmp2=invcdfnorm(dble(k  )*prob,mu,sqrt(sigma),1,0)

         if(tmp1.ge.tmp2)then
            call rexit("Error in the limits")
         end if  

         ainf=.false.
         binf=.false.
         linf=tmp1
         lsup=tmp2
         theta=rtnorm(mu,sqrt(sigma),linf,lsup,ainf,binf)
      end if
      
      return
      end

c=======================================================================
      subroutine locationupt(x,m,n,loca,mu,sigma)
c=======================================================================
c     function that return the succesive location of x across the 
c     finite univariate tree of length m. This is based on a 
c     normal(mu,sigma) distribution.
c
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer m,n,loca(n)
      integer i,k,k1,k2,nint
      double precision invcdfnorm
      double precision prob,quan
      double precision mu,sigma
      double precision x
      
      if(m.gt.n)then
        call rexit("Error in 'locationpt'")
      end if

      quan=mu
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
         quan=invcdfnorm(dble(k1)*prob,mu,sqrt(sigma),1,0)
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
      subroutine accumcountupt(nsubject,mu,sigma,m,b,maxm,loca,
     &                         ntotals,nhash,maxnzr,tmp,counts,nr,
     &                         ia,ja,a)
c=======================================================================
c     accumulate the number of subjects.  
c
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer i,j
      integer m,nsubject,maxm
      integer loca(maxm)
      double precision mu,sigma
      double precision b(nsubject),theta,one
      parameter(one=1.d0)
      
      integer nhash,maxnzr,nr
      double precision counts(nhash,3)
      integer ntotals,ia(ntotals+1),ja(maxnzr),tmp(ntotals)
      double precision a(maxnzr)

      nr=0
      do i=1,nsubject
         theta=b(i)
         call locationupt(theta,m,maxm,loca,mu,sigma)         
         do j=1,m
            call hashm(one,loca(j),j,counts,nhash,nr)
         end do
      end do
      
      call hashiajaa(counts,nhash,ntotals,ia,ja,a,maxnzr,tmp)
      
      return
      end


c=======================================================================
      subroutine sprobupt(mdzero,cpar,m,ntotals,maxnzr,
     &                    ia,ja,a,tmp1,tmp2)
c=======================================================================
c     generate probabilities for a finite univariate PT.  
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer mdzero,i,j,k,l,ista,iend
      integer m,narea,ntotals,maxnzr
      integer ia(ntotals+1),ja(maxnzr)
      double precision a(maxnzr)
      double precision cpar
      double precision mass(2)
      double precision rvecs(2)
      double precision tmp1(ntotals),tmp2(ntotals)

      integer ngroup,nvar
      double precision je2
      
      narea=2 
      nvar=1
      
c++++ First level probabilities
      ngroup=1
      do k=1,narea
         mass(k)=cpar
      end do

      ista=1
      iend=2
      
      if(mdzero.ne.0)then
         do k=ista,iend
            do l=ia(k),ia(k+1)-1
               if(ja(l).eq.1)mass(k-ista+1)=mass(k-ista+1)+a(l)
            end do
         end do
         call dirichlet(mass,2,2,rvecs)
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
            call dirichlet(mass,2,2,rvecs)
            
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
      subroutine setsupt(m,ntotals,liminf1,liminf2,limsup1,limsup2)
c=======================================================================
c     computes the limtis of the sets in a finite univariate PT
c
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer m,ntotals
      integer pattern(1)
      double precision liminf1(ntotals),liminf2(ntotals)
      double precision limsup1(ntotals),limsup2(ntotals)
      integer i,j,k,ngroup,ista,narea,nvar
      double precision linf,lsup,quan

c++++ Algorithm
      nvar=1
      narea=2**nvar  
      
c++++ First level sets

      do i=1,narea
         call binaryrepinv(nvar,i,pattern)
         if(pattern(1).eq.0)then
            liminf1(i)=-999999.d0
            limsup1(i)=0.d0
            liminf2(i)=-999999.d0
            limsup2(i)=0.d0
          else
            liminf1(i)=0.d0
            limsup1(i)=+999999.d0
            liminf2(i)=0.d0
            limsup2(i)=+999999.d0
         end if 
      end do

c++++ Sets for other levels

      do i=2,m
         ngroup=2**(nvar*(i-1))
         do j=1,ngroup

            ista=(2**nvar)*(j-1)+1

            do k=1,narea
               call binaryrepinv(nvar,k,pattern)

               linf=liminf1(j)
               lsup=limsup1(j)
               call quandetpt(linf,lsup,quan)
                  
               if(pattern(1).eq.0)then
                  liminf2(ista+k-1)=linf
                  limsup2(ista+k-1)=quan
                else
                  liminf2(ista+k-1)=quan
                  limsup2(ista+k-1)=lsup
               end if 
               
            end do
         end do
         
         if(i.lt.m)then
            do j=1,2**(nvar*i)
               liminf1(j)=liminf2(j)
               limsup1(j)=limsup2(j)
            end do
         end if   
      end do
      return
      end


c=======================================================================
      subroutine momupt(m,ntotals,liminf,limsup,probs,
     &                  means,covs)
c=======================================================================
c     computes the moments of a finite univariate PT
c
c     Alejandro Jara, 2007
c=======================================================================
      implicit none
      integer m,ntotals
      double precision liminf(ntotals)
      double precision limsup(ntotals)
      double precision probs(ntotals)
      double precision means,covs
      
      integer i,nsets,nvar
      double precision linf,lsup,muwork1,muwork2
      double precision cdfnorm,dnrm
      double precision tmp1,tmp2,tmp3,tmp4
      double precision tmass

c++++ Algorithm
      nvar=1
      nsets=2**(nvar*m)  
      means=0.d0 
      covs=0.d0 
      
      if(nsets.gt.ntotals)then
         call rexit("Error in 'mompt': nsets > ntotals")            
      end if
      
      tmass=0.d0
      do i=1,nsets
         linf=liminf(i)
         lsup=limsup(i)
         tmp1=dnrm(linf,0.d0,1.d0,0)
         tmp2=dnrm(lsup,0.d0,1.d0,0)
         tmp3=cdfnorm(linf,0.d0,1.d0,1,0)
         tmp4=cdfnorm(lsup,0.d0,1.d0,1,0)
         muwork1= - (tmp2-tmp1)/(tmp4-tmp3)
         means=means+probs(i)*muwork1
            
         muwork2=1.d0 - (lsup*tmp2-linf*tmp1)/(tmp4-tmp3) 
         covs=covs+probs(i)*muwork2

         tmass=tmass+probs(i)
      end do

      covs=covs-means*means   

      return
      end


c=======================================================================
      subroutine samplefuncupt(mdzero,m,nsubject,cpar,b,mu,sigma,
     &                         means,covs)
c=======================================================================
      implicit none
      integer i,j
      integer mdzero,m,nsubject
      double precision cpar,b(nsubject)
      double precision mu,sigma
      double precision means,covs

c++++ parameters
      integer maxm,ntotals,ntotalp,nhash,maxnzr
      parameter(maxm=20)
      parameter(ntotals=2**(2*6),ntotalp=5460)
      parameter(nhash=ntotalp,maxnzr=ntotalp)
      integer loca(maxm),mwork,nr
      integer ia(ntotals+1),ja(maxnzr),tmp(ntotals)

      double precision a(maxnzr)    
      double precision counts(nhash,3)
      double precision probw1(ntotals),probw2(ntotals)
      double precision liminf1(ntotals),liminf2(ntotals)
      double precision limsup1(ntotals),limsup2(ntotals)
      
      mwork=m
      do while(ntotals.lt.(2**mwork))
         mwork=mwork-1
      end do

      if(mwork.lt.1)then
        call rexit("Error in 'samplefuncupt': maxm < 1")      
      end if
      
      if(mwork.gt.maxm)then
        call rexit("Error in 'samplefuncupt': increase maxm")      
      end if

      nr=0
      do i=1,nhash
         do j=1,3
            counts(i,j)=0.d0
         end do   
      end do

c++++ set the MPT parameters      

      call accumcountupt(nsubject,mu,sigma,m,b,maxm,loca,
     &                   ntotals,nhash,maxnzr,tmp,counts,nr,
     &                   ia,ja,a)

c++++ sampling probabilities

      call sprobupt(mdzero,cpar,mwork,ntotals,maxnzr,
     &              ia,ja,a,probw1,probw2)

c++++ finding sets

      call setsupt(mwork,ntotals,liminf1,liminf2,limsup1,limsup2)

c++++ computing the means for each set

      call momupt(mwork,ntotals,liminf2,limsup2,probw2,
     &            means,covs)

      return
      end


c=======================================================================                  
      subroutine gridupptprior(theta,maxm,mdzero,nsubject,alpha,mu,
     &                         sigma,b,
     &                         whicho,whichn,logprioro)
c======================================================================= 
c     This subroutine evaluate the log-contional prior distribution,
c     arising in a marginal univariate partially specified PT, 
c     for a value in a grid 'theta'.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer maxm,mdzero,nsubject
      double precision alpha,mu,sigma,b(nsubject),theta

c++++ Working Externals
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internals
      integer countero,countern
      integer j,je2,k,k1,k2,l,nint,parti
      integer ok
      double precision dnrm,invcdfnorm
      double precision prob
      double precision quan
      double precision tmp1,tmp2

c++++ Output
      double precision logprioro
      
      logprioro=0.d0
      
      quan=mu
      countero=0
      if(theta.le.quan) then
          parti=1
          do l=1,nsubject
             if(b(l).le.quan)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
        else
          parti=2
          do l=1,nsubject
             if(b(l).gt.quan)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
      end if  

      if(mdzero.ne.0)then 
         logprioro=logprioro+
     &    log(2.d0*alpha+dble(2*countero))-
     &    log(2.d0*alpha+dble(nsubject))
      end if 

      if(countero.eq.0)go to 1

      ok=1
      j=2
      do while(ok.eq.1.and.j.le.maxm)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
        
         k1=2*(parti-1)+1
         k2=2*(parti-1)+2
         quan=invcdfnorm(dble(k1)*prob,mu,sqrt(sigma),1,0)
      
         if(theta.le.quan)then
           parti=k1
           k=k1
          else
           parti=k2
           k=k2
         end if  
         
         countern=0

         if(k.eq.1)then
            do l=1,countero
               if(b(whicho(l)).le.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else if(k.eq.nint)then
            quan=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0) 
            do l=1,countero
               if(b(whicho(l)).gt.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else
            tmp1=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0)
            tmp2=invcdfnorm(dble(k  )/dble(nint),mu,
     &                      sqrt(sigma),1,0)

            if(tmp1.ge.tmp2)then
              call rexit("Error in the limits")
            end if  
         
            do l=1,countero
               if(b(whicho(l)).gt.tmp1.and.
     &            b(whicho(l)).le.tmp2)then
                 countern=countern+1
                 whichn(countern)=whicho(l)
               end if
            end do
         end if
        
         logprioro=logprioro+
     &       log(2.d0)+                
     &       log(     alpha*dble(je2)+dble(  countern))-
     &       log(2.d0*alpha*dble(je2)+dble(  countero))

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

      logprioro=logprioro+dnrm(theta,mu,sqrt(sigma),1)
      
      return
      end


c=======================================================================                  
      subroutine gridufptprior(theta,mdzero,nsubject,alpha,mu,
     &                         sigma,b,
     &                         whicho,whichn,logprioro)
c======================================================================= 
c     This subroutine evaluate the log-contional prior distribution,
c     arising in a marginal univariate fully specified PT, 
c     for a value in a grid 'theta'.
c
c     Alejandro Jara, 2007
c======================================================================= 
      implicit none

c++++ Input
      integer mdzero,nsubject
      double precision alpha,mu,sigma,b(nsubject),theta

c++++ Working Externals
      integer whicho(nsubject),whichn(nsubject)

c++++ Working Internals
      integer countero,countern
      integer j,je2,k,k1,k2,l,nint,parti
      integer ok
      double precision dnrm,invcdfnorm
      double precision prob
      double precision quan
      double precision tmp1,tmp2

c++++ Output
      double precision logprioro
      
      logprioro=0.d0
      
      quan=mu
      countero=0
      if(theta.le.quan) then
          parti=1
          do l=1,nsubject
             if(b(l).le.quan)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
        else
          parti=2
          do l=1,nsubject
             if(b(l).gt.quan)then
                countero=countero+1
                whicho(countero)=l
             end if   
          end do
      end if  

      if(mdzero.ne.0)then 
         logprioro=logprioro+
     &    log(2.d0*alpha+dble(2*countero))-
     &    log(2.d0*alpha+dble(nsubject))
      end if 

      if(countero.eq.0)go to 1

      ok=1
      j=2
      do while(ok.eq.1)
         nint=2**j
         je2=j**2
         prob=1.d0/dble(nint)
        
         k1=2*(parti-1)+1
         k2=2*(parti-1)+2
         quan=invcdfnorm(dble(k1)*prob,mu,sqrt(sigma),1,0)
      
         if(theta.le.quan)then
           parti=k1
           k=k1
          else
           parti=k2
           k=k2
         end if  
         
         countern=0

         if(k.eq.1)then
            do l=1,countero
               if(b(whicho(l)).le.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else if(k.eq.nint)then
            quan=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0) 
            do l=1,countero
               if(b(whicho(l)).gt.quan)then
                  countern=countern+1
                  whichn(countern)=whicho(l)
               end if   
            end do
          else
            tmp1=invcdfnorm(dble(k-1)/dble(nint),mu,
     &                      sqrt(sigma),1,0)
            tmp2=invcdfnorm(dble(k  )/dble(nint),mu,
     &                      sqrt(sigma),1,0)

            if(tmp1.ge.tmp2)then
              call rexit("Error in the limits")
            end if  
         
            do l=1,countero
               if(b(whicho(l)).gt.tmp1.and.
     &            b(whicho(l)).le.tmp2)then
                 countern=countern+1
                 whichn(countern)=whicho(l)
               end if
            end do
         end if
        
         logprioro=logprioro+
     &       log(2.d0)+                
     &       log(     alpha*dble(je2)+dble(  countern))-
     &       log(2.d0*alpha*dble(je2)+dble(  countero))

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

      logprioro=logprioro+dnrm(theta,mu,sqrt(sigma),1)
      
      return
      end

