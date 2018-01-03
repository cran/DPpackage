c=======================================================================                      
      subroutine ldpddoublyintsba(                                            
     &                  nsubject,nvar,p,q,ngrid,npred,                  
     &                  x,llower,lupper,grid,xpred,                     
     &                  maxm,a0b0,nu,tinv1,smu,psiinv,tinv2,            
     &                  mcmc,nsave,                                       
     &                  ncluster,ss,alpha,b,                            
     &                  sigma,sigmainv,mub,sigmab,sigmabinv,            
     &                  y,    
     &                  acrate,f,h,thetasave,randsave,           
     &                  seed,                                           
     &                  model,possi,varind,                                    
     &                  cstrt,ccluster,prob,prob2,weights,              
     &                  iflagc,theta,workmc,workmc2,workmc3,                 
     &                  workmhc,workmhc2,workvc,                         
     &                  iflagn,workmn,workmn2,workmn3,                  
     &                  workmhn,workmhn2,workvn,workvn2,workvn3,        
     &                  ztz,zty,fw,fw2)                                 
c=======================================================================                      
c
c     # of arguments = 65.
c
c     Subroutine `ldpddoublyintsba' to run a Markov chain for a 
c     LDPD model for doubly-interval-censored data using
c     a stick-breaking representation.
c
c     Copyright: Alejandro Jara, 2008-2010.
c
c     Version 1.0: 
c
c     Last modification: 23-01-2010.
c     
c     This program is free software; you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published by
c     the Free Software Foundation; either version 2 of the License, or (at
c     your option) any later version.
c
c     This program is distributed in the hope that it will be useful, but
c     WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     General Public License for more details.
c
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
c
c     The author's contact information:
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

      implicit none

c+++++Data
      integer nsubject,nvar,p,q,ngrid,npred
      double precision x(nsubject,nvar*(p+q))
      double precision xpred(npred,nvar*(p+q))
      double precision llower(nsubject,2*nvar)
      double precision lupper(nsubject,2*nvar)
      double precision grid(2*nvar,ngrid)

c+++++Prior 
      integer maxm,nu01,nu02,nu(2) 
      double precision aa0,ab0,mubp,sigmabp,qq,a0b0(7) 
      double precision tinv1(2*nvar,2*nvar)
      double precision smu(nvar*(p+q)),psiinv(nvar*(p+q),nvar*(p+q))
      double precision tinv2(nvar*(p+q),nvar*(p+q))

c+++++MCMC parameters
      integer mcmc(3),nburn,nskip,nsave,ndisplay
      double precision tune1,tune2

c+++++Current values of the parameters
      integer ncluster,ss(nsubject)
      double precision alpha(2),ap,bp
      double precision b(maxm,nvar*(p+q))
      double precision sigma(2*nvar,2*nvar),sigmainv(2*nvar,2*nvar)
      double precision mub(nvar*(p+q))
      double precision sigmab(nvar*(p+q),nvar*(p+q))
      double precision sigmabinv(nvar*(p+q),nvar*(p+q))
      double precision y(nsubject,2*nvar)

c+++++Output
      double precision acrate(2)
      double precision f(npred*2*nvar,ngrid)
      double precision h(npred*2*nvar,ngrid)
      double precision thetasave(nsave,2*nvar*(2*nvar+1)/2+ 
     &                 nvar*(p+q)+nvar*(p+q)*(nvar*(p+q)+1)/2+3)

      double precision randsave(nsave,npred*2*nvar)

c+++++Seeds
      integer seed(2),seed1,seed2

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++External working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++Model Specification
      integer model(2*nvar,p+q)
      integer possi(2*nvar,p+q)
      integer varind(nvar*(p+q))

c+++++Clusters (maxm)
      integer cstrt(maxm,nsubject)
      integer ccluster(maxm)
      double precision prob(maxm)
      double precision prob2(maxm)
      double precision weights(maxm)

c+++++Locations
      integer iflagc(nvar*(p+q))
      double precision theta(nvar*(p+q))
      double precision workmc(nvar*(p+q),nvar*(p+q))
      double precision workmc2(nvar*(p+q),nvar*(p+q))
      double precision workmc3(nvar*(p+q),nvar*(p+q))
      double precision workmhc(nvar*(p+q)*(nvar*(p+q)+1)/2)
      double precision workmhc2(nvar*(p+q)*(nvar*(p+q)+1)/2)
      double precision workvc(nvar*(p+q))

      integer iflagn(2*nvar)
      double precision workmn(2*nvar,2*nvar)      
      double precision workmn2(2*nvar,2*nvar)      
      double precision workmn3(2*nvar,2*nvar)      
      double precision workmhn(2*nvar*(2*nvar+1)/2)
      double precision workmhn2(2*nvar*(2*nvar+1)/2)
      double precision workvn(2*nvar)
      double precision workvn2(2*nvar)
      double precision workvn3(2*nvar)
      double precision ztz(nvar*(p+q),nvar*(p+q))
      double precision zty(nvar*(p+q))
      
      double precision fw(ngrid)
      double precision fw2(ngrid)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++Internal working space
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c+++++General
      integer i,ii,j,k,l,sprint
      double precision tmp1,tmp2,tmp3,tmp4,ratio
      double precision astar

      double precision acan,bcan,uni
      double precision logpn,logpo,logcgkn,logcgko

c+++++Latent variables
      double precision slow,supp  
      logical ainf,binf

c+++++DP
      double precision lsweight

c+++++MCMC
      integer dispcount,isave,iscan,nscan,skipcount 

c+++++RNG and distributions
      double precision dnrm,dlnrm,cdfnorm
      double precision cdflnorm,rtnorm,rbeta
      real runif

c+++++CPU time
      double precision sec00,sec0,sec1,sec

c++++ parameters
      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)

      aa0=a0b0(1)
      ab0=a0b0(2)
      qq=a0b0(3)
      mubp=a0b0(4)
      sigmabp=a0b0(5)
      tune1=a0b0(6)
      tune2=a0b0(7)

      nu01=nu(1)
      nu02=nu(2)
      
      ap=alpha(1)
      bp=alpha(2)

c+++++Model Specification
      do i=1,2*nvar
         do j=1,p+q
            model(i,j)=0
         end do
      end do

      k=0      
      do i=1,nvar
         do j=1,p
            k=k+1 
            model(i,j)=1
            possi(i,j)=k
         end do
      end do

      do i=nvar+1,2*nvar
         do j=p+1,p+q
            k=k+1 
            model(i,j)=1
            possi(i,j)=k
         end do
      end do

c++++ set random number generator
      seed1=seed(1)
      seed2=seed(2)
      call setall(seed1,seed2)

c++++ opening files

      open(unit=1,file='out_dppackage1.txt',status='unknown',
     &     form='unformatted')

      open(unit=2,file='out_dppackage2.txt',status='unknown',
     &     form='unformatted')

      open(unit=3,file='out_dppackage3.txt',status='unknown',
     &     form='unformatted')

      open(unit=4,file='out_dppackage4.txt',status='unknown',
     &     form='unformatted')

      open(unit=5,file='out_dppackage5.txt',status='unknown',
     &     form='unformatted')

c++++ set configurations
      do i=1,nsubject
         if(ss(i).gt.maxm)then
            call rexit("number of clusters bigger than 'maxm'")
          else
            ccluster(ss(i))=ccluster(ss(i))+1
            cstrt(ss(i),ccluster(ss(i)))=i
         end if   
      end do
      
      astar=rbeta(aa0,ab0)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*nsave

      call cpu_time(sec0)
      sec00=0.d0
      
      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Latent data
c+++++++++++++++++++++++++++++++++++++++++++++++++

c         call dblepr("Latent",-1,dble(iscan),1)
         
         do i=1,nsubject

            do k=1,2*nvar
               tmp1=0.d0
               do l=1,p+q
                  if(model(k,l).eq.1)then
                     tmp1=tmp1+x(i,possi(k,l))*b(ss(i),possi(k,l))
                  end if   
               end do
               workvn(k)=tmp1
            end do   
            
            do j=1,2*nvar

               if(llower(i,j).ne.lupper(i,j))then

                  call condmvn(j,sigma,2*nvar,workmn,workmn2)
                  tmp1=sqrt(workmn(1,1))

                  do k=1,2*nvar
                     tmp2=0.d0
                     do l=1,2*nvar
                        tmp2=tmp2+workmn2(k,l)*workvn(l) 
                     end do
                     workvn2(k)=tmp2

                     tmp2=0.d0
                     do l=1,2*nvar
                        tmp2=tmp2+workmn2(k,l)*y(i,l) 
                     end do
                     workvn3(k)=tmp2
                  end do

                  tmp2=workvn(j)
                  do k=2,2*nvar
                     tmp2=tmp2-workmn(1,k)*(workvn3(k)-workvn2(k))
                  end do            

                  ainf=.false.
                  binf=.false. 
        
                  if(j.le.nvar)then

                     tmp3=(
     &                     max(
     &                         llower(i,j),
     &                         llower(i,nvar+j)-
     &                         exp(y(i,nvar+j))
     &                        )
     &                 )

                     tmp4=(
     &                     min(
     &                         lupper(i,j),
     &                         lupper(i,nvar+j)-
     &                         exp(y(i,nvar+j))
     &                        )
     &                 )


                     slow=log(tmp3)
                     supp=log(tmp4)
 
                     if(slow.gt.supp)then
                        call dblepr("slow",-1,slow,1) 
                        call dblepr("supp",-1,supp,1) 
                        call rexit("error in limits")
                     end if

                   else
                  
                     if(llower(i,j).eq.llower(i,j-nvar))then
                        tmp3=0.d0
                        slow=0.d0
                        ainf=.true.
                        tmp4=lupper(i,j) - exp(y(i,j-nvar)) 
                        supp=log(tmp4)
                       else
                        tmp3=llower(i,j) - exp(y(i,j-nvar))
                        tmp4=lupper(i,j) - exp(y(i,j-nvar)) 
                        slow=log(tmp3)
                        supp=log(tmp4)

                       if(slow.gt.supp)then
                          call dblepr("slow",-1,slow,1) 
                          call dblepr("supp",-1,supp,1) 
                          call rexit("error in limits")
                       end if
                     end if   
                  end if

                  y(i,j)=rtnorm(tmp2,tmp1,slow,supp,ainf,binf)

c++++++++++++++++ testing

                  if(j.le.nvar)then
                 
                  if(exp(y(i,j)).gt.lupper(i,j).or.
     &               exp(y(i,j)).lt.llower(i,j))then
                     call dblepr("expY",-1,exp(y(i,j)),1) 
                     call dblepr("low",-1,llower(i,j),1) 
                     call dblepr("upp",-1,lupper(i,j),1) 
               call rexit("error in latent variable generation 1")
                  end if
                 
                  else

               if(exp(y(i,j))+exp(y(i,j-nvar)).gt.lupper(i,j).or.
     &            exp(y(i,j))+exp(y(i,j-nvar)).lt.llower(i,j))then
                  call dblepr("expY",-1,exp(y(i,j))+exp(y(i,j-nvar)),1) 
                  call dblepr("low",-1,llower(i,j),1) 
                  call dblepr("upp",-1,lupper(i,j),1) 
               call rexit("error in latent variable generation 2")
                  end if

                  end if

               end if

            end do
         end do

c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Probabilities
c+++++++++++++++++++++++++++++++++++++++++++++++++

c         call dblepr("Probabilities",-1,dble(iscan),1)

         call spyprob_sba(maxm,ccluster,prob,weights,
     &                    ap,bp,lsweight)
 
   
c+++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++ Confingurations
c+++++++++++++++++++++++++++++++++++++++++++++++++

c         call dblepr("Config",-1,dble(iscan),1)

         call spyconf_sba(maxm,nsubject,nvar,p,q,
     &                    ncluster,ccluster,cstrt,ss,
     &                    model,possi,y,x,b,sigma,
     &                    prob,iflagn,prob2,
     &                    workmn,workvn)

c         call intpr("ncluster",-1,ncluster,1)
c         call intpr("ccluster",-1,ccluster,maxm)
c         call dblepr("prob",-1,prob,maxm)
c         call dblepr("weights",-1,weights,maxm)

c++++++++++++++++++++++++++++++
c+++++++ Sampling locations
c++++++++++++++++++++++++++++++

c         call dblepr("Location",-1,dble(iscan),1)
         
c         call spylocnc_sba(nsubject,nvar,p,q,
c     &                     maxm,ccluster,cstrt,
c     &                     model,possi,y,x,sigmainv,
c     &                     mub,sigmab,sigmabinv,b,
c     &                     iflagc,theta,
c     &                     workvc,workmhc,
c     &                     ztz,zty)

          call spylocnc_sba2(nsubject,nvar,p,q,
     &                       maxm,ccluster,cstrt,
     &                       varind,y,x,sigmainv,
     &                       mub,sigmab,sigmabinv,b,
     &                       iflagc,theta,
     &                       workvc,workmhc,
     &                       ztz,zty)


c++++++++++++++++++++++++++++++
c+++++++ Kernel Covariance
c++++++++++++++++++++++++++++++

c         call dblepr("Kernel",-1,dble(iscan),1)

         call spykcnc(maxm,nsubject,nvar,p,q,ss,
     &                model,possi,
     &                y,x,b,nu01,tinv1,sigma,sigmainv,
     &                iflagn,workvn,workmn,workmn2,workmn3,
     &                workmhn,workmhn2)

c++++++++++++++++++++++++++++++++++         
c+++++++ Base line distribution
c++++++++++++++++++++++++++++++++++

c         call dblepr("Mean",-1,dble(iscan),1)

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,nvar*(p+q)
            zty(i)=smu(i)
            do j=1,nvar*(p+q)
               ztz(i,j)=(sigmabinv(i,j)*dble(maxm))+psiinv(i,j)
            end do
         end do

         call inverse(ztz,nvar*(p+q),iflagc) 

         do i=1,maxm
            do j=1,nvar*(p+q)
               tmp1=0.d0
               do k=1,nvar*(p+q)
                  tmp1=tmp1+sigmabinv(j,k)*b(i,k)
               end do
               zty(j)=zty(j)+tmp1
            end do
         end do

         do i=1,nvar*(p+q)
            tmp1=0.d0
            do j=1,nvar*(p+q)
               tmp1=tmp1+ztz(i,j)*zty(j)
            end do
            workvc(i)=tmp1
         end do

         call rmvnorm(nvar*(p+q),workvc,ztz,workmhc,zty,theta)

c         call dblepr("mub",-1,mub,nvar*(p+q))

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c         call dblepr("Covariance",-1,dble(iscan),1)

         do i=1,nvar*(p+q)
            mub(i)=theta(i)
            do j=1,nvar*(p+q)
               workmc(i,j)=0.d0
            end do
         end do

         do i=1,maxm
            do j=1,nvar*(p+q)
               do k=1,nvar*(p+q)
                  workmc(j,k)=workmc(j,k)+               
     &                 (b(i,j)-mub(j))*(b(i,k)-mub(k))
               end do
            end do
         end do

         do i=1,nvar*(p+q)
            do j=1,nvar*(p+q)
               workmc(i,j)=workmc(i,j)+tinv2(i,j)
            end do
         end do

         call riwishart(nvar*(p+q),nu02+maxm,workmc,workmc2,workmc3,
     &                  workvc,workmhc,workmhc2,iflagc)

         do i=1,nvar*(p+q)
            do j=1,nvar*(p+q)
               sigmab(i,j)=workmc(i,j)
               sigmabinv(i,j)=workmc2(i,j)
            end do
         end do

c         call dblepr("sigmab",-1,sigmab,nvar*(p+q)*nvar*(p+q))


c+++++++++++++++++++++++++++++++++++++++++++++++++         
c+++++++ Parameters of the PD process
c+++++++++++++++++++++++++++++++++++++++++++++++++

c         call dblepr("Precision",-1,dble(iscan),1)

         if(aa0.gt.0.d0)then

c++++++++++ generates the candidate

            uni=dble(runif())
            if(uni.le.0.5d0)then
               acan=0.d0
             else
               ainf=.false.
               binf=.false. 
               acan=rtnorm(alpha(1),tune1,0.d0,1.d0,ainf,binf)
            end if
 
            ainf=.false.
            binf=.true. 
            bcan=rtnorm(alpha(2),tune2,-acan,0.d0,ainf,binf)

c++++++++++ evaluates the candidate generating distribution

            tmp1=0.d0  
            if(alpha(1).eq.0.d0)then
               tmp1=tmp1+0.5d0
            end if

            tmp2=dnrm(alpha(1),acan,tune1,1)-log(
     &           cdfnorm(1.d0,acan,tune1,1,0)-
     &           cdfnorm(0.d0,acan,tune1,1,0))
            logcgkn=log(tmp1+0.5d0*exp(tmp2))

            logcgkn=logcgkn+dnrm(alpha(2),bcan,tune2,1)-
     &           cdfnorm(-acan,bcan,tune2,0,1)


            tmp1=0.d0  
            if(acan.eq.0.d0)then
               tmp1=tmp1+0.5d0
            end if

            tmp2=dnrm(acan,alpha(1),tune1,1)-log(
     &           cdfnorm(1.d0,alpha(1),tune1,1,0)-
     &           cdfnorm(0.d0,alpha(1),tune1,1,0))
            logcgko=log(tmp1+0.5d0*exp(tmp2))

            logcgko=logcgko+dnrm(bcan,alpha(2),tune2,1)-
     &           cdfnorm(-alpha(1),alpha(2),tune2,0,1)

c++++++++++ evaluates the log-posterior

            call spydens_sbapi(maxm,weights,qq,acan,bcan,
     &                         mubp,sigmabp,aa0,ab0,logpn)

            call spydens_sbapi(maxm,weights,qq,alpha(1),alpha(2),
     &                         mubp,sigmabp,aa0,ab0,logpo)

c++++++++++ acceptance step

            ratio=logpn-logpo+logcgkn-logcgko
            
            if(log((runif())).le.ratio)then
               alpha(1)=acan
               alpha(2)=bcan
               acrate(1)=acrate(1)+1.d0
            end if
         end if
         
c         call dblepr("alpha",-1,alpha,2)
c         call dblepr("acan",-1,acan,1)
c         call dblepr("bcan",-1,bcan,1)

c++++++++++++++++++++++++++++++++++         
c+++++++ save samples
c++++++++++++++++++++++++++++++++++         
         
         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then

               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ random effects (in a file)
               do i=1,maxm
                  write(1) (b(i,j),j=1,nvar*(p+q))
               end do   

c+++++++++++++ configurations (in a file)
               write(2) (ss(i),i=1,nsubject) 

c+++++++++++++ predictive survival information

               do i=1,npred
                  do j=1,2*nvar

                     tmp2=sqrt(sigma(j,j))
                     ii=(i-1)*2*nvar+j
                     
                     do k=1,ngrid
                        fw(k)=0.d0
                        fw2(k)=0.d0
                     end do
                     
                     do k=1,maxm
                        if(prob(k).gt.0.d0)then
                           tmp1=0.d0
                           do l=1,p+q
                              if(model(j,l).eq.1)then
                                 tmp1=tmp1+xpred(i,possi(j,l))*
     &                                b(k,possi(j,l))
                              end if   
                           end do
                           
                           do l=1,ngrid
                              fw(l)=fw(l)+prob(k)*
     &                              cdflnorm(grid(j,l),tmp1,tmp2,1,0)

                              fw2(l)=fw2(l)+prob(k)*
     &                              dlnrm(grid(j,l),tmp1,tmp2,0)
                           end do
                        end if
                     end do

                     do k=1,ngrid
                        fw(k)=1.d0-fw(k)
                        fw2(k)=fw2(k)/fw(k)
                        f(ii,k)=f(ii,k)+fw(k)
                        h(ii,k)=h(ii,k)+fw2(k)
                     end do

                     call fmedian_sba(j,nvar,ngrid,fw,grid,tmp1)
                     randsave(isave,ii)=tmp1

                     write(3) (fw(k),k=1,ngrid)
                     write(4) (fw2(k),k=1,ngrid)

                  end do
               end do   
               
c+++++++++++++ kernel variance
               k=0
               do i=1,2*nvar
                  do j=i,2*nvar
                     k=k+1
                     thetasave(isave,k)=sigma(i,j)
                  end do
               end do
               
c+++++++++++++ baseline mean
               do i=1,nvar*(p+q)
                  k=k+1 
                  thetasave(isave,k)=mub(i) 
               end do

c+++++++++++++ baseline covariance
               do i=1,nvar*(p+q)
                  do j=i,nvar*(p+q)
                     k=k+1
                     thetasave(isave,k)=sigmab(i,j)
                  end do
               end do

c+++++++++++++ PY informartion
               k=k+1
               thetasave(isave,k)=ncluster

               k=k+1
               thetasave(isave,k)=alpha(1)

               k=k+1
               thetasave(isave,k)=alpha(2)

               if(alpha(1).eq.0.d0)then
                  acrate(2)=acrate(2)+1.d0
               end if

c+++++++++++++ parameters (in a file)
               write(5)(thetasave(isave,j),j=1,k)

c+++++++++++++ print
               skipcount = 0
               if(dispcount.ge.ndisplay)then
                  call cpu_time(sec1)
                  sec00=sec00+(sec1-sec0)
                  sec=sec00
                  sec0=sec1
                  tmp1=sprint(isave,nsave,sec)
                  dispcount=0
               end if   
            end if
         end if   

      end do

      do i=1,npred*2*nvar
         do j=1,ngrid
            f(i,j)=f(i,j)/dble(nsave)
            h(i,j)=h(i,j)/dble(nsave)
         end do
      end do   
      
      acrate(2)=acrate(2)/dble(nsave)
      acrate(1)=acrate(1)/dble(nscan)    

      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=5)
      
      return
      end


c=======================================================================
      subroutine spydens_sbapi(maxm,weights,qq,a,b,mub,sigmab,aa0,
     &                         ab0,logp)
c=======================================================================      
c
c     Alejandro Jara, 2007 
c=======================================================================            

      implicit none

c+++++Parameters
      integer maxm
      double precision weights(maxm)
      double precision mub,sigmab
      double precision qq,a,b,logp
      double precision aa0,ab0

c+++++Working internal
      integer i
      double precision tmp1,tmp2,tmp3
      
      double precision dbet,dnrm,cdfnorm

c+++++algorithm
     
      logp=0.d0 
     
      do i=1,(maxm-1)
         call rchkusr()

         tmp1=1.d0-a
         tmp2=b+dble(i)*a
         tmp3=dbet(weights(i),tmp1,tmp2,1)

         logp=logp+tmp3
      end do
      
      logp=logp+dnrm(b,mub,sqrt(sigmab),1)-
     &          cdfnorm(-a,mub,sqrt(sigmab),0,1)


      if(a.eq.0.d0)then
         logp=logp+qq
        else
         logp=logp+log(1.d0-qq)+dbet(a,aa0,ab0,1)
      end if 
      
      return
      end



c=======================================================================
      subroutine spylocnc_sba2(nsubject,nvar,p,q,
     &                         maxm,ccluster,cstrt,
     &                         varind,y,x,sigmainv,
     &                         mub,sigmab,sigmabinv,b,
     &                         iflagc,theta,
     &                         workvc,workmhc,
     &                         ztz,zty)
c=======================================================================      
c
c     Alejandro Jara, 2007 
c=======================================================================            
      implicit none
c+++++Parameters
      integer nsubject,nvar,p,q
      integer maxm
      integer ccluster(maxm)
      integer cstrt(maxm,nsubject)
      integer varind(nvar*(p+q))
      
      double precision y(nsubject,2*nvar)
      double precision x(nsubject,nvar*(p+q))
      
      double precision b(maxm,nvar*(p+q))
      double precision sigmainv(2*nvar,2*nvar)            
      double precision mub(nvar*(p+q))
      double precision sigmab(nvar*(p+q),nvar*(p+q))
      double precision sigmabinv(nvar*(p+q),nvar*(p+q))

c+++++Working external
      integer iflagc(nvar*(p+q))
      double precision theta(nvar*(p+q))
      double precision workvc(nvar*(p+q))
      double precision workmhc(nvar*(p+q)*(nvar*(p+q)+1)/2)
      double precision ztz(nvar*(p+q),nvar*(p+q))
      double precision zty(nvar*(p+q))

c+++++Working internal
      integer ii,j,jj,k,l
      integer ns
      integer pos1,pos2
      double precision tmp1

c+++++Algorithm

      do j=1,maxm

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         ns=ccluster(j)
         
         if(ns.gt.0)then

            do k=1,nvar*(p+q)
               tmp1=0.d0
               do l=1,nvar*(p+q)
                  ztz(k,l)=sigmabinv(k,l)
                  tmp1=tmp1+sigmabinv(k,l)*mub(l)
               end do
               zty(k)=tmp1
            end do

            do k=1,ns
               do ii=1,nvar*(p+q)
                  pos1=varind(ii)
                  do jj=1,nvar*(p+q)
                     pos2=varind(jj)
                     ztz(ii,jj)=ztz(ii,jj)+
     &                          x(cstrt(j,k),ii)*sigmainv(pos1,pos2)* 
     &                          x(cstrt(j,k),jj)
                  end do

                  do jj=1,2*nvar
                     zty(ii)=zty(ii)+
     &                        x(cstrt(j,k),ii)*sigmainv(pos1,jj)* 
     &                        y(cstrt(j,k),jj)
                  end do
               end do
            end do 

            call inverse(ztz,nvar*(p+q),iflagc) 
               
            do k=1,nvar*(p+q)
               tmp1=0.d0
               do l=1,nvar*(p+q)
                  tmp1=tmp1+ztz(k,l)*zty(l)                     
               end do
               workvc(k)=tmp1
            end do

            call rmvnorm(nvar*(p+q),workvc,ztz,workmhc,zty,theta)

            do k=1,nvar*(p+q)
               b(j,k)=theta(k)
            end do
          
          else

            call rmvnorm(nvar*(p+q),mub,sigmab,workmhc,workvc,theta)          

            do k=1,nvar*(p+q)
               b(j,k)=theta(k)
            end do

         end if
      end do

      return
      end


c=======================================================================
      subroutine fmedian_sba(ii,nvar,ngrid,f,grid,meds)
c=======================================================================
c     find the median in the survival curve
c
c     Alejandro Jara, 2007 
c=======================================================================            
      implicit none

c+++++Parameters 
      integer ii,nvar,ngrid
      double precision grid(2*nvar,ngrid)
      double precision f(ngrid)
      double precision meds

c+++++Working internal
      integer i
      double precision tmp1,tmp2

c+++++Algorithm

      if(f(ngrid).gt.0.5d0)then
         tmp1=f(ngrid-1)
         tmp2=(f(ngrid)-f(ngrid-1))/(grid(ii,ngrid)-grid(ii,ngrid-1))
         meds=grid(ii,ngrid-1)+(0.5d0-tmp1)/tmp2
        else
         i=1
         tmp1=f(1)
         do while(tmp1.gt.0.5d0.and.i.lt.ngrid)
            i=i+1     
            tmp1=f(i)
         end do
         
         if(tmp1.eq.0.5d0)then
            meds=grid(ii,i)
          else
            tmp1=f(i-1)
            tmp2=(f(i)-f(i-1))/(grid(ii,i)-grid(ii,i-1))
            meds=grid(ii,i-1)+(0.5d0-tmp1)/tmp2
         end if   
      end if
      return
      end


c=======================================================================
      subroutine spydens_sba(maxm,weights,a,b,mub,sigmab,loglik)
c=======================================================================      
c     evaluate the weights density in a steak-breaking approximation. 
c
c     Alejandro Jara, 2007 
c=======================================================================            

      implicit none

c+++++Parameters
      integer maxm
      double precision weights(maxm)
      double precision mub,sigmab
      double precision a,b,loglik

c+++++Working internal
      integer i
      double precision tmp1,tmp2,tmp3
      
      double precision dbet,dnrm,cdfnorm

c+++++algorithm
     
      loglik=0.d0 
     
      do i=1,(maxm-1)
         call rchkusr()

         tmp1=1.d0-a
         tmp2=b+dble(i)*a
         tmp3=dbet(weights(i),tmp1,tmp2,1)

         loglik=loglik+tmp3
      end do
      
      loglik=loglik+dnrm(b,mub,sqrt(sigmab),1)-
     &              cdfnorm(-a,mub,sqrt(sigmab),0,1)
      
      return
      end


c=======================================================================
      subroutine spydens_sba2(maxm,weights,a,b,loglik)
c=======================================================================      
c     evaluate the weights density in a steak-breaking approximation. 
c
c     Alejandro Jara, 2007 
c=======================================================================            

      implicit none

c+++++Parameters
      integer maxm
      double precision weights(maxm)
      double precision a,b,loglik

c+++++Working internal
      integer i
      double precision tmp1,tmp2,tmp3
      double precision dbet

c+++++algorithm
     
      loglik=0.d0 
     
      do i=1,(maxm-1)
         call rchkusr()

         tmp1=1.d0-a
         tmp2=b+dble(i)*a
         tmp3=dbet(weights(i),tmp1,tmp2,1)

         loglik=loglik+tmp3
      end do
      
      return
      end



c=======================================================================
      subroutine spyprob_sba(maxm,ccluster,prob,weights,a,b,lsweight)
c=======================================================================      
c     sampling weights using a steak-breaking approximation in a
c     dependent survival model for doubly-interval-censored data. 
c
c     Alejandro Jara, 2007 
c=======================================================================            

      implicit none

c+++++Parameters
      integer maxm
      integer ccluster(maxm)
      double precision prob(maxm)
      double precision weights(maxm)
      double precision a,b,lsweight

c+++++Working internal
      integer i,j
      double precision tmp1,tmp2,tmp3
      double precision lweight1
      
      double precision rbeta

c+++++algorithm
     
      lweight1=0.d0
      lsweight=0.d0

      do i=1,(maxm-1)
         call rchkusr()

         tmp1=1.d0-a+dble(ccluster(i))

         tmp2=b+dble(i)*a
         do j=i+1,maxm
            tmp2=tmp2+dble(ccluster(j)) 
         end do

         tmp3=rbeta(tmp1,tmp2)
         weights(i)=tmp3
            
         prob(i)=lweight1+log(tmp3)
         lweight1=lweight1+log(1.d0-tmp3)

         if((1.d0-weights(i)).lt.1.d-300)then
             lsweight=lsweight+0.d0
           else 
             lsweight=lsweight+log(1.d0-weights(i))
         end if   

      end do
      weights(maxm)=1.d0
      prob(maxm)=lweight1
      
      do i=1,maxm
         prob(i)=dexp(prob(i))
      end do
      
      return
      end
 

c=======================================================================
      subroutine spyconf_sba(maxm,nsubject,nvar,p,q,
     &                       ncluster,ccluster,cstrt,ss,
     &                       model,possi,y,x,b,sigma,
     &                       prob,iflagn,prob2,
     &                       workmn,workvn)
c=======================================================================      
c     Steak-breaking approximation for the DP prior is the
c     dependent survival model for doubly-interval-censored data. 
c
c     Alejandro Jara, 2007 
c=======================================================================            
      implicit none
c+++++Parameters
      integer maxm,nsubject,nvar,p,q
      integer ncluster
      integer ccluster(maxm)
      integer cstrt(maxm,nsubject)
      integer ss(nsubject)
      integer model(2*nvar,p+q)
      integer possi(2*nvar,p+q)
      
      double precision y(nsubject,2*nvar)
      double precision x(nsubject,nvar*(p+q))
      double precision b(maxm,nvar*(p+q))
      double precision sigma(2*nvar,2*nvar)      
      double precision prob(maxm)

c+++++Working external
      integer iflagn(2*nvar)
      double precision prob2(maxm)
      double precision workmn(2*nvar,2*nvar)      
      double precision workvn(2*nvar)

c+++++Working internal
      integer i,j,k,l,evali
      double precision detlog,tmp1
      double precision tpi,work1,work2,work3,sse

c+++++Algorithm

      do i=1,maxm
         ccluster(i)=0
      end do

      tpi=6.283185307179586476925286766559d0
      work1=-(dble(2*nvar)*log(tpi))

      do i=1,2*nvar
         do j=1,2*nvar
            workmn(i,j)=sigma(i,j) 
         end do
      end do

      call inversedet(workmn,2*nvar,iflagn,detlog)
      work2=detlog

      do i=1,nsubject

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do j=1,maxm 
            
            if(prob(j).gt.0.d0)then
               do k=1,2*nvar
                  tmp1=0.d0
                  do l=1,p+q
                     if(model(k,l).eq.1)then
                        tmp1=tmp1+x(i,possi(k,l))*b(j,possi(k,l))
                     end if   
                  end do
                  workvn(k)=y(i,k)-tmp1
               end do   
            
               sse=0.d0
               do k=1,2*nvar
                  do l=1,2*nvar
                     sse=sse+workvn(k)*workmn(k,l)*workvn(l)
                  end do
               end do
               work3=sse
            
               tmp1=(work1-work2-work3)/2.d0
               if(tmp1.gt.700.d0)tmp1=700.d0
               prob2(j)=prob(j)*exp(tmp1)
             else
               prob2(j)=0.d0
            end if   
         end do
         
         call simdisc(prob2,maxm,maxm,evali)         
         
         ss(i)=evali
         ccluster(evali)=ccluster(evali)+1
         cstrt(evali,ccluster(evali))=i
      end do   
      
      ncluster=0
      do i=1,maxm
         if(ccluster(i).gt.0)ncluster=ncluster+1
      end do
      
      return
      end
      

c=======================================================================
      subroutine spylocnc_sba(nsubject,nvar,p,q,
     &                        maxm,ccluster,cstrt,
     &                        model,possi,y,x,sigmainv,
     &                        mub,sigmab,sigmabinv,b,
     &                        iflagc,theta,
     &                        workvc,workmhc,
     &                        ztz,zty)
c=======================================================================      
c
c     Alejandro Jara, 2007 
c=======================================================================            
      implicit none
c+++++Parameters
      integer nsubject,nvar,p,q
      integer maxm
      integer ccluster(maxm)
      integer cstrt(maxm,nsubject)
      integer model(2*nvar,p+q)
      integer possi(2*nvar,p+q)
      
      double precision y(nsubject,2*nvar)
      double precision x(nsubject,nvar*(p+q))
      
      double precision b(maxm,nvar*(p+q))
      double precision sigmainv(2*nvar,2*nvar)            
      double precision mub(nvar*(p+q))
      double precision sigmab(nvar*(p+q),nvar*(p+q))
      double precision sigmabinv(nvar*(p+q),nvar*(p+q))

c+++++Working external
      integer iflagc(nvar*(p+q))
      double precision theta(nvar*(p+q))
      double precision workvc(nvar*(p+q))
      double precision workmhc(nvar*(p+q)*(nvar*(p+q)+1)/2)
      double precision ztz(nvar*(p+q),nvar*(p+q))
      double precision zty(nvar*(p+q))

c+++++Working internal
      integer ii,j,jj,k,l,m
      integer ns
      integer pos1,pos2
      double precision tmp1

c+++++Algorithm

      do j=1,maxm

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         ns=ccluster(j)
         
         if(ns.gt.0)then

            do k=1,nvar*(p+q)
               tmp1=0.d0
               do l=1,nvar*(p+q)
                  ztz(k,l)=sigmabinv(k,l)
                  tmp1=tmp1+sigmabinv(k,l)*mub(l)
               end do
               zty(k)=tmp1
            end do

            do k=1,ns
               do ii=1,2*nvar
                  do l=1,p+q
                     if(model(ii,l).eq.1)then
                        pos1=possi(ii,l)
                        do jj=1,2*nvar 
                           do m=1,p+q   
                              if(model(jj,m).eq.1)then
                                 pos2=possi(jj,m)
                                 ztz(pos1,pos2)=ztz(pos1,pos2)+
     &                             x(cstrt(j,k),pos1)*sigmainv(ii,jj)* 
     &                             x(cstrt(j,k),pos2)
                              end if
                           end do
                           zty(pos1)=zty(pos1)+
     &                       x(cstrt(j,k),pos1)*sigmainv(ii,jj)* 
     &                       y(cstrt(j,k),jj)
                        end do
                     end if   
                  end do
               end do
            end do 

            call inverse(ztz,nvar*(p+q),iflagc) 
               
            do k=1,nvar*(p+q)
               tmp1=0.d0
               do l=1,nvar*(p+q)
                  tmp1=tmp1+ztz(k,l)*zty(l)                     
               end do
               workvc(k)=tmp1
            end do

            call rmvnorm(nvar*(p+q),workvc,ztz,workmhc,zty,theta)

            do k=1,nvar*(p+q)
               b(j,k)=theta(k)
            end do
          
          else

            call rmvnorm(nvar*(p+q),mub,sigmab,workmhc,workvc,theta)          

            do k=1,nvar*(p+q)
               b(j,k)=theta(k)
            end do

         end if
      end do

      return
      end


c=======================================================================
      subroutine spykcnc(maxm,nsubject,nvar,p,q,ss,
     &                   model,possi,
     &                   y,x,b,nu01,tinv1,sigma,sigmainv,
     &                   iflagn,workvn,workmn,workmn2,workmn3,
     &                   workmhn,workmhn2)
c=======================================================================
c     Updating the Kernel covariance
c
c     Alejandro Jara, 2007
c=======================================================================

      implicit none
c+++++Parameters
      integer maxm   
      integer nsubject,nvar,p,q
      integer ss(nsubject)
      integer model(2*nvar,p+q)
      integer possi(2*nvar,p+q)
      integer nu01

      double precision y(nsubject,2*nvar)
      double precision x(nsubject,nvar*(p+q))
      double precision b(maxm,nvar*(p+q))
      double precision tinv1(2*nvar,2*nvar)
      double precision sigma(2*nvar,2*nvar)      
      double precision sigmainv(2*nvar,2*nvar)            

c+++++Working external
      integer iflagn(2*nvar)
      double precision workvn(2*nvar)      
      double precision workmn(2*nvar,2*nvar)      
      double precision workmn2(2*nvar,2*nvar)      
      double precision workmn3(2*nvar,2*nvar)      
      double precision workmhn(2*nvar*(2*nvar+1)/2)
      double precision workmhn2(2*nvar*(2*nvar+1)/2)

c+++++Working internal
      integer i,j,k,l
      double precision tmp1

c++++ check if the user has requested an interrupt
      call rchkusr()

      do i=1,2*nvar
         iflagn(i)=0 
         workvn(i)=0.d0 
         do j=1,2*nvar
            workmn(i,j)=0.d0
            workmn2(i,j)=0.d0
            workmn3(i,j)=0.d0
         end do
      end do
      
      do i=1,2*nvar*(2*nvar+1)/2
         workmhn(i)=0.d0
         workmhn2(i)=0.d0
      end do
      
      do i=1,nsubject

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do k=1,2*nvar
            tmp1=0.d0
            do l=1,p+q
               if(model(k,l).eq.1)then
                  tmp1=tmp1+x(i,possi(k,l))*b(ss(i),possi(k,l))
               end if   
            end do
            workvn(k)=y(i,k)-tmp1
         end do   
      
         do j=1,2*nvar
            do k=1,2*nvar
               workmn(j,k)=workmn(j,k)+workvn(j)*workvn(k)
            end do
         end do            
      end do
      
      do i=1,2*nvar
         do j=1,2*nvar
            workmn(i,j)=workmn(i,j)+tinv1(i,j)
         end do
      end do

      call riwishart(2*nvar,nu01+nsubject,workmn,workmn2,workmn3,
     &               workvn,workmhn,workmhn2,iflagn)

      do i=1,2*nvar
         do j=1,2*nvar
            sigma(i,j)=workmn(i,j)
            sigmainv(i,j)=workmn2(i,j)
         end do
      end do

      return
      end
      
      
c=======================================================================      
      subroutine hpdspy(nsave,npred,nvar,ngrid,alpha,tint,
     &                  workv1,workv2,llower,lupper,llower2,lupper2)
c=======================================================================
c     Compute CI for survival curves.
c
c     Alejandro Jara, 2007
c=======================================================================
      implicit none 
c+++++External parameters
      integer tint
      integer nsave,npred,nvar,ngrid 
      double precision alpha

c+++++External working
      double precision workv1(nsave),workv2(ngrid)

c+++++Output      
      double precision llower(npred*2*nvar,ngrid)
      double precision lupper(npred*2*nvar,ngrid)
      double precision llower2(npred*2*nvar,ngrid)
      double precision lupper2(npred*2*nvar,ngrid)

c+++++Internal parameters
      integer maxnsave,maxngrid
      parameter(maxnsave=30000,maxngrid=300)
      double precision aupp(2),alow(2)
      double precision workm(maxnsave,maxngrid)
      double precision workm2(maxnsave,maxngrid)

c+++++Internal working
      integer i,ii,j,jj,k,l   


c+++++algorithm

      if(maxnsave.lt.nsave)then
         call rexit("Increase 'maxnsave' in 'hpdspy'")
      end if   

      if(maxngrid.lt.ngrid)then
         call rexit("Increase 'maxngrid' in 'hpdspy'")
      end if   

      open(unit=1,file='out_dppackage3.txt',status='old',
     &     form='unformatted')

      open(unit=2,file='out_dppackage4.txt',status='old',
     &     form='unformatted')

      do ii=1,npred*2*nvar
         do i=1,nsave 
            call rchkusr()
            do j=1,npred
               do k=1,2*nvar
                  jj=(j-1)*2*nvar+k
                  read(1) (workv2(l),l=1,ngrid)
                  if(ii.eq.jj)then
                     do l=1,ngrid
                        workm(i,l)=workv2(l) 
                     end do
                  end if

                  read(2) (workv2(l),l=1,ngrid)
                  if(ii.eq.jj)then
                     do l=1,ngrid
                        workm2(i,l)=workv2(l) 
                     end do
                  end if

               end do   
            end do
          end do  
          rewind(unit=1)
          rewind(unit=2)
          
          do i=1,ngrid
             do j=1,nsave
                workv1(j)=workm(j,i) 
             end do
          
             call hpd(nsave,alpha,workv1,alow,aupp)
          
             if(tint.eq.1)then
c               (alow(1),aupp(1)): 100(1-alpha)% HPD interval
c 
                llower(ii,i)=alow(1)
                lupper(ii,i)=aupp(1)

              else
c              (alow(2),aupp(2)): 100(1-alpha)% Bayesian credible 
c              interval

                llower(ii,i)=alow(2)
                lupper(ii,i)=aupp(2)
             end if


             do j=1,nsave
                workv1(j)=workm2(j,i) 
             end do
          
             call hpd(nsave,alpha,workv1,alow,aupp)
          
             if(tint.eq.1)then
c               (alow(1),aupp(1)): 100(1-alpha)% HPD interval
c 
                llower2(ii,i)=alow(1)
                lupper2(ii,i)=aupp(1)

              else
c              (alow(2),aupp(2)): 100(1-alpha)% Bayesian credible 
c              interval

                llower2(ii,i)=alow(2)
                lupper2(ii,i)=aupp(2)
             end if

          end do

      end do      

      close(unit=1)      
      return
      end
