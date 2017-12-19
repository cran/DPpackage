c=======================================================================
      subroutine hdpmdensity (nrec,nvar,nstudy,
     &                        study,y,
     &                        ngrid,grid,
     &                        nu,nub,
     &                        a0b0,tinv,m0,prec0,tbinv,aebe,pieps,
     &                        sameps,sammu,samsb,
     &                        ncluster,ss,sc,alpha,muclus,sigma,mub,
     &                        sigmab,eps,
     &                        mcmc,nsave,seed,
     &                        densmc,densms,thetasave,randsave,
     &                        cstrt,ccluster,scstrt,sccluster,
     &                        iflagp,ywork,muwork,sigmainv,sigmabinv,
     &                        sigmawork,prob,quadf,workm1,
     &                        workmh1,workmh2,workv1)
c=======================================================================
c     # 47
c
c     Subroutine `hdpmdensity' to run a Markov chain in a  
c     hierarchical DPM of normals model for density
c     estimation. 
c
c     Copyright: Alejandro Jara, 2009-2010.
c
c     Version 2.0
c
c     Last modification: 14-05-2010.
c
c     Changes and Bug fixes: 
c
c     Version 1.0 to Version 2.0:
c          - Bug in computing mean vector for new cluster. Thanks to
c            Ana Paula Sales for reporting the bug.  
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
c++++ data
      integer nrec,nvar,nstudy
      integer study(nrec)
      double precision y(nrec,nvar)

c++++ prediction
      integer ngrid
      double precision grid(ngrid,nvar)

c++++ prior
      integer nu,nub,sameps,sammu,samsb
      double precision a0b0(nstudy+1,2)
      double precision tinv(nvar,nvar)
      double precision m0(nvar),prec0(nvar,nvar)
      double precision tbinv(nvar,nvar)
      double precision aebe(2)
      double precision pieps(2)

c++++ current value
      integer ncluster
      integer ss(nrec)
      integer sc(nrec)
      double precision alpha(nstudy+1)
      double precision muclus(nrec+100,nvar)
      double precision sigma(nvar,nvar)
      double precision mub(nvar)
      double precision sigmab(nvar,nvar)
      double precision eps

c++++ mcmc
      integer mcmc(3),nburn,nskip,nsave,ndisplay

c++++ output
      double precision densmc(nstudy+1,nvar*ngrid)
      double precision densms(nstudy,nvar*ngrid)
      double precision thetasave(nsave,nvar*(nvar+1)+1+
     &                 nstudy+1+nvar)
      double precision randsave(nsave,nstudy)
  
c++++ seeds
      integer seed1,seed2,seed(2)

c++++ external working space
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      integer scstrt(nstudy+1,nrec)
      integer sccluster(nstudy+1)
      integer iflagp(nvar)
      double precision ywork(nvar)
      double precision muwork(nvar)
      double precision sigmainv(nvar,nvar)
      double precision sigmabinv(nvar,nvar)
      double precision sigmawork(nvar,nvar)
      double precision prob(nrec+100)
      double precision quadf(nvar,nvar)
      double precision workm1(nvar,nvar)
      double precision workmh1(nvar*(nvar+1)/2)
      double precision workmh2(nvar*(nvar+1)/2)
      double precision workv1(nvar)

c++++ internal working space - General
      integer i,ii,j,jj,k,kk,l,ok
      integer dispcount
      integer evali
      integer iscan
      integer isave
      integer ns
      integer nscan
      integer since
      integer skipcount
      integer sprint
      double precision dnrm
      double precision tmp1,tmp2,tmp3,tmp4,tmp5
      double precision rbeta,lbetaf

c++++ internal working space - WDP
      integer studyind
      integer nj,n0,n1 
      integer ncj,nc0

c++++ CPU time
      double precision sec00,sec0,sec1,sec

c++++++++++++++++++++++++++
c     initialize variables
c++++++++++++++++++++++++++

c++++ mcmc, priors and "zipped"

      nburn=mcmc(1)
      nskip=mcmc(2)
      ndisplay=mcmc(3)
   
      seed1=seed(1)
      seed2=seed(2)

c++++ set random number generator

      call setall(seed1,seed2)
      
c++++ start the MCMC algorithm

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*nsave

      call cpu_time(sec0)
      sec00=0.0

c++++ cluster structure

      do i=1,nrec
         ccluster(ss(i))=ccluster(ss(i))+1
         cstrt(ss(i),ccluster(ss(i)))=i
      end do

      do i=1,ncluster
         sccluster(sc(i))=sccluster(sc(i))+1
         scstrt(sc(i),sccluster(sc(i)))=i
      end do

c++++ inverse of the kernel covariance
      do i=1,nvar
         do j=1,nvar
            sigmainv(i,j)=sigma(i,j)
            sigmabinv(i,j)=sigmab(i,j)
         end do
      end do
      call inverse(sigmainv,nvar,iflagp) 
      call inverse(sigmabinv,nvar,iflagp) 

c      call dblepr("sigmainv",-1,sigmainv,nvar*nvar)
c      call dblepr("sigmabinv",-1,sigmabinv,nvar*nvar)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ start the MCMC algorithm
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      isave=0
      skipcount=0
      dispcount=0
      nscan=nburn+(nskip+1)*(nsave)

      call cpu_time(sec0)
      sec00=0.d0
      
      do iscan=1,nscan

c+++++++ check if the user has requested an interrupt
         call rchkusr()

c++++++++++++++++++++++++++++++++++++++++
c+++++++ clustering structure         +++
c++++++++++++++++++++++++++++++++++++++++

         do i=1,nrec

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            do j=1,nvar
               ywork(j)=y(i,j)
            end do
            studyind=study(i)

            ns=ccluster(ss(i))
            
            if(ns.gt.1)then
               ccluster(ss(i))=ccluster(ss(i))-1 
               j=1
               ok=0
               do while(ok.eq.0.and.j.le.ns)
                  if(cstrt(ss(i),j).eq.i)ok=j
                  j=j+1
               end do
               if(ok.lt.ns)then
                  do j=ok,ns-1
                     cstrt(ss(i),j)=cstrt(ss(i),j+1)
                  end do
               end if
               cstrt(ss(i),ns)=0
             else

               since=ss(i)

               if(since.lt.ncluster)then
                  call relabelwdpm(i,since,nrec,nvar,nstudy,
     &                             ncluster,
     &                             ss,cstrt,ccluster,
     &                             sc,muclus)

                else
                  sc(ncluster)=0
                  do k=1,nvar
                     muclus(ncluster,k)=0.d0
                  end do
                  ccluster(ncluster)=0 
                  ncluster=ncluster-1
               end if   

               do k=1,nstudy+1
                  sccluster(k)=0
               end do     

               do k=1,ncluster
                  sccluster(sc(k))=sccluster(sc(k))+1
                  scstrt(sc(k),sccluster(sc(k)))=k
               end do

            end if

            nj=0
            ncj=sccluster(studyind)
            if(ncj.gt.0)then
               do j=1,ncj
                  ii=scstrt(studyind,j)
                  nj=nj+ccluster(ii)
               end do
            end if

            n0=0 
            nc0=sccluster(nstudy+1)
            if(nc0.gt.0)then
               do j=1,nc0
                  ii=scstrt(nstudy+1,j)
                  n0=n0+ccluster(ii)
               end do
            end if

            if(ncj.gt.0)then
               do j=1,ncj
                  ii=scstrt(studyind,j)
                  do k=1,nvar
                     muwork(k)=muclus(ii,k)
                     do l=1,nvar
                        sigmawork(k,l)=sigma(k,l)
                     end do 
                  end do                
                  call dmvnd(nvar,ywork,muwork,sigmawork,tmp2,iflagp)
                  prob(j)=exp(tmp2)*(1.d0-eps)*
     &                 dble(ccluster(ii))/(alpha(studyind)+dble(nj))
               end do
            end if

            if(nc0.gt.0)then
               do j=1,nc0
                  ii=scstrt(nstudy+1,j)
                  do k=1,nvar
                     muwork(k)=muclus(ii,k)
                     do l=1,nvar
                        sigmawork(k,l)=sigma(k,l)
                     end do 
                  end do                
                  call dmvnd(nvar,ywork,muwork,sigmawork,tmp2,iflagp)
                  prob(ncj+j)=exp(tmp2)*eps*
     &                dble(ccluster(ii))/(alpha(nstudy+1)+dble(n0))
               end do
            end if

            do k=1,nvar
               muwork(k)=mub(k)
               do l=1,nvar
                  sigmawork(k,l)=sigma(k,l)+sigmab(k,l)
               end do 
            end do                
            call dmvnd(nvar,ywork,muwork,sigmawork,tmp2,iflagp)
            prob(ncj+nc0+1)=exp(tmp2)*(1.d0-eps)*
     &                      alpha(studyind)/(alpha(studyind)+dble(nj))

            do k=1,nvar
               muwork(k)=mub(k)
               do l=1,nvar
                  sigmawork(k,l)=sigma(k,l)+sigmab(k,l)
               end do 
            end do                
            call dmvnd(nvar,ywork,muwork,sigmawork,tmp2,iflagp)
            prob(ncj+nc0+2)=exp(tmp2)*eps*
     &                      alpha(nstudy+1)/(alpha(nstudy+1)+dble(n0))

            call simdisc(prob,nrec+100,ncj+nc0+2,evali)

            if(evali.le.ncj)then
               ii=scstrt(studyind,evali)
               ss(i)=ii
               ccluster(ii)=ccluster(ii)+1
               cstrt(ii,ccluster(ii))=i

            else if(evali.gt.ncj.and.evali.le.(ncj+nc0))then
               ii=scstrt(nstudy+1,evali-ncj)
               ss(i)=ii
               ccluster(ii)=ccluster(ii)+1
               cstrt(ii,ccluster(ii))=i

            else
               ncluster=ncluster+1

               do k=1,nvar
                  tmp1=0.d0
                  tmp2=0.d0
                  do l=1,nvar
                     workm1(k,l)=sigmabinv(k,l)+sigmainv(k,l)
                     tmp1=tmp1+sigmainv(k,l)*ywork(l)
                     tmp2=tmp2+sigmabinv(k,l)*mub(l)
                  end do 
                  workv1(k)=tmp1+tmp2
               end do                
               call inverse(workm1,nvar,iflagp)
               do k=1,nvar
                  tmp1=0.d0
                  do l=1,nvar
                     tmp1=tmp1+workm1(k,l)*workv1(l)
                  end do
                  muwork(k)=tmp1
               end do   
  
               call rmvnorm(nvar,muwork,workm1,workmh1,workv1,ywork)
 
               do k=1,nvar
                  muclus(ncluster,k)=ywork(k)
               end do

               ss(i)=ncluster
               ccluster(ncluster)=1
               cstrt(ncluster,ccluster(ncluster))=i

               if(evali.eq.(ncj+nc0+1))then
                  sc(ncluster)=studyind
                  sccluster(studyind)=sccluster(studyind)+1
                  scstrt(studyind,sccluster(studyind))=ncluster
               else
                  sc(ncluster)=nstudy+1
                  sccluster(nstudy+1)=sccluster(nstudy+1)+1
                  scstrt(nstudy+1,sccluster(nstudy+1))=ncluster
               end if
            end if
      
         end do

c         call intpr("ncluster",-1,ncluster,1)
c         call intpr("ss",-1,ss,nrec)
c         call intpr("sc",-1,sc,ncluster)
c         call intpr("ccluster",-1,ccluster,ncluster)
c         call intpr("sccluster",-1,sccluster,nstudy+1)


c++++++++++++++++++++++++++++++++++++++++
c+++++++ resampling cluster locations +++
c++++++++++++++++++++++++++++++++++++++++

         do ii=1,nstudy+1

c++++++++++ check if the user has requested an interrupt
            call rchkusr()

            ns=sccluster(ii)

c            call intpr("DP",-1,ii,1)
c            call intpr("ns",-1,ns,1)

            if(ns.gt.0)then

               do jj=1,ns

c++++++++++++++++ check if the user has requested an interrupt
                  call rchkusr()

                  kk=scstrt(ii,jj)

c                  call intpr("cluster",-1,kk,1)

                  n0=ccluster(kk)

                  do i=1,nvar
                     tmp1=0.d0
                     do j=1,nvar
                        quadf(i,j)=sigmabinv(i,j)+
     &                             dble(n0)*sigmainv(i,j)
                        tmp1=tmp1+sigmabinv(i,j)*mub(j)
                     end do
                     workv1(i)=tmp1
                  end do 
                  call inverse(quadf,nvar,iflagp) 

                  do i=1,n0 
                     do j=1,nvar
                        tmp1=0.d0 
                        do k=1,nvar
                           tmp1=tmp1+sigmainv(j,k)*y(cstrt(kk,i),k)
                        end do
                        workv1(j)=workv1(j)+tmp1
                     end do
                  end do
                  
                  do i=1,nvar
                     tmp1=0.d0
                     do j=1,nvar
                        tmp1=tmp1+quadf(i,j)*workv1(j)
                     end do
                     ywork(i)=tmp1
                  end do

                  call rmvnorm(nvar,ywork,quadf,workmh1,workv1,muwork)

                  do i=1,nvar
                     muclus(kk,i)=muwork(i)
                  end do

c                  call dblepr("mu",-1,muwork,nvar)

               end do
            end if
         end do

 
c++++++++++++++++++++++++++++++++++++++++
c+++++++ kernel covariance matrix     +++
c++++++++++++++++++++++++++++++++++++++++
    
c+++++++ check if the user has requested an interrupt
         call rchkusr()

         if(nu.gt.0)then
            do i=1,nvar
               do j=1,nvar
                  quadf(i,j)=0.d0
               end do
            end do
         
            do i=1,nrec
               do j=1,nvar
                  do k=1,nvar
                     quadf(j,k)=quadf(j,k)+               
     &                   (y(i,j)-muclus(ss(i),j))*
     &                   (y(i,k)-muclus(ss(i),k))
                  end do
              end do            
           end do

           do i=1,nvar
              do j=1,nvar
                 quadf(i,j)=quadf(i,j)+tinv(i,j)
              end do
           end do

           call riwishart(nvar,nu+nrec,quadf,sigmawork,workm1,
     &                    workv1,
     &                    workmh1,workmh2,iflagp)
           do i=1,nvar
              do j=1,nvar
                 sigma(i,j)=quadf(i,j)
                 sigmainv(i,j)=sigmawork(i,j)
              end do
           end do

c           call dblepr("sigma",-1,sigma,nvar*nvar)
         end if

c++++++++++++++++++++++++++++++++++++++++
c+++++++ baseline mean                +++
c++++++++++++++++++++++++++++++++++++++++
        
         if(sammu.eq.1)then
            do i=1,nvar
               tmp1=0.d0
               do j=1,nvar
                  quadf(i,j)=prec0(i,j)+dble(ncluster)*sigmabinv(i,j)
                  tmp1=tmp1+prec0(i,j)*m0(j)
               end do
               workv1(i)=tmp1
            end do

            call inverse(quadf,nvar,iflagp) 

            do ii=1,ncluster

c+++++++++++++ check if the user has requested an interrupt
               call rchkusr()

               do i=1,nvar
                  tmp1=0.d0
                  do j=1,nvar
                     tmp1=tmp1+sigmabinv(i,j)*muclus(ii,j)
                  end do
                  workv1(i)=workv1(i)+tmp1  
               end do
            end do 

            do i=1,nvar
               tmp1=0.d0
               do j=1,nvar
                  tmp1=tmp1+quadf(i,j)*workv1(j)
               end do
               ywork(i)=tmp1
            end do

            call rmvnorm(nvar,ywork,quadf,workmh1,workv1,mub)

c            call dblepr("mub",-1,mub,nvar)
         end if


c++++++++++++++++++++++++++++++++++++++++
c+++++++ baseline covariance          +++
c++++++++++++++++++++++++++++++++++++++++

         if(samsb.eq.1)then
            do i=1,nvar
               do j=1,nvar
                  quadf(i,j)=0.d0
               end do
            end do
         
            do i=1,ncluster
               do j=1,nvar
                  do k=1,nvar
                     quadf(j,k)=quadf(j,k)+               
     &               (muclus(i,j)-mub(j))*
     &               (muclus(i,k)-mub(k))
                  end do
              end do            
            end do

            do i=1,nvar
               do j=1,nvar
                  quadf(i,j)=quadf(i,j)+tbinv(i,j)
               end do
            end do

            call riwishart(nvar,nub+ncluster,quadf,sigmawork,workm1,
     &                  workv1,
     &                  workmh1,workmh2,iflagp)
            do i=1,nvar
               do j=1,nvar
                  sigmab(i,j)=quadf(i,j)
                  sigmabinv(i,j)=sigmawork(i,j)
               end do
            end do

c            call dblepr("sigmab",-1,sigmab,nvar*nvar)
         end if
 
c++++++++++++++++++++++++++++++++++++++++
c+++++++ precision parameters         +++
c++++++++++++++++++++++++++++++++++++++++

c+++++++ check if the user has requested an interrupt
         call rchkusr()

         do i=1,nstudy+1
            if(a0b0(i,1).gt.0.0)then
               nj=0
               ncj=sccluster(i)
               do j=1,ncj
                  ii=scstrt(i,j)
                  nj=nj+ccluster(ii)
               end do
               if(nj.gt.0)then
                  call samalph(tmp1,a0b0(i,1),a0b0(i,2),ncj,nj)
                else
                  tmp1=rbeta(a0b0(i,1),a0b0(i,2))
               end if
               alpha(i)=tmp1
            end if
         end do
         
c         call dblepr("alpha",-1,alpha,nstudy+1)

c++++++++++++++++++++++++++++++++++++++++
c+++++++ sampling eps                 +++
c++++++++++++++++++++++++++++++++++++++++
         if(sameps.eq.1)then

            n0=0 
            nc0=sccluster(nstudy+1)
            if(nc0.gt.0)then
               do j=1,nc0
                  ii=scstrt(nstudy+1,j)
                  n0=n0+ccluster(ii)
               end do
            end if
            n1=nrec-n0

            tmp1=aebe(1)+dble(n0)
            tmp2=aebe(2)+dble(n1)
            tmp3=aebe(1)
            tmp4=aebe(2)
            tmp5=exp(lbetaf(tmp1,tmp2)-lbetaf(tmp3,tmp4))

            prob(1)=0.d0
            prob(2)=0.d0
            if(n0.eq.nrec)prob(1)=pieps(1)
            if(n1.eq.nrec)prob(2)=pieps(2)
            prob(3)=(1.d0-pieps(1)-pieps(2))*tmp5

            call simdisc(prob,nrec+100,3,evali)

            if(evali.eq.1)then
               eps=0.d0
              else if(evali.eq.2)then
               eps=1.d0
              else
               eps=rbeta(tmp1,tmp2)
            end if 

c            call dblepr("eps",-1,eps,1)
         end if
 
c++++++++++++++++++++++++++++++++++++++++
c+++++++ saving samples               +++
c++++++++++++++++++++++++++++++++++++++++

         if(iscan.gt.nburn)then
            skipcount=skipcount+1
            if(skipcount.gt.nskip)then
               isave=isave+1
               dispcount=dispcount+1

c+++++++++++++ kernel covariance

               ii=0
               do i=1,nvar
                  do j=i,nvar
                     ii=ii+1
                     thetasave(isave,ii)=sigma(i,j)
                  end do
               end do

c+++++++++++++ eps
               ii=ii+1
               thetasave(isave,ii)=eps

c+++++++++++++ precision parameters
               do i=1,nstudy+1
                  ii=ii+1
                  thetasave(isave,ii)=alpha(i)
               end do 

c+++++++++++++ baseline means

               do i=1,nvar
                  ii=ii+1 
                  thetasave(isave,ii)=mub(i)
               end do

c+++++++++++++ baseline covariance

               do i=1,nvar
                  do j=i,nvar
                     ii=ii+1
                     thetasave(isave,ii)=sigmab(i,j)
                  end do
               end do
               
c+++++++++++++ density of each component in the mixture

               do i=1,nstudy+1
                  nj=0
                  ncj=sccluster(i)

                  if(ncj.gt.0)then
                     do j=1,ncj
                        ii=scstrt(i,j)
                        nj=nj+ccluster(ii)
                     end do

                     do j=1,ncj
                        ii=scstrt(i,j)
                        prob(j)=dble(ccluster(ii))/(alpha(i)+dble(nj))
                     end do
                  end if
                  prob(ncj+1)=alpha(i)/(alpha(i)+dble(nj)) 

                  if(ncj.gt.0)then
                     do j=1,ncj
                        ii=scstrt(i,j)
                        kk=0
                        do k=1,nvar
                           tmp1=muclus(ii,k)
                           tmp2=sigma(k,k)
                           do l=1,ngrid
                              kk=kk+1
                            densmc(i,kk)=densmc(i,kk)+
     &                      prob(j)*dnrm(grid(l,k),tmp1,sqrt(tmp2),0)   
                           end do 
                        end do 
                     end do
                  end if

                  kk=0
                  do k=1,nvar
                     tmp1=mub(k)
                     tmp2=sigma(k,k)+sigmab(k,k)
                     do l=1,ngrid
                        kk=kk+1
                        densmc(i,kk)=densmc(i,kk)+
     &                  prob(ncj+1)*dnrm(grid(l,k),tmp1,sqrt(tmp2),0)   
                     end do 
                  end do 
               end do

c+++++++++++++ density for each study

               do i=1,nstudy

                  nj=0
                  ncj=sccluster(i)
                  if(ncj.gt.0)then
                     do j=1,ncj
                        ii=scstrt(i,j)
                        nj=nj+ccluster(ii)
                     end do
                  end if

                  n0=0 
                  nc0=sccluster(nstudy+1)
                  if(nc0.gt.0)then
                     do j=1,nc0
                        ii=scstrt(nstudy+1,j)
                        n0=n0+ccluster(ii)
                     end do
                  end if

                  tmp3=0.d0
                  if(ncj.gt.0)then
                     do j=1,ncj
                        ii=scstrt(i,j)
                        prob(j)=(1.d0-eps)*
     &                     dble(ccluster(ii))/(alpha(i)+dble(nj))
                        tmp3=tmp3+prob(j)
                     end do
                  end if
                 
                  if(nc0.gt.0)then
                     do j=1,nc0
                        ii=scstrt(nstudy+1,j)
                        prob(ncj+j)=eps*
     &                    dble(ccluster(ii))/(alpha(nstudy+1)+dble(n0))
                        tmp3=tmp3+prob(ncj+j)
                     end do
                  end if

                  prob(ncj+nc0+1)=(1.d0-eps)*alpha(i)/
     &                            (alpha(i)+dble(nj))
                  tmp3=tmp3+prob(ncj+nc0+1)
 
                  prob(ncj+nc0+2)=eps*alpha(nstudy+1)/
     &                            (alpha(nstudy+1)+dble(n0)) 
                  tmp3=tmp3+prob(ncj+nc0+2)

                  if(ncj.gt.0)then
                     do j=1,ncj
                        ii=scstrt(i,j)
                        kk=0
                        do k=1,nvar
                           tmp1=muclus(ii,k)
                           tmp2=sigma(k,k)
                           do l=1,ngrid
                              kk=kk+1
                              densms(i,kk)=densms(i,kk)+
     &                        prob(j)*dnrm(grid(l,k),tmp1,sqrt(tmp2),0)/
     &                        tmp3   
                           end do 
                        end do 
                     end do
                  end if
 
                  if(nc0.gt.0)then
                     do j=1,nc0
                        ii=scstrt(nstudy+1,j)
                        kk=0
                        do k=1,nvar
                           tmp1=muclus(ii,k)
                           tmp2=sigma(k,k)
                           do l=1,ngrid
                              kk=kk+1
                              densms(i,kk)=densms(i,kk)+
     &                        prob(ncj+j)*
     &                        dnrm(grid(l,k),tmp1,sqrt(tmp2),0)/
     &                        tmp3   
                           end do 
                        end do 
                     end do
                  end if

                  kk=0
                  do k=1,nvar
                     tmp1=mub(k)
                     tmp2=sigma(k,k)+sigmab(k,k)
                     do l=1,ngrid
                        kk=kk+1
                        densms(i,kk)=densms(i,kk)+
     &                  prob(ncj+nc0+1)*
     &                  dnrm(grid(l,k),tmp1,sqrt(tmp2),0)/tmp3   
                     end do 
                  end do 

                  kk=0
                  do k=1,nvar
                     tmp1=mub(k)
                     tmp2=sigma(k,k)+sigmab(k,k)
                     do l=1,ngrid
                        kk=kk+1
                        densms(i,kk)=densms(i,kk)+
     &                  prob(ncj+nc0+2)*
     &                  dnrm(grid(l,k),tmp1,sqrt(tmp2),0)/tmp3   
                     end do 
                  end do 

               end do
 
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


      do i=1,nstudy
         do j=1,ngrid*nvar
            densms(i,j)=densms(i,j)/dble(nsave)
         end do
      end do

      do i=1,nstudy+1
         do j=1,ngrid*nvar
            densmc(i,j)=densmc(i,j)/dble(nsave)
         end do
      end do

      return
      end


c=======================================================================
      subroutine relabelwdpm(ind,since,nrec,nvar,nstudy,
     &                       ncluster,
     &                       ss,cstrt,ccluster,
     &                       sc,muclus)
c=======================================================================
      implicit none
      integer ind,since,nrec,nvar,nstudy,ncluster
      integer ss(nrec)
      integer cstrt(nrec,nrec)
      integer ccluster(nrec)
      integer sc(nrec)
      double precision muclus(nrec+100,nvar)

      integer i,j,ns,ii

      do i=since+1,ncluster

         ns=ccluster(i)    
         
         do j=1,ns
c++++++++++ check if the user has requested an interrupt
            call rchkusr()
            ii=cstrt(i,j) 
            ss(ii)=i-1
         end do

         do j=1,ns
            cstrt(i-1,j)=cstrt(i,j) 
         end do

         do j=1,nvar
            muclus(i-1,j)=muclus(i,j)
         end do
         ccluster(i-1)=ccluster(i)

         sc(i-1)=sc(i)
      end do

      ss(ind)=0

      sc(ncluster)=0
      do i=1,nvar
         muclus(ncluster,i)=0.d0
      end do
      ccluster(ncluster)=0
      ncluster=ncluster-1      

      return
      end
      

