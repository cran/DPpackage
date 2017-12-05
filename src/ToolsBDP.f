c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR BDP Models
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
      integer function ifloor(x)
c=======================================================================      
      implicit none
      double precision x
      ifloor=x
      if(ifloor.gt.x) ifloor=ifloor-1
      end  
      
c=======================================================================      
      subroutine nrmu(nsteps,nsubject,betar,b,z,sigma2,kk,mu0,s0,
     &                mu,meancg,varcg)
c=======================================================================      
c     Newton-Raphson for the full conditional of the mean in an
c     extended BDP prior.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none
      integer i,j,l
      integer nsteps,nsubject,kk
      double precision betar,b(nsubject),z(nsubject),sigma2,mu,mu0,s0
      double precision meancg,varcg,tmp1,tt1,tt2,q1,q2,q3,q4,q5,q6
      double precision first,second
      double precision cdfnorm,dnrm
      
      tmp1=mu

      do i=1,nsteps
      
         q1=0.d0
         q2=0.d0
         q3=0.d0
         q4=0.d0
         q5=0.d0
         q6=0.d0
         
         do j=1,nsubject
            call jcomponentbd(z(j),kk,l)
            
            tt1=dnrm((b(j)-tmp1)/sqrt(sigma2),0.d0,1.d0,0)
            tt2=cdfnorm(b(j),tmp1,sqrt(sigma2),1,0)
            q1=q1+dble(l-1)*tt1/tt2 
            q2=q2+dble(kk-l)*tt1/(1.d0-tt2) 
            
            q3=q3-((b(j)-tmp1)/sqrt(sigma2))*dble(l-1)*tt1/tt2 
            q4=q4+dble(l-1)*(tt1*tt1)/(tt2*tt2)
            q5=q5-((b(j)-tmp1)/sqrt(sigma2))*dble(kk-l)*tt1/(1.d0-tt2) 
            q6=q6+dble(kk-l)*(tt1*tt1)/((1.d0-tt2)*(1.d0-tt2)) 
         end do
         
         first=betar-dble(nsubject)*tmp1-q1+q2
         first=first/sigma2
         first=first-(tmp1/s0)+(mu0/s0)

         second= -(dble(nsubject)/sigma2) +
     &            (1.d0/(sigma2*sqrt(sigma2)))*q3 -   
     &            (1.d0/(sigma2*sqrt(sigma2)))*q4 -   
     &            (1.d0/(sigma2*sqrt(sigma2)))*q5 -   
     &            (1.d0/(sigma2*sqrt(sigma2)))*q6 -   
     &            (1.d0/s0) 
     
         tmp1=tmp1-first/second 
      end do


      q3=0.d0
      q4=0.d0
      q5=0.d0
      q6=0.d0
      
      do j=1,nsubject
         call jcomponentbd(z(j),kk,l)

         tt1=dnrm((b(j)-tmp1)/sqrt(sigma2),0.d0,1.d0,0)
         tt2=cdfnorm(b(j),tmp1,sqrt(sigma2),1,0)

         q3=q3-((b(j)-tmp1)/sqrt(sigma2))*dble(l-1)*tt1/tt2 
         q4=q4+dble(l-1)*(tt1*tt1)/(tt2*tt2)
         q5=q5-((b(j)-tmp1)/sqrt(sigma2))*dble(kk-l)*tt1/(1.d0-tt2) 
         q6=q6+dble(kk-l)*(tt1*tt1)/((1.d0-tt2)*(1.d0-tt2)) 
      end do
      
      second= -(dble(nsubject)/sigma2) +
     &         (1.d0/(sigma2*sqrt(sigma2)))*q3 -   
     &         (1.d0/(sigma2*sqrt(sigma2)))*q4 -   
     &         (1.d0/(sigma2*sqrt(sigma2)))*q5 -   
     &         (1.d0/(sigma2*sqrt(sigma2)))*q6 -   
     &         (1.d0/s0) 

      
      meancg=tmp1
      varcg=-1.d0/second

      if(varcg.lt.0.d0)then
        call rexit("Error Newton-Raphson approximation for mu")
      end if
      
      return
      end


c=======================================================================
       subroutine jcomponentbd(y,k,j)
c=======================================================================
c      return the component j in {1,...,k} corresponding to the y
c      in [0,1] random variable in a Bernstein-Dirichlet prior.
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none
       integer iceil,k,j
       double precision y
c       j=int(y*dble(k))+1
       
       j=iceil(y*k)
       
c       if(y.gt.dble(j)/dble(k))then
c          call dblepr("y",-1,y,1)
c          call intpr("k",-1,k,1)
c          call intpr("j",-1,j,1)          
c          call rexit("Error in 'jcomponent'")      
c       end if

c       if(y.lt.dble(j-1)/dble(k))then
c          call dblepr("y",-1,y,1)
c          call intpr("k",-1,k,1)
c          call intpr("j",-1,j,1)                    
c          call rexit("Error in 'jcomponent'")      
c       end if
      
       return
       end
       
c=======================================================================
       subroutine baseevalbd(x,k,a0,b0,eval)
c=======================================================================
c      evaluates the Bernstein polinomial at the 
c      baseline distribution, Beta(a0,b0), in a Bernstein-Dirichlet 
c      prior.
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none
       integer i,k
       double precision a0,b0,cdfbetas,dbet,x,eval

       eval=0.d0
       do i=1,k
          eval=eval+
     &        (cdfbetas(dble(i)/dble(k),a0,b0,1,0)-
     &         cdfbetas(dble(i-1)/dble(k),a0,b0,1,0))*
     &         dbet(x,dble(i),dble(k-i+1),0) 
       end do 
       return
       end
       
c=======================================================================
       subroutine clustevalbd(x,k,y,eval)
c=======================================================================
c      evaluates the cluster contribution for the cluster
c      defined by y in a Bernstein-Dirichlet prior.
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none
       integer j,k
       double precision dbet,x,y,eval

       call jcomponentbd(y,k,j)       
       eval=dbet(x,dble(j),dble(k-j+1),0) 

       return
       end

c=======================================================================
       subroutine samplekbd(nrec,x,y,prob,kmax,k)
c=======================================================================
c      generate k from its conditional distribution in a 
c      Bernstein-Dirichlet prior. This function assumes that the prior 
c      for k is a uniform U(1,kmax).
c
c      Alejandro Jara, 2007.
c=======================================================================
       implicit none
       integer i,j,kmax,k,nrec
       double precision eval,prob(kmax),y(nrec),x(nrec),tmp1

       do i=1,kmax
          tmp1=0.d0         
          do j=1,nrec
             call clustevalbd(x(j),i,y(j),eval)
             tmp1=tmp1+log(eval)
          end do
          prob(i)=exp(tmp1)*(1.d0/dble(kmax))
       end do

       call simdisc(prob,kmax,kmax,k)

       return
       end


c=======================================================================
       subroutine sampleybd(x,kmax,prob,a0,b0,k,y)
c=======================================================================
c      generate y from the baseline in a 
c      Bernstein-Dirichlet prior. 
c
c      Alejandro Jara, 2007.
c=======================================================================
       implicit none
       integer i,j,kmax,k,status
       double precision a0,b0,bound,dbet,prob(kmax),tmp1,tmp2,tmp3,x
       double precision y,y2
       double precision tt1,tt2,tt3,tt4
       double precision cdfbetas
       real runif

       do i=1,k
          tmp1=       
     &        (cdfbetas(dble(i)/dble(k),a0,b0,1,0)-
     &         cdfbetas(dble(i-1)/dble(k),a0,b0,1,0))*
     &         dbet(x,dble(i),dble(k-i+1),0)   
          prob(i)=tmp1
       end do

       call simdisc(prob,kmax,k,j)

       if(a0.eq.1.d0.and.b0.eq.1.d0)then
          y=(dble(j-1)+dble(runif()))/dble(k)
        else
          tt3=dble(j-1)/dble(k)
          tt4=1.d0-dble(j-1)/dble(k)
          call cdfbet(1,tt1,tt2,tt3,tt4,a0,b0,status,bound)
          if(status.ne.0)then
             call rexit("Error in 'sampleybd'")      
          end if
          tmp1=tt1
       
          tt3=dble(j)/dble(k)
          tt4=1.d0-dble(j)/dble(k)
          call cdfbet(1,tt1,tt2,tt3,tt4,a0,b0,status,bound)
          if(status.ne.0)then
             call rexit("Error in 'sampleybd'")      
          end if
          tmp2=tt1
 
          tmp3=tmp1+dble(runif())*(tmp2-tmp1) 
       
          call cdfbet(2,tmp3,1.d0-tmp3,y,y2,a0,b0,status,bound)
          if(status.ne.0)then
             call rexit("Error in 'sampleybd'")      
          end if

       end if
        
       if(y.gt.dble(j)/dble(k))then
          call rexit("Error in 'sampleybd'")      
       end if

       if(y.le.dble(j-1)/dble(k))then
          call rexit("Error in 'sampleybd'")      
       end if

       if(y.eq.0.d0)then
          call rexit("Error in 'sampleybd' (0)")      
       end if


       return
       end


c=======================================================================           
      subroutine logpsigma(nsubject,kk,b,z,mu,sigma,loglikn)
c=======================================================================      
c     Compute the full conditional of sigma in an extended BDP prior
c     when a uniform prior on sigma is used.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none
      integer i
      integer nsubject,kk,l
      double precision b(nsubject),z(nsubject),mu,sigma,loglikn
      double precision dnrm,cdfnorm
      
      loglikn=0.d0
      do i=1,nsubject
         call jcomponentbd(z(i),kk,l)       
         loglikn=loglikn+dble(l-1)*cdfnorm(b(i),mu,sigma,1,1)  
         loglikn=loglikn+dble(kk-l)*cdfnorm(b(i),mu,sigma,0,1)  
         loglikn=loglikn+dnrm(b(i),mu,sigma,1)
      end do   

      return
      end

    
c=======================================================================           
      subroutine logpsigma2(nsubject,kk,b,z,mu,sigma2,tau1,tau2,
     &                      loglikn)
c=======================================================================      
c     Compute the full conditional of sigma2 in an extended BDP prior
c     in the case an inverse gamma prior on sigma2 is used.
c
c     Alejandro Jara, 2007
c=======================================================================      
      implicit none
      integer i
      integer nsubject,kk,l
      double precision b(nsubject),z(nsubject),mu,sigma2,loglikn
      double precision dnrm,cdfnorm,sigma,tau1,tau2
      
      sigma=sqrt(sigma2)
      loglikn=0.d0
      do i=1,nsubject
         call jcomponentbd(z(i),kk,l)       
         loglikn=loglikn+dble(l-1)*cdfnorm(b(i),mu,sigma,1,1)  
         loglikn=loglikn+dble(kk-l)*cdfnorm(b(i),mu,sigma,0,1)  
         loglikn=loglikn+dnrm(b(i),mu,sigma,1)
      end do   

      loglikn=loglikn-(0.5d0*tau1+1.d0)*log(sigma2)-0.5d0*tau2/sigma2

      return
      end
                 
