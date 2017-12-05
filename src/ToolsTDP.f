c=======================================================================                  
c=======================================================================                  
c     SUBROUTINES FOR TDP Models
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
       subroutine jcomponenttd(y,k,j)
c=======================================================================
c      return the component j in {0,...,k} corresponding to the y
c      in [0,1] random variable in a Triangular-Dirichlet prior.
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none
       integer iceil,k,j
       double precision y

       j=iceil(y*dble(2*k))
       j=iceil(dble(j-1)/2.d0)
       
c       if(j.eq.0)then
c          if(y.gt.1.d0/dble(2*k))then
c             call dblepr("y",-1,y,1)
c             call intpr("k",-1,k,1)
c             call intpr("j",-1,j,1)          
c             call rexit("Error in 'jcomponentt'")      
c          end if
c        else if(j.eq.k)then
c          if(y.gt.1.d0)then
c             call dblepr("y",-1,y,1)
c             call intpr("k",-1,k,1)
c             call intpr("j",-1,j,1)          
c             call rexit("Error in 'jcomponentt'")      
c          end if
c          if(y.lt.dble(2*j-1)/dble(2*k))then
c             call dblepr("y",-1,y,1)
c             call intpr("k",-1,k,1)
c             call intpr("j",-1,j,1)                    
c             call rexit("Error in 'jcomponentt'")      
c          end if
c        else
c          if(y.lt.dble(2*j-1)/dble(2*k))then
c             call dblepr("y",-1,y,1)
c             call intpr("k",-1,k,1)
c             call intpr("j",-1,j,1)                    
c             call rexit("Error in 'jcomponentt'")      
c          end if
c          if(y.gt.dble(2*j+1)/dble(2*k))then
c             call dblepr("y",-1,y,1)
c             call intpr("k",-1,k,1)
c             call intpr("j",-1,j,1)                    
c             call rexit("Error in 'jcomponentt'")      
c          end if
c       end if  
      
       return
       end
       
c=======================================================================
       subroutine baseevaltd(x,k,a0,b0,eval)
c=======================================================================
c      evaluates the Mixture of Triangular Distributions at the 
c      baseline distribution, Beta(a0,b0), in a Triangular-Dirichlet 
c      prior.
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none
       integer i,k
       double precision a0,b0,cdfbetas,x,eval
       double precision a,b,c,tmp1

       eval=0.d0
       do i=0,k

          if(i.eq.0)then
             a=0.d0
             b=dble(i+1)/dble(k)
             c=0.d0

             call dtriang(x,a,b,c,tmp1)    

             eval=eval+
     &            cdfbetas(dble(i+1)/dble(2*k),a0,b0,1,0)*
     &            tmp1

           else if(i.eq.k)then  
             a=dble(i-1)/dble(k)
             b=1.d0
             c=1.d0
   
             call dtriang(x,a,b,c,tmp1)    

             eval=eval+
     &           (1.d0-
     &            cdfbetas(dble(2*i-1)/dble(2*k),a0,b0,1,0))*
     &            tmp1

           else
             a=dble(i-1)/dble(k)
             b=dble(i+1)/dble(k)
             c=dble(i)/dble(k)

             call dtriang(x,a,b,c,tmp1)    

             eval=eval+
     &           (cdfbetas(dble(2*i+1)/dble(2*k),a0,b0,1,0)-
     &            cdfbetas(dble(2*i-1)/dble(2*k),a0,b0,1,0))*
     &            tmp1
             
          end if  


       end do 
       return
       end
       
c=======================================================================
       subroutine clustevaltd(x,k,y,eval)
c=======================================================================
c      evaluates the cluster contribution for the cluster
c      defined by y in a Triangular-Dirichlet prior.
c
c      Alejandro Jara, 2007 
c=======================================================================
       implicit none
       integer j,k
       double precision x,y,eval
       double precision a,b,c,tmp1

       call jcomponenttd(y,k,j)       

       if(j.eq.0)then
          a=0.d0
          b=dble(j+1)/dble(k)
          c=0.d0
        else if(j.eq.k)then  
          a=dble(j-1)/dble(k)
          b=1.d0
          c=1.d0
        else
          a=dble(j-1)/dble(k)
          b=dble(j+1)/dble(k)
          c=dble(j)/dble(k)
       end if  

       call dtriang(x,a,b,c,tmp1)   
       eval=tmp1

       return
       end

c=======================================================================
       subroutine samplektd(nrec,x,y,prob,kmax,k)
c=======================================================================
c      generate k from its conditional distribution in a 
c      Triangular-Dirichlet prior. This function assumes that the prior 
c      for k is a uniform U(1,kmax).
c
c      Alejandro Jara, 2007.
c=======================================================================
       implicit none
       integer i,j,kmax,k,nrec
       double precision eval,prob(kmax+1),y(nrec),x(nrec),tmp1

       do i=1,kmax
          tmp1=0.d0         
          do j=1,nrec
             call clustevaltd(x(j),i,y(j),eval)
             if(eval.eq.0.d0)then
                eval=0.0001d0
             end if
             tmp1=tmp1+log(eval)
          end do
          
          if(tmp1.gt.700.d0)then
             tmp1=700.d0
          end if
          prob(i)=exp(tmp1-log(dble(kmax)))
       end do

       call simdisc(prob,kmax+1,kmax,k)

       return
       end


c=======================================================================
       subroutine sampleytd(x,kmax,prob,a0,b0,k,y)
c=======================================================================
c      generate y from the baseline in a 
c      Triangular-Dirichlet prior. 
c
c      Alejandro Jara, 2007.
c=======================================================================
       implicit none
       integer i,j,kmax,k,status
       double precision a0,b0,bound,prob(kmax+1),tmp1,tmp2,tmp3,x
       double precision y,y2
       double precision tt1,tt2,tt3,tt4
       double precision cdfbetas
       double precision a,b,c
       real runif

       do i=0,k
          if(i.eq.0)then
             a=0.d0
             b=dble(i+1)/dble(k)
             c=0.d0

             call dtriang(x,a,b,c,tmp2)    

             tmp1=
     &            cdfbetas(dble(i+1)/dble(2*k),a0,b0,1,0)*
     &            tmp2
             
           else if(i.eq.k)then  
             a=dble(i-1)/dble(k)
             b=1.d0
             c=1.d0

             call dtriang(x,a,b,c,tmp2)    

             tmp1=
     &           (cdfbetas(dble(2*i)/dble(2*k),a0,b0,1,0)-
     &            cdfbetas(dble(2*i-1)/dble(2*k),a0,b0,1,0))*
     &            tmp2
             
           else
             a=dble(i-1)/dble(k)
             b=dble(i+1)/dble(k)
             c=dble(i)/dble(k)

             call dtriang(x,a,b,c,tmp2)    

             tmp1=
     &           (cdfbetas(dble(2*i+1)/dble(2*k),a0,b0,1,0)-
     &            cdfbetas(dble(2*i-1)/dble(2*k),a0,b0,1,0))*
     &            tmp2
             
          end if  

          prob(i+1)=tmp1
       end do

       call simdisc(prob,kmax+1,k+1,j)
       j=j-1

       if(a0.eq.1.d0.and.b0.eq.1.d0)then
          if(j.eq.0)then
             y=dble(runif())/dble(2*k)
           else if(j.eq.k)then
             y=(dble(2*j-1)/dble(2*k))+dble(runif())/dble(2*k)
           else
             y=(dble(2*j-1)/dble(2*k))+dble(runif())/dble(k) 
          end if  

        else

          if(j.eq.0)then
             tt3=0.d0
           else if(j.eq.k)then
             tt3=dble(2*j-1)/dble(2*k)
           else
             tt3=dble(2*j-1)/dble(2*k)
          end if   
          tt4=1.d0-tt3
          call cdfbet(1,tt1,tt2,tt3,tt4,a0,b0,status,bound)
          if(status.ne.0)then
             call rexit("Error in 'sampleytd'")      
          end if
          tmp1=tt1

          if(j.eq.0)then
             tt3=1.d0/dble(2*k)
           else if(j.eq.k)then
             tt3=1.d0
           else
             tt3=dble(2*j+1)/dble(2*k)
          end if   
          tt4=1.d0-tt3
          call cdfbet(1,tt1,tt2,tt3,tt4,a0,b0,status,bound)
          if(status.ne.0)then
             call rexit("Error in 'sampleytd'")      
          end if
          tmp2=tt1
 
          tmp3=tmp1+dble(runif())*(tmp2-tmp1) 
       
          call cdfbet(2,tmp3,1.d0-tmp3,y,y2,a0,b0,status,bound)
          if(status.ne.0)then
             call rexit("Error in 'sampleytd'")      
          end if

       end if

       if(j.eq.0)then
          if(y.gt.1.d0/dble(2*k))then
             call dblepr("y",-1,y,1)
             call intpr("k",-1,k,1)
             call intpr("j",-1,j,1)          
             call rexit("Error in 'sampleytd'")      
          end if
        else if(j.eq.k)then
          if(y.gt.1.d0)then
             call dblepr("y",-1,y,1)
             call intpr("k",-1,k,1)
             call intpr("j",-1,j,1)          
             call rexit("Error in 'sampleytd'")      
          end if
          if(y.lt.dble(2*j-1)/dble(2*k))then
             call dblepr("y",-1,y,1)
             call intpr("k",-1,k,1)
             call intpr("j",-1,j,1)                    
             call rexit("Error in 'sampleytd'")      
          end if
        else
          if(y.lt.dble(2*j-1)/dble(2*k))then
             call dblepr("y",-1,y,1)
             call intpr("k",-1,k,1)
             call intpr("j",-1,j,1)                    
             call rexit("Error in 'sampleytd'")      
          end if
          if(y.gt.dble(2*j+1)/dble(2*k))then
             call dblepr("y",-1,y,1)
             call intpr("k",-1,k,1)
             call intpr("j",-1,j,1)                    
             call rexit("Error in 'sampleytd'")      
          end if
       end if  
      
       return
       end
     
