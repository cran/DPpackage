
c=======================================================================      
       subroutine hpd(n,alpha,x,alow,aupp)
c=======================================================================       
c      computing 100(1-alpha)% HPD and credible intervals for x 
c      use Chen-Shao HPD Estimation Algorithm 
c      (see page 219 of Monte Carlo Methods in Bayesian Computation, 
c      Springer-Verlag, 2000) 
c      ming-hui chen
c      july 23, 2001 at wpi
c      input:
c            alpha: confidence level,  0 < alpha < 1
c            n = mcmc sample size
c            x(n): a univariate vector of mcmc sample 
c      output: 
c            (alow(1),aupp(1)): 100(1-alpha)% HPD interval
c            (alow(2),aupp(2)): 100(1-alpha)% Bayesian credible interval
c
       implicit double precision (a-h,o-z)
       double precision x(n)
       double precision aupp(2),alow(2)

       whpd=0.d0
       aupp1=0.d0
       alow1=0.d0

       q1=(alpha/2.0d0)*float(n)
       q2=(1.0d0-alpha/2.0d0)*float(n)
       nq1=nint(q1)
       nq2=nint(q2)
       nq=nq2-nq1
       do 100 i=1,n-1
          do 110 j=i+1,n
             if (x(i) .gt. x(j)) then
                temp=x(i)
                x(i)=x(j)
                x(j)=temp
             end if
 110   continue
 100   continue
       do 120 j=1,n-nq
              pdiff1=x(j)
              pdiff2=x(j+nq)
              wb=pdiff2-pdiff1
              if (j .eq. 1) then
                 whpd=wb
                 aupp1=pdiff2
                 alow1=pdiff1
               else
                 if (whpd .gt. wb) then
                    whpd=wb
                    aupp1=pdiff2
                    alow1=pdiff1
                 end if
              end if
 120   continue
       alow(1)=alow1       
       aupp(1)=aupp1       
       alow(2)=x(nq1)
       aupp(2)=x(nq2)

       return
       end
      
 
