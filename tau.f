      !****************************************************
      ! get time step
      !****************************************************
      subroutine get_tau(timer,tend)
        use mesh
        use machine
        use omp_lib
        implicit none

        ! input timer and maximal time
        double precision, intent(in) :: timer, tend

        ! variables 
        double precision dmaxv, dmaxvgl, dsnd, densval
        double precision dpressval, dvelx, dvely, dvelz
        integer i, k, l, m

        ! calculate maximal velocity
        dmaxv   = 0.d0
        dmaxvgl = 0.d0
        do m=1,Ndepth
!$omp parallel do reduction(max : dmaxv)
!$omp& private(i,k,l)
!$omp& private(dpressval,dvelx,dvely,dvelz,dsnd, densval) 
!$omp& shared(m,dmaxvgl)
!$omp& shared(rho,velx,vely,velz,press)
!$omp& schedule(dynamic)
!$omp& num_threads(MIC_NUM_THREADS)
         do l=1,N+2
          do k=1,N+2
           do i=1,N+2
       
            ! get physics variables
            densval   = rho(i,k,l,m).vl
            dpressval = press(i,k,l,m).vl
            dvelx     = velx(i,k,l,m).vl
            dvely     = vely(i,k,l,m).vl
            dvelz     = velz(i,k,l,m).vl
            dsnd      = dsqrt(dgamma * dpressval / densval)

            ! get maximal velocity
            if(dabs(dvelx)+dsnd > dmaxv) dmaxv = dabs(dvelx)+dsnd 
            if(dabs(dvely)+dsnd > dmaxv) dmaxv = dabs(dvely)+dsnd 
            if(dabs(dvelz)+dsnd > dmaxv) dmaxv = dabs(dvelz)+dsnd 
 				
           enddo   
          enddo
         enddo
!$omp end parallel do 
         if(dmaxv > dmaxvgl) then
          dmaxvgl = dmaxv
         endif

        enddo

        ! calculate time step
        tau = cfl * (h/(2.d0**(Ndepth-1))) / dmaxvgl

        ! reset to maximal time
        if(timer+tau >= tend) tau = tend - timer

      end subroutine get_tau
      !**************************************************
