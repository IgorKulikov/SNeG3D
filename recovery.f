      !****************************************************
      ! recovery primitive variables
      !****************************************************
      subroutine primitive
        use mesh
        use omp_lib
        implicit none
        double precision, external :: eos

        ! variables 
        double precision densval, dimpx, dimpy, dimpz
        integer i, k, l, m

        ! recovery primitive variables
        do m=1,Ndepth
!$omp parallel do 
!$omp& private(i,k,l)
!$omp& private(densval,dimpx,dimpy,dimpz) 
!$omp& shared(m,density,dimplsx,dimplsy,dimplsz)
!$omp& shared(rho,velx,vely,velz,press)
!$omp& schedule(dynamic)
!$omp& num_threads(MIC_NUM_THREADS)
         do l=1,N+2
          do k=1,N+2
           do i=1,N+2
       
            ! get conservative variables
            densval = density(i,k,l,m)
            dimpx   = dimplsx(i,k,l,m)
            dimpy   = dimplsy(i,k,l,m)
            dimpz   = dimplsz(i,k,l,m)

            ! set physics variables
            rho(i,k,l,m).vl   = densval
            velx(i,k,l,m).vl  = dimpx/densval
            vely(i,k,l,m).vl  = dimpy/densval
            velz(i,k,l,m).vl  = dimpz/densval
            press(i,k,l,m).vl = eos(densval)

           enddo   
          enddo
         enddo
!$omp end parallel do 
        enddo

      end subroutine primitive
      !**************************************************
