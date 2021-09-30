      !****************************************************
      ! projection m to m-1 meshes
      !****************************************************
      subroutine projection(m)
        use mesh
        use omp_lib
        implicit none 

        ! input mesh level
        integer, intent(in) :: m

        ! variables 
        double precision dmdens, dmimpx, dmimpy, dmimpz
        integer i, k, l, ib, kb, lb

!$omp parallel do 
!$omp& private(i,k,l,ib,kb,lb)
!$omp& private(dmdens,dmimpx,dmimpy,dmimpz) 
!$omp& shared(m,density,dimplsx,dimplsy,dimplsz)
!$omp& schedule(dynamic)
!$omp& num_threads(MIC_NUM_THREADS)
        do l=2,N+1,2
         do k=2,N+1,2
          do i=2,N+1,2

            ! calculate median density and momentum
            dmdens = ( density(i+1,k+1,l+1,m) + 
     =                 density(i+1,k+1,l  ,m) + 
     =                 density(i+1,k  ,l+1,m) + 
     =                 density(i+1,k  ,l  ,m) +
     =                 density(i  ,k+1,l+1,m) + 
     =                 density(i  ,k+1,l  ,m) +
     =                 density(i  ,k  ,l+1,m) + 
     =                 density(i  ,k  ,l  ,m) ) / 8.d0
            dmimpx = ( dimplsx(i+1,k+1,l+1,m) + 
     =                 dimplsx(i+1,k+1,l  ,m) + 
     =                 dimplsx(i+1,k  ,l+1,m) + 
     =                 dimplsx(i+1,k  ,l  ,m) +
     =                 dimplsx(i  ,k+1,l+1,m) + 
     =                 dimplsx(i  ,k+1,l  ,m) +
     =                 dimplsx(i  ,k  ,l+1,m) + 
     =                 dimplsx(i  ,k  ,l  ,m) ) / 8.d0
            dmimpy = ( dimplsy(i+1,k+1,l+1,m) + 
     =                 dimplsy(i+1,k+1,l  ,m) + 
     =                 dimplsy(i+1,k  ,l+1,m) + 
     =                 dimplsy(i+1,k  ,l  ,m) +
     =                 dimplsy(i  ,k+1,l+1,m) + 
     =                 dimplsy(i  ,k+1,l  ,m) +
     =                 dimplsy(i  ,k  ,l+1,m) + 
     =                 dimplsy(i  ,k  ,l  ,m) ) / 8.d0
            dmimpz = ( dimplsz(i+1,k+1,l+1,m) + 
     =                 dimplsz(i+1,k+1,l  ,m) + 
     =                 dimplsz(i+1,k  ,l+1,m) + 
     =                 dimplsz(i+1,k  ,l  ,m) +
     =                 dimplsz(i  ,k+1,l+1,m) + 
     =                 dimplsz(i  ,k+1,l  ,m) +
     =                 dimplsz(i  ,k  ,l+1,m) + 
     =                 dimplsz(i  ,k  ,l  ,m) ) / 8.d0

            ib = N/4 + i/2 + 1
            kb = N/4 + k/2 + 1
            lb = N/4 + l/2 + 1

            density(ib,kb,lb,m-1) = dmdens
            dimplsx(ib,kb,lb,m-1) = dmimpx
            dimplsy(ib,kb,lb,m-1) = dmimpy
            dimplsz(ib,kb,lb,m-1) = dmimpz

          enddo   
         enddo
        enddo
!$omp end parallel do

      end subroutine projection
      !**************************************************
