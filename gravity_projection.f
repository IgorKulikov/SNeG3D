      !****************************************************
      ! gravity projection for Gauss theorem
      !****************************************************
      subroutine gravity_projection
        use mesh
        implicit none

        ! variables 
        integer i, k, l, m, ism, ksm, lsm
        logical is_empty, is_xp, is_xm, is_yp, is_ym, is_zp, is_zm
        double precision dgrad, hloc, hlocd2

        ! godunov method
        do m=1,Ndepth
         ! calculate actual h
         hloc = h/(2.d0**(m-1))
         hlocd2 = hloc/2.d0

         do l=2,N+1
          do k=2,N+1
           do i=2,N+1

            ! current cell is empty
            is_empty = (m < Ndepth) 
     =                 .and. (i > N/4+1)
     =                 .and. (i < 3*N/4+2) 
     =                 .and. (k > N/4+1) 
     =                 .and. (k < 3*N/4+2)
     =                 .and. (l > N/4+1) 
     =                 .and. (l < 3*N/4+2)

            ! current cell is not empty
            if( .not. is_empty  ) then

             ! internal boundary
             is_xp    = (m < Ndepth) 
     =                  .and. (i == N/4+1) 
     =                  .and. (k > N/4+1) 
     =                  .and. (k < 3*N/4+2)
     =                  .and. (l > N/4+1) 
     =                  .and. (l < 3*N/4+2)
             is_xm    = (m < Ndepth) 
     =                  .and. (i == 3*N/4+2) 
     =                  .and. (k > N/4+1) 
     =                  .and. (k < 3*N/4+2)
     =                  .and. (l > N/4+1) 
     =                  .and. (l < 3*N/4+2)
             is_yp    = (m < Ndepth) 
     =                  .and. (i > N/4+1) 
     =                  .and. (i < 3*N/4+2) 
     =                  .and. (k == N/4+1) 
     =                  .and. (l > N/4+1) 
     =                  .and. (l < 3*N/4+2)
             is_ym    = (m < Ndepth) 
     =                  .and. (i > N/4+1) 
     =                  .and. (i < 3*N/4+2) 
     =                  .and. (k == 3*N/4+2)
     =                  .and. (l > N/4+1) 
     =                  .and. (l < 3*N/4+2)
             is_zp    = (m < Ndepth) 
     =                  .and. (i > N/4+1) 
     =                  .and. (i < 3*N/4+2) 
     =                  .and. (k > N/4+1) 
     =                  .and. (k < 3*N/4+2)
     =                  .and. (l == N/4+1)
             is_zm    = (m < Ndepth) 
     =                  .and. (i > N/4+1) 
     =                  .and. (i < 3*N/4+2) 
     =                  .and. (k > N/4+1) 
     =                  .and. (k < 3*N/4+2)
     =                  .and. (l == 3*N/4+2)
 
             ! x plus flux
             if( is_xp ) then
               ism = 2
               ksm = 2*(k - N/4 - 1)
               lsm = 2*(l - N/4 - 1)

               dgrad = ( (gravity(ism,ksm,lsm,m+1)   - 
     =                     gravity(ism-1,ksm,lsm,m+1))  /hlocd2 + 
     =                   (gravity(ism,ksm+1,lsm,m+1) - 
     =                     gravity(ism-1,ksm+1,lsm,m+1))/hlocd2 + 
     =                   (gravity(ism,ksm,lsm+1,m+1) - 
     =                     gravity(ism-1,ksm,lsm+1,m+1))/hlocd2 + 
     =                   (gravity(ism,ksm+1,lsm+1,m+1) - 
     =                     gravity(ism-1,ksm+1,lsm+1,m+1))/hlocd2 )/4.d0

               gravity(i+1,k,l,m) = gravity(i,k,l,m) + dgrad * hloc

             endif

             ! x minus flux
             if( is_xm ) then
               ism = N + 1
               ksm = 2*(k - N/4 - 1)
               lsm = 2*(l - N/4 - 1)

               dgrad = ( (gravity(ism+1,ksm,lsm,m+1)   - 
     =                     gravity(ism,ksm,lsm,m+1))  /hlocd2 + 
     =                   (gravity(ism+1,ksm+1,lsm,m+1) - 
     =                     gravity(ism,ksm+1,lsm,m+1))/hlocd2 + 
     =                   (gravity(ism+1,ksm,lsm+1,m+1) - 
     =                     gravity(ism,ksm,lsm+1,m+1))/hlocd2 + 
     =                   (gravity(ism+1,ksm+1,lsm+1,m+1) - 
     =                     gravity(ism,ksm+1,lsm+1,m+1))/hlocd2 )/4.d0

               gravity(i-1,k,l,m) = gravity(i,k,l,m) - dgrad * hloc

             endif

             ! y plus flux
             if( is_yp ) then
               ism = 2*(i - N/4 - 1)
               ksm = 2
               lsm = 2*(l - N/4 - 1)

               dgrad = ( (gravity(ism,ksm,lsm,m+1)   - 
     =                     gravity(ism,ksm-1,lsm,m+1))  /hlocd2 + 
     =                   (gravity(ism+1,ksm,lsm,m+1) - 
     =                     gravity(ism+1,ksm-1,lsm,m+1))/hlocd2 + 
     =                   (gravity(ism,ksm,lsm+1,m+1) - 
     =                     gravity(ism,ksm-1,lsm+1,m+1))/hlocd2 + 
     =                   (gravity(ism+1,ksm,lsm+1,m+1) - 
     =                     gravity(ism+1,ksm-1,lsm+1,m+1))/hlocd2 )/4.d0

               gravity(i,k+1,l,m) = gravity(i,k,l,m) + dgrad * hloc

             endif

             ! y minus flux
             if( is_ym )  then
               ism = 2*(i - N/4 - 1)
               ksm = N + 1
               lsm = 2*(l - N/4 - 1)

               dgrad = ( (gravity(ism,ksm+1,lsm,m+1)   - 
     =                     gravity(ism,ksm,lsm,m+1))  /hlocd2 + 
     =                   (gravity(ism+1,ksm+1,lsm,m+1) - 
     =                     gravity(ism+1,ksm,lsm,m+1))/hlocd2 + 
     =                   (gravity(ism,ksm+1,lsm+1,m+1) - 
     =                     gravity(ism,ksm,lsm+1,m+1))/hlocd2 + 
     =                   (gravity(ism+1,ksm+1,lsm+1,m+1) - 
     =                     gravity(ism+1,ksm,lsm+1,m+1))/hlocd2 )/4.d0

               gravity(i,k-1,l,m) = gravity(i,k,l,m) - dgrad * hloc

             endif 

             ! z plus flux
             if( is_zp ) then
               ism = 2*(i - N/4 - 1)
               ksm = 2*(k - N/4 - 1)
               lsm = 2

               dgrad = ( (gravity(ism,ksm,lsm,m+1)   - 
     =                     gravity(ism,ksm,lsm-1,m+1))  /hlocd2 + 
     =                   (gravity(ism+1,ksm,lsm,m+1) - 
     =                     gravity(ism+1,ksm,lsm-1,m+1))/hlocd2 + 
     =                   (gravity(ism,ksm+1,lsm,m+1) - 
     =                     gravity(ism,ksm+1,lsm-1,m+1))/hlocd2 + 
     =                   (gravity(ism+1,ksm+1,lsm,m+1) - 
     =                     gravity(ism+1,ksm+1,lsm-1,m+1))/hlocd2 )/4.d0

               gravity(i,k,l+1,m) = gravity(i,k,l,m) + dgrad * hloc

             endif

             ! z minus flux
             if( is_zm ) then
               ism = 2*(i - N/4 - 1)
               ksm = 2*(k - N/4 - 1)
               lsm = N + 1

               dgrad = ( (gravity(ism,ksm,lsm+1,m+1)   - 
     =                     gravity(ism,ksm,lsm,m+1))  /hlocd2 + 
     =                   (gravity(ism+1,ksm,lsm+1,m+1) - 
     =                     gravity(ism+1,ksm,lsm,m+1))/hlocd2 + 
     =                   (gravity(ism,ksm+1,lsm+1,m+1) - 
     =                     gravity(ism,ksm+1,lsm,m+1))/hlocd2 + 
     =                   (gravity(ism+1,ksm+1,lsm+1,m+1) - 
     =                     gravity(ism+1,ksm+1,lsm,m+1))/hlocd2 )/4.d0

               gravity(i,k,l-1,m) = gravity(i,k,l,m) - dgrad * hloc
             endif

            endif

           enddo   
          enddo
         enddo

        enddo

      end subroutine gravity_projection
      !**************************************************
