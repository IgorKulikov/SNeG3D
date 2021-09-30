      !****************************************************
      ! limiter function
      !****************************************************
      double precision function dlimiter(dminus, dzero, dplus)
        use machine
        use mesh 
        implicit none 

        ! input index of cell
        double precision, intent(in) :: dminus, dzero, dplus

        ! calculate limiter
        dlimiter = dmin1( dmax1(dzero-dminus,0.d0), 
     =                    dmax1(dplus-dzero,0.d0) ) + 
     =                    dmax1(dmin1(dzero-dminus,0.d0), 
     =                    dmin1(dplus-dzero,0.d0) )

      end function dlimiter

      !****************************************************
      ! piecewise linear reconstruction subroutine
      !****************************************************
      subroutine reconstruction(m)
        use mesh
        implicit none
        double precision, external :: dlimiter

        ! input mesh level
        integer, intent(in) :: m

        ! variables 
        integer i, k, l, ib, kb, lb

        ! reset linear reconstruction
        do l=1,N+2
         do k=1,N+2
           i = 1        ! left boundary reconstruction

           if(m == 1) then    ! zero reconstruction for basis mesh
             rho(i,k,l,m).dx   = 0.d0
             rho(i,k,l,m).dy   = 0.d0
             rho(i,k,l,m).dz   = 0.d0
             press(i,k,l,m).dx = 0.d0 
             press(i,k,l,m).dy = 0.d0 
             press(i,k,l,m).dz = 0.d0
             velx(i,k,l,m).dx  = 0.d0 
             velx(i,k,l,m).dy  = 0.d0 
             velx(i,k,l,m).dz  = 0.d0
             vely(i,k,l,m).dx  = 0.d0 
             vely(i,k,l,m).dy  = 0.d0 
             vely(i,k,l,m).dz  = 0.d0
             velz(i,k,l,m).dx  = 0.d0 
             velz(i,k,l,m).dy  = 0.d0 
             velz(i,k,l,m).dz  = 0.d0
           else               ! internal reconstruction
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             rho(i,k,l,m).dx   = rho(ib,kb,lb,m-1).dx   
             rho(i,k,l,m).dy   = rho(ib,kb,lb,m-1).dy   
             rho(i,k,l,m).dz   = rho(ib,kb,lb,m-1).dz   
             press(i,k,l,m).dx = press(ib,kb,lb,m-1).dx  
             press(i,k,l,m).dy = press(ib,kb,lb,m-1).dy  
             press(i,k,l,m).dz = press(ib,kb,lb,m-1).dz 
             velx(i,k,l,m).dx  = velx(ib,kb,lb,m-1).dx   
             velx(i,k,l,m).dy  = velx(ib,kb,lb,m-1).dy   
             velx(i,k,l,m).dz  = velx(ib,kb,lb,m-1).dz  
             vely(i,k,l,m).dx  = vely(ib,kb,lb,m-1).dx   
             vely(i,k,l,m).dy  = vely(ib,kb,lb,m-1).dy   
             vely(i,k,l,m).dz  = vely(ib,kb,lb,m-1).dz  
             velz(i,k,l,m).dx  = velz(ib,kb,lb,m-1).dx   
             velz(i,k,l,m).dy  = velz(ib,kb,lb,m-1).dy   
             velz(i,k,l,m).dz  = velz(ib,kb,lb,m-1).dz  
           endif

           i = N+2      ! right boundary condition

           if(m == 1) then    ! zero reconstruction for basis mesh
             rho(i,k,l,m).dx   = 0.d0
             rho(i,k,l,m).dy   = 0.d0
             rho(i,k,l,m).dz   = 0.d0
             press(i,k,l,m).dx = 0.d0 
             press(i,k,l,m).dy = 0.d0 
             press(i,k,l,m).dz = 0.d0
             velx(i,k,l,m).dx  = 0.d0 
             velx(i,k,l,m).dy  = 0.d0 
             velx(i,k,l,m).dz  = 0.d0
             vely(i,k,l,m).dx  = 0.d0 
             vely(i,k,l,m).dy  = 0.d0 
             vely(i,k,l,m).dz  = 0.d0
             velz(i,k,l,m).dx  = 0.d0 
             velz(i,k,l,m).dy  = 0.d0 
             velz(i,k,l,m).dz  = 0.d0
           else               ! internal reconstruction
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             rho(i,k,l,m).dx   = rho(ib,kb,lb,m-1).dx   
             rho(i,k,l,m).dy   = rho(ib,kb,lb,m-1).dy   
             rho(i,k,l,m).dz   = rho(ib,kb,lb,m-1).dz   
             press(i,k,l,m).dx = press(ib,kb,lb,m-1).dx  
             press(i,k,l,m).dy = press(ib,kb,lb,m-1).dy  
             press(i,k,l,m).dz = press(ib,kb,lb,m-1).dz 
             velx(i,k,l,m).dx  = velx(ib,kb,lb,m-1).dx   
             velx(i,k,l,m).dy  = velx(ib,kb,lb,m-1).dy   
             velx(i,k,l,m).dz  = velx(ib,kb,lb,m-1).dz  
             vely(i,k,l,m).dx  = vely(ib,kb,lb,m-1).dx   
             vely(i,k,l,m).dy  = vely(ib,kb,lb,m-1).dy   
             vely(i,k,l,m).dz  = vely(ib,kb,lb,m-1).dz  
             velz(i,k,l,m).dx  = velz(ib,kb,lb,m-1).dx   
             velz(i,k,l,m).dy  = velz(ib,kb,lb,m-1).dy   
             velz(i,k,l,m).dz  = velz(ib,kb,lb,m-1).dz  
           endif

         enddo
        enddo

        do l=1,N+2
         do i=1,N+2
           k = 1        ! top boundary reconstruction

           if(m == 1) then    ! zero reconstruction for basis mesh
             rho(i,k,l,m).dx   = 0.d0
             rho(i,k,l,m).dy   = 0.d0
             rho(i,k,l,m).dz   = 0.d0
             press(i,k,l,m).dx = 0.d0 
             press(i,k,l,m).dy = 0.d0 
             press(i,k,l,m).dz = 0.d0
             velx(i,k,l,m).dx  = 0.d0 
             velx(i,k,l,m).dy  = 0.d0 
             velx(i,k,l,m).dz  = 0.d0
             vely(i,k,l,m).dx  = 0.d0 
             vely(i,k,l,m).dy  = 0.d0 
             vely(i,k,l,m).dz  = 0.d0
             velz(i,k,l,m).dx  = 0.d0 
             velz(i,k,l,m).dy  = 0.d0 
             velz(i,k,l,m).dz  = 0.d0
           else               ! internal reconstruction
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             rho(i,k,l,m).dx   = rho(ib,kb,lb,m-1).dx   
             rho(i,k,l,m).dy   = rho(ib,kb,lb,m-1).dy   
             rho(i,k,l,m).dz   = rho(ib,kb,lb,m-1).dz   
             press(i,k,l,m).dx = press(ib,kb,lb,m-1).dx  
             press(i,k,l,m).dy = press(ib,kb,lb,m-1).dy  
             press(i,k,l,m).dz = press(ib,kb,lb,m-1).dz 
             velx(i,k,l,m).dx  = velx(ib,kb,lb,m-1).dx   
             velx(i,k,l,m).dy  = velx(ib,kb,lb,m-1).dy   
             velx(i,k,l,m).dz  = velx(ib,kb,lb,m-1).dz  
             vely(i,k,l,m).dx  = vely(ib,kb,lb,m-1).dx   
             vely(i,k,l,m).dy  = vely(ib,kb,lb,m-1).dy   
             vely(i,k,l,m).dz  = vely(ib,kb,lb,m-1).dz  
             velz(i,k,l,m).dx  = velz(ib,kb,lb,m-1).dx   
             velz(i,k,l,m).dy  = velz(ib,kb,lb,m-1).dy   
             velz(i,k,l,m).dz  = velz(ib,kb,lb,m-1).dz  
           endif

           k = N+2      ! bottom boundary reconstruction

           if(m == 1) then    ! zero reconstruction for basis mesh
             rho(i,k,l,m).dx   = 0.d0
             rho(i,k,l,m).dy   = 0.d0
             rho(i,k,l,m).dz   = 0.d0
             press(i,k,l,m).dx = 0.d0 
             press(i,k,l,m).dy = 0.d0 
             press(i,k,l,m).dz = 0.d0
             velx(i,k,l,m).dx  = 0.d0 
             velx(i,k,l,m).dy  = 0.d0 
             velx(i,k,l,m).dz  = 0.d0
             vely(i,k,l,m).dx  = 0.d0 
             vely(i,k,l,m).dy  = 0.d0 
             vely(i,k,l,m).dz  = 0.d0
             velz(i,k,l,m).dx  = 0.d0 
             velz(i,k,l,m).dy  = 0.d0 
             velz(i,k,l,m).dz  = 0.d0
           else               ! internal reconstruction
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             rho(i,k,l,m).dx   = rho(ib,kb,lb,m-1).dx   
             rho(i,k,l,m).dy   = rho(ib,kb,lb,m-1).dy   
             rho(i,k,l,m).dz   = rho(ib,kb,lb,m-1).dz   
             press(i,k,l,m).dx = press(ib,kb,lb,m-1).dx  
             press(i,k,l,m).dy = press(ib,kb,lb,m-1).dy  
             press(i,k,l,m).dz = press(ib,kb,lb,m-1).dz 
             velx(i,k,l,m).dx  = velx(ib,kb,lb,m-1).dx   
             velx(i,k,l,m).dy  = velx(ib,kb,lb,m-1).dy   
             velx(i,k,l,m).dz  = velx(ib,kb,lb,m-1).dz  
             vely(i,k,l,m).dx  = vely(ib,kb,lb,m-1).dx   
             vely(i,k,l,m).dy  = vely(ib,kb,lb,m-1).dy   
             vely(i,k,l,m).dz  = vely(ib,kb,lb,m-1).dz  
             velz(i,k,l,m).dx  = velz(ib,kb,lb,m-1).dx   
             velz(i,k,l,m).dy  = velz(ib,kb,lb,m-1).dy   
             velz(i,k,l,m).dz  = velz(ib,kb,lb,m-1).dz  
           endif

         enddo
        enddo

        do k=1,N+2
         do i=1,N+2
           l = 1        ! down boundary reconstruction

           if(m == 1) then    ! zero reconstruction for basis mesh
             rho(i,k,l,m).dx   = 0.d0
             rho(i,k,l,m).dy   = 0.d0
             rho(i,k,l,m).dz   = 0.d0
             press(i,k,l,m).dx = 0.d0 
             press(i,k,l,m).dy = 0.d0 
             press(i,k,l,m).dz = 0.d0
             velx(i,k,l,m).dx  = 0.d0 
             velx(i,k,l,m).dy  = 0.d0 
             velx(i,k,l,m).dz  = 0.d0
             vely(i,k,l,m).dx  = 0.d0 
             vely(i,k,l,m).dy  = 0.d0 
             vely(i,k,l,m).dz  = 0.d0
             velz(i,k,l,m).dx  = 0.d0 
             velz(i,k,l,m).dy  = 0.d0 
             velz(i,k,l,m).dz  = 0.d0
           else               ! internal reconstruction
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             rho(i,k,l,m).dx   = rho(ib,kb,lb,m-1).dx   
             rho(i,k,l,m).dy   = rho(ib,kb,lb,m-1).dy   
             rho(i,k,l,m).dz   = rho(ib,kb,lb,m-1).dz   
             press(i,k,l,m).dx = press(ib,kb,lb,m-1).dx  
             press(i,k,l,m).dy = press(ib,kb,lb,m-1).dy  
             press(i,k,l,m).dz = press(ib,kb,lb,m-1).dz 
             velx(i,k,l,m).dx  = velx(ib,kb,lb,m-1).dx   
             velx(i,k,l,m).dy  = velx(ib,kb,lb,m-1).dy   
             velx(i,k,l,m).dz  = velx(ib,kb,lb,m-1).dz  
             vely(i,k,l,m).dx  = vely(ib,kb,lb,m-1).dx   
             vely(i,k,l,m).dy  = vely(ib,kb,lb,m-1).dy   
             vely(i,k,l,m).dz  = vely(ib,kb,lb,m-1).dz  
             velz(i,k,l,m).dx  = velz(ib,kb,lb,m-1).dx   
             velz(i,k,l,m).dy  = velz(ib,kb,lb,m-1).dy   
             velz(i,k,l,m).dz  = velz(ib,kb,lb,m-1).dz  
           endif

           l = N+2      ! up boundary reconstruction

           if(m == 1) then    ! zero reconstruction for basis mesh
             rho(i,k,l,m).dx   = 0.d0
             rho(i,k,l,m).dy   = 0.d0
             rho(i,k,l,m).dz   = 0.d0
             press(i,k,l,m).dx = 0.d0 
             press(i,k,l,m).dy = 0.d0 
             press(i,k,l,m).dz = 0.d0
             velx(i,k,l,m).dx  = 0.d0 
             velx(i,k,l,m).dy  = 0.d0 
             velx(i,k,l,m).dz  = 0.d0
             vely(i,k,l,m).dx  = 0.d0 
             vely(i,k,l,m).dy  = 0.d0 
             vely(i,k,l,m).dz  = 0.d0
             velz(i,k,l,m).dx  = 0.d0 
             velz(i,k,l,m).dy  = 0.d0 
             velz(i,k,l,m).dz  = 0.d0
           else               ! internal reconstruction
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             rho(i,k,l,m).dx   = rho(ib,kb,lb,m-1).dx   
             rho(i,k,l,m).dy   = rho(ib,kb,lb,m-1).dy   
             rho(i,k,l,m).dz   = rho(ib,kb,lb,m-1).dz   
             press(i,k,l,m).dx = press(ib,kb,lb,m-1).dx  
             press(i,k,l,m).dy = press(ib,kb,lb,m-1).dy  
             press(i,k,l,m).dz = press(ib,kb,lb,m-1).dz 
             velx(i,k,l,m).dx  = velx(ib,kb,lb,m-1).dx   
             velx(i,k,l,m).dy  = velx(ib,kb,lb,m-1).dy   
             velx(i,k,l,m).dz  = velx(ib,kb,lb,m-1).dz  
             vely(i,k,l,m).dx  = vely(ib,kb,lb,m-1).dx   
             vely(i,k,l,m).dy  = vely(ib,kb,lb,m-1).dy   
             vely(i,k,l,m).dz  = vely(ib,kb,lb,m-1).dz  
             velz(i,k,l,m).dx  = velz(ib,kb,lb,m-1).dx   
             velz(i,k,l,m).dy  = velz(ib,kb,lb,m-1).dy   
             velz(i,k,l,m).dz  = velz(ib,kb,lb,m-1).dz  
           endif

         enddo
        enddo

        ! reconstruction
!$omp parallel do 
!$omp& private(i,k,l)
!$omp& shared(m,rho,velx,vely,velz,press)
!$omp& schedule(dynamic)
!$omp& num_threads(MIC_NUM_THREADS)
        do l=2,N+1
         do k=2,N+1
          do i=2,N+1

            ! density
            rho(i,k,l,m).dx = dlimiter(rho(i-1,k,l,m).vl, 
     =                                 rho(i,k,l,m).vl, 
     =                                 rho(i+1,k,l,m).vl)
            rho(i,k,l,m).dy = dlimiter(rho(i,k-1,l,m).vl, 
     =                                 rho(i,k,l,m).vl, 
     =                                 rho(i,k+1,l,m).vl)
            rho(i,k,l,m).dz = dlimiter(rho(i,k,l-1,m).vl, 
     =                                 rho(i,k,l,m).vl, 
     =                                 rho(i,k,l+1,m).vl)

            ! pressure 
            press(i,k,l,m).dx = dlimiter(press(i-1,k,l,m).vl, 
     =                                   press(i,k,l,m).vl, 
     =                                   press(i+1,k,l,m).vl)
            press(i,k,l,m).dy = dlimiter(press(i,k-1,l,m).vl, 
     =                                   press(i,k,l,m).vl, 
     =                                   press(i,k+1,l,m).vl)
            press(i,k,l,m).dz = dlimiter(press(i,k,l-1,m).vl, 
     =                                   press(i,k,l,m).vl, 
     =                                   press(i,k,l+1,m).vl)
 
            ! velocity
            velx(i,k,l,m).dx = dlimiter(velx(i-1,k,l,m).vl, 
     =                                  velx(i,k,l,m).vl, 
     =                                  velx(i+1,k,l,m).vl)
            velx(i,k,l,m).dy = dlimiter(velx(i,k-1,l,m).vl, 
     =                                  velx(i,k,l,m).vl, 
     =                                  velx(i,k+1,l,m).vl)
            velx(i,k,l,m).dz = dlimiter(velx(i,k,l-1,m).vl, 
     =                                  velx(i,k,l,m).vl, 
     =                                  velx(i,k,l+1,m).vl)

            vely(i,k,l,m).dx = dlimiter(vely(i-1,k,l,m).vl, 
     =                                  vely(i,k,l,m).vl, 
     =                                  vely(i+1,k,l,m).vl)
            vely(i,k,l,m).dy = dlimiter(vely(i,k-1,l,m).vl, 
     =                                  vely(i,k,l,m).vl, 
     =                                  vely(i,k+1,l,m).vl)
            vely(i,k,l,m).dz = dlimiter(vely(i,k,l-1,m).vl, 
     =                                  vely(i,k,l,m).vl, 
     =                                  vely(i,k,l+1,m).vl)

            velz(i,k,l,m).dx = dlimiter(velz(i-1,k,l,m).vl, 
     =                                  velz(i,k,l,m).vl, 
     =                                  velz(i+1,k,l,m).vl)
            velz(i,k,l,m).dy = dlimiter(velz(i,k-1,l,m).vl, 
     =                                  velz(i,k,l,m).vl, 
     =                                  velz(i,k+1,l,m).vl)
            velz(i,k,l,m).dz = dlimiter(velz(i,k,l-1,m).vl, 
     =                                  velz(i,k,l,m).vl, 
     =                                  velz(i,k,l+1,m).vl)

          enddo   
         enddo
        enddo
!$omp end parallel do

      end subroutine reconstruction
      !**************************************************

