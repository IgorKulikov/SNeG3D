      !****************************************************
      ! Poisson solver
      !****************************************************
      program code3d
        ! using modules
        use machine
        use mesh
        implicit none

        double precision timer, tend
        integer i, m
        integer begin, end, rate

        !****************************************************
        ! Main code
        !****************************************************
        timer = 0.d0                   ! reset timer   
        call load                      ! load problem
        call save(timer)               ! save init density

        do i = 1,Ntime                 ! begin cycle by global time 
          tend = prntau * i            ! calculate local max time
          do while(timer < tend)       ! begin cycle by local time

            call system_clock(begin,rate)
            do m = Ndepth, 2, -1       ! projection from m to m-1 
             call projection(m)
            enddo
            call system_clock(end)
            !print *,"Projection ",real(end-begin)/real(rate)

            do m = 1, Ndepth           ! hydrodynamics boundary conditions
             call hydro_boundary(m)    ! internal boundary conditions     
            enddo

            call system_clock(begin,rate)
            call primitive             ! physics variables recovery
            call system_clock(end)
            !print *,"Primitive ",real(end-begin)/real(rate)

            call system_clock(begin,rate)
            do m = 1, Ndepth           ! piecewise-linear recosntruction 
             call reconstruction(m) 
            enddo
            call system_clock(end)
            !print *,"Reconstruction ",real(end-begin)/real(rate)

            call system_clock(begin,rate)
            call get_tau(timer,tend)   ! calculate tau
            timer = timer + tau        ! increment timer
            call system_clock(end)
            !print *,"CFL ",real(end-begin)/real(rate)

            call system_clock(begin,rate)
            call multigrid             ! multigrid poisson solver
            call system_clock(end)
            !print *,"Poisson solver ",real(end-begin)/real(rate)

            call system_clock(begin,rate)
            call godunov               ! godunov method
            call system_clock(end)
            !print *,"Godunov solver ",real(end-begin)/real(rate)
            
            print *,"Time = ",timer    ! info print of current timer
          enddo                        ! end cycle by local time
          call save(tend)              ! save density
        enddo                          ! end cycle by global time 

      end program code3d

