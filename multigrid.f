      !****************************************************
      ! multigrid poisson solver
      !****************************************************
      subroutine multigrid
        use mesh
        implicit none

        ! variables 
        integer m

        do m = 1, Ndepth
          call poisson_boundary(m)  ! gravity boundary conditions
          call poisson_solver(m)    ! poisson solver
        enddo

        call gravity_projection     ! theorem Gauss for gravity force

      end subroutine multigrid
      !**************************************************



