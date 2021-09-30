      !****************************************************
      ! module for poisson solver
      !****************************************************
      module poisson
        use machine
        use mesh
        implicit none 

        !**************************************************
        ! special variables for poisson solver
        !**************************************************
        double precision psn_x(N+2, N+2, N+2)
        double precision psn_q(N+2, N+2, N+2)
        double precision psn_r(N+2, N+2, N+2)
        double precision psn_f(N+2, N+2, N+2)
        double precision psn_z(N+2, N+2, N+2)
        !**************************************************

      end module poisson