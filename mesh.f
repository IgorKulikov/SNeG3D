      !****************************************************
      ! module for mesh description
      !****************************************************
      module mesh
        implicit none 

        !**************************************************
        ! structure for physics variables
        !**************************************************
        structure /physics/
          double precision vl
          double precision dx, dy, dz
        end structure

        !**************************************************
        ! mesh constants & variables
        !**************************************************
        double precision, parameter :: length = 3.2d0    ! length of domain
        double precision, parameter :: prntau = 1.0d0    ! time step for print
        double precision center                          ! center of domain
        integer, parameter :: N = 64                     ! length of mesh
        integer, parameter :: Ndepth = 3                 ! depth of mesh
        integer, parameter :: Ntime = 10                 ! number of print 
        integer, parameter :: MIC_NUM_THREADS = 4        ! number of OpenMP threads
        double precision h                               ! length of cell
        double precision tau                             ! time step
        double precision, parameter :: cfl = 0.2         ! Courant number
        double precision, parameter :: dgamma = 1.d0     ! adiabatic index
        double precision, parameter :: temper = 0.1d0    ! temperature
        !**************************************************

        !**************************************************
        ! conservative and physics variables on mesh
        !**************************************************
        double precision density(N+2, N+2, N+2, Ndepth)  ! conservative density
        double precision dimplsx(N+2, N+2, N+2, Ndepth)  ! conservative momentum by x
        double precision dimplsy(N+2, N+2, N+2, Ndepth)  ! conservative momentum by y
        double precision dimplsz(N+2, N+2, N+2, Ndepth)  ! conservative momentum by z
        double precision gravity(N+2, N+2, N+2, Ndepth)  ! gravity
        record /physics/ rho(N+2, N+2, N+2, Ndepth)      ! physics density
        record /physics/ press(N+2, N+2, N+2, Ndepth)    ! physics pressure
        record /physics/ velx(N+2, N+2, N+2, Ndepth)     ! physics velocity by x
        record /physics/ vely(N+2, N+2, N+2, Ndepth)     ! physics velocity by y
        record /physics/ velz(N+2, N+2, N+2, Ndepth)     ! physics velocity by z
        !**************************************************

      end module mesh
