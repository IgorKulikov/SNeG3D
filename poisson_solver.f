      !****************************************************
      ! poisson solver subroutine
      !****************************************************
      subroutine poisson_solver(m)
        use machine
        use mesh
        use poisson
        implicit none

        ! input mesh level
        integer, intent(in) :: m

        ! variables 
        integer i, k, l, NumIter
        double precision NormRight, NormR, Residual, hloc
        double precision alpha, alphach, alphazn, beta, betach, betazn

        ! conjugate gradient method
        hloc = h/(2.d0**(m-1))

        psn_f = gravity(:,:,:,m)                  ! right part formation
        forall(i = 2:N+1, k = 2:N+1, l = 2:N+1)
            psn_f(i,k,l) = 4.d0 * pi * density(i,k,l,m) * hloc * hloc
        end forall
        
        NormRight = dsqrt(sum(psn_f * psn_f))     ! |f|
        psn_x = gravity(:,:,:,m)                  ! x = x0

        psn_q = psn_x                             ! q = Ax

        forall(i = 2:N+1, k = 2:N+1, l = 2:N+1)
            psn_q(i,k,l) = -38.d0/9.d0 * psn_x(i,k,l) +
     =         4.d0/9.d0  * ( psn_x(i+1,k,l) + psn_x(i-1,k,l) + 
     =                        psn_x(i,k+1,l) + psn_x(i,k-1,l) +
     =                        psn_x(i,k,l+1) + psn_x(i,k,l-1) ) + 
     =         1.d0/9.d0  * ( psn_x(i,k+1,l+1) + psn_x(i,k-1,l+1) + 
     =                        psn_x(i,k+1,l-1) + psn_x(i,k-1,l-1) +
     =                        psn_x(i+1,k,l+1) + psn_x(i-1,k,l+1) + 
     =                        psn_x(i+1,k,l-1) + psn_x(i-1,k,l-1) +
     =                        psn_x(i+1,k+1,l) + psn_x(i-1,k+1,l) + 
     =                        psn_x(i+1,k-1,l) + psn_x(i-1,k-1,l) ) + 
     =         1.d0/36.d0 * ( psn_x(i+1,k+1,l+1) + psn_x(i-1,k+1,l+1) +
     =                        psn_x(i+1,k+1,l-1) + psn_x(i-1,k+1,l-1) +
     =                        psn_x(i+1,k-1,l+1) + psn_x(i-1,k-1,l+1) +
     =                        psn_x(i+1,k-1,l-1) + psn_x(i-1,k-1,l-1) )
        end forall

        psn_r = psn_f - psn_q                     ! r = f - q
        psn_q = psn_r                             ! q = r
        psn_z = psn_q                             ! z = q
	
        NormR = dsqrt(sum(psn_r * psn_r))         ! |r|
	Residual = NormR/NormRight                ! |r|/|f|

        alphach = sum(psn_q * psn_r)              ! alch = (q,r)

        NumIter = 0
        do while(Residual > EPS .and. NumIter < 1000)
          psn_q = 0.d0                            ! q = Az
                        
          forall(i = 2:N+1, k = 2:N+1, l = 2:N+1)
              psn_q(i,k,l) = -38.d0/9.d0 * psn_z(i,k,l) +
     =         4.d0/9.d0  * ( psn_z(i+1,k,l) + psn_z(i-1,k,l) + 
     =                        psn_z(i,k+1,l) + psn_z(i,k-1,l) +
     =                        psn_z(i,k,l+1) + psn_z(i,k,l-1) ) + 
     =         1.d0/9.d0  * ( psn_z(i,k+1,l+1) + psn_z(i,k-1,l+1) + 
     =                        psn_z(i,k+1,l-1) + psn_z(i,k-1,l-1) +
     =                        psn_z(i+1,k,l+1) + psn_z(i-1,k,l+1) + 
     =                        psn_z(i+1,k,l-1) + psn_z(i-1,k,l-1) +
     =                        psn_z(i+1,k+1,l) + psn_z(i-1,k+1,l) + 
     =                        psn_z(i+1,k-1,l) + psn_z(i-1,k-1,l) ) + 
     =         1.d0/36.d0 * ( psn_z(i+1,k+1,l+1) + psn_z(i-1,k+1,l+1) +
     =                        psn_z(i+1,k+1,l-1) + psn_z(i-1,k+1,l-1) +
     =                        psn_z(i+1,k-1,l+1) + psn_z(i-1,k-1,l+1) +
     =                        psn_z(i+1,k-1,l-1) + psn_z(i-1,k-1,l-1) )
          end forall

          alphazn = sum(psn_q * psn_z)            ! alzn = (q,z)
          alpha = alphach/alphazn                 ! al = alch/alzn
          psn_x = psn_x + alpha * psn_z           ! x = x + al * z
          psn_r = psn_r - alpha * psn_q           ! r = r - al * q
          psn_q = psn_r                           ! q = r
          betach = sum(psn_q * psn_r)             ! btch = (q,r)
          betazn = alphach                        ! btzn = alch
          alphach = betach                        ! alch = betach
          beta = betach/betazn                    ! bt = btch/btzn
          psn_z = psn_q + beta * psn_z            ! z = q + bt * z	
          NormR = dsqrt(sum(psn_r * psn_r))       ! |r|
         
          Residual = NormR/NormRight
          NumIter = NumIter + 1

        enddo

        !print *,NumIter,Residual

        gravity(:,:,:,m) = psn_x

      end subroutine poisson_solver
      !**************************************************

