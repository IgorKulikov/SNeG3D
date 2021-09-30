      !**************************************************
      ! initial density profiles
      !**************************************************
      double precision function density_profile(i,k,l,m)
        use machine
        use mesh 
        implicit none 

        ! input index of cell
        integer, intent(in) :: i, k, l, m 

        ! variables 
        double precision xcord, ycord, zcord, radius
        double precision hloc, centerloc

        ! compute actual h and center for local mesh
        hloc      = h     /(2.d0**(m-1))
        centerloc = length/(2.d0**m)

        ! compute coordinates
        xcord = (i-1)*hloc - centerloc - hloc/2.d0
        ycord = (k-1)*hloc - centerloc - hloc/2.d0
        zcord = (l-1)*hloc - centerloc - hloc/2.d0
        radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)

        ! return exact density
        if(radius <= 1.d0) then
          density_profile = 1.d0
!          density_profile = 2.d0*(radius**3.d0) - 
!     =                      3.d0*(radius**2.d0) + 1.d0 + eps
        else
          density_profile = eps
        endif

      end function density_profile
      !**************************************************

