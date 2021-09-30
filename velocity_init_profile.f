      !**************************************************
      ! initial velocity profiles
      !**************************************************
      double precision function velocity_profile(i,k,l,d,m)
        use machine
        use mesh 
        implicit none 

        ! input index of cell and d = 0(x), 1(y), 2(z)
        integer, intent(in) :: i, k, l, d, m

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
        velocity_profile = 0.d0

      end function velocity_profile
      !**************************************************

