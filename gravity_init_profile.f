      !**************************************************
      ! initial gravity profiles
      !**************************************************
      double precision function gravity_profile(i,k,l,m)
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

        ! return exact gravity
        if(radius <= 1.d0) then
          gravity_profile = 4.d0*PI*(radius**5.d0)/15.d0 - 
     =                      3.d0*PI*(radius**4.d0)/5.d0 + 
     =                      2.d0*PI*(radius**2.d0)/3.d0 - 3.d0*PI/5.d0
        else
          gravity_profile = -4.d0*PI/radius/15.d0;
        endif

      end function gravity_profile
      !**************************************************
