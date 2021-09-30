      !**************************************************
      ! equation of state
      !**************************************************
      double precision function eos(dens)
        use machine
        use mesh 
        implicit none 

        ! input density
        double precision, intent(in) :: dens 

        eos = temper * (dens ** dgamma)

      end function eos
      !**************************************************

