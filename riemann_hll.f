      !****************************************************
      ! riemann solver
      !****************************************************
      subroutine riemann(frho, frvx, frvy, frvz,
     =                   drhol, drhor, dpressl, dpressr,
     =                   dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)
        use mesh
        use machine
        implicit none

        ! output flux solution
        double precision, intent(out) :: frho, frvx, frvy, frvz

        ! input timer and maximal time
        double precision, intent(in) :: drhol, dpressl, dvxl, dvyl, dvzl
        double precision, intent(in) :: drhor, dpressr, dvxr, dvyr, dvzr

        ! variables 
        double precision dcrhol, dcrvxl, dcrvyl, dcrvzl
        double precision dcrhor, dcrvxr, dcrvyr, dcrvzr
        double precision dfrhol, dfrvxl, dfrvyl, dfrvzl
        double precision dfrhor, dfrvxr, dfrvyr, dfrvzr
        double precision dcl, dcr, dlambdal, dlambdar
        double precision dinvlambda

        ! left and right conservative vector formation
        dcrhol = drhol
        dcrvxl = drhol * dvxl
        dcrvyl = drhol * dvyl
        dcrvzl = drhol * dvzl

        dcrhor = drhor
        dcrvxr = drhor * dvxr
        dcrvyr = drhor * dvyr
        dcrvzr = drhor * dvzr

        ! left and right fluxes formation
        dfrhol = drhol * dvxl
        dfrvxl = drhol * dvxl * dvxl + dpressl
        dfrvyl = drhol * dvyl * dvxl
        dfrvzl = drhol * dvzl * dvxl

        dfrhor = drhor * dvxr
        dfrvxr = drhor * dvxr * dvxr + dpressr
        dfrvyr = drhor * dvyr * dvxr
        dfrvzr = drhor * dvzr * dvxr

        ! calculate sound velocities
        dcl = dsqrt(dgamma * dpressl / drhol)
        dcr = dsqrt(dgamma * dpressr / drhor)
        
        ! calculate lambda
        dlambdal = dmin1(dvxl - dcl, 0.d0)
        dlambdar = dmax1(dvxr + dcr, 0.d0)

        ! calculate fluxes
        if( dlambdal .eq. 0.d0 .and. dlambdar .eq. 0.d0 ) then
          frho = 0.5d0 * (dfrhol + dfrhor) 
          frvx = 0.5d0 * (dfrvxl + dfrvxr) 
          frvy = 0.5d0 * (dfrvyl + dfrvyr) 
          frvz = 0.5d0 * (dfrvzl + dfrvzr) 
        else
          dinvlambda = 1.d0/(dlambdar - dlambdal)
          frho = (dlambdar * dfrhol - dlambdal * dfrhor + 
     =      dlambdal * dlambdar * (dcrhor - dcrhol)) * dinvlambda
          frvx = (dlambdar * dfrvxl - dlambdal * dfrvxr + 
     =      dlambdal * dlambdar * (dcrvxr - dcrvxl)) * dinvlambda
          frvy = (dlambdar * dfrvyl - dlambdal * dfrvyr + 
     =      dlambdal * dlambdar * (dcrvyr - dcrvyl)) * dinvlambda
          frvz = (dlambdar * dfrvzl - dlambdal * dfrvzr + 
     =      dlambdal * dlambdar * (dcrvzr - dcrvzl)) * dinvlambda
        endif

      end subroutine riemann
      !**************************************************
 


