      !****************************************************
      ! riemann solver (Rusanov 24 august 2021 solver)
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
        double precision dcl, dcr, dlambda, dlambdal, dlambdar

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

        ! calculate lambda
        dlambda = dmax1(dabs(dlambdal),dabs(dlambdar))

        ! calculate fluxes
        frho = 0.5d0 * (dfrhol + dfrhor) + 
     =         0.5d0 * dlambda * (dcrhol - dcrhor)
        frvx = 0.5d0 * (dfrvxl + dfrvxr) + 
     =         0.5d0 * dlambda * (dcrvxl - dcrvxr)
        frvy = 0.5d0 * (dfrvyl + dfrvyr) + 
     =         0.5d0 * dlambda * (dcrvyl - dcrvyr)
        frvz = 0.5d0 * (dfrvzl + dfrvzr) + 
     =         0.5d0 * dlambda * (dcrvzl - dcrvzr)

      end subroutine riemann
      !**************************************************
