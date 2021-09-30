      !****************************************************
      ! riemann solver
      !****************************************************
      subroutine riemann(frho, frvx, frvy, frvz,
     =                   drhol, drhor, dpressl, dpressr,
     =                   dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)
        use mesh
        use machine
        implicit none
        double precision, external :: eos

        ! output flux solution
        double precision, intent(out) :: frho, frvx, frvy, frvz

        ! input timer and maximal time
        double precision, intent(in) :: drhol, dpressl, dvxl, dvyl, dvzl
        double precision, intent(in) :: drhor, dpressr, dvxr, dvyr, dvzr
  
        ! variables
        double precision dmedrho, dmedu, dmedp, dmedc
        double precision dw1, dw2, dw3, dw4, dlamupc, dlamumc, dlamu
        double precision dbrho, dbux, dbuy, dbuz, dbpress 

        ! Roe average
        dmedrho = 0.5d0 * (drhol + drhor)
        dmedu   = 0.5d0 * (dvxl  + dvxr)
        dmedp   = 0.5d0 * (dpressl + dpressr)
        dmedc   = dsqrt(dgamma * dmedp / dmedrho)
        
        ! Eigenvalues
        dlamupc = dmedu + dmedc
        dlamumc = dmedu - dmedc
        dlamu   = dmedu

        ! Characteristic
        if(dlamupc > 0.d0) then
         dw1 = dvxl + drhol * dmedc/dmedrho
        else
         dw1 = dvxr + drhor * dmedc/dmedrho
        endif

        if(dlamumc > 0.d0) then
         dw2 = dvxl - drhol * dmedc/dmedrho
        else
         dw2 = dvxr - drhor * dmedc/dmedrho
        endif 

        if(dlamu > 0.d0) then
         dw3 = dvyl
         dw4 = dvzl
        else
         dw3 = dvyr
         dw4 = dvzr
        endif

        ! Riemann primitive variables
        dbrho = (dw1 - dw2) / (2.d0 * dmedc/dmedrho)
        dbux  = (dw1 + dw2) / 2.d0
        dbuy  = dw3
        dbuz  = dw4
        dbpress = eos(dbrho)

        ! Calculate fluxes
        frho = dbrho * dbux
        frvx = dbrho * dbux * dbux + dbpress
        frvy = dbrho * dbuy * dbux
        frvz = dbrho * dbuz * dbux

      end subroutine riemann
      !**************************************************
 


