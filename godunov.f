      !****************************************************
      ! godunov method
      !****************************************************
      subroutine godunov
        use mesh
        use omp_lib
        implicit none

        ! variables 
        integer i, k, l, m, ism, ksm, lsm
        logical is_empty, is_xp, is_xm, is_yp, is_ym, is_zp, is_zm
        double precision drhol, drhor, dpressl, dpressr
        double precision dvxl, dvxr, dvyl, dvyr, dvzl, dvzr
        double precision frho1, frvx1, frvy1, frvz1
        double precision frho2, frvx2, frvy2, frvz2
        double precision frho3, frvx3, frvy3, frvz3
        double precision frho4, frvx4, frvy4, frvz4
        double precision fxplusrho, fxplusrvx, fxplusrvy, fxplusrvz
        double precision fxminusrho, fxminusrvx, fxminusrvy, fxminusrvz
        double precision fyplusrho, fyplusrvx, fyplusrvy, fyplusrvz
        double precision fyminusrho, fyminusrvx, fyminusrvy, fyminusrvz
        double precision fzplusrho, fzplusrvx, fzplusrvy, fzplusrvz
        double precision fzminusrho, fzminusrvx, fzminusrvy, fzminusrvz
        double precision dfix, dfiy, dfiz, drho, hloc

        ! godunov method
        do m=1,Ndepth
         ! calculate actual h
         hloc = h/(2.d0**(m-1))
!$omp parallel do 
!$omp& private(i,k,l,ism,ksm,lsm)
!$omp& private(is_empty,is_xp,is_xm,is_yp,is_ym,is_zp,is_zm)
!$omp& private(drhol,drhor,dpressl,dpressr) 
!$omp& private(dvxl,dvxr,dvyl,dvyr,dvzl,dvzr)
!$omp& private(frho1,frvx1,frvy1,frvz1)
!$omp& private(frho2,frvx2,frvy2,frvz2)
!$omp& private(frho3,frvx3,frvy3,frvz3)
!$omp& private(frho4,frvx4,frvy4,frvz4)
!$omp& private(fxplusrho,fxplusrvx,fxplusrvy,fxplusrvz)
!$omp& private(fxminusrho,fxminusrvx,fxminusrvy,fxminusrvz)
!$omp& private(fyplusrho,fyplusrvx,fyplusrvy,fyplusrvz)
!$omp& private(fyminusrho,fyminusrvx,fyminusrvy,fyminusrvz)
!$omp& private(fzplusrho,fzplusrvx,fzplusrvy,fzplusrvz)
!$omp& private(fzminusrho,fzminusrvx,fzminusrvy,fzminusrvz)
!$omp& private(dfix,dfiy,dfiz,drho)
!$omp& shared(m,hloc,tau)
!$omp& shared(density,dimplsx,dimplsy,dimplsz)
!$omp& shared(rho,velx,vely,velz,press,gravity)
!$omp& schedule(dynamic)
!$omp& num_threads(MIC_NUM_THREADS)
         do l=2,N+1
          do k=2,N+1
           do i=2,N+1

            ! current cell is empty
            is_empty = (m < Ndepth) 
     =                 .and. (i > N/4+1)
     =                 .and. (i < 3*N/4+2) 
     =                 .and. (k > N/4+1) 
     =                 .and. (k < 3*N/4+2)
     =                 .and. (l > N/4+1) 
     =                 .and. (l < 3*N/4+2)

            ! current cell is not empty
            if( .not. is_empty  ) then

             ! internal boundary
             is_xp    = (m < Ndepth) 
     =                  .and. (i == N/4+1) 
     =                  .and. (k > N/4+1) 
     =                  .and. (k < 3*N/4+2)
     =                  .and. (l > N/4+1) 
     =                  .and. (l < 3*N/4+2)
             is_xm    = (m < Ndepth) 
     =                  .and. (i == 3*N/4+2) 
     =                  .and. (k > N/4+1) 
     =                  .and. (k < 3*N/4+2)
     =                  .and. (l > N/4+1) 
     =                  .and. (l < 3*N/4+2)
             is_yp    = (m < Ndepth) 
     =                  .and. (i > N/4+1) 
     =                  .and. (i < 3*N/4+2) 
     =                  .and. (k == N/4+1) 
     =                  .and. (l > N/4+1) 
     =                  .and. (l < 3*N/4+2)
             is_ym    = (m < Ndepth) 
     =                  .and. (i > N/4+1) 
     =                  .and. (i < 3*N/4+2) 
     =                  .and. (k == 3*N/4+2)
     =                  .and. (l > N/4+1) 
     =                  .and. (l < 3*N/4+2)
             is_zp    = (m < Ndepth) 
     =                  .and. (i > N/4+1) 
     =                  .and. (i < 3*N/4+2) 
     =                  .and. (k > N/4+1) 
     =                  .and. (k < 3*N/4+2)
     =                  .and. (l == N/4+1)
             is_zm    = (m < Ndepth) 
     =                  .and. (i > N/4+1) 
     =                  .and. (i < 3*N/4+2) 
     =                  .and. (k > N/4+1) 
     =                  .and. (k < 3*N/4+2)
     =                  .and. (l == 3*N/4+2)
 
 
             ! x plus flux
             drhol   = rho(i,k,l,m).vl   + 0.5d0*rho(i,k,l,m).dx
             dpressl = press(i,k,l,m).vl + 0.5d0*press(i,k,l,m).dx
             dvxl    = velx(i,k,l,m).vl  + 0.5d0*velx(i,k,l,m).dx
             dvyl    = vely(i,k,l,m).vl  + 0.5d0*vely(i,k,l,m).dx
             dvzl    = velz(i,k,l,m).vl  + 0.5d0*velz(i,k,l,m).dx

             if( is_xp ) then
               ism = 2
               ksm = 2*(k - N/4 - 1)
               lsm = 2*(l - N/4 - 1)

               drhor   = rho(ism,ksm,lsm,m+1).vl   - 
     =                    0.5d0*rho(ism,ksm,lsm,m+1).dx
               dpressr = press(ism,ksm,lsm,m+1).vl - 
     =                    0.5d0*press(ism,ksm,lsm,m+1).dx
               dvxr    = velx(ism,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velx(ism,ksm,lsm,m+1).dx
               dvyr    = vely(ism,ksm,lsm,m+1).vl  - 
     =                    0.5d0*vely(ism,ksm,lsm,m+1).dx
               dvzr    = velz(ism,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velz(ism,ksm,lsm,m+1).dx
               call riemann(frho1, frvx1, frvy1, frvz1,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)

               drhor   = rho(ism,ksm+1,lsm,m+1).vl   - 
     =                    0.5d0*rho(ism,ksm+1,lsm,m+1).dx
               dpressr = press(ism,ksm+1,lsm,m+1).vl - 
     =                    0.5d0*press(ism,ksm+1,lsm,m+1).dx
               dvxr    = velx(ism,ksm+1,lsm,m+1).vl  - 
     =                    0.5d0*velx(ism,ksm+1,lsm,m+1).dx
               dvyr    = vely(ism,ksm+1,lsm,m+1).vl  - 
     =                    0.5d0*vely(ism,ksm+1,lsm,m+1).dx
               dvzr    = velz(ism,ksm+1,lsm,m+1).vl  - 
     =                    0.5d0*velz(ism,ksm+1,lsm,m+1).dx
               call riemann(frho2, frvx2, frvy2, frvz2,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)

               drhor   = rho(ism,ksm,lsm+1,m+1).vl   - 
     =                    0.5d0*rho(ism,ksm,lsm+1,m+1).dx
               dpressr = press(ism,ksm,lsm+1,m+1).vl - 
     =                    0.5d0*press(ism,ksm,lsm+1,m+1).dx
               dvxr    = velx(ism,ksm,lsm+1,m+1).vl  - 
     =                    0.5d0*velx(ism,ksm,lsm+1,m+1).dx
               dvyr    = vely(ism,ksm,lsm+1,m+1).vl  - 
     =                    0.5d0*vely(ism,ksm,lsm+1,m+1).dx
               dvzr    = velz(ism,ksm,lsm+1,m+1).vl  - 
     =                    0.5d0*velz(ism,ksm,lsm+1,m+1).dx
               call riemann(frho3, frvx3, frvy3, frvz3,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)

               drhor   = rho(ism,ksm+1,lsm+1,m+1).vl   - 
     =                    0.5d0*rho(ism,ksm+1,lsm+1,m+1).dx
               dpressr = press(ism,ksm+1,lsm+1,m+1).vl - 
     =                    0.5d0*press(ism,ksm+1,lsm+1,m+1).dx
               dvxr    = velx(ism,ksm+1,lsm+1,m+1).vl  - 
     =                    0.5d0*velx(ism,ksm+1,lsm+1,m+1).dx
               dvyr    = vely(ism,ksm+1,lsm+1,m+1).vl  - 
     =                    0.5d0*vely(ism,ksm+1,lsm+1,m+1).dx
               dvzr    = velz(ism,ksm+1,lsm+1,m+1).vl  - 
     =                    0.5d0*velz(ism,ksm+1,lsm+1,m+1).dx
               call riemann(frho4, frvx4, frvy4, frvz4,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)

               fxplusrho = (frho1 + frho2 + frho3 + frho4)/4.d0
               fxplusrvx = (frvx1 + frvx2 + frvx3 + frvx4)/4.d0
               fxplusrvy = (frvy1 + frvy2 + frvy3 + frvy4)/4.d0
               fxplusrvz = (frvz1 + frvz2 + frvz3 + frvz4)/4.d0
             else
               drhor   = rho(i+1,k,l,m).vl   - 0.5d0*rho(i+1,k,l,m).dx
               dpressr = press(i+1,k,l,m).vl - 0.5d0*press(i+1,k,l,m).dx
               dvxr    = velx(i+1,k,l,m).vl  - 0.5d0*velx(i+1,k,l,m).dx
               dvyr    = vely(i+1,k,l,m).vl  - 0.5d0*vely(i+1,k,l,m).dx
               dvzr    = velz(i+1,k,l,m).vl  - 0.5d0*velz(i+1,k,l,m).dx
               call riemann(fxplusrho, fxplusrvx, fxplusrvy, fxplusrvz,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)
             endif

             ! x minus flux
             drhor   = rho(i,k,l,m).vl   - 0.5d0*rho(i,k,l,m).dx
             dpressr = press(i,k,l,m).vl - 0.5d0*press(i,k,l,m).dx
             dvxr    = velx(i,k,l,m).vl  - 0.5d0*velx(i,k,l,m).dx
             dvyr    = vely(i,k,l,m).vl  - 0.5d0*vely(i,k,l,m).dx
             dvzr    = velz(i,k,l,m).vl  - 0.5d0*velz(i,k,l,m).dx

             if( is_xm ) then
               ism = N + 1
               ksm = 2*(k - N/4 - 1)
               lsm = 2*(l - N/4 - 1)

               drhol   = rho(ism,ksm,lsm,m+1).vl   + 
     =                    0.5d0*rho(ism,ksm,lsm,m+1).dx
               dpressl = press(ism,ksm,lsm,m+1).vl + 
     =                    0.5d0*press(ism,ksm,lsm,m+1).dx
               dvxl    = velx(ism,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velx(ism,ksm,lsm,m+1).dx
               dvyl    = vely(ism,ksm,lsm,m+1).vl  + 
     =                    0.5d0*vely(ism,ksm,lsm,m+1).dx
               dvzl    = velz(ism,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velz(ism,ksm,lsm,m+1).dx
               call riemann(frho1, frvx1, frvy1, frvz1,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)

               drhol   = rho(ism,ksm+1,lsm,m+1).vl   + 
     =                    0.5d0*rho(ism,ksm+1,lsm,m+1).dx
               dpressl = press(ism,ksm+1,lsm,m+1).vl + 
     =                    0.5d0*press(ism,ksm+1,lsm,m+1).dx
               dvxl    = velx(ism,ksm+1,lsm,m+1).vl  + 
     =                    0.5d0*velx(ism,ksm+1,lsm,m+1).dx
               dvyl    = vely(ism,ksm+1,lsm,m+1).vl  + 
     =                    0.5d0*vely(ism,ksm+1,lsm,m+1).dx
               dvzl    = velz(ism,ksm+1,lsm,m+1).vl  + 
     =                    0.5d0*velz(ism,ksm+1,lsm,m+1).dx
               call riemann(frho2, frvx2, frvy2, frvz2,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)

               drhol   = rho(ism,ksm,lsm+1,m+1).vl   + 
     =                    0.5d0*rho(ism,ksm,lsm+1,m+1).dx
               dpressl = press(ism,ksm,lsm+1,m+1).vl + 
     =                    0.5d0*press(ism,ksm,lsm+1,m+1).dx
               dvxl    = velx(ism,ksm,lsm+1,m+1).vl  + 
     =                    0.5d0*velx(ism,ksm,lsm+1,m+1).dx
               dvyl    = vely(ism,ksm,lsm+1,m+1).vl  + 
     =                    0.5d0*vely(ism,ksm,lsm+1,m+1).dx
               dvzl    = velz(ism,ksm,lsm+1,m+1).vl  + 
     =                    0.5d0*velz(ism,ksm,lsm+1,m+1).dx
               call riemann(frho3, frvx3, frvy3, frvz3,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)

               drhol   = rho(ism,ksm+1,lsm+1,m+1).vl   + 
     =                    0.5d0*rho(ism,ksm+1,lsm+1,m+1).dx
               dpressl = press(ism,ksm+1,lsm+1,m+1).vl + 
     =                    0.5d0*press(ism,ksm+1,lsm+1,m+1).dx
               dvxl    = velx(ism,ksm+1,lsm+1,m+1).vl  + 
     =                    0.5d0*velx(ism,ksm+1,lsm+1,m+1).dx
               dvyl    = vely(ism,ksm+1,lsm+1,m+1).vl  + 
     =                    0.5d0*vely(ism,ksm+1,lsm+1,m+1).dx
               dvzl    = velz(ism,ksm+1,lsm+1,m+1).vl  + 
     =                    0.5d0*velz(ism,ksm+1,lsm+1,m+1).dx
               call riemann(frho4, frvx4, frvy4, frvz4,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)

               fxminusrho = (frho1 + frho2 + frho3 + frho4)/4.d0
               fxminusrvx = (frvx1 + frvx2 + frvx3 + frvx4)/4.d0
               fxminusrvy = (frvy1 + frvy2 + frvy3 + frvy4)/4.d0
               fxminusrvz = (frvz1 + frvz2 + frvz3 + frvz4)/4.d0
             else
               drhol   = rho(i-1,k,l,m).vl   + 0.5d0*rho(i-1,k,l,m).dx
               dpressl = press(i-1,k,l,m).vl + 0.5d0*press(i-1,k,l,m).dx
               dvxl    = velx(i-1,k,l,m).vl  + 0.5d0*velx(i-1,k,l,m).dx
               dvyl    = vely(i-1,k,l,m).vl  + 0.5d0*vely(i-1,k,l,m).dx
               dvzl    = velz(i-1,k,l,m).vl  + 0.5d0*velz(i-1,k,l,m).dx
               call riemann(fxminusrho,fxminusrvx,fxminusrvy,fxminusrvz, 
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvxl, dvxr, dvyl, dvyr, dvzl, dvzr)
             endif

             ! y plus flux
             drhol   = rho(i,k,l,m).vl   + 0.5d0*rho(i,k,l,m).dy
             dpressl = press(i,k,l,m).vl + 0.5d0*press(i,k,l,m).dy
             dvxl    = velx(i,k,l,m).vl  + 0.5d0*velx(i,k,l,m).dy
             dvyl    = vely(i,k,l,m).vl  + 0.5d0*vely(i,k,l,m).dy
             dvzl    = velz(i,k,l,m).vl  + 0.5d0*velz(i,k,l,m).dy

             if( is_yp ) then
               ism = 2*(i - N/4 - 1)
               ksm = 2
               lsm = 2*(l - N/4 - 1)

               drhor   = rho(ism,ksm,lsm,m+1).vl   - 
     =                    0.5d0*rho(ism,ksm,lsm,m+1).dy
               dpressr = press(ism,ksm,lsm,m+1).vl - 
     =                    0.5d0*press(ism,ksm,lsm,m+1).dy
               dvxr    = velx(ism,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velx(ism,ksm,lsm,m+1).dy
               dvyr    = vely(ism,ksm,lsm,m+1).vl  - 
     =                    0.5d0*vely(ism,ksm,lsm,m+1).dy
               dvzr    = velz(ism,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velz(ism,ksm,lsm,m+1).dy
               call riemann(frho1, frvy1, frvz1, frvx1,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)

               drhor   = rho(ism+1,ksm,lsm,m+1).vl   - 
     =                    0.5d0*rho(ism+1,ksm,lsm,m+1).dy
               dpressr = press(ism+1,ksm,lsm,m+1).vl - 
     =                    0.5d0*press(ism+1,ksm,lsm,m+1).dy
               dvxr    = velx(ism+1,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velx(ism+1,ksm,lsm,m+1).dy
               dvyr    = vely(ism+1,ksm,lsm,m+1).vl  - 
     =                    0.5d0*vely(ism+1,ksm,lsm,m+1).dy
               dvzr    = velz(ism+1,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velz(ism+1,ksm,lsm,m+1).dy
               call riemann(frho2, frvy2, frvz2, frvx2,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)

               drhor   = rho(ism,ksm,lsm+1,m+1).vl   - 
     =                    0.5d0*rho(ism,ksm,lsm+1,m+1).dy
               dpressr = press(ism,ksm,lsm+1,m+1).vl - 
     =                    0.5d0*press(ism,ksm,lsm+1,m+1).dy
               dvxr    = velx(ism,ksm,lsm+1,m+1).vl  - 
     =                    0.5d0*velx(ism,ksm,lsm+1,m+1).dy
               dvyr    = vely(ism,ksm,lsm+1,m+1).vl  - 
     =                    0.5d0*vely(ism,ksm,lsm+1,m+1).dy
               dvzr    = velz(ism,ksm,lsm+1,m+1).vl  - 
     =                    0.5d0*velz(ism,ksm,lsm+1,m+1).dy
               call riemann(frho3, frvy3, frvz3, frvx3,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)

               drhor   = rho(ism+1,ksm,lsm+1,m+1).vl   - 
     =                    0.5d0*rho(ism+1,ksm,lsm+1,m+1).dy
               dpressr = press(ism+1,ksm,lsm+1,m+1).vl - 
     =                    0.5d0*press(ism+1,ksm,lsm+1,m+1).dy
               dvxr    = velx(ism+1,ksm,lsm+1,m+1).vl  - 
     =                    0.5d0*velx(ism+1,ksm,lsm+1,m+1).dy
               dvyr    = vely(ism+1,ksm,lsm+1,m+1).vl  - 
     =                    0.5d0*vely(ism+1,ksm,lsm+1,m+1).dy
               dvzr    = velz(ism+1,ksm,lsm+1,m+1).vl  - 
     =                    0.5d0*velz(ism+1,ksm,lsm+1,m+1).dy
               call riemann(frho4, frvy4, frvz4, frvx4,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)

               fyplusrho = (frho1 + frho2 + frho3 + frho4)/4.d0
               fyplusrvx = (frvx1 + frvx2 + frvx3 + frvx4)/4.d0
               fyplusrvy = (frvy1 + frvy2 + frvy3 + frvy4)/4.d0
               fyplusrvz = (frvz1 + frvz2 + frvz3 + frvz4)/4.d0
             else
               drhor   = rho(i,k+1,l,m).vl   - 0.5d0*rho(i,k+1,l,m).dy
               dpressr = press(i,k+1,l,m).vl - 0.5d0*press(i,k+1,l,m).dy
               dvxr    = velx(i,k+1,l,m).vl  - 0.5d0*velx(i,k+1,l,m).dy
               dvyr    = vely(i,k+1,l,m).vl  - 0.5d0*vely(i,k+1,l,m).dy
               dvzr    = velz(i,k+1,l,m).vl  - 0.5d0*velz(i,k+1,l,m).dy
               call riemann(fyplusrho, fyplusrvy, fyplusrvz, fyplusrvx,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)
             endif

             ! y minus flux
             drhor   = rho(i,k,l,m).vl   - 0.5d0*rho(i,k,l,m).dy
             dpressr = press(i,k,l,m).vl - 0.5d0*press(i,k,l,m).dy
             dvxr    = velx(i,k,l,m).vl  - 0.5d0*velx(i,k,l,m).dy
             dvyr    = vely(i,k,l,m).vl  - 0.5d0*vely(i,k,l,m).dy
             dvzr    = velz(i,k,l,m).vl  - 0.5d0*velz(i,k,l,m).dy

             if( is_ym )  then
               ism = 2*(i - N/4 - 1)
               ksm = N + 1
               lsm = 2*(l - N/4 - 1)

               drhol   = rho(ism,ksm,lsm,m+1).vl   + 
     =                    0.5d0*rho(ism,ksm,lsm,m+1).dy
               dpressl = press(ism,ksm,lsm,m+1).vl + 
     =                    0.5d0*press(ism,ksm,lsm,m+1).dy
               dvxl    = velx(ism,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velx(ism,ksm,lsm,m+1).dy
               dvyl    = vely(ism,ksm,lsm,m+1).vl  + 
     =                    0.5d0*vely(ism,ksm,lsm,m+1).dy
               dvzl    = velz(ism,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velz(ism,ksm,lsm,m+1).dy
               call riemann(frho1, frvy1, frvz1, frvx1,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)

               drhol   = rho(ism+1,ksm,lsm,m+1).vl   + 
     =                    0.5d0*rho(ism+1,ksm,lsm,m+1).dy
               dpressl = press(ism+1,ksm,lsm,m+1).vl + 
     =                    0.5d0*press(ism+1,ksm,lsm,m+1).dy
               dvxl    = velx(ism+1,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velx(ism+1,ksm,lsm,m+1).dy
               dvyl    = vely(ism+1,ksm,lsm,m+1).vl  + 
     =                    0.5d0*vely(ism+1,ksm,lsm,m+1).dy
               dvzl    = velz(ism+1,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velz(ism+1,ksm,lsm,m+1).dy
               call riemann(frho2, frvy2, frvz2, frvx2,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)

               drhol   = rho(ism,ksm,lsm+1,m+1).vl   + 
     =                    0.5d0*rho(ism,ksm,lsm+1,m+1).dy
               dpressl = press(ism,ksm,lsm+1,m+1).vl + 
     =                    0.5d0*press(ism,ksm,lsm+1,m+1).dy
               dvxl    = velx(ism,ksm,lsm+1,m+1).vl  + 
     =                    0.5d0*velx(ism,ksm,lsm+1,m+1).dy
               dvyl    = vely(ism,ksm,lsm+1,m+1).vl  + 
     =                    0.5d0*vely(ism,ksm,lsm+1,m+1).dy
               dvzl    = velz(ism,ksm,lsm+1,m+1).vl  + 
     =                    0.5d0*velz(ism,ksm,lsm+1,m+1).dy
               call riemann(frho3, frvy3, frvz3, frvx3,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)

               drhol   = rho(ism+1,ksm,lsm+1,m+1).vl   + 
     =                    0.5d0*rho(ism+1,ksm,lsm+1,m+1).dy
               dpressl = press(ism+1,ksm,lsm+1,m+1).vl + 
     =                    0.5d0*press(ism+1,ksm,lsm+1,m+1).dy
               dvxl    = velx(ism+1,ksm,lsm+1,m+1).vl  + 
     =                    0.5d0*velx(ism+1,ksm,lsm+1,m+1).dy
               dvyl    = vely(ism+1,ksm,lsm+1,m+1).vl  + 
     =                    0.5d0*vely(ism+1,ksm,lsm+1,m+1).dy
               dvzl    = velz(ism+1,ksm,lsm+1,m+1).vl  + 
     =                    0.5d0*velz(ism+1,ksm,lsm+1,m+1).dy
               call riemann(frho4, frvy4, frvz4, frvx4,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)

               fyminusrho = (frho1 + frho2 + frho3 + frho4)/4.d0
               fyminusrvx = (frvx1 + frvx2 + frvx3 + frvx4)/4.d0
               fyminusrvy = (frvy1 + frvy2 + frvy3 + frvy4)/4.d0
               fyminusrvz = (frvz1 + frvz2 + frvz3 + frvz4)/4.d0
             else
               drhol   = rho(i,k-1,l,m).vl   + 0.5d0*rho(i,k-1,l,m).dy
               dpressl = press(i,k-1,l,m).vl + 0.5d0*press(i,k-1,l,m).dy
               dvxl    = velx(i,k-1,l,m).vl  + 0.5d0*velx(i,k-1,l,m).dy
               dvyl    = vely(i,k-1,l,m).vl  + 0.5d0*vely(i,k-1,l,m).dy
               dvzl    = velz(i,k-1,l,m).vl  + 0.5d0*velz(i,k-1,l,m).dy
               call riemann(fyminusrho,fyminusrvy,fyminusrvz,fyminusrvx,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvyl, dvyr, dvzl, dvzr, dvxl, dvxr)
             endif 

             ! z plus flux
             drhol   = rho(i,k,l,m).vl   + 0.5d0*rho(i,k,l,m).dz
             dpressl = press(i,k,l,m).vl + 0.5d0*press(i,k,l,m).dz
             dvxl    = velx(i,k,l,m).vl  + 0.5d0*velx(i,k,l,m).dz
             dvyl    = vely(i,k,l,m).vl  + 0.5d0*vely(i,k,l,m).dz
             dvzl    = velz(i,k,l,m).vl  + 0.5d0*velz(i,k,l,m).dz

             if( is_zp ) then
               ism = 2*(i - N/4 - 1)
               ksm = 2*(k - N/4 - 1)
               lsm = 2

               drhor   = rho(ism,ksm,lsm,m+1).vl   - 
     =                    0.5d0*rho(ism,ksm,lsm,m+1).dz
               dpressr = press(ism,ksm,lsm,m+1).vl - 
     =                    0.5d0*press(ism,ksm,lsm,m+1).dz
               dvxr    = velx(ism,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velx(ism,ksm,lsm,m+1).dz
               dvyr    = vely(ism,ksm,lsm,m+1).vl  - 
     =                    0.5d0*vely(ism,ksm,lsm,m+1).dz
               dvzr    = velz(ism,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velz(ism,ksm,lsm,m+1).dz
               call riemann(frho1, frvz1, frvx1, frvy1,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)

               drhor   = rho(ism+1,ksm,lsm,m+1).vl   - 
     =                    0.5d0*rho(ism+1,ksm,lsm,m+1).dz
               dpressr = press(ism+1,ksm,lsm,m+1).vl - 
     =                    0.5d0*press(ism+1,ksm,lsm,m+1).dz
               dvxr    = velx(ism+1,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velx(ism+1,ksm,lsm,m+1).dz
               dvyr    = vely(ism+1,ksm,lsm,m+1).vl  - 
     =                    0.5d0*vely(ism+1,ksm,lsm,m+1).dz
               dvzr    = velz(ism+1,ksm,lsm,m+1).vl  - 
     =                    0.5d0*velz(ism+1,ksm,lsm,m+1).dz
               call riemann(frho2, frvz2, frvx2, frvy2,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)

               drhor   = rho(ism,ksm+1,lsm,m+1).vl   - 
     =                    0.5d0*rho(ism,ksm+1,lsm,m+1).dz
               dpressr = press(ism,ksm+1,lsm,m+1).vl - 
     =                    0.5d0*press(ism,ksm+1,lsm,m+1).dz
               dvxr    = velx(ism,ksm+1,lsm,m+1).vl  - 
     =                    0.5d0*velx(ism,ksm+1,lsm,m+1).dz
               dvyr    = vely(ism,ksm+1,lsm,m+1).vl  - 
     =                    0.5d0*vely(ism,ksm+1,lsm,m+1).dz
               dvzr    = velz(ism,ksm+1,lsm,m+1).vl  - 
     =                    0.5d0*velz(ism,ksm+1,lsm,m+1).dz
               call riemann(frho3, frvz3, frvx3, frvy3,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)

               drhor   = rho(ism+1,ksm+1,lsm,m+1).vl   - 
     =                    0.5d0*rho(ism+1,ksm+1,lsm,m+1).dz
               dpressr = press(ism+1,ksm+1,lsm,m+1).vl - 
     =                    0.5d0*press(ism+1,ksm+1,lsm,m+1).dz
               dvxr    = velx(ism+1,ksm+1,lsm,m+1).vl  - 
     =                    0.5d0*velx(ism+1,ksm+1,lsm,m+1).dz
               dvyr    = vely(ism+1,ksm+1,lsm,m+1).vl  - 
     =                    0.5d0*vely(ism+1,ksm+1,lsm,m+1).dz
               dvzr    = velz(ism+1,ksm+1,lsm,m+1).vl  - 
     =                    0.5d0*velz(ism+1,ksm+1,lsm,m+1).dz
               call riemann(frho4, frvz4, frvx4, frvy4,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)

               fzplusrho = (frho1 + frho2 + frho3 + frho4)/4.d0
               fzplusrvx = (frvx1 + frvx2 + frvx3 + frvx4)/4.d0
               fzplusrvy = (frvy1 + frvy2 + frvy3 + frvy4)/4.d0
               fzplusrvz = (frvz1 + frvz2 + frvz3 + frvz4)/4.d0
             else
               drhor   = rho(i,k,l+1,m).vl   - 0.5d0*rho(i,k,l+1,m).dz
               dpressr = press(i,k,l+1,m).vl - 0.5d0*press(i,k,l+1,m).dz
               dvxr    = velx(i,k,l+1,m).vl  - 0.5d0*velx(i,k,l+1,m).dz
               dvyr    = vely(i,k,l+1,m).vl  - 0.5d0*vely(i,k,l+1,m).dz
               dvzr    = velz(i,k,l+1,m).vl  - 0.5d0*velz(i,k,l+1,m).dz
               call riemann(fzplusrho, fzplusrvz, fzplusrvx, fzplusrvy,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)
             endif

             ! z minus flux
             drhor   = rho(i,k,l,m).vl   - 0.5d0*rho(i,k,l,m).dz
             dpressr = press(i,k,l,m).vl - 0.5d0*press(i,k,l,m).dz
             dvxr    = velx(i,k,l,m).vl  - 0.5d0*velx(i,k,l,m).dz
             dvyr    = vely(i,k,l,m).vl  - 0.5d0*vely(i,k,l,m).dz
             dvzr    = velz(i,k,l,m).vl  - 0.5d0*velz(i,k,l,m).dz

             if( is_zm ) then
               ism = 2*(i - N/4 - 1)
               ksm = 2*(k - N/4 - 1)
               lsm = N + 1

               drhol   = rho(ism,ksm,lsm,m+1).vl   + 
     =                    0.5d0*rho(ism,ksm,lsm,m+1).dz
               dpressl = press(ism,ksm,lsm,m+1).vl + 
     =                    0.5d0*press(ism,ksm,lsm,m+1).dz
               dvxl    = velx(ism,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velx(ism,ksm,lsm,m+1).dz
               dvyl    = vely(ism,ksm,lsm,m+1).vl  + 
     =                    0.5d0*vely(ism,ksm,lsm,m+1).dz
               dvzl    = velz(ism,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velz(ism,ksm,lsm,m+1).dz
               call riemann(frho1, frvz1, frvx1, frvy1,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)

               drhol   = rho(ism+1,ksm,lsm,m+1).vl   + 
     =                    0.5d0*rho(ism+1,ksm,lsm,m+1).dz
               dpressl = press(ism+1,ksm,lsm,m+1).vl + 
     =                    0.5d0*press(ism+1,ksm,lsm,m+1).dz
               dvxl    = velx(ism+1,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velx(ism+1,ksm,lsm,m+1).dz
               dvyl    = vely(ism+1,ksm,lsm,m+1).vl  + 
     =                    0.5d0*vely(ism+1,ksm,lsm,m+1).dz
               dvzl    = velz(ism+1,ksm,lsm,m+1).vl  + 
     =                    0.5d0*velz(ism+1,ksm,lsm,m+1).dz
               call riemann(frho2, frvz2, frvx2, frvy2,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)

               drhol   = rho(ism,ksm+1,lsm,m+1).vl   + 
     =                    0.5d0*rho(ism,ksm+1,lsm,m+1).dz
               dpressl = press(ism,ksm+1,lsm,m+1).vl + 
     =                    0.5d0*press(ism,ksm+1,lsm,m+1).dz
               dvxl    = velx(ism,ksm+1,lsm,m+1).vl  + 
     =                    0.5d0*velx(ism,ksm+1,lsm,m+1).dz
               dvyl    = vely(ism,ksm+1,lsm,m+1).vl  + 
     =                    0.5d0*vely(ism,ksm+1,lsm,m+1).dz
               dvzl    = velz(ism,ksm+1,lsm,m+1).vl  + 
     =                    0.5d0*velz(ism,ksm+1,lsm,m+1).dz
               call riemann(frho3, frvz3, frvx3, frvy3,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)

               drhol   = rho(ism+1,ksm+1,lsm,m+1).vl   + 
     =                    0.5d0*rho(ism+1,ksm+1,lsm,m+1).dz
               dpressl = press(ism+1,ksm+1,lsm,m+1).vl + 
     =                    0.5d0*press(ism+1,ksm+1,lsm,m+1).dz
               dvxl    = velx(ism+1,ksm+1,lsm,m+1).vl  + 
     =                    0.5d0*velx(ism+1,ksm+1,lsm,m+1).dz
               dvyl    = vely(ism+1,ksm+1,lsm,m+1).vl  + 
     =                    0.5d0*vely(ism+1,ksm+1,lsm,m+1).dz
               dvzl    = velz(ism+1,ksm+1,lsm,m+1).vl  + 
     =                    0.5d0*velz(ism+1,ksm+1,lsm,m+1).dz
               call riemann(frho4, frvz4, frvx4, frvy4,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)

               fzminusrho = (frho1 + frho2 + frho3 + frho4)/4.d0
               fzminusrvx = (frvx1 + frvx2 + frvx3 + frvx4)/4.d0
               fzminusrvy = (frvy1 + frvy2 + frvy3 + frvy4)/4.d0
               fzminusrvz = (frvz1 + frvz2 + frvz3 + frvz4)/4.d0
             else
               drhol   = rho(i,k,l-1,m).vl   + 0.5d0*rho(i,k,l-1,m).dz
               dpressl = press(i,k,l-1,m).vl + 0.5d0*press(i,k,l-1,m).dz
               dvxl    = velx(i,k,l-1,m).vl  + 0.5d0*velx(i,k,l-1,m).dz
               dvyl    = vely(i,k,l-1,m).vl  + 0.5d0*vely(i,k,l-1,m).dz
               dvzl    = velz(i,k,l-1,m).vl  + 0.5d0*velz(i,k,l-1,m).dz
               call riemann(fzminusrho,fzminusrvz,fzminusrvx,fzminusrvy,
     =                      drhol, drhor, dpressl, dpressr,
     =                      dvzl, dvzr, dvxl, dvxr, dvyl, dvyr)
             endif

             ! gravity term
             drho = density(i,k,l,m)
             dfix = (gravity(i+1,k,l,m) - gravity(i-1,k,l,m))/2.d0/hloc
             dfiy = (gravity(i,k+1,l,m) - gravity(i,k-1,l,m))/2.d0/hloc
             dfiz = (gravity(i,k,l+1,m) - gravity(i,k,l-1,m))/2.d0/hloc

             if(dfix > hloc/tau/tau) then
               print *,"Strong X signal"
               stop
             endif

             if(dfiy > hloc/tau/tau) then
               print *,"Strong Y signal"
               stop
             endif

             if(dfiz > hloc/tau/tau) then
               print *,"Strong Z signal"
               stop
             endif

             ! Godunov method
             density(i,k,l,m) = density(i,k,l,m) - tau * ( 
     =                          (fxplusrho - fxminusrho)/hloc + 
     =                          (fyplusrho - fyminusrho)/hloc + 
     =                          (fzplusrho - fzminusrho)/hloc )
             dimplsx(i,k,l,m) = dimplsx(i,k,l,m) - tau * ( 
     =                          (fxplusrvx - fxminusrvx)/hloc + 
     =                          (fyplusrvx - fyminusrvx)/hloc + 
     =                          (fzplusrvx - fzminusrvx)/hloc ) -
     =                          tau * drho * dfix
             dimplsy(i,k,l,m) = dimplsy(i,k,l,m) - tau * ( 
     =                          (fxplusrvy - fxminusrvy)/hloc + 
     =                          (fyplusrvy - fyminusrvy)/hloc + 
     =                          (fzplusrvy - fzminusrvy)/hloc ) -
     =                          tau * drho * dfiy
             dimplsz(i,k,l,m) = dimplsz(i,k,l,m) - tau * ( 
     =                          (fxplusrvz - fxminusrvz)/hloc + 
     =                          (fyplusrvz - fyminusrvz)/hloc + 
     =                          (fzplusrvz - fzminusrvz)/hloc ) -
     =                          tau * drho * dfiz

            endif

           enddo   
          enddo
         enddo
!$omp end parallel do 
        enddo

      end subroutine godunov
      !**************************************************

