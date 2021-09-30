      !****************************************************
      ! boundary condition for hydrodynamics
      !****************************************************
      subroutine hydro_boundary(m)
        use mesh
        implicit none 

        ! input mesh level
        integer, intent(in) :: m

        ! variables 
        double precision densval, vxval, vyval, vzval
        integer i, k, l, ib, kb, lb
                        
        ! set boundary conditions
        do l=1,N+2
         do k=1,N+2
           i = 1        ! left boundary condition

           if(m == 1) then    ! neuman for basic mesh
             density(i,k,l,m) = density(i+1,k,l,m)
             dimplsx(i,k,l,m) = dimplsx(i+1,k,l,m)
             dimplsy(i,k,l,m) = dimplsy(i+1,k,l,m)
             dimplsz(i,k,l,m) = dimplsz(i+1,k,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             density(i,k,l,m) = density(ib,kb,lb,m-1)
             dimplsx(i,k,l,m) = dimplsx(ib,kb,lb,m-1)
             dimplsy(i,k,l,m) = dimplsy(ib,kb,lb,m-1)
             dimplsz(i,k,l,m) = dimplsz(ib,kb,lb,m-1)
           endif

           i = N+2      ! right boundary condition

           if(m == 1) then    ! neuman for basic mesh
             density(i,k,l,m) = density(i-1,k,l,m)
             dimplsx(i,k,l,m) = dimplsx(i-1,k,l,m)
             dimplsy(i,k,l,m) = dimplsy(i-1,k,l,m)
             dimplsz(i,k,l,m) = dimplsz(i-1,k,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             density(i,k,l,m) = density(ib,kb,lb,m-1)
             dimplsx(i,k,l,m) = dimplsx(ib,kb,lb,m-1)
             dimplsy(i,k,l,m) = dimplsy(ib,kb,lb,m-1)
             dimplsz(i,k,l,m) = dimplsz(ib,kb,lb,m-1)
           endif

         enddo
        enddo

        do l=1,N+2
         do i=1,N+2
           k = 1        ! top boundary condition

           if(m == 1) then    ! neuman for basic mesh
             density(i,k,l,m) = density(i,k+1,l,m)
             dimplsx(i,k,l,m) = dimplsx(i,k+1,l,m)
             dimplsy(i,k,l,m) = dimplsy(i,k+1,l,m)
             dimplsz(i,k,l,m) = dimplsz(i,k+1,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             density(i,k,l,m) = density(ib,kb,lb,m-1)
             dimplsx(i,k,l,m) = dimplsx(ib,kb,lb,m-1)
             dimplsy(i,k,l,m) = dimplsy(ib,kb,lb,m-1)
             dimplsz(i,k,l,m) = dimplsz(ib,kb,lb,m-1)
           endif

           k = N+2      ! bottom boundary condition

           if(m == 1) then    ! neuman for basic mesh
             density(i,k,l,m) = density(i,k-1,l,m)
             dimplsx(i,k,l,m) = dimplsx(i,k-1,l,m)
             dimplsy(i,k,l,m) = dimplsy(i,k-1,l,m)
             dimplsz(i,k,l,m) = dimplsz(i,k-1,l,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             density(i,k,l,m) = density(ib,kb,lb,m-1)
             dimplsx(i,k,l,m) = dimplsx(ib,kb,lb,m-1)
             dimplsy(i,k,l,m) = dimplsy(ib,kb,lb,m-1)
             dimplsz(i,k,l,m) = dimplsz(ib,kb,lb,m-1)
           endif

         enddo
        enddo

        do k=1,N+2
         do i=1,N+2
           l = 1        ! down boundary condition

           if(m == 1) then    ! neuman for basic mesh
             density(i,k,l,m) = density(i,k,l+1,m)
             dimplsx(i,k,l,m) = dimplsx(i,k,l+1,m)
             dimplsy(i,k,l,m) = dimplsy(i,k,l+1,m)
             dimplsz(i,k,l,m) = dimplsz(i,k,l+1,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             density(i,k,l,m) = density(ib,kb,lb,m-1)
             dimplsx(i,k,l,m) = dimplsx(ib,kb,lb,m-1)
             dimplsy(i,k,l,m) = dimplsy(ib,kb,lb,m-1)
             dimplsz(i,k,l,m) = dimplsz(ib,kb,lb,m-1)
           endif

           l = N+2      ! up boundary condition

           if(m == 1) then    ! neuman for basic mesh
             density(i,k,l,m) = density(i,k,l-1,m)
             dimplsx(i,k,l,m) = dimplsx(i,k,l-1,m)
             dimplsy(i,k,l,m) = dimplsy(i,k,l-1,m)
             dimplsz(i,k,l,m) = dimplsz(i,k,l-1,m)
           else               ! internal mesh
             ib = N/4 + i/2 + 1         
             kb = N/4 + k/2 + 1
             lb = N/4 + l/2 + 1
             density(i,k,l,m) = density(ib,kb,lb,m-1)
             dimplsx(i,k,l,m) = dimplsx(ib,kb,lb,m-1)
             dimplsy(i,k,l,m) = dimplsy(ib,kb,lb,m-1)
             dimplsz(i,k,l,m) = dimplsz(ib,kb,lb,m-1)
           endif

         enddo
        enddo

      end subroutine hydro_boundary
      !**************************************************
