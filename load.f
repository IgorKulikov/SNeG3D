      !****************************************************
      ! load physics problem subroutine
      !****************************************************
      subroutine load
        use mesh
        implicit none
        double precision, external :: density_profile 
        double precision, external :: gravity_profile
        double precision, external :: velocity_profile

        ! variables 
        double precision densval
        integer i, k, l, m

        ! compute length of basis cell
        h = length/N

        ! compute center of domain
	center = length/2.d0

        ! load init profiles
        do m=1,Ndepth
         do l=2,N+1
          do k=2,N+1
           do i=2,N+1

            ! set density
            densval = density_profile(i,k,l,m)
            density(i,k,l,m) = densval

            ! set momentum of impulse
            dimplsx(i,k,l,m) = densval * velocity_profile(i,k,l,0,m)
            dimplsy(i,k,l,m) = densval * velocity_profile(i,k,l,1,m)
            dimplsz(i,k,l,m) = densval * velocity_profile(i,k,l,2,m)

            ! set gravity
            gravity(i,k,l,m) = gravity_profile(i,k,l,m)

           enddo   
          enddo
         enddo
        enddo

      end subroutine load
      !**************************************************

