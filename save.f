      !****************************************************
      ! save results
      !****************************************************
      subroutine save(timer)
        use machine
        use mesh 
        implicit none 

        ! input timer
        double precision, intent(in) :: timer

        ! variables 
        integer i, k, l, m
        logical is_empty
        double precision xcord, ycord, zcord, radius
        double precision hloc, centerloc
        character*32 filenamexy, filenamexz

201     format(f21.14,1x,f21.14,1x,f21.14)

        ! save results
        do m=1,Ndepth
         ! compute actual h and center for local mesh
         hloc      = h     /(2.d0**(m-1))
         centerloc = length/(2.d0**m)

         ! create file for any mesh
         write(filenamexy, 101) timer,m
101      format(f6.3,"_xy_m",i5.5,".dat")
         open(501, file = filenamexy)

         write(filenamexz, 102) timer,m
102      format(f6.3,"_xz_m",i5.5,".dat")
         open(502, file = filenamexz)

         ! save XY surface
         l = N/2 + 1
         do k=2,N+1
          do i=2,N+1
           xcord = (i-1)*hloc - centerloc - hloc/2.d0
           ycord = (k-1)*hloc - centerloc - hloc/2.d0
           zcord = (l-1)*hloc - centerloc - hloc/2.d0
           radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)
           write(501,201) xcord,ycord,density(i,k,l,m)
          enddo 
         enddo

         ! save XZ surface
         k = N/2 + 1
         do l=2,N+1
          do i=2,N+1
           xcord = (i-1)*hloc - centerloc - hloc/2.d0
           ycord = (k-1)*hloc - centerloc - hloc/2.d0
           zcord = (l-1)*hloc - centerloc - hloc/2.d0
           radius = dsqrt(xcord*xcord + ycord*ycord + zcord*zcord)
           write(502,201) xcord,zcord,density(i,k,l,m)
          enddo 
         enddo

         close(501)
         close(502)

        enddo


      end subroutine save
      !**************************************************

