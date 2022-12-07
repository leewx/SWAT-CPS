      subroutine read_matrix
      !! produce by lwx 2022-4-21 
      use parm
      
      integer kkii,count,mmii,stater
      count = 0
      matrix_hru = 0

      open(1008613, file = 'matrix.txt', status = 'old')
      do kkii = 1,9574
          read(1008613, *) (matrix_hru(kkii,mmii),mmii = 1,9574)
!!          do mmii = 1,9574
!!          if (matrix_hru(kkii,mmii) == 1) then 
!!              print *,mmii
!!          end if 
!!          end do
      end do
      close (1008613)
      open(1008614, file = 'slope.txt', status = 'old')
      read (1008614,*) !!ÌøÍ·
      do kkii = 1,9574
          read (1008614,*) slope_hru(kkii)
      end do
      open(1008615, file = 'dem_high.txt',status = 'old')
      read (1008615,*) !!Ìø¹ı
      do kkii = 1,9574
          read (1008615,*) high_hru(kkii)
      end do
      print *,"read_matrix have been processed"
      end