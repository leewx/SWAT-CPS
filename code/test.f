      program main
      integer kkii,count,mmii,stater
      integer,dimension(9574,9574)::matrix_hru
      count = 0
      matrix_hru = 0
      close(1008611)
      open(1008611, file = 'matrix.txt', status = 'old')
      do kkii = 1,9574
          read(1008611, *) (matrix_hru(kkii,mmii),mmii = 1,9574)
!!          do mmii = 1,9574
!!          if (matrix_hru(kkii,mmii) == 1) then 
!!              print *,mmii
!!          end if 
!!          end do
      print *, kkii
      end do 
      end