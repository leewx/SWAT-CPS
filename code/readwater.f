      subroutine readwater
      !! produce by lwx 2022-4-17 
      use parm
      integer :: conut, iiii, stater
      count = 0 
      count = 1825
      open(1008612, file = 'water_high.high', status = 'old')
      read(1008612, *) !����head
      do iiii = 1 , count 
          read(1008612 , * ) water_high_c(iiii)
      end do
      close (1008612)
      
      do iiii = 1, count 
          if (water_high_c(iiii) > 34.4) then
              water_high_warn(iiii) = 1 !!ˮλ����34.4�͸���ID
          end if 
      end do 
      !!���ҽ���
      return
      end