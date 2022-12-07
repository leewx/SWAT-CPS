    subroutine readsol_K
    use parm
    integer :: i_lwx_soil 
    i_lwx_soil  = 0
    !!!!!
      open(1008616, file = 'theta_r.txt', status = 'old')
      read(1008616, *) !跳过head
      do i_lwx_soil = 1, 49
          read(1008616, *) lwx_theta_r(i_lwx_soil)
      end do
      close(1008616)
      !!!!!
      open(1008617, file = 'theta_s.txt', status = 'old')
      read(1008617, *) !跳过head
      do i_lwx_soil = 1, 49
          read(1008617, *) lwx_theta_s(i_lwx_soil)
      end do
      close(1008617)
      !!!!!
      open(1008618, file = 'alpha.txt', status = 'old')
      read(1008618, *) !跳过head
      do i_lwx_soil = 1, 49
          read(1008618, *) lwx_alpha(i_lwx_soil)
      end do
      close(1008618)
      !!!!!
      open(1008619, file = 'npar.txt', status = 'old')
      read(1008619, *) !跳过head
      do i_lwx_soil = 1, 49
          read(1008619, *) lwx_npar(i_lwx_soil)
      end do
      close(1008619)
      !!!!! 有49种土壤，把参数存入
      open(1234567,file = 'temp.txt', status = 'old')
      !!建立映射关系
      do i_lwx_soil = 1, mhru 
      if (snam(i_lwx_soil) == 'LEPTOSOLS       ') then 
          soil_k_code(i_lwx_soil) = 1
      end if 
      if (snam(i_lwx_soil) == 'Ferric Lixisols ') then 
          soil_k_code(i_lwx_soil) = 2
      end if 
      if (snam(i_lwx_soil) == 'Rendzic Leptosol') then 
          soil_k_code(i_lwx_soil) = 3
      end if 
      if (snam(i_lwx_soil) == 'Haplic Luvisols2') then 
          soil_k_code(i_lwx_soil) = 4
      end if 
      if (snam(i_lwx_soil) == 'Chromic Luvisols') then 
          soil_k_code(i_lwx_soil) = 5
      end if 
      if (snam(i_lwx_soil) == 'Dystric Cambisol') then 
          soil_k_code(i_lwx_soil) = 6
      end if 
      if (snam(i_lwx_soil) == 'Dystric Cambisols2') then 
          soil_k_code(i_lwx_soil) = 7
      end if 
      if (snam(i_lwx_soil) == 'Calcaric Regosol') then 
          soil_k_code(i_lwx_soil) = 8
      end if 
      if (snam(i_lwx_soil) == 'Dystric Leptosol') then 
          soil_k_code(i_lwx_soil) = 9
      end if 
      if (snam(i_lwx_soil) == 'Humic Cambisols ') then 
          soil_k_code(i_lwx_soil) = 10
      end if 
      if (snam(i_lwx_soil) == 'Calcaric Fluviso') then 
          soil_k_code(i_lwx_soil) = 11
      end if 
      if (snam(i_lwx_soil) == 'Eutric Fluvisols') then 
          soil_k_code(i_lwx_soil) = 12
      end if 
      if (snam(i_lwx_soil) == 'Eutric Fluvisols2') then 
          soil_k_code(i_lwx_soil) = 13
      end if 
      if (snam(i_lwx_soil) == 'Calcaric Fluviso') then 
          soil_k_code(i_lwx_soil) = 14
      end if 
      if (snam(i_lwx_soil) == 'Mollic Gleysols ') then 
          soil_k_code(i_lwx_soil) = 15
      end if 
      if (snam(i_lwx_soil) == 'Cumulic Anthroso') then 
          soil_k_code(i_lwx_soil) = 16
      end if 
      if (snam(i_lwx_soil) == 'Cumulic Anthrosols2') then 
          soil_k_code(i_lwx_soil) = 17
      end if 
      if (snam(i_lwx_soil) == 'Cumulic Anthrosols3') then 
          soil_k_code(i_lwx_soil) = 18
      end if 
      if (snam(i_lwx_soil) == 'Cumulic Anthrosols4') then 
          soil_k_code(i_lwx_soil) = 19
      end if 
      if (snam(i_lwx_soil) == 'Cumulic Anthrosols5') then 
          soil_k_code(i_lwx_soil) = 20
      end if 
      if (snam(i_lwx_soil) == 'Cumulic Anthrosols6') then 
          soil_k_code(i_lwx_soil) = 21
      end if 
      if (snam(i_lwx_soil) == 'Cumulic Anthrosols7') then 
          soil_k_code(i_lwx_soil) = 22
      end if 
      if (snam(i_lwx_soil) == 'Cumulic Anthrosols8') then 
          soil_k_code(i_lwx_soil) = 23
      end if 
      if (snam(i_lwx_soil) == 'Eutric Gleysols ') then 
          soil_k_code(i_lwx_soil) = 24
      end if 
      if (snam(i_lwx_soil) == 'Haplic Acrisols1') then 
          soil_k_code(i_lwx_soil) = 25
      end if 
      if (snam(i_lwx_soil) == 'Haplic Acrisols2') then 
          soil_k_code(i_lwx_soil) = 26
      end if 
      if (snam(i_lwx_soil) == 'Haplic Acrisols3') then 
          soil_k_code(i_lwx_soil) = 27
      end if 
      if (snam(i_lwx_soil) == 'Haplic Acrisols4') then 
          soil_k_code(i_lwx_soil) = 28
      end if 
      if (snam(i_lwx_soil) == 'Haplic Acrisols5') then 
          soil_k_code(i_lwx_soil) = 29
      end if 
      if (snam(i_lwx_soil) == 'Haplic Acrisols6') then 
          soil_k_code(i_lwx_soil) = 30
      end if 
      if (snam(i_lwx_soil) == 'Humic Acrisols1 ') then 
          soil_k_code(i_lwx_soil) = 31
      end if 
      if (snam(i_lwx_soil) == 'Humic Acrisols2 ') then 
          soil_k_code(i_lwx_soil) = 32
      end if 
      if (snam(i_lwx_soil) == 'Humic Acrisols3 ') then 
          soil_k_code(i_lwx_soil) = 33
      end if 
      if (snam(i_lwx_soil) == 'Humic Acrisols4 ') then 
          soil_k_code(i_lwx_soil) = 34
      end if 
      if (snam(i_lwx_soil) == 'Ferric Alisols  ') then 
          soil_k_code(i_lwx_soil) = 35
      end if 
      if (snam(i_lwx_soil) == 'Chromic Cambisol') then 
          soil_k_code(i_lwx_soil) = 36
      end if 
      if (snam(i_lwx_soil) == 'Ferralic Cambiso') then 
          soil_k_code(i_lwx_soil) = 37
      end if 
      if (snam(i_lwx_soil) == 'Ferralic Cambisols2') then 
          soil_k_code(i_lwx_soil) = 38
      end if 
      if (snam(i_lwx_soil) == 'Haplic Alisols1 ') then 
          soil_k_code(i_lwx_soil) = 39
      end if 
      if (snam(i_lwx_soil) == 'Haplic Alisols2 ') then 
          soil_k_code(i_lwx_soil) = 40
      end if 
      if (snam(i_lwx_soil) == 'Haplic Alisols3 ') then 
          soil_k_code(i_lwx_soil) = 41
      end if 
      if (snam(i_lwx_soil) == 'Dystric Cambisol') then 
          soil_k_code(i_lwx_soil) = 42
      end if 
      if (snam(i_lwx_soil) == 'Haplic Luvisols1') then 
          soil_k_code(i_lwx_soil) = 43
      end if 
      if (snam(i_lwx_soil) == 'Haplic Luvisols2') then 
          soil_k_code(i_lwx_soil) = 44
      end if 
      if (snam(i_lwx_soil) == 'Haplic Luvisols3') then 
          soil_k_code(i_lwx_soil) = 45
      end if
      if (snam(i_lwx_soil) == 'Haplic Luvisols ') then 
          soil_k_code(i_lwx_soil) = 45
      end if
      if (snam(i_lwx_soil) == 'Urban, mining, e') then 
          soil_k_code(i_lwx_soil) = 46
      end if 
      if (snam(i_lwx_soil) == 'Water bodies1   ') then 
          soil_k_code(i_lwx_soil) = 47
      end if 
      if (snam(i_lwx_soil) == 'Water bodies2   ') then 
          soil_k_code(i_lwx_soil) = 48
      end if 
      if (snam(i_lwx_soil) == 'Dunes & shift.sa') then 
          soil_k_code(i_lwx_soil) = 49 
      end if
      write(1234567,*) snam(i_lwx_soil)
      end do
      close(1234567)
      do i_lwx_soil =1, mhru
        lwx_theta_r_data(i_lwx_soil) = lwx_theta_r(soil_k_code(i_lwx_soil))
        lwx_theta_s_data(i_lwx_soil) = lwx_theta_s(soil_k_code(i_lwx_soil))
        lwx_alpha_data(i_lwx_soil) = lwx_alpha(soil_k_code(i_lwx_soil))
        lwx_npar_data(i_lwx_soil) = lwx_npar(soil_k_code(i_lwx_soil))
      end do 
      print *,"readsol_K have been processed"

      return
      end