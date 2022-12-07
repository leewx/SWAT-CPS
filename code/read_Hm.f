      subroutine read_Hm
      !! produce by lwx 2022-4-14 
      use parm
      integer dao_i 
      integer dao_N
      real, dimension (:), allocatable :: Init_Month_read
      real, dimension (:), allocatable :: Init_day_read
      real, dimension (:), allocatable :: End_Month_read
      real, dimension (:), allocatable :: End_Day_read
      real, dimension (:), allocatable :: H_pot_max_read
      real, dimension (:), allocatable :: H_pot_min_read
      real, dimension (:), allocatable :: Init_Month
      real, dimension (:), allocatable :: Init_day
      real, dimension (:), allocatable :: End_Month
      real, dimension (:), allocatable :: End_Day
      real, dimension (:), allocatable :: map_month
      integer map_month_sum,temp_day1, temp2_month   !!ע�⵽ȫ������Ϊreal��
      dao_N = 19 !! �ı��ļ�������
      temp2_month = 1
      map_month_sum = 0
      temp_day1 = 0
        if (curyr == nbyr .and. idal > 0) then
          idlst = idal
        else
          idlst = 366 - leapyr
        end if
      allocate (map_month(12))
      allocate (Init_Month(idlst))
      allocate (Init_day(idlst))
      allocate (End_Month(idlst))
      allocate (End_Day(idlst))
      allocate (Init_Month_read(dao_N))
      allocate (Init_day_read(dao_N))
      allocate (End_Month_read(dao_N))
      allocate (End_Day_read(dao_N))
      allocate (H_pot_max_read(dao_N))
      allocate (H_pot_min_read(dao_N))
      map_month(1) = 31
      if (leapyr == 1) then
          map_month(2) = 28
      else
          map_month(2) = 29
      end if
      map_month(3) = 31
      map_month(4) = 30
      map_month(5) = 31
      map_month(6) = 30
      map_month(7) = 31
      map_month(8) = 31
      map_month(9) = 30
      map_month(10) = 31
      map_month(11) = 30
      map_month(12) = 31 !!!�����·ݵ�ӳ��ͼ
      H_pot_max_idl = 0 
      H_pot_min_idl = 0
      open(1008611, file = 'water.rice', status = 'old')
      read(1008611, *) !����head
      do dao_i = 1,dao_N
      read(1008611, *) Init_month_read(dao_i),Init_day_read(dao_i),
     &End_month_read(dao_i),End_day_read(dao_i),H_pot_max_read(dao_i),
     &H_pot_min_read(dao_i)
      end do
      do dao_j =1 ,12 
      if (dao_j < Init_month_read(1)) then !! ����1-12 �£���δ���ý�ˮ���·���Ϊ0��������ˮ
        map_month_sum = map_month_sum + map_month(dao_j)
        temp2_month = temp2_month + 1
      end if 
      end do
      do dao_i =1,map_month_sum 
          H_pot_max_idl(dao_i) = 0 
          H_pot_min_idl(dao_i) = 0
      end do
      temp_day1 = map_month_sum 
      do dao_i = map_month_sum, temp_day1
          H_pot_max_idl(dao_i) = 0 
          H_pot_min_idl(dao_i) = 0
      end do
      do dao_i = temp_day1 + 1, temp_day1 + Init_day_read(1) -1
          H_pot_max_idl(dao_i) = 0 
          H_pot_min_idl(dao_i) = 0          
      end do!! ���϶���������
      temp_day1 = temp_day1 + Init_day_read(1) -1
      do dao_i =1, dao_N !! ����ÿһ�У���max �� min ֵ����idl��
          if (Init_month_read(dao_i) == End_month_read(dao_i)) then 
          do dao_j = temp_day1 +1 , temp_day1 + End_day_read(dao_i) - 
     &       Init_day_read(dao_i) + 1
                      H_pot_max_idl(dao_j) = H_pot_max_read(dao_i) 
                      H_pot_min_idl(dao_j) = H_pot_min_read(dao_i)
          end do
              temp_day1 = temp_day1 + End_day_read(dao_i) - 
     &        Init_day_read(dao_i) + 1
          end if
          if (Init_month_read(dao_i) /= End_month_read(dao_i)) then !!ע��ע��!!�ֽ��ܱ����ؽ��·ݷ�Ϊ��ͬһ�µĲ�����������һ���µĲ�������һ������֧�ֿ�������!!
           if (Init_month_read(dao_i) +1  == End_month_read(dao_i)) then
               do dao_j = temp_day1 +1, temp_day1 +1 +
     &            map_month(Init_month_read(dao_i))-
     &            Init_day_read(dao_i)+ End_day_read(dao_i)
                  H_pot_max_idl(dao_j) = H_pot_max_read(dao_i) 
                  H_pot_min_idl(dao_j) = H_pot_min_read(dao_i)    
               end do
               temp_day1 = temp_day1 +1 +
     &         map_month(Init_month_read(dao_i))-
     &         Init_day_read(dao_i)+ End_day_read(dao_i)
           end if 
          end if 
      end do
      close (1008611)
      return !!��ģ����������
      end