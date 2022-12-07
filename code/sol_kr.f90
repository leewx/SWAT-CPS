      subroutine sol_kr
      !! produce by lwx 2022-5-12 ����sol_kr��ֵ
      use parm
        real*8 :: tempkr1,tempkr2, tempkr3, tempkr4, tempkr5, tempkr6, Ep_lwx
        tempkr1 = 0.
        tempkr2 = 0.
        tempkr3 = 0.
        tempkr4 = 0.
        tempkr5 = 0.
        tempkr6 = 0.
        Ep_lwx = 0.
        z_lwx = 300
        j = 0
        j = ihru
                    !����Ա���ˮ�������ȵĵĽ�����
                    ! Se = (��-��r)/(��s - ��r)=
                    !1/(1+ABS(��h)^n)^m   h<0    hΪѹ��ˮͷ
                    ! �� = ��r +����s -��r ����1+ABS(��h)^n��^-m
                    ! K = sol_k *Kr
                    ! Kr = Se^0.5*(1-(1-Se^1/m)^m)^2
                    ! �˴�ΪK������
                    !����ROSETTA����ÿһ��������vanģ�Ͳ���
        lwx_m(j) = 1 - (1/lwx_npar_data(j))
        if (lwx_theta_r_data(j) < 0 ) then
            lwx_kr(j) = 1
            return
        end if
        if (days_pot == 1 ) then 
                theta_sw(days_pot,j) = lwx_theta_s_data(j)
                es_day_last(days_pot,j) = es_day
        end if 
        if (days_pot > 1) then  !!��ֹ��һ��
        if (H_pot_max_idl(days_pot) < 0.2 .and. H_pot_min_idl(days_pot) < 0.2 .and. H_pot_max_idl(days_pot) > 0 .and. H_pot_min_idl(days_pot) > 0) then
            if (H_pot_max_idl(days_pot-1) > 0.2 .and. H_pot_min_idl(days_pot-1) > 0.2) then
                theta_sw(days_pot,j) = lwx_theta_s_data(j)
                es_day_last(days_pot,j) = es_day
            end if
!!�ڱ�㣬ֲ����ˮΪ0��ֻ������
!!            Tp = pet_day_last * (1 - exp(-k_lwx*laiday(j-1)))
!            Ep = pet_day_last * (exp(-k_lwx * laiday(j-1)))
!            pet_day_last = pet_day
        end if
        if (H_pot_max_idl(days_pot) < 0.2 .and. H_pot_min_idl(days_pot) < 0.2 .and. H_pot_max_idl(days_pot) > 0 .and. H_pot_min_idl(days_pot) > 0) then 
            if(H_pot_max_idl(days_pot-1) == 0.1 .and. H_pot_min_idl(days_pot-1) == 0.1) then 
                Ep_lwx = es_day_last(days_pot-1,j)*(z_lwx/(z_lwx+exp(2.374-0.00713*z_Lwx)))
                theta_sw(days_pot,j) = (theta_sw(days_pot-1,j)*z_lwx - Ep_lwx)/z_lwx
                es_day_last(days_pot,j) = es_day 
            end if 
        end if
        end if 
        lwx_se(j) = (theta_sw(days_pot,j) - lwx_theta_r_data(j))/(lwx_theta_s_data(j) - lwx_theta_r_data(j))
        if (lwx_se(j) < 0 )then 
            lwx_se(j) = 0.0001
        end if 
        if (lwx_se(j) > 1) then 
            lwx_se(j) = 1
        end if 
        tempkr1 = 1/(lwx_m(j))
        tempkr2 = lwx_se(j) ** tempkr1
        tempkr3 = 1 - tempkr2
        tempkr4 = tempkr3 ** lwx_m(j)
        tempkr5 = 1 - tempkr4
        tempkr6 = tempkr5 ** 2
        lwx_kr(j) = tempkr6 * lwx_se(j)**0.5
        if (lwx_kr(j) > 1 )then 
            lwx_kr(j) = 1
        end if 
        if (lwx_kr(j) < 0 )then 
            lwx_kr(j) = 0.000001
        end if 
    return 
    end 