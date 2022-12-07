      subroutine irr_rch

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine performs the irrigation operation when the water
!!    source is a reach

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name            |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    aird(:)         |mm H2O        |amount of water applied to HRU on current
!!                                   |day
!!    auto_wstr(:)    |none or mm    |water stress factor which triggers auto
!!                                   |irrigation
!!    divmax(:)       |mm H2O or     |maximum daily irrigation diversion from
!!                    |  10^4 m^3 H2O|the reach (when IRR=1): when value is
!!                                   |positive the units are mm H2O; when the
!!                                   |value is negative, the units are (10**4
!!                                   |m^3 H2O
!!    flag                           |1 = manual 2 = auto
!!    flowfr(:)       |none          |fraction of available flow in reach that
!!                                   |is allowed to be applied to the HRU
!!    flowmin(:)      |m**3/s        |minimum instream flow for irrigation
!!                                   |diversions when IRR=1, irrigation water
!!                                   |will be diverted only when streamflow is
!!                                   |at or above FLOWMIN.
!!    iida            |julian date   |day being simulated (current julian date)
!!    wstrs_id(:)     |none          |water stress identifier:
!!                                   |1 plant water demand
!!                                   |2 soil water deficit
!!    inum1           |none          |reach number
!!    ipot(:)         |none          |number of HRU (in subbasin) that is ponding
!!                                   |water--the HRU that the surface runoff from
!!                                   |current HRU drains into. This variable is
!!                                   |used only for rice paddys or closed
!!                                   |depressional areas
!!    irramt(:)       |mm H2O        |depth of irrigation water applied to
!!                                   |HRU
!!    irrno(:)        |none          |irrigation source location
!!                                   |if IRR=1, IRRNO is the number of the
!!                                   |          reach
!!                                   |if IRR=2, IRRNO is the number of the
!!                                   |          reservoir
!!                                   |if IRR=3, IRRNO is the number of the
!!                                   |          subbasin
!!                                   |if IRR=4, IRRNO is the number of the
!!                                   |          subbasin
!!                                   |if IRR=5, not used
!!    irrsc(:)        |none          |irrigation source code:
!!                                   |1 divert water from reach
!!                                   |2 divert water from reservoir
!!                                   |3 divert water from shallow aquifer
!!                                   |4 divert water from deep aquifer
!!                                   |5 divert water from source outside
!!                                   |  watershed
!!    nair(:)         |none          |sequence number of auto-irrigation
!!                                   |application within the year
!!    nhru            |none          |number of HRUs in watershed
!!    nirr(:)         |none          |sequence number of irrigation application
!!                                   |within the year
!!    nro(:)          |none          |sequence number of year in rotation
!!    phuacc(:)       |none          |fraction of plant heat units accumulated
!!    pot_vol(:)      |m**3 H2O      |current volume of water stored in the
!!                                   |depression/impounded area
!!    rtwtr           |m^3 H2O       |water leaving reach on day
!!    sedrch          |metric tons   |sediment transported out of reach on day
!!    sol_sumfc(:)    |mm H2O        |amount of water held in the soil profile
!!                                   |at field capacity
!!    sol_sw(:)       |mm H2O        |amount of water stored in soil profile on any
!!                                   |given day
!!    strsw(:)        |none          |fraction of potential plant growth achieved
!!                                   |on the day where the reduction is caused by
!!                                   |water stress
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    nirr(:)     |none          |sequence number of irrigation application
!!                               |within the year
!!    pot_vol(:)  |mm            |current volume of water stored in the
!!                               |depression/impounded area
!!    rtwtr       |m^3 H2O       |water leaving reach on day
!!    sedrch      |metric tons   |sediment transported out of reach on day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    cnv         |none          |conversion factor (mm => m^3)
!!    flag        |none          |irrigation flag:
!!                               |0 no irrigation operation on current day
!!                               |1 scheduled irrigation
!!                               |2 auto irrigation
!!    jrch        |none          |reach number
!!    k           |none          |HRU number
!!    vminmm      |mm H2O        |maximum amount of water available for
!!                               |irrigation from reach
!!    vmm         |mm H2O        |depth of irrigation water over HRU
!!    vmxi        |mm H2O        |amount of water specified in irrigation
!!                               |operation
!!    vol         |m^3 H2O       |volume of water applied in irrigation 
!!                               |operation
!!    wtrin       |m^3 H2O       |water outflow from reach prior to subtracting
!!                               |irrigation diversions
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: abs, Min
!!    SWAT: irrigate

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: jrch, k, flag, ii, ID
      real*8 :: high, ggxx, pot_runoff
      real*8 :: cnv, vmm, vminmm, vol, wtrin
      real*8 :: pppp
      integer :: release_day 
      release_day = 10 
      pot_runoff = 0.
      ggxx=0
      jrch = 0
      jrch = inum1
      
      wtrin = 0.
      wtrin = rtwtr + rchstor(jrch)
      do ii = 1, mhru
          pot_vol_temp(ii) = pot_volxmm(ii) !!����ƽ�������
      end do
      do k = 1, nhru
          !! check for timing of irrigation operation
          flag = 0
          flag = irr_flag(k)
          if (auto_wstr(k) > 0.) then
            if (wstrs_id(k) == 1 .and. strsw(k) < auto_wstr(k)) flag = 2
            if (wstrs_id(k) == 2 .and. sol_sumfc(k) - sol_sw(k) >       
     &             auto_wstr(k)) flag = 2
          end if

            !! Set parameters based on manual or auto irrigation
			if (flag == 1) then
			  sq_rto = irrsq(k)
			  irrsc(k) = irr_sc(k)                              
			  irrno(k) = irr_no(k)
			else
			  sq_rto = irr_asq(k)
			  irrsc(k) = irr_sca(k)
			  irrno(k) = irr_noa(k)                   
			endif

        if (pot_fr(k) > 0.001) then
          aird(k) = 0.                                          

          if (flag > 0) then
            !!irrigate only if flow is greater than minimum flow
            if (rtwtr > flowmin(k) * 86400.) then
              cnv = 0.
              cnv = hru_ha(k) * 10.

              vmm = 0.
              vminmm = 0.
              !! compute maximum amount of water allowed in HRU
              if (divmax(k) < 0.) then
                !!divmax units are 10^4 m^3
                vmm = abs(divmax(k)) * 10000. / cnv
              else
                !! divmax units are mm H2O
                vmm = divmax(k)
              endif
              !! compute maximum amount of water available for irrigation
              !! from reach
              wtr_avail = rtwtr + rchstor(jrch)
              vminmm = (wtr_avail - flowmin(k) * 86400.) * flowfr(k)/cnv
              vmm = Min(vminmm, vmm)

              !! check available against set amount in scheduled operation
              if (flag == 1) then
                vmxi = 0.
                vmxi = irramt(k)                       
                if (vmxi < 1.e-6) vmxi = sol_sumfc(k)
                if (vmm > vmxi) vmm = vmxi
              end if
              if (flag == 2) then
                vmxi = 0.
                vmxi = irr_mx(k)
                if (vmm > vmxi) vmm = vmxi
              end if

              if (vmm > 0.) then
                vol = 0.
                vol = vmm * cnv
              ID = 1
              high = 0
              end if 
       !! ���������
       !! 2022-4-18 ����һ����Χ���䣬Hmax,Hmin��ʹ�����ʵʩ��ȴ�ʩ����ˮ��С��Hminʱ�򲹳�ˮ��(Hmin+Hmax)/2��ˮ������Hmax��ˮ��(Hmin+Hmax)/2
              H_pot_max(k) = H_pot_max_idl(days_pot)
              H_pot_min(k) = H_pot_min_idl(days_pot)
              if (water_high_warn(days_pot) == 0) then !!���ڿ��Ƹ�ģ���Ƿ�����
              if(pot_volxmm(k) > pot_vol_temp(k) .and. count_fl>1) then  !!���֮ǰ�޸Ĺ������ˮ,�ֶ��ͷ�ˮ��ƽ����ˮ
                      pot_volxmm(k) = pot_volxmm(k) - (pot_volxmm(k) -
     &                pot_vol_temp(k))/release_day
                      count_fl = count_fl - 1
                  end if   
          if (pot_volxmm(k) == pot_vol_temp(k) .and. count_fl == 1) then               
                      pot_volxmm(k) = pot_vol_temp(k) !!���ƽʱ״̬����������Ϊ����ˮ
                    if (H_pot_max(k) > pot_volxmm(k)) then 
                      H_pot_max(k) = pot_volxmm(k)
                    end if 
                    if (H_pot_min(k) < 0) then 
                      H_pot_min(k) = 0 
                    end if
                    H_pot_avg(k) = (H_pot_min(k)+H_pot_max(k))/2
                    !!��Χ����,�����������max��Ҫ����С���ָ��� 
                    if (pot_vol(k)< H_pot_min(k)) then !!���pot�ĺ�ˮ��ˮ����С��ˮ��ҪС������Havg(k) - pot_vol(k)�ĸ߶�
                        vol=10 * (H_pot_avg(k) - pot_vol(k)) * hru_ha(k)
                        pot_vol(k) = pot_vol(k) + vol / (10. * hru_ha(k))
                        pppp = (H_pot_avg(k) - pot_vol(k))
                        irr(k) = pppp
                        call irrigate(k,pppp)
                    end if
         if (pot_vol(k)> H_pot_min(k).and.pot_vol(k)< H_pot_max(k)) then
                        vol=0
                    end if !!���pot��ˮ��Hmin �� Hmax֮�䣬�Ͳ����
                    if(pot_vol(k) > H_pot_max(k)) then
                      vol = 0 
                      pot_runoff = pot_vol(k) - H_pot_avg(k)  !!���potˮ̫�࣬������������ޣ�����ˮ����ˮ���ƽ��ֵ����֤ˮ������
                      pot_runoff = 10 * pot_runoff * hru_ha(k)
                      pot_vol(k) = H_pot_avg(k)
                    end if 
              end if
              end if
                !!�޸Ľ��� !!����Ϊˮ����ʽ
                if (water_high_warn(days_pot) == 20) then !!�����������
                    pot_volxmm(k) = 2500 !!�����ˮ�������������ṩ��ˮ��
                    count_fl = release_day !! 10�����ſ�
                    if (pot_vol(k) < pot_volxmm(k)) then 
                        high = pot_volxmm(k) - pot_vol(k) !!����
                        vol=10 * high * hru_ha(k)
                    end if 
                    ggxx = (wtr_avail - flowmin(k) * 86400.) * flowfr(k)
                    vol = Min (vol,ggxx)
                    if (pot_fr(k) > 1.e-6) then !!�����ˮ���͸���ˮ��
                      pot_vol(k) = pot_vol(k) + vol / (10. * hru_ha(k))
                    else
                    call irrigate(k,vmm)
                    end if    
                end if
                !!����Ϊ����ģʽ����һ�׶�����
                if (pot_fr(k) > 1.e-6) then
                  vol = 0.
                  vol = aird(k) * cnv
                end if
                ggxx = (wtr_avail - flowmin(k) * 86400.) * flowfr(k)
                vol = Min (vol,ggxx)
                if (pot_runoff > 0) then  !!�����Ҫ��ˮ������ºӵ�ˮ��
                    rchstor(jrch) = rchstor(jrch) + pot_runoff
                    rtwtr = rtwtr + pot_runoff
                    rtwtr = amax1(0., rtwtr)
                end if 
                if (vol > 0) then
                if (vol > rchstor(jrch)) then
                  vol = vol - rchstor(jrch)                                         !!��ˮ
                  rchstor(jrch) = 0.
                else
                  rchstor(jrch) = rchstor(jrch) - vol
                  vol = 0.
                end if
                if (vol > 0.) then
                  rtwtr = rtwtr - vol
                  rtwtr = dmax1(0., rtwtr)
                end if
                end if
                if (rchstor(jrch) < 0) then
                    rchstor(jrch)=0
                end if
                !!�޸Ľ���
       !!         if (ipot(k) == k) then
                !if (pot_fr(k) > 1.e-6) then
                !  pot_vol(k) = pot_vol(k) + vol / (10. * potsa(k))
                !else

                !end if

                !! subtract irrigation from reach outflow
           !!     if (ipot(k) /= k) then
                if (pot_fr(k) > 1.e-6) then
                  vol = 0.
                  vol = aird(k) * cnv
                end if
                if (ievent > 0) then
                  do ii = 1, nstep
                    hrtwtr(ii) = hrtwtr(ii) - vol * hrtwtr(ii) / rtwtr
                    if (hrtwtr(ii) < 0.) hrtwtr(ii) = 0.
                  end do
                end if
!!                xx = vol     							                           !! BN: replaced "wtrin" with "vol"
!!                vol = vol / irr_eff(k)   !! BN: inserted to account for irr. efficiency                                             
!                xx = (wtr_avail - flowmin(k) * 86400.) * flowfr(k)                 !! BN: inserted: xx = available/allowed amount in m3/s
!                xx = Min(xx, vol)                                                  !! BN: inserted dabstracted water cannot be more than allowed/available amount
!                if (xx > rchstor(jrch)) then
!                  xx = vol - rchstor(jrch)                                         !! BN: replaced "wtrin" with "vol"
!                  rchstor(jrch) = 0.
!                else
!                  rchstor(jrch) = rchstor(jrch) - xx
!                  xx = 0.
!                end if
!                if (xx > 0.) then
!                  rtwtr = rtwtr - xx
!                  rtwtr = dmax1(0., rtwtr)
!                end if

                !! advance irrigation operation number
                if (flag == 1) then
                  nirr(k) = nirr(k) + 1
                end if
            
            
20001  format (a5,a4,1x,2i5,9f10.2)
          end if
          end if
          end if
          
      end do

      if (wtrin /= rtwtr .and. wtrin > 0.01) then
        sedrch = sedrch * rtwtr / wtrin

        rch_san = rch_san * rtwtr / wtrin
        rch_sil = rch_sil * rtwtr / wtrin
        rch_cla = rch_cla * rtwtr / wtrin
        rch_sag = rch_sag * rtwtr / wtrin
        rch_lag = rch_lag * rtwtr / wtrin
        rch_gra = rch_gra * rtwtr / wtrin

        if (sedrch  < 1.e-6) then
	    sedrch = 0.
	    rch_san = 0.
          rch_sil = 0.
          rch_cla = 0.
          rch_sag = 0.
          rch_lag = 0.
          rch_gra = 0.
	  end if

        if (ievent > 0) then
          do ii = 1, nstep
            hsedyld(ii) = hsedyld(ii) * rtwtr / wtrin
          end do
        end if
      end if

      return
      end