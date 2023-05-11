      subroutine subbasin
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine controls the simulation of the land phase of the 
!!    hydrologic cycle

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name           |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    auto_wstr(:)   |none          |water stress factor which triggers auto
!!                                  |irrigation
!!    bio_e(:)       |(kg/ha)/      |biomass-energy ratio
!!                   |     (MJ/m**2)|The potential (unstressed) growth rate per
!!                                  |unit of intercepted photosynthetically
!!                                  |active radiation.
!!    canev          |mm H2O        |amount of water evaporated from canopy
!!                                  |storage
!!    ep_day         |mm H2O        |actual amount of transpiration that occurs
!!                                  |on day in HRU
!!    es_day         |mm H2O        |actual amount of evaporation (soil et) that
!!                                  |occurs on day in HRU
!!    gw_q(:)        |mm H2O        |groundwater contribution to streamflow from
!!                                  |HRU on current day
!!    hru_ra(:)      |MJ/m^2        |solar radiation for the day in HRU
!!    iida           |julian date   |day being simulated (current julian date)
!!    idplt(:)       |none          |land cover code from crop.dat
!!    igro(:)        |none          |land cover status code
!!                                  |0 no land cover currently growing
!!                                  |1 land cover growing
!!    inum1          |none          |subbasin number
!!    imp_trig(:)    |none          |release/impound action code:
!!                                  |0 begin impounding water
!!                                  |1 release impounded water
!!    irrsc(:)       |none          |irrigation source code:
!!                                  |1 divert water from reach
!!                                  |2 divert water from reservoir
!!                                  |3 divert water from shallow aquifer
!!                                  |4 divert water from deep aquifer
!!                                  |5 divert water from source outside
!!                                  |  watershed
!!    iurban(:)      |none          |urban simulation code:
!!                                  |0  no urban sections in HRU
!!                                  |1  urban sections in HRU, simulate using
!!                                  |   USGS regression equations
!!                                  |2  urban sections in HRU, simulate using
!!                                  |   build up/wash off algorithm
!!    latq(:)        |mm H2O        |total lateral flow in soil profile for the
!!                                  |day in HRU
!!    nafert(:)      |none          |sequence number of auto-fert application
!!                                  |within the year
!!    nair(:)        |none          |sequence number of auto-irrigation 
!!                                  |application within the year
!!    nfert(:)       |none          |sequence number of fertilizer application
!!                                  |within the year
!!    nirr(:)        |none          |sequence number of irrigation application
!!                                  |within the year
!!    nrelease(:)    |none          |sequence number of impound/release
!!                                  |operation within the year
!!    nro(:)         |none          |sequence number of year in rotation
!!    peakr          |m^3/s         |peak runoff rate
!!    pet_day        |mm H2O        |potential evapotranspiration on current
!!                                  |day in HRU
!!    phuacc(:)      |none          |fraction of plant heat units accumulated
!!    phubase(:)     |heat units    |base zero total heat units (used when no
!!                                  |land cover is growing)
!!                                  |pesticide application occurs
!!    pot_fr(:)      |km2/km2       |fraction of HRU area that drains into
!!                                  |pothole
!!    pot_vol(:)     |m**3 H2O      |current volume of water stored in the 
!!                                  |depression/impounded area
!!    precipday      |mm H2O        |precipitation for the day in HRU
!!    qday           |mm H2O        |surface runoff loading to main channel from
!!                                  |HRU for day
!!    qtile          |mm H2O        |drainage tile flow in soil layer for the 
!!                                  |day
!!    sci(:)         |none          |retention coefficient for CN method based
!!                                  |on plant ET
!!    sedyld(:)      |metric tons   |soil loss for day in HRU
!!    smx(:)         |none          |retention coefficient for CN method based
!!                                  |on soil moisture
!!    surfq(:)       |mm H2O        |surface runoff generated on day in HRU
!!    tmn(:)         |deg C         |minimum temperature for the day in HRU
!!    tmpav(:)       |deg C         |average temperature for the day in HRU
!!    tmx(:)         |deg C         |maximum temperature for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    albday      |none          |albedo, the fraction of the solar radiation
!!                               |reflected at the soil surface back into
!!                               |space
!!    etday       |mm H2O        |actual evapotranspiration occuring on day
!!                               |in HRU
!!    ihru        |none          |HRU number
!!    inflpcp     |mm H2O        |amount of precipitation that infiltrates
!!                               |into soil (enters soil)
!!    nafert(:)   |none          |sequence number of auto-fert application
!!                               |within the year
!!    nair(:)     |none          |sequence number of auto-irrigation 
!!                               |application within the year
!!    qdfr        |none          |fraction of water yield that is surface
!!                               |runoff
!!    qdr(:)      |mm H2O        |total amount of water entering main channel
!!                               |for day from HRU
!!    sci(:)      |none          |retention coefficient for CN method based
!!                               |on plant ET
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    d           |
!!    gma         |kPa/deg C     |psychrometric constant
!!    ho          |              |net radiation
!!    j           |none          |HRU number
!!    pet_alpha   |none          |alpha factor in Priestley-Taylor ET 
!!                               |equation
!!    tmpk        |deg K         |average temperature for the day in the HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Max
!!    SWAT: varinit, albedo, solt, surface, percmain, etpot, etact, fert
!!    SWAT: confert, graze, plantmod, nminrl, nitvol, pminrl, gwmod, apply, gwmod_deep
!!    SWAT: washp, decay, pestlch, enrsb, pesty, orgn, psed, nrain, nlch
!!    SWAT: solp, subwq, bacteria, urban, pothole, latsed, surfstor
!!    SWAT: substor, wetland, hrupond, irrsub, autoirr, watuse, watbal
!!    SWAT: sumv, virtual

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j,sb,kk
      real*8 :: tmpk, d, gma, ho, pet_alpha, aphu, phuop
      real*8 :: pot_sum,sqrt_tempi, sqrt_tempj,Q_seep_lwx_i,Q_seep_lwx_j
      real*8 :: H_P_lwx_1, H_P_lwx_2, Q_ban_i, Q_ban_j
      real*8 :: sub_pot_num, temp1,temp3, sort_flag, pot_vol_avg,h_sub
      real*8 :: sro, ssfnlyr, percnlyr, vv, vno3, co,Q_stream_flow_iihru,iino3,spillo_lwx
      real*8 :: xxorgn, wt1orgn, erorgn1,seepage_overflow,seepage_overflow_sed,sedorgn_seep1
      sedorgn_seep1 =0.
      spillo_lwx = 0.
      Q_stream_flow_iihru = 0
      h_sub = 0
      Q_ban_i = 0
      Q_ban_j = 0
      Q_seep_lwx_i = 0
      Q_seep_lwx_j = 0
      Q_seep_lwx = 0
      H_P_lwx_1 = 0
      H_P_lwx_2 = 0
      sqrt_tempi = 0
      sqrt_tempj = 0
      pot_sum = 0.
      pot_vol_avg = 0.
      sub_pot_num = 0
      temp1 = 0
      temp2 = 0
      temp3 = 0
      sort_flag = 1
      landuse_hru_num = 0
      landuse_hru_avg = 0
      do temp1 = temp2 + 1, temp2 + hrutot(inum1)
          if (pot_fr(temp1) > 0.000001) then
              sub_pot_num = sub_pot_num + 1
              pot_flag(temp1) = 1
          end if
      end do
      ! by lwx 2022-5-9
      ihru = 0
      ihru = hru1(inum1) 

      call sub_subbasin

      do iihru = 1, hrutot(inum1)

      j = 0
      j = ihru


      !!by zhang DSSAT tillage
      !!======================
      !!    deptil(:)   |mm  |depth of mixing caused by tillage operation
      !jj is hru number
      if (cswat == 2) then
          if (tillage_switch(ihru) .eq. 1) then
              if (tillage_days(ihru) .ge. 30) then
                    tillage_switch(ihru) = 0
                    tillage_days(ihru) = 0
              else
                    tillage_days(ihru) = tillage_days(ihru) + 1
              end if                
              !tillage_depth(ihru) = dtil
              !tillage_switch(ihru) = .TRUE. 
          end if
      end if
      !!by zhang DSSAT tillage  
      !!====================== 



      call varinit
      if (icr(j) <= 0) icr(j) = 1
      
      i_wtrhru = 0
      idplrot(icr(j),ihru) = idplt(j)
      if (idplt(j) /= 0) then
          if (cpnm(idplt(j)) == "WATR") then
              i_wtrhru = 1
          end if
      endif
      
	if (i_wtrhru == 1) then
         call water_hru
      else 

        !! Simulate land covers other than water

        !! update base zero total heat units
        if (tmpav(j) > 0. .and. phutot(hru_sub(j)) > 0.01) then
           phubase(j) = phubase(j) + tmpav(j) / phutot(hru_sub(j))
        end if
        
        call schedule_ops

        !! calculate albedo for day
        call albedo

        !! calculate soil temperature for soil layers
        call solt
!        j) /= j .and. imp_trig(nro(j),nrelease(j),j)==1)       &  Srini pothole
!
!     &        then             
          !! calculate surface runoff if HRU is not impounded or an 
          !! undrained depression--
          call surface

          !! add surface flow that was routed across the landscape on the previous day
       !!   qday = qday + surfq_ru(j)
       !!   surfq_ru(j) = 0.
          
          !! compute effective rainfall (amount that percs into soil)
          inflpcp = Max(0.,precipday - surfq(j))
!        end if
         
        !! perform management operations
        if (yr_skip(j) == 0) call operatn
          
        if (auto_wstr(j) > 1.e-6 .and. irrsc(j) > 2) call autoirr       
        
        !! perform soil water routing
        call percmain

        !! compute evapotranspiration
        call etpot
!        if (pot_vol(j) < 1.e-6) call etact
        call etact

        !! compute water table depth using climate drivers
        call wattable

        !! new CN method
        if (icn == 1) then 
        sci(j) = sci(j) + pet_day*exp(-cncoef_sub(hru_sub(j))*sci(j)/   
     &    smx(j)) - precipday + qday + qtile + latq(j) + sepbtm(j)
        else if (icn == 2) then 
        sci(j) = sci(j) + pet_day*exp(-cncoef_sub(hru_sub(j))*sci(j)/   
     &    smx(j)) - precipday + qday + latq(j) + sepbtm(j) + qtile
        sci(j) = dmin1(sci(j),smxco * smx(j))
        end if 
        
        !! apply fertilizer/manure in continuous fert operation
        if (icfrt(j) == 1) then
          ndcfrt(j) = ndcfrt(j) + 1
          call confert
        end if
        
        !! apply pesticide in continuous pest operation
        if (icpst(j) == 1) then 
          ndcpst(j) = ndcpst(j) + 1
          call conapply
        end if 
        
        !! remove biomass from grazing and apply manure
        if (igrz(j) == 1) then
          ndeat(j) = ndeat(j) + 1
          call graze
        end if
       
        !! compute crop growth
        call plantmod
        
        !! check for dormancy
        if (igro(j) == 1) call dormant
        !! compute actual ET for day in HRU
        etday = ep_day + es_day + canev

        !! write daily air and soil temperature file
        !! can be uncommmented if needed by user and also in readfile.f

!      write (120,12112) i,j,tmx(j),tmn(j),(sol_tmp(k,j),k=1,sol_nly(j))
!12112  format (2i4,12f8.2)

        !! compute nitrogen and phosphorus mineralization 

      if (cswat == 0) then
        call nminrl
	end if
	if (cswat == 1) then
		call carbon
	end if
	
	!! Add by zhang
	!!=================
	if (cswat == 2) then
	  call carbon_zhang2
	end if
	!! Add by zhang
	!!=================	

        call nitvol
        if (sol_P_model == 1) then
            call pminrl
        else
            call pminrl2
        end if

!!    compute biozone processes in septic HRUs
!!    if 1)current is septic hru and 2)  soil temperature is above zero
	  if (isep_opt(j)/=0.and.iyr>=isep_iyr(j)) then
	   if (sol_tmp(i_sep(j),j) > 0.) call biozone     
	  endif

        !! compute ground water contribution
        call gwmod
        call gwmod_deep

        !! compute pesticide washoff   
        if (precipday >= 2.54) call washp

        !! compute pesticide degradation
        call decay

        !! compute pesticide movement in soil
        call pestlch

        if (surfq(j) > 0. .and. peakr > 1.e-6) then
          if (precipday > 0.) then
            call enrsb(0)
            if (sedyld(j) > 0.) call pesty(0)

		  if (cswat == 0) then
			call orgn(0)
	    end if
	    if (cswat == 1) then
	    
		    call orgncswat(0)
		  end if
		  
		  !! Add by zhang
		  !! ====================
		  if (cswat == 2) then
		    call orgncswat2(0)
		  end if
		  !! Add by zhang
		  !! ====================

            call psed(0)
          end if
        end if

        !! add nitrate in rainfall to soil profile
        call nrain

        !! compute nitrate movement leaching
        call nlch

        !! compute phosphorus movement
        call solp

        !! compute bacteria transport
        call bacteria

        !! compute loadings from urban areas
        if (urblu(j) > 0) then
	     if(ievent == 0) then
	        call urban ! daily simulation
	     else
		     call urbanhr ! subdaily simulation J.Jeong 4/20/2009
	     endif
	  endif
	  
!! Srini Pothole
        !! compute undrained depression/impounded area (eg rice) processes
!        if (pot_fr(j) > 0.) then
!           if (ievent == 0) then   
!          call pothole
!           else
!              call potholehr
!           endif
!        endif
        
        !! compute sediment loading in lateral flow and add to sedyld
        call latsed

        !! compute nutrient loading in groundwater flow
        call gwnutr
        call gw_no3

        !! lag nutrients and sediment in surface runoff
        call surfstor

        !! lag subsurface flow and nitrate in subsurface flow

        call substor

        !! add lateral flow that was routed across the landscape on the previous day
      !!  latq(j) = latq(j) + latq_ru(j)
      !!  latq_ru(j) = 0.
        
        !! compute reduction in pollutants due to edge-of-field filter strip
        if (vfsi(j) >0.)then
          call filter
          if (filterw(j) > 0.) call buffer
        end if
              if (vfsi(j) == 0. .and. filterw(j) > 0.) then 
                call filtw
                call buffer
              end if

	 !! compute reduction in pollutants due to in field grass waterway
         if (grwat_i(j) == 1) then
          call grass_wway
        end if

	 !! compute reduction in pollutants due to in fixed BMP eff
	 !  if (bmp_flag(j) == 1) then
       !   call bmpfixed
       ! end if


        !! compute water yield for HRU
        qdr(j) = qday + latq(j) + gw_q(j) + qtile + gw_qdeep(j)
        if (qdr(j) < 0.) qdr(j) = 0.
        if (qdr(j) > 0.) then
          qdfr = qday / qdr(j)
        else
          qdfr = 0.
        end if

        !! compute wetland processes
        call wetlan

        !! compute pond processes
        if (ievent == 0) then
           call hrupond
        else
           call hrupondhr
        endif
        
!       Srini pothole        
        if (pot_fr(j) > 0.) call pothole
                
        xx = sed_con(j)+soln_con(j)+solp_con(j)+orgn_con(j)+orgp_con(j)
        if (xx > 1.e-6) then
          call urb_bmp
      end if
      !!计算水力传导
        call sol_kr
         !!consumptive water use (ponds, shallow aquifer, deep aquifer)
        call watuse

        !! perform water balance
        call watbal
        
        !! compute chl-a, CBOD and dissolved oxygen loadings
        call subwq

        !! qdayout is surface runoff leaving the hru - after wetlands, ponds, and potholes
        qdayout(j) = qday

      endif

      !! perform output summarization
      call sumv

      !! summarize output for multiple HRUs per subbasin
      !! store reach loadings for new fig method
      call virtual
      aird(j) = 0.
      
      ihru = ihru + 1
      end do
      
      !! 坡度数据读到 slope_hru(ihru) 中，read_matrix(ihru)储存一个邻接矩阵，i,j从0开始读取矩阵
      !! 
      !! 从这里开始，对Pot进行平均，刚刚的计算已经将所有的pot_vol全部计算完毕，现需要一均值 
      if (icode == 1 .and. icodes(idum + 1) /= 1) then !!在计算完所有子流域参数之后，进行平均,计算subbasin 即将进入主河道汇流计算
          !! 对K进行声明
        Each_iihru_flow_go = 0
        Each_iihru_flow_accept = 0
        Each_iihru_no3_go = 0
        Each_iihru_orgn_go = 0
          do iihru = 1,mhru  !!start loop
        sro = 0.
        vv = 0.
        vno3 = 0.
        co = 0.
        jj = 1
        Q_seep_lwx_i = 0
        Q_seep_lwx_j = 0
        Q_stream_flow_iihru = 0
        Q_stream_flow_jjhru = 0
        iino3 = 0
        !!土壤硝酸根参数计算
              do jjhru = 1,mhru   !!对matrix中的每一个hru进行遍历
                  if (pot_fr(iihru)>0 .and.pot_fr(jjhru)>0) then
                  if(matrix_hru(iihru, jjhru) == 1) then
                      sqrt_tempi =1+slope_hru(iihru)**2
                      sqrt_tempj =1+slope_hru(jjhru)**2
                      Lmin1 = AMIN1(SQRT(hru_km(iihru) * 1000000), 
     &                SQRT(hru_km(jjhru)*1000000))  !!米
                      
          H_P_lwx_1 = pot_vol(iihru)*0.001
          H_P_lwx_2 = pot_vol(jjhru)*0.001
                  !!!对SOL_K的矫正 
                  !SOL_K = exp(1.9582 + 0.0308 *SAND -0.6142*SOL_BD-0.01566*SOL_CBN)
                  !对相对饱和水力传导度的的矫正：
                  ! Se = (θ-θr)/(θs - θr)=
                  !1/(1+ABS(αh)^k)^β   h<0    h为压力水头
                  ! θ = θr +（θs -θr ）（1+ABS(αh)^k）^-β= sol_awc(mhru)
                  ! K = sol_k *Kr
                  ! Kr = Se^0.5*(1-(1-Se^1/β)^β)^2
                  ! 此处为K的修正
                  !利用ROSETTA计算每一种土壤的van模型参数
                  !晒田期取消水面的蒸发
                     K_lwx(iihru) = sol_k(1,iihru) *1.74*24
                     K_lwx_nowater(iihru) = K_lwx(iihru) * lwx_kr(iihru)
                 if (H_pot_max_idl(days_pot) < 0.2 .and. H_pot_min_idl(days_pot) < 0.2 
     &           .and. H_pot_max_idl(days_pot) > 0 .and. H_pot_min_idl(days_pot) > 0) then 
                  K_lwx(iihru) = K_lwx_nowater(iihru) !!如果晒田，则开始非饱和
                 end if
                      if (high_hru(iihru) > high_hru(jjhru) .and. pot_vol(iihru) > pot_vol(jjhru)) then !!如果高
                 !!如果非饱和，附加基质势，由于上下界面不连通，取消压力势
                 !!h = 1/alpha*(((Se)^-1/m)-1)^1/n
                    if (lwx_npar_data(iihru) > 0 ) then
                    h_sub = 1/lwx_alpha_data(iihru)* ((lwx_se(iihru)**(-1/lwx_m(iihru)))-1)**lwx_npar_data(iihru)
                    else
                    h_sub = 0
                    end if!! 修正数据缺失地区的基质势等
                    !!!溢流开始
                    if (pot_vol(iihru) > pot_volxmm(iihru)) then
                    Q_stream_flow = pot_vol(iihru) - pot_volxmm(iihru)
                    Q_stream_flow_iihru = Q_stream_flow
                    if (Q_stream_flow < 0) then 
                        Q_stream_flow = 0
                    end if
                    pot_vol(iihru) = pot_vol(iihru) - Q_stream_flow
                    Q_stream_flow = Q_stream_flow * 0.001*hru_km(iihru)
     &              *1000000
              Q_stream_flow_jjhru=(Q_stream_flow*1000)/(hru_km(jjhru)
     &              *1000000)
                    pot_vol(jjhru) = pot_vol(jjhru)+Q_stream_flow_jjhru
                    end if 
                    !!以上是溢流的水量
                  H_lwx_1 = pot_vol(jjhru)*0.001
     &            /SQRT(1+slope_hru(jjhru)**2) !! 米
     
                  H_lwx_2 = (slope_hru(jjhru)*(SQRT(hru_km(jjhru) 
     &            *1000000) - pot_vol(jjhru)*0.001*slope_hru(jjhru)))/
     &            (SQRT (1+ slope_hru(jjhru)**2)+ pot_vol(jjhru)*0.001* !!米
     &            SQRT(1+slope_hru(jjhru)**2))
     
                  H_lwx_jj = (H_lwx_1+H_lwx_2)/2 
                  
                 H_lwx_y =(SQRT(hru_km(jjhru)*1000000)*slope_hru(jjhru))
     &           /(SQRT(1+slope_hru(jjhru)**2))
     
                  H_lwx_3 = ((H_lwx_y *SQRT(1+slope_hru(iihru)**2) +
     &             pot_vol(iihru)*0.001))/(SQRT(1+slope_hru(iihru)**2))
     
                  H_lwx_y31 = (slope_hru(jjhru)*
     &             SQRT((hru_km(jjhru))*1000000))
     &             / (SQRT(1+slope_hru(jjhru)**2))
     
                  H_lwx_K = (SQRT(hru_km(iihru)*1000000)) - 
     &             (slope_hru(iihru)*pot_vol(iihru)*0.001)
     
                  H_lwx_y32 = (H_lwx_K*slope_hru(iihru))
     &            /(SQRT(1+slope_hru(iihru)**2))
     
                  H_lwx_y33 = (pot_vol(iihru) * 0.001*
     &             SQRT(1+slope_hru(iihru)**2))
     
                  H_lwx_4 = H_lwx_y31 + H_lwx_y32 + H_lwx_y33
                  
                  H_lwx_ii = (H_lwx_3+H_lwx_4)/2
                  Lpot = 0.9
                  !!如果晒田，附加基质势,去除压力势
                  if (K_lwx(iihru) > 80) then 
                      K_lwx(iihru) = 80
                  end if  !!对异常地形范围限制
                  tempss1 = K_lwx(iihru)*R7(iihru)*0.001*24
                  tempss2 = Lmin1*pot_vol(iihru)*0.001
                  tempss3 = H_lwx_ii-H_lwx_jj+H_P_lwx_1-H_P_lwx_2
                 Q_seep_lwx = (tempss1*tempss2*tempss3)/Lpot
                  
                 if (H_pot_max_idl(days_pot) < 0.2 .and. H_pot_min_idl(days_pot) < 0.2 
     &           .and. H_pot_max_idl(days_pot) > 0 .and. H_pot_min_idl(days_pot) > 0) then
                 Q_seep_lwx = (R7(iihru)*K_lwx(iihru)*0.001*24*Lmin1*
     &           pot_vol(iihru)*0.001*
     &           (H_lwx_ii-H_lwx_jj+h_sub))/Lpot !!m^3,注意Lmin单位 
                 end if 
                      if (Q_seep_lwx < 0) then
                      
                      print *,'warning'
                      print *, pot_vol(iihru)
                      
                      end if 
                      if(H_lwx_ii<H_lwx_jj) then
                          print *,'big warning'
                      end if 
                 Q_seep_lwx_i = (Q_seep_lwx*1000) 
     &           /(hru_km(iihru) * 1000000)
                 Q_seep_lwx_j = (Q_seep_lwx*1000) 
     &           /(hru_km(jjhru) * 1000000)
                 if (Q_seep_lwx_i > pot_vol(iihru)) then
                     Q_seep_lwx_i = pot_vol(iihru)
                     pot_vol(iihru) = 0
                     Q_seep_lwx = Q_seep_lwx_i*(hru_km(iihru) * 1000000)
     &               /1000
                 Q_seep_lwx_j = Q_seep_lwx*1000 
     &           /(hru_km(jjhru) * 1000000)
                 end if
                 if (Q_seep_lwx_i < pot_vol(iihru)) then
                     if (pot_vol(jjhru)+Q_seep_lwx_j < pot_vol(iihru)-Q_seep_lwx_i) then 
                 pot_vol(iihru) = pot_vol(iihru)-Q_seep_lwx_i
                 pot_vol(jjhru) = pot_vol(jjhru)+Q_seep_lwx_j
                     else
                     avg = (pot_vol(jjhru)*hru_km(jjhru) * 1000000 + pot_vol(iihru)*hru_km(iihru) * 1000000)/
     &                ((hru_km(iihru)+hru_km(jjhru))* 1000000)
                      Q_seep_lwx_j = avg - pot_vol(jjhru)
                      Q_seep_lwx_i = pot_vol(iihru) - avg
                      if (Q_seep_lwx_j < 0 )then 
                          Q_seep_lwx_j = 0
                      end if 
                      if (Q_seep_lwx_i < 0 )then 
                          Q_seep_lwx_i = 0
                      end if 
                      !!传递的水
                      pot_vol(iihru) = pot_vol(iihru) - Q_seep_lwx_i
                      pot_vol(jjhru) = pot_vol(jjhru) + Q_seep_lwx_j
                     end if 
                 end if
                 !!!以上是渗流的水量
                 !!以上是溢流的水量
              end if
              !!如果某一池子的水过高，就排水
              if(pot_vol(jjhru) > H_pot_max_idl(days_pot)) then
                  spillo_lwx = pot_vol(jjhru) - 0.5 * (H_pot_min_idl(days_pot) + H_pot_max_idl(days_pot))
                  pot_vol(jjhru) = 0.5*(H_pot_min_idl(days_Pot) + H_pot_max_idl(days_pot))
                  qdr(jjhru) = qdr(jjhru) + spillo_lwx
              end if 
                  if (high_hru(jjhru) == high_hru(iihru)
     &                .and. iihru /= jjhru) then !! 如果i与j的高度相等，高流向低，则水头损失 = pot高 -pot低
                          Lpot = 0.9
                  Q_seep_lwx = K_lwx(iihru)*0.001*Lmin1*pot_vol(iihru)*0.001*24*((pot_vol(iihru) - pot_vol(jjhru))*0.001)+                         !! K mm/h       
     &           (H_P_lwx_1 - H_P_lwx_2)/Lpot
                 if (H_pot_max_idl(days_pot) < 0.2 .and. H_pot_min_idl(days_pot) < 0.2 
     &           .and. H_pot_max_idl(days_pot) > 0 .and. H_pot_min_idl(days_pot) > 0) then
                 Q_seep_lwx = (K_lwx(iihru)*0.001*24*Lmin1*
     &           pot_vol(iihru)*0.001*
     &           (H_lwx_ii-H_lwx_jj+h_sub))/Lpot !!m^3 /D ,注意Lmin单位 
                 end if 
                 Q_seep_lwx = ABS(Q_seep_lwx)
                 if (pot_vol(iihru)>pot_vol(jjhru)) then
                 Q_seep_lwx_i = Q_seep_lwx*1000 
     &           /(hru_km(iihru) * 1000000)
                 Q_seep_lwx_j = Q_seep_lwx*1000 
     &           /(hru_km(jjhru) * 1000000)
                 if (Q_seep_lwx_i > pot_vol(iihru)) then
                     Q_seep_lwx_i = pot_vol(iihru)
                     pot_vol(iihru) = 0
                     Q_seep_lwx = Q_seep_lwx_i*(hru_km(iihru) * 1000000)
     &               /1000
                 Q_seep_lwx_j = Q_seep_lwx*1000 
     &           /(hru_km(jjhru) * 1000000)
                 end if
                 if (Q_seep_lwx_i < pot_vol(iihru)) then
                 pot_vol(iihru) = pot_vol(iihru)-Q_seep_lwx_i
                 pot_vol(jjhru) = pot_vol(jjhru)+Q_seep_lwx_j
                 end if
                 end if
                 if (pot_vol(iihru)<pot_vol(jjhru)) then
                 Q_seep_lwx_i = Q_seep_lwx*1000 
     &           /(hru_km(iihru) * 1000000)
                 Q_seep_lwx_j = Q_seep_lwx*1000 
     &           /(hru_km(jjhru) * 1000000)
                 if (Q_seep_lwx_j > pot_vol(jjhru)) then
                     Q_seep_lwx_j = pot_vol(jjhru)
                     pot_vol(jjhru) = 0
                     Q_seep_lwx = Q_seep_lwx_j*(hru_km(jjhru) * 1000000)
     &               /1000
                 Q_seep_lwx_i = Q_seep_lwx*1000 
     &           /(hru_km(iihru) * 1000000)
                 end if
                 if (Q_seep_lwx_j < pot_vol(jjhru)) then
                 pot_vol(iihru) = pot_vol(iihru)+Q_seep_lwx_i
                 pot_vol(jjhru) = pot_vol(jjhru)-Q_seep_lwx_j
                 end if
                 end if 
                      end if
                    !!计算iihru 的总流量
                    vv = Q_seep_lwx_i + Q_stream_flow_iihru + 1.e-10
                    ww = -vv / ((1. - anion_excl(iihru)) * sol_ul(jj,iihru))
                    vno3 = sol_no3(jj,iihru) * (1. - Exp(ww))
                    if (vv > 1.e-10)  co = Max(vno3 / vv, 0.)
                    !!计算iihru 的总流量
                    iino3 = co * (Q_seep_lwx_i + Q_stream_flow_iihru)  !!KG  透过土壤带走的NO3量 Kg/hm^2
                    iino3l = 0.2*iino3 *hru_km(iihru)*100  !!换算成kg 0.2为渗流系数
                    if (pot_vol(iihru) /= 0 ) then 
                    iino3l = iino3l + (pot_no3(iihru)/pot_vol(iihru))
     &              *(Q_seep_lwx_i + Q_stream_flow_iihru) !!kg/mm^3 水流中包含的NO3 Kg
                    end if 
                    iino3 = Min(iino3, sol_no3(1,iihru)) !!土壤中no3 单位是kg/hm^2
                    sol_no3(1,iihru) = sol_no3(1,iihru) - iino3
                    pot_no3(jjhru) = pot_no3(jjhru) + iino3l !!pot单位是kg
                    !!以下计算有机氮
                    xxorgn = 0.
                    wt1orgn = 0.
                    erorgn1= 1.
                    xxorgn = sol_orgn(1,iihru) + sol_aorgn(1,iihru) + sol_fon(1,iihru)
                    wt1orgn = sol_bd(1,iihru) * sol_z(1,iihru) / 100.
                    erorgn1 = erorgn(iihru)
                    conc = 0.
                    conc = xxorgn * erorgn1 / wt1orgn
                    sand_usle(iihru) = usle_cfac(j) * usle_mult(j)
                    seepage_overflow = 1000*(Q_seep_lwx_i + Q_stream_flow_iihru+ 1.e-10)/(Lmin1*hru_km(iihru) * 1000000) !!runoff Q mm
                    seepage_overflow_sed = (seepage_overflow * peakr )** .56 * sand_usle(iihru)
                    sedorgn_seep(iihru) = .001 * conc * seepage_overflow_sed / hru_ha(iihru)  !!kg/ha
                    sedorgn_seep(iihru) = bmp_pn(iihru) * sedorgn_seep(iihru)
                    sedorgn(iihru) = sedorgn_seep(iihru) +sedorgn(iihru)
                    sedorgn_seep1 = sedorgn_seep(iihru) *hru_km(iihru)*100  !!转移的有机氮 kg
                  if (iwave <= 0 .and. xxorgn > 1.e-6) then
                   sol_aorgn(1,iihru) = sol_aorgn(1,iihru) - sedorgn_seep(iihru) *(sol_aorgn(1,iihru) / xxorgn)                              
                   sol_orgn(1,iihru) = sol_orgn(1,iihru) - sedorgn_seep(iihru) * (sol_orgn(1,iihru) / xxorgn)
                   sol_fon(1,iihru) = sol_fon(1,iihru) - sedorgn_seep(iihru) * (sol_fon(1,iihru) / xxorgn)

                   if (sol_aorgn(1,iihru) < 0.) then
                     sedorgn(iihru) = sedorgn_seep(iihru) + sol_aorgn(1,iihru)
                     sol_aorgn(1,iihru) = 0.
                   end if

                   if (sol_orgn(1,iihru) < 0.) then
                     sedorgn(iihru) = sedorgn_seep(iihru) + sol_orgn(1,iihru)
                     sol_orgn(1,iihru) = 0.
                   end if

                   if (sol_fon(1,iihru) < 0.) then
                     sedorgn(iihru) = sedorgn_seep(iihru) + sol_fon(1,iihru)
                     sol_fon(1,iihru) = 0.
                   end if
                  end if
                  !!总氮结束
                    Each_iihru_flow_go(iihru) = Each_iihru_flow_go(iihru) + Q_seep_lwx_i + Q_stream_flow_iihru !!释放的水
                    Each_iihru_flow_accept(jjhru) = Each_iihru_flow_accept(jjhru) + Q_seep_lwx_j + Q_stream_flow_jjhru
                    Each_iihru_no3_go(iihru) = Each_iihru_no3_go(iihru) + iino3l
                    Each_iihru_orgn_go(iihru) = Each_iihru_orgn_go(iihru) + sedorgn_seep1
                      Q_ban_i = Q_ban_i + Q_seep_lwx_i
                      Q_ban_j = Q_ban_j + Q_seep_lwx_j
                  end if !!如果有联通的话，计算两个之间的流量传递率 Q = A*K*(水头损失)/L   !!  A = Lmin1 * min(pot_vol(hruii),pot_vol(hrujj)) 5-19修改R10矫正参数
                  end if
              end do
          end do
          do iihru = 1,mhru 
             if (pot_fr(iihru) > 0) then 
             write (10086, 20001) subnum(iihru), hruno(iihru),i, iyr, iida,Each_iihru_flow_go(iihru),Each_iihru_flow_accept(iihru),
     *       Each_iihru_no3_go(iihru),Each_iihru_orgn_go(iihru)
             end if 
          end do 
      print *, Q_ban_i
      print *, Q_ban_j

      end if
20001 format (a5,a4,1x,2i5,9f10.2) 
      !! by lwx 2022-5-9
      if (icode == 1 .and. icodes(idum + 1) /= 1) then !!在计算完所有子流域参数之后，进行平均,计算subbasin 即将进入主河道汇流计算
      do j = 1, mhru
          if (pot_fr(j) > 0) then 
              write (10086001,2000) subnum(j), hruno(j), i, iyr, pot_no3(j)
2000      format (a5,a4,1x,2i5,9f10.2)
          end if 
      end do 
      end if
      !! route 2 landscape units
      if (ils2flag(inum1) > 0) then
      isub = inum1                        ! save the subbasin number
      
      !! calculate outputs from hillslope
      ihout1 = mhyd_bsn + (inum1 - 1) * 4 ! first outflow hyd number
      ihout = ihout1                      ! outflow hyd number
      inum1 = 1                           ! landscape unit number
      inum2 = isub                        ! subbasin number
      call routeunit                      ! hillslope unit
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout
      
      !! calculate outputs from valley bottom
      inum1 = 2                           ! landscape unit number
      ihout = ihout + 1                   ! outflow hyd number
      sumdaru = 0.
      do j = 1, hrutot(isub)
        sumdaru = sumdaru + hru_km(j)
      end do 
      daru_km(inum2,inum1) = sumdaru
      call routeunit                      ! valley bottom unit
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout
      
      !! route output from hillslope across valley bottom
      ihout = ihout + 1                   ! outflow hyd number
      inum1 = 2                           ! valley bottom landscape unit
      inum2 = ihout1                      ! inflow hyd=outlfow from hillslope
      inum3 = isub                        ! subbasin number
      rnum1 = 1.                          ! fraction overland flow
      iru_sub = 1                         ! route across landscape unit
      !! compute weighted K factor for sediment transport capacity
      sumk = 0.
      ovsl = 0.
      ovs = 0.
      do j = 1, hrutot(isub)
        sumk = sumk + usle_k(j) * hru_rufr(inum1,j)
        ovsl = ovsl + slsubbsn(j)
        ovs = ovs + hru_slp(j)
      end do 
      ovsl = ovsl / hrutot(isub)
      ovs = ovs / hrutot(isub)
      ru_k(isub,inum1) = sumk
      ru_ovsl(isub,inum1) = ovsl
      ru_ovs(isub,inum1) = ovs
      ru_ktc(isub,inum1) = 50.
      ru_a(isub,inum1) = daru_km(isub,1) / ru_ovsl(isub,inum1)
      call routels(iru_sub)               ! route across valley bottom
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      inum3s(ihout) = inum3
      ihouts(ihout) = ihout
      
      !! add routed with valley bottom loading
      inum1 = ihout                       ! hyd from routed 
      inum2 = ihout - 1                   ! hyd from loading
      ihout = ihout + 1                   ! outflow hyd number
      call addh                           ! add hyd's
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout
      
      !! save landscape routed output in place of subbasin output for routing
      varoute(isub,:) = varoute(ihout,:)
      end if
      
 1000 format(4i10,a10)
      return
      end