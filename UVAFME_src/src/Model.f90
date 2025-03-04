module Model

!*******************************************************************************
  !
  ! This contains the five main subroutines of UVAFME:
  !
  ! BioGeoClimate: computes daily and yearly site- and plot-level weather and
  !                soil dynamics
  !
  ! Canopy:        computes the plot-level LAI and light availability
  !
  ! Growth:        computes annual tree growth and branch thinning
  !
  ! Mortality:     determines which trees die and adds their components to the
  !                soil
  !
  ! Renewal        updates the seed and seedling banks for each species and
  !                regenerates new trees
  !
!*******************************************************************************


    use Parameters
    use Constants
    use Soil
    use Site
    use Species
    use Tree
    use Random
    use Climate
    use Input
    use Output
    use Manage

    implicit none

    ! Data dictionary: constants
    real, parameter :: MIN_HT = 1.83 ! Minimum height for "canopy" fuels (m)

contains

    !:.........................................................................:

    subroutine BioGeoClimate(site, year)
        !
        !  Computes daily weather data and annaul sums of weather and soil
        !  characteristics
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong          Original Code
        !    01/01/12     K. Holcomb           Updated to OOP structure
        !    10/10/16     A. C. Foster         Updated for soil/plot overhaul
        !                                       and permafrost updates
        !    05/01/17     A. C. Foster         Updated for moss and nutrient
        !                                       updates
        !

        ! Data dictionary: constants
        real,    parameter :: MIN_GROW_TEMP = 5.0   ! Minimum temperature for growing season (degC)
        real,    parameter :: MAX_DRY_PARM = 1.0001 ! Threshold for below wilting point
        real,    parameter :: DRY_THRESH = 0.6      ! Threshold for droughty conditions
        real,    parameter :: MIN_FLOOD_PARM = 0.8  ! Threshold for flooded conditions
        real,    parameter :: MIN_MOIST_PARM = 0.6  ! Moisture index threshold
        real,    parameter :: FIRE_THRESH = 50.0    ! Threshold intensity for self-sustaining fire (kW/m)
        integer, parameter :: N_MEM = 3             ! Number of days to calculate temperature "memory"

        ! Last Julian Day in each month
        integer, dimension(12), parameter :: MODAYS = [31, 59, 90, 120, 151,   &
            181, 212, 243, 273, 304, 334, 365]


        ! Data dictionary: calling arguments
        integer,        intent(in)    :: year ! Calendar year
        type(SiteData), intent(inout) :: site ! Site object

        ! Data dictionary: local variables
        real,    dimension(NTEMPS, 2)     :: tdd            ! Thawing degree-days (>0degC)
        real,    dimension(NTEMPS, 2)     :: fdd            ! Freezing degree-days (<0degC)
        real,    dimension(NTEMPS)        :: tmin           ! Monthly minimum temperature (degC)
        real,    dimension(NTEMPS)        :: tmax           ! Monthly maximum temperature (degC)
        real,    dimension(NTEMPS)        :: prcp           ! Monthly precipitation (cm)
        real,    dimension(NTEMPS)        :: tmean          ! Monthly average temperature (deg)
        real,    dimension(NTEMPS)        :: cld            ! Monthly cloudiness (tenths of sky covered)
        real,    dimension(NTEMPS)        :: rh             ! Monthly relative humidity (%)
        real,    dimension(NTEMPS)        :: strikes        ! Monthly lightning (strikes/km2/day)
        real,    dimension(NTEMPS)        :: tmptmin        ! Temporary variable for calculating actual tmin (degC)
        real,    dimension(NTEMPS)        :: tmptmax        ! Temporary variable for calculating actual tmax (degC)
        real,    dimension(NTEMPS)        :: tmpprec        ! Temporary variable for calculating actual prcp (cm)
        real,    dimension(NTEMPS)        :: tmpcld         ! Temporary variable for calculating actual cld (tenths of sky)
        real,    dimension(NTEMPS)        :: tmprh          ! Temporary variable for calculating actual rh (%)
        real,    dimension(NTEMPS)        :: tmpwind        ! Temporary variable for calculating actual wind (m/s)
        real,    dimension(NTEMPS)        :: tmpstrikes     ! Temporary variable for calculating actual
                                                            ! lightning (strikes/km2/day)
        real,    dimension(NTEMPS)        :: wind           ! Monthly wind speed (m/s)
        real,    dimension(NTEMPS)        :: solrad         ! Monthly solar radiation
        real,    dimension(NTEMPS)        :: soilW          ! Monthly soil moisture
        real,    dimension(NTEMPS)        :: soilI          ! Monthly soil ice content
        real,    dimension(NTEMPS)        :: alt            ! Monthly active layer thickness
        real,    dimension(DAYS_PER_YEAR) :: daywind        ! Daily wind speed (m/s)
        real,    dimension(DAYS_PER_YEAR) :: dayRH          ! Daily relative humidity (%)
        real,    dimension(DAYS_PER_YEAR) :: daystrikes     ! Daily relative lightning (strikes/km2/day)
        real,    dimension(DAYS_PER_YEAR) :: ffmc           ! Daily fine fuel moisture code
        real,    dimension(DAYS_PER_YEAR) :: dmc            ! Daily duff moisture code
        real,    dimension(DAYS_PER_YEAR) :: daytemp        ! Daily temperature (degC)
        real,    dimension(DAYS_PER_YEAR) :: daytemp_min    ! Daily minimum temperature (degC)
        real,    dimension(DAYS_PER_YEAR) :: daytemp_max    ! Daily maximum temperature (degC)
        real,    dimension(DAYS_PER_YEAR) :: daycld         ! Daily cloud cover (tenths of sky covered)
        real,    dimension(DAYS_PER_YEAR) :: dayprecip      ! Daily precipitation (cm)
        real,    dimension(DAYS_PER_YEAR) :: sun            ! Surface solar radiation (cal/cm2/day)
        real,    dimension(DAYS_PER_YEAR) :: st             ! Horizontal surface solar radiation (cal/cm2/day)
        real,    dimension(DAYS_PER_YEAR) :: exrad          ! Top of atmosphere solar radiation (cal/cm2/day)
        real,    dimension(DAYS_PER_YEAR) :: pot_ev_day     ! Potential evapotranspiration (cm)
        real,    dimension(DAYS_PER_YEAR) :: solar_rad_day  ! Solar radiation (cal/cm2/day)
        real,    dimension(DAYS_PER_YEAR) :: soilW_day      ! Daily soil liquid moisture (volumetric)
        real,    dimension(DAYS_PER_YEAR) :: soilI_day      ! Daily soil ice content (volumetric)
        real,    dimension(DAYS_PER_YEAR) :: alt_day        ! Daily active layer thickness (m)
        character(len = MAX_CHAR)         :: message        ! Error message
        real                              :: rain           ! Annual precipitation (cm)
        real                              :: rain_n         ! Annual N deposition (tN)
        real                              :: temp_f         ! Factor for creating temperature randomness
        real                              :: prcp_f         ! Factor for creating precipitation randomness
        real                              :: cld_f          ! Factor for creating cloud cover randomness
        real                              :: rh_f           ! Factor for creating relative humidity randomness
        real                              :: wnd_f          ! Factor for creating wind speed randomness
        real                              :: strike_f       ! Factor for creating lightning randomness
        real                              :: temp_max       ! Maximum temperature of warmest month (degC)
        real                              :: temp_min       ! Mininum temperature of warmest month (degC)
        real                              :: daytemp_mem    ! Average temperature over last N_MEM days
        real                              :: tmean_max      ! Maximum average temperature - for finding warmest month
        real                              :: n_avail        ! Plant-available nitrogen (tN/ha)
        real                              :: pet            ! Potential evapotranspiration (cm)
        real                              :: e1             ! Saturation vapor presstion at tmin of warmest month
        real                              :: e2             ! Saturation vapor presstion at tmax of warmest month
        real                              :: aet            ! Annual actual evapotranspiration (cm)
        real                              :: aet_mm         ! Annual actual evapotranspiration (mm)
        real                              :: growdays       ! Growing season length (days)
        real                              :: soildays       ! Soil degree-days (>0degC)
        real                              :: moistdays      ! Moisture index
        real                              :: wpdays         ! Proportion of growing season below wilting point
        real                              :: drydays        ! Proportion of growing season with drought conditions
        real                              :: flooddays      ! Proportion of growing season with flooded conditiosn
        real                              :: degday         ! Growing degree-days (>5degC)
        real                              :: atm_days       ! Proportion of growing season with atm. demand > supply
        real                              :: outwater       ! Runoff (cm)
        real                              :: tot_sun        ! Annual surface solar radiation (cal/cm2/day)
        real                              :: tot_st         ! Annual horizontal surface solar radiation (cal/cm2/day)
        real                              :: cfs            ! Ratio of surface:horizontal surface solar radiation
        real                              :: act_ev_day     ! Actual evapotranspiration (cm)
        real                              :: tcum           ! Cumulative thawing degree-days (>0degC)
        real                              :: fcum           ! Cumulative freezing degree-days (<0degC)
        real                              :: amlt           ! Last year's active layer depth (m)
        real                              :: xmlt           ! Thaw depth (m)
        real                              :: xfrz           ! Freezing depth (m)
        real                              :: zh             ! Soil layer depth
        real                              :: alff           ! Available light on the forest floor (0-1)
        real                              :: pc_germ        ! Effect of temperature on germination (0-1)
        real                              :: aow0_ByMin     ! Organic layer moisture scaled by wilting point
        real                              :: saw0_ByFC      ! Mineral layer moisture scaled by field capacity
        real                              :: saw0_BySAT     ! Mineral layer moisture scaled by saturation capacity
        real                              :: saw0_ByWP      ! Mineral layer moisture scaled by wilting point
        real                              :: saw0_ByFC_sum  ! Sum of mineral layer moisture scaled by field capacity
        real                              :: aow0_ByMin_sum ! Sum of organic layer moisture scaled by wilting point
        real                              :: saw0_BySAT_sum ! Sum of mineral layer moisture scaled by saturation capacity
        real                              :: tmpstep1       ! Temporary variable for implementing linear climate change
        real                              :: tmpstep2       ! Temporary variable for implementing linear climate change
        real                              :: tmp            ! Temporary variable for implementing linear climate change
        real                              :: FDI            ! Fire Danger Index (0-1)
        real, dimension(numplots)         :: max_FDI          ! Maximum annual Fire Danger Index (0-1)
        real, dimension(numplots)         :: cum_FDI          ! Cumulative annual Fire Danger Index (0-1)
        real        :: max_FDI_mn          ! Mean maximum annual Fire Danger Index (0-1)
        real        :: cum_FDI_mn          ! Mean maximum annual Fire Danger Index (0-1)
        real        :: max_FDI_sd          ! SD maximum annual Fire Danger Index (0-1)
        real        :: cum_FDI_sd          ! SD maximum annual Fire Danger Index (0-1)
        ! real, dimension(DAYS_PER_YEAR)         :: dailyFDI          ! Maximum annual Fire Danger Index (0-1)
        ! real, dimension(numplots)         :: mean_FDI          ! Maximum annual Fire Danger Index (0-1)
        ! real        :: mean_FDI_sd          ! SD maximum annual Fire Danger Index (0-1)
        real                              :: rosf           ! Rate of forward spread of fire (m/min)
        real                              :: I_surf         ! Surface fire intensity (kW/m)
        real                              :: tau_l          ! Residence time of fire (min)
        real                              :: ign_prob       ! Probability of an ignition event
        real                              :: rand           ! Random number (uniform)
        real                              :: tmp_sum        ! Temporary variable for calculating atm. demand
        real                              :: cadd           ! Available thawing days for soil
        real                              :: cft            ! Correction for soil thawing degree days
        real, dimension(:),  allocatable  :: I_r             ! Reaction intensity for all fires that occur
        real, dimension(:),  allocatable  :: I_r_temp       ! Holder variable
        real                              :: I_r_mn ! Mean reaction intensity for all fires that occur
        real                              :: I_r_sd ! SD reaction intensity for all fires that occur
        integer                           :: Nfire ! Number of fires per site per year
        integer                           :: siteid         ! Site ID
        integer                           :: warmest_month  ! Warmest month
        integer                           :: hrise          ! Hour of sunrise
        integer                           :: i, j, m, ip    ! Looping indices
        integer                           :: l
        logical                           :: strike         ! Do we have a strike event?


        ! Initialize accumulators
        rain = 0.0
        rain_n = 0.0
        tmean_max = RNVALID
        site%plots(:)%fire = 0
        max_FDI = 0
        cum_FDI = 0
        Nfire = 0

        ! Set site ID - in case we need to warn user
        siteid = site%site_id

        ! Check for and implement climate chnage
        ! The user is expected to input decr_by values as positive
        if (linear_cc) then
            ! Using linear climate change
            if (year >= start_gcm) then
                site%accum_tmin = site%accum_tmin + tmin_change
                site%accum_tmax = site%accum_tmax + tmax_change
                do m = 1, NTEMPS
                    tmpstep1 = site%precip(m) + site%accum_precip(m)
                    tmpstep2 = tmpstep1 * precip_change
                    site%accum_precip(m) = site%accum_precip(m) + tmpstep2
                end do
            endif
        else if (use_gcm) then
            ! Using climate change from input file - figure out which year
            ! we are in the file
            if (year >= start_gcm) then

                ! Read in climate change data
                call read_gcm_climate(site%site_id, year, start_gcm, tmin, &
                    tmax, prcp)

                call read_gcm_lightning(site%site_id, year, start_gcm,    &
                   strikes)

                if (site%site_id /= INVALID) then
                    ! Update climate values
                    site%tmin = tmin
                    site%tmax = tmax
                    site%precip = prcp*MM_TO_CM
                  !  site%rh = rh
                    site%strikes = strikes

                    ! Readjust for altitude if needed
                    if (adjust_altitude) then
                        call adjustForAltitude(site)
                    end if
                else
                    ! Problem reading in climate data - tell user, but don't
                    ! update
                    write(message, '(A, I6, A)') "Bad climate data for site ", &
                        siteid, " - using historical climate."
                    call warning(message)
                end if
            endif
        else if (use_monthly_clim) then
            call read_monthly_climate(site%site_id, year, tmean, tmin, tmax,   &
                prcp, rh, wind, solrad, soilW, soilI, alt)
            site%tmin = tmin
            site%tmax = tmax
            site%precip = prcp*mm_to_cm
            site%rh = rh
            site%wind = wind
            solrad = solrad/4.1868*60*60*24/10000
		endif

        ! Generate current year's weather from distributions of input climate
        ! data
        do i = 1, NTEMPS
            if (linear_cc) then
                ! Adjust for linear climate change
                tmptmin(i) = site%tmin(i) + site%accum_tmin
                tmptmax(i) = site%tmax(i) + site%accum_tmax
                tmpprec(i) = site%precip(i) + site%accum_precip(i)
                tmpcld(i) = site%cld(i)
                tmprh(i)   = site%rh(i)
                tmpwind(i) = site%wind(i)
                tmpstrikes(i) = site%strikes(i)
            else
                tmptmin(i) = site%tmin(i)
                tmptmax(i) = site%tmax(i)
                tmpprec(i) = site%precip(i)
                tmpcld(i)  = site%cld(i)
                tmprh(i)   = site%rh(i)
                tmpwind(i) = site%wind(i)
                tmpstrikes(i) = site%strikes(i)
            endif

            if (.not. use_monthly_clim) then

                ! Calculate climate fluctuations
                temp_f = clim_nrand(0.0, 1.0)
                prcp_f = clim_nrand(0.0, 1.0)
                cld_f = clim_nrand(0.0, 1.0)
                rh_f  = clim_nrand(0.0, 1.0)
                wnd_f = clim_nrand(0.0, 1.0)
                strike_f = clim_nrand(0.0, 1.0)

                ! Adjust
                prcp_f = max(-1.0, min(prcp_f, 1.0))
                temp_f = max(-1.0, min(temp_f, 1.0))
                cld_f = max(-1.0, min(cld_f, 1.0))
                rh_f = max(-1.0, min(rh_f, 1.0))
                wnd_f = max(-1.0, min(wnd_f, 1.0))
                strike_f = max(-1.0, min(strike_f, 1.0))

                ! Adjust monthly climate vars with std and random numbers
                if (use_gcm .and. year >= start_gcm ) then

                    ! Just use input values for that month for precip, temperature,
                    ! relative humidity, and lightning
                    tmin(i) = tmptmin(i)
                    tmax(i) = tmptmax(i)
                    prcp(i) = max(tmpprec(i), 0.0)
                    rh(i) = max(tmprh(i) + rh_f*site%rh_std(i), 0.0)
                    strikes(i) = max(strikes(i), 0.0)
                    wind(i) = max(tmpwind(i), 0.0)
                    cld(i) = max(tmpcld(i) + cld_f*site%cld_std(i), 0.0)
                else
                    ! Adjust for std
                    tmin(i) = tmptmin(i) + temp_f*site%tmin_std(i)
                    tmax(i)= tmptmax(i) + temp_f*site%tmax_std(i)

                    ! Can't be less than 0.0
                    cld(i) = max(tmpcld(i) + cld_f*site%cld_std(i), 0.0)
                    prcp(i) = max(tmpprec(i) + prcp_f*site%precip_std(i), 0.0)
                    rh(i) = max(tmprh(i) + rh_f*site%rh_std(i), 0.0)
                    wind(i) = max(tmpwind(i), 0.0)
                    strikes(i) = max(tmpstrikes(i) +                           &
                        strike_f*site%strikes_std(i), 0.0)
                end if

                ! Get mean monthly temperature for warmest month calculation
                ! tmean(i) = (site%tmin(i) + site%tmax(i))/2.0
                tmean(i) = (tmin(i) + tmax(i))/2.0
            else

                cld_f = clim_nrand(0.0, 1.0)
                strike_f = clim_nrand(0.0, 1.0)

                cld_f = max(-1.0, min(cld_f, 1.0))
                strike_f = max(-1.0, min(strike_f, 1.0))

                ! Just use input values for that month
                tmin(i) = tmptmin(i)
                tmax(i) = tmptmax(i)
                prcp(i) = max(tmpprec(i), 0.0)
                rh(i) = max(tmprh(i), 0.0)
                wind(i) = max(tmpwind(i), 0.0)

                cld(i) = max(tmpcld(i) + cld_f*site%cld_std(i), 0.0)
                strikes(i) = max(tmpstrikes(i) + strike_f*site%strikes_std(i), &
                    0.0)

            end if

            ! Accumulate precipitation and N deposition
            rain = rain + prcp(i)
            rain_n = rain_n + prcp(i)*PRCP_N

        end do

        ! Find warmest month of this year
        do i = 1, NTEMPS
            tmean_max = max(tmean_max, tmean(i))
            if (tmean_max == tmean(i)) warmest_month = i
        end do

        ! Get tmax and tmin of warmest month
        temp_max = site%tmax(warmest_month)
        temp_min = site%tmin(warmest_month)

        ! Calculate e2 and e1 (used for PET calculation)
        e1 = esat(temp_min)
        e2 = esat(temp_max)

        ! Convert monthly weather data into daily weather data
        call cov365_state(tmin, daytemp_min)
        call cov365_state(tmax, daytemp_max)
        call cov365_integr(prcp, dayprecip)
        call cov365_state(cld, daycld)
        call cov365_state(rh, dayRH)
        call cov365_state(wind, daywind)
        call cov365_state(strikes, daystrikes)

        if (use_monthly_clim) then
            call cov365_state(tmean, daytemp)
            call cov365_state(solrad, solar_rad_day)
            call cov365_state(soilW, soilW_day)
            call cov365_state(soilI, soilI_day)
            call cov365_state(alt, alt_day)
        end if

        if (fire_on .or. (site%prescribed_burn .and. year == site%year_management)) then

            ! Loop to calculate FFMC and DMC
            ! FFMC and DMC initialized at beginning of year
            ffmc(1) = 85.0
            dmc(1) = 6.0
            m = 1
            do j = 2, DAYS_PER_YEAR

                if (j > MODAYS(m)) m = m + 1 ! Get correct month

                ! Calculate FFMC
                call calc_ffmc(daytemp_max(j), dayRH(j), daywind(j),           &
                    dayprecip(j), ffmc(j-1), ffmc(j))

                ! Calculate DMC
                call calc_dmc(daytemp_max(j), dayRH(j), dayprecip(j),          &
                    dmc(j-1), m, site%latitude, dmc(j))

                ! Calculate daily fuel conditions, rate of spread, and fire
                ! intensity
            end do

            ! Loop to calculate fire conditions for each plot
            plotloop: do ip = 1, site%numplots

                ! Update fuels for fresh litter
                call update_fuels(site%plots(ip)%soil)

                ! Reset fire
                site%plots(ip)%fire = 0

                dayloop: do j = 1, DAYS_PER_YEAR

                    !check for testing conditions
                    if (fire_testing .and. year == site%fire_year .and.        &
                        j == site%fire_day) then
                        strike = .true.
                        ffmc(j) = site%fire_ffmc
                        dmc(j) = site%fire_dmc
                        daywind(j) = site%fire_wind
                    else if (.not. fire_testing .and.                          &
                        site%prescribed_burn .and.                             &
                        year == site%year_management .and.                     &
                         j == 180) then
                        ! j == site%fire_day) then
                        strike = .true.
                    !else if (.not. fire_testing .and.                          &
                    !        .not. site%prescribed_burn .and.                   &
                    !        year == 100 .and.                                  &
                    !        j == site%fire_day) then
                    !    strike = .true.
                    !    ffmc(j) = site%fire_ffmc
                    !    dmc(j) = site%fire_dmc
                    !    daywind(j) = site%fire_wind
                    else if (.not. fire_testing .and.                          &
                        .not. site%prescribed_burn .and.                       &
                        daystrikes(j) > 0.0)  then
                        strike = .true.
                    else
                        strike = .false.
                    end if

                    strikecheck: if (strike) then

                        ! Get fuel conditions for each day
                        call fuel_conditions(site%plots(ip)%soil,              &
                            ffmc(j), dmc(j), FDI, site%site_id,  &
                            ip, year, j)

                        max_FDI(ip) = max(max_FDI(ip), FDI)
                        cum_FDI(ip) = cum_FDI(ip) + FDI
                        ! daily_FDI(j) = FDI
                        ! call stddev(daily_FDI, mean_FDI(ip), mean_FDI_sd, 0)
                        ! Probility of an ignition event
                        if (fire_testing .or. site%prescribed_burn) then
                            ign_prob = 1.0
                        !else if (.not. fire_testing .and. year == 100) then
                        !    ign_prob = 1.0
                        else
                            ign_prob = daystrikes(j)*FDI*0.1
                        end if
                        rand = urand()

                        ! Check for ignition event
                        ignite: if (rand < ign_prob) then
                            ! Rate of spread (m/min)
                            call rate_of_spread(site%plots(ip)%soil,           &
                                daywind(j), ffmc(j), FDI, rosf, site%site_id,  &
                                    ip, year, j, I_r)

                            ! Surface fire intensity (kW/m)
                            call fire_intensity(site%plots(ip)%soil, rosf,     &
                                I_surf)

                            if (fire_testing .or. I_surf >= 50.0 .or.          &
                                site%prescribed_burn) then
                                ! Fire is enough to ignite and spread
                                site%plots(ip)%fire = 1
                                site%plots(ip)%fire_day = j
                          !      site%plots(ip)%soil%dmc_fire = site%fire_dmc
                                site%plots(ip)%soil%dmc_fire = dmc(j)

                                Nfire = Nfire + 1
                                exit dayloop
                            else 
                                ! Need to remove last I_r from list because fire didn't start
                                allocate(I_r_temp(size(I_r)))
                                I_r_temp = I_r
                                if(allocated(I_r)) deallocate(I_r)
                                if(size(I_r_temp) > 1) then 
                                    allocate(I_r(size(I_r_temp)-1))
                                    I_r = I_r_temp(1:size(I_r_temp) - 1)
                                end if
                                if(allocated(I_r_temp)) deallocate(I_r_temp)
                            end if
                       end if ignite
                    end if strikecheck
                end do dayloop
            end do plotloop
            if (mgmt_testing) then 
                call stddev(max_FDI(:), max_FDI_mn, max_FDI_sd, RNVALID)   
                call stddev(cum_FDI(:), cum_FDI_mn, cum_FDI_sd, RNVALID)   
                if(allocated(I_r)) then 
                    call stddev(I_r(:), I_r_mn, I_r_sd, RNVALID)   
                else 
                    I_r_mn = RNVALID
                    I_r_sd = RNVALID
                end if
                call csv_write(sitefire, site%site_id, .false.)
                call csv_write(sitefire, year, .false.)
                call csv_write(sitefire, max_FDI_mn, .false.)
                call csv_write(sitefire, cum_FDI_mn, .false.)
                call csv_write(sitefire, I_r_mn, .false.)
                call csv_write(sitefire, I_r_sd, .false.)
                call csv_write(sitefire, Nfire, .true.)
            end if 
        end if

        if(allocated(I_r)) deallocate (I_r)

        ! Initialize accumulators
        pet = 0.0
        tot_sun = 0.0
        tot_st = 0.0
        site%pc_germ = 0.0
        m = 1

        ! Calculate mean daily temperature, solar radiation, and PET
        do i = 1, DAYS_PER_YEAR

            if (.not. use_monthly_clim) then
                ! Mean daily temperature (degC)
                daytemp(i) = 0.5*(daytemp_min(i) + daytemp_max(i))
            end if

            ! Save these values for fire weather subroutines
            site%maxdaytemp = daytemp_max(i)
            site%dayprecip = dayprecip(i)

            ! Calculate temperature impact on regeneration
            if (daytemp(i) > MIN_GROW_TEMP) then ! During growing season

                ! Calculate temperature "memory"
                daytemp_mem = sum(daytemp(i-N_MEM+1:i))/float(N_MEM)

                ! Calculate effect on germination
                if (daytemp_mem <= 15.0) then
                    pc_germ = min(1.0, max(0.0, daytemp_mem*0.1 - 0.5))
                else
                    pc_germ = 1.0
                endif

                ! Set to max, so we only neet to hit it once
                site%pc_germ = max(site%pc_germ, pc_germ)
            end if

            ! Calculate solar radiation (cal/cm2/day)
            call solar_rad(i, site%latitude, site%slope, site%aspect,          &
                daycld(i), exrad(i), sun(i), st(i), hrise)

            if (use_monthly_clim) then
                sun(i) = solar_rad_day(i)
            end if

             ! Accumulate surface and horizontal surface radiation
             tot_sun = tot_sun + sun(i) ! Actual surface
             tot_st = tot_st + st(i)    ! Horizontal surface

            ! Calculate PET (cm)
            pot_ev_day(i) = pot_evap(daytemp(i), sun(i), site%altitude, e2, e1)

            ! Accumulate PET (cm)
            pet = pet + pot_ev_day(i)
        end do

        ! Calculate ratio of actual surface to horizontal surface radiation
        cfs = tot_sun/tot_st

        ! Calculate freezing and thawing degree days for permafrost subroutine
        tdd = 0.0
        fdd = 0.0
        m = 1
        do j = 1, DAYS_PER_YEAR
            if (j > MODAYS(m)) m = m + 1
            if (tmean(m) > epsilon(1.0) .and. daytemp(j) > epsilon(1.0)) then
                 tdd(m, 1) = tdd(m, 1) + daytemp(j)
            end if
            if (tmean(m) <= epsilon(1.0) .and. daytemp(j) <= epsilon(1.0)) then
                 fdd(m, 1) = fdd(m, 1) + abs(daytemp(j))
            end if
        end do

        ! Calculate cumulative freezing and thawing degree days
        tcum = 0.0
        fcum = 0.0
        do m = 12, 1, -1
             if (fdd(m, 1) > 0.0) fcum = fcum + fdd(m, 1)
             if (fdd(m, 1) == 0.0) exit
        end do
        do m = 1, 12
             if (tdd(m, 1) == 0.0) tcum = 0.0
             tcum = tcum + tdd(m, 1)
             tdd(m, 2) = tcum
             if (fdd(m, 1) == 0.0) fcum = 0.0
             fcum = fcum + fdd(m, 1)
             fdd(m, 2) = fcum
        end do

        ! Loop through each plot to calculate soil dynamics
        do ip = 1, site%numplots

            ! Initialize accumulators
            aet = 0.0
            degday = 0.0
            growdays = 0.0
            soildays = 0.0
            drydays = 0.0
            flooddays = 0.0
            moistdays = 0.0
            outwater = 0.0
            tmp_sum = 0.0
            wpdays = 0.0
            aow0_ByMin_sum = 0.0
            saw0_ByFC_sum = 0.0
            saw0_BySAT_sum = 0.0

            ! Store depth of thaw from previous year (m)
            amlt = min(site%plots(ip)%soil%active, site%plots(ip)%soil%A_depth)
            site%plots(ip)%amlt = amlt

            ! Reset
            site%plots(ip)%soil%active = 0.0
            site%plots(ip)%soil%z_freeze = 0.0

            ! Initialize depths freeze and thaw
            xmlt = 0.0
            xfrz = site%plots(ip)%soil%M_depth +                               &
                site%plots(ip)%soil%O_depth + site%plots(ip)%soil%A_depth

            ! Calculate light on the forest floor
            alff = 1.0*exp(-0.25*site%plots(ip)%cla/plotsize)

            do l = 1, 2
               ! Calculate drainage conditions
                site%plots(ip)%soil%z_drain(l) =                               &
                    (site%plots(ip)%soil%sat(l)*(1.0 - amlt) +                 &
                    site%plots(ip)%soil%fc(l)*(amlt - 0.32))/(1.0 - 0.32)

                ! Must be between field capacity and saturation capacity
                site%plots(ip)%soil%z_drain(l) =                               &
                    min(site%plots(ip)%soil%z_drain(l),                        &
                    site%plots(ip)%soil%sat(l))
                site%plots(ip)%soil%z_drain(l) =                               &
                    max(site%plots(ip)%soil%z_drain(l),                        &
                     site%plots(ip)%soil%fc(l))

                ! Set soil to fully saturated at onset of year
                if (l == 1) then
                     site%plots(ip)%soil%wc(l) =                               &
                        site%plots(ip)%soil%z_drain(l)*H2O_D/                  &
                         site%plots(ip)%soil%O_bulk_dens
                     zh = site%plots(ip)%soil%O_depth +                        &
                        site%plots(ip)%soil%M_depth
                else
                     site%plots(ip)%soil%wc(l) =                               &
                        site%plots(ip)%soil%z_drain(l)*H2O_D/                  &
                         site%plots(ip)%soil%A_bulk_dens
                     zh = site%plots(ip)%soil%A_depth
                end if

                ! Soil is completely frozen at start of year
                site%plots(ip)%soil%H2Oice(l) =                                &
                    site%plots(ip)%soil%z_drain(l)*zh
                 site%plots(ip)%soil%water(l) = 0.0
                 site%plots(ip)%soil%d_melt(l) = 0.0
                 site%plots(ip)%soil%d_freeze(l) = zh
            end do

            ! Loop on days - thaw depths are calculated monthly so have to
            ! increment the months as well
            m = 1
            do j = 1, DAYS_PER_YEAR

                if (j > MODAYS(m)) m = m + 1

                ! Calculate freeze/thaw depths (xfrz, xmlt) and maximum depths
                ! of freeze and thaw
                call permf(site%plots(ip)%soil, m, 1, alff, tdd, fdd, cfs, xfrz)
                call permf(site%plots(ip)%soil, m, 2, alff, tdd, fdd, cfs, xmlt)
                if (use_monthly_clim) then
                    xmlt = min(alt_day(j), 3.0)
                end if

                ! Update maximum depths (z_freeze and active)
                site%plots(ip)%soil%z_freeze = max((xfrz -                     &
                    site%plots(ip)%soil%M_depth -                              &
                    site%plots(ip)%soil%O_depth), site%plots(ip)%soil%z_freeze)
                site%plots(ip)%soil%active = max((xmlt -                       &
                    site%plots(ip)%soil%M_depth -                              &
                    site%plots(ip)%soil%O_depth), site%plots(ip)%soil%active)

                 ! Calculate soil water dynamics for the day
                if (use_monthly_clim) then
                    call moist(site%plots(ip)%soil, site%site_id, ip, year, j, &
                        daytemp(j), dayprecip(j), pot_ev_day(j),               &
                        site%leaf_area_ind, site%slope, amlt, xmlt, xfrz, tdd, &
                        m, act_ev_day, site%flow, aow0_ByMin, saw0_ByFC,       &
                        saw0_ByWP, saw0_BySAT, soilW_day(j), soilI_day(j))
                else
                    call moist(site%plots(ip)%soil, site%site_id, ip, year, j, &
                        daytemp(j), dayprecip(j), pot_ev_day(j),               &
                        site%leaf_area_ind, site%slope, amlt, xmlt, xfrz, tdd, &
                        m, act_ev_day, site%flow, aow0_ByMin, saw0_ByFC,       &
                        saw0_ByWP, saw0_BySAT)
                end if

                ! Accumulate variables
                outwater = outwater + site%plots(ip)%soil%runoff
                aet = act_ev_day + aet
                tmp = 1.0 - max(min(dayprecip(j)/pot_ev_day(j), 1.0),          &
                    min(act_ev_day/pot_ev_day(j), 1.0))

                ! Compute degday, dry days, flood days, and growing season
                ! length (days)
                if (daytemp(j) >= MIN_GROW_TEMP) then

                    ! Growing degree-days
                    degday = degday + (daytemp(j) - MIN_GROW_TEMP)

                    ! Growing season length
                    growdays = growdays + 1.0

                    ! For averageing values
                    saw0_ByFC_sum = saw0_ByFC_sum + saw0_ByFC
                    saw0_BySAT_sum = saw0_BySAT_sum + saw0_BySAT
                    aow0_ByMin_sum = aow0_ByMin_sum + aow0_ByMin

                    if (use_monthly_clim) then
                        if (saw0_ByWP < 1.01) then
                            drydays = drydays + 1.0
                        end if
                        if (saw0_BySAT > 0.6) then
                             flooddays = flooddays + 1.0
                             tmp = 0.0
                        endif
                    else
                        if (saw0_BySAT > MIN_FLOOD_PARM) then
                             flooddays = flooddays + 1.0
                             tmp = 0.0
                        else if (saw0_ByFC < DRY_THRESH) then
                             drydays = drydays + 1.0
                        end if

                    end if

                    if (aow0_ByMin < MAX_DRY_PARM) then
                         wpdays = wpdays + 1.0
                    end if

                    tmp_sum = tmp_sum + tmp
                end if

                ! Accumulate soil degree-days
                if (daytemp(j) >= 0.0) then
                    soildays = soildays + (daytemp(j) - 0.0)
                end if

                if (saw0_BySAT > MIN_MOIST_PARM) then
                     moistdays = moistdays + 1.0
                endif

            end do
            moistdays = moistdays/DAYS_PER_YEAR

            ! Convert drydays, flooddays, and wpdays to proportion of growing
            ! season
            if (growdays == 0) then
                 drydays = 0.0
                 flooddays = 0.0
                 wpdays = 0.0
                 atm_days = 0.0
            else
                 atm_days = tmp_sum/growdays
                 tmp = max(min(rain/pet, 1.0), min(aet/pet, 1.0))
                 !drydays = ((drydays/growdays) + atm_days)/2.0
                 drydays = ((drydays/growdays) + (1.0 - tmp))/2.0
                 flooddays = flooddays/growdays
                 wpdays = wpdays/growdays
            endif

            if (flooddays >= 0.1) then
                drydays = 0.0
            end if

            ! Convert aet to mm for decomposition
            aet_mm = aet*10.0
            
                ! if (year == 0 .and. init_on .and. &
                ! site%plots(ip)%soil%cohorts(1,2) == site%plots(ip)%soil%cohorts(1,1)*0.006) then
                !   ! this is an initialized black spruce site; let's initialize some moss!
                !   site%plots(ip)%soil%moss_biom = 300.0
                !   ! site%plots(ip)%soil%litter(IMOSS) = site%plots(ip)%soil%litter(IMOSS)/plotsize*HEC_TO_M2*KG_TO_T
                !   ! site%plots(ip)%soil%M_depth = site%plots(ip)%soil%moss_biom/plotsize/BULK_MOSS
                ! end if

            call moss(site%plots(ip)%soil, alff, site%plots(ip)%cla,           &
                site%plots(ip)%soil%dec_fuel, drydays, site%site_id, ip, year)


            call soiln(site%plots(ip)%soil, aet_mm, site%plots(ip)%cla,        &
                soildays, moistdays, n_avail)
            ! Set fan to 0.0 (used last year's fan for this year's soiln
            ! calculation)
            site%plots(ip)%soil%fan = 0.0

            ! Set plot-level attributes for year-end values
            site%plots(ip)%soil%avail_N = n_avail + rain_n
            site%plots(ip)%saw0_ByFC = saw0_ByFC_sum/growdays
            site%plots(ip)%aow0_ByMin = aow0_ByMin_sum/growdays
            site%plots(ip)%saw0_BySAT = saw0_BySAT_sum/growdays
            site%plots(ip)%soil%runoff = outwater
            site%plots(ip)%act_evap_day = aet
            site%plots(ip)%flood_days = flooddays
            site%plots(ip)%dry_days = drydays
            site%plots(ip)%wilt_days = wpdays
        end do

        ! Set site-level attributes to yearly sums of climate values
        site%deg_days = degday
        site%grow_days = growdays
        site%pot_evap_day = pet
        site%solar = tot_sun
        site%rain = rain

        ! Summer vapor pressure deficit
        site%summerVPD = ((6.1078 * 10.0**((7.5*tmean(7))/(237.3+tmean(7)))) - (rh(7)/100.0) + &
            (6.1078 * 10.0**((7.5*tmean(8))/(237.3+tmean(8)))) - (rh(8)/100.0))/2

        ! Reducing moss at initialization...
        ! if(year .eq. start_year) then
        !   do ip = 1, numplots
        !     site%plots(ip)%soil%M_depth = site%plots(ip)%soil%M_depth/4
        !   end do
        ! end if

    end subroutine BioGeoClimate

    !:.........................................................................:

    subroutine Canopy(site)
        !
        !  Calculates plot-level leaf area, LAI, and light environment
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong          Original Code
        !    01/01/12     K. Holcomb           Updated to OOP structure

        ! Data dictionary: constants
        real, parameter :: XT = -0.40      ! Light extinction coefficient
        real, parameter :: MIN_CBD = 0.011 ! Minimum canopy bulk density (kg/m3)

        ! Data dictionary: calling arguments
        type(SiteData), intent(inout) :: site ! Site object

        ! Data dictionary: local variables
        real, dimension(maxheight) :: la_dec      ! Leaf area experienced by deciduous plants (m2)
        real, dimension(maxheight) :: la_con      ! Leaf area experienced by evergreen plants (m2)
        real, dimension(maxheight) :: cla_dec     ! Cumulative leaf area experienced by deciduous plants (m2)
        real, dimension(maxheight) :: cla_con     ! Cumulative leaf area experienced by evergreen plants (m2)
        real, dimension(maxheight) :: cbiom       ! Canopy biomass (kg/m2)
        real, dimension(maxheight) :: treecanopy       ! Tree vertial LA (m2)
        real                       :: forht       ! Tree height (m)
        real                       :: canht       ! Clear branch bole height (m)
        real                       :: canopy_bh   ! Canopy base height
        real                       :: tla         ! Tree leaf area (m2)
        real                       :: tla_adj     ! Leaf area per 1-m height bin (m2)
        real                       :: cbiom_adj   ! Canopy biomass per 1-m height bin (kg/m2)
        real                       :: sum_bd      ! Sum of canopy bulk density (kg/m2)
        real                       :: max_ht      ! Maximum tree height (m)
        real                       :: canopy_biom ! Tree canopy biomass (kg)
        integer                    :: ntrees      ! Number of trees on plot
        integer                    :: iht         ! Tree height (m)
        integer                    :: cl          ! Canopy length (m)
        integer                    :: i, ip       ! Looping indices
        integer                    :: ih, it, is, ist
        integer                    :: m

        ! Initialize accumulators
        site%leaf_area_ind = 0.0
        site%lai_array = 0.0
        site%la_spec = 0.0

        ! Loop through plots to calculate leaf area of each tree and LAI of
        ! stand
        do ip = 1, site%numplots

            ! Get number of trees on plot
            ntrees = site%plots(ip)%numtrees

            if (ntrees == 0) then

                ! Full light conditions and no nutrient pressure
                site%plots(ip)%con_light = 1.0
                site%plots(ip)%dec_light = 1.0
                site%plots(ip)%fc_nutr = 1.0

                ! No canopy biomass
                site%plots(ip)%canopy_bd = 0.0
                site%plots(ip)%canopy_bh = 0.0
                sum_bd = 0.0
                max_ht = 0.0
                canopy_bh = 0.0
                canopy_biom = 0.0

            else

                ! Initialize leaf area and canopy biomass arrays
                la_dec = 0.0
                la_con = 0.0
                cla_dec = 0.0
                cla_con = 0.0
                cbiom = 0.0
                max_ht = 0.0
                sum_bd = 0.0
                canopy_bh = 0.0

                do it = 1, ntrees

                    if (site%plots(ip)%trees(it)%form <= 3) then

                        ! Tree or tall/erect shrub

                        ! Total tree height (m)
                        forht = max(site%plots(ip)%trees(it)%forska_ht, 1.0)

                        ! Integer of tree height (m)
                        iht = min(int(forht), maxheight)

                        ! Clear branch bole height (m)
                        canht = max(site%plots(ip)%trees(it)%canopy_ht, 1.0)

                    else

                        ! Erect shrub
                        forht = 1.0
                        iht = min(int(forht), maxheight)
                        canht = 1.0

                    end if

                    ! Get tallest tree/shrub
                    if (site%plots(ip)%trees(it)%forska_ht >= MIN_HT) then
                        max_ht = max(max_ht, forht)
                    end if

                    ! Calculate leaf area (m2)
                    tla = leaf_area(site%plots(ip)%trees(it))

                    ! Calculate 1-hr canopy biomass (kg)
                    if (site%plots(ip)%trees(it)%forska_ht >= MIN_HT) then
                        canopy_biom = (site%plots(ip)%trees(it)%leaf_bm +      &
                            site%plots(ip)%trees(it)%branchC*                  &
                            PERC_BRANCHES(1))*T_TO_KG/B_TO_C
                        else
                            canopy_biom = 0.0
                    end if

                    ! Accumulate site leaf area
                    site%leaf_area_ind = site%leaf_area_ind + tla

                    !Calculate canopy depth and divide leaf area and biomass
                    ! into 1m sections
                    cl = max(int(forht) - int(canht) + 1, 1)
                    tla_adj = tla/float(cl)
                    cbiom_adj = canopy_biom/float(cl)

                    ! Calculating (summer) vertical 1-m LA by species
                    ist = 0
                    is = 1
                    do is = 1, size(site%species)
                      if(site%species(is)%form == 1) then
                      ist = ist + 1
                      end if
                      if (site%plots(ip)%trees(it)%spec_ptr%unique_id ==   &
                       site%species(is)%unique_id) then
                       exit
                      end if
                    end do
                    if(site%plots(ip)%trees(it)%spec_ptr%form /= 1) then
                      ist = 0
                    end if
                    if (ist >=1) then
                      do ih = int(canht), int(forht)
                        site%la_spec(ist, ih) = site%la_spec(ist, ih) + tla_adj
                        ! site%la_spec(ist,:) = site%la_spec(ist,:) + treecanopy
                      end do
                    end if

                    ! Fill temporary arrays with leaf area/biomass
                    if (site%plots(ip)%trees(it)%conifer) then
                        ! Leaf experienced by evergreens reduced for deciduous
                        ! plants. This accounts for part of each year without
                        ! decid. leaves
                        do ih = int(canht), int(forht)
                            la_dec(ih) = la_dec(ih) + tla_adj
                            la_con(ih) = la_con(ih) + tla_adj
                            cbiom(ih) = cbiom(ih) + cbiom_adj
                        end do
                        ! la_dec = la_dec + treecanopy
                        ! la_con = la_con + treecanopy
                        ! cbiom = cbiom + canopy_biom/treecanopy
                    else
                        do ih = int(canht), int(forht)
                            la_dec(ih) = la_dec(ih) + tla_adj
                            la_con(ih) = la_con(ih) + tla_adj*0.8
                            cbiom(ih) = cbiom(ih) + cbiom_adj
                        end do
                        ! la_dec = la_dec + treecanopy
                        ! la_con = la_con + treecanopy*0.8
                        ! cbiom = cbiom + canopy_biom/treecanopy
                    end if

                end do

                ! Calculate cumulative leaf area from top down
                cla_dec(maxheight) = la_dec(maxheight)
                cla_con(maxheight) = la_con(maxheight)

                do ih = 1, maxheight - 1
                    ! Reduced leaf area from deciduous plants means higher
                    ! light environment for evergreens
                    ! Leaf area for deciduous plants normal b/c evergreen and
                    ! decid both present when decid has leaves
                    cla_dec(maxheight - ih) = cla_dec(maxheight - ih + 1) +    &
                        la_dec(maxheight - ih)
                    cla_con(maxheight - ih) = cla_con(maxheight - ih + 1) +    &
                        la_con(maxheight - ih)
                end do

                ! Calculate light environment
                do ih = 1, maxheight - 1
                    site%plots(ip)%con_light(ih) = exp(XT*cla_con(ih + 1)/     &
                        plotsize)
                    site%plots(ip)%dec_light(ih) = exp(XT*cla_dec(ih + 1)/     &
                        plotsize)
                end do
                site%plots(ip)%con_light(maxheight) = 1.0
                site%plots(ip)%dec_light(maxheight) = 1.0

                ! Convert canopy biomass to bulk density (kg/m3)
                cbiom(:) = cbiom(:)/plotsize

                ! Find canopy base height
                canopy_bh = 0.0
                do ih = 1, maxheight
                    if (cbiom(ih) >= MIN_CBD) then
                        canopy_bh = float(ih)
                        exit
                    end if
                end do

                if (canopy_bh >= 1.0) then
                    ! Sum up bulk density
                    sum_bd = sum(cbiom(int(canopy_bh):maxheight))
                else
                    sum_bd = 0.0
                end if

            end if ! end if any trees

            ! Save plot attributes
            site%plots(ip)%canopy_bd = sum_bd/(max_ht - canopy_bh)
            site%plots(ip)%canopy_bh = canopy_bh
            site%plots(ip)%canopy_biom = canopy_biom
            site%plots(ip)%cla = cla_dec(1)
            site%lai_array = site%lai_array + cla_dec

        end do !end plot loop

        ! Get average LAI (m2/m2) for site
        site%leaf_area_ind = site%leaf_area_ind/float(site%numplots)/plotsize
        site%lai_array = site%lai_array/float(site%numplots)/plotsize
        site%la_spec = site%la_spec/float(site%numplots)/plotsize

    end subroutine Canopy

    !:.........................................................................:

    subroutine Growth(site, year)
        !
        !  Calculates annual growth and branch thinning of each tree
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong          Original Code
        !    01/01/12     K. Holcomb           Updated to OOP structure
        !    10/10/16     A. C. Foster         Updated for soil/plot overhaul
        !                                       and permafrost updates
        !

        ! Data dictionary: constants
        integer, parameter :: MCOUNT = 2

        ! Data dictionary: calling arguments
        type(SiteData), intent(inout) :: site ! Site object
        integer,        intent(in)    :: year        ! Calendar year


        ! Data dictionary: local variables
        real,    dimension(size(site%species)) :: recr_trees  ! Number of reproductively active trees
        real,    dimension(maxcells*maxcells)  :: diam        ! Tree diameter (dbh)
        real,    dimension(maxcells*maxcells)  :: shade       ! Shade stress at top of tree's canopy (0-1)
        real,    dimension(maxcells*maxcells)  :: can_shade   ! Shade stress at bottom of tree's canopy (0-1)
        real,    dimension(maxcells*maxcells)  :: biomC       ! Woody biomass (tC)
        real,    dimension(maxcells*maxcells)  :: bleaf       ! Leaf biomass (tC)
        real                                   :: ht          ! Tree height (m)
        real                                   :: canht       ! Clear branch bole height(m)
        real                                   :: envstress   ! Growth stress factor (0-1)
        real                                   :: totbiomC    ! Total biomass on plot (tC)
        real                                   :: NPP         ! Net primary production
        real                                   :: d_leafb     ! Change in leaf biomass (tC)
        real                                   :: N_used      ! Nitrogen used by plant growth (tN)
        real                                   :: N_req       ! Nitrogen required for plant growth (tN)
        real                                   :: Npavail     ! Percent N available of required N
        real                                   :: dt          ! Diameter increment (cm)
        real                                   :: mindt       ! Minimum diameter increment before "stressed" (cm)
        real                                   :: dbiomC      ! Change in woody biomass from previous year (tC)
        real                                   :: leafbm      ! Leaf biomass (tC)
        real                                   :: bct         ! Aboveground woody biomass (tC)
        real                                   :: bcr         ! Branch biomass (tC)
        real                                   :: bcs         ! Stem biomass (tC)
        real                                   :: bcbr        ! Total branch biomass (tC)
        real                                   :: d_bc        ! Woody biomass lost to branch thinning (tC)
        real                                   :: d_bcs       ! Stem biomass lost to branch thinning (tC)
        real                                   :: d_bcr       ! Root biomass lost to branch thinning (tC)
        real                                   :: d_bcbr      ! Total branch biomass lost to branch thinning (tC)
        real                                   :: d_bctw      ! Twig biomass lost to branch thinning (tC)
        real                                   :: d_bcsb      ! Small branch biomass lost to branch thinning (tC)
        real                                   :: d_bclb      ! Large branch biomass lost to branch thinning (tC)
        integer                                :: hc          ! Canopy height (m)
        integer                                :: h           ! Tree height (m)
        integer                                :: ntrees      ! Number of trees on plot
        integer                                :: num_species ! Number of species in site
        integer                                :: it, is, ip  ! Looping indices
        integer                                :: lc          ! Litter class
        integer :: height !loop Index
        real :: tlight !Light available for individual tree
        logical :: conifer ! 0 is deciduous, 1 is conifer
        real              :: PAR    ! Annual photosynthetically active radiation at top of forest (MJ/m2/year)
        real,              dimension(maxheight) :: APAR_con   ! Annual available PAR by height bin (MJ/m2/year)
        real,              dimension(maxheight) :: APAR_dec    ! Annual available PAR by height bin (MJ/m2/year)
        real,              dimension(maxheight) :: treecanopy    ! Vertial LA distribution (m2)
        real :: forht, tla, tla_adj


        ! Get number of species at site
        num_species = size(site%species)

        ! Initialize recruiting trees accumulator
        recr_trees = 0.0

        plot: do ip = 1, site%numplots
            ! Initialize accumulators
            N_used = 0.0
            N_req = 0.0
            NPP = 0.0
            totbiomC = 0.0
            site%plots(ip)%mature(:) = 0
            site%plots(ip)%avail_spec = 0.0

            ! Calculate species-level response to drought, over-saturation,
            ! temperature, and permafrost
            do is = 1, num_species
                call temp_rsp(site%species(is), site%deg_days,                 &
                    site%plots(ip)%fc_gdd(is))
                call drought_rsp(site%species(is), site%plots(ip)%dry_days,    &
                    site%plots(ip)%fc_drought(is))
                call flood_rsp(site%species(is), site%plots(ip)%flood_days,    &
                    site%plots(ip)%fc_flood(is))
                call perm_rsp(site%species(is)%perm_tol,                       &
                    site%plots(ip)%soil%active, site%plots(ip)%fc_perm(is))
            end do

            ! Get number of trees
            ntrees = site%plots(ip)%numtrees

            numtrees: if (ntrees > 0) then

                stress: do it = 1, ntrees

                    ! Get species index and update tree
                    is = site%plots(ip)%trees(it)%species_index

                    ! Save diameter here
                    diam(it) = site%plots(ip)%trees(it)%diam_bht

                    ! Convenience variables to reduce table lookups
                    canht = site%plots(ip)%trees(it)%canopy_ht
                    ht = site%plots(ip)%trees(it)%forska_ht

                    ! Calculate if species is able to regenerate
                    ! site%plots(ip)%avail_spec(is) = max(kron(diam(it) -        &
                    !     site%species(is)%max_diam*site%species(is)%dbh_min),   &
                    !     site%plots(ip)%avail_spec(is))

                    site%plots(ip)%avail_spec(is) = max(kron(diam(it) -        &
                        site%species(is)%min_recr_dbh),   &
                        site%plots(ip)%avail_spec(is))

                    ! Calculate the number of reproductively active trees
                    if (site%plots(ip)%trees(it)%diam_bht >                    &
                        site%species(is)%min_recr_dbh) then ! .and.                    &
                        ! site%plots(ip)%trees(it)%tree_age >=                   &
                        ! site%species(is)%recr_age) then

                        site%plots(ip)%mature(is) =                            &
                            site%plots(ip)%mature(is) + 1
                        recr_trees(is) = recr_trees(is) + 1.0
                    end if

                    ! Get canopy and total tree height as integers
                    h = max(int(ht), 1)
                    hc = max(int(canht), 1)

                    ! Get leaf biomass and maximum possible DBH growth
                    call leaf_biomass_c(site%plots(ip)%trees(it))
                    call max_growth(site%plots(ip)%trees(it))

                    ! Save current value of leaf biomass
                    bleaf(it) = site%plots(ip)%trees(it)%leaf_bm

                    ! Calculate shading effect on tree
                    if (site%plots(ip)%trees(it)%conifer) then
                        shade(it) = light_rsp(site%species(is),                &
                            site%plots(ip)%con_light(h))
                        can_shade(it) = light_rsp(site%species(is),            &
                            site%plots(ip)%con_light(hc)) !can_shade isn't used in this version

                        ! Trying something new for calculating light availability:
                        tlight = 0.0
                        do height = hc, h
                          tlight = tlight + site%plots(ip)%con_light(height)
                        end do
                        if(h /= hc) then
                          tlight = tlight/(h-hc)
                        end if
                        shade(it) = light_rsp(site%species(is), tlight)
                    else
                        shade(it) = light_rsp(site%species(is),                &
                            site%plots(ip)%dec_light(h))
                        can_shade(it) = light_rsp(site%species(is),            &
                            site%plots(ip)%dec_light(hc))

                            tlight = 0.0
                            do height = hc, h
                              tlight = tlight + site%plots(ip)%dec_light(height)
                            end do
                            if(h /= hc) then
                              tlight = tlight/(h-hc)
                            end if
                            shade(it) = light_rsp(site%species(is), tlight)
                    end if

                    ! Compute new value
                    call biomass_c(site%plots(ip)%trees(it))
                    call biomass_n(site%plots(ip)%trees(it))

                    ! Calculate environmental stress (excluding nutrients)
                    call env_stress(site%plots(ip)%trees(it), shade(it),       &
                        site%plots(ip)%fc_gdd(is),                             &
                        site%plots(ip)%fc_drought(is),                         &
                        site%plots(ip)%fc_perm(is),                            &
                        site%plots(ip)%fc_flood(is), envstress)

                    ! Increment tree diameter using potential DBH growth
                    site%plots(ip)%trees(it)%diam_bht =                        &
                        diam(it) +                    &
                        max(0.0,site%plots(ip)%trees(it)%diam_opt*envstress)

                    ! Compute total height for new diameter
                    call forska_height(site%plots(ip)%trees(it))

                    ! Update clear branch bole height and leaf biomass with
                    ! new height
                    call stem_shape(site%plots(ip)%trees(it))
                    call leaf_biomass_c(site%plots(ip)%trees(it))

                    ! Calculate leaf and fine root growth N requirement
                    if (site%species(is)%conifer) then
                        N_req = N_req +                                        &
                            (site%plots(ip)%trees(it)%leaf_bm -                &
                            bleaf(it)*(1.0 - CON_LEAF_RATIO))/CON_LEAF_C_N
                    else
                        N_req = N_req +                                        &
                            site%plots(ip)%trees(it)%leaf_bm/DEC_LEAF_C_N
                    end if

                    ! Store old value of woody biomass
                    biomC(it) = site%plots(ip)%trees(it)%biomC +               &
                        site%plots(ip)%trees(it)%rootC

                    ! Compute new value
                    call biomass_c(site%plots(ip)%trees(it))
                    call biomass_n(site%plots(ip)%trees(it))

                    ! Calculate woody growth N requirement
                    N_req = N_req + ((site%plots(ip)%trees(it)%biomC +         &
                        site%plots(ip)%trees(it)%rootC) - biomC(it))/STEM_C_N

                end do stress

                ! Convert N_req tonnes N/ha
                N_req = max(N_req*HEC_TO_M2/plotsize, epsilon(1.0))

                ! Calculate percent available N
                Npavail = site%plots(ip)%soil%avail_N/N_req

                ! Calculate species-level response to available N
                do is = 1, num_species
                    site%plots(ip)%fc_nutr(is) = poor_soil_rsp(Npavail,        &
                        site%species(is)%lownutr_tol)
                end do

                ! Calculate actual DBH growth and N used
                grow: do it = 1, ntrees

                    ! Get species index
                    is = site%plots(ip)%trees(it)%species_index

                    ! Calculate environmental stress - including nutrients
                    call env_stress(site%plots(ip)%trees(it), shade(it),       &
                        site%plots(ip)%fc_gdd(is),                             &
                        site%plots(ip)%fc_drought(is),                         &
                        site%plots(ip)%fc_perm(is),                            &
                        site%plots(ip)%fc_flood(is),                           &
                        envstress,  site%plots(ip)%fc_nutr(is))

                    ! Calculate actual diameter increment growth
                    dt = max(0.0, site%plots(ip)%trees(it)%diam_opt*envstress)

                    ! Increment old diameter
                    site%plots(ip)%trees(it)%diam_bht = diam(it) + dt

                    ! Get check value for age and growth-related mortality
                    mindt = min(site%species(is)%max_diam/                     &
                        site%species(is)%max_age*0.1, site%species(is)%dbh_min)

                    ! Check for possible mortality age/growth stress mortality
                    if (dt <= site%species(is)%dbh_min) then

                        ! Diameter growth is below minimum level, increment
                        ! mortality counter
                        site%plots(ip)%trees(it)%mort_count =                  &
                            site%plots(ip)%trees(it)%mort_count + 1

                        if (site%plots(ip)%trees(it)%mort_count >= MCOUNT) then
                            ! Tree has been stressed for too many years,
                            ! turn on mortality flag
                            site%plots(ip)%trees(it)%mort_marker = .true.
                        ! else
                        !     ! Still possible to live
                        !     site%plots(ip)%trees(it)%mort_count = 0
                        endif
                    else
                        ! Rest mortality counter and set flag to false
                        site%plots(ip)%trees(it)%mort_count = 0
                        site%plots(ip)%trees(it)%mort_marker = .false.
                    endif

                    ! Compute actual new height and diameter
                    call forska_height(site%plots(ip)%trees(it))
                    call stem_shape(site%plots(ip)%trees(it))

                    ! Update biomass, saving leaf biomass into convenience
                    ! variable
                    call leaf_biomass_c(site%plots(ip)%trees(it))
                    call biomass_c(site%plots(ip)%trees(it))
                    call biomass_n(site%plots(ip)%trees(it))
                    leafbm = site%plots(ip)%trees(it)%leaf_bm

                    ! Calculate change in biomass from previous year
                    dbiomC = (site%plots(ip)%trees(it)%biomC +                 &
                        site%plots(ip)%trees(it)%rootC) - biomC(it)

                    ! Calculate C and N used
                    NPP = NPP + dbiomC
                    N_used = N_used + dbiomC/STEM_C_N

                    ! Update NPP and N_used from leaf biomass
                    if (site%plots(ip)%trees(it)%conifer) then

                        ! Conifers don't have to put all of their leaf biomass
                        ! back on
                        NPP = NPP + leafbm - bleaf(it)*(1.0 - CON_LEAF_RATIO)
                        N_used = N_used + (leafbm -                           &
                            bleaf(it)*(1.0 - CON_LEAF_RATIO))/CON_LEAF_C_N

                        ! Accumulate total biomass
                        totbiomC = totbiomC + site%plots(ip)%trees(it)%biomC + &
                            site%plots(ip)%trees(it)%rootC + leafbm
                    else

                        NPP = NPP + leafbm
                        N_used = N_used + leafbm/DEC_LEAF_C_N

                        ! Accumulate total woody biomass (no leaves)
                        totbiomC = totbiomC + site%plots(ip)%trees(it)%biomC + &
                            site%plots(ip)%trees(it)%rootC

                    end if

                    ! Calculate stand age as average tree age
                    site%plots(ip)%stand_age = site%plots(ip)%stand_age +      &
                        site%plots(ip)%trees(it)%tree_age

                    ! Get updated tree height and clear branch bole height
                    ht = site%plots(ip)%trees(it)%forska_ht
                    canht = site%plots(ip)%trees(it)%canopy_ht
                    hc = max(int(canht), 1)
                    conifer = site%plots(ip)%trees(it)%conifer

                    ! Check for lower branch thinning
                    ! This will increase clear branch bole height
                    ! branchfall: if (dt <= site%species(is)%dbh_min             &
                    !     .and. site%plots(ip)%trees(it)%form <= 2) then
                    ! branchfall: if (((conifer == .true. .and. site%plots(ip)%con_light(hc) < 0.2) .or. &
                    ! (conifer == .false. .and. site%plots(ip)%dec_light(hc) < 0.2)) &
                    !     .and. site%plots(ip)%trees(it)%form <= 2) then

                    branchfall: if (((conifer == .true. .and. light_rsp(site%plots(ip)%trees(it)%spec_ptr, &
                    site%plots(ip)%con_light(hc)) < 0.2) .or. &
                    (conifer == .false. .and. light_rsp(site%plots(ip)%trees(it)%spec_ptr, &
                    site%plots(ip)%dec_light(hc)) < 0.2)) &
                        .and. site%plots(ip)%trees(it)%form <= 2) then


                        ! Tree form and will drop some branches
                        ! Increment clear branch bole height
                        ! hc = hc + 1

                        ! Find height where there is enough light
                        do height = hc, ht
                          if(conifer == .true.) then
                            if(site%plots(ip)%con_light(height) >= 0.2) then
                              hc = height
                            exit
                            end if
                          else if (conifer == .false.) then
                            if(site%plots(ip)%dec_light(height) >= 0.2) then
                              hc = height
                              exit
                            end if
                          end if
                        end do

                        htcheck: if (hc < int(ht)) then

                            ! Clear branch bole height is still less than tree
                            ! height
                            ! So we can update it without any issues
                            site%plots(ip)%trees(it)%canopy_ht = float(hc) ! +  &
                              !  0.01

                            ! Update diameter at clear branch bole height
                            call stem_shape(site%plots(ip)%trees(it))

                            ! Save old woody biomass
                            bct = site%plots(ip)%trees(it)%biomC +             &
                                site%plots(ip)%trees(it)%rootC

                            ! Save old root, twig, stem C biomass
                            bcr = site%plots(ip)%trees(it)%rootC
                            bcbr = site%plots(ip)%trees(it)%branchC
                            bcs = site%plots(ip)%trees(it)%stemC

                            ! Update biomass C and N given new clear branch
                            ! bole height
                            call biomass_c(site%plots(ip)%trees(it))
                            call biomass_n(site%plots(ip)%trees(it))

                            ! How much wood litter did we lose?
                            d_bc = bct - (site%plots(ip)%trees(it)%biomC +     &
                                site%plots(ip)%trees(it)%rootC)
                            d_bcr = bcr - site%plots(ip)%trees(it)%rootC
                            d_bcbr = bcbr - site%plots(ip)%trees(it)%branchC
                            d_bcs = bcs -  site%plots(ip)%trees(it)%stemC

                            ! Divide branch biomass into twigs, large branches,
                            ! and small branches (Thonicke et al. 2010)
                            d_bctw = d_bcbr*PERC_BRANCHES(1)
                            d_bcsb = d_bcbr*PERC_BRANCHES(2)
                            d_bclb = d_bcbr*PERC_BRANCHES(3)

                            ! Add litter loss to litter pools
                            ! Here we convert to dry biomass because soiln
                            ! subroutine calculates weight loss not C loss

                            ! Roots
                            site%plots(ip)%soil%litter(IROOT) =                &
                                site%plots(ip)%soil%litter(IROOT) + d_bcr/B_TO_C

                            ! Twigs
                            site%plots(ip)%soil%litter(ITW) =                   &
                                site%plots(ip)%soil%litter(ITW) + d_bctw/B_TO_C

                            ! Small branches
                            site%plots(ip)%soil%litter(ISBR) =                 &
                                site%plots(ip)%soil%litter(ISBR) + d_bcsb/B_TO_C

                            !large branches
                            site%plots(ip)%soil%litter(ILBR) =                 &
                                site%plots(ip)%soil%litter(ILBR) + d_bclb/B_TO_C

                            ! Tree DBH < 10 cm go into smallwood, otherwise into
                            ! large wood
                            if (site%plots(ip)%trees(it)%diam_bht  >        &
                                10.0) then
                                ! If salvage is on, large boles don't enter litterpool
                                if(site%mgmtflag /= 3) then 
                                site%plots(ip)%soil%litter(ILBL) =             &
                                    site%plots(ip)%soil%litter(ILBL) +         &
                                    d_bcs/B_TO_C
                                end if 
                            else
                                site%plots(ip)%soil%litter(ISBL) =             &
                                    site%plots(ip)%soil%litter(ISBL) +         &
                                    d_bcs/B_TO_C
                            end if

                            ! Save previous value of leaf biomass
                            leafbm = site%plots(ip)%trees(it)%leaf_bm

                            ! Update leaf bm and get difference
                            call leaf_biomass_c(site%plots(ip)%trees(it))
                            d_leafb = leafbm - site%plots(ip)%trees(it)%leaf_bm

                            ! Add that litter to correct leaf litter class
                            lc = site%species(is)%litter_class
                            site%plots(ip)%soil%litter(lc) =                   &
                                site%plots(ip)%soil%litter(lc) + d_leafb/B_TO_C

                        end if htcheck

                    end if branchfall

                end do grow

                site%plots(ip)%stand_age = site%plots(ip)%stand_age/site%plots(ip)%numtrees

            end if numtrees

            ! Update plot-level soil characteristics
            N_used = N_used*HEC_TO_M2/plotsize ! tN/ha
            site%plots(ip)%NPP = NPP*HEC_TO_M2/plotsize ! tC/ha
            site%plots(ip)%soil%N_used = N_used
            site%plots(ip)%soil%avail_N = max(0.0,                             &
                site%plots(ip)%soil%avail_N - site%plots(ip)%soil%N_used)


        end do plot

        ! Update reproductively active trees to t/ha
        site%recr_trees = recr_trees/plotsize*HEC_TO_M2/numplots

        ! Update whether or not species is native
        ! Must have more than 5 trees/ha of repr. active trees
        do is = 1, num_species
            if (site%recr_trees(is) >= 5.0) then
                site%species(is)%native = .true.
              else
                site%species(is)%native = .false.
            end if
        end do
        !
        ! print *, PAR
        ! print *, APAR_con
        ! print *, APAR_dec

        call classify_site(site)
    end subroutine Growth

    !:.........................................................................:

    subroutine Mortality(site, year)
        !
        !  Determines which trees die by age, stress, or disturbances, and adds
        !  their biomass to the appropriate soil litter pools
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    01/01/12     K. Holcomb          Updated to OOP structure
        !    10/10/16     A. C. Foster        Updated for soil/plot overhaul
        !                                       and permafrost updates
        !    03/01/17     A. C. Foster        Updated for fire updates
        !    06/10/19     A. C. Foster        Updated for SPITFIRE
        !    10/10/20     A. C. Foster        Updated for harvest routines
        !

        ! Data dictionary: constants
        real, parameter :: SBR_RT = 0.5 ! Rate of consumption of live small branches
        real, parameter :: LBR_RT = 0.1  ! Rate of consumption of live large branches
        real, parameter :: BL_RT = 0.05  ! Rate of consumption of live boles

        ! Data dictionary: calling argumnets
        type(SiteData),        intent(inout) :: site ! Site object
        integer,               intent(in)    :: year ! Calendar year

        ! Data dictionary: local variables
        integer, dimension(:), allocatable :: dbh_ind          ! Index of sorted tree DBH array
        real,    dimension(:), allocatable :: dbh              ! Tree DBH (cm)
        real                               :: totbiomC         ! Plot-wide biomass (tC)
        real                               :: NPP_loss         ! Loss of NPP from mortality
        real                               :: fan              ! N volatilized by fires (tN)
        real                               :: consRoot         ! Proportion roots consumed by fire (0-1)
        real                               :: wind_prob        ! Random number for windthrow
        real                               :: biomC            ! Total tree biomass (tC)
        real                               :: leaf_bm          ! Leaf biomass (tC)
        real                               :: bcr              ! Root biomass (tC)
        real                               :: N_cons           ! Proportion of N consumed by fire (0-1)
        real                               :: burn             ! Amount of tree burned (tC)
        real                               :: bctw             ! Twig biomass (tC)
        real                               :: bcs              ! Stem biomass (tC)
        real                               :: bcbr             ! Total branch biomass (tC)
        real                               :: bcsb             ! Small branch biomass (tC)
        real                               :: bclb             ! Large branch biomass (tC)
        real                               :: ab_combust       ! Aboveground combustion (tC)
        real                               :: bgr_combust      ! Root combustion (tC)
        real                               :: not_burn         ! Didn't burn (tC)
        real                               :: R_a              ! Critical active rate of spread (m/min)
        real                               :: rosf_active      ! Active rate of spread (m/min)
        real                               :: CFB              ! Crown fraction burnt (0-1)
        real                               :: R_final          ! Final rate of spread (surface + crown) (m/min)
        real                               :: I_final          ! Final fire intensity (surface + crown) (kW/m)
        real, dimension(14)                :: pactvmgmt        ! Probability of active management at the site
        real                               :: tpactvmgmt        ! Probability of active management at the site
        integer                            :: num_species      ! Number of species on site
        integer                            :: it, ip, iu           ! Looping indices
        integer                            :: dt               ! Counter for dead trees
        integer                            :: lt               ! Counter for live trees
        integer                            :: ind              ! Tree array index
        integer                            :: thin_num         ! Number of trees to thin
        integer                            :: is               ! Species index
        integer                            :: trow             ! Row of tree
        integer                            :: tcol             ! Column of tree
        logical                            :: active_crowning  ! Active crown fire?
        logical                            :: passive_crowning ! Passive crown fire?
        logical                            :: age_survive      ! Does tree survive age check?
        logical                            :: growth_survive   ! Does tree survive stress check?
        logical                            :: newplant_survive ! Is the tree a sad little nursery seedling?
        logical                            :: fire_survive     ! Does tree survive fire?
        logical                            :: wind_survive     ! Does tree survive windthrow?
        logical                            :: cutting_survive  ! Does tree survive selective harvest?
        logical                            :: sel_cut          ! Are we selective cutting?
        logical                            :: newplant         ! Less than 5 years since  planting?
        integer                            :: snum             ! Tree growth stressor
        integer                            :: lc               ! Litter class
        integer                            :: i      ! Looping index

        ! Get number of species on site
        num_species = size(site%species)

        site%buldozed = .false.

        if(active_management) then
            call write_management_data(site, year, mgmt)
          active_treatments: if(site%mgmtflag == 1) then !harvest
            print *, "Harvesting", site%site_id
              call Harvest(site) ! add argument for cut threshold
          else if (site%mgmtflag == 2) then active_treatments ! thin 30%
            print *, "Thinning", site%site_id
            call Thin(site, 0.30)
            site%mgmtflag = 0
          else if (site%mgmtflag == 3) then active_treatments ! salvage
            print *, "Salvage", site%site_id
            ! This "salvages" a site after a significant mortality event. 
            ! Large boles for dead trees are removed. The rest goes into litter/ dead pools 
            ! like normal mortality. This process plays out within Mortality() below
            ! site%mgmtflag = 0
          else if (site%mgmtflag == 4) then active_treatments ! shearblade
            print *, "Shearblading", site%site_id
            call Shearblade(site)
            site%mgmtflag = 0
          else if (site%mgmtflag == 5) then active_treatments ! prune
            print *, "Pruning", site%site_id
            call Prune(site)
            site%mgmtflag = 0

          end if active_treatments

        else
            ! if(site%sel_cutting .and. site%planting /= "" .and. sel_cut == .false. &
            ! .and. ((year > site%year_management) .and. (mod((year - site%yedayfirear_management),     &
            ! site%rotation_time) <= 5))) then
            ! ! runs selective cut at management year and at rotations afterward
            !     newplant = .true.
            ! end if
              treatments: if (site%management .and. site%thinning .and.          &
                  year == site%year_management) then
                  ! We are thinning the stand (DOD)
                  site%mgmtflag = 2
                  call Thin(site, site%thin_perc)

              else if (site%management .and. site%shearblading .and.             &
                  year == site%year_management) then treatments
    
                  ! We are shearblading
                  site%mgmtflag = 4
                  call Shearblade(site)

                ! else if (site%management .and. site%sel_cutting .and. (year == site%year_management                   &
                !     .or. (year > site%year_management .and. mod((year - site%year_management), site%rotation_time == 0)))) then treatments
                else if (site%management .and. site%sel_cutting .and. (year == site%year_management)) then treatments
                site%mgmtflag = 1
                call Harvest(site)


              else if(site%management .and. site%pruning .and. &
                year == site%year_management) then treatments 
                ! We are pruning
                site%mgmtflag = 5
                call Prune(site)
    
              end if treatments
    
              call write_management_data(site, year, mgmt)
            
              !A lazy way to test planting after shearblading. 
              ! Eventually I may want to replace mgmtflag as the mechanism to initiate planting
              if(site%mgmtflag == 4 .or. site%mgmtflag == 1) then 
                  site%mgmtflag = 1
              else 
                  site%mgmtflag = 0
              end if 
          end if
  

            ! if (site%already_cut) then
            !     if (site%management .and. site%sel_cutting .and.               &
            !         site%plots(ip)%stand_age >= float(site%rotation_time)) then
            !         sel_cut = .true.
            !         print *, "stand age", site%plots(ip)%stand_age, "rotation time", site%rotation_time
            !     end if
            ! else
            !     if (site%management .and. site%sel_cutting .and.               &
            !         site%plots(ip)%stand_age >= float(site%rotation_time)) then
            !         site%already_cut = .true.
            !         sel_cut = .true.
            !
            !     end if
            ! end if

          plot: do ip = 1, site%numplots

            ! Initialize accumulators
            site%plots(ip)%num_dead = 0
            totbiomC = 0.0
            NPP_loss = 0.0
            fan = 0.0
            site%plots(ip)%d_type = 0.0
            active_crowning = .false.
            passive_crowning = .false.
            site%plots(ip)%soil%shrubLitter = 0.0
            ab_combust = 0.0
            bgr_combust = 0.0

            ! Set wind to 0
            site%plots(ip)%wind = 0

            ! Get random number for wind throw
            wind_prob = urand()

            if ((fire_on .or. site%prescribed_burn) .and. site%plots(ip)%fire == 1) then

                ! Check for active or passive crown fire
                call active_passive_check(site%plots(ip)%soil,                 &
                    site%plots(ip)%canopy_bh, site%plots(ip)%canopy_bd,        &
                    site%plots(ip)%canopy_biom, site%plots(ip)%soil%wind_fire, &
                    site%plots(ip)%soil%ffmc_fire, R_a, rosf_active, CFB,      &
                    R_final, I_final, passive_crowning, active_crowning)

                ! Calculate fuel consumption
                call fuel_consumption(site%site_id, ip, year,                  &
                    site%plots(ip)%soil, consRoot, N_cons)

                fntrees: if (site%plots(ip)%numtrees > 0) then

                    ! Kill trees that died by fire, age, or low growth - only
                    ! copy surviving trees, rest go into soil or burn

                    ! Initialize dead and live tree counters
                    dt = 0
                    lt = 0

                    fireloop: do it = 1, site%plots(ip)%numtrees

                        ! Get species index
                        is = site%plots(ip)%trees(it)%species_index

                        ! Get leaf biomass
                        call leaf_biomass_c(site%plots(ip)%trees(it))
                        leaf_bm = site%plots(ip)%trees(it)%leaf_bm

                        ! Check for growth and age survival
                        call growth_survival(site%plots(ip)%trees(it),         &
                            growth_survive)
                        call age_survival(site%plots(ip)%trees(it),            &
                            age_survive)

                        ! Check for fire survival
                        call fire_survival(site%plots(ip)%trees(it),           &
                            site%plots(ip)%soil%I_surf,                        &
                            site%plots(ip)%soil%tau_l, active_crowning,        &
                            CFB, fire_survive)

                        fdeathcheck: if (growth_survive .and.                  &
                            age_survive .and. fire_survive) then

                            ! Tree survives

                            lt = lt + 1 ! Increment live tree counter

                            ! Increment tree age
                            site%plots(ip)%trees(it)%tree_age =                &
                            site%plots(ip)%trees(it)%tree_age + 1

                            ! Copy tree to front of list
                            call copy_tree(site%plots(ip)%trees(lt),           &
                                site%plots(ip)%trees(it))

                            if (site%plots(ip)%trees(it)%form > 0 .and.        &
                                site%plots(ip)%trees(it)%forska_ht < MIN_HT) then
                                ! Add up shrub live foliage and fine twigs
                                site%plots(ip)%soil%shrubLitter =              &
                                    site%plots(ip)%soil%shrubLitter + leaf_bm + &
                                    site%plots(ip)%trees(it)%branchC*         &
                                    PERC_BRANCHES(1)

                            end if

                            ! Calculate leaf litter
                            if (site%species(is)%conifer) then
                                ! Conifers only drop some of their needles

                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm*CON_LEAF_RATIO/B_TO_C

                                ! Accumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm*CON_LEAF_RATIO

                            else

                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm/B_TO_C

                                ! Acumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm

                            end if

                            if (site%plots(ip)%trees(it)%CK > 0.0) then

                                ! Tree survived but was still damaged

                                !TODO - add in damage to live trees
                                burn = site%plots(ip)%trees(it)%leaf_bm*       &
                                    site%plots(ip)%trees(it)%CK +              &
                                    site%plots(ip)%trees(it)%branchC*          &
                                    PERC_BRANCHES(1)*                          &
                                    site%plots(ip)%trees(it)%CK +              &
                                    site%plots(ip)%trees(it)%branchC*          &
                                    PERC_BRANCHES(2)*                          &
                                    site%plots(ip)%trees(it)%CK*SBR_RT +       &
                                    site%plots(ip)%trees(it)%branchC*          &
                                    PERC_BRANCHES(3)*                          &
                                    site%plots(ip)%trees(it)%CK*LBR_RT +       &
                                    site%plots(ip)%trees(it)%stemC*            &
                                    site%plots(ip)%trees(it)%CK*BL_RT

                                ab_combust = ab_combust + burn
                                not_burn = not_burn +                          &
                                    max((site%plots(ip)%trees(it)%biomC +      &
                                    site%plots(ip)%trees(it)%leaf_bm) - burn,  &
                                    0.0)

                            end if

                        else fdeathcheck

                            ! Tree dies from something (need to check which)

                            dt = dt + 1 ! Increment dead tree counter

                            ! Copy to dead tree list for output
                            call copy_tree(site%plots(ip)%deadtrees(dt),       &
                                site%plots(ip)%trees(it))

                            ! Set cells of that tree to unfilled
                            trow = site%plots(ip)%trees(it)%row
                            tcol = site%plots(ip)%trees(it)%col
                            site%plots(ip)%cells(trow, tcol) = 0

                            ! Roots - fire consumption from Bonan (1989)
                            bcr = site%plots(ip)%trees(it)%rootC
                            burn = bcr*consRoot
                            bgr_combust = bgr_combust + burn
                            bcr = bcr - burn
                            site%plots(ip)%soil%litter(IROOT) =                &
                                site%plots(ip)%soil%litter(IROOT) + bcr/B_TO_C

                            ! Accumulate volatilized N
                            fan = fan + burn*litter_params(IROOT, 2)*          &
                                (1.0 - N_cons)

                            ! Branch biomass
                            bcbr = site%plots(ip)%trees(it)%branchC

                            ! Convert branch litter into twigs, small
                            ! branches, and large branches (Thonicke et al. 2010)
                            bctw = bcbr*PERC_BRANCHES(1)
                            bcsb = bcbr*PERC_BRANCHES(2)
                            bclb = bcbr*PERC_BRANCHES(3)

                            ! Twigs
                            burn = bctw*(site%plots(ip)%trees(it)%CK)
                            ab_combust = ab_combust + burn
                            not_burn = not_burn + bctw - burn
                            bctw = bctw - burn
                            site%plots(ip)%soil%litter(ITW) =                  &
                                site%plots(ip)%soil%litter(ITW) + bctw/B_TO_C

                            ! Accumulate volatilized N
                            fan = fan + burn*litter_params(ITW, 2)*            &
                                (1.0 - N_cons)

                            ! Small branches
                            burn = bcsb*(SBR_RT*site%plots(ip)%trees(it)%CK)
                            ab_combust = ab_combust + burn
                            not_burn = not_burn + bcsb - burn
                            bcsb = bcsb - burn
                            site%plots(ip)%soil%litter(ISBR) =                 &
                                site%plots(ip)%soil%litter(ISBR) + bcsb/B_TO_C

                            ! Accumulate volatilized N
                            fan = fan + burn*litter_params(ISBR, 2)*           &
                                (1.0 - N_cons)

                            ! Large branches
                            burn = bclb*(LBR_RT*site%plots(ip)%trees(it)%CK)
                            ab_combust = ab_combust + burn
                            not_burn = not_burn + bclb - burn
                            bclb = bclb - burn
                            site%plots(ip)%soil%litter(ILBR) =                 &
                                site%plots(ip)%soil%litter(ILBR) + bclb/B_TO_C

                            ! Accumulate volatilized N
                            fan = fan + burn*litter_params(ILBR, 2)*           &
                                (1.0 - N_cons)

                            ! Stems
                            bcs = site%plots(ip)%trees(it)%stemC
                            burn = bcs*(BL_RT*site%plots(ip)%trees(it)%CK)
                            ab_combust = ab_combust + burn
                            not_burn = not_burn + bcs - burn
                            bcs = bcs - burn

                            ! Small boles (DBH < 10) vs. large boles
                            if (site%plots(ip)%trees(it)%diam_bht > 10.0) then
                                ! If salvage is on, large boles don't enter litterpool
                                if(site%mgmtflag /= 3) then 
                                    site%plots(ip)%soil%litter(ILBL) =             &
                                    site%plots(ip)%soil%litter(ILBL) + bcs/B_TO_C
                                ! Accumulate volatilized N
                                fan = fan + burn*litter_params(ILBL, 2)*       &
                                    (1.0 - N_cons)
                                end if 
                            else
                                site%plots(ip)%soil%litter(ISBL) =             &
                                    site%plots(ip)%soil%litter(ISBL) +         &
                                    bcs/B_TO_C

                                ! Accumulate volatilized N
                                fan = fan + burn*litter_params(ISBL, 2)*       &
                                    (1.0 - N_cons)
                            end if

                            ! Leaves
                            lc = site%species(is)%litter_class
                            burn = leaf_bm*(site%plots(ip)%trees(it)%CK)
                            ab_combust = ab_combust + burn
                            not_burn = not_burn + leaf_bm - burn
                            leaf_bm = leaf_bm - burn

                            site%plots(ip)%soil%litter(lc) =                   &
                            site%plots(ip)%soil%litter(lc) + leaf_bm/B_TO_C

                            ! Accumulate volatilized N
                            fan = fan + burn*litter_params(lc, 2)*             &
                                (1.0 - N_cons)

                            ! Get most limiting growth factor
                            snum = site%plots(ip)%trees(it)%stressor

                            firecheck: if (.not. fire_survive) then

                                ! Died by fire, fire consumes some of each
                                !litter class. also calculate N volatilized in
                                ! tree burning

                                ! Add biomass to array of dead biomass
                                site%plots(ip)%d_type(IFIRE) =                 &
                                    site%plots(ip)%d_type(IFIRE) +             &
                                    site%plots(ip)%trees(it)%biomC + leaf_bm

                                ! Update "stressor"
                                site%plots(ip)%deadtrees(dt)%stressor = IFIRE

                                ! Acumulate NPP losses
                                NPP_loss = NPP_loss + biomC + leaf_bm

                            else if (.not. growth_survive .or.                 &
                                .not. age_survive) then firecheck

                                ! Died from growth or age-related stress, all
                                ! litter goes into soil

                                ! Add biomass to array of dead biomass
                                site%plots(ip)%d_type(snum) =                  &
                                    site%plots(ip)%d_type(snum) +              &
                                    site%plots(ip)%trees(it)%biomC + leaf_bm


                                ! Acumulate NPP losses
                                NPP_loss = NPP_loss + biomC + leaf_bm

                            end if firecheck

                        end if  fdeathcheck

                    end do fireloop

                    ! Set number of live and dead trees
                    site%plots(ip)%numtrees = lt
                    site%plots(ip)%num_dead = dt

                end if fntrees

                ! Convert to kg/m2
                ab_combust = ab_combust/B_TO_C*T_TO_KG/plotsize
                bgr_combust = bgr_combust/B_TO_C*T_TO_KG/plotsize
                not_burn = not_burn/B_TO_C*T_TO_KG/plotsize

                if (testing) then
                    call csv_write(cons_out, not_burn, .false.)
                    call csv_write(cons_out, site%plots(ip)%canopy_bd, .false.)
                    call csv_write(cons_out, site%plots(ip)%canopy_bh, .false.)
                    call csv_write(cons_out, site%plots(ip)%canopy_biom, .false.)
                    call csv_write(cons_out, R_a, .false.)
                    call csv_write(cons_out, rosf_active, .false.)
                    call csv_write(cons_out, CFB, .false.)
                    call csv_write(cons_out, R_final, .false.)
                    call csv_write(cons_out, I_final, .false.)
                    call csv_write(cons_out, ab_combust, .false.)
                    call csv_write(cons_out, bgr_combust, .true.)
                end if

            else if (wind_prob < site%wind_prob) then

                ! We have a windthrow event

                ! Set fire to 0 and wind to 1
                site%plots(ip)%fire = 0
                site%plots(ip)%wind = 1

                ! Set wind counter for regeneration
                site%plots(ip)%windCount = 3

                wntrees: if (site%plots(ip)%numtrees > 0) then

                    ! Initialize counters for live and dead trees
                    lt = 0
                    dt = 0

                    windloop: do it = 1, site%plots(ip)%numtrees

                        ! Get species index and update tree
                        is = site%plots(ip)%trees(it)%species_index

                        ! Get leaf biomass
                        call leaf_biomass_c(site%plots(ip)%trees(it))
                        leaf_bm = site%plots(ip)%trees(it)%leaf_bm

                        ! Check for growth and age mortality
                        call growth_survival(site%plots(ip)%trees(it),         &
                            growth_survive)
                        call age_survival(site%plots(ip)%trees(it),            &
                            age_survive)

                        ! Check for wind survival
                        call wind_survival(site%plots(ip)%trees(it),           &
                            wind_survive)

                        wdeathcheck: if (growth_survive .and.                  &
                            age_survive .and. wind_survive) then

                            ! Tree survives

                            lt = lt + 1 ! Increment live tree counter

                            ! Increment tree age
                            site%plots(ip)%trees(it)%tree_age =                &
                                site%plots(ip)%trees(it)%tree_age + 1

                            ! Copy tree object to top of list
                            call copy_tree(site%plots(ip)%trees(lt),           &
                                site%plots(ip)%trees(it))

                            if (site%plots(ip)%trees(it)%form > 0 .and.        &
                                site%plots(ip)%trees(it)%forska_ht < MIN_HT) then
                                ! Add up shrub live foliage and fine twigs
                                site%plots(ip)%soil%shrubLitter =              &
                                    site%plots(ip)%soil%shrubLitter + leaf_bm + &
                                    site%plots(ip)%trees(it)%branchC*         &
                                    PERC_BRANCHES(1)

                            end if

                            ! Calculate leaf litter
                            if (site%species(is)%conifer) then

                                ! Evergreens only lose some of their leaves
                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm*CON_LEAF_RATIO/B_TO_C

                                ! Acumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm*CON_LEAF_RATIO

                            else

                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm/B_TO_C

                                ! Accumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm
                            end if

                        else wdeathcheck

                            ! Tree dies
                            dt = dt + 1 ! Increment dead tree counter

                            ! Copy to list of dead trees for output
                            call copy_tree(site%plots(ip)%deadtrees(dt),       &
                                site%plots(ip)%trees(it))

                            ! Set cells of that tree to unfilled
                            trow = site%plots(ip)%trees(it)%row
                            tcol = site%plots(ip)%trees(it)%col
                            site%plots(ip)%cells(trow, tcol) = 0

                            ! Get most limiting growth factor
                            snum = site%plots(ip)%trees(it)%stressor

                            windcheck: if (.not. growth_survive  .or.          &
                                .not. age_survive) then

                                ! Tree died from low growth or age

                                ! Add biomass to mortality array
                                site%plots(ip)%d_type(snum) =                  &
                                    site%plots(ip)%d_type(snum) +              &
                                    site%plots(ip)%trees(it)%biomC + leaf_bm

                            else windcheck

                                ! Died by windthrow

                                ! Add biomass to mortality array
                                site%plots(ip)%d_type(IWIND) =                 &
                                    site%plots(ip)%d_type(IWIND) +             &
                                    site%plots(ip)%trees(it)%biomC + leaf_bm

                                site%plots(ip)%deadtrees(dt)%stressor = IWIND

                            end if windcheck

                            ! Add total tree litter components to all
                            ! litter categories

                            ! Roots
                            bcr = site%plots(ip)%trees(it)%rootC
                            site%plots(ip)%soil%litter(IROOT) =                &
                                site%plots(ip)%soil%litter(IROOT) + bcr/B_TO_C

                            ! Branches
                            bcbr = site%plots(ip)%trees(it)%branchC
                            ! Convert branch litter into twigs, small
                            ! branches, and large branches
                            bctw = bcbr*PERC_BRANCHES(1)
                            bcsb = bcbr*PERC_BRANCHES(2)
                            bclb = bcbr*PERC_BRANCHES(3)

                            ! Twigs
                            site%plots(ip)%soil%litter(ITW) =                  &
                                site%plots(ip)%soil%litter(ITW) + bctw/B_TO_C

                            ! Small branches
                            site%plots(ip)%soil%litter(ISBR) =                 &
                                site%plots(ip)%soil%litter(ISBR) + bcsb/B_TO_C

                            ! Large branches
                            site%plots(ip)%soil%litter(ILBR) =                 &
                                site%plots(ip)%soil%litter(ILBR) + bclb/B_TO_C

                            ! Stems
                            bcs = site%plots(ip)%trees(it)%stemC

                            ! Small boles (DBH < 10) vs. large boles
                            if (site%plots(ip)%trees(it)%diam_bht > 10.0) then
                                ! If salvage is on, large boles don't enter litterpool
                                if(site%mgmtflag /= 3) then 
                                site%plots(ip)%soil%litter(ILBL) =             &
                                    site%plots(ip)%soil%litter(ILBL) +         &
                                    bcs/B_TO_C
                                end if 
                            else
                                site%plots(ip)%soil%litter(ISBL) =             &
                                    site%plots(ip)%soil%litter(ISBL) +         &
                                    bcs/B_TO_C
                            end if

                            ! Leaves
                            lc = site%species(is)%litter_class
                            site%plots(ip)%soil%litter(lc) =                   &
                                site%plots(ip)%soil%litter(lc) + leaf_bm/B_TO_C

                            ! Acumulate NPP losses
                            NPP_loss = NPP_loss +                              &
                                site%plots(ip)%trees(it)%biomC  + leaf_bm

                        end if wdeathcheck
                    end do windloop

                        ! Set number of live and trees
                        site%plots(ip)%numtrees = lt
                        site%plots(ip)%num_dead = dt

                end if wntrees

            else

                ! No disturbances - just check for age and growth stress
                ! mortality

                ! Set fire and wind to 0
                site%plots(ip)%fire = 0
                site%plots(ip)%wind = 0

                gntrees: if (site%plots(ip)%numtrees > 0) then

                    ! Initialize counters for live and dead trees
                    lt = 0
                    dt = 0

                    deathloop: do it = 1, site%plots(ip)%numtrees

                        ! Get species index
                        is = site%plots(ip)%trees(it)%species_index


                        ! Get leaf biomass
                        call leaf_biomass_c(site%plots(ip)%trees(it))
                        leaf_bm = site%plots(ip)%trees(it)%leaf_bm


                        ! Check for age and growth survival
                        call growth_survival(site%plots(ip)%trees(it),         &
                            growth_survive)

                        call age_survival(site%plots(ip)%trees(it),            &
                            age_survive)

                        newplant_survive = .true.
                        if(newplant .and.                                      &
                        site%plots(ip)%trees(it)%planted == 1) then
                          call newplant_survival(site%plots(ip)%trees(it),     &
                            site%viability, newplant_survive)
                        end if

                        deathcheck: if (growth_survive .and. age_survive .and. &
                         newplant_survive) then

                            ! Tree survives

                            lt = lt + 1 ! Increment live tree counter

                            ! Increment tree age
                            site%plots(ip)%trees(it)%tree_age=                 &
                                site%plots(ip)%trees(it)%tree_age + 1

                            ! Copy tree to top of list
                            call copy_tree(site%plots(ip)%trees(lt),           &
                                site%plots(ip)%trees(it))

                            if (site%plots(ip)%trees(it)%form > 0 .and.        &
                                site%plots(ip)%trees(it)%forska_ht < MIN_HT) then
                                ! Add up shrub live foliage and fine twigs
                                site%plots(ip)%soil%shrubLitter =              &
                                    site%plots(ip)%soil%shrubLitter + leaf_bm + &
                                    site%plots(ip)%trees(it)%branchC*         &
                                    PERC_BRANCHES(1)

                            end if

                            ! Calculate litterfall
                            if (site%species(is)%conifer) then

                                ! Conifers only drop some of their leaves
                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm*CON_LEAF_RATIO/B_TO_C

                                ! Accumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm*CON_LEAF_RATIO

                            else

                                lc = site%species(is)%litter_class
                                site%plots(ip)%soil%litter(lc) =               &
                                    site%plots(ip)%soil%litter(lc) +           &
                                    leaf_bm/B_TO_C

                                ! Accumulate NPP losses
                                NPP_loss = NPP_loss + leaf_bm

                            end if

                        else deathcheck

                            ! Tree dies

                            dt = dt + 1

                            ! Copy to list of dead trees for output
                            call copy_tree(site%plots(ip)%deadtrees(dt),       &
                                site%plots(ip)%trees(it))

                            ! Set cells of that tree to unfilled
                            trow = site%plots(ip)%trees(it)%row
                            tcol = site%plots(ip)%trees(it)%col
                            site%plots(ip)%cells(trow, tcol) = 0


                            ! Get most limiting growth factor
                            snum = site%plots(ip)%trees(it)%stressor

                            ! Add biomass to mortality array
                            site%plots(ip)%d_type(snum) =                      &
                                site%plots(ip)%d_type(snum) +                  &
                                site%plots(ip)%trees(it)%biomC + leaf_bm

                            ! Add total tree litter components to all
                            ! litter categories

                            ! Roots
                            bcr = site%plots(ip)%trees(it)%rootC
                            site%plots(ip)%soil%litter(IROOT) =                &
                                site%plots(ip)%soil%litter(IROOT) + bcr/B_TO_C

                            ! Branches
                            bcbr = site%plots(ip)%trees(it)%branchC
                            ! Convert branch litter into twigs, small
                            ! branches, and large branches
                            bctw = bcbr*PERC_BRANCHES(1)
                            bcsb = bcbr*PERC_BRANCHES(2)
                            bclb = bcbr*PERC_BRANCHES(3)

                            ! Twigs
                            site%plots(ip)%soil%litter(ITW) =                  &
                                site%plots(ip)%soil%litter(ITW) + bctw/B_TO_C

                            ! Small branches
                            site%plots(ip)%soil%litter(ISBR) =                 &
                                site%plots(ip)%soil%litter(ISBR) + bcsb/B_TO_C

                            ! Large branches
                            site%plots(ip)%soil%litter(ILBR) =                 &
                                site%plots(ip)%soil%litter(ILBR) + bclb/B_TO_C

                            ! Stems
                            bcs = site%plots(ip)%trees(it)%stemC

                            ! Small boles (DBH < 10) vs. large boles
                            if (site%plots(ip)%trees(it)%diam_bht > 10.0) then
                                ! If salvage is on, large boles don't enter litterpool
                                if(site%mgmtflag /= 3) then 
                                site%plots(ip)%soil%litter(ILBL) =             &
                                    site%plots(ip)%soil%litter(ILBL) +         &
                                    bcs/B_TO_C
                                end if 
                            else
                                site%plots(ip)%soil%litter(ISBL) =             &
                                    site%plots(ip)%soil%litter(ISBL) +         &
                                    bcs/B_TO_C
                            end if

                            ! Leaves
                            lc = site%species(is)%litter_class
                            site%plots(ip)%soil%litter(lc) =                   &
                                site%plots(ip)%soil%litter(lc) + leaf_bm/B_TO_C

                            ! Acumulate NPP losses
                            NPP_loss = NPP_loss +                              &
                                site%plots(ip)%trees(it)%biomC  + leaf_bm

                        end if deathcheck

                    end do deathloop

                    ! Set number of live and dead trees
                    site%plots(ip)%num_dead = dt
                    site%plots(ip)%numtrees = lt

                end if gntrees
            end if

            ! Update N volatilized N (tN/ha)
            site%plots(ip)%soil%fan = fan/plotsize/M2_TO_HEC

            ! Update NPP (tC/ha)
            site%plots(ip)%NPP = site%plots(ip)%NPP -                          &
                NPP_loss/plotsize/M2_TO_HEC

        end do plot

        if (site%mgmtflag == 3) site%mgmtflag = 0 !Done with salvage

    end subroutine Mortality

    !:.........................................................................:

    subroutine Renewal(site, year, tree_id)
        !
        !  Calculates the seedling and seed banks of each species and
        !  regenerates new trees and shrubs
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    01/01/12     K. Holcomb          Updated to OOP structure
        !    10/10/16     A. C. Foster        Updated for soil/plot overhaul
        !                                       and permafrost updates
        !    06/10/20     A. C. Foster        Updated for shrubs and seed rain
        !

        ! Data dictionary: calling arguments
        type(SiteData), intent(inout) :: site ! Site object
        integer,        intent(in)    :: year ! Calendar year 
        integer,        intent(inout) :: tree_id !Unique tree counter

        ! Data dictionary: constants
        real, parameter :: REGMIN = 0.01 ! Minimum threshold for regeneration

        ! Data dictionary: local variables
        real,    dimension(size(site%species)) :: regstress   ! Regeneration stress factor (0-1)
        real, dimension(numplots, size(site%species)) :: regstress_all ! Regeneration stress factor (0-1) among all plots
        real, dimension(size(site%species)) :: regstress_mn ! Regeneration stress factor (0-1) mean for output
        real, dimension(size(site%species)) :: regstress_sd ! Regeneration stress factor (0-1) sd for output 
        real,    dimension(size(site%species)) :: recr_trees  ! Number of reproductively active trees
        real,    dimension(site%num_trees)     :: probt       ! Probability of regeneration (trees)
        real,    dimension(site%num_shrubs)    :: probs       ! Probability of regeneration (shrubs)
        integer, dimension(site%num_trees)     :: tree_sp     ! Id locations of tree species
        integer, dimension(site%num_shrubs)    :: shrub_sp    ! Id locations of shrub species
        integer, dimension(maxcells*maxcells)  :: r_empty     ! Rows of empty cells
        integer, dimension(maxcells*maxcells)  :: c_empty     ! Columns of empty cells
        integer, dimension(:), allocatable     :: locs        ! Locations of empty cells
        real                                   :: NPP         ! Net primary productivity (tC/ha)
        real                                   :: N_used      ! N used in regeneration (tN/ha)
        real                                   :: tregmax     ! Maximum growth stress factor for trees
        real                                   :: sregmax     ! Maximum growth stress fractor for shrubs
        real                                   :: fc_temp     ! Effect of temperature on germination
        real                                   :: fc_org      ! Effect of organic layer depth on regeneration
        real                                   :: germinants  ! Number of seedlings regenerating (#/m2)
        real                                   :: germinants1  ! Number of seedlings regenerating (#/m2)
        real                                   :: probtsum    ! Sum of probability of regeneration for trees
        real                                   :: probssum    ! Sum of probability of regeneration for shrubs
        real                                   :: rand        ! Random number for determining species (uniform)
        real                                   :: dbh         ! Tree dbh (cm)
        real                                   :: leafbm      ! Leaf biomass (tC)
        real                                   :: shade       ! Shade stress (0-1)
        real                                   :: envstress   ! Growth stress factor (0-1)
        real                                   :: layers !How much layering is allowed?
        integer                                :: num_species ! Number of species at site
        integer                                :: new_trees   ! Number of trees regenerated
        integer                                :: new_shrubs  ! Number of shrubs regenerated
        integer                                :: numtrees    ! Number of trees on plot
        integer                                :: numshrubs   ! Number of shrubs on plot
        integer                                :: max_trenew  ! Max number of trees to renew
        integer                                :: max_srenew  ! Max number of shrubs to renew
        integer                                :: ntrenew     ! Number of trees to renew
        integer                                :: nsrenew     ! Number of shrubs to renew
        integer                                :: org_tol     ! Ability to regenerate on deep soil
        integer                                :: n_empty     ! Number of empty cells on plot
        integer                                :: ht          ! Tree height (m)
        integer                                :: lc          ! Litter class
        integer                                :: is, ip, it  ! Looping indices
        integer                                :: stp, ssp    ! Species counters
        integer                                :: t, r, c     ! Looping indices
        integer                                :: irenew, i   ! Looping indices

        ! Get number of species at site
        num_species = size(site%species)

        plot: do ip = 1, site%numplots
        ! Initialize accumulators 
        new_trees = 0
        new_shrubs = 0
        NPP = 0.0
        N_used = 0.0
        tregmax = 0.0
        sregmax = 0.0

        ! Get number of trees and shrubs
        numtrees = count(site%plots(ip)%trees(:)%form == 1)
        numshrubs = count(site%plots(ip)%trees(:)%form > 1)


            ! Check to make sure we aren't still waiting on wind counter
            windcheck: if (site%plots(ip)%windCount == 0) then
                ! Calculate species-level environmental stressors for
                ! the plot
                do is = 1, num_species

                    ! First get minimum of gdd, nutrient, and drought
                    ! stress
                    regstress(is) = min(site%plots(ip)%fc_gdd(is),         &
                        site%plots(ip)%fc_nutr(is),                        &
                        site%plots(ip)%fc_drought(is))

                    if (site%plots(ip)%numtrees == 0) then
                        ! We don't have to worry about shade stress
                        regstress(is) = regstress(is)*                     &
                            site%plots(ip)%fc_flood(is)*                   &
                            site%plots(ip)%fc_perm(is)
                            ! if (is==1 .and. ip < 5) write(sitef, *) ("All the trees are dead"), site%site_id, year
                    else
                        ! Need to consider shading
                        if (site%species(is)%conifer) then

                            regstress(is) = min(regstress(is),             &
                                light_rsp(site%species(is),                &
                                site%plots(ip)%con_light(1)))*             &
                                site%plots(ip)%fc_flood(is)*               &
                                site%plots(ip)%fc_perm(is)
                        else
                            regstress(is) = min(regstress(is),             &
                                light_rsp(site%species(is),                &
                                site%plots(ip)%dec_light(1)))*             &
                                site%plots(ip)%fc_flood(is)*               &
                                site%plots(ip)%fc_perm(is)

                        end if
                    end if

                    ! print *, site%plots(ip)%con_light(1)

                    ! Check for enough mature trees
                    if (site%plots(ip)%mature(is) <= 5 .and. site%species(is)%conifer) then
                        regstress(is) = regstress(is)*0.5
                        !site%plots(ip)%avail_spec(is) = 0.0
                    end if

                    ! Can't regenerate if below minimum
                    if (site%species(is)%unique_id == 'PICEmari') then
                        if (regstress(is) < 0.01) then
                            regstress(is) = 0.0
                        end if
                    else
                        if (regstress(is) < REGMIN) then
                            regstress(is) = 0.0
                        end if
                    end if

                    regstress_all(ip,:) = regstress

                    ! Calculate maximum regrowth capacity across all species
                    if (site%species(is)%form == 1) then
                        tregmax = max(tregmax, regstress(is))
                    else
                        sregmax = max(sregmax, regstress(is))
                    end if

                end do

                ! Compute the max renew number trees and shrubs
                    max_trenew = max(min(int(maxtrees*tregmax) - numtrees,     &
                        maxtrees), 0)
                max_srenew = max(min(int(maxshrubs*sregmax) - numshrubs,   &
                    maxshrubs), 0)

                ! Compute actual number of renewable trees and shrubs
                ! if(ip==1) print *, numtrees, "numtrees"
                ntrenew = min(max_trenew, maxtrees - numtrees)
                nsrenew = min(max_srenew, maxshrubs - numshrubs)

            seedbank: do is = 1, num_species

            ! Get species-level regeneration response to fire
            call fire_rsp(site%species(is), site%plots(ip)%fire,   &
                site%plots(ip)%fc_fire(is))

            ! Check for enough mature trees
            if (site%plots(ip)%mature(is) <= 5) then
                site%plots(ip)%fc_fire(is) = min(1.0,            &
                    site%plots(ip)%fc_fire(is))
            end if
            
            if (site%species(is)%native) then
                ! 'native' species start with/or have achieved the
                ! number, size, and age for regenration
                site%plots(ip)%seedbank(is) =                      &
                    site%plots(ip)%seedbank(is) +                  &
                (site%species(is)%invader +            &
                    site%species(is)%seed_num*            &
                    site%plots(ip)%avail_spec(is)*                 &
                    site%plots(ip)%fc_fire(is))
            else if(seedrain) then
                ! Only get some seedrain from adjacent sites
                site%plots(ip)%seedbank(is) =                      &
                    site%plots(ip)%seedbank(is) +                  &
                    site%species(is)%fc_recr
            else
            ! "Native" has different meaning with init_stand()
            site%plots(ip)%seedbank(is) =                      &
                site%plots(ip)%seedbank(is) +                  &
                (site%species(is)%invader)
            end if
            ! if(init_on == .false. .and. year - start_year < 3) then !boosting invaders at start of spinup
            !   site%plots(ip)%seedbank(is) = site%plots(ip)%seedbank(is) + &
            !   site%species(is)%invader*10
            ! end if

            ! if(ip==1) print *, is, site%plots(ip)%seedbank(is)*plotsize, "seedbank per plot"

            ! We don't allow seedling regeneration the first
            ! year of a fire
            firecheck: if (site%plots(ip)%fire == 0) then


                ! Put seeds into seedling bank if envstress is
                ! high enough
             !  seedcheck: if (regstress(is) > site%species(is)%dbh_min*5) then
                seedcheck: if (regstress(is) > 0.0) then

                    ! Black spruce limited by temperature
                    if (site%species(is)%unique_id == 'PICEmari') then
                        fc_temp = site%pc_germ
                    else
                        fc_temp = 1.0
                    end if
                    if (fc_temp > 1.0) fc_temp = 1.0
                    if (fc_temp < 0.0) fc_temp = 0.0

                    ! Ability to reproduce on moss-covered soil
                    org_tol = site%species(is)%org_tol
                    ! fc_org = 1
                    fc_org = max(1.0,0.1 + exp(ORG_GF(org_tol)*            &
                        (site%plots(ip)%soil%O_depth +             &
                        site%plots(ip)%soil%M_depth)))

                    ! Calculate the number of new seedlings
                    germinants = site%plots(ip)%seedbank(is)*      &
                        fc_org*fc_temp

                    ! Remove new seedlings from seedbank
                    site%plots(ip)%seedbank(is) =                  &
                        site%plots(ip)%seedbank(is) - germinants

                    ! Add germinants to seedling bank
                    site%plots(ip)%seedling(is) =                  &
                        site%plots(ip)%seedling(is) + germinants + 1! +10 ! this was getting scaled *500...

                    ! Decrease seedbank for survival
                    site%plots(ip)%seedbank(is) =                  &
                        max(0.0, site%plots(ip)%seedbank(is)*               &
                        site%species(is)%seed_surv)

                else seedcheck

                    ! No new seedlings from seedbank

                    ! Decrease seedbank for survival
                    site%plots(ip)%seedbank(is) =                  &
                        max(0.0, site%plots(ip)%seedbank(is)*               &
                        site%species(is)%seed_surv)

                endif seedcheck

                ! Add seedlings from layering
                if (site%species(is)%layering) then
                    ! if ((site%plots(ip)%soil%M_depth +             &
                    !     site%plots(ip)%soil%O_depth) > 0.1) then
                    !     site%plots(ip)%seedling(is) =              &
                    !         min(site%plots(ip)%seedling(is)*5/plotsize,1000.0/plotsize)
                    ! end if
                    if ((site%plots(ip)%soil%M_depth +             &
                        site%plots(ip)%soil%O_depth) > 0.1) then
                        layers = min(2.0, 20.0*(site%plots(ip)%soil%M_depth +             &
                            site%plots(ip)%soil%O_depth))
                        site%plots(ip)%seedling(is) = &
                        min(site%plots(ip)%seedling(is)*layers,100.0)
                end if
                end if

                ! Add seedlings from sprouts
                site%plots(ip)%seedling(is) =                      &
                    site%plots(ip)%seedling(is) +                  &
                    site%species(is)%sprout_num*                   &
                    site%plots(ip)%avail_spec(is)

            end if firecheck
            ! Convert seedling bank to #/plot
            site%plots(ip)%seedling(is) =                      &
            site%plots(ip)%seedling(is)*plotsize
        end do seedbank

          ! if(init_on .and. year - start_year < 5) then
          !   if(ip==1) print *, "no trees added in first 5 years when initialization is on"
          !   exit plot
          ! end if
           
                    ! if(ip==1) print *, year, ntrenew, "number of max new trees"
                    ! Update the seed bank size and seedling bank size (#/m2)


                    ! Calculate probability of regeneration
                    ! Need to calculate shrub and tree probabilities
                    ! separately
                    probtsum = 0.0
                    probssum = 0.0
                    stp = 1
                    ssp = 1
                    do is = 1, num_species
                        if (site%species(is)%form == 1) then
                          probt(stp) = site%plots(ip)%seedling(is)*          &
                                regstress(is)

                            probtsum = probtsum + probt(stp)
                            tree_sp(stp) = is
                            stp = stp + 1
                        else
                            probs(ssp) = site%plots(ip)%seedling(is)*          &
                                regstress(is)
                            probssum = probssum + probs(ssp)
                            shrub_sp(ssp) = is
                            ssp = ssp + 1
                        end if
                    end do

                else windcheck

                    ! We are still waiting after a windthrow event
                    ! Decrease counter
                    site%plots(ip)%windCount = max(0,                          &
                        site%plots(ip)%windCount - 1)

                    ! Set probsums to 0
                    probtsum = 0.0
                    probssum = 0.0

                end if windcheck

                ! After setting seed and seedling banks
                ! Calculate cumulative probability of regeneration for trees
                if (probtsum > epsilon(1.0)) then
                    do is = 1, site%num_trees
                        probt(is) = probt(is)/probtsum
                    end do
                    do is = 2, site%num_trees
                        probt(is) = probt(is - 1) + probt(is)
                    end do
                else
                    ntrenew = 0
                end if

                ! Calculate cumulative probability of regeneration for shrubs
                if (probssum > epsilon(1.0)) then
                    do is = 1, site%num_shrubs
                        probs(is) = probs(is)/probssum
                    end do
                    do is = 2, site%num_shrubs
                        probs(is) = probs(is - 1) + probs(is)
                    end do
                else
                    nsrenew = 0
                end if


                ! Get current number of trees on plot
                it = site%plots(ip)%numtrees

                planting: if(site%mgmtflag == 1) then

                  ! We initialize planting routine for a species
                  is = 1
                  do while(is <= num_species)
                    if(site%species(is)%unique_id == site%planting) exit
                    is = is + 1
                  end do
                  if(is > num_species) then 
                    site%mgmtflag = 0
                    exit planting 
                  end if 
                  if(ip == 1) write(sitef, *) "Initializing planting. Calendar year", year, site%species(is)%unique_id
                  if (site%species(is)%plant_dens > 0) then
                    print *, is, site%species(is)%unique_id

                    !
                !   !  Add seedlings from planting
                !     if(ip == 1) print *, "before", site%plots(ip)%seedling(is)
                !     site%plots(ip)%seedling(is) = site%plots(ip)%seedling(is) + &
                !     site%species(is)%plant_dens/HEC_TO_M2*plotsize ! seedlings/plot
                !     if(ip == 1) print *, "after", site%plots(ip)%seedling(is)
                !     do is = 1, num_species
                !       if(ip == 1) print *, site%species(is)%unique_id, site%plots(ip)%seedling(is)
                !     end do
                !
                !   end if
                ! end if planting
                !  ! jk, adding seedlings doesn't have enough impact

                    n_empty = count(site%plots(ip)%cells(:,:) == 0)
                    plantempty: if (n_empty > 0) then

                        ! Some cells are unfilled - can place trees
                        allocate(locs(n_empty))

                        ! Loop through whole rows and columns and fill
                        ! r_empty and c_empty with empty cell indices -
                        ! keeping them together with the same 't' index
                        t = 1
                        do while (t <= n_empty)
                            do r = 1, maxcells
                                do c = 1, maxcells
                                    if (site%plots(ip)%cells(r, c) == 0) then
                                        r_empty(t) = r
                                        c_empty(t) = c
                                        locs(t) = t
                                        t = t + 1
                                    end if
                                end do
                            end do
                          end do

                          ! Shuffle locations array so we can randomly place
                          ! new trees
                          call shuffle(locs)

                        it = 0
                        plant_renew: do i = 1, int(site%species(is)%plant_dens/HEC_TO_M2*plotsize)
                              ! Increment number of plants and
                              ! initialize new tree
                              new_trees = new_trees + 1
                              it = it + 1
                              call initialize_tree(site%plots(ip)%trees(it), &
                                  site%species(is), is, tree_id, 1)
                                  tree_id = tree_id + 1
                              ! Grab the r and c values of that index
                              r = r_empty(locs(i))
                              c = c_empty(locs(i))

                              ! Set tree location to that value
                              site%plots(ip)%trees(it)%row = r
                              site%plots(ip)%trees(it)%col = c

                              ! Set this to filled in cells array
                              site%plots(ip)%cells(r,c) = 1
                              !print *, "placed"
                              ! Get dbh value of new tree
                              ! (must be between 0.5 and 2.5 cm)
                              dbh = nrand(0.2, 0.1)
                              if (dbh >= 0.5) dbh = 0.5
                              if (dbh <= 0.01) dbh = 0.01
                              site%plots(ip)%trees(it)%diam_bht = dbh

                              ! Set clear branch bole height
                              site%plots(ip)%trees(it)%canopy_ht = 1

                              ! Get other characteristics
                              call forska_height(site%plots(ip)%trees(it))
                              call stem_shape(site%plots(ip)%trees(it))
                              call biomass_c(site%plots(ip)%trees(it))
                              call biomass_n(site%plots(ip)%trees(it))
                              call leaf_biomass_c(site%plots(ip)%trees(it))

                              ! Leaf biomass
                              leafbm = site%plots(ip)%trees(it)%leaf_bm

                              ! Tree height
                              ht = max(int(site%plots(ip)%trees(it)%forska_ht), 1)

                              !calculate shading effect on tree
                              if (site%plots(ip)%trees(it)%conifer) then
                                  shade = light_rsp(site%species(is),        &
                                      site%plots(ip)%con_light(ht))
                              else
                                  shade = light_rsp(site%species(is),        &
                                      site%plots(ip)%dec_light(ht))
                              end if
                              !calculate environmental stressors
                              call env_stress(site%plots(ip)%trees(it),      &
                                  shade, site%plots(ip)%fc_gdd(is),          &
                                  site%plots(ip)%fc_drought(is),             &
                                  site%plots(ip)%fc_perm(is),                &
                                  site%plots(ip)%fc_flood(is), envstress,    &
                                  site%plots(ip)%fc_nutr(is))
                              ! Add leaf litter to soil and update NPP and
                              ! N used
                              if (site%species(is)%conifer) then

                                  NPP = NPP + leafbm*                        &
                                  (1.0 - CON_LEAF_RATIO) +                   &
                                  site%plots(ip)%trees(it)%biomC +           &
                                  site%plots(ip)%trees(it)%rootC

                                  N_used = N_used + leafbm/CON_LEAF_C_N +    &
                                      site%plots(ip)%trees(it)%biomN

                                  lc = site%species(is)%litter_class
                                  site%plots(ip)%soil%litter(lc) =           &
                                      site%plots(ip)%soil%litter(lc) +       &
                                      leafbm*CON_LEAF_RATIO/B_TO_C
                              else
                                  NPP = NPP + leafbm +                       &
                                  site%plots(ip)%trees(it)%biomC +           &
                                  site%plots(ip)%trees(it)%rootC

                                  N_used = N_used +                          &
                                      site%plots(ip)%trees(it)%biomN +       &
                                      leafbm/DEC_LEAF_C_N

                                  lc = site%species(is)%litter_class
                                  site%plots(ip)%soil%litter(lc) =           &
                                      site%plots(ip)%soil%litter(lc) +       &
                                      leafbm/B_TO_C
                              end if
                          end do plant_renew

                      end if plantempty

                  else
                    if(ip == 1) call warning("Skipping planting: planting density is 0 or species info not provided")
                  end if
                if (allocated(locs)) deallocate(locs)
                if(ip == numplots) site%mgmtflag = 0
               !  end if planting
                else 

                navail: if (site%plots(ip)%soil%avail_N > epsilon(1.0)) then
                renew: if (ntrenew >= 1 .or. nsrenew >= 1) then
                    ! Count number of unfilled cells
                    n_empty = count(site%plots(ip)%cells(:,:) == 0)

                    nempty: if (n_empty > 0) then

                        ! Some cells are unfilled - can place trees
                        allocate(locs(n_empty))


                        ! Loop through whole rows and columns and fill
                        ! r_empty and c_empty with empty cell indices -
                        ! keeping them together with the same 't' index
                        t = 1
                        do while (t <= n_empty)
                            do r = 1, maxcells
                                do c = 1, maxcells
                                    if (site%plots(ip)%cells(r, c) == 0) then
                                        r_empty(t) = r
                                        c_empty(t) = c
                                        locs(t) = t
                                        t = t + 1
                                    end if
                                end do
                            end do
                        end do

                        ! Shuffle locations array so we can randomly place
                        ! new trees
                        call shuffle(locs)

                        ntrees: if (ntrenew >= 1) then
                            ! Renew trees
                            tree_renew: do irenew = 1, ntrenew
                            if(irenew > n_empty) exit

                              if(N_used*HEC_TO_M2/plotsize >= site%plots(ip)%soil%avail_N) then
                                exit ! Stop adding trees if N is gone
                              end if

                              if(irenew > n_empty) exit

                                ! Determine species of new tree
                                rand = urand()
                                is = 1
                                do while (rand > probt(is))
                                    is = is + 1
                                    if (is > site%num_trees) then
                                        is = 1 + int(urand(0.0,                &
                                            real(site%num_trees)))
                                        rand = urand()
                                    endif
                                end do

                                ! Increment new tree counter
                                new_trees = new_trees + 1
                                is = tree_sp(is)

                                ! Can't add tree if there are no seedlings left
                                grow_seedling: if (site%plots(ip)%seedling(is) > epsilon(1.0)) then
                                ! Decrement seedling bank of species
                                site%plots(ip)%seedling(is) =                  &
                                    max(0.0, site%plots(ip)%seedling(is) - 1.0)

                                ! Increment number of plants and
                                ! initialize new tree
                                it = it + 1
                                call initialize_tree(site%plots(ip)%trees(it), &
                                    site%species(is), is, tree_id, 0)
                                    tree_id = tree_id + 1

                                ! Grab the r and c values of that index
                                r = r_empty(locs(irenew))
                                c = c_empty(locs(irenew))

                                ! Set tree location to that value
                                site%plots(ip)%trees(it)%row = r
                                site%plots(ip)%trees(it)%col = c

                                ! Set this to filled in cells array
                                site%plots(ip)%cells(r,c) = 1

                                ! Get dbh value of new tree
                                ! (must be between 0.5 and 2.5 cm)
                                dbh = 1.5 + nrand(0.0, 1.0)
                                if (dbh >= 2.5) dbh = 2.5
                                if (dbh <= 0.5) dbh = 0.5
                                site%plots(ip)%trees(it)%diam_bht = dbh

                                ! Set clear branch bole height
                                site%plots(ip)%trees(it)%canopy_ht = 1.0


                                ! Get other characteristics
                                call forska_height(site%plots(ip)%trees(it))
                                call stem_shape(site%plots(ip)%trees(it))
                                call biomass_c(site%plots(ip)%trees(it))
                                call biomass_n(site%plots(ip)%trees(it))
                                call leaf_biomass_c(site%plots(ip)%trees(it))


                                ! Leaf biomass
                                leafbm = site%plots(ip)%trees(it)%leaf_bm

                                ! Tree height
                                ht = max(int(site%plots(ip)%trees(it)%forska_ht), 1)

                                !calculate shading effect on tree
                                if (site%plots(ip)%trees(it)%conifer) then
                                    shade = light_rsp(site%species(is),        &
                                        site%plots(ip)%con_light(ht))
                                else
                                    shade = light_rsp(site%species(is),        &
                                        site%plots(ip)%dec_light(ht))
                                end if

                                !calculate environmental stressors
                                call env_stress(site%plots(ip)%trees(it),      &
                                    shade, site%plots(ip)%fc_gdd(is),          &
                                    site%plots(ip)%fc_drought(is),             &
                                    site%plots(ip)%fc_perm(is),                &
                                    site%plots(ip)%fc_flood(is), envstress,    &
                                    site%plots(ip)%fc_nutr(is))

                                ! Add leaf litter to soil and update NPP and
                                ! N used
                                if (site%species(is)%conifer) then

                                    NPP = NPP + leafbm*                        &
                                    (1.0 - CON_LEAF_RATIO) +                   &
                                    site%plots(ip)%trees(it)%biomC +           &
                                    site%plots(ip)%trees(it)%rootC

                                    N_used = N_used + leafbm/CON_LEAF_C_N +    &
                                        site%plots(ip)%trees(it)%biomN

                                    lc = site%species(is)%litter_class
                                    site%plots(ip)%soil%litter(lc) =           &
                                        site%plots(ip)%soil%litter(lc) +       &
                                        leafbm*CON_LEAF_RATIO/B_TO_C
                                else
                                    NPP = NPP + leafbm +                       &
                                    site%plots(ip)%trees(it)%biomC +           &
                                    site%plots(ip)%trees(it)%rootC

                                    N_used = N_used +                          &
                                        site%plots(ip)%trees(it)%biomN +       &
                                        leafbm/DEC_LEAF_C_N

                                    lc = site%species(is)%litter_class
                                    site%plots(ip)%soil%litter(lc) =           &
                                        site%plots(ip)%soil%litter(lc) +       &
                                        leafbm/B_TO_C
                                end if

                              end if grow_seedling

                            end do tree_renew

                            if (site%plots(ip)%trees(it)%form > 0 .and.        &
                                site%plots(ip)%trees(it)%forska_ht < MIN_HT) then
                                ! Add up shrub live foliage and fine twigs
                                site%plots(ip)%soil%shrubLitter =              &
                                    site%plots(ip)%soil%shrubLitter + leafbm + &
                                    site%plots(ip)%trees(it)%branchC*         &
                                    PERC_BRANCHES(1)

                            end if

                        end if ntrees

                        nshrubs: if (nsrenew >= 1) then

                            
                            shrub_renew: do irenew = 1, nsrenew
                              if(N_used*HEC_TO_M2/plotsize >= site%plots(ip)%soil%avail_N) then
                                exit
                              end if

                              if(irenew + new_trees > n_empty) exit

                                ! Determine species of new shrub
                                rand = urand()
                                is = 1
                                do while (rand > probs(is))
                                    is = is + 1
                                    if (is > site%num_shrubs) then
                                        is = 1 + int(urand(0.0,                &
                                            real(site%num_shrubs)))
                                        rand = urand()
                                    endif
                                end do

                                ! Increment new individual counter
                                new_shrubs = new_shrubs + 1
                                is = shrub_sp(is)

                                grow_seedling2: if (site%plots(ip)%seedling(is) > epsilon(1.0)) then

                                ! Decrement seedling bank of species
                                site%plots(ip)%seedling(is) =                  &
                                    max(0.0, site%plots(ip)%seedling(is) - 1.0)

                                ! Increment number of plants and
                                ! initialize new tree
                                it = it + 1
                                call initialize_tree(site%plots(ip)%trees(it), &
                                    site%species(is), is, tree_id, 0)
                                    tree_id = tree_id + 1

                                ! Grab the r and c values of that index
                                r = r_empty(locs(irenew + new_trees))
                                c = c_empty(locs(irenew + new_trees))

                                ! Set tree location to that value
                                site%plots(ip)%trees(it)%row = r
                                site%plots(ip)%trees(it)%col = c

                                ! Set this to filled in cells array
                                site%plots(ip)%cells(r,c) = 1

                                ! Get dbh value of new shrub
                                ! Must be between 0.1 and 0.5 cm
                                dbh = 0.3 + nrand(0.0, 1.0)
                                if (dbh >= 0.5) dbh = 0.5
                                if (dbh <= 0.1) dbh = 0.1
                                site%plots(ip)%trees(it)%diam_bht = dbh

                                ! Set clear branch bole height
                                site%plots(ip)%trees(it)%canopy_ht = 0.0

                                ! Get other characteristics
                                call forska_height(site%plots(ip)%trees(it))
                                call stem_shape(site%plots(ip)%trees(it))
                                call biomass_c(site%plots(ip)%trees(it))
                                call biomass_n(site%plots(ip)%trees(it))
                                call leaf_biomass_c(site%plots(ip)%trees(it))
                                ! Leaf biomass
                                leafbm = site%plots(ip)%trees(it)%leaf_bm

                                ! Tree height
                                ht = max(int(site%plots(ip)%trees(it)%forska_ht), 1)

                                !calculate shading effect on tree
                                if (site%plots(ip)%trees(it)%conifer) then
                                    shade = light_rsp(site%species(is),        &
                                        site%plots(ip)%con_light(ht))
                                else
                                    shade = light_rsp(site%species(is),        &
                                        site%plots(ip)%dec_light(ht))
                                end if
                                !calculate environmental stressors
                                call env_stress(site%plots(ip)%trees(it),      &
                                    shade, site%plots(ip)%fc_gdd(is),          &
                                    site%plots(ip)%fc_drought(is),             &
                                    site%plots(ip)%fc_perm(is),                &
                                    site%plots(ip)%fc_flood(is), envstress,    &
                                    site%plots(ip)%fc_nutr(is))
                                ! Add leaf litter to soil and update NPP and
                                ! N used
                                if (site%species(is)%conifer) then

                                    NPP = NPP + leafbm*                        &
                                    (1.0 - CON_LEAF_RATIO) +                   &
                                    site%plots(ip)%trees(it)%biomC +           &
                                    site%plots(ip)%trees(it)%rootC

                                    N_used = N_used + leafbm/CON_LEAF_C_N +    &
                                        site%plots(ip)%trees(it)%biomN

                                    lc = site%species(is)%litter_class
                                    site%plots(ip)%soil%litter(lc) =           &
                                        site%plots(ip)%soil%litter(lc) +       &
                                        leafbm*CON_LEAF_RATIO/B_TO_C
                                else
                                    NPP = NPP + leafbm +                       &
                                    site%plots(ip)%trees(it)%biomC +           &
                                    site%plots(ip)%trees(it)%rootC

                                    N_used = N_used +                          &
                                        site%plots(ip)%trees(it)%biomN +       &
                                        leafbm/DEC_LEAF_C_N

                                    lc = site%species(is)%litter_class
                                    site%plots(ip)%soil%litter(lc) =           &
                                        site%plots(ip)%soil%litter(lc) +       &
                                        leafbm/B_TO_C
                                end if
                            if (site%plots(ip)%trees(it)%form > 0 .and.        &
                                site%plots(ip)%trees(it)%forska_ht < MIN_HT) then
                                ! Add up shrub live foliage and fine twigs
                                site%plots(ip)%soil%shrubLitter =              &
                                    site%plots(ip)%soil%shrubLitter + leafbm + &
                                    site%plots(ip)%trees(it)%branchC*         &
                                    PERC_BRANCHES(1)

                            end if
                          end if grow_seedling2

                            end do shrub_renew

                        end if nshrubs

                    end if nempty

                end if renew

                ! Set new number of trees
                site%plots(ip)%numtrees = it

            end if navail
        end if planting

            
            site%plots(ip)%seedling = site%plots(ip)%seedling + 2.0*plotsize ! adding random germinants without the massive scaler
            ! if(ip==1) print *, site%plots(ip)%seedling, "seedlings"

            ! Decrease seedling bank for survivability
            ! Also convert to #/m2
            do is = 1, num_species
                site%plots(ip)%seedling(is) = (site%plots(ip)%seedling(is)* &
                    ! regstress(is)* &
                    site%species(is)%seedling_surv)/plotsize
            end do

            ! if(init_on == .false. .and. year - start_year < 5) then !boosting invader establishment early on
            !   site%plots(ip)%seedling = site%plots(ip)%seedling*10
            !   site%plots(ip)%seedbank = site%plots(ip)%seedbank*10
            ! end if


            ! Update site and soil variables
            N_used = N_used*HEC_TO_M2/plotsize
            site%plots(ip)%NPP = site%plots(ip)%NPP + NPP*HEC_TO_M2/plotsize
            site%plots(ip)%soil%shrubLitter =                                  &
                site%plots(ip)%soil%shrubLitter/B_TO_C/plotsize
            ! if(ip==1) print*, "N used ", N_used, " avail ", site%plots(ip)%soil%avail_N
            ! Convert to tonnes/ha
            do i = 1, 18
                site%plots(ip)%soil%litter(i) =                                &
                    site%plots(ip)%soil%litter(i)*HEC_TO_M2/plotsize
            end do

            ! Deallocate locs array
            if (allocated(locs)) deallocate(locs)

            ! Write out regeneration data if we are testing
            if (reg_testing) then
                do is = 1, num_species
                    call csv_write(spec_regen, site%site_id, .false.)
                    call csv_write(spec_regen, ip, .false.)
                    call csv_write(spec_regen, year, .false.)
                    call csv_write(spec_regen, probtsum, .false.)
                    call csv_write(spec_regen, probssum, .false.)
                    call csv_write(spec_regen, site%species(is)%unique_id,     &
                        .false.)
                    call csv_write(spec_regen, site%plots(ip)%seedling(is),    &
                        .false.)
                    call csv_write(spec_regen, site%plots(ip)%seedbank(is),    &
                        .false.)
                    call csv_write(spec_regen, site%plots(ip)%fc_gdd(is), .false.)
                    call csv_write(spec_regen, site%plots(ip)%fc_nutr(is), .false.)
                    call csv_write(spec_regen, site%plots(ip)%fc_drought(is), .false.)
                    call csv_write(spec_regen, site%plots(ip)%fc_flood(is), .false.)
                    call csv_write(spec_regen, site%plots(ip)%fc_perm(is), .false.)
                    call csv_write(spec_regen, site%plots(ip)%con_light(1), .false.)
                    call csv_write(spec_regen, site%plots(ip)%dec_light(1), .false.)
                    call csv_write(spec_regen, regstress(is), .true.)
                end do
            end if

        end do plot
        
        ! Summarize regeneration data if testing
        if (mgmt_testing) then 
            do is = 1, num_species
               call stddev(regstress_all(:, is), regstress_mn(is), regstress_sd(is),  &
               RNVALID)   
                   call csv_write(siteregen, site%site_id, .false.)
                   call csv_write(siteregen, year, .false.)
                   call csv_write(siteregen, site%species(is)%unique_id,     &
                       .false.)
                   call csv_write(siteregen, regstress_mn(is), .false.)
                   call csv_write(siteregen, regstress_sd(is), .true.)
               end do 
           end if 
    end subroutine Renewal

    !:.........................................................................:

  	subroutine Init_Stand(site, year, tree_id)

  		type(SiteData),                      intent(inout) :: site
  		integer,                             intent(in)    :: year
      integer,                             intent(inout) :: tree_id !Unique tree counter

  		real,    dimension(size(site%species_ids))             :: perc_spec ! Convert species counts to percentages
  		real,    dimension(size(site%species_ids))             :: cumm_perc ! Cumulative percentage vector
  		real,    dimension(size(site%species_ids))             :: tot_spec ! Numner of species on site
  		real,    dimension(size(site%species_ids), 10)         :: dbh_bins ! Trees in each DBH bin
  		real,    dimension(size(site%species_ids), 10)         :: perc_dbh ! Convert DBH counts to percentages
  		real,    dimension(size(site%species_ids), 10)         :: cumm_dbh ! Cumulative dbh vector

  		integer, dimension(:),    allocatable              :: r_empty !empty rows
  		integer, dimension(:),    allocatable              :: c_empty ! empty columns

  		real,    dimension(11)                             :: dbh_lims ! cutpoints for DBH bins
  		real                                               :: leaf_b ! leaf biomass
  		real                                               :: N_used ! Nitrogen used
  		real                                               :: net_prim_prod !Net primary productivity   
  		real                                               :: q0, qdbh ! To get dbh value for new trees
  		real                                               :: rand_t ! Holder to choose location of new trees
  		real                                               :: z, zz !holders
  		real                                               :: uconvert ! Conversion value
  		real                                               :: canopy_shade ! To calculate light environment after initializing trees
  		real                                               :: ilitter ! initial litter
  		real                                               :: avtt !Nitrogen holder
      real                                               :: envstress !holder
  		integer                                            :: num_species ! Number of species on site
  		integer                                            :: new ! Counter for trees
  		integer                                            :: ip ! Looping (plots)
  		integer                                            :: n_init ! Number of trees to initialize
  		integer                                            :: i_init ! Looping index
  		integer                                            :: n_empty ! Number of empty "cells" for trees
  		integer                                            :: s, is, it ! Looping indices
  		integer                                            :: k, t, r, c ! Looping indices
  		integer                                            :: kh, d, lc ! Temporary values 
  		integer                                            :: lit, hc ! Temporary values 
      real                                               :: bct, bcr, bcs, bcbr, d_bc, d_bcs, d_bcr, d_bcbr, d_bctw, d_bcsb, d_bclb, leafbm, d_leafb ! Branch stuff (tC)
      real                                               :: con_la, dec_la, tot_la, m, hum_input, hum_mod ! for calculating species dominance, humus depth
      real                                               :: net_N_into_A0 !A0 is new N in humus
      real                                               :: temp_BA  !counter; m2/ha

      site%species%native = .false.

      con_la = 0.0
      dec_la = 0.0
  		!dbh bins
  		data dbh_lims /2.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0,   &
  						70.0, 80.0, 90.0/

  		!calculate conifer leaf ratio
  		leaf_b = 1.0 + con_leaf_ratio

  		!get number of species in initialization files
  		num_species = size(site%species_ids)

  		perc_spec = site%percent_sp
  		dbh_bins = site%dbh_means

  		!calculate cumulative percent species
  		cumm_perc(1) = perc_spec(1)
  		do s = 2, num_species
  			cumm_perc(s) = perc_spec(s) + cumm_perc(s - 1)
  		end do

  		!calculate total, percent, and cumulative percent stems
  		!for each species
      cumm_dbh = 0.0
  		do s = 1, num_species
  			tot_spec(s) = sum(dbh_bins(s, :))
  			if (tot_spec(s) .le. 0.0) then
  				perc_dbh(s, :) = 0.0
  			else
  				perc_dbh(s,:) = dbh_bins(s, :)/tot_spec(s)
  			end if
  			cumm_dbh(s, 1) = perc_dbh(s, 1)
  			do d = 2, 10
  				cumm_dbh(s, d) = perc_dbh(s, d) + cumm_dbh(s, d - 1)
  			end do

  		end do

  		do ip = 1, site%numplots

  			!initialize new trees, NPP, N_used
  			new = 0
  			net_prim_prod = 0.0
  			N_used = 0.0

  			if (site%plots(ip)%soil%avail_N .gt. 0.0) then
  				!recalculate species-level response to soil moisture
  				!and temperature
  				do is = 1, size(site%species)
  					call temp_rsp(site%species(is),                     &
  						site%deg_days, site%plots(ip)%fc_gdd(is))
  					call drought_rsp(site%species(is),                  &
  						site%plots(ip)%dry_days,           &
  						site%plots(ip)%fc_drought(is))
  					call flood_rsp(site%species(is),                    &
  						site%plots(ip)%flood_days, site%plots(ip)%fc_flood(is))
  					call perm_rsp(site%species(is)%perm_tol,            &
  						site%plots(ip)%soil%active, site%plots(ip)%fc_perm(is))
  				end do

    				it = site%plots(ip)%numtrees

  				!find out how many total trees to initialize
  				n_init = int((site%total_stems/HEC_TO_M2)*plotsize)

  				!have to initialize at least 1 tree
  				if (n_init .ge. 1) then

  					!initialize trees
                    do i_init = 1, n_init

                        !determine species of tree
                        q0 = urand()
                        do k = 1, num_species
                            if (q0 .lt. cumm_perc(k)) then
                                exit
                            end if
                        end do

                        do is = 1, size(site%species)
                            if(site%species_ids(k) .eq. site%species(is)%unique_id) then
                            if(site%species(is)%native .eq. .false.) then
                                site%species(is)%native = .true.
                            end if
                            exit
                            end if
                        end do
                         ! k is index for init file inputs, "is" is index for species from rangelist

  						!increment new tree counter
  						new = new + 1

  						!increment number of trees and initialize new tree
  						it = it + 1
  						call initialize_tree(site%plots(ip)%trees(it), &
  								site%species(is), is, tree_id, 2)
                         tree_id = tree_id + 1
  						!get tree location ----

  						!count number of unfilled cells
  						n_empty = count(site%plots(ip)%cells(:,:) == 0)

  						if (n_empty .gt. 0) then

  							!some cells are unfilled - can place new trees

  							!allocate lists of empty rows and cols to n_empty
  							allocate(r_empty(n_empty))
  							allocate(c_empty(n_empty))

  							!loop through whole rows and cols and fill
  							! r_emtpy and c_empty with empty cell
  							!indices - keeping them together with the
  							!same 't' index
  							t = 1
  							do while (t .le. n_empty)
  								do r = 1, maxcells
  									do c = 1, maxcells
  										if (site%plots(ip)%cells(r,c)  &
  												.eq. .false.) then
  											r_empty(t) = r
  											c_empty(t) = c
  											t = t + 1
  										end if
  									end do
  								end do
  							end do

  							!pull random number between 1 and n_empty
  							!convert to integer
  							rand_t = int(urand(1.0, float(n_empty)))

  							!grab the r and c values of that index
  							r = r_empty(rand_t)
  							c = c_empty(rand_t)

  							!set tree location to that value
  							site%plots(ip)%trees(it)%row = r
  							site%plots(ip)%trees(it)%col = c

  							!set this to filled in cells array
  							site%plots(ip)%cells(r,c) = .true.

  							!deallocate r and c empty lists
  							if (allocated(r_empty)) deallocate(r_empty)
  							if (allocated(c_empty)) deallocate(c_empty)

  							!get dbh value of new tree
  							qdbh = urand()
  							do d = 1, 10
  								if (qdbh .lt. cumm_dbh(k, d)) then
  									exit
  								end if
  							end do

  							z = urand(dbh_lims(d), dbh_lims(d + 1))
  							if (z .le. dbh_lims(d)) z = dbh_lims(d)
  							if (z .ge. dbh_lims(d + 1)) z = dbh_lims(d + 1)
  							site%plots(ip)%trees(it)%diam_bht = z

  							!get height, canopy diameter, and biomass
  							call forska_height(site%plots(ip)%trees(it))
                            !set canopy height
                            site%plots(ip)%trees(it)%canopy_ht = max(1.0, &
                            site%plots(ip)%trees(it)%forska_ht*0.4)

  							call stem_shape(site%plots(ip)%trees(it))
  							call biomass_c(site%plots(ip)%trees(it))
  							call biomass_n(site%plots(ip)%trees(it))
  							call leaf_biomass_c(site%plots(ip)%trees(it))

                            ! set tree age (assuming optimal growth)
                            ! call init_age(site%plots(ip)%trees(it))

                            !working on soil spinup
                            ! figure out what to keep wrt litter. And read soil stuff closely

  							zz = leaf_area(site%plots(ip)%trees(it))*  &
  								site%species(is)%leafarea_c*2.0


  							!add leaf litter to soil and update N_used
  							!and NPP
  							if (site%species(is)%conifer) then
  								net_prim_prod = net_prim_prod +        &
  									zz*leaf_b +                        &
  									site%plots(ip)%trees(it)%biomC

  								N_used = N_used + zz/CON_LEAF_C_N +    &
  									site%plots(ip)%trees(it)%biomN

  								! site%plots(ip)%soil%C_into_A0 =        &
  								! 	site%plots(ip)%soil%C_into_A0 +    &
  								! 	zz*(leaf_b - 1.0)

  								lc = site%species(is)%litter_class
  								site%plots(ip)%soil%litter(lc) =       &
  									site%plots(ip)%soil%litter(lc) +   &
  									zz*(leaf_b - 1.0)/B_TO_C

                                 ! Don't know if this is the right place for leaf litter N
  								! site%plots(ip)%soil%cohorts(1,2) =        &
  								! 	site%plots(ip)%soil%cohorts(1,2) +    &
  								! 	zz*(leaf_b - 1.0)/con_leaf_c_n

                                 con_la = con_la + leaf_area(site%plots(ip)%trees(it))

  							else
  								net_prim_prod = net_prim_prod +        &
  									site%plots(ip)%trees(it)%biomC + zz

  								N_used = N_used +                      &
  									site%plots(ip)%trees(it)%biomN +   &
  									zz/DEC_LEAF_C_N

  							! 	site%plots(ip)%soil%C_into_A0 =        &
  							! 		site%plots(ip)%soil%C_into_A0 + zz
                !
  								lc = site%species(is)%litter_class
  								site%plots(ip)%soil%litter(lc) =       &
  									site%plots(ip)%soil%litter(lc) +   &
  									zz/B_TO_C

  								! site%plots(ip)%soil%cohorts(1,2) =        &
  								! 	site%plots(ip)%soil%cohorts(1,2) +    &
  								! 	zz/dec_leaf_c_n

                dec_la = dec_la + leaf_area(site%plots(ip)%trees(it))

                             end if

  						end if !end if empty cells

  					end do !end i_init loop

                else if(n_init .lt. 1 .and. site%init_BA .gt. 0) then ! initializing to basal area 

                    !initialize trees
                    temp_BA = 0.0
                    do while(temp_BA < (site%init_BA)/HEC_TO_M2*plotsize)

                        !determine species of tree
                        q0 = urand()
                        do k = 1, num_species
                            if (q0 .lt. cumm_perc(k)) then
                                exit
                            end if
                        end do

                        do is = 1, size(site%species)
                            if(site%species_ids(k) .eq. site%species(is)%unique_id) then
                            if(site%species(is)%native .eq. .false.) then
                                site%species(is)%native = .true.
                            end if
                            exit
                            end if
                        end do
                         ! k is index for init file inputs, "is" is index for species from rangelist

  						!increment new tree counter
  						new = new + 1

  						!increment number of trees and initialize new tree
  						it = it + 1
  						call initialize_tree(site%plots(ip)%trees(it), &
  								site%species(is), is, tree_id, 2)
                         tree_id = tree_id + 1
  						!get tree location ----

  						!count number of unfilled cells
  						n_empty = count(site%plots(ip)%cells(:,:) == 0)

  						if (n_empty .gt. 0) then

  							!some cells are unfilled - can place new trees

  							!allocate lists of empty rows and cols to n_empty
  							allocate(r_empty(n_empty))
  							allocate(c_empty(n_empty))

  							!loop through whole rows and cols and fill
  							! r_emtpy and c_empty with empty cell
  							!indices - keeping them together with the
  							!same 't' index
  							t = 1
  							do while (t .le. n_empty)
  								do r = 1, maxcells
  									do c = 1, maxcells
  										if (site%plots(ip)%cells(r,c)  &
  												.eq. .false.) then
  											r_empty(t) = r
  											c_empty(t) = c
  											t = t + 1
  										end if
  									end do
  								end do
  							end do

  							!pull random number between 1 and n_empty
  							!convert to integer
  							rand_t = int(urand(1.0, float(n_empty)))

  							!grab the r and c values of that index
  							r = r_empty(rand_t)
  							c = c_empty(rand_t)

  							!set tree location to that value
  							site%plots(ip)%trees(it)%row = r
  							site%plots(ip)%trees(it)%col = c

  							!set this to filled in cells array
  							site%plots(ip)%cells(r,c) = .true.

  							!deallocate r and c empty lists
  							if (allocated(r_empty)) deallocate(r_empty)
  							if (allocated(c_empty)) deallocate(c_empty)

  							!get dbh value of new tree
  							qdbh = urand()
  							do d = 1, 10
  								if (qdbh .lt. cumm_dbh(k, d)) then
  									exit
  								end if
  							end do

  							z = urand(dbh_lims(d), dbh_lims(d + 1))
  							if (z .le. dbh_lims(d)) z = dbh_lims(d)
  							if (z .ge. dbh_lims(d + 1)) z = dbh_lims(d + 1)
  							site%plots(ip)%trees(it)%diam_bht = z

                            temp_BA = temp_BA + pi*(z*CM_TO_M/2)**2 !convert z to m and radius

  							!get height, canopy diameter, and biomass
  							call forska_height(site%plots(ip)%trees(it))
                            !set canopy height
                            site%plots(ip)%trees(it)%canopy_ht = max(1.0, &
                            site%plots(ip)%trees(it)%forska_ht*0.4)

  							call stem_shape(site%plots(ip)%trees(it))
  							call biomass_c(site%plots(ip)%trees(it))
  							call biomass_n(site%plots(ip)%trees(it))
  							call leaf_biomass_c(site%plots(ip)%trees(it))

                            ! set tree age (assuming optimal growth)
                            ! call init_age(site%plots(ip)%trees(it))

                            !working on soil spinup
                            ! figure out what to keep wrt litter. And read soil stuff closely

  							zz = leaf_area(site%plots(ip)%trees(it))*  &
  								site%species(is)%leafarea_c*2.0


  							!add leaf litter to soil and update N_used
  							!and NPP
  							if (site%species(is)%conifer) then
  								net_prim_prod = net_prim_prod +        &
  									zz*leaf_b +                        &
  									site%plots(ip)%trees(it)%biomC

  								N_used = N_used + zz/CON_LEAF_C_N +    &
  									site%plots(ip)%trees(it)%biomN

  								! site%plots(ip)%soil%C_into_A0 =        &
  								! 	site%plots(ip)%soil%C_into_A0 +    &
  								! 	zz*(leaf_b - 1.0)

  								lc = site%species(is)%litter_class
  								site%plots(ip)%soil%litter(lc) =       &
  									site%plots(ip)%soil%litter(lc) +   &
  									zz*(leaf_b - 1.0)/B_TO_C

                                 ! Don't know if this is the right place for leaf litter N
  								! site%plots(ip)%soil%cohorts(1,2) =        &
  								! 	site%plots(ip)%soil%cohorts(1,2) +    &
  								! 	zz*(leaf_b - 1.0)/con_leaf_c_n

                                 con_la = con_la + leaf_area(site%plots(ip)%trees(it))

  							else
  								net_prim_prod = net_prim_prod +        &
  									site%plots(ip)%trees(it)%biomC + zz

  								N_used = N_used +                      &
  									site%plots(ip)%trees(it)%biomN +   &
  									zz/DEC_LEAF_C_N

  							! 	site%plots(ip)%soil%C_into_A0 =        &
  							! 		site%plots(ip)%soil%C_into_A0 + zz
                !
  								lc = site%species(is)%litter_class
  								site%plots(ip)%soil%litter(lc) =       &
  									site%plots(ip)%soil%litter(lc) +   &
  									zz/B_TO_C

  								! site%plots(ip)%soil%cohorts(1,2) =        &
  								! 	site%plots(ip)%soil%cohorts(1,2) +    &
  								! 	zz/dec_leaf_c_n

                dec_la = dec_la + leaf_area(site%plots(ip)%trees(it))

                             end if

  						end if !end if empty cells

  					end do !end i_init loop

  				end if !end if trees to initialize

  				site%plots(ip)%numtrees = it

                
  			end if !end if enough available N

  			!update site and soil variables
  			uconvert = HEC_TO_M2/plotsize
  			N_used = N_used*uconvert
  			net_prim_prod = net_prim_prod*uconvert
  			avtt = site%plots(ip)%soil%avail_N - N_used

        ! Nitrogen in soil
        ! print *, "Nused", N_used, "avtt", avtt
        ! site%plots(ip)%soil%avail_N = max(0.0, avtt)

  			! if (avtt .gt. 0.0) then
  			! 	net_N_into_A0 = avtt*min(site%plots(ip)%soil%runoff/1000.0,0.1)
  			! 	site%plots(ip)%soil%cohorts(1,2) = site%plots(ip)%soil%cohorts(1,2) +  &
  			! 		avtt - net_N_into_A0
  			! else
  			! 	site%plots(ip)%soil%cohorts(1,2) = site%plots(ip)%soil%cohorts(1,2) +  &
  			! 		avtt
  			! 	net_N_into_A0 = 0.0
  			! end if

        if(avtt>epsilon(1.0)) then
          site%plots(ip)%soil%avail_N = avtt
        else
          site%plots(ip)%soil%avail_N = 0.0
        end if

        ! print *, site%plots(ip)%soil%avail_N

  			! site%plots(ip)%soil%cohorts(1,2) = site%plots(ip)%soil%cohorts(1,2) -      &
  			! 	0.00002*site%plots(ip)%soil%runoff
  			! ! site%plots(ip)%soil%A_c0 = site%plots(ip)%soil%A_c0 -      &
  			! ! 	net_N_into_A0*uconvert
  			! ! site%plots(ip)%soil%BL_c0 = site%plots(ip)%soil%BL_c0 +    &
  			! ! 	net_N_into_A0*uconvert
  			! site%plots(ip)%soil%min_Nmin = site%plots(ip)%soil%min_Nmin +    &
  			! 	net_N_into_A0
  			! site%plots(ip)%soil%C_into_A0 =                            &
  			! 	(site%plots(ip)%soil%C_into_A0)*uconvert
          ! site%plots(ip)%soil%cohorts =                            &
          !   (site%plots(ip)%soil%cohorts)*uconvert
  			 site%plots(ip)%soil%N_used = N_used
  			! site%plots(ip)%soil%net_prim_prodC =                       &
  			! site%plots(ip)%soil%net_prim_prodC + net_prim_prod
  			! site%plots(ip)%soil%net_prim_prodN =                       &
  			! 	site%plots(ip)%soil%net_prim_prodN + N_used
  			! site%plots(ip)%soil%new_growth = int(float(new)*uconvert)

  		end do !end plot loop
      ! do ip = 1, numplots
      !   do it = 1, site%plots(ip)%numtrees
      !       ! Adding coarse woody debris
      !       ! Idea: for each tree, figure out previous droppings and "decay rate" for how far to go back dropping branches
      !       ! Trees drop branches when stressed. Will use envstress to calc who drops ... nvm, trees aren't stressed at initialization
      !
      !       !Tree is big - will remove some of lower branches and add *part* of this to current cwd
      !       branchfall: if(site%plots(ip)%trees(it)%diam_bht >= 10.0) then
      !         hc = max(int(site%plots(ip)%trees(it)%forska_ht*0.2), 1) ! remove branches from bottom 1/5 of tree
      !         htcheck: if (hc < int(site%plots(ip)%trees(it)%forska_ht)) then
      !           site%plots(ip)%trees(it)%canopy_ht = float(hc)
      !           ! Update diameter at clear branch bole height
      !           call stem_shape(site%plots(ip)%trees(it))
      !
      !           ! Save old woody biomass
      !           bct = site%plots(ip)%trees(it)%biomC +             &
      !               site%plots(ip)%trees(it)%rootC
      !
      !           ! Save old root, twig, stem C biomass
      !           bcr = site%plots(ip)%trees(it)%rootC
      !           bcbr = site%plots(ip)%trees(it)%branchC
      !           bcs = site%plots(ip)%trees(it)%stemC
      !
      !           ! Update biomass C and N given new clear branch
      !           ! bole height
      !           call biomass_c(site%plots(ip)%trees(it))
      !           call biomass_n(site%plots(ip)%trees(it))
      !
      !           ! How much wood litter did we lose?
      !           d_bc = bct - (site%plots(ip)%trees(it)%biomC +     &
      !               site%plots(ip)%trees(it)%rootC)
      !           d_bcr = bcr - site%plots(ip)%trees(it)%rootC
      !           d_bcbr = bcbr - site%plots(ip)%trees(it)%branchC
      !           d_bcs = bcs -  site%plots(ip)%trees(it)%stemC
      !
      !           ! Reduce branch biomass to account for gradual senescence and decay
      !           d_bcbr = d_bcbr*0.2
      !           d_bcs = d_bcs*0.2
      !
      !           ! Divide branch biomass into twigs, large branches,
      !           ! and small branches (Thonicke et al. 2010)
      !           d_bctw = d_bcbr*PERC_BRANCHES(1)
      !           d_bcsb = d_bcbr*PERC_BRANCHES(2)
      !           d_bclb = d_bcbr*PERC_BRANCHES(3)
      !
      !           ! Add litter loss to litter pools
      !           ! Here we convert to dry biomass because soiln
      !           ! subroutine calculates weight loss not C loss
      !
      !           ! if(ip==1) print *, d_bcbr
      !           ! Roots
      !           site%plots(ip)%soil%litter(IROOT) =                &
      !               site%plots(ip)%soil%litter(IROOT) + d_bcr/B_TO_C
      !
      !           ! Twigs
      !           site%plots(ip)%soil%litter(ITW) =                   &
      !               site%plots(ip)%soil%litter(ITW) + d_bctw/B_TO_C
      !
      !           ! Small branches
      !           site%plots(ip)%soil%litter(ISBR) =                 &
      !               site%plots(ip)%soil%litter(ISBR) + d_bcsb/B_TO_C
      !
      !           !large branches
      !           site%plots(ip)%soil%litter(ILBR) =                 &
      !               site%plots(ip)%soil%litter(ILBR) + d_bclb/B_TO_C
      !
      !           ! large wood
      !             site%plots(ip)%soil%litter(ISBL) =             &
      !                 site%plots(ip)%soil%litter(ISBL) +         &
      !                 d_bcs/B_TO_C
      !
      !           ! Save previous value of leaf biomass
      !           leafbm = site%plots(ip)%trees(it)%leaf_bm
      !
      !           ! Update leaf bm and get difference
      !           call leaf_biomass_c(site%plots(ip)%trees(it))
      !           d_leafb = leafbm - site%plots(ip)%trees(it)%leaf_bm
      !
      !           d_leafb = d_leafb * 0.1 ! Reducing to account for gradual decay
      !
      !           ! Add that litter to correct leaf litter class
      !           lc = site%species(is)%litter_class
      !           site%plots(ip)%soil%litter(lc) =                   &
      !               site%plots(ip)%soil%litter(lc) + d_leafb/B_TO_C
      !
      !         end if htcheck
      !
      !       end if branchfall
      !
      !     end do
      !   end do

            ! Now that all trees are initialized, calculate light environment and envstress
            call Canopy(site)

            do ip = 1, numplots
              do it = 1, site%plots(ip)%numtrees

                kh = int(site%plots(ip)%trees(it)%forska_ht)

                !calculate shading effect on tree
                if (site%plots(ip)%trees(it)%conifer) then
                  canopy_shade = light_rsp(site%species(is),  &
                    site%plots(ip)%con_light(kh))
                else
                  canopy_shade = light_rsp(site%species(is),  &
                    site%plots(ip)%dec_light(kh))
                end if

              do is = 1, size(site%species)
                if(site%plots(ip)%trees(it)%spec_ptr%unique_id .eq. site%species(is)%unique_id) then
                  exit
                end if
              end do

              !calculate environmental stress
              call env_stress(site%plots(ip)%trees(it),      &
                  canopy_shade, site%plots(ip)%fc_gdd(is),   &
                  site%plots(ip)%fc_drought(is),             &
                  site%plots(ip)%fc_perm(is),                &
                  site%plots(ip)%fc_flood(is), envstress,    &
                  site%plots(ip)%fc_nutr(is))

          end do ! tree loop
          !convert to tonnes/ha
          do lit = 1, 20
            !if (lit .ne. 17) then
              ilitter = site%plots(ip)%soil%litter(lit)
              if (ilitter .gt. 0.00001) then
                site%plots(ip)%soil%litter(lit) =              &
                  site%plots(ip)%soil%litter(lit)*uconvert
              else
                site%plots(ip)%soil%litter(lit) = 0.0
              end if
            !end if
          end do

      end do !plot loop

  	end subroutine Init_Stand

  	!:.................................................................:

end module Model
