module Site

!*******************************************************************************
  !
  ! The Site module contains the definition of the SiteData type, which
  ! holds attributes (physical and geographical) of a particular site,
  ! as well as procedures directly operating on site variables.
  !
!*******************************************************************************

    use Constants
    use Input
    use Species
    use Soil
    use Plot
    use Utilities
    use lists

    implicit none

    ! Define SiteData types
    type SiteData
        type(PlotData),    dimension(:),     allocatable :: plots           ! Array of plot objects
        type(SpeciesData), dimension(:),     allocatable :: species         ! Array of species objects
        real,              dimension(:),     allocatable :: recr_trees      ! Reproductively active trees (t/ha)
        real,              dimension(:),     allocatable :: lai_array       ! LAI with height (m2/m2)
        real,              dimension(:,:),   allocatable :: la_spec        ! LAI with height by species (m2/m2)
        real,              dimension(NTEMPS)             :: temp_lapse_r    ! Temperature lapse rate (degC/km)
        real,              dimension(NTEMPS)             :: precip_lapse_r  ! Precipitation lapse rate (mm/km)
        real,              dimension(NTEMPS)             :: tmin            ! Mean monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: tmax            ! Mean monthly maximum temperature (degC)
        real,              dimension(NTEMPS)             :: precip          ! Mean monthly precipitation (cm)
        real,              dimension(NTEMPS)             :: tmin_std        ! SD of monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: tmax_std        ! SD of monthly maximum temperature (degC)
        real,              dimension(NTEMPS)             :: precip_std      ! SD of monthly precipitation (cm)
        real,              dimension(NTEMPS)             :: cld             ! Mean monthly cloudiness (tenths of sky covered)
        real,              dimension(NTEMPS)             :: cld_std         ! SD of monthly cloudiness (tenths of sky covered)
        real,              dimension(NTEMPS)             :: rh              ! Mean monthly relative humidity (%)
        real,              dimension(NTEMPS)             :: rh_std          ! SD of monthly relative humidity (%)
        real,              dimension(NTEMPS)             :: wind            ! Mean monthly wind speed (m/s)
        real,              dimension(NTEMPS)             :: strikes         ! Mean monthly lightning (strikes/km2/day)
        real,              dimension(NTEMPS)             :: strikes_std     ! SD of monthly lightning (strikes/km2/day)
        real,              dimension(NTEMPS)             :: accum_precip    ! Accumulated precipitation for linear clim change (cm)
        real,              dimension(DAYS_PER_YEAR)      :: maxdaytemp      ! Daily maximum temperature (degC)
        real,              dimension(DAYS_PER_YEAR)      :: dayprecip       ! Daily precipitation (cm)
        character(len = MAX_NLEN)                        :: region          ! Site region
        character(len = MAX_NLEN)                        :: site_name       ! Site name
        integer                                          :: site_id         ! Site ID
        integer                                          :: runID           ! Run ID
        integer                                          :: numplots        ! Number of plots
        integer                                          :: fire_year       ! Year to burn plots
        integer                                          :: fire_day        ! Day to burn plots
        integer                                          :: num_trees       ! Number of tree species
        integer                                          :: num_shrubs      ! Number of shrub species
        real                                             :: latitude        ! Site latitude (degrees)
        real                                             :: longitude       ! Site longitude (degrees)
        real                                             :: elevation       ! Site elevation (m)
        real                                             :: altitude        ! Site altitude (for altitude climate adjustments) (m)
        real                                             :: slope           ! Site slope (degrees)
        real                                             :: aspect          ! Site aspect (degrees)
        real                                             :: leaf_area_ind   ! Leaf area index (m2/m2)
        real                                             :: rain            ! Annual rainfall (cm)
        real                                             :: pot_evap_day    ! Annual potential evapotranspiration (cm)
        real                                             :: solar           ! Annual solar radiation (cal/cm2/year)
        real                                             :: grow_days       ! Growing season (>5degC) length (days)
        real                                             :: deg_days        ! Growing degree-days (>5degC)
        real                                             :: wind_prob       ! Windthrow probability (0-1)
        real                                             :: accum_tmax      ! Accumulated minimum temperature for linear climate
                                                                            ! change (degC)
        real                                             :: accum_tmin      ! Accumulated maximum temperature for linear climate
                                                                            ! change (degC)
        real                                             :: flow            ! Moisture input from overland flow (mm)
        real                                             :: alff            ! Available light on the forst floor (0-1)
        real                                             :: pc_germ         ! Impact of temperature on germination
        real                                             :: fire_wind       ! Wind speed for forced fire (m/s)
        real                                             :: fire_ffmc       ! FFMC for forced fire
        real                                             :: fire_dmc        ! DMC for forced fire
        integer                                          :: ind             ! Site index (maps to row/column)
        integer                                          :: row, col        ! Site row and column
        integer                                          :: n_adj           ! Number of adjacent sites contributing to seed rain
        logical                                          :: management      ! Are we conducting management?
        logical                                          :: thinning        ! Are we thinning?
        logical                                          :: shearblading    ! Are we shearblading?
        logical                                          :: pruning    ! Are we shearblading?
        logical                                          :: sel_cutting     ! Are we selectively cutting?
        logical                                          :: buldozed        ! Did we just buldoze soil?
        logical                                          :: prescribed_burn ! Are we burning stand?
        logical                                          :: already_cut     ! Did we already cut?
        integer                                          :: rotation_time   ! How often to cut?
        character (len = 8)                              :: planting          ! If cutting, what species will be planted?
        real                                             :: viability         ! If planting, what proportion of seedlings survive?
        integer                                          :: year_management ! Year to conduct management eg. 1950, not 100
        real                                             :: thin_perc       ! Percent to thin in management (0-1)
        real                                             :: dbh_thresh      ! DBH threshold for selective cutting (cm)
        real                                             :: psolrem         ! Proportion to remove soil in shearblading (0-1)
        real                                             :: summerVPD         !Vapor pressure deficit in Jul/Aug (millibars)
        real                                             :: total_stems     ! For site initialization (per ha)  
        real                                             :: init_BA     ! For site initialization (m2/ per ha)
        real, dimension(:), allocatable                  :: percent_sp     ! For site initialization
        character(len = 8), dimension(:), allocatable    :: species_ids    ! For site initialization
        real, dimension(:, :), allocatable               :: dbh_means     ! For site initialization
        character(len=3) :: unit ! Unit designation for active management
        integer :: init_year ! What year is the stand being initialized in? (In case inventory has different years)
        character(len = 4) :: classification ! for active management
        integer :: mgmtflag ! Flag for active management; 0 = nothing, 1 = harvest (clearcut), 2 = thin, 3 = salvage, 4 = shearblade, 5 = prune


    end type SiteData

contains

    !:...........................................................................:

    subroutine initialize_site(self, site_vals, species_data, species_ids, sndx)
        !
        !  Inititalizes a site with input site parameters and starting values.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb          Original Code
        !    01/26/21     A. C. Foster        Separated attach_species and
        !                                       initialize_plot calls and skip
        !                                       site if bad rangelist data
        !

        ! Data dictionary: calling arguments
        class(SiteData),                    intent(inout) :: self         ! Site object
        type(SpeciesData),  dimension(:),   intent(inout) :: species_data ! Array of species objects
        integer,            dimension(:,:), intent(inout) :: species_ids  ! Array of species present at each site
        real,               dimension(:),   intent(in)    :: site_vals    ! Site-level runtime parameters
        integer,                            intent(in)    :: sndx         ! Site index in site array

        ! Data dictionary: local variables
        character(len = 8), dimension(:),    allocatable :: range_species_ids ! Unique IDs of species at this site
        character(len = 8), dimension(:),    allocatable :: init_species_ids ! Unique IDs of species at this site
        real,              dimension(NTEMPS)             :: tmin              ! Mean monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: tmax              ! Mean monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: prcp              ! Mean monthly precipitation (mm)
        real,              dimension(NTEMPS)             :: tmin_std          ! SD of monthly minimum temperature (degC)
        real,              dimension(NTEMPS)             :: tmax_std          ! SD of monthly maximum temperature (degC)
        real,              dimension(NTEMPS)             :: prcp_std          ! SD of monthly precipitation (mm)
        real,              dimension(NTEMPS)             :: cld               ! Mean monthly cloud cover (tenths of sky covered)
        real,              dimension(NTEMPS)             :: cld_std           ! SD monthly cloud cover (tenths of sky covered)
        real,              dimension(NTEMPS)             :: rh                ! Mean monthly relative humidity (%)
        real,              dimension(NTEMPS)             :: rh_std            ! SD of monthly relative humidity (%)
        real,              dimension(NTEMPS)             :: wind              ! Mean monthly wind speed (m/s)
        real,              dimension(NTEMPS)             :: strikes           ! Mean monthly lightning (strikes/km2/day)
        real,              dimension(NTEMPS)             :: strikes_sd        ! SD of monthly lightning (strikes/km2/day)
        character(len = MAX_NLEN)                        :: sitename          ! Name of site
        character(len = MAX_NLEN)                        :: siteregion        ! Site region
        character(len = MAX_CHAR)                        :: message           ! Warning message
        real                                             :: fire_year         ! Year to burn sites
        real                                             :: fire_day          ! Day to burn sites
        real                                             :: lat               ! Latitude (degrees)
        real                                             :: long              ! Longitude (degrees)
        real                                             :: elevation         ! Elevation (m)
        real                                             :: slope             ! Slope (degrees)
        real                                             :: aspect            ! Aspect (degrees)
        real                                             :: altitude          ! Site altitude (for altitude climate adjustments) (m)
        real                                             :: wind_prob         ! Number of windthrow events in 1000 years
        real                                             :: a_sat             ! A-layer saturation capacity (volumetric)
        real                                             :: a_fc              ! A-layer field capacity (volumetric)
        real                                             :: a_pwp             ! A-layer permanent wilting point (volumetric)
        real                                             :: hum_input         ! Initial humus amount (t/ha)
        real                                             :: o_sat             ! Organic layer saturation capacity (volumetric)
        real                                             :: o_fc              ! Organic layer field capacity (volumetric)
        real                                             :: o_pwp             ! Organic permanent wilting point (volumetric)
        real                                             :: o_bd              ! Organic layer bulk density (kg/m3)
        real                                             :: a_bd              ! A-layer bulk density (kg/m3)
        real                                             :: flow              ! Moisture input from overland flow (mm)
        real                                             :: fire_wind         ! Wind speed for forced fire (m/s)
        real                                             :: fire_ffmc         ! FFMC for forced fire
        real                                             :: fire_dmc          ! DMC for forced fire
        real                                             :: A_depth           ! A-layer depth (m)
        integer                                          :: itxt              ! Soil texture (0: very coarse; 1: coarse; 2: fine)
        integer                                          :: management        ! 0: no management; 1: management
        integer                                          :: thinning          ! 0: no thinning; 1: thinning
        integer                                          :: prescribed_burn   ! 0: no prescribed burn; 1: prescribed burn
        integer                                          :: shearblading      ! 0: no shearblading; 1: shearblading
        integer                                          :: pruning      ! 0: no shearblading; 1: shearblading
        integer                                          :: sel_cutting       ! 0: no; 1: yes
        character (len = 8)                              :: planting          ! If cutting, what species will be planted?
        real                                             :: viability         ! If planting, what proportion of seedlings survive?
        integer                                          :: year_management   ! Year to conduct management
        integer                                          :: rotation_time     ! How often to conduct management (years)
        real                                             :: thin_perc         ! Percent to thin in thinning treatments (0-1)
        integer                                          :: min_age           ! Minimum age for harvesting
        real                                             :: dbh_thresh        ! DBH threshold for selective cutting (cm)
        real                                             :: total_stems     ! For site initialization (per ha)
        real, dimension(:), allocatable                  :: percent_sp     ! For site initialization
        ! character(len = 8), dimension(:), allocatable    :: species_ids    ! For site initialization
        real, dimension(:, :), allocatable               :: dbh_means     ! For site initialization
        real                                             :: ind               ! Site index (maps to row/column)
        real                                             :: row, col          ! Site row/column
        integer                                          :: ip, is                ! Looping index
        integer                                          :: k,l                 ! Looping index
        logical                                          :: match                   ! Holding value
        real :: PCTN ! Humus nitrogen content
        character (len = 3) :: unit ! which management subunit is this site located in (for active management)
        integer :: init_year ! What year is the stand being initialized in? (In case inventory has different years)

        ! Initialize properties from the sitelist file
        self%site_id = int(site_vals(1))
        self%runID = int(site_vals(2))
        altitude = site_vals(3)
        management =  0
        sel_cutting = 0
        planting = ""
        viability = 0.0
        rotation_time = 0
        dbh_thresh = 0
        fire_year = 0
        fire_day = 0
        fire_wind = 0
        fire_ffmc = 0
        fire_dmc = 0

        thinning = 0
        year_management = 0
        thin_perc = 0
        shearblading = 0

        ! Read in site file
        call read_site(self%site_id, sitename, siteregion, lat, long,          &
            elevation, slope, aspect, a_sat, a_fc, a_pwp, o_sat, o_fc, o_pwp,  &
            o_bd, a_bd, itxt, hum_input, A_depth, wind_prob,       &
            flow, row, col, ind, management, sel_cutting, planting,  &
            viability, rotation_time, dbh_thresh, thinning, year_management,   &
            thin_perc, shearblading, pruning, unit, init_year)

        ! Initialize values
        self%site_name = sitename
        self%region = siteregion
        self%latitude = lat
        self%longitude = long
        self%elevation = elevation
        self%slope = slope
        self%aspect = aspect
        self%wind_prob = wind_prob
        self%flow = flow
        self%ind = int(ind)
        self%row = int(row)
        self%col = int(col)
        self%fire_year = fire_year
        self%fire_day = fire_day
        self%fire_wind = fire_wind
        self%fire_ffmc = fire_ffmc
        self%fire_dmc = fire_dmc

        ! Standard adjustments and settings
        self%wind_prob = self%wind_prob/1000.0
        self%leaf_area_ind = 1.0
        allocate(self%lai_array(maxheight))
        self%lai_array = 0.0
        self%pc_germ = 0.0
        self%psolrem = 0.15
        self%n_adj = 0
        self%mgmtflag = 0

        ! Management
        if (management == INVALID .or. management == 0) then

            ! Set all to false
            self%management = .false.
            self%thinning = .false.
            self%shearblading = .false.
            self%sel_cutting = .false.
            self%pruning = .false.
            self%planting = ""
            self%viability = 0.0
            self%prescribed_burn = .false.

            ! Set all to 0
            self%rotation_time = 0
            self%dbh_thresh = 0.0
            self%year_management = 0
            self%thin_perc = 0.0

        else

            ! We are conducting management
            self%management = .true.

            if (year_management == INVALID) year_management = 0
            if (shearblading == INVALID) shearblading = 0
            if (prescribed_burn == INVALID) prescribed_burn = 0
            if (thinning == INVALID) thinning = 0
            if (thin_perc == RNVALID) thin_perc = 0.0
            if (sel_cutting == INVALID) sel_cutting = 0
            if (pruning == INVALID) pruning = 0
            if (rotation_time == INVALID) rotation_time = 0
            if (dbh_thresh == RNVALID) dbh_thresh = 0.0

            if (year_management == 0 .and. sel_cutting == 0 .and.              &
                shearblading == 0 .and. thinning == 0 .and.                    &
                prescribed_burn == 0) then

                ! Something wrong here - warn user
                call warning("Management turned on but no type chosen")

                ! Set all to false
                self%management = .false.
                self%thinning = .false.
                self%shearblading = .false.
                self%sel_cutting = .false.
                self%prescribed_burn = .false.
                self%pruning = .false.

                ! Set all to 0
                self%planting = ""
                self%viability = 0.0
                self%rotation_time = 0
                self%dbh_thresh = 0.0
                self%year_management = 0
                self%thin_perc = 0.0

            else

                if (shearblading == 1) then
                    if (year_management == 0) then

                        ! Something wrong - warn user
                        call warning("Shearblading turned on but no year chosen")

                        ! Set all to false
                        self%management = .false.
                        self%thinning = .false.
                        self%shearblading = .false.
                        self%sel_cutting = .false.
                        self%prescribed_burn = .false.
                        self%pruning = .false.

                        ! Set all to 0
                        self%rotation_time = 0
                        self%dbh_thresh = 0.0
                        self%year_management = 0
                        self%thin_perc = 0.0
                    else

                        ! We are shearblading
                        self%sel_cutting = .false.
                        self%prescribed_burn = .false.
                        self%thinning = .false.
                        self%shearblading = .true.
                        self%pruning = .false.
                        self%year_management = year_management
                        self%thin_perc = 0.0
                        self%rotation_time = 0
                        self%dbh_thresh = 0.0
                    end if

                else if (thinning == 1) then
                    ! We are thinning
                    if (thin_perc <= 0.0 .or. year_management == 0) then
                        ! Something wrong here - warn user
                        call warning("Thinning chosen but invalid choices")

                        ! Set all to false
                        self%management = .false.
                        self%thinning = .false.
                        self%shearblading = .false.
                        self%sel_cutting = .false.
                        self%prescribed_burn = .false.
                        self%pruning = .false.

                        ! Set all to 0
                        self%rotation_time = 0
                        self%dbh_thresh = 0.0
                        self%sel_cutting = .false.
                        self%planting = planting
                        self%viability = viability
                        self%year_management = 0
                        self%thin_perc = 0.0
                        self%rotation_time = 0
                    else

                        ! We are thinning
                        self%thinning = .true.
                        self%shearblading = .false.
                        self%sel_cutting = .false.
                        self%pruning = .false.

                        self%dbh_thresh = 0.0
                        self%rotation_time = 0
                        self%year_management = year_management
                        self%thin_perc = thin_perc
                    end if

                else if (sel_cutting == 1) then

                    ! We are selectively cutting
                    self%thinning = .false.
                    self%shearblading = .false.
                    self%sel_cutting = .true.
                    self%planting = planting
                    self%viability = viability
                    self%year_management = year_management
                    self%thin_perc = 0.0
                    self%rotation_time = rotation_time
                    self%dbh_thresh = dbh_thresh
                    self%already_cut = .false.
                    self%pruning = .false.

                else if (pruning == 1) then
                  if(year_management == 0) then 
                  
                  ! Something wrong - warn user
                  call warning("Pruning turned on but no year chosen")

                  ! Set all to false
                  self%management = .false.
                  self%thinning = .false.
                  self%shearblading = .false.
                  self%sel_cutting = .false.
                  self%prescribed_burn = .false.
                  self%pruning = .false.

                  ! Set all to 0
                  self%rotation_time = 0
                  self%dbh_thresh = 0.0
                  self%year_management = 0
                  self%thin_perc = 0.0
                  else

                  ! We are pruning
                  self%sel_cutting = .false.
                  self%prescribed_burn = .false.
                  self%thinning = .false.
                  self%shearblading = .false.
                  self%pruning = .true.
                  self%year_management = year_management
                  self%thin_perc = 0.0
                  self%rotation_time = 0
                  self%dbh_thresh = 0.0
                  end if
                end if 
            end if
        endif

        self%planting = planting
        self%unit = unit
        self%init_year = init_year

        ! Read in climate data
        call read_climate(self%site_id, tmin, tmax, prcp, cld, rh, wind)
        call read_lightning(self%site_id, strikes, strikes_sd)

        ! Adjust climate data for correct units
        self%tmin = tmin
        self%tmax = tmax
        self%precip = prcp*MM_TO_CM
        self%cld = cld/10.0
        self%rh = rh
        self%wind = wind
        self%strikes = strikes
        self%strikes_std = strikes_sd
        self%accum_tmax = 0.0
        self%accum_tmin = 0.0
        self%accum_precip = 0.0

        ! Adjust temperature and precipitation for altitude if necessary
        if (altitude /= RNVALID) then
            adjust_altitude = .true.
            self%altitude = altitude
            call adjustForAltitude(self)
        else
            adjust_altitude = .false.
            self%altitude = self%elevation
        endif

        ! Read in climate stddev data
        if (self%site_id /= INVALID) then
            if (use_climstd) then
                call read_climate_stds(self%site_id, tmin_std, tmax_std,       &
                    prcp_std, cld_std, rh_std)
                self%tmin_std = tmin_std
                self%tmax_std = tmax_std
                self%precip_std = prcp_std*MM_TO_CM
                self%cld_std = cld_std/10.0
                self%rh_std = rh_std
            end if
        endif

        ! Read in rangelist for site and add to site, also initialize plot-level
        ! soil information
        if (self%site_id /= INVALID) then
            if  (use_rangelist) then
                ! We are using a rangelist - read in
                call read_rangelist(self%site_id, range_species_ids)

                ! Add range information to site
                if (size(range_species_ids) /= 0) then
                    ! Read was successfull - attach some species
                    call attach_species(self, species_data, species_ids, sndx, &
                        range_species_ids)
                else
                    ! Read was not successful or something weird happened
                    ! Tell user
                    write(message, '(A,A,I9)') "Bad read of rangelist for ",   &
                        " site ", self%site_id
                    call warning(message)
                    self%site_id = INVALID
                endif
            else
                ! Not using rangelist
                ! All species present in this site
                call attach_species(self, species_data, species_ids, sndx)
            endif

            if (init_on) then
              call read_init_sp(self%site_id, self%total_stems, self%init_BA, self%percent_sp, init_species_ids)
              ! Add range information to site
              ! if (size(range_species_ids) /= 0) then
              !     ! Read was successfull - attach some species
              !     call attach_species(self, species_data, species_ids, sndx, &
              !         range_species_ids)
              !
              ! else
              !     ! Read was not successful or something weird happened
              !     ! Tell user
              !     write(message, '(A,A,I9)') "Bad read of rangelist for ",   &
              !         " site ", self%site_id
              !     call warning(message)
              !     self%site_id = INVALID
              ! endif
              call read_init_dbh(self%site_id, self%species_ids, self%dbh_means)

            !make sure species in init_sp are attached
            do l = 1, size(init_species_ids)
              match = .false.
              do k = 1, size(range_species_ids)
                if(init_species_ids(l) .eq. range_species_ids(k)) then
                  match = .true.
                  exit
                end if
              end do
              if(match .eq. .false.) then
                stop "Error in initialization - species in input file is not in rangelist"
              end if
            end do

          end if !init_on
        end if !site not invalid

        allocate(self%la_spec(self%num_trees, maxheight))
        self%la_spec = 0.0

        ! Now that we have the species info, initialize plots
        if (self%site_id /= INVALID) then
            self%numplots = numplots
          allocate(self%plots(self%numplots))
            if(init_on) then
              do is = 1, size(self%species_ids)
                if(self%species_ids(is) == "PICEmari") exit
              end do
              if(self%percent_sp(is) > 0.7) then
                ! Reducing nitrogen for black spruce sites
                PCTN = 0.006      ! Initial N % in humus
                 print *, "Black spruce site initialized"
              else
                PCTN = 0.008      ! Initial N % in humus
              end if
            end if
            do ip = 1, self%numplots
                call initialize_plot(self%plots(ip), size(self%species),       &
                    a_sat, a_fc, a_pwp, o_sat, o_fc, o_pwp, o_bd, a_bd, &
                    itxt, hum_input, A_depth, PCTN)
            enddo
        end if
    end subroutine initialize_site

    !:.........................................................................:

    subroutine attach_species(self, species_data, species_ids, sndx,           &
         range_species_ids)
         !
         !  Attaches correct species list for current site, depending on the
         !  rangelist for the site.
         !
         !  Record of revisions:
         !      Date       Programmer          Description of change
         !      ====       ==========          =====================
         !    01/01/12     K. Holcomb          Original Code
         !    10/10/19     A. C. Foster        Updated for seed rain so that we can
         !                                     have 'native' and 'non-native'
         !                                     species
         !    01/26/21     A. C. Foster        Moved initialize plot outside of this
         !                                       subroutine
         !

        ! Data dictionary: calling arguments
        class(SiteData),                              intent(inout) :: self              ! Site object
        type(SpeciesData),            dimension(:),   intent(inout) :: species_data      ! Array of species objects
        character(len = *), optional, dimension(:),   intent(in)    :: range_species_ids ! Unique IDs for species present this site
        integer,                      dimension(:,:), intent(inout) :: species_ids       ! Site IDs x species IDs
        integer,                                      intent(in)    :: sndx              ! Index of site in site array

        ! Data dictionary: local variables
        integer :: num_all_species   ! Number of total species in input file
        integer :: num_site_species  ! Number of total species in site
        integer :: num_range_species ! Number of total species in rangelist file
        integer :: n, nn, t          ! Looping indices

        if (present(range_species_ids)) then
            ! Total number of species in specieslist
            num_all_species = size(species_data)

            ! Number of species in rangelist file
            num_range_species = size(range_species_ids)

            ! Number of species native to this site
            num_site_species = count(range_species_ids /= 'NP')

            ! Loop through and add native species to species object
            if (num_site_species == 0) then
                allocate(self%species(0)) ! Nothing to allocate, no native species
            else
                ! Loop through and add species data for each native species
                ! Also add species id to species_id array for specific site
                t = 1
                do n = 1, num_all_species
                    do nn = 1, num_range_species
                        if (species_data(n)%unique_id ==                     &
                            range_species_ids(nn)) then
                            call appendspecies(self%species, species_data(n))
                            species_ids(sndx, t) = n
                            t = t + 1
                        endif
                    enddo
                enddo
            endif
        else
            t = 1
            ! No range list, all species in this site
            do n = 1, num_all_species
                call appendspecies(self%species, species_data(n))
                species_ids(sndx, t) = n
                t = t + 1
            enddo

        endif
        self%species%native = .true. ! Set these species to native

        if (seedrain) then

            ! Now add in non-native species if doing seed rain
            t = num_site_species + 1
            do n = 1, num_all_species
                if (.not. any(self%species%unique_id ==                      &
                    species_data(n)%unique_id)) then
                    call appendspecies(self%species, species_data(n))
                    species_ids(sndx, t) = n
                    t = t + 1
                end if
            end do
        end if

        ! Allocate the recruiting trees array and initialize
        allocate(self%recr_trees(size(self%species)))
        self%recr_trees = 0.0
        ! Count number of shrub and tree species
        self%num_trees = count(self%species(:)%form == 1)
        self%num_shrubs = count(self%species(:)%form /= 1)

    end subroutine attach_species

    !:.........................................................................:

    subroutine adjustForAltitude(self)
        !
        !  Adjusts precipitaiton and temperature data for altitude
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb          Original Code
        !

        ! Data dictionary: calling arguments
        class(SiteData), intent(inout) :: self ! Site object

        ! Data dictionary: local variables
        integer :: z ! Looping index

        ! Loop through and adjust each month
        ! Lapse rates are in degC/km and mm/km, elevation is in m
        if (self%altitude /= RNVALID .and. adjust_altitude) then
            do  z = 1, 12
                self%tmax(z) = self%tmax(z) -                                  &
                    (self%altitude - self%elevation)*self%temp_lapse_r(z)*0.01
                self%tmin(z) = self%tmin(z) -                                  &
                    (self%altitude - self%elevation)*self%temp_lapse_r(z)*0.01
                self%precip(z) = (max(self%precip(z) +                         &
                    (self%altitude - self%elevation)*self%precip_lapse_r(z)*   &
                    0.001, 0.0))
            end do
        endif

    end subroutine adjustForAltitude

    !:.........................................................................:

    subroutine delete_site(self)
        !
        !  Frees memory associated with a site
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb          Original Code
        !

        ! Data dictionary: calling arguments
        class(SiteData), intent(inout) :: self ! Site object

        ! Data dictionary: local variables
        integer :: ip ! Looping index

        ! Free memory from plots
        do ip = 1, numplots
            call delete_plot(self%plots(ip))
        enddo

        ! Deallocate site arrays
        if (allocated(self%plots)) deallocate(self%plots)
        if (allocated(self%species)) deallocate(self%species)
        if (allocated(self%recr_trees)) deallocate(self%recr_trees)
        if (allocated(self%lai_array)) deallocate(self%lai_array)
if (allocated(self%la_spec)) deallocate(self%la_spec)

    end subroutine delete_site

    !:.........................................................................:

    subroutine classify_site(self)
      !
      !  Classifies site to describe stand
      !
      !  Record of revisions:
      !      Date       Programmer          Description of change
      !      ====       ==========          =====================
      !    08/30/22     S. Sundquist          Original Code
      !
      class(SiteData), intent(inout) :: self ! Site object
      integer :: ip, it, is, i
      real, dimension(self%num_trees) :: BA ! basal area by species (m2/ha)
      real, dimension(self%num_trees) :: PBA ! percent basal area by species
      real :: sumBA !m2/ha
      character(len = 8), dimension(self%num_trees) :: species !tree species ids
      integer :: bsp, wsp, asp, bir, bal, lpn ! species locations in "species"
      integer :: mixcount
      real :: counter, total, qmeandbh

      BA = 0
      sumBA = 0
      self%classification = "UNC" !unclassified
      mixcount = 0

      ! count basal area per tree species
      i = 1
      do is = 1, size(self%species)
        if (self%species(is)%form == 1) then
          species(i) = self%species(is)%unique_id
          i = i+1
        end if
      end do

      do ip = 1, self%numplots
        do it = 1, self%plots(ip)%numtrees
          if(self%plots(ip)%trees(it)%spec_ptr%form ==1) then
            do is = 1, self%num_trees
              if(self%plots(ip)%trees(it)%spec_ptr%unique_id == species(is)) then
                BA(is) = BA(is) + CM_TO_M*CM_TO_M*0.25*pi*(self%plots(ip)%trees(it)%diam_bht**2)
                exit
              end if
            end do
          end if
        end do
      end do

      BA = BA/float(numplots)/plotsize*HEC_TO_M2

      do is = 1, self%num_trees
        sumBA = sumBA + BA(is)
      end do

      if(sumBA > 0) then 
        PBA = BA/sumBA
      else 
        PBA = 0.01 
      end if 

      bsp = -1
      wsp = -1
      asp = -1
      bir = -1
      bal = -1
      lpn = -1
      ! now classify forest type
      do is = 1, size(species)
        if(species(is) == "PICEglau") wsp = is
        if(species(is) == "PICEmari") bsp = is
        if(species(is) == "POPUtrem") asp = is
        if(species(is) == "POPUbals") bal = is
        if(species(is) == "BETUneoa") bir = is
        if(species(is) == "PINUcont") lpn = is
      end do

      ! white spruce

      
      if(wsp > 0) then 
        if(PBA(wsp) > 0.6) then
          self%classification = "WSPP"
          counter = 0.0
          total = 0.0
          do ip = 1, self%numplots
            do it = 1, self%plots(ip)%numtrees
              if(self%plots(ip)%trees(it)%spec_ptr%unique_id == "PICEglau") then
                counter = counter + 1.0
                total = total + self%plots(ip)%trees(it)%diam_bht**2
              end if
            end do
          end do
          qmeandbh = total/counter
          if(qmeandbh > 450) then
            self%classification = "WSPS"
          end if
        end if 
      end if
      if(wsp > 0 .and. bir  > 0) then 
        if((PBA(wsp) + PBA(bir) > 0.75) .and. &
        PBA(wsp) > 0.2 .and. PBA(bir) > 0.2) then
        self%classification = "WSB"
        end if
      end if 
      if (bsp > 0) then 
        if(PBA(bsp) > 0.33) then 
          self%classification = "BSP"
          if (wsp > 0) then 
            if(BA(wsp) >= 18) then
              self%classification = "WSPP"
              counter = 0.0
              total = 0.0
              do ip = 1, self%numplots
                do it = 1, self%plots(ip)%numtrees
                  if(self%plots(ip)%trees(it)%spec_ptr%unique_id == "PICEglau") then
                    counter = counter + 1.0
                    total = total + self%plots(ip)%trees(it)%diam_bht**2
                  end if
                end do
              end do
              qmeandbh = total/counter
              if(qmeandbh > 450) then
                self%classification = "WSPS"
              end if
            end if
         end if
        end if 
      end if 
      if(wsp > 0 .and. bsp > 0) then 
        if(PBA(bsp) + PBA(wsp) > 0.85 .and. BA(wsp) > 18) then 
          self%classification = "WSPP" 
          counter = 0.0
          total = 0.0
          do ip = 1, self%numplots
            do it = 1, self%plots(ip)%numtrees
              if(self%plots(ip)%trees(it)%spec_ptr%unique_id == "PICEglau") then
                counter = counter + 1.0
                total = total + self%plots(ip)%trees(it)%diam_bht**2
              end if
            end do
          end do
          qmeandbh = total/counter
          if(qmeandbh > 450) then
            self%classification = "WSPS"
          end if
        end if 
      end if 
      if (asp > 0) then 
        if (PBA(asp) > 0.6) then ! this would override bsp but that's okay... these won't co-occur for long
          self%classification = "ASP"
        end if 
      end if 
      if (bir > 0) then 
        if(PBA(bir) > 0.6) then ! this would override bsp but that's okay... these won't co-occur for long
          self%classification = "BIR"
        end if 
      end if 
      if (bir > 0 .and. asp > 0) then 
        if((PBA(asp) + PBA(bir) > 0.8) .and. &
          PBA(asp) > 0.3 .and. PBA(bir) > 0.3) then
          self%classification = "ABI"
        end if 
      end if 
      if(lpn > 0) then 
        if(PBA(lpn) > 0.6) then 
          self%classification = "LPIN"
        end if
      end if 
      if(asp > 0 .and. lpn > 0) then 
        if(PBA(lpn) > 0.3 .and. PBA(asp) > 0.3 .and. &
        PBA(lpn) + PBA(asp) > 0.8) then 
          self%classification = "APIN"
        end if 
      end if
      if(bsp > 0)  then 
        if(sumBA < 5 .and. PBA(bsp) < 0.33) then !if still unclassified... 
            self%classification = "REP"
        end if 
      end if
      if (self%classification == "UNC") then 
        do is = 1, size(species)
          if(PBA(is) > 0.10) mixcount = mixcount + 1
        end do
        if(mixcount >= 3) then
          self%classification = "MIX"
        else
        self%classification = "OTH"
        end if
      end if


    end subroutine classify_site

    !:.........................................................................:

    subroutine countldtrees(self, lt, dt)

      ! Counts live and dead trees for site in a given year to determine if a 
      ! site is salvage-eligible. This only considers growth (inc. flood, 
      ! permafrost, drought) and fire deaths

      ! By S. Sundquist 12/19/2022

      class(SiteData), intent(inout) :: self
      integer, intent(inout) :: lt, dt !Live and dead tree counters

      integer :: ip, it ! looping indices
      logical :: fire_survive, growth_survive
      logical :: active_crowning, passive_crowning ! crown fire?
      real                               :: R_a              ! Critical active rate of spread (m/min)
      real                               :: rosf_active      ! Active rate of spread (m/min)
      real                               :: CFB              ! Crown fraction burnt (0-1)
      real                               :: R_final          ! Final rate of spread (surface + crown) (m/min)
      real                               :: I_final          ! Final fire intensity (surface + crown) (kW/m)

      ! Reset counters for live and dead trees
      lt = 0
      dt = 0
      do ip = 1, self%numplots
        if(self%plots(ip)%fire == 1) then !if there's a fire...
          do it = 1, self%plots(ip)%numtrees
            if (self%plots(ip)%trees(it)%diam_bht>= 12.7) then
              call growth_survival(self%plots(ip)%trees(it),         &
                growth_survive)

              call active_passive_check(self%plots(ip)%soil,                 &
              self%plots(ip)%canopy_bh, self%plots(ip)%canopy_bd,        &
              self%plots(ip)%canopy_biom, self%plots(ip)%soil%wind_fire, &
              self%plots(ip)%soil%ffmc_fire, R_a, rosf_active, CFB,      &
              R_final, I_final, passive_crowning, active_crowning)

              call fire_survival(self%plots(ip)%trees(it),           &
              self%plots(ip)%soil%I_surf,                        &
              self%plots(ip)%soil%tau_l, active_crowning,        &
              CFB, fire_survive)

              if(growth_survive .and. fire_survive) then 
                lt = lt + 1
              else 
                dt = dt + 1
              end if         
            end if 
          end do
          else 
            do it = 1, self%plots(ip)%numtrees
              if(self%plots(ip)%trees(it)%diam_bht >= 12.7) then 
                call growth_survival(self%plots(ip)%trees(it),         &
                  growth_survive)
                if(growth_survive) then 
                  lt = lt + 1
                else 
                  dt = dt + 1
                end if      
              end if    
            end do
        end if 
      end do

      if(lt + dt == 0) lt = 1 ! avoiding error from division by 0


    end subroutine countldtrees

    !:.........................................................................:

    subroutine write_site_csv(self, site_unit)
        !
        !  Helper function for writing climate & site data to the output climate
        !   file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb          Original Code
        !
        use csv_file

        ! Data dictionary: calling arguments
        class(SiteData), intent(in) :: self      ! Site object
        integer,         intent(in) :: site_unit ! File unit for output climate

        ! Data dictionary: local variables
        real, dimension(numplots) :: aet           ! AET (cm)
        real, dimension(numplots) :: drydays       ! Drought index
        real, dimension(numplots) :: flooddays     ! Flooding index
        real, dimension(numplots) :: active        ! Active layer depth (cm)
        real, dimension(numplots) :: org           ! Organic layer depth (cm)
        real, dimension(numplots) :: availn        ! Plant-available N (kgN/ha)
        real, dimension(numplots) :: wilt_days     ! Wilting-point index
        real, dimension(numplots) :: saw0_ByFC     ! A-layer moisture scaled by field capacity
        real, dimension(numplots) :: aow0_ByMin    ! Organic-layer moisture scaled by wilting point
        real, dimension(numplots) :: saw0_BySAT    ! A-layer moisture scaled by saturation capacity
        real                      :: aet_mn        ! Average AET (cm)
        real                      :: aet_sd        ! SD of AET (cm)
        real                      :: drydays_mn    ! Average drought index
        real                      :: drydays_sd    ! SD of drought index
        real                      :: flooddays_mn  ! Average Flooding index
        real                      :: flooddays_sd  ! SD of flooding index
        real                      :: active_mn     ! Average active layer depth (cm)
        real                      :: active_sd     ! SD of active layer depth (cm)
        real                      :: org_mn        ! Average organic layer depth (cm)
        real                      :: org_sd        ! SD of organic layer depth (cm)
        real                      :: availn_mn     ! Average plant-available N (kgN/ha)
        real                      :: availn_sd     ! SD of plant-available N (tkgN/ha)
        real                      :: wilt_days_mn  ! Average wilting point index
        real                      :: wilt_days_sd  ! SD of wilting point index
        real                      :: saw0_ByFC_mn  ! Average A-layer moisture scaled by field capacity
        real                      :: saw0_ByFC_sd  ! SD of A-layer moisture scaled by field capacity
        real                      :: saw0_BySAT_mn ! Average A-layer moisture scaled by saturation capacity
        real                      :: saw0_BySAT_sd ! SD of A-layer moisture scaled by saturation capacity
        real                      :: aow0_ByMin_mn ! Average organic layer moisture scaled by wilting point
        real                      :: aow0_ByMin_sd ! SD of organic layer moisture scaled by wilting point
        integer                   :: ip            ! Looping index

        do ip = 1, numplots

            ! Read in relevant variables
            aet(ip) = self%plots(ip)%act_evap_day
            drydays(ip) = self%plots(ip)%dry_days
            flooddays(ip) = self%plots(ip)%flood_days
            active(ip) = self%plots(ip)%soil%active*M_TO_CM
            org(ip) = self%plots(ip)%soil%O_depth*M_TO_CM
            availn(ip) = self%plots(ip)%soil%avail_N*T_TO_KG
            wilt_days(ip) = self%plots(ip)%wilt_days
            saw0_ByFC(ip) = self%plots(ip)%saw0_ByFC
            saw0_BySAT(ip) = self%plots(ip)%saw0_BySAT
            aow0_ByMin(ip) = self%plots(ip)%aow0_ByMin
        end do

        ! Get mean and sd
        call stddev(aet, aet_mn, aet_sd, RNVALID)
        call stddev(drydays, drydays_mn, drydays_sd, RNVALID)
        call stddev(flooddays, flooddays_mn, flooddays_sd, RNVALID)
        call stddev(active, active_mn, active_sd, RNVALID)
        call stddev(org, org_mn, org_sd, RNVALID)
        call stddev(availn, availn_mn, availn_sd, RNVALID)
        call stddev(wilt_days, wilt_days_mn, wilt_days_sd, RNVALID)
        call stddev(saw0_ByFC, saw0_ByFC_mn, saw0_ByFC_sd, RNVALID)
        call stddev(saw0_BySAT, saw0_BySAT_mn, saw0_BySAT_sd, RNVALID)
        call stddev(aow0_ByMin, aow0_ByMin_mn, aow0_ByMin_sd, RNVALID)

        ! Write to file
        call csv_write(site_unit, self%rain, .false.)
        call csv_write(site_unit, self%pot_evap_day, .false.)
        call csv_write(site_unit, self%solar, .false.)
        call csv_write(site_unit, active_mn, .false.)
        call csv_write(site_unit, org_mn, .false.)
        call csv_write(site_unit, availn_mn, .false.)
        call csv_write(site_unit, aet_mn, .false.)
        call csv_write(site_unit, self%grow_days, .false.)
        call csv_write(site_unit, self%pc_germ, .false.)
        call csv_write(site_unit, self%deg_days, .false.)
        call csv_write(site_unit, self%summerVPD, .false.)
        call csv_write(site_unit, drydays_mn, .false.)
        call csv_write(site_unit, saw0_ByFC_mn, .false.)
        call csv_write(site_unit, saw0_BySAT_mn, .false.)
        call csv_write(site_unit, aow0_ByMin_mn, .false.)
        call csv_write(site_unit, wilt_days_mn, .false.)
        call csv_write(site_unit, flooddays_mn, .true.)

    end subroutine write_site_csv

    !:.........................................................................:

   end module Site
