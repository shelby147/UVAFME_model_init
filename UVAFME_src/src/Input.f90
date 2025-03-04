module Input

!*******************************************************************************
  !
  ! This module contains subroutines to open and read input files
  ! and to initialize each site, plot, species, and tree list with
  ! initial input parameters/variables.
  !
!*******************************************************************************

    use Constants
    use Parameters
    use IO
    use FileUtils
    use Species

    implicit none

contains

    !:.........................................................................:

    subroutine initialize_inputFiles(filelist)
        !
        !  Calls subroutines to open and read filelist, runtime file,
        !   initializes runtime and litter parameters, and opens all other
        !   input files
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !    10/10/16     A. C. Foster         Added litter parameters
        !

        ! Data dictionary: calling arguments
        character(len = *), intent(in) :: filelist ! File list file name
        ! character(len = MAX_FILE) :: filename ! File name
        character(len = MAX_PATH) :: pathname ! Full path to file

        ! Read the file list file
        call read_file_list(filelist)
        ! Open runtime file
        call read_runFile
        ! Read runtime file and initialize runtime parameters
        ! Logfile (must be opened as an "input" file so other input can write
        ! to it)
        call build_pathname(outputdir, 'log.txt', pathname)
        logf = open_file(pathname)
        call initialize_parameters
        ! Open all other input files
        call open_inputFiles
        ! Read and initialize litter parameters file
        call read_litterpars(litter_params)

    end subroutine initialize_inputFiles

    !:.........................................................................:

    subroutine initialize_parameters
        !
        !  Reads runtime file, and initializes runtime parameters or sets to
        !  default values
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !    10/25/21     A. C. Foster         Deleted deprecated parameters
        !

        ! Data dictionary: local variables
        integer                   :: ios     ! I/O status
        character(len = MAX_CHAR) :: message ! Error message
        character(len = 80)       :: msg     ! I/O Error message

        ! Namelist of runtime parameters
        namelist /uvafme/ numyears, numplots, plotsize, year_print_interval,   &
                        maxtrees, maxcells, maxshrubs, fixed_seed,             &
                        debug, use_rangelist, use_climstd,                     &
                        testing, plot_level_data, tree_level_data,             &
                        with_clim_change, linear_cc, use_gcm, incr_tmin_by,    &
                       incr_tmax_by, incr_precip_by, decr_tmin_by,            &
                       decr_tmax_by, decr_precip_by, &
                       start_year, end_year, start_gcm,   &
                       incr_or_decr_prcp, incr_or_decr_temp, seedrain,        &
                        pixel_size, nrows, ncols, min_seedRain_prob, fire_on,  &
                        conds_testing, use_monthly_clim,              &
                        init_on, active_management, management_scenario, &
                        management_scale, fire_testing, reg_testing, mgmt_testing

        ! First initialize to default values
        init_on = .false.        
        active_management = .false.
        management_scenario = ""
        management_scale = 1.0

        ! Basic parameters
        numplots = 200
        plotsize = 500.0
        year_print_interval = 10
        maxtrees = 1000
        maxcells = 100
        maxshrubs = 9000
        maxheight = 60.0

        ! RNG parameters
        fixed_seed = .true.
        debug = .false.

        ! Settings and options
        fire_on = .false.
        fire_testing = .false.
        adjust_altitude = .false.
        use_rangelist = .true.
        use_monthly_clim = .false.
        use_climstd = .true.
        testing = .false.
        conds_testing = .false.
        plot_level_data = .false.
        tree_level_data = .false.
        reg_testing = .false.
        mgmt_testing = .false.

        ! Climate change
        with_clim_change = .false.
        linear_cc = .false.
        use_gcm = .false.
        incr_tmin_by = 0.00
        incr_tmax_by = 0.00
        incr_precip_by = 0.0000
        decr_tmin_by = 0.00
        decr_tmax_by = 0.00
        decr_precip_by = 0.0000
        start_year = 0
        end_year = 100
        start_gcm = 100
        start_amgmt = 2022
        numyears = end_year - start_year
        gcm_duration = end_year - start_gcm
        incr_or_decr_prcp = 'decr'
        incr_or_decr_temp = 'incr'

        ! Seed rain
        seedrain = .false.
        min_seedRain_prob = 0.001
        pixel_size = 1000
        nrows = 0
        ncols = 0

        ! Now read parameters namelist.
        ! Any variables not set in this file will take the default values above.

        if (rt_file /= INVALID) then
            read(rt_file, uvafme, iostat = ios, iomsg = msg)
            if (ios /= 0) then
                ! Problem reading file - tell user.
                write(message, '(A, I6, A)') "Error reading run time file",    &
                    ios, "IOMSG: ", msg
                call warning(message)
                call warning("Using default parameters")
            end if
        end if

        numyears = end_year - start_year
        gcm_duration = end_year - start_gcm

        ! Create climate change variables
        if (testing) then
          write(logf, *) 'Testing is on'
        end if

        if(fire_on) then 
            write(logf, *) "Fire is on"
        else 
            write(logf, *) "Fire is off"
        end if 

        if (with_clim_change) then
          write(logf, *)  "Climate change ON"
            ! Using climate change - check for which kind

            if (.not. use_gcm) then
                ! Using linear cc
                write(logf, *) "Linear climate change"
                linear_cc = .true.

                ! Check for linear cc inputs to make sure they make sense,
                ! if so, create yearly changes in temperature/precipitation
                if (gcm_duration == 0) then

                    ! Bad inputs - tell user
                    write(logf, *) "Inconsistent input for climate change: ",  &
                        "0 years duration"
                    stop 9
                else if (incr_tmin_by == 0.0 .and.                             &
                    decr_tmax_by == 0.0 .and.                                  &
                    decr_tmin_by == 0.0 .and.                                  &
                    decr_tmax_by == 0.0 .and.                                  &
                    incr_precip_by == 0.0 .and.                                &
                    decr_precip_by == 0.0) then

                    ! More bad inputs - tell user
                    write(logf, *) "Inconsistent climate change: no temp or ", &
                        "precip changes"
                    stop 9

                else if (incr_or_decr_temp == 'incr' .and.                     &
                    incr_or_decr_prcp == 'incr') then

                    ! We are increasing temperature and precipitation

                    if (incr_tmin_by <= 0.0 .and. incr_tmax_by <= 0.0          &
                            .and. incr_precip_by <= 0.0) then

                        ! All negative?
                        write(logf, *) "Inconsistent climate change: ",        &
                            "bad incr input"
                        stop 9
                    else
                        ! Calculate annual change
                        tmin_change = incr_tmin_by/(gcm_duration + 1)
                        tmax_change = incr_tmax_by/(gcm_duration + 1)
                        precip_change = incr_precip_by/(gcm_duration + 1)
                    end if

                else if (incr_or_decr_temp == 'decr' .and.                     &
                        incr_or_decr_prcp == 'decr') then

                    ! We are decreasing temperature and precipi

                    if (decr_tmin_by == 0.0 .and. decr_tmax_by == 0.0          &
                            .and. decr_precip_by == 0.0) then

                        ! We allow it to be positive or negative - but not 0
                        write(logf, *) "Inconsistent climate change: ",        &
                            "bad decr data"
                        stop 9
                    else
                        ! Input for decr should be positive, but assume user
                        ! made mistake if negative
                        if (decr_tmin_by < 0.0) then
                            write(logf, *) "Assuming decrease intended"
                            decr_tmin_by = -decr_tmin_by
                        endif
                        if (decr_tmax_by < 0.0) then
                            write(logf, *) "Assuming decrease intended"
                            decr_tmax_by = -decr_tmax_by
                        endif
                        if (decr_precip_by < 0.0) then
                            write(logf,*) "Assuming decrease intended"
                            decr_precip_by = -decr_precip_by
                        endif

                        ! Calculate annual change
                        tmin_change = -decr_tmin_by/(gcm_duration + 1)
                        tmax_change = -decr_tmax_by/(gcm_duration + 1)
                        precip_change = -decr_precip_by/(gcm_duration + 1)

                    end if

                else if (incr_or_decr_temp == 'incr' .and.                     &
                        incr_or_decr_prcp == 'decr') then

                    ! We are increasing temperature and decreasing precip
                    if (incr_tmin_by <= 0.0 .and. incr_tmax_by <= 0.0          &
                            .and. decr_precip_by == 0.0) then

                        ! Bad data
                        write(logf, *) "Bad incr or decr data"
                        stop 9
                    else
                        if (decr_precip_by < 0.0) then
                            ! Input for decr should be positive, but assume user
                            ! made mistake if negative
                            write(logf,*) "Assuming decrease intended"
                            decr_precip_by = -decr_precip_by
                        end if

                        ! Calculate annual change
                        tmin_change = incr_tmin_by/(gcm_duration + 1)
                        tmax_change = incr_tmax_by/(gcm_duration + 1)
                        precip_change = -decr_precip_by/(gcm_duration + 1)
                    end if

                else if (incr_or_decr_temp == 'decr' .and.                     &
                        incr_or_decr_prcp == 'incr') then

                    ! We are decreasing temperature and increasing precip
                    if (decr_tmin_by == 0.0 .and. decr_tmax_by == 0.0          &
                            .and. incr_precip_by <= 0.0) then

                        ! Bad data
                        write(logf, *) "Inconsistent climate change data"
                        stop 9
                    else
                        ! Input for decr should be positive, but assume user
                        ! made mistake if negative
                        if (decr_tmax_by < 0.0) then
                            write(logf,*) "Assuming decrease intended"
                            decr_tmax_by = -decr_tmax_by
                        end if
                        if (decr_tmin_by < 0.0) then
                            write(logf,*) "Asuming decrease intended"
                            decr_tmin_by = -decr_tmin_by
                        end if

                        ! Calculate annual change
                        tmin_change = -decr_tmin_by/(gcm_duration + 1)
                        tmax_change = -decr_tmax_by/(gcm_duration + 1)
                        precip_change = incr_precip_by/(gcm_duration + 1)
                    end if
                end if
                write(logf, *) "Temperature ", incr_or_decr_temp, "tmin diff: ", incr_tmin_by
                write(logf, *) "Precip ", incr_or_decr_prcp, "prcp increase: ", incr_precip_by
            else
              write(logf, *) "Climate from GCM"
                ! Climate change is from input file
                linear_cc = .false.
                tmin_change = 0.0
                tmax_change = 0.0
                precip_change = 0.0

                if (gcm_duration == 0) then
                    ! Bad data - need to check
                    write(logf, *) "Inconsistent input for climate change: ",  &
                        "0 years duration"
                    stop 9
                else
                    ! We start counting at 0
                    write(logf, *) "Start year: ", start_year, "End year: ", end_year
                    write(logf, *) "Start GCM: ", start_gcm, "GCM duration: ", gcm_duration
                endif
            endif
          else
            write(logf, *) "Climate change OFF"
        endif

        if(active_management) then
            write(logf, *) "Active management is ON. Scenario: ", management_scenario, "Grid resolution (m2)", pixel_size
            if(management_scenario /= "BAU" .and. management_scenario /= "PRJ") then
              call warning("Inconsistent active management scenario: not BAU or PRJ")
            end if
            if(management_scale < 0.0) then 
              call warning("Negative management scaler. No management will be triggered")
            end if 
          end if

        ! If we are debugging, fixed_seed automatically true
        if (debug) then
            fixed_seed = .true.
        endif

        ! Close runtime file
        close(rt_file)
    end subroutine initialize_parameters

    !:.........................................................................:

    subroutine read_litterpars(litter_pars)
        !
        !  Reads in parameters for litter decomposition
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/10/18     A. C. Foster         Original Code
        !

        ! Data dictionary: calling arguments
        real, dimension(LIT_LEVS, 12), intent(out) :: litter_pars ! Litter parameters array

        ! Data dictionary: local variables
        character(len = MAX_NLEN) :: litter_name ! Litter name
        character(len = MAX_NLEN) :: header      ! File header
        integer                   :: lc, m, i    ! Looping indices
        integer                   :: ios         ! I/O status

        ! Set to INVALID at first
        litter_pars = RNVALID

        ! We don't care about the header
        read(litterfile, '(A)') header

        ! Read litter parameters and set to litter_pars array
        lc = 0
        do while (lc <= LIT_LEVS)
            lc = lc + 1
            read(litterfile, *, iostat = ios, end = 10) litter_name,           &
                (litter_pars(lc, m), m = 2, 12)
        enddo

10      continue

        ! Set name column to id
        do i = 1, LIT_LEVS
            litter_pars(i, 1) = i
        end do

        !check if any are INVALID
        if (any(litter_pars == RNVALID)) then
            call fatal_error('Error in litter parameters file, RNVALID')
        endif

        ! Close file
        close(litterfile)

    end subroutine read_litterpars

    !:.........................................................................:

    subroutine read_sitelist(site_vals)
        !
        !  Reads sitelist file - runtime parameters for individual sites
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/14     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        real, allocatable, dimension(:,:), intent(inout) :: site_vals ! Array of site parameters

        ! Data dictionary: local variables
        integer                                       :: site_id    ! Site ID of site
        integer                                       :: numoptions ! Number of parameters to change
        character(len = MAX_LINE)                     :: row        ! Row of sitelist file
        character(len = :), allocatable, dimension(:) :: fields     ! Header broken into fields
        character(len = :), allocatable, dimension(:) :: rowvals    ! Row split into fields
        integer                                       :: numsites   ! Number of sites to run
        integer                                       :: ival       ! Value to read in
        integer                                       :: ios        ! I/O status
        integer                                       :: lc, l      ! Looping indices

        ! Count rows in sitelist file to get number of sites to run
        numsites = count_records(slist, 1)
        if (numsites <= 0) then
            call fatal_error("Error in site list file - no records")
        endif

        ! Read header to obtain the number of site variables to modify
        read(slist, '(A)') row
        call split(row, fields, ',')
        numoptions = size(fields) - 1

        ! Allocate size of site_vals to number of sites x number of variables to
        ! change
        allocate(site_vals(numsites, numoptions + 1))

        ios = 0
        lc = 0
        ! Set all site_vals values (except first column) to -999
        site_vals(:,2:) = RNVALID
        do
            ! Read each row and increment line counter
            read(slist, '(a)', end = 10) row
            lc = lc + 1

            ! Split row by comma
            call split(row, rowvals, ',')

            ! If nothing there, skip site
            if (size(rowvals) == 0) then
                cycle
            else if (len(rowvals(1)) == 0 .or. rowvals(1) == ' ') then
                cycle
            else
                ! Set first column to site_id and add to site_vals array
                read(rowvals(1), '(I9)') site_id
                site_vals(lc, 1) = site_id

                ! Grab all other site_vals data
                if (size(rowvals) > 1) then
                    do l = 2, size(rowvals)
                        if (len(rowvals(l)) == 0 .or. rowvals(l) == ' ') then
                            cycle
                        else
                            ! Check if it is a real or integer, still read in
                            ! as a real
                            if (scan(rowvals(l), '.') == 0) then
                                read(rowvals(l), '(I8)') ival
                                site_vals(lc, l) = real(ival)
                            else
                                read(rowvals(l), '(F15.7)') site_vals(lc, l)
                            endif
                        endif
                    enddo
                endif
            endif
        enddo

    10      continue

        ! Tell user how many sites we are running
        write(logf, *) 'Site data initialized. Total read in: ', numsites

    end subroutine read_sitelist

    !:.........................................................................:

    subroutine read_rangelist(site_id, species_present)
        !
        !  Reads in species rangelist for current site
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        integer,                                       intent(inout) :: site_id         ! Site ID of site
        character(len = 8), dimension(:), allocatable, intent(out)   :: species_present ! Unique IDs of species present

        ! Data dictionary: local variables
        character(len = MAX_LONG_LINE)                        :: line         ! Header
        character(len = :),        dimension(:),  allocatable :: fields       ! Header broken into fields
        character(len = 8),        dimension(:),  allocatable :: species_ids  ! All species in rangelist file
        integer,                   dimension(:),  allocatable :: species_pres ! Species present at site
        real                                                  :: lat          ! Latitude of site
        real                                                  :: long         ! Longitude of site
        integer                                               :: siteid       ! Site ID of site
        integer                                               :: max_numspec  ! Maximum number of species
        integer                                               :: lc, n        ! Looping indices
        integer                                               :: ios          ! I/O Status

        ! First line has the list of species.  Parse the line.
        read(rglist, '(a)', iostat = ios) line
        if (ios == 0) then
            call split(line, fields, ',')
            ! Number of species in the file is header length minus 3
            ! (siteid, lat, lon)
            max_numspec = size(fields) - 3
        else
            allocate(species_present(0))
            call fatal_error("Problem in rangelist file.")
            return
        endif

        ! Allocate arrays with max number of species
        allocate(species_ids(max_numspec))
        allocate(species_pres(max_numspec))
        allocate(species_present(max_numspec))

        ! Set up the list of species
        do n = 1, max_numspec
            species_ids(n) = fields(n + 3)
        enddo

        ! Set all species to not present ('NP') to start
        species_present = 'NP'

        ! Now read the file until we find the correct siteid
        lc = 0
        do
            read(rglist, *, iostat = ios, end = 10) siteid, lat, long,         &
                species_pres
            lc = lc + 1
            if (site_id == siteid) then
                ! Set to present where 1
                where (species_pres /= 0)
                    species_present = species_ids
                endwhere
                exit
            endif
        end do

10      continue

        ! Rewind the file
        rewind(rglist)

    end subroutine read_rangelist

    !:.........................................................................:

  	subroutine read_init_sp(site_id, tot_stems, init_BA, perc_spec, species_ids)
  		!reads total stems and percent species from inventory data

  		integer,                         intent(inout) :: site_id
  		real, dimension(:), allocatable, intent(out)   :: perc_spec
  		real,                            intent(out)   :: tot_stems ! trees/ha
  		real,                            intent(out)   :: init_BA ! m2/ha
  		character(len = 8), dimension(:), allocatable,                 &
  										 intent(out)   :: species_ids

      character(len = 8), dimension(:), allocatable  :: all_species_ids
  		character(len = MAX_LONG_LINE)                 :: line
  		character(len = :), dimension(:), allocatable  :: fields
  		integer                                        :: siteid
  		integer                                        :: max_numspec
  		integer                                        :: ios, i, n, lc

  		!tot_stems should be in trees/ha
  		!perc_spec is out of 1

  		!first line has the list of species.  parse the line.
  		read(initspec, '(a)', iostat = ios) line
  		if (ios .eq. 0) then
  			call split(line, fields, ',')
  			max_numspec = size(fields) - 3
  		else
  			call warning("Problem in species initialization file.")
  			return
  		endif
  		!allocate arrays with number of species from file
  		allocate(all_species_ids(max_numspec))
  		allocate(perc_spec(max_numspec))

  		!set up list of species
  		do n = 1, max_numspec
  			all_species_ids(n) = fields(n + 3)
  		end do
  		perc_spec = rnvalid; tot_stems = rnvalid

  		!now read the species list, per site
  		!do while doesn't really work quite as we'd like
  		!(we want repeat until)
  		lc = 0
  		do
  			read(initspec, *, iostat = ios, end = 10) siteid,          &
  					tot_stems, init_BA, perc_spec

  			if (site_id .eq. siteid) then
  				exit
  			end if
  		end do
  10      continue
  		rewind(initspec)

      if (site_id .ne. siteid) then
        stop "Error in initialization - site number invalid"
      end if

  		if ((tot_stems .eq. rnvalid .and. init_BA .eq. rnvalid)  .or.                                &
  				any(perc_spec .eq. rnvalid)) then

  			write(*, *) 'No species initilization data in site', site_id
  				site_id = invalid
  		end if

      ! get species actually present on initialized site
      n = 0
      do i = 1, max_numspec
        if(perc_spec(i)>0.0) n = n + 1
      end do
      allocate(species_ids(n))
      n = 1
      do i = 1, max_numspec
        if(perc_spec(i)>0.0) then
          species_ids(n) = all_species_ids(i)
          n = n + 1
        end if
      end do
      if(allocated(all_species_ids)) deallocate(all_species_ids)

  	end subroutine read_init_sp

  	!:.................................................................:

  	subroutine read_init_dbh(site_id, spec_ids, dbh_means)
  		!reads in dbh distributions for each species at site from
  		!inventory data

        integer,                             intent(inout) :: site_id
        character(len = 8), dimension(:), allocatable,                 &
                                             intent(out)   :: spec_ids
        real, dimension(:, :), allocatable, intent(out)    :: dbh_means

        character(len = 2500)                          :: line
        character(len = :), dimension(:), allocatable  :: fields, cols
        integer                                        :: siteid
        integer                                        :: max_numspec
        integer                                        :: ios, i, n, lc
        integer                                        :: s, d
        !dbh_means is in units of trees/ha

  		!first line has the list of species and dbh.  parse the line.
  		read(initdbh, '(a)', iostat = ios) line
  		if (ios .eq. 0) then
  			call split(line, fields, ',')
  			max_numspec = (size(fields) - 1)/10
  		else
  			call warning("Problem in species initialization file.")
  			return
  		endif

  		!allocate arrays with number of species from file
  		allocate(spec_ids(max_numspec))
  		allocate(dbh_means(max_numspec, 10))

  		!set up list of species
  		do n = 0, max_numspec - 1
  			call split(fields(n*10 + 2), cols, '_')
  			spec_ids(n+1) = cols(1)
  		end do

  		dbh_means = rnvalid

  		!now read the species list, per site
  		!do while doesn't really work quite as we'd like
  		!(we want repeat until)
  		lc = 0
  		do
  			read(initdbh, *, iostat = ios, end = 10) siteid,           &
  					(dbh_means(s,:), s = 1, max_numspec)

  			if (site_id .eq. siteid) then
  				exit
  			end if
  		end do
  10      continue
  		rewind(initdbh)

  		if (any(dbh_means .eq. rnvalid)) then

  			write(*, *) 'No dbh initilization data in site', site_id
  				site_id = invalid
  		end if
  	end subroutine read_init_dbh
    
    !:.........................................................................:
      subroutine read_active_management(unit, mgmt)
        character(len=3), intent(in) :: unit
        real, dimension(14), intent(out) :: mgmt
  
        ! real :: bau_wsp, bau_mixed, bau_birch, bau_thin, bau_salv ! these are ha/year
        ! real :: prj_wsp, prj_mixed, prj_birch, prj_thin, prj_salv
        ! real :: prj_bspharv, prj_bspthin, prj_bspprune, prj_bspshear
        character(len = 3) :: iunit
        character(MAX_LINE) :: header ! File header
        integer             :: ios    ! I/O status
  
        iunit = ""
  
        ! We don't care about the header
        read(amfile, '(A)') header
  
        ! Read until we find the correct site ID
        do
          ! read(amfile, *, iostat = ios, end = 10) iunit, &
          ! bau_wsp, bau_mixed, bau_birch, bau_thin, bau_salv, &
          ! prj_wsp, prj_mixed, prj_birch, prj_thin, prj_salv, &
          ! prj_bspharv, prj_bspthin, prj_bspprune, prj_bspshear
          read(amfile, *, iostat = ios, end = 10) iunit, mgmt
  
        ! Note: fortran breaks if there are empty columns in the csv (ie ,,,,, vs 0,0,0,0,0)
            if (unit == iunit) then
                exit
            endif
  
        enddo
  
        10      continue
  
        rewind(amfile)
  
        if (unit /= iunit) then
            ! Will be skipped
            call warning("Unit not found in active management data file")
        endif
  
      end subroutine read_active_management

    !:.........................................................................:

    subroutine read_climate(site_id, tmin, tmax, prcp, cld, rh, wind)
        !
        !  Reads in climate data for a site
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !    01/25/21     A. C. Foster         Merged with "extra" climate
        !                                       subroutine

        ! Data dictionary: calling arguments
        integer,                    intent(inout) :: site_id ! Site id of site
        real,    dimension(NTEMPS), intent(out)   :: tmin    ! Mean monthly minimum temperature (degC)
        real,    dimension(NTEMPS), intent(out)   :: tmax    ! Mean monthly maximum temperature (degC)
        real,    dimension(NTEMPS), intent(out)   :: prcp    ! Mean monthly precipitation (mm)
        real,    dimension(NTEMPS), intent(out)   :: cld     ! Mean monthly cloudiness (%)
        real,    dimension(NTEMPS), intent(out)   :: rh      ! Mean monthly relative humidity (%)
        real,    dimension(NTEMPS), intent(out)   :: wind    ! Mean monthly wind speed (m/s)

        ! Data dictionary: local variables
        integer                        :: siteid ! Site id of site
        real                           :: lat    ! Latitude of site
        real                           :: long   ! Longitude of site
        character(len = MAX_LONG_LINE) :: header ! File header
        integer                        :: ios    ! I/O status
        integer                        :: m      ! Looping index

        ! Set output variables to RNVALID to start
        tmin = RNVALID; tmax = RNVALID; prcp = RNVALID;
        cld = RNVALID; rh = RNVALID; wind = RNVALID

        ! Read in first climate file ---
        ! We don't care about the header
        read(cfile, '(A)') header

        ! Read until we find the correct site id
        do
            read(cfile, *, iostat = ios, end = 10) siteid, lat, long,          &
                (tmin(m), m = 1, NTEMPS), (tmax(m), m = 1, NTEMPS),            &
                (prcp(m), m = 1, NTEMPS)

            if (site_id == siteid) then
                exit
            endif

        enddo

10      continue

        ! Rewind the climate file
        rewind(cfile)

        ! Now read in the next climate file ----
        ! We don't care about the header
        read(cexfile, '(A)') header

        ! Now read the file until we find the correct siteid
        do
            read(cexfile, *, iostat = ios, end = 11) siteid, lat, long,        &
                (cld(m), m = 1, NTEMPS), (rh(m), m = 1, NTEMPS),               &
                (wind(m), m = 1, NTEMPS)

            if (site_id == siteid) then
                exit
            end if

        end do

11      continue

        ! Rewind the second climate file
        rewind(cexfile)

        if (any(tmin == RNVALID) .or. any(tmax == RNVALID) .or.                &
            any(prcp == RNVALID) .or. any(cld == RNVALID) .or.                 &
            any(wind == RNVALID) .or. any(rh == RNVALID)) then

            ! Can't find site - skip
            write(logf, *) 'No climate data for site number ', site_id

            site_id = INVALID

        endif

    end subroutine read_climate

    !:.........................................................................:

    subroutine read_climate_stds(site_id, tmin_std, tmax_std, prcp_std,        &
            cld_std, rh_std)
        !
        !  Reads in standard devation climate data for a site
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !    01/25/21     A. C. Foster         Merged with "extra" climate
        !                                       stddev subroutine

        ! Data dictionary: calling arguments
        integer,                    intent(inout) :: site_id  ! Site ID of site
        real,    dimension(NTEMPS), intent(out)   :: tmin_std ! SD of monthly minimum temperature (degC)
        real,    dimension(NTEMPS), intent(out)   :: tmax_std ! SD of monthly maximum temperature (degC)
        real,    dimension(NTEMPS), intent(out)   :: prcp_std ! SD of monthly precipitation (mm)
        real,    dimension(NTEMPS), intent(out)   :: cld_std  ! SD of monthly cloudiness (%)
        real,    dimension(NTEMPS), intent(out)   :: rh_std   ! SD of monthly relative humidity (%)

        ! Data dictionary: local variables
        integer                                :: siteid ! Site ID
        real                                   :: lat    ! Site latitude
        real                                   :: long   ! Site longitude
        character(len = MAX_LONG_LINE)         :: header ! File header
         integer                                :: ios    ! I/O status
        integer                                :: m      ! Looping index

        ! Set standard deviations to RNVALID
        tmin_std = RNVALID; tmax_std = RNVALID; prcp_std = RNVALID
        cld_std = RNVALID; rh_std = RNVALID

        ! Read first file ---
        ! We don't care about the header
        read(cstdfile, '(A)') header

        ! Read the file until we find the correct siteid
        do
            read(cstdfile, *, iostat = ios, end = 10) siteid, lat, long,       &
                (tmin_std(m), m = 1, NTEMPS), (tmax_std(m), m = 1, NTEMPS),    &
                (prcp_std(m), m = 1, NTEMPS)
            if (site_id == siteid) then
                exit
            endif
        enddo

10      continue

        ! Rewind the first file
        rewind(cstdfile)

        ! Now read in the second file
        ! We don't care about the header
        read(cexstdfile, '(A)') header

        ! Read the file until we find the correct siteid
        do
            read(cexstdfile, *, iostat = ios, end = 11) siteid, lat, long,     &
                (cld_std(m), m = 1, NTEMPS), (rh_std(m), m = 1, NTEMPS)

            if (site_id == siteid) then
                exit
            end if

        end do

11      continue

        ! Rewind the second file
        rewind(cexstdfile)

        if (any(tmin_std == RNVALID) .or. any(tmax_std == RNVALID) .or.        &
            any(prcp_std  == RNVALID) .or. any(cld_std  == RNVALID) .or.       &
            any(rh_std == RNVALID)) then

                ! Bad data - will be skipped
                write(logf, *) 'No climate stdev data for site number ', site_id
                site_id = INVALID
        endif

    end subroutine read_climate_stds

    !:.........................................................................:

    subroutine read_monthly_climate(site_id, simyear, tmean, tmin, tmax, prcp, &
            rh, wind, solar_rad, soilW, soilI, alt)
        !
        !  Reads in optional monthly forcing data
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    11/12/20     A. C. Foster        Original Code
        !
        ! Data dictionary: calling arguments
        integer,                 intent(inout) :: site_id   ! Site Id
        integer,                 intent(in)    :: simyear   ! Simulation year
        real, dimension(NTEMPS), intent(out)   :: tmean     ! Mean monthly temperature (degC)
        real, dimension(NTEMPS), intent(out)   :: tmin      ! Mean monthly minimum temperature (degC)
        real, dimension(NTEMPS), intent(out)   :: tmax      ! Mean monthly maximum temperature (degC)
        real, dimension(NTEMPS), intent(out)   :: prcp      ! Mean monthly precipitation (mm)
        real, dimension(NTEMPS), intent(out)   :: rh        ! Mean monthly relative humidity (%)
        real, dimension(NTEMPS), intent(out)   :: wind      ! Mean monthly wind speed (m/s)
        real, dimension(NTEMPS), intent(out)   :: solar_rad ! Mean monthly solar radiation
        real, dimension(NTEMPS), intent(out)   :: soilW     ! Mean monthly liquid soil moisture (volumetric)
        real, dimension(NTEMPS), intent(out)   :: soilI     ! Mean monthly frozen soil moisture (volumetric)
        real, dimension(NTEMPS), intent(out)   :: alt       ! Mean monthly active layer thickness

        ! Data dictionary: local variables
        integer                        :: siteid    ! Site ID
        real                           :: year      ! Year of simualtion
        real                           :: latitude  ! Site latitude
        real                           :: longitude ! Site longitude
        character(len = MAX_LONG_LINE) :: header    ! File header
        integer                        :: m         ! Looping index

        print *, "reading monthly climate"

        ! Set all to rnvalid initially
        tmean = RNVALID
        tmin = RNVALID
        tmax = RNVALID
        prcp = RNVALID
        rh = RNVALID
        wind = RNVALID
        solar_rad = RNVALID
        soilW = RNVALID
        soilI = RNVALID
        alt = RNVALID
        if (simyear == 0) then

            ! New site.  We don't assume numerical ordering in the site file.
            rewind(mclimfile)

            ! We don't care about the header
            read(mclimfile, '(A)') header

        endif

        do
            read(mclimfile, *, end = 10) siteid, latitude, longitude, year,    &
                (tmean(m), m = 1, NTEMPS), (tmin(m), m = 1, NTEMPS),           &
                (tmax(m), m = 1, NTEMPS), (prcp(m), m = 1, NTEMPS),            &
                (wind(m), m = 1, NTEMPS), (rh(m), m = 1, NTEMPS),              &
                (solar_rad(m), m = 1, NTEMPS), (soilW(m), m = 1, NTEMPS),      &
                (soilI(m), m = 1, NTEMPS), (alt(m), m = 1, NTEMPS)


            if (site_id == siteid .and. int(year) == simyear) then
                exit
            endif
        enddo

    10  	continue

        if (any(tmean == rnvalid) .or. any(tmin == rnvalid) .or.               &
            any(tmax == rnvalid) .or. any(prcp == rnvalid) .or.                &
            any(wind == rnvalid) .or. any(rh == rnvalid) .or.                  &
            any(solar_rad == rnvalid) .or. any(soilW == rnvalid) .or.          &
            any(soilI == rnvalid) .or. any(alt == rnvalid)) then
                write(logf, *) 'No monthly climate data for site number ',     &
                    site_id
            site_id = invalid
        endif

        if (int(year) /= simyear) then
            write(logf, *) 'Incomplete monthly climate data for site number ', &
                site_id
            site_id = invalid
        end if

    end subroutine read_monthly_climate

    !:.........................................................................:

    subroutine read_lightning(site_id, strikes, strikes_sd)
        !
        !  Reads monthly lightning strike data
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/10/19     A. C. Foster        Original Code

        ! Data dictionary: calling arguments
        integer,                    intent(inout) :: site_id    ! Site ID of site
        real,    dimension(NTEMPS), intent(out)   :: strikes    ! Mean monthly lightning (strikes/km2/day)
        real,    dimension(NTEMPS), intent(out)   :: strikes_sd ! SD of monthly lightning (strikes/km2/day)

        ! Data dictionary: local variables
        integer                        :: siteid ! Site id
        real                           :: lat    ! Latitude
        real                           :: long   ! Longitude
        character(len = MAX_LONG_LINE) :: header ! File header
        integer                        :: ios    ! I/O status
        integer                        :: m      ! Looping index

        ! Set outputs to RNVALID to start
        strikes = RNVALID; strikes_sd = RNVALID

        ! We don't care about the header
        read(lightfile, '(A)') header

        ! Read until we find the correct site id
        do
            read(lightfile, *, iostat = ios, end = 10) siteid, lat, long,      &
                (strikes(m), m = 1, NTEMPS), (strikes_sd(m), m = 1, NTEMPS)

            if (site_id == siteid) then
                exit
            endif

        enddo

10      continue

        ! Rewind the file
        rewind(lightfile)

        if (any(strikes == RNVALID) .or. any(strikes_sd == RNVALID)) then
            ! Bad data - will be skipped
            write(logf, *) 'No lightning data for site number ', site_id
            site_id = INVALID
        endif
    end subroutine read_lightning

    !:.........................................................................:

    subroutine read_gcm_climate(site_id, ygcm, startyear, tmin, tmax, prcp)
        !
        !  Reads in optional GCM climate change file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb          Original Code
        !    10/10/20     A. C. Foster        Added relative humidity

        ! Data dictionary: calling arguments
        integer,                 intent(inout) :: site_id   ! Site ID of site
        integer,                 intent(in)    :: ygcm      ! Current year of climate change
        integer,                 intent(in)    :: startyear ! Year to start climate change
        real, dimension(NTEMPS), intent(out)   :: tmin      ! Mean monthly minimum temperature (degC)
        real, dimension(NTEMPS), intent(out)   :: tmax      ! Mean monthly maximum temperature (degC)
        real, dimension(NTEMPS), intent(out)   :: prcp      ! Mean monthly precipitation (mm)
        !real, dimension(NTEMPS), intent(out)   :: rh        ! Mean monthly relative humidity (%)

        ! Data dictionary: local variables
        integer                        :: siteid ! Site ID
        real                           :: lat    ! Latitude
        real                           :: long   ! Longitude
        real                           :: year   ! Year
        character(len = MAX_LONG_LINE) :: header ! File header
        integer                        :: m      ! Looping index

        ! Set values to RNVALID to start
        tmin = RNVALID; tmax = RNVALID; prcp = RNVALID !; rh = RNVALID

      !  if (ygcm == startyear) then
            ! New site.  We don't assume numerical ordering in the site file.
            rewind(cgcmfile)

            ! We don't care about the header
            read(cgcmfile, '(A)') header

      !  endif

        ! Loop until we get the correct siteID and year
        do
            read(cgcmfile, *, end = 10) siteid, lat, long, year,               &
                (tmin(m), m = 1, NTEMPS), (tmax(m), m = 1, NTEMPS),            &
                (prcp(m), m = 1, NTEMPS) !, (rh(m), m = 1, NTEMPS)

            if (site_id == siteid .and. int(year) == ygcm) then
                exit
            endif

        enddo

10      continue

        if (any(tmin == RNVALID) .or. any(tmax == RNVALID) .or.                &
                any(prcp == RNVALID)) then
                write(logf, *) 'No climate change data for site number ',      &
                    site_id, 'Check input years for consistency'
                print *, 'Buggy climate change data for this site. Check year inputs.'
            site_id = INVALID
        endif

    end subroutine read_gcm_climate

    !:.........................................................................:

    subroutine read_gcm_lightning(site_id, ygcm, startyear, strmn)
        !
        !  Reads in optional lightning strike climate change fileC
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/01/20     A.C. Foster          Original Code
        !

        ! Data dictionary: calling arguments
        integer,                 intent(inout) :: site_id   ! Site ID of site
        integer,                 intent(in)    :: ygcm      ! Current year of climate change
        integer,                 intent(in)    :: startyear ! Start year of climate change
        real, dimension(NTEMPS), intent(out)   :: strmn     ! Mean monthly lightning (strikes/km2/day)

        ! Data dictionary: local variables
        integer                        :: siteid ! Site ID
        real                           :: lat    ! Latitude
        real                           :: long   ! Longitude
        real                           :: year   ! Year
        character(len = MAX_LONG_LINE) :: header ! File header
        integer                        :: m      ! Looping index

        ! Set to RNVALID to start
        strmn = RNVALID

        ! if (ygcm == startyear) then
            ! New site.  We don't assume numerical ordering in the site file.
            rewind(cgclightfile)

            ! We don't care about the header
            read(cgclightfile, '(A)') header
        ! endif

        ! Read until we find the correct site and year
        do
            read(cgclightfile, *, end = 10) siteid, lat, long, year,           &
                (strmn(m), m = 1, NTEMPS)

            if (site_id == siteid .and. int(year) == ygcm) then
                exit
            endif

        enddo

10      continue

        if (any(strmn == RNVALID)) then
            ! Bad data - skip site
            write(logf, *) 'No climate change lightning data for site ',       &
                ' number ', site_id
            site_id = INVALID
        endif

    end subroutine read_gcm_lightning

    !:.........................................................................:

    subroutine read_daily_climate(site_id, tmean, tmin, tmax, prcp)
        !
        !  Reads in optional daily climate file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/01/19     A.C. Foster          Original Code
        !

        ! Data dictionary: calling arguments
        integer,                        intent(inout) :: site_id ! Site ID of site
        real, dimension(DAYS_PER_YEAR), intent(out)   :: tmean   ! Mean daily temperature (degC)
        real, dimension(DAYS_PER_YEAR), intent(out)   :: tmin    ! Minimum daily temperature (degC)
        real, dimension(DAYS_PER_YEAR), intent(out)   :: tmax    ! Maximum daily temperature (degC)
        real, dimension(DAYS_PER_YEAR), intent(out)   :: prcp    ! Daily precipitation (mm)

        ! Data dictionary: local variables
        integer                        :: siteid ! Site ID of site
        real                           :: lat    ! Latitude
        real                           :: long   ! Longitude
        real                           :: year   ! Year
        real                           :: jd     ! Looping index
        real                           :: tdmean ! Mean temperature (degC)
        real                           :: tdmax  ! Max temperature (degC)
        real                           :: tdmin  ! Min temperature (degC)
        real                           :: dprcp  ! Precip (mm)
        character(len = MAX_LONG_LINE) :: header ! File header
        integer                        :: j      ! Looping index

        ! Set to RNVALID to start
        tmean = RNVALID; tmin = RNVALID; tmax = RNVALID; prcp = RNVALID

        ! New site.  We don't assume numerical ordering in the site file.
        rewind(dailyclim)

        ! We don't care about the header
        read(dailyclim, '(A)') header

        ! Loop until we find the correct site ID and jd = 1
        do
            read(dailyclim, *, end = 10) siteid, lat, long, year, jd, tdmean,  &
                tdmax, tdmin, dprcp

            if (site_id == siteid .and. int(jd) == 1) then

                tmean(1) = tdmean
                tmax(1) = tdmax
                tmin(1) = tdmin
                prcp(1) = dprcp

                ! Now read in the rest of the values
                do j = 2, DAYS_PER_YEAR
                    read(dailyclim, *, end = 10) siteid, lat, long, year, jd,  &
                        tmean(j), tmax(j), tmin(j), prcp(j)
                end do
                exit

            endif

        enddo

10      continue

        if (any(tmean == RNVALID) .or. any(tmin == RNVALID) .or.               &
            any(tmax == RNVALID) .or. any(prcp == RNVALID) ) then
                ! Bad data - will skip
                write(logf, *) 'No climate change data for site number ',      &
                    site_id
            site_id = INVALID
        endif

    end subroutine read_daily_climate

    !:.........................................................................:

    subroutine read_site(site_id, sitename, siteregion, lat, long, elevation,  &
        slope, aspect, asat, afc, apwp, osat, ofc, opwp, obd, abd, itxt,       &
        hum_input, A_depth, wprob, flow, row, col, ind,   &
        management, sel_cutting, planting, viability, rotation_time, dbh_thresh,&
        thinning,year_management, thin_perc, shearblading, pruning, unit, init_year)
    !
    !  Reads in site and soil input for a site
    !
    !  Record of revisions:
    !      Date       Programmer          Description of change
    !      ====       ==========          =====================
    !    01/01/12     K. Holcomb           Original Code
    !

        ! Data dictionary: calling arguments
        character(len = MAX_NLEN), intent(out)   :: sitename   ! Name of site
        character(len = MAX_NLEN), intent(out)   :: siteregion ! Site region
        integer,                   intent(inout) :: site_id    ! SiteID
        real,                      intent(out)   :: lat        ! Latitude (degrees)
        real,                      intent(out)   :: long       ! Longitude (degrees)
        real,                      intent(out)   :: elevation  ! Elevation (m)
        real,                      intent(out)   :: slope      ! Slope (degrees)
        real,                      intent(out)   :: aspect     ! Aspect (degrees)
        real,                      intent(out)   :: asat       ! A-layer saturation capacity (volumetric)
        real,                      intent(out)   :: afc        ! A-layer field capacity (volumetric), itxt, apwp
        real,                      intent(out)   :: apwp       ! A-layer permanent wilting poitn (volumetric)
        real,                      intent(out)   :: osat       ! Organic layer saturation capacity (volumetric)
        real,                      intent(out)   :: ofc        ! Organic layer field capacity (volumetric)
        real,                      intent(out)   :: opwp       ! Organic layer permanent wilting pint (volumetric)
        real,                      intent(out)   :: obd        ! Organic layer bulk density (kg/m3)
        real,                      intent(out)   :: abd        ! A-layer bulk density (kg/m3)
        real,                      intent(out)   :: hum_input  ! Initial humus amount (t/ha)
        real,                      intent(out)   :: A_depth    ! A-layer soil depth (m)
        real,                      intent(out)   :: wprob      ! Windthrow events in 1000 years
        real,                      intent(out)   :: flow       ! Moisture input from overland flow (mm)
        real,                      intent(out)   :: row, col   ! Row/column of site
        real,                      intent(out)   :: ind        ! Site index - maps to row/column
        integer,                   intent(out)   :: itxt       ! Soil texture (0: very coarse; 1: coarse; 2: fine)
        integer,                   intent(out)   :: management        ! 0: no management; 1: management
        integer,                   intent(out)   :: thinning          ! 0: no thinning; 1: thinning
        integer,                   intent(out)   :: shearblading      ! 0: no shearblading; 1: shearblading
        integer,                   intent(out)   :: pruning       ! 0: no pruning; 1: pruning
        integer,                   intent(out)   :: sel_cutting       ! 0: no; 1: yes
        character(len = 8),        intent(out)   :: planting          ! If cutting, what species will be planted?
        real,                      intent(out)   :: viability         ! What percent of seedlings survive annually for 5 years after planting?
        integer,                   intent(out)   :: year_management   ! Year to conduct management
        integer,                   intent(out)   :: rotation_time     ! How often to conduct management (years)
        real,                      intent(out)   :: thin_perc         ! Percent to thin in thinning treatments (0-1)
        real,                      intent(out)   :: dbh_thresh        ! DBH threshold for selective cutting (cm)
        character(len = 3), intent(out) :: unit
        integer, intent(out) ::init_year

        ! Data dictionary: local variables
        integer             :: siteid ! Site ID
        character(MAX_LINE) :: header ! File header
        integer             :: ios    ! I/O status

        ! Set siteID to INVALID
        siteid = INVALID

        ! We don't care about the header
        read(sfile, '(A)') header

        ! Read until we find the correct site ID
        do
          read(sfile, *, iostat = ios, end = 10) siteid, lat, long,          &
                sitename, siteregion, elevation, slope, aspect, asat, afc,     &
                apwp, osat, ofc, opwp, abd, obd, itxt, hum_input,              &
                A_depth, wprob, flow, row, col, ind,      &
                management, sel_cutting, planting, viability, rotation_time,   &
                dbh_thresh, thinning, year_management, thin_perc, shearblading, &
                pruning, unit, init_year 

        ! Note: fortran breaks if there are empty columns in the csv (ie ,,,,, vs 0,0,0,0,0)
            if (site_id == siteid) then
                exit
            endif

        enddo

10      continue

        rewind(sfile)

        if (siteid == INVALID .or. siteid /= site_id) then
            print *, ("Site not found in site data file: "), site_id
            ! Will be skipped
            call warning("Site not found in site data file ")
            site_id = INVALID
        endif

        if(init_year < start_year) then 
            call warning("Initialization year is less than start year")
        end if 


    end subroutine read_site

    !:.........................................................................:

    subroutine read_speciesdata(species_data)
        !
        !  Reads in species parameters and intializes an array of species
        !   objects
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        type(SpeciesData), dimension(:), allocatable, intent(out) :: species_data ! Array of species objects

        ! Data dictionary: local variables
        integer                        :: ios             ! I/O status
        integer                        :: lc              ! Looping index
        character(len = MAX_LONG_LINE) :: line            ! Header of file
        character(len = MAX_NLEN)      :: genus_name      ! Genus
        character(len = MAX_NLEN)      :: taxonomic_name  ! Scientific name
        character(len = 8)             :: unique_id       ! Unique 8-char ID of species
        character(len = MAX_NLEN)      :: common_name     ! Common name
        integer                        :: leaf_type       ! 0: deciduous; 1: evergreen
        integer                        :: layer_type      ! 0: no layering; 1: layering
        integer                        :: genus_sort      ! ID in genus list
        integer                        :: genus_id        ! ID in species/genus list
        integer                        :: form            ! 1: tree; 2: tree-like shrub;
                                                          ! 3: erect shrub; 4: prostrate shrub
        integer                        :: shade_tol       ! Relative shade tolerance (1-5; 5 = least tolerant)
        integer                        :: lownutr_tol     ! Relative low nutrient tolerance (0-3; 3 = least tolerant)
        integer                        :: stress_tol      ! Relative stress tolerance (1-5; 5 = least tolerant)
        integer                        :: age_tol         ! Probability of surviving to max age (1-3; 3 = least likely)
        integer                        :: dry_tol         ! Relative drought tolerance (1-6; 6 = least tolerant)
        integer                        :: flood_tol       ! Relative flood tolernace (1-7; 7 = least tolerant)
        integer                        :: perm_tol        ! Relative permafrost tolerane (1-2; 2 = least tolerant)
        integer                        :: org_tol         ! Relative ability to germinate on deep soil (1-3; 3 = least able)
        integer                        :: recr_age        ! Minimum age for reproduction
        integer                        :: fire_regen      ! Relative ability to reproduce following fire
                                                          ! (1-2: fire beneficial; 3: 4-6 fire detrimental)
        integer                        :: litter_class    ! Litter class
        integer                        :: host            ! Pests: 0: None 1: Aspen Leaf Miner
        integer                        :: plant_dens      ! Post-harvest recommended p  lanting density
        real                           :: max_age         ! Average maximum age (years)
        real                           :: max_diam        ! Average maximum DBH (cm)
        real                           :: max_ht          ! Average maximum height (m)
        real                           :: rootdepth       ! Average rooting depth (m)
        real                           :: wood_bulk_dens  ! Wood bulk density (tonnes/m3)
        real                           :: leafdiam_a      ! Scalar for leaf area-diameter relationship
        real                           :: leafdiam_exp    ! Exponent for leaf area-diameter relationship
        real                           :: leafarea_c      ! Scalar for leaf area:biomass relationship
        real                           :: deg_day_min     ! Minimum growing degree-days (>5degC)
        real                           :: deg_day_opt     ! Optimum growing degree-days (>5degC)
        real                           :: deg_day_max     ! Maximum growing degree-days (>5degC)
        real                           :: seed_surv       ! Proportion of seebank lost annually (0-1)
        real                           :: seedling_surv   ! Proportion of seedling bank lost annually (0-1)
        real                           :: invader         ! Seed rain from outside the plot (seeds/m2)
        real                           :: bark_thick      ! Bark thickness (cm bark/cm DBH)
        real                           :: F_i             ! Fire scorch parameter (Thonicke et al. 2010)
        real                           :: seed_num        ! Seed rain from within the plot (seeds/m2)
        real                           :: sprout_num      ! Regeneration from sprouting (sprouts/m2)
        real                           :: s, g            ! Growth parameters
        real                           :: beta            ! Growth parameter
        real                           :: dbh_min         ! Mininum DBH growth (cm)
        real                           :: min_recr_dbh    ! Minimum DBH for reproduction (cm)
        real                           :: seedRain_a      ! Seed rain parameter
        real                           :: seedRain_b      ! Seed rain parameter
        logical                        :: conifer         ! Is species a conifer?
        logical                        :: layering        ! Can species reproduce by layering?
        integer                        :: num_all_species ! Total number of species in file

        ! Count the records to check how many species there will be
        num_all_species = count_records(splist, 1)
        if (num_all_species < 0) then
            call fatal_error('Error in species-list file - no rows')
        endif

        ! Allocate output array to the number of species
        allocate(species_data(num_all_species))

        ! Read and discard the first line
        read(splist, '(A500)', iostat = ios) line
        if (ios /= 0) then
            call fatal_error('Unable to read the specieslist file')
        endif

        lc = 0
        ! Read the file
        do while (lc < num_all_species)
            lc = lc + 1
            read(splist, *, end = 100)                                         &
            genus_id,                                                          &
            genus_name,                                                        &
            genus_sort,                                                        &
            taxonomic_name,                                                    &
            common_name,                                                       &
            form,                                                              &
            max_age,                                                           &
            max_diam,                                                          &
            max_ht,                                                            &
            rootdepth,                                                         &
            s,                                                                 &
            g,                                                                 &
            beta,                                                              &
            wood_bulk_dens,                                                    &
            leafdiam_a,                                                        &
            leafarea_c,                                                        &
            deg_day_min,                                                       &
            deg_day_opt,                                                       &
            deg_day_max,                                                       &
            shade_tol,                                                         &
            dry_tol,                                                           &
            flood_tol,                                                         &
            perm_tol,                                                          &
            lownutr_tol,                                                       &
            bark_thick,                                                        &
            F_i,                                                               &
            fire_regen,                                                        &
            stress_tol,                                                        &
            age_tol,                                                           &
            dbh_min,                                                           &
            leaf_type,                                                         &
            litter_class,                                                      &
            invader,                                                           &
            seed_num,                                                          &
            sprout_num,                                                        &
            layer_type,                                                        &
            org_tol,                                                           &
            recr_age,                                                          &
            min_recr_dbh,                                                      &
            seedRain_a,                                                        &
            seedRain_b,                                                        &
            seed_surv,                                                         &
            seedling_surv,                                                     &
            unique_id,                                                         &
            host,                                                         &
            plant_dens

            ! Set conifer and layering
            if (leaf_type == 0) then
                conifer = .false.
            else
                conifer = .true.
            endif

            if (layer_type == 0) then
                layering = .false.
            else
                layering = .true.
            end if
            ! Initialize the species object
            call initialize_species(species_data(lc), genus_name,              &
                taxonomic_name, unique_id, common_name, form, shade_tol,       &
                lownutr_tol, stress_tol, age_tol, dry_tol, flood_tol,          &
                perm_tol, org_tol, bark_thick, F_i, fire_regen, max_age,       &
                max_diam, max_ht, wood_bulk_dens, rootdepth, leafdiam_a,       &
                leafarea_c, deg_day_min, deg_day_opt, deg_day_max,             &
                seedling_surv, invader, seed_num, sprout_num, layering,        &
                seed_surv, s, g, beta, conifer, litter_class, recr_age,        &
                dbh_min, min_recr_dbh, seedRain_a, seedRain_b, host, plant_dens)
        enddo



100     continue

        ! Tell user how many species we read in
        write(logf, *) 'Species data initialized. total read in: ',            &
            num_all_species

    end subroutine read_speciesdata

    !:.........................................................................:

end module Input
