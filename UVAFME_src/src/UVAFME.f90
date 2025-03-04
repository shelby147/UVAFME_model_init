program UVAFME

!*******************************************************************************
  !
  ! This is the main UVAFME program which runs the model and prints out the main
  ! banner and progress messages to the screen
  ! Current version: UVAFME v. 3.0 October 2018
  !
!*******************************************************************************

    use Constants
    use Parameters
    use Input
    use Output
    use Species
    use Site
    use GenusGroups
    use Model
    use data_dicts
    use dictionary
    use Landscape

    implicit none

    ! Data dictionary: constants
    character(len=8), parameter :: GENFIELD = 'genus'    ! Field name for genus output
    character(len=8), parameter :: SPECFIELD = 'species' ! Field name for species output

    ! Data dictonary: global variables
    type(SiteData),    dimension(:),    allocatable  :: sites           ! Array of site objects
    type(SpeciesData), dimension(:),    allocatable  :: species_data    ! Array of species data objects
    real,              dimension(:, :), allocatable  :: all_site_vals   ! Array of site runtime parameters
    integer,           dimension(:, :), allocatable  :: species_ids     ! Array of species ids for each site
    type(dict_struct), pointer                       :: dict            ! Dictionary of site ids and locations in landscape
    type(SiteData)                                   :: current_site    ! Site object
    type(Groups)                                     :: species_present ! Species/genus names
    type(dict_dat)                                   :: dat             ! Data for dictionary
    real                                             :: start_time      ! Start time of run
    real                                             :: total_time      ! Time for run
    real                                             :: init_time       ! Time for intitalization
    integer                                          :: sndx            ! Species array index
    integer                                          :: year            ! Calendar year
    integer                                          :: nargs           ! Number of command-line arguments
    integer                                          :: numsites        ! Number of sites to run
    integer, dimension(:), allocatable               :: tree_id         ! Number of trees initialized; each has unique id
    character(len=80)                                :: filelist        ! File list file name
    character(len=80)                                :: indc            ! Site index - maps to row/column
    character(len=80)                                :: sndxc           ! Site array index
    real, dimension(14) :: pactvmgmt ! Tracker in case management rolls over to next year 
    integer :: i


    interface

        subroutine drawBanner(numsites, species_present)
            !
            !  Writes the UVAFME banner and runtime parameters currently in use
            !
            !  Record of revisions:
            !      Date       Programmer          Description of change
            !      ====       ==========          =====================
            !    07/26/12    K. Holcomb           Original Code
            !

            use Parameters
            use GenusGroups
            implicit none

            type(Groups) ::  species_present ! Species/genus names
            integer      ::  numsites        ! Number of sites to run

        end subroutine drawBanner

        subroutine showProgress(asite)
            !
            !  Prints out current site being run and some site parameters
            !
            !  Record of revisions:
            !      Date       Programmer          Description of change
            !      ====       ==========          =====================
            !    07/26/12    K. Holcomb           Original Code
            !

            use Constants
            use Parameters
            use Site
            implicit none

            type(SiteData), intent(in) :: asite ! Site object

        end subroutine showProgress

    end interface

    ! Get command-line filelist name argument
    nargs = command_argument_count()
    if (nargs /= 1) then
        filelist = ''
    else
        call get_command_argument(1, filelist)
    endif

    ! Open and read file list & runtime file, initialize runtime parameters,
    ! open all other input files, and read and initialize litter parameters file
    call initialize_inputFiles(filelist)

    ! Read site list file and count how many sites to run
    call read_sitelist(all_site_vals)
    numsites = size(all_site_vals(:, 1))

    ! Read in species data file and initialize groups of genera and species
    call read_speciesdata(species_data)
    call initialize_genus_groups(species_present, species_data)

    ! Open output files and write headers
    call initialize_outputFiles

    ! Write runtime vars to screen
    call drawBanner(numsites, species_present)
    ! Start timing
    call cpu_time(start_time)

    ! Allocate species_ids to the number of site and species objects
    allocate(species_ids(numsites, size(species_data)))

    version: if (seedrain) then
      write(logf, *) "Model is running with seedrain on"

        ! We are considering seed rain - loop on years and sites and
        ! calculate landscape-level seed rain effects

        ! Set the random number generator seed
        call set_site_rng_seed(fixed_seed)

        ! Allocate site array with the number of sites
        allocate(sites(numsites))
        allocate(tree_id(numsites))
        tree_id = 0

        ! Loop through and initialize sites and site dictionary
        do sndx = 1, numsites

            ! Initialize current site
            call initialize_site(sites(sndx), all_site_vals(sndx, :),          &
                species_data, species_ids, sndx)

            ! Make sure the site exists and has valid climate data
            if (sites(sndx)%site_id == INVALID) then
                 write(logf, *) 'No valid site or climate data for site ',     &
                    sites(sndx)%site_name
                 write(logf, *) 'Skipping site ', sites(sndx)%site_name
                 write(*, *) '                          '

                ! Free the memory used by this site
                call delete_site(sites(sndx))
                cycle
             endif

            ! Also skip this site if no species are present
            if (size(sites(sndx)%species) == 0)  then
                write(logf, *) 'No species present in site ',                  &
                    sites(sndx)%site_id
                write(logf, *) 'Skipping site ', sites(sndx)%site_name
                 write(*, *) '             '

                ! Free the memory used by this site
                call delete_site(sites(sndx))
                 cycle
             endif

            ! Create dictionary and add sites to it
             ! key - landscape-wide index (maps to row, column)
             ! value - site array index
             write(sndxc, *) sndx
             write(indc, *) sites(sndx)%ind
             dat%string = sndxc
             if (sndx == 1) then
                call dict_create(dict, indc, dat)
            else
                call dict_add_key(dict, indc, dat)
             endif
        end do

        ! Tell user how long it took to initialize all sites
        call cpu_time(init_time)
         write(*, *) '  Initialization time : ', init_time - start_time
        write(*, '(A80)')                                                      &
     '============================================================================='
        ! Loop through each year
        do year = start_year, end_year
             do sndx = 1, numsites


                ! Initialize landscape-level parameters for this year
                sites(sndx)%n_adj = 0.0
                sites(sndx)%species(:)%fc_recr = 0.0

                 ! Calculate current weather/site variables for the year
                 call BioGeoClimate(sites(sndx), year)

                 ! Write climate and soil data
          !      if (mod(year, year_print_interval) == 0 .or.                   &
          !          year == numyears) then
                    call write_site_data(sites(sndx), year)
                    call write_soiln_data(sites(sndx), year)
          !      end if

                 ! Calculate current leaf area and light environment
                 call Canopy(sites(sndx))

                 ! Run annual growth and mortality
                 call Growth(sites(sndx), year)
                 call Mortality(sites(sndx), year)

             end do

            do sndx = 1, numsites
                ! Now seed other sites if possible
                if (any(sites(sndx)%recr_trees > 0.0)) then
                    call Seed_Rain(sites, sndx, dict, species_ids)
                end if
            end do

            do sndx = 1, numsites

                ! Average recruitment factors
                if (sites(sndx)%n_adj >= 1) then
                    sites(sndx)%species(:)%fc_recr =                           &
                        (sites(sndx)%species(:)%fc_recr)/sites(sndx)%n_adj
                end if

                ! Renew trees
                 call Renewal(sites(sndx), year, tree_id(sndx))

                 ! Print output
                 if (mod(year, year_print_interval) == 0 .or. year ==          &
                    numyears) then

                     ! Across-species attributes
                     call total_plot_values(sites(sndx), year)

                    ! Genus-and species-level attributes
                    call write_genus_or_species_data(sites(sndx),              &
                        species_present, year, SPECFIELD,                      &
                        species_present%numspecies, biom_by_s, pl_biom_by_s,   &
                        species_ids, sndx)

                    call write_genus_or_species_data(sites(sndx),               &
                        species_present, year, GENFIELD,                       &
                        species_present%numgenera, biom_by_g, pl_biom_by_g)

                     ! Tree-level attributes
                     if (tree_level_data) then
                         call write_tree_data(sites(sndx), year)
                     endif
                 end if
             end do

            ! Tell user we finished the year
            write(sitef, *) 'Finished year ', year

        end do

        do sndx = 1, numsites
              ! Free the memory used by each site
             call delete_site(sites(sndx))
        end do

        ! Deallocate arrays and dictionary
        if (allocated(sites)) deallocate(sites)
        call dict_destroy(dict)
        if (allocated(species_ids)) deallocate(species_ids)
        if(allocated(tree_id)) deallocate(tree_id)


    else if (init_on .and. active_management == .false.) then version
        write(logf, *) "Running model with initialized stands"
      ! We are initializing sites with mature trees

      allocate(tree_id(1))

      do sndx = 1, numsites
        tree_id = 0

        ! Set the random number generator seed
        call set_site_rng_seed(fixed_seed)

        call initialize_site(current_site, all_site_vals(sndx, :),        &
           species_data, species_ids, sndx)



        ! Make sure the site exists and has valid climate data
        if (current_site%site_id == INVALID) then
             write(logf, *) 'No valid site or climate data for site ',     &
                current_site%site_name
             write(logf, *) 'Skipping site ', current_site%site_name
             write(*, *) '                          '

            ! Free the memory used by this site
            call delete_site(current_site)
             cycle
         endif

         ! Also skip this site if no species are present
        if (size(current_site%species) == 0 .or.                           &
            current_site%site_id == INVALID)  then
            write(logf, *) 'No species present in site ',                  &
                current_site%site_id
            write(logf, *) 'Skipping site ', current_site%site_name
            write(*, *) '             '

            ! Free the memory used by this site
            call delete_site(current_site)
             cycle
         endif
         ! Print out current site and site vars to screen
         call showProgress(current_site)

         ! Run the model
         do year = start_year, end_year

             ! Calculate current weather/site variables for the year
             call BioGeoClimate(current_site, year)

                if (year == current_site%init_year) then 
                    call Init_Stand(current_site, year, tree_id(1))
                    call BioGeoClimate(current_site, year) ! recalculating some values

                    write(sitef, *) "Initialized stand on site", current_site%site_id
                else if (year > current_site%init_year) then 
                    ! Write climate and soil data
                   if (mod(year, year_print_interval) == 0 .or.                   &
                   year == end_year) then
                        call write_site_data(current_site, year)
                        call write_soiln_data(current_site, year)
                    end if

                    call Canopy(current_site)
                    call Growth(current_site, year)
      				call Mortality(current_site, year)
      				call Renewal(current_site, year, tree_id(1))
                    !don't run mortality or renewal on first year if init_on
      			end if 
                
            if (year >= current_site%init_year) then
                ! Print output
                if (mod(year, year_print_interval) == 0 .or.                  &
                    year == end_year) then

                    ! Across-species attributes
                    call total_plot_values(current_site, year)
                    ! Genus-and species-level attributes
                    call write_genus_or_species_data(current_site,             &
                        species_present, year, SPECFIELD,                      &
                        species_present%numspecies, biom_by_s, pl_biom_by_s,   &
                        species_ids, sndx)

                    call write_genus_or_species_data(current_site,             &
                        species_present, year, GENFIELD,                       &
                        species_present%numgenera, biom_by_g, pl_biom_by_g)

                    call write_dead_genus_or_species_data(current_site,        &
                        species_present, year, SPECFIELD,                      &
                        species_present%numspecies, dead_s, dead_ps)

                    call write_dead_genus_or_species_data(current_site,        &
                        species_present, year, GENFIELD,                       &
                        species_present%numgenera, dead_g, dead_pg)

                    ! Tree-level attributes
                    if (tree_level_data) then
                        call write_tree_data(current_site, year)
                    end if

                    if (mod(year, 5) == 0 .or. year == end_year) then
                        call write_forestry_data(current_site, species_present,      &
                        year, forestry)
                    end if 
                end if 
            end if 
        end do

        ! Free the memory used by this site and tell user we finished the
        ! site
        call delete_site(current_site)
        write(sitef, *) 'Finished site ', current_site%site_id

        ! Print current run time
        call cpu_time(total_time)
        write(*, *) '  Cumulative time : ', total_time - start_time
        write(*, '(A80)')                                                  &
    '============================================================================='

      end do
      if(allocated(tree_id)) deallocate(tree_id)


    else if(init_on .and. active_management) then version 


        write(logf, *) "Running model with initialized stands and active landscape management"
  
        call set_site_rng_seed(fixed_seed)
        allocate(sites(numsites))
        allocate(tree_id(numsites))
        tree_id = 0
  
  
        do sndx = 1, numsites
          ! thinking to split inputs by unit/ subunit - will have to make different script from A's
          call initialize_site(sites(sndx), all_site_vals(sndx, :),          &
              species_data, species_ids, sndx)
  
          ! Make sure the site exists and has valid climate data
          if (sites(sndx)%site_id == INVALID) then
               write(logf, *) 'No valid site or climate data for site ',     &
                  sites(sndx)%site_name
               write(logf, *) 'Skipping site ', sites(sndx)%site_name
               write(*, *) '                          '
  
              ! Free the memory used by this site
              call delete_site(sites(sndx))
               cycle
           endif
  
           ! Also skip this site if no species are present
           if (size(sites(sndx)%species) == 0)  then
               write(logf, *) 'No species present in site ',                  &
                   sites(sndx)%site_id
               write(logf, *) 'Skipping site ', sites(sndx)%site_name
                write(*, *) '             '
  
               ! Free the memory used by this site
               call delete_site(sites(sndx))
                cycle
            endif
          end do
  
            ! Tell user how long it took to initialize all sites
            call cpu_time(init_time)
             write(*, *) '  Initialization time : ', init_time - start_time
            write(*, '(A80)')                                                      &
          '============================================================================='
          ! Loop through each year
          do year = start_year, end_year
            ! note: all sites must run for same number of years; set in runtime
            do sndx = 1, numsites
  
                if (year .gt. sites(sndx)%init_year) then
                 call BioGeoClimate(sites(sndx), year)
  
                 ! Write climate and soil data
                   !  if (mod(year, year_print_interval) == 0 .or.                   &
                   !  year == numyears) then
                    call write_site_data(sites(sndx), year)
                    call write_soiln_data(sites(sndx), year)
                  !   end if
  
                   call Canopy(sites(sndx))
                   
                   call Growth(sites(sndx), year)
                   if (sndx == 1 .and. year .gt. start_amgmt) call Management(sites, numsites, year, pactvmgmt)
                   
                   if(sites(sndx)%mgmtflag > 0) then 
                      ! Printing before/ after forestry data for managed sites
                      call write_forestry_data(sites(sndx), species_present,      &
                      year, forestry) 
                      call Mortality(sites(sndx), year)
                      ! call write_forestry_data(sites(sndx), species_present,      &
                      ! year, forestry) 
                   else 
                      call Mortality(sites(sndx), year)
                      if (mod(year, year_print_interval) == 0 .or. year == end_year) then
                          call write_forestry_data(sites(sndx), species_present,      &
                          year, forestry)
                        end if
                   end if 
  
                   call Renewal(sites(sndx), year, tree_id(sndx))
                  
                   ! Print output
                   if (mod(year, year_print_interval) == 0 .or. year ==          &
                      numyears) then
  
                       ! Across-species attributes
                       call total_plot_values(sites(sndx), year)
  
                      ! Genus-and species-level attributes
                      call write_genus_or_species_data(sites(sndx),              &
                          species_present, year, SPECFIELD,                      &
                          species_present%numspecies, biom_by_s, pl_biom_by_s,   &
                          species_ids, sndx)
  
                      call write_genus_or_species_data(sites(sndx),               &
                          species_present, year, GENFIELD,                       &
                          species_present%numgenera, biom_by_g, pl_biom_by_g)
  
                      call write_dead_genus_or_species_data(sites(sndx),        &
                          species_present, year, SPECFIELD,                      &
                          species_present%numspecies, dead_s, dead_ps)
  
                      call write_dead_genus_or_species_data(sites(sndx),        &
                          species_present, year, GENFIELD,                       &
                          species_present%numgenera, dead_g, dead_pg)
  
  
                       ! Tree-level attributes
                       if (tree_level_data) then
                           call write_tree_data(sites(sndx), year)
                       endif
  
                   end if
                  
                  else if (year .eq. sites(sndx)%init_year) then 
                  call BioGeoClimate(sites(sndx), year)
                   call Init_Stand(sites(sndx), year, tree_id(sndx))
                   write(sitef, *) "Initialized stand on site", sites(sndx)%site_id
                   call BioGeoClimate(sites(sndx), year)
                 end if
  
  
                 end do
              ! Tell user we finished the year
              write(sitef, *) 'Finished year ', year
  
         ! Print current run time
          call cpu_time(total_time)
          write(*, *) '  Cumulative time : ', total_time - start_time, '    year  ', year
          write(*, '(A80)')                                                  &
      '============================================================================='
  
          end do
  
          do i = 1, 14
              if(pactvmgmt(i) > 0.0) then 
                  write(logf, *) "Model was not able to run active management item ", i, &
                  " on ", pactvmgmt(i) ," sites as intended"
              end if 
          end do 
  
          !Free up space
          do sndx = 1, numsites
                ! Free the memory used by each site
               call delete_site(sites(sndx))
          end do
          ! Deallocate arrays and dictionary
          if (allocated(sites)) deallocate(sites)
          if (allocated(species_ids)) deallocate(species_ids)
          if(allocated(tree_id)) deallocate(tree_id)
          
      else version 
        !No seedrain, no initialization, no active management
      ! We are not considering seed rain - loop on sites and free memory
      ! associated with each site after it is finished
        write(logf, *) "Classic UVAFME running"

        allocate(tree_id(1))

        do sndx = 1, numsites
            tree_id = 0
          ! Set the random number generator seed
          call set_site_rng_seed(fixed_seed)

           ! Initialize current site
           call initialize_site(current_site, all_site_vals(sndx, :),        &
              species_data, species_ids, sndx)

          ! Make sure the site exists and has valid climate data
          if (current_site%site_id == INVALID) then
               write(logf, *) 'No valid site or climate data for site ',     &
                  current_site%site_name
               write(logf, *) 'Skipping site ', current_site%site_name
               write(*, *) '                          '

              ! Free the memory used by this site
              call delete_site(current_site)
               cycle
           endif

           ! Also skip this site if no species are present
          if (size(current_site%species) == 0 .or.                           &
              current_site%site_id == INVALID)  then
              write(logf, *) 'No species present in site ',                  &
                  current_site%site_id
              write(logf, *) 'Skipping site ', current_site%site_name
              write(*, *) '             '

              ! Free the memory used by this site
              call delete_site(current_site)
               cycle
           endif

               ! Print out current site and site vars to screen
               call showProgress(current_site)

               ! Run the model
                 do year = start_year, end_year
               !do year = 0, 5

                   ! Calculate current weather/site variables for the year
                   call BioGeoClimate(current_site, year)

                  ! Write climate and soil data
          !                if (mod(year, year_print_interval) == 0 .or.                   &
          !                year == end_year) then
                      call write_site_data(current_site, year)
                      call write_soiln_data(current_site, year)
          !            end if

                  ! Calculate current LAI and light environment
                  call Canopy(current_site)

                  ! Run annual growth, mortality, and renewal
                  call Growth(current_site, year)

                  call Mortality(current_site, year)

                  call Renewal(current_site, year, tree_id(1))
                   ! Print output
                   if (mod(year, year_print_interval) == 0 .or.                  &
                      year == end_year) then

                      ! Across-species attributes
                      call total_plot_values(current_site, year)
                      ! Genus-and species-level attributes
                      call write_genus_or_species_data(current_site,             &
                          species_present, year, SPECFIELD,                      &
                          species_present%numspecies, biom_by_s, pl_biom_by_s,   &
                          species_ids, sndx)

                      call write_genus_or_species_data(current_site,             &
                          species_present, year, GENFIELD,                       &
                          species_present%numgenera, biom_by_g, pl_biom_by_g)

                      call write_dead_genus_or_species_data(current_site,        &
                          species_present, year, SPECFIELD,                      &
                          species_present%numspecies, dead_s, dead_ps)

                      call write_dead_genus_or_species_data(current_site,        &
                          species_present, year, GENFIELD,                       &
                          species_present%numgenera, dead_g, dead_pg)

                      ! Tree-level attributes
                     if (tree_level_data) then
                         call write_tree_data(current_site, year)
                     end if
                   end if
                   if (mod(year, 5) == 0 .or. end_year) then
                    call write_forestry_data(current_site, species_present,      &
                    year, forestry)

                  end if


          !          call write_tree_data(current_site, year)
              end do

              ! Free the memory used by this site and tell user we finished the
              ! site
              call delete_site(current_site)
              write(sitef, *) 'Finished site ', current_site%site_id

              ! Print current run time
              call cpu_time(total_time)
              write(*, *) '  Cumulative time : ', total_time - start_time
              write(*, '(A80)')                                                  &
          '============================================================================='

           end do
           if(allocated(tree_id)) deallocate(tree_id)
      end if version 

  !  end if !end initialization

     ! Close output files created
     call close_outputFiles

end program UVAFME

!:.............................................................................:

subroutine drawBanner(numsites, species_present)
    !
    !  Writes the UVAFME banner and runtime parameters currently in use
    !
    !  Record of revisions:
    !      Date       Programmer          Description of change
    !      ====       ==========          =====================
    !    07/26/12    K. Holcomb           Original Code
    !

    use Parameters
    use GenusGroups
    implicit none

    ! Data dictionary: calling arguments
    type(Groups), intent(in) :: species_present ! Species/genus names
    integer,      intent(in) :: numsites        ! Number of sites to run


    write(*, 500)                                                              &
'============================================================================='
    write(*, 500)                                                              &
'                       UVA Forest Model Enhanced                      '
    write(*, 500)                                                              &
'                           2018 Version 3.0                           '
    write(*, 500)                                                              &
'============================================================================='

    write(*, *) '  Running with parameters:'
    write(*, 400) 'Number of sites:', numsites
    write(*, 400) 'Number of years:', numyears
    write(*, 400) 'Number of plots:', numplots
    write(*, 400) 'Number of species:', species_present%numspecies
    write(*, 400) 'Maximum number of trees:', maxtrees
    write(*, 400) 'Maximum number of shrubs:', maxshrubs
    write(*, 400) 'Maximum number of cells:', maxcells*maxcells
    write(*, 401) 'Plotsize:', plotsize

    if (with_clim_change) then

        write(*, *) 'Running with climate change'
        write(*, 400) 'Duration in years:', gcm_duration

        if (linear_cc) then

            write(*, *) 'Running with linear cc'

            if (incr_or_decr_temp == 'incr' .and.                            &
                incr_or_decr_prcp == 'incr') then
                write(*, 401) 'Total tmin increase', incr_tmin_by
                write(*, 401) 'Total tmax increase', incr_tmax_by
                write(*, 401) 'Total precip increase', incr_precip_by
            else if (incr_or_decr_temp == 'decr' .and.                       &
                incr_or_decr_prcp == 'decr') then
                write(*, 401) 'Total tmin decrease', decr_tmin_by
                write(*, 401) 'Total tmax decrease', decr_tmax_by
                write(*, 401) 'Total precip decrease', decr_precip_by
            else if (incr_or_decr_temp == 'incr' .and.                       &
                incr_or_decr_prcp == 'decr') then
                write(*, 401) 'Total tmin increase', incr_tmin_by
                write(*, 401) 'Total tmax increase', incr_tmax_by
                write(*, 401) 'Total precip decrease', decr_precip_by
            else if (incr_or_decr_temp == 'decr' .and.                       &
                incr_or_decr_prcp == 'incr') then
                write(*, 401) 'Total tmin decrease', decr_tmin_by
                write(*, 401) 'Total tmax decrease', decr_tmin_by
                write(*, 401) 'Total precip increase', incr_precip_by
            end if

        else if (use_gcm) then

            write(*, *) 'Using GCM data:'
            write(*, 400) 'Model start year ', start_year
            write(*, 400) 'GCM start year ', start_gcm
            write(*, 400) 'GCM end year ', end_year
        end if
    end if

    write(*, 400) 'Printing interval in years:', year_print_interval
    write(*, 500)                                                              &
'============================================================================='
    write(*, *)

    400 format(A30, I10)
    401 format(A30, F10.3)
    402 format(A30, A)
    500 format(A80)

end subroutine drawBanner

!:.............................................................................:

subroutine showProgress(asite)
    !
    !  Prints out current site being run and some site parameters
    !
    !  Record of revisions:
    !      Date       Programmer          Description of change
    !      ====       ==========          =====================
    !    07/26/12    K. Holcomb           Original Code
    !

    use Parameters
    use Site
    implicit none

    ! Data dictionary: calling arguments
    type(SiteData), intent(in) :: asite ! Site object

    ! Data dictionary: local variables
    integer ::  num_site_species ! Number of species present at site

    ! Get number of species at site
    num_site_species = size(asite%species)

    write(*, 500) 'Running for site ',  asite%site_id, asite%site_name
    write(*,*) '                    '
    write(*, 501) 'Number of species present: ', num_site_species

    ! Get altitude adjustment (if present), and print out
    if (adjust_altitude .and. asite%altitude /= RNVALID) then
        write(*, *) '             Site altitude adjustment ', asite%altitude
    endif

    write(*,*)

    ! Write some other site parameters
    write(*, 502) asite%elevation, asite%slope, asite%aspect, asite%wind_prob

    500 format(14X, A, I10, 4X, a)
    501 format(14X, A, I8)
    502 format(7X, 'Site parameters: elevation ', F9.3, '   slope     ', F7.3, &
        /  23X, ' aspect      ', F7.3, '   wind/1000 ', F7.3)

end subroutine showProgress

!:.............................................................................:
