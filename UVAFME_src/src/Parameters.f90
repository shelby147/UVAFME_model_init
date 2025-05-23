module Parameters
!*******************************************************************************
  !
  ! This module sets parameters used throughout the model
  !
!*******************************************************************************

  use Constants
  implicit none

!:.............................................................................:

  ! Basic parameters
  integer             :: numyears            ! Number of years to run simulation (for seedrain only)
  integer             :: numplots            ! Number of plots to run per site
  real                :: plotsize            ! Area of plots (m2)
  integer             :: year_print_interval ! Interval to print output
  integer             :: maxtrees            ! Maximum number of trees per plot
  integer             :: maxcells            ! Maximum number of rowsxcolumns per plot
  integer             :: maxshrubs           ! Maximum number of shrubs per plot
  integer             :: maxheight           ! Maximum height of trees (m)

  ! Variables determining whether to use explicit seeds for RGNs
  logical             :: fixed_seed ! Fixed seed for RNG applications
  logical             :: debug      ! Are we debugging?

  ! Optional files and settings
  logical             :: fire_on           ! Are fire dynamics on?
  logical             :: adjust_altitude   ! Whether or not to adjust for altitude
  logical             :: use_rangelist     ! Whether or not to use rangelist
  logical             :: use_climstd       ! Whether or not to use stddevs adjustments for climate
  logical             :: testing           ! Whether or not to print out testing files
  logical             :: plot_level_data   ! Whether or not to print out plot-level data
  logical             :: tree_level_data   ! Whether or not to print out tree-level data
  logical             :: conds_testing     ! Whether or not to print out daily fuel conditions data
  logical             :: use_monthly_clim  ! Whether or not to read in monthly climate
  logical             :: init_on           ! Are we initializing the site with grown trees?
  logical             :: active_management ! Will stands be managed based on modeled stand conditions?
  character(len = 3)  :: management_scenario ! Business as usual or projected management?
  integer             :: management_scale  ! Multiply management vectors by this 
  logical             :: fire_testing      ! Are we forcing a fire event?
  logical             :: reg_testing       ! Whether or not to print out regeneration data 
  logical             :: mgmt_testing       ! Whether or not to print out fire and regen summary data

  ! Climate change
  logical             :: with_clim_change  ! Whether or not we are simulating climate change
  logical             :: use_gcm           ! Whether or not we are using an input climate change file
  logical             :: linear_cc         ! Whether or not we are conducting linear climate change
  real                :: incr_tmin_by      ! How much to increase minimum temperature for linear climate change (degC)
  real                :: incr_tmax_by      ! How much to increase maximum temperature for linear climate change (degC)
  real                :: incr_precip_by    ! How much to increase preciptation for linear climate change (proportion)
  real                :: decr_tmin_by      ! How much to decrease minimum temperature for linear climate change (degC)
  real                :: decr_tmax_by      ! How much to decrease maximum temperature for linear climate change (degC)
  real                :: decr_precip_by    ! How much to decrease precipitation for linear climate change (proportion)
  real                :: tmin_change       ! Annual minimum temperature change for linear climate change (degC)
  real                :: tmax_change       ! Annual maximum temperature change for linear climate change (degC)
  real                :: precip_change     ! Annual precipitation change for linear climate change (proportion)
  integer             :: start_year         ! Start year
  integer             :: end_year           ! End year
  integer             :: start_gcm          ! Year to start using annual input climate data
  integer             :: start_amgmt       ! Year when active management will start
  integer             :: gcm_duration      ! How long to simulate climate change
  character(len=4)    :: incr_or_decr_prcp ! Increasing ('incr') or deacreasing ('decr') precipitation
  character(len=4)    :: incr_or_decr_temp ! Increasing ('incr') or deacreasing ('decr') temperature

  ! Seed rain
  logical             :: seedrain          ! Simulating seed rain?
  integer             :: pixel_size        ! Pixel size for gridcell communication
  integer             :: nrows             ! Number of total rows of gridcells
  integer             :: ncols             ! Number of total columns of gridcells
  real                :: min_seedRain_prob ! Minimum seed rain probability to consider (0-1)

  ! Litter parameters
  real, dimension(LIT_LEVS, 12) :: litter_params

end module Parameters
