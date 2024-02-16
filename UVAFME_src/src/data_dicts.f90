module dict_dat_module

!*******************************************************************************
  !
  ! This module sets up type definition for the data for dictionaries.
  !
!*******************************************************************************

    implicit none

    ! Define dictionary data type
    type dict_dat
        character(len=20) :: string ! The data that will be stored with each key
    end type dict_dat

    ! Data dictionary: constants
    integer,        parameter :: DICT_KEY_LENGTH = 20       ! The length of the keys
    type(dict_dat), parameter :: DICT_NULL = dict_dat( '' ) ! The "null" value for these data

end module dict_dat_module

module data_dicts

!*******************************************************************************
  !
  ! Helper module for dictionary data
  !
!*******************************************************************************

    use dict_dat_module
    implicit none

end module data_dicts
