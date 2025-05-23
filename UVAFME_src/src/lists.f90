module lists

!*******************************************************************************
  !
  ! This module contains subroutines for appending and deleting items
  !    from SpeciesData arrays
  !
  ! So far we have only implemented append and delete, but that is all we
  !    need for most programs
  !
!*******************************************************************************

use Species

implicit none


contains

    subroutine appendspecies(list, new_item)
        !
        !  Appends a new Species object to an array of Species Objects
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        type(SpeciesData), dimension(:), allocatable, intent(inout) :: list     ! List to append
        type(SpeciesData),                            intent(inout) :: new_item ! Item to append

        ! Data dictionary: local variables
        type(SpeciesData), dimension(:), allocatable :: temp_list ! Temporary list
        integer                                      :: n         ! Loopign index
        integer                                      :: list_len  ! Current length of list

        ! If list is empty, just use new item
        if (.not. allocated(list)) then
            allocate(list(1))
            list(1) = new_item
        else
            ! Otherwise, get current length of list and allocate temporary list
            list_len = size(list)
            allocate(temp_list(list_len))

            ! Set values in old list to new
            do n = 1, list_len
                temp_list(n) = list(n)
            enddo

            ! Reallocate input list
            deallocate(list)
            allocate(list(list_len + 1))

            ! Fill with old values plus new one
            do n = 1, list_len
                list(n) = temp_list(n)
            enddo
            list(list_len + 1) = new_item

        endif

    end subroutine appendspecies

    !:.........................................................................:

    subroutine deletespecies(list, list_index)
        !
        !  Deletes a Species object from an array of Species Objects
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        type(SpeciesData), dimension(:), allocatable, intent(inout) :: list       ! Input list
        integer,                                      intent(in)    :: list_index ! Index to delete

        ! Data dictionary: local variables
        type(SpeciesData), dimension(:), allocatable :: temp_list ! Temporary list
        integer                                      :: list_len  ! Current length of list
        integer                                      :: n         ! Looping index
        integer                                      :: ncount    ! Looping index

        ! If empty list, can't allocate anything
        if (.not. allocated(list) .or. size(list) < 1) then
            write(*, *) "Empty list, cannot remove item"
            return
        endif

        ! Get current size of list
        list_len = size(list)

        ! If outside of range, can't delete
        if (list_index > list_len .or. list_index < 0) then
            write(*, *) "List index not in range of list"
            return
        endif

        ! Allocate temporary list and fill with old values
        allocate(temp_list(list_len))
        do n = 1, list_len
            temp_list(n) = list(n)
        enddo

        ! Reallocate list
        deallocate(list)
        allocate(list(list_len - 1))

        ! Loop through and fill with old values, skipping index to delete
        ncount = 1
        do n = 1, list_len
            if (n /= list_index) then
                list(ncount) = temp_list(n)
                ncount = ncount + 1
            endif
        enddo

    end subroutine deletespecies

    !:.........................................................................:

    subroutine appendix(list, newix)
      integer, dimension(:), allocatable, intent(inout) :: list
      integer, intent(in) :: newix

      integer, dimension(:), allocatable :: temp_list
      integer :: i, l

      if (.not. allocated(list)) then
          allocate(list(1))
          list(1) = newix
      else
          ! Otherwise, get current length of list and allocate temporary list
          l = size(list)
          allocate(temp_list(l))

          ! Set values in old list to new
          do i = 1, l
              temp_list(i) = list(i)
          enddo

          ! Reallocate input list
          deallocate(list)
          allocate(list(l + 1))

          ! Fill with old values plus new one
          do i = 1, l
              list(i) = temp_list(i)
          enddo
          list(l + 1) = newix

      endif

    end subroutine appendix

    !:.........................................................................:

    subroutine deleteix(list, removeix)
      integer, dimension(:), allocatable, intent(inout) :: list
      integer, intent(in) :: removeix

      integer, dimension(:), allocatable :: temp_list
      integer :: i,j,l

      if (.not. allocated(list) .or. size(list) < 1) then
          write(*, *) "Empty list, cannot remove item"
          return
      endif

      l = size(list)

      ! If outside of range, can't delete
      if (removeix > l .or. removeix < 0) then
          write(*, *) "List index not in range of list"
          return
      endif

      ! Allocate temporary list and fill with old values
      allocate(temp_list(l))
      do i = 1, l
          temp_list(i) = list(i)
      enddo

      ! Reallocate list
      deallocate(list)
      allocate(list(l - 1))

      ! Loop through and fill with old values, skipping index to delete
      j = 1
      do i = 1, l
          if (i /= removeix) then
              list(j) = temp_list(i)
              j = j + 1
          endif
      enddo

    end subroutine deleteix

    !:.........................................................................:

end module lists
