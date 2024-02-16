module dictionary

!*******************************************************************************
  !
  ! This module sets up types and subroutines for dictionaries:
  !    a mapping of strings to some data.
  !    Note:
  !       Use is made of a hash table. This should speed up most operations.
  !       The algorithm for determining the hashkey is taken from
  !       Kernighan and Pike: The Practice of Programming
  !
!*******************************************************************************

    use data_dicts
    use FileUtils

    implicit none

    integer, parameter, private :: HASH_SIZE  = 4993 ! Size of hashes
    integer, parameter, private :: MULTIPLIER = 31   ! Multiplier for hash

    ! Define list_dat type
    type list_dat
        character(len=DICT_KEY_LENGTH) :: key
        type(dict_dat)                 :: value
    end type list_dat

    ! Define hash_lit type
    type hash_list
        type(linked_list), pointer :: list
    end type hash_list

    ! Define dict_struct type
    type dict_struct
        private
        type(hash_list), pointer, dimension(:) :: table
    end type dict_struct

     ! Define linked_list type
    type linked_list
        type(linked_list), pointer :: next
        type(list_dat)             :: data
    end type linked_list

    ! We do not want everything to be public
    private :: list_dat
    private :: hash_list
    private :: linked_list
    private :: ll_create
    private :: ll_destroy
    private :: ll_count
    private :: ll_next
    private :: ll_insert
    private :: ll_insert_head
    private :: ll_delete_element
    private :: ll_get_data
    private :: ll_put_data
    private :: dict_get_elem
    private :: dict_hashkey

    contains

    !:.........................................................................:

    subroutine dict_create(dict, key, value)
        !
        !  Creates and initializes a dictionary.
        !
        !  Note:
        !     This verision assumes a shallow copy is enough (i.e. no
        !     pointers within the data to be stored). Also assumes the argument
        !     list does not already refer to a list
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(dict_struct), pointer    :: dict  ! Pointer to new dictionary
        character(len=*),  intent(in) :: key   ! Key for first element
        type(dict_dat),    intent(in) :: value ! Value for first element

        ! Data dictionary: local variables
        type(list_dat) :: data ! Local list
        integer        :: i    ! Looping index
        integer        :: hash ! Hash for key

        ! Ceate pointer and allocate table to HASH_SIZE
        allocate(dict)
        allocate(dict%table(HASH_SIZE))

        ! Set all pointers to null
        do i = 1, HASH_SIZE
            dict%table(i)%list => null()
        enddo

        ! Set to key, value pair
        data%key = key
        data%value = value

        ! Hash the key and create the dictionary
        hash = dict_hashkey(trim(key))
        call ll_create(dict%table(hash)%list, data)

    end subroutine dict_create

    !:.........................................................................:

    subroutine dict_destroy(dict)
        !
        !  Destroys an entire dictionary.
        !
        !  Note:
        !     This verision assumes that there are no pointers within the data
        !      that need deallocation
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(dict_struct), pointer :: dict ! Input dictionary

        ! Data dictionary: local variables
        integer :: i ! Looping index

        ! Destroy all associated lists
        do i = 1, size(dict%table)
            if (associated(dict%table(i)%list)) then
                call ll_destroy(dict%table(i)%list)
            endif
        enddo

        ! Deallocate
        deallocate(dict%table)
        deallocate(dict)

    end subroutine dict_destroy

    !:.........................................................................:

    subroutine dict_add_key(dict, key, value)
        !
        !  Adds a key to an existing dictionary
        !
        !  Note:
        !     If the key already exists, the key's value is simply replaced.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(dict_struct), pointer    :: dict  ! Input dictionary
        character(len=*),  intent(in) :: key   ! Key
        type(dict_dat),    intent(in) :: value ! Value

        ! Data dictionary: calling arguments
        type(list_dat)             :: data ! Local list
        type(linked_list), pointer :: elem ! Pointer to key (if it exists already)
        integer                    :: hash ! Hash for key

        ! Try to get element from dictionary
        elem => dict_get_elem(dict, key)

        ! Replace with new value if already exists, otherwise create new
        ! key/value pair
        if (associated(elem)) then
            write(logf, *) "Value replaced in dictionary for key ", key
            elem%data%value = value
        else
            data%key = key
            data%value = value

            ! Hash the key and insert the value, or create a list
            !!TODO!! Need to check for hash collision
            hash = dict_hashkey(trim(key))
            if (associated(dict%table(hash)%list)) then
                call ll_insert(dict%table(hash)%list, data)
            else
                call ll_create(dict%table(hash)%list, data)
            endif
        endif

    end subroutine dict_add_key

    !:.........................................................................:

    subroutine dict_delete_key(dict, key)
        !
        !  Deletes a key/value pair from a dictionary
        !
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: local variables
        type(dict_struct), pointer    :: dict ! Input dictionary
        character(len=*),  intent(in) :: key  ! Key to be removed

        ! Data dictionary: local variables
        type(linked_list), pointer :: elem ! Pointer to key/value pair
        integer                    :: hash ! Hash for key

        ! Get element
        elem => dict_get_elem(dict, key)

        ! Delete if it exists
        if (associated(elem)) then
            hash = dict_hashkey(trim(key))
            call ll_delete_element(dict%table(hash)%list, elem)
        endif

    end subroutine dict_delete_key

    !:.........................................................................:

    function dict_get_key(dict, key) result(value)
        !
        !  Gets the value belonging to a key
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(dict_struct), pointer    :: dict ! Input dictionary
        character(len=*),  intent(in) :: key  ! Key to retrieve
        type(dict_dat)                :: value ! Value

        ! Data dictionary: local variables
        type(linked_list), pointer :: elem  ! Pointer to key/value

        ! Point to key/value pair
        elem => dict_get_elem(dict, key)

        ! Return value if it exists
        if (associated(elem)) then
            value = elem%data%value
        else
            write(logf, *) "Can't find value for key", key
            value = DICT_NULL
        endif
    end function dict_get_key

    !:.........................................................................:

    function dict_has_key(dict, key) result(has)
        !
        !  Checks if the dictionary has a particular key
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(dict_struct), pointer    :: dict ! Input dictionary
        character(len=*),  intent(in) :: key  ! Key to check
        logical                       :: has  ! Does the dictionary have that key?

        ! Data dictionary: local variables
        type(linked_list), pointer :: elem ! Pointer to key/value

        ! Try to point to key/value
        elem => dict_get_elem(dict, key)

        ! Does it exist?
        has = associated(elem)

    end function dict_has_key

    !:.........................................................................:

    function dict_get_elem(dict, key) result(elem)
        !
        !  Finds the element with a particular key
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        !finds the element with a particular key
        !Arguments:
        !    dict: pointer to the dictionary
        !    key:  key to be sought

        type(dict_struct), pointer    :: dict ! Input dictionary
        character(len=*),  intent(in) :: key  ! Key to grab

        ! Data dictionary: local variables
        type(linked_list), pointer :: elem ! Pointer to list
        integer                    :: hash ! Hash for key

        ! Get the hash
        hash = dict_hashkey(trim(key))

        ! Get the list
        elem => dict%table(hash)%list

        ! Find the key/value pair
        do while(associated(elem))
            if (elem%data%key == key) then
                exit
            else
                elem => ll_next(elem)
            endif
        enddo

    end function dict_get_elem

    !:.........................................................................:

    integer function dict_hashkey(key)
        !
        !  Determines the hash value from an input string
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling argumnets
        character(len=*), intent(in) :: key ! Input string

        ! Data dictionary: local variables
        integer :: i ! Looping index

        ! Initialize to 0
        dict_hashkey = 0

        ! Loop through and update for each character in key
        do i = 1, len(key)
            dict_hashkey = MULTIPLIER*dict_hashkey + ichar(key(i:i))
        enddo

        ! Modulo to get final hash value
        dict_hashkey = 1 + modulo(dict_hashkey-1, HASH_SIZE)

    end function dict_hashkey

    !:.........................................................................:

    subroutine ll_create(list, data)
        !
        !  Creates and initializes a linked list.
        !
        !  Note:
        !     This version assumes a shallow copy is enough (i.e. no pointers
        !       within the data to be stored). It also assumes the argument list
        !      does not already refer to a list.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(linked_list), pointer    :: list ! List
        type(list_dat),    intent(in) :: data ! Data for list

        ! Allocate and initialize list
        allocate(list)
        list%next => null()
        list%data =  data

    end subroutine ll_create

    !:.........................................................................:

    subroutine ll_destroy(list)
        !
        !  Destroys an entire linked list.
        !
        !  Note:
        !     This version assumes that there are no pointers within the data
        !     that need deallocation
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(linked_list), pointer  :: list ! Input linked list

        ! Data dictionary: local variables
        type(linked_list), pointer :: current ! Current item in list
        type(linked_list), pointer :: next    ! Next item list

        ! Destroy all buckets within the list
        current => list
        do while (associated(current%next))
            next => current%next
            deallocate(current)
            current => next
        enddo

    end subroutine ll_destroy

    !:.........................................................................:

    integer function ll_count(list)
        !
        !  Counts the number of items in a linked list
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(linked_list), pointer :: list ! Input list

        ! Data dictionary: local variables
        type(linked_list), pointer :: current ! Current item in list

        ! Count all buckets in the list
        if (associated(list)) then
            ll_count = 1
            current => list
            do while (associated(current%next))
                current => current%next
                ll_count = ll_count + 1
            enddo
        else
            ll_count = 0
        endif

    end function ll_count

    !:.........................................................................:

    function ll_next(elem) result(next)
        !
        !  Returns the next element in a linked list (if any)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(linked_list), pointer :: elem ! Current element in list
        type(linked_list), pointer :: next ! Next element in list

        next => elem%next

    end function ll_next

    !:.........................................................................:

    subroutine ll_insert(elem, data)
        !
        !  Inserts a new item in a linked list after "elem"
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(linked_list), pointer    :: elem ! Element after which to insert
        type(list_dat),    intent(in) :: data ! Data for item

        ! Data dictionary: local variables
        type(linked_list), pointer :: next ! Next item in list

        ! Allocate new pointer
        allocate(next)

        ! Set the "next" to point to elem's next
        next%next => elem%next

        ! Set elem's "next" to point to new elem
        elem%next => next

        ! Set data for new item
        next%data =  data

    end subroutine ll_insert

    !:.........................................................................:

    subroutine ll_insert_head(list, data)
        !
        !  Inserts a new item in a linked list before the first element in the
        !  list
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(linked_list), pointer    :: list ! Head of linked list
        type(list_dat),    intent(in) :: data ! Data to insert

        ! Data dictionary: local variables
        type(linked_list), pointer :: elem ! New item in list

        ! Allocate and set new item
        allocate(elem)
        elem%data = data

        ! Set the new item to the current head, point to previous head
        elem%next => list
        list      => elem

    end subroutine ll_insert_head

    !:.........................................................................:

    subroutine ll_delete_element(list, elem)
        !
        !  Deletes an item from a linked list
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary:
        type(linked_list), pointer :: list ! Input linked list
        type(linked_list), pointer :: elem ! Item to delete

        ! Data dictionary: local variables
        type(linked_list), pointer :: current ! Current item in list
        type(linked_list), pointer :: prev    ! Previous item in list

        ! Deallocate the list
        if (associated(list, elem)) then
            list => elem%next
            deallocate(elem)
        else
            current => list
            prev    => list
            do while (associated(current))
                if (associated(current, elem)) then
                    prev%next => current%next
                    deallocate(current) ! Is also "elem"
                    exit
                endif
                prev    => current
                current => current%next
            enddo
        endif

    end subroutine ll_delete_element

    !:.........................................................................:

    function ll_get_data(elem) result(data)
        !
        !  Gets data stored in a linked list element
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(linked_list), pointer :: elem ! Element to get
        type(list_dat)             :: data ! Data for element

        data = elem%data

    end function ll_get_data

    !:.........................................................................:

    subroutine ll_put_data(elem, data)
        !
        !  Stores new data with a linked list element
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/26/07     Arjen Markus        Original Code
        !

        ! Data dictionary: calling arguments
        type(linked_list), pointer    :: elem ! Element to get
        type(list_dat),    intent(in) :: data ! Data to set

        elem%data = data

    end subroutine ll_put_data

    !:.........................................................................:

end module dictionary
