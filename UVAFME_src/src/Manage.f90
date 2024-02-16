module Manage

!*******************************************************************************
  !
  ! This contains the management routines of UVAFME:
  !
  ! Harvest
  !
  ! Selective harvest
  !
  ! Thin
  !
  ! Shearblade
  !
  ! Salvage
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

    implicit none

contains


    subroutine Thin(site, thin_perc)
      ! Data dictionary: calling argumnets
      type(SiteData),        intent(inout) :: site ! Site object
      real, intent(in) :: thin_perc ! What percentage of stems to remove

      ! Data dictionary: local variables
      integer, dimension(:), allocatable :: dbh_ind          ! Index of sorted tree DBH array
      real,    dimension(:), allocatable :: dbh              ! Tree DBH (cm)
      real                               :: totbiomC         ! Plot-wide biomass (tC)
      real                               :: NPP_loss         ! Loss of NPP from mortality
      real                               :: biomC            ! Total tree biomass (tC)
      real                               :: leaf_bm          ! Leaf biomass (tC)
      real                               :: bcr              ! Root biomass (tC)
      integer                            :: it, ip, iu           ! Looping indices
      integer                            :: dt               ! Counter for dead trees
      integer                            :: lt               ! Counter for live trees
      integer                            :: ind              ! Tree array index
      integer                            :: thin_num         ! Number of trees to thin
      integer :: numtreeslarge !Number of trees >5cm dbh... not counting smaller than that in what I'm leaving because then everything becomes seedlings
      integer                            :: is               ! Species index
      integer                            :: trow             ! Row of tree
      integer                            :: tcol             ! Column of tree
      integer                            :: lc               ! Litter class
      integer                            :: i      ! Looping index

      numtreeslarge = 0

      ! We are thinning the stand (DOD)
      do ip = 1, site%numplots
        ! Set fire and wind to 0
        site%plots(ip)%fire = 0
        site%plots(ip)%wind = 0

                    ! Initialize counters for dead and live trees
        dt = 0
        lt = 0
        
        hntrees: if (site%plots(ip)%numtrees > 0) then


            ! Allocate dbh arrays
            allocate(dbh_ind(site%plots(ip)%numtrees))
            allocate(dbh(site%plots(ip)%numtrees))

            do it = 1, site%plots(ip)%numtrees
                if (site%plots(ip)%trees(it)%form == 1) then
                    ! Tree - add to list
                    dbh_ind(it) = it
                    dbh(it) = site%plots(ip)%trees(it)%diam_bht
                    if(site%plots(ip)%trees(it)%diam_bht > 5.0) numtreeslarge = numtreeslarge + 1 
                end if
            end do

            if (site%plots(ip)%numtrees > 1) then
                ! Sort the array so we can harvest from largest to
                ! smallest
                call merge_index(dbh, dbh_ind)
            end if

            ! Get number of trees to thin
            ! thin_num = int(float(site%plots(ip)%numtrees)*             &
            !     site%thin_perc)
            thin_num = int(float(numtreeslarge)*             &
                site%thin_perc)

            thinloop: do it = 1, site%plots(ip)%numtrees

                ! Index of tree to thin from tree object array
                ind = dbh_ind(it)

                ! Get species index
                is = site%plots(ip)%trees(it)%species_index

                ! Get leaf biomass
                call leaf_biomass_c(site%plots(ip)%trees(it))
                leaf_bm = site%plots(ip)%trees(it)%leaf_bm

               thinned: if (it > thin_num) then

                    ! Tree is large enough to keep

                    ! We have already thinned the appropriate amount
                    ! Tree survives

                    lt = lt + 1 ! Increment live tree counter

                    ! Increment tree age
                    site%plots(ip)%trees(it)%tree_age=                &
                        site%plots(ip)%trees(it)%tree_age + 1

                    ! Copy tree to new location in array
                    call copy_tree(site%plots(ip)%trees(lt),           &
                        site%plots(ip)%trees(it))

                    if (site%plots(ip)%trees(it)%form > 0 .and.       &
                        site%plots(ip)%trees(it)%forska_ht < 1.83) then
                        ! Add up shrub live foliage and fine twigs
                        site%plots(ip)%soil%shrubLitter =              &
                            site%plots(ip)%soil%shrubLitter + leaf_bm + &
                            site%plots(ip)%trees(it)%branchC*         &
                            PERC_BRANCHES(1)

                    end if

                    ! Calculate litterfall
                    ! Here we convert to dry biomass because soiln
                    ! subroutine calculates weight loss, not C loss
                    if (site%species(is)%conifer) then

                        ! Conifers don't drop all their needles
                        lc = site%species(is)%litter_class
                        site%plots(ip)%soil%litter(lc) =               &
                            site%plots(ip)%soil%litter(lc) +           &
                            leaf_bm*(CON_LEAF_RATIO)/B_TO_C
                        NPP_loss = NPP_loss + leaf_bm*CON_LEAF_RATIO

                    else

                        lc = site%species(is)%litter_class
                        site%plots(ip)%soil%litter(lc) =               &
                            site%plots(ip)%soil%litter(lc) +           &
                            leaf_bm/B_TO_C
                        NPP_loss = NPP_loss + leaf_bm

                    end if

                else thinned
                    ! We harvest the tree

                    dt = dt + 1 ! Increment dead tree counter

                    ! Copy to list of dead trees for output
                    call copy_tree(site%plots(ip)%deadtrees(dt),       &
                        site%plots(ip)%trees(it))

                    ! Set cells of that tree to unfilled
                    trow = site%plots(ip)%trees(it)%row
                    tcol = site%plots(ip)%trees(it)%col
                    site%plots(ip)%cells(trow, tcol) = 0

                    ! Add biomass to array of dead biomass
                    site%plots(ip)%d_type(IHARV) =                     &
                        site%plots(ip)%d_type(IHARV) +                 &
                        site%plots(ip)%trees(it)%biomC + leaf_bm

                    ! Add roots to litter (we are removing the rest
                    ! of the slash off the plot)
                    bcr = site%plots(ip)%trees(it)%rootC
                    site%plots(ip)%soil%litter(IROOT) =                &
                        site%plots(ip)%soil%litter(IROOT) + bcr/B_TO_C

                    ! Get woody biomass
                    biomC = site%plots(ip)%trees(it)%biomC +          &
                        site%plots(ip)%trees(it)%rootC

                    ! Decrease NPP
                    if (site%species(is)%conifer) then
                        NPP_loss = NPP_loss + biomC + leaf_bm
                    else
                        NPP_loss = NPP_loss + biomC + leaf_bm
                    end if

                end if thinned

            end do thinloop

            ! Deallocate dbh arrays
            if (allocated(dbh)) deallocate(dbh)
            if (allocated(dbh_ind)) deallocate(dbh_ind)

        end if hntrees

        ! Set number of live and dead trees
        site%plots(ip)%numtrees = lt
        if(site%site_id .eq. 183100)   print *, lt
        site%plots(ip)%num_dead = dt
      end do

    end subroutine Thin



    ! subroutine Thin(site)
    !   ! Data dictionary: calling argumnets
    !   type(SiteData),        intent(inout) :: site ! Site object

    !   ! Data dictionary: local variables
    !   integer, dimension(:), allocatable :: dbh_ind          ! Index of sorted tree DBH array
    !   real,    dimension(:), allocatable :: dbh              ! Tree DBH (cm)
    !   real                               :: totbiomC         ! Plot-wide biomass (tC)
    !   real                               :: NPP_loss         ! Loss of NPP from mortality
    !   real                               :: biomC            ! Total tree biomass (tC)
    !   real                               :: leaf_bm          ! Leaf biomass (tC)
    !   real                               :: bcr              ! Root biomass (tC)
    !   integer                            :: it, ip, iu           ! Looping indices
    !   integer                            :: dt               ! Counter for dead trees
    !   integer                            :: lt               ! Counter for live trees
    !   integer                            :: ind              ! Tree array index
    !   integer                            :: thin_num         ! Number of trees to thin
    !   integer                            :: is               ! Species index
    !   integer                            :: trow             ! Row of tree
    !   integer                            :: tcol             ! Column of tree
    !   integer                            :: lc               ! Litter class
    !   integer                            :: i      ! Looping index
    !   real :: thin_thresh ! All trees below this dbh (cm) will be thinned

    !     thin_thresh = 12.7 !5 in 

    !   ! We are thinning the stand (DOD)
    !   do ip = 1, site%numplots
    !     ! Set fire and wind to 0
    !     site%plots(ip)%fire = 0
    !     site%plots(ip)%wind = 0

    !     hntrees: if (site%plots(ip)%numtrees > 0) then

    !         ! Initialize counters for dead and live trees
    !         dt = 0
    !         lt = 0

    !         ! Allocate dbh arrays
    !         ! allocate(dbh_ind(site%plots(ip)%numtrees))
    !         ! allocate(dbh(site%plots(ip)%numtrees))

    !         ! do it = 1, site%plots(ip)%numtrees
    !         !     if (site%plots(ip)%trees(it)%form == 1) then
    !         !         ! Tree - add to list
    !         !         dbh_ind(it) = it
    !         !         dbh(it) = site%plots(ip)%trees(it)%diam_bht
    !         !     end if
    !         ! end do

    !         ! if (site%plots(ip)%numtrees > 1) then
    !         !     ! Sort the array so we can harvest from largest to
    !         !     ! smallest
    !         !     call merge_index(dbh, dbh_ind)
    !         ! end if

    !         ! Get number of trees to thin
    !         ! thin_num = int(float(site%plots(ip)%numtrees)*             &
    !         !     site%thin_perc)

    !         thinloop: do it = 1, site%plots(ip)%numtrees

    !             ! Index of tree to thin from tree object array
    !             ! ind = dbh_ind(it)

    !             ! Get species index
    !             is = site%plots(ip)%trees(it)%species_index

    !             ! Get leaf biomass
    !             call leaf_biomass_c(site%plots(ip)%trees(it))
    !             leaf_bm = site%plots(ip)%trees(it)%leaf_bm

    !         !    thinned: if (it > thin_num) then
    !             thinned: if(site%plots(ip)%trees(it)%diam_bht > thin_thresh) then 

    !                 ! Tree is large enough to keep

    !                 ! We have already thinned the appropriate amount
    !                 ! Tree survives

    !                 lt = lt + 1 ! Increment live tree counter

    !                 ! Increment tree age
    !                 site%plots(ip)%trees(it)%tree_age=                &
    !                     site%plots(ip)%trees(it)%tree_age + 1

    !                 ! Copy tree to new location in array
    !                 call copy_tree(site%plots(ip)%trees(lt),           &
    !                     site%plots(ip)%trees(it))

    !                 if (site%plots(ip)%trees(it)%form > 0 .and.       &
    !                     site%plots(ip)%trees(it)%forska_ht < 1.83) then
    !                     ! Add up shrub live foliage and fine twigs
    !                     site%plots(ip)%soil%shrubLitter =              &
    !                         site%plots(ip)%soil%shrubLitter + leaf_bm + &
    !                         site%plots(ip)%trees(it)%branchC*         &
    !                         PERC_BRANCHES(1)

    !                 end if

    !                 ! Calculate litterfall
    !                 ! Here we convert to dry biomass because soiln
    !                 ! subroutine calculates weight loss, not C loss
    !                 if (site%species(is)%conifer) then

    !                     ! Conifers don't drop all their needles
    !                     lc = site%species(is)%litter_class
    !                     site%plots(ip)%soil%litter(lc) =               &
    !                         site%plots(ip)%soil%litter(lc) +           &
    !                         leaf_bm*(CON_LEAF_RATIO)/B_TO_C
    !                     NPP_loss = NPP_loss + leaf_bm*CON_LEAF_RATIO

    !                 else

    !                     lc = site%species(is)%litter_class
    !                     site%plots(ip)%soil%litter(lc) =               &
    !                         site%plots(ip)%soil%litter(lc) +           &
    !                         leaf_bm/B_TO_C
    !                     NPP_loss = NPP_loss + leaf_bm

    !                 end if

    !             else thinned
    !                 ! We harvest the tree

    !                 dt = dt + 1 ! Increment dead tree counter

    !                 ! Copy to list of dead trees for output
    !                 call copy_tree(site%plots(ip)%deadtrees(dt),       &
    !                     site%plots(ip)%trees(it))

    !                 ! Set cells of that tree to unfilled
    !                 trow = site%plots(ip)%trees(it)%row
    !                 tcol = site%plots(ip)%trees(it)%col
    !                 site%plots(ip)%cells(trow, tcol) = 0

    !                 ! Add biomass to array of dead biomass
    !                 site%plots(ip)%d_type(IHARV) =                     &
    !                     site%plots(ip)%d_type(IHARV) +                 &
    !                     site%plots(ip)%trees(it)%biomC + leaf_bm

    !                 ! Add roots to litter (we are removing the rest
    !                 ! of the slash off the plot)
    !                 bcr = site%plots(ip)%trees(it)%rootC
    !                 site%plots(ip)%soil%litter(IROOT) =                &
    !                     site%plots(ip)%soil%litter(IROOT) + bcr/B_TO_C

    !                 ! Get woody biomass
    !                 biomC = site%plots(ip)%trees(it)%biomC +          &
    !                     site%plots(ip)%trees(it)%rootC

    !                 ! Decrease NPP
    !                 if (site%species(is)%conifer) then
    !                     NPP_loss = NPP_loss + biomC + leaf_bm
    !                 else
    !                     NPP_loss = NPP_loss + biomC + leaf_bm
    !                 end if

    !             end if thinned

    !         end do thinloop

    !         ! Deallocate dbh arrays
    !         ! if (allocated(dbh)) deallocate(dbh)
    !         ! if (allocated(dbh_ind)) deallocate(dbh_ind)

    !     end if hntrees

    !     ! Set number of live and dead trees
    !     site%plots(ip)%numtrees = lt
    !     site%plots(ip)%num_dead = dt
    !   end do

    ! end subroutine Thin



    subroutine Shearblade(site)

      ! Data dictionary: calling argumnets
      type(SiteData),        intent(inout) :: site ! Site object

      ! Data dictionary: local variables
      integer, dimension(:), allocatable :: dbh_ind          ! Index of sorted tree DBH array
      real,    dimension(:), allocatable :: dbh              ! Tree DBH (cm)
      real                               :: totbiomC         ! Plot-wide biomass (tC)
      real                               :: NPP_loss         ! Loss of NPP from mortality
      real                               :: consRoot         ! Proportion roots consumed by fire (0-1)
      real                               :: wind_prob        ! Random number for windthrow
      real                               :: biomC            ! Total tree biomass (tC)
      real                               :: leaf_bm          ! Leaf biomass (tC)
      real                               :: bcr              ! Root biomass (tC)
      real                               :: bctw             ! Twig biomass (tC)
      real                               :: bcs              ! Stem biomass (tC)
      real                               :: bcbr             ! Total branch biomass (tC)
      real                               :: bcsb             ! Small branch biomass (tC)
      real                               :: bclb             ! Large branch biomass (tC)
      integer                            :: num_species      ! Number of species on site
      integer                            :: it, ip, iu           ! Looping indices
      integer                            :: dt               ! Counter for dead trees
      integer                            :: lt               ! Counter for live trees
      integer                            :: ind              ! Tree array index
      integer                            :: is               ! Species index
      integer                            :: trow             ! Row of tree
      integer                            :: tcol             ! Column of tree
      logical                            :: sel_cut          ! Are we selective cutting?
      integer                            :: snum             ! Tree growth stressor
      integer                            :: lc               ! Litter class
      integer                            :: i      ! Looping index

      ! We are shearblading

      ! Bulldoze the soil
      do ip = 1, site%numplots
        call bulldoze_soil(site%plots(ip)%soil, site%psolrem)

        sntrees: if (site%plots(ip)%numtrees > 0) then

            ! Initialize counters for dead and live trees
            dt = 0
            lt = 0

            ! Set fire and wind to 0
            site%plots(ip)%fire = 0
            site%plots(ip)%wind = 0

            shearloop: do it = 1, site%plots(ip)%numtrees

                ! Get species index
                is = site%plots(ip)%trees(it)%species_index

                ! Get leaf biomass
                call leaf_biomass_c(site%plots(ip)%trees(it))
                leaf_bm = site%plots(ip)%trees(it)%leaf_bm

                dt = dt + 1 ! Increment dead tree counter

                ! Copy tree to list of deadtrees
                call copy_tree(site%plots(ip)%deadtrees(dt),           &
                    site%plots(ip)%trees(it))

                ! Set cells of that tree to unfilled
                trow = site%plots(ip)%trees(it)%row
                tcol = site%plots(ip)%trees(it)%col
                site%plots(ip)%cells(trow, tcol) = 0

                ! Add biomass to array of dead biomass
                site%plots(ip)%d_type(IHARV) =                         &
                    site%plots(ip)%d_type(IHARV) +                     &
                    site%plots(ip)%trees(it)%biomC + leaf_bm

                ! Update "stressor"
                site%plots(ip)%deadtrees(dt)%stressor = IHARV

                    ! Add roots to litter (we are removing the rest
                    ! of the slash off the plot)
                    bcr = site%plots(ip)%trees(it)%rootC
                    site%plots(ip)%soil%litter(IROOT) =                &
                        site%plots(ip)%soil%litter(IROOT) + bcr/B_TO_C

                    ! Get woody biomass
                    biomC = site%plots(ip)%trees(it)%biomC +          &
                        site%plots(ip)%trees(it)%rootC

                    ! Decrease NPP
                    if (site%species(is)%conifer) then
                        NPP_loss = NPP_loss + biomC + leaf_bm
                    else
                        NPP_loss = NPP_loss + biomC + leaf_bm
                    end if

            end do shearloop

        end if sntrees

        ! Set number of live and trees
        site%plots(ip)%numtrees = lt
        site%plots(ip)%num_dead = dt
      end do

    end subroutine Shearblade

    subroutine Harvest(site)
    ! selective cut
    type(SiteData),        intent(inout) :: site ! Site object

    ! Data dictionary: local variables
    real                               :: totbiomC         ! Plot-wide biomass (tC)
    real                               :: NPP_loss         ! Loss of NPP from mortality
    real                               :: biomC            ! Total tree biomass (tC)
    real                               :: leaf_bm          ! Leaf biomass (tC)
    real                               :: bcr              ! Root biomass (tC)
    real                               :: bctw             ! Twig biomass (tC)
    real                               :: bcs              ! Stem biomass (tC)
    real                               :: bcbr             ! Total branch biomass (tC)
    real                               :: bcsb             ! Small branch biomass (tC)
    real                               :: bclb             ! Large branch biomass (tC)
    integer                            :: num_species      ! Number of species on site
    integer                            :: it, ip, iu           ! Looping indices
    integer                            :: dt               ! Counter for dead trees
    integer                            :: lt               ! Counter for live trees
    integer                            :: ind              ! Tree array index
    integer                            :: is               ! Species index
    integer                            :: trow             ! Row of tree
    integer                            :: tcol             ! Column of tree
    logical                            :: cutting_survive  ! Does tree survive selective harvest?
    integer                            :: snum             ! Tree growth stressor
    integer                            :: lc               ! Litter class
    integer                            :: i      ! Looping index

    ! turn off flag if not white spruce; this won't be replanted
    if(site%classification /= "WSPS" .and. site%classification /= "WSPP") site%mgmtflag = 0


    do ip = 1, site%numplots
      site%plots(ip)%fire = 0
      site%plots(ip)%wind = 0

      ! Initialize counters for dead and live trees
      dt = 0
      lt = 0

      ! Kill (trample?) lots of seedlings
      do is = 1, num_species
        site%plots(ip)%seedling(is) = site%plots(ip)%seedling(is)*0.01
        site%plots(ip)%seedbank(is) = site%plots(ip)%seedbank(is)*0.02
      end do
      ! Doesn't seem to do anything for succession

      cntrees: if (site%plots(ip)%numtrees > 0) then

          cutloop: do it = 1, site%plots(ip)%numtrees

              ! Get species index
              is = site%plots(ip)%trees(it)%species_index

              ! Get leaf biomass
              call leaf_biomass_c(site%plots(ip)%trees(it))
              leaf_bm = site%plots(ip)%trees(it)%leaf_bm

              ! Do we cut the tree?
              call selective_cut(site%plots(ip)%trees(it),           &
                  0.0, cutting_survive)

              cut: if (cutting_survive) then

                  ! Tree does not get harvested and does not
                  ! die from age/growth stress

                  lt = lt + 1 ! Increment live tree counter

                  ! Increment tree age
                  site%plots(ip)%trees(it)%tree_age =                &
                      site%plots(ip)%trees(it)%tree_age + 1

                  ! Copy tree to new location in array
                  call copy_tree(site%plots(ip)%trees(lt),           &
                      site%plots(ip)%trees(it))

                  if (site%plots(ip)%trees(it)%form > 0 .and.        &
                      site%plots(ip)%trees(it)%forska_ht < 1.83) then
                      ! Add up shrub live foliage and fine twigs
                      site%plots(ip)%soil%shrubLitter =              &
                          site%plots(ip)%soil%shrubLitter + leaf_bm + &
                          site%plots(ip)%trees(it)%branchC*          &
                          PERC_BRANCHES(1)

                  end if

                  ! Calculate litterfall
                  ! Here we convert to dry biomass because soiln
                  ! subroutine calculates weight loss, not C loss
                  if (site%species(is)%conifer) then

                      ! Conifers don't drop all their needles
                      lc = site%species(is)%litter_class
                      site%plots(ip)%soil%litter(lc) =               &
                          site%plots(ip)%soil%litter(lc) +           &
                          leaf_bm*(CON_LEAF_RATIO)/B_TO_C
                      NPP_loss = NPP_loss + leaf_bm*CON_LEAF_RATIO

                  else

                      lc = site%species(is)%litter_class
                      site%plots(ip)%soil%litter(lc) =               &
                          site%plots(ip)%soil%litter(lc) +           &
                          leaf_bm/B_TO_C
                      NPP_loss = NPP_loss + leaf_bm

                  end if

              else cut

                  ! We potentially harvest the tree

                  dt = dt + 1 ! Increment dead tree counter

                  ! Copy to list of dead trees for output
                  call copy_tree(site%plots(ip)%deadtrees(dt),       &
                      site%plots(ip)%trees(it))
                  site%plots(ip)%deadtrees(dt)%stressor = IHARV

                  ! Set cells of that tree to unfilled
                  trow = site%plots(ip)%trees(it)%row
                  tcol = site%plots(ip)%trees(it)%col
                  site%plots(ip)%cells(trow, tcol) = 0

                  ! Get most limiting growth factor
                  snum = site%plots(ip)%trees(it)%stressor

                  cutcheck: if (.not. cutting_survive) then

                      ! Add biomass to array of dead biomass
                      site%plots(ip)%d_type(IHARV) =                 &
                          site%plots(ip)%d_type(IHARV) +             &
                          site%plots(ip)%trees(it)%biomC + leaf_bm

                      ! Add roots to litter
                      bcr = site%plots(ip)%trees(it)%rootC
                      site%plots(ip)%soil%litter(IROOT) =            &
                          site%plots(ip)%soil%litter(IROOT) +        &
                          bcr/B_TO_C

                      if (site%plots(ip)%trees(it)%diam_bht <       &
                          site%dbh_thresh) then

                          ! These remain as slash

                          ! Get woody biomass
                          biomC = site%plots(ip)%trees(it)%biomC +   &
                              site%plots(ip)%trees(it)%rootC

                          ! Branches
                          bcbr = site%plots(ip)%trees(it)%branchC
                          ! Convert branch litter into twigs, small
                          ! branches, and large branches
                          bctw = bcbr*PERC_BRANCHES(1)
                          bcsb = bcbr*PERC_BRANCHES(2)
                          bclb = bcbr*PERC_BRANCHES(3)

                          ! Twigs
                          site%plots(ip)%soil%litter(ITW) =          &
                              site%plots(ip)%soil%litter(ITW) +      &
                              bctw/B_TO_C

                          ! Small branches
                          site%plots(ip)%soil%litter(ISBR) =         &
                              site%plots(ip)%soil%litter(ISBR) +     &
                              bcsb/B_TO_C

                          ! Large branches
                          site%plots(ip)%soil%litter(ILBR) =         &
                              site%plots(ip)%soil%litter(ILBR) +     &
                              bclb/B_TO_C

                          ! Stems
                          bcs = site%plots(ip)%trees(it)%stemC

                          ! Small boles (DBH < 10) vs. large boles
                          if (site%plots(ip)%trees(it)%diam_bht      &
                              > 10.0) then
                              site%plots(ip)%soil%litter(ILBL) =     &
                                  site%plots(ip)%soil%litter(ILBL) + &
                                  bcs/B_TO_C
                          else
                              site%plots(ip)%soil%litter(ISBL) =     &
                                  site%plots(ip)%soil%litter(ISBL) + &
                                  bcs/B_TO_C
                          end if

                          ! Leaves
                          lc = site%species(is)%litter_class
                          site%plots(ip)%soil%litter(lc) =           &
                                site%plots(ip)%soil%litter(lc) +       &
                              leaf_bm/B_TO_C

                          ! Decrease NPP
                          if (site%species(is)%conifer) then
                              NPP_loss = NPP_loss + biomC + leaf_bm
                          else
                              NPP_loss = NPP_loss + biomC + leaf_bm
                          end if

                      end if

                  end if cutcheck

              end if cut

          end do cutloop

      end if cntrees

      ! Set number of live and dead trees
      site%plots(ip)%numtrees = lt
      site%plots(ip)%num_dead = dt

    end do
    end subroutine Harvest


    subroutine Prune(site)
      !Prune low branches and remove small trees
      type(SiteData), intent(inout) :: site

      !local variables
      integer :: ip, it, is ! indices
      integer :: dt, lt !live and dead trees
      integer :: trow, tcol
      real :: leaf_bm, NPP_loss, biomC
      real :: bcbr, d_bcbr, bcr, d_bcr


      do ip = 1, site%numplots
        pntrees: if (site%plots(ip)%numtrees > 0) then
          lt = 0
          dt = 0

          site%plots(ip)%fire = 0
          site%plots(ip)%wind = 0

          do it = 1, site%plots(ip)%numtrees
            pruneloop: if(site%plots(ip)%trees(it)%forska_ht <= 2.0) then
              ! kill tree and remove AGB
              dt = dt + 1
              ! Get species index
              is = site%plots(ip)%trees(it)%species_index
              ! Get leaf biomass
              call leaf_biomass_c(site%plots(ip)%trees(it))
              leaf_bm = site%plots(ip)%trees(it)%leaf_bm
              ! Copy tree to list of deadtrees
              call copy_tree(site%plots(ip)%deadtrees(dt),           &
                  site%plots(ip)%trees(it))
                  ! Set cells of that tree to unfilled
              trow = site%plots(ip)%trees(it)%row
              tcol = site%plots(ip)%trees(it)%col
              site%plots(ip)%cells(trow, tcol) = 0
              ! Add biomass to array of dead biomass
              site%plots(ip)%d_type(IHARV) =                         &
                  site%plots(ip)%d_type(IHARV) +                     &
                  site%plots(ip)%trees(it)%biomC + leaf_bm

                ! Update "stressor"
                site%plots(ip)%deadtrees(dt)%stressor = IHARV
                ! Add roots to litter (we are removing the rest
                ! of the slash off the plot)
                bcr = site%plots(ip)%trees(it)%rootC
                site%plots(ip)%soil%litter(IROOT) =                &
                site%plots(ip)%soil%litter(IROOT) + bcr/B_TO_C

                ! Get woody biomass
                biomC = site%plots(ip)%trees(it)%biomC +          &
                    site%plots(ip)%trees(it)%rootC

                    NPP_loss = NPP_loss + biomC + leaf_bm

              else pruneloop
                ! remove branches below 2m
                if(site%plots(ip)%trees(it)%canopy_ht<= 2.0) then
                  site%plots(ip)%trees(it)%canopy_ht = 2.0
                  call stem_shape(site%plots(ip)%trees(it))
                  leaf_bm = site%plots(ip)%trees(it)%leaf_bm
                  call leaf_biomass_c(site%plots(ip)%trees(it))

                  bcbr = site%plots(ip)%trees(it)%branchC

                  ! Update biomass C and N given new clear branch
                  ! bole height
                  call biomass_c(site%plots(ip)%trees(it))
                  call biomass_n(site%plots(ip)%trees(it))

                  ! How much wood litter did we lose?
                  d_bcbr = bcbr - site%plots(ip)%trees(it)%branchC

                  ! Here we convert to dry biomass because soiln
                  ! subroutine calculates weight loss not C loss

                  NPP_loss = d_bcr/B_TO_C + (leaf_bm - site%plots(ip)%trees(it)%leaf_bm)
                end if

            end if pruneloop

              site%plots(ip)%NPP = site%plots(ip)%NPP -                          &
                  NPP_loss/plotsize/M2_TO_HEC

          end do
        end if pntrees

        ! Set number of live and dead trees
        site%plots(ip)%numtrees = lt
        site%plots(ip)%num_dead = dt

      end do


    end subroutine Prune

end module Manage
