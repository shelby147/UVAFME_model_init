module Landscape

!*******************************************************************************
  !
  ! This module deals with site-site calculations
  !
!*******************************************************************************

  use Parameters
  use Constants
  use Site
  use Species
  use data_dicts
  use dictionary
  use Input
  use lists

  implicit none


 contains

    !:.........................................................................:

    subroutine Seed_Rain(all_sites, sndx, dict, species_ids)
        !
        !  Computes seed rain to adjacent sites from a site
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/10/20    A. C. Foster         Original Code
        !

        ! Data dictionary: calling arguments
        type(SiteData), dimension(:),   intent(inout) :: all_sites   ! Array of site objects
        type(dict_struct),              pointer       :: dict        ! Dictionary for site ids and location
        integer,                        intent(in)    :: sndx        ! Index of current site
        integer,        dimension(:,:), intent(in)    :: species_ids ! Species index locations for each site

        ! Data dictionary: local variables
        type(dict_dat)        :: dat          ! Array position of site
        integer, dimension(1) :: locs         ! Location of species in array
        real                  :: r_tree       ! Factor for determining seed rain
        real                  :: a, b         ! Species-specific seed rain parameters
        real                  :: maxdist_m    ! Maximum distance for seed rain (m)
        real                  :: pdist        ! Seed dispersal probability
        real                  :: dist         ! Distance to adjacent site (m)
        real                  :: pixel_length ! Diagonal length of pixels (m)
        integer               :: r, c         ! Row and column of site
        integer               :: maxdist_p    ! Maximum distance for seed rain (pixels)
        integer               :: num_species  ! Number of species at site
        integer               :: rmin, cmin   ! Minimum row/columns to check
        integer               :: rmax, cmax   ! Maximum row/columns to check
        integer               :: rr, cc, is   ! Looping indices
        integer               :: ind          ! Index position mapped from row/column
        integer               :: isite        ! Array position of site
        integer               :: site_is      ! species_id for site sndx
        integer               :: isite_is     ! species_id for site isite
        character(len=80)     :: indc         ! Index position mapped from row/column

        ! Calculate pixel length (diagonal)
        pixel_length = sqrt(float(pixel_size**2) + float(pixel_size**2))

        ! Get row/column value of site
        r = all_sites(sndx)%row
        c = all_sites(sndx)%col

        ! Get number of species at site
        num_species = size(all_sites(sndx)%species)

        specloop: do is = 1, num_species

            ! Get seed rain factor
            r_tree = min(1.0, all_sites(sndx)%recr_trees(is)/5.0)

            ! If enough reproductive trees of this species
            minprob: if (r_tree >= min_seedRain_prob) then

                a = all_sites(sndx)%species(is)%seedRain_a
                b = all_sites(sndx)%species(is)%seedRain_b

                ! Get maximum distance to check for seed rain
                maxdist_m = -1.0*(log(min_seedRain_prob)/b)
                maxdist_p = 1 + int(maxdist_m/pixel_size)

                ! Get minimum row/column values
                rmin = max(1, (r - maxdist_p))
                rmax = min(nrows, (r + maxdist_p))
                cmin = max(1, (c - maxdist_p))
                cmax = min(ncols, (c + maxdist_p))

                ! Loop through landscape to compute distance & probabilities
                rloop: do rr = rmin, rmax
                    cloop: do cc = cmin, cmax

                        ! Get index position in landscape-wide grid
                        ind = ((rr-1)*nrows) + cc

                        ! Write to character and check if it is in the
                        ! dictionary
                        write(indc, *) ind
                        haskey: if (dict_has_key(dict, indc)) then

                            ! Calculate distance to row, col (meters)
                            dist = max((float(pixel_size)*                     &
                                sqrt((float(abs(rr -r))*float(abs(rr-r))) +    &
                                (float(abs(cc-c))*float(abs(cc-c))))) -        &
                                pixel_length, 0.0)

                            ! Calculate seed dispersal probability
                            pdist = (min(1.0, exp(a - b*dist)/(exp(a))))*r_tree

                            ! Get array position of adjacent site
                            dat = dict_get_key(dict, indc)
                            read(dat%string, *) isite

                            ! is is the species_id for site sndx
                            site_is = species_ids(sndx, is)
                            locs = findloc(species_ids(isite,:), value = site_is)
                            isite_is = locs(1)

                            ! Add that probability to loop pixel
                            all_sites(isite)%species(isite_is)%fc_recr =       &
                                all_sites(isite)%species(isite_is)%fc_recr +   &
                                pdist
                            all_sites(isite)%n_adj = all_sites(isite)%n_adj + 1

                        end if haskey

                    end do cloop
                end do rloop

            end if minprob

        end do specloop

    end subroutine Seed_Rain

    !:.........................................................................:

    subroutine Management(all_sites, numsites, year, pactvmgmt)
        !
        !  Computes seed rain to adjacent sites from a site
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    08/31/22    S. Sundquist         Original Code
        !
  
        ! Data dictionary
        type(SiteData), dimension(:),   intent(inout) :: all_sites   ! Array of site objects
        integer, intent(in) :: numsites ! number of sites
        integer, intent(in) :: year
        real, dimension(14), intent(inout) :: pactvmgmt ! active management likelihood per site (sites treated per year)
  
        ! Local variables
        character (len = 3) :: unit ! which management subunit is this site located in (for active management)
        real, dimension(14) :: actvmgmt ! active management input data (ha treated per year)
        integer :: sndx
        integer :: sitearea
        integer, dimension(:), allocatable :: REP, WSPP, WSB, BSP, ASP, BIR, ABI, MIX, OTH, WSPS ! site indexes for each forest class
        real, dimension(10) :: clengths
        real :: rand
        integer :: irand
        integer :: salv_sndx !site ID for salvage-eligible site 
        integer :: lt, dt ! Live and dead tree counters for salvage 
        integer :: i !looping indices
  
        active_mgmt: if(active_management) then
  
        ! Determine unit
        unit = all_sites(1)%unit
        do sndx = 1, numsites
          if(unit /= all_sites(sndx)%unit) then
            call warning("Warning: not all sites are in same unit")
            write(logf, *) "Applying management for unit", unit
            exit
          end if
        end do
  
        ! Determine classifications for each site, add to list
        clengths = 0.0
        do sndx = 1, numsites
          if(all_sites(sndx)%classification == "REP") then
            call appendix(REP, sndx)
            clengths(1) = clengths(1) + 1.0
          else if(all_sites(sndx)%classification == "WSPP") then
            call appendix(WSPP, sndx)
            clengths(2) = clengths(2) + 1.0
          else if(all_sites(sndx)%classification == "WSB") then
            call appendix(WSB, sndx)
            clengths(3) = clengths(3) + 1.0
          else if(all_sites(sndx)%classification == "BSP") then
            call appendix(BSP, sndx)
            clengths(4) = clengths(4) + 1.0
          else if(all_sites(sndx)%classification == "ASP") then
            call appendix(ASP, sndx)
            clengths(5) = clengths(5) + 1.0
          else if(all_sites(sndx)%classification == "BIR") then
            call appendix(BIR, sndx)
            clengths(6) = clengths(6) + 1.0
          else if(all_sites(sndx)%classification == "ABI") then
            call appendix(ABI, sndx)
            clengths(7) = clengths(7) + 1.0
          else if(all_sites(sndx)%classification == "MIX") then
            call appendix(MIX, sndx)
            clengths(8) = clengths(8) + 1.0
          else if(all_sites(sndx)%classification == "OTH") then
            call appendix(OTH, sndx)
            clengths(9) = clengths(9) + 1.0
          else if(all_sites(sndx)%classification == "WSPS") then
            call appendix(WSPS, sndx)
            clengths(10) = clengths(10) + 1.0
          end if
        end do
  
        call read_active_management(unit, actvmgmt)
  
        ! Scale management to scaler from runtime (default 1.0)
        actvmgmt = actvmgmt*management_scale
  
        sitearea = pixel_size**2*M2_TO_HEC ! area in ha per site
        pactvmgmt = actvmgmt/sitearea
  
        ! make lists(?) of sites with each classification type
  
          !mgmt columns 1-5: bau_wsp, bau_mixed, bau_birch, bau_thin, bau_salvage
          !columns 6-10: prj_wsp, prj_mixed, prj_birch, prj_thin, prj_salvage
          !columns 11-14: prj_bsp, prj_thinbsp, prj_prunebsp, prj_shearbsp
  
          scenario: if(management_scenario == "BAU") then
            ! spruce harvest?
            rand = urand()
            do while (rand <= pactvmgmt(1))
              print *, "Harvesting spruce"
              ! flag a wsp site for harvest
              if(clengths(10) == 0) then
                ! write(logf, *) "No forest available for selected treatment: year", year, "WSP sawtimber harvest"
                pactvmgmt(1) = 1.0 + (int(pactvmgmt(1)))
                exit
              end if
              irand = int(urand(1.0, clengths(10))) ! select index to treat from WSP classiied list
              sndx = WSPS(irand) ! identify site index to treat
              call deleteix(WSPS, irand)
              clengths(10) = clengths(10) - 1
              all_sites(sndx)%mgmtflag = 1
              ! reduce probability of another site being selected (in case p > 1)
              pactvmgmt(1) = pactvmgmt(1) - 1.0
              rand = urand()
            end do
  
            ! mixed harvest?
            rand = urand()
            do while (rand <= pactvmgmt(2))
              print *, "Harvesting mixed forest"
              if(clengths(8) == 0) then
                pactvmgmt(2) = 1.0 + (int(pactvmgmt(2)))
                exit
              end if
              irand = int(urand(1.0, clengths(8)))
              sndx = MIX(irand)
              call deleteix(MIX, irand)
              clengths(8) = clengths(8) - 1
              all_sites(sndx)%mgmtflag = 1
              pactvmgmt(2) = pactvmgmt(2) - 1.0
              rand = urand()
            end do
  
            ! birch harvest?
            rand = urand()
            do while (rand <= pactvmgmt(3))
              print *, "Harvesting birch"
              if(clengths(6) == 0) then
                pactvmgmt(3) = 1.0 + (int(pactvmgmt(3)))
                exit
              end if
              irand = int(urand(1.0, clengths(6)))
              sndx = BIR(irand)
              call deleteix(BIR, irand)
              clengths(6) = clengths(6) - 1
              all_sites(sndx)%mgmtflag = 1
              pactvmgmt(3) = pactvmgmt(3) - 1.0
              rand = urand()
            end do
            rand = urand()
  
            ! thinning?
            rand = urand()
            do while (rand <= pactvmgmt(4))
              print *, "Thinning pole spruce"
              if(clengths(2) == 0) then
                pactvmgmt(4) = 1.0 + (int(pactvmgmt(4)))
                exit
              end if
              irand = int(urand(1.0, clengths(2))) ! select index to treat from WSP classiied list
              sndx = WSPP(irand) ! identify site index to treat
              call deleteix(WSPP, irand)
              clengths(2) = clengths(2) - 1
              all_sites(sndx)%mgmtflag = 2
              ! reduce probability of another site being selected (in case p > 1)
              pactvmgmt(4) = pactvmgmt(4) - 1.0
              rand = urand()
            end do
  
            ! ! salvage?
            rand = urand()
            do while(rand <= pactvmgmt(5))
              print *, "Salvaging"
              ! flag a salvage-eligible site for salvage 
              salv_sndx = -1 !haven't found one
  
              ! First, look for WSPS site with mortality event 
              do i = 1, clengths(10)
                if(salv_sndx > 0 ) exit 
                sndx = WSPS(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
  
              ! Next, WSPP
              do i = 1, clengths(2)
                if(salv_sndx > 0 ) exit 
                sndx = WSPP(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, WSB
              do i = 1, clengths(3)
                if(salv_sndx > 0 ) exit 
                sndx = WSB(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, BIR
              do i = 1, clengths(6)
                if(salv_sndx > 0 ) exit 
                sndx = BIR(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, MIX
              do i = 1, clengths(8)
                if(salv_sndx > 0 ) exit 
                sndx = MIX(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, ABI
              do i = 1, clengths(7)
                if(salv_sndx > 0 ) exit 
                sndx = ABI(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, ASP
              do i = 1, clengths(5)
                if(salv_sndx > 0 ) exit 
                sndx = ASP(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              if(salv_sndx > 0) then 
                all_sites(salv_sndx)%mgmtflag = 3
                pactvmgmt(5) = pactvmgmt(5) - 1.0
                rand = urand()
              else if (salv_sndx == -1) then ! can't find site to salvage
                pactvmgmt(5) = (int(pactvmgmt(5))) + 1.0
                exit 
              end if 
              
            end do
  
  
  
          else if(management_scenario == "PRJ") then scenario
            ! spruce harvest?
            rand = urand()
            do while (rand <= pactvmgmt(6))
              print *, "Harvesting spruce"
              ! flag a wsp site for harvest
              if(clengths(10) == 0) then
                pactvmgmt(6) = 1.0 + (int(pactvmgmt(6)))
                exit
              end if
              irand = int(urand(1.0, clengths(10))) 
              sndx = WSPS(irand)
              call deleteix(WSPS, irand)
              clengths(10) = clengths(10) - 1
              all_sites(sndx)%mgmtflag = 1 
              ! reduce probability of another site being selected (in case p > 1)
              pactvmgmt(6) = pactvmgmt(6) - 1.0
              rand = urand()
            end do
  
            ! mixed harvest?
            rand = urand()
            do while (rand <= pactvmgmt(7))
              print *, "Harvesting mixed forest"
              if(clengths(8) == 0) then
                pactvmgmt(7) = 1.0 + (int(pactvmgmt(7)))
                exit
              end if
              irand = int(urand(1.0, clengths(8)))
              sndx = MIX(irand)
              call deleteix(MIX, irand)
              clengths(8) = clengths(8) - 1
              all_sites(sndx)%mgmtflag = 1
              pactvmgmt(7) = pactvmgmt(7) - 1.0
              rand = urand()
            end do
  
            ! birch harvest?
            rand = urand()
            do while (rand <= pactvmgmt(8))
              print *, "Harvesting birch"
              if(clengths(6) == 0) then
                pactvmgmt(8) = 1.0 + (int(pactvmgmt(8)))
                exit
              end if
              irand = int(urand(1.0, clengths(6)))
              sndx = BIR(irand)
              call deleteix(BIR, irand)
              clengths(6) = clengths(6) - 1
              all_sites(sndx)%mgmtflag = 1
              pactvmgmt(8) = pactvmgmt(8) - 1.0
              rand = urand()
            end do
            rand = urand()
  
            ! thinning?
            rand = urand()
            do while (rand <= pactvmgmt(9))
              print *, "Thinning pole spruce"
              if(clengths(2) == 0) then
                pactvmgmt(9) = 1.0 + (int(pactvmgmt(9)))
                exit
              end if
              irand = int(urand(1.0, clengths(2))) 
              sndx = WSPP(irand)
              call deleteix(WSPP, irand)
              clengths(2) = clengths(2) - 1
              all_sites(sndx)%mgmtflag = 2 
              ! reduce probability of another site being selected (in case p > 1)
              pactvmgmt(9) = pactvmgmt(9) - 1.0
              rand = urand()
            end do
  
            ! salvage? 
            rand = urand()
            do while(rand <= pactvmgmt(10))
              print *, "Salvaging"
              ! flag a salvage-eligible site for salvage 
              salv_sndx = -1 !haven't found one
  
              ! First, look for WSPS site with mortality event 
              do i = 1, clengths(10)
                if(salv_sndx > 0 ) exit 
                sndx = WSPS(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
  
              ! Next, WSPP
              do i = 1, clengths(2)
                if(salv_sndx > 0 ) exit 
                sndx = WSPP(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, WSB
              do i = 1, clengths(3)
                if(salv_sndx > 0 ) exit 
                sndx = WSB(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, BIR
              do i = 1, clengths(6)
                if(salv_sndx > 0 ) exit 
                sndx = BIR(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, MIX
              do i = 1, clengths(8)
                if(salv_sndx > 0 ) exit 
                sndx = MIX(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, ABI
              do i = 1, clengths(7)
                if(salv_sndx > 0 ) exit 
                sndx = ABI(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
              
              ! Next, ASP
              do i = 1, clengths(5)
                if(salv_sndx > 0 ) exit 
                sndx = ASP(i)
                call countldtrees(all_sites(sndx), lt, dt)
                if (dt/(lt+dt) > 0.33) salv_sndx = sndx ! more than 1/3 of trees on plot dead
              end do 
  
              if(salv_sndx > 0) then 
                all_sites(salv_sndx)%mgmtflag = 3
                pactvmgmt(5) = pactvmgmt(5) - 1.0
                rand = urand()
              else if (salv_sndx == -1) then ! can't find site to salvage
                pactvmgmt(5) = (int(pactvmgmt(5))) + 1.0
                exit 
              end if 
              
            end do ! End of salvage
  
            ! black spruce harvest?
            rand = urand()
            do while (rand <= pactvmgmt(11))
              print *, "Harvesting black spruce"
              if(clengths(4) == 0) then
                pactvmgmt(11) = 1.0 + (int(pactvmgmt(11)))
                exit
              end if
              irand = int(urand(1.0, clengths(4))) 
              sndx = BSP(irand) 
              call deleteix(BSP, irand)
              clengths(4) = clengths(4) - 1
              all_sites(sndx)%mgmtflag = 1
              ! reduce probability of another site being selected (in case p > 1)
              pactvmgmt(11) = pactvmgmt(11) - 1.0
              rand = urand()
            end do
  
            ! black spruce thinning?
            rand = urand()
            do while (rand <= pactvmgmt(12))
              print *, "Thinning black spruce"
              if(clengths(4) == 0) then
                pactvmgmt(12) = 1.0 + (int(pactvmgmt(12)))
                exit
              end if
              irand = int(urand(1.0, clengths(4))) 
              sndx = BSP(irand) 
              call deleteix(BSP, irand)
              clengths(4) = clengths(4) - 1
              all_sites(sndx)%mgmtflag = 2
              ! reduce probability of another site being selected (in case p > 1)
              pactvmgmt(12) = pactvmgmt(12) - 1.0
              rand = urand()
            end do
  
            ! black spruce prunning?
            rand = urand()
            do while (rand <= pactvmgmt(13))
              print *, "Pruning black spruce"
              if(clengths(4) == 0) then
                pactvmgmt(13) = 1.0 + (int(pactvmgmt(13)))
                exit
              end if
              irand = int(urand(1.0, clengths(4))) 
              sndx = BSP(irand) 
              call deleteix(BSP, irand)
              clengths(4) = clengths(4) - 1
              all_sites(sndx)%mgmtflag = 5 
              ! reduce probability of another site being selected (in case p > 1)
              pactvmgmt(13) = pactvmgmt(13) - 1.0
              rand = urand()
            end do
  
            ! black spruce shearblading?
            rand = urand()
            do while (rand <= pactvmgmt(14))
              print *, "Shearblading black spruce"
              if(clengths(4) == 0) then
                pactvmgmt(14) = 1.0 + (int(pactvmgmt(14)))
                exit
              end if
              irand = int(urand(1.0, clengths(4))) 
              sndx = BSP(irand)
              call deleteix(BSP, irand)
              clengths(4) = clengths(4) - 1
              all_sites(sndx)%mgmtflag = 4
              ! reduce probability of another site being selected (in case p > 1)
              pactvmgmt(14) = pactvmgmt(14) - 1.0
              rand = urand()
            end do 
  
          end if scenario
  
          do i = 1, 14
            if (pactvmgmt(i) < 1.0) pactvmgmt(i) = 0.0 
            ! reset to 0 unless we're waiting on a site(s) to be eligible for management 
          end do
  
        end if active_mgmt
  
        if(allocated(REP)) deallocate(REP)
        if(allocated(WSPP)) deallocate(WSPP)
        if(allocated(WSB)) deallocate(WSB)
        if(allocated(BSP)) deallocate(BSP)
        if(allocated(ASP)) deallocate(ASP)
        if(allocated(BIR)) deallocate(BIR)
        if(allocated(ABI)) deallocate(ABI)
        if(allocated(MIX)) deallocate(MIX)
        if(allocated(OTH)) deallocate(OTH)
        if(allocated(WSPS)) deallocate(WSPS)
  
      end subroutine Management
  
      !:.........................................................................:
  end module Landscape
