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

end module Landscape
