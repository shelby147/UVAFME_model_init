module Tree
!*******************************************************************************
  !
  ! The Tree module contains attributes and procedures relevant to tree growth
  ! and mortality
  !
!*******************************************************************************

    use Constants
    use Species
    use Random

    implicit none

    ! Define tree data type
    type                            :: TreeData
        class(SpeciesData), pointer :: spec_ptr => null() ! Pointer to species data
        real, dimension(FC_NUM)     :: env_resp           ! Species growth response to stressors (0-1)
        real                        :: diam_opt           ! Optimum diameter increment growth (cm)
        real                        :: diam_bht           ! Diameter at breast height (cm)
        real                        :: diam_canht         ! Diameter at clear branch bole height (cm)
        real                        :: canopy_ht          ! Clear branch bole height (m)
        real                        :: forska_ht          ! Tree height (m)
        real                        :: leaf_bm            ! Leaf biomass (tC)
        real                        :: biomC              ! Aboveground woody biomass (tC)
        real                        :: stemC              ! Stem biomass (tC)
        real                        :: branchC            ! Branch biomass (tC)
        real                        :: rootC              ! Root biomass (tC)
        real                        :: biomN              ! N content (tC)
        real                        :: CK                 ! Proportion of crown scorched by fire (0-1)
        real                        :: pm_tau             ! Effect of cambial burning (0-1)
        real                        :: pm_CK              ! Effect of scorch from fires (0-1)
        real                        :: row                ! Row location
        real                        :: col                ! Column location
        integer                     :: mort_count         ! How many years tree has experiences low growth
        integer                     :: tree_age           ! Tree age (years)
        integer                     :: stressor           ! Which factor is most stressful
        logical                     :: conifer            ! 1: conifer; 0: deciduous
        integer                     :: form               ! Form: 1: tree; 2: tree-like shrub;
                                                          ! 3: erect shrub; 4: prostrate shrub
        integer                     :: species_index      ! Species index - points to SpeciesData array location
        logical                     :: mort_marker        ! Is tree marked for death?
        integer                     :: id                 ! Unique id to track tree growth
        integer                     :: planted            ! Marker for whether tree was planted (management) 0/1; 2 if tree was initialized
    end type TreeData

contains


    !:.........................................................................:

    subroutine initialize_tree(self, tree_species, is, tree_id, planted)
        !
        !  Initializes a tree object with initial conditions
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        class(TreeData),           intent(inout) :: self         ! Tree object
        type(SpeciesData), target, intent(in)    :: tree_species ! Species object
        integer,                   intent(in)    :: is           ! Integer for species array
        integer,                   intent(in)    :: tree_id
        integer,                   intent(in)    :: planted     ! 0/1/2 (1 is planted, 2 is stand initialized, 0 tree is from annual renewal routine)

        ! Constructor
        self%env_resp = RNVALID
        self%diam_opt = 0.0
        self%diam_bht = 0.0
        self%diam_canht = 0.0
        self%canopy_ht = STD_HT
        self%forska_ht = STD_HT
        self%leaf_bm = 0.0
        self%biomC = 0.0
        self%stemC = 0.0
        self%branchC = 0.0
        self%rootC = 0.0
        self%biomN = 0.0
        self%CK = 0.0
        self%pm_tau = 0.0
        self%pm_CK = 0.0
        self%row = 0
        self%col = 0
        self%mort_count = 0
        self%tree_age = 0
        self%stressor = 0
        self%mort_marker = .false.
        self%species_index = is

        ! Attach species data and species index
        if (associated(self%spec_ptr)) nullify(self%spec_ptr)
        self%spec_ptr => tree_species

        ! Grab these from species data for easier use
        self%form = self%spec_ptr%form
        self%conifer = self%spec_ptr%conifer

        self%id = tree_id
        self%planted = planted

    end subroutine initialize_tree

    !:.........................................................................:

    subroutine copy_tree(self, tree)
        !
        !  Copies attributes from one tree object to another
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: self
        class(TreeData), intent(in)    :: tree
        ! Copy attributes
        self%env_resp       = tree%env_resp
        !if (associated(self%spec_ptr)) nullify(self%spec_ptr)
        self%spec_ptr       => tree%spec_ptr
        self%species_index  = tree%species_index
        self%diam_opt       = tree%diam_opt
        self%diam_bht       = tree%diam_bht
        self%diam_canht     = tree%diam_canht
        self%canopy_ht      = tree%canopy_ht
        self%forska_ht      = tree%forska_ht
        self%leaf_bm        = tree%leaf_bm
        self%biomC          = tree%biomC
        self%stemC          = tree%stemC
        self%branchC        = tree%branchC
        self%rootC          = tree%rootC
        self%biomN          = tree%biomN
        self%CK             = tree%CK
        self%pm_tau         = tree%pm_tau
        self%pm_CK          = tree%pm_CK
        self%row            = tree%row
        self%col            = tree%col
        self%mort_count     = tree%mort_count
        self%tree_age       = tree%tree_age
        self%stressor       = tree%stressor
        self%mort_marker    = tree%mort_marker
        self%form           = tree%form
        self%conifer        = tree%conifer
        self%id             = tree%id
        self%planted        = tree%planted

    end subroutine copy_tree

    !:.........................................................................:

    subroutine update_tree(self, tree_species)
        !
        !  A helper function required to have more flexible object structure
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/12     K. Holcomb           Original Code
        !

        ! Data dictionary: calling arguments
        class(TreeData),   intent(inout) :: self         ! Tree object
        type(SpeciesData), intent(in)    :: tree_species ! Species object


    end subroutine update_tree

    !:.........................................................................:

    subroutine env_stress(self, shade, fc_gdd, fc_drought, fc_perm, fc_flood, &
        envstress, fc_nutr)
        !
        !  Calculates the environmental stress on a tree using Liebig's Law of
        !   the Minimum, and sets what the most stressful factor is
        !
        !    stressor:  which factor was most stressful
        !                1: temperature
        !                2: drought
        !                3: shade
        !                4: permafrost
        !               5: flooding
        !                6: nutrients
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    01/01/16     A. C. Foster        Original Code
        !    01/27/21     A. C. Foster        Updated to read in effects rather
        !                                       than use attached species
        !                                       attributes

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout)        :: self       ! Tree object
        real,            intent(in)           :: shade      ! Growth response to shading (0-1)
        real,            intent(in)           :: fc_gdd     ! Growth response to growing degree-days (0-1)
        real,            intent(in)           :: fc_drought ! Growth response to drought (0-1)
        real,            intent(in)           :: fc_perm    ! Growth response to permafrost (0-1)
        real,            intent(in)           :: fc_flood   ! Growth response to flooding (0-1)
        real,            intent(out)          :: envstress  ! Growth response to stress
        real,            intent(in), optional :: fc_nutr

        ! Data dictionary: local variables
        real    :: minstress  ! Minimum stress from gdd, drought, shade, and (optionally) nutrients
        integer :: stressor   ! Which factor had the highest impact

        ! Get the minimum value
        if (present(fc_nutr)) then

            ! Calculate using nutrient stress
            minstress = min(fc_gdd, fc_drought, shade, fc_nutr)

            ! Find the most limiting factor
            if (minstress == fc_gdd) then
                stressor = 1
            else if (minstress == fc_drought) then
                stressor = 2
            else if (minstress == shade) then
                stressor = 3
            else if (minstress == fc_nutr) then
                stressor = 6
            end if

            if (minstress > fc_perm) then
                stressor = 4
            else if (minstress > fc_flood .and. fc_perm > fc_flood) then
                stressor = 5
            endif

            ! Output is minimum of shade, gdd, nutrient, and drought times flood
            ! and permafrost
            envstress = minstress*fc_flood*fc_perm
            !envstress = fc_gdd*fc_drought*fc_shade*fc_nutr*fc_flood*fc_perm

        else
            ! Don't use nutrient stress
            ! This is so we can calculate the potential growth and then
            ! decrease it based on available N
            minstress = min(fc_gdd, fc_drought, shade)

            ! Find the most limiting factor
            if (minstress == fc_gdd) then
                stressor = 1
            else if (minstress == fc_drought) then
                stressor = 2
            else if (minstress == shade) then
                stressor = 3
            end if

            if (minstress > fc_perm) then
                stressor = 4
            else if (minstress > fc_flood .and. fc_perm > fc_flood) then
                stressor = 5
            endif

            ! Output is minimum of shade, gdd, and drought times flood and
            ! permafrost
            envstress = minstress*fc_flood*fc_perm
            !envstress = fc_gdd*fc_drought*fc_shade*fc_flood*fc_perm

        end if

        ! Set the environmental stressors
        self%stressor = stressor
        self%env_resp(1) = fc_gdd
        self%env_resp(2) = fc_drought
        self%env_resp(3) = shade
        self%env_resp(4) = fc_perm
        self%env_resp(5) = fc_flood
        if(present(fc_nutr)) self%env_resp(6) = fc_nutr

    end subroutine env_stress

    !:.........................................................................:

    subroutine stem_shape(tree)
        !
        !  Calculates diameter at clear branch bole height (cm)
        !   Notes from original code:
        !       basal diameter = dhshape(0.0, h, dbh)
        !       diameter @ ccb = dhshape(hc, h, sdbh)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    10/10/20     A. C. Foster        Update for shrub forms
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        ! Data dictionary: local variables
        real :: hcbb ! Clear branch bole height (m)
        real :: ht   ! Tree height (m)
        real :: dbh  ! Diameter at breast height (cm)
        real :: dcbb ! Diameter at clear branch bole height (cm)
        real :: beta ! Stem shape parameter

        if (tree%form <= 2) then
            ! Tree or tree-like shrub form

            hcbb = tree%canopy_ht
            ht = tree%forska_ht
            dbh = tree%diam_bht
            beta = tree%spec_ptr%beta

            if (ht <= hcbb .or. ht <= STD_HT) then
                ! Just set to dbh
                dcbb = dbh
            else
                ! Calculate
                dcbb = ((ht - hcbb)/(ht - STD_HT))**(1.0/beta)*dbh
            end if
            ! Set
            tree%diam_canht = dcbb
        else
            ! Shrub form
            ! Just set to dbh
            tree%diam_canht = tree%diam_bht
        end if

    end subroutine stem_shape

    !:.........................................................................:

    subroutine forska_height(tree)
        !
        !  Calculates total tree height (m)
        !  Equation adapted from FORSKA model Leemans & Prentice 1989
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    10/10/20     A. C. Foster        Update for shrub forms
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        ! Data dictionary: local variables
        real :: dbh      ! Diameter at breast height (cm)
        real :: s        ! Growth parameter
        real :: hmax     ! Average maximum height (m)
        real :: delta_ht ! Temporary variable

        dbh = tree%diam_bht
        hmax = tree%spec_ptr%max_ht
        s = tree%spec_ptr%s

        if (tree%form <= 2) then
            ! Tree or tree-like shrub form
            delta_ht = hmax - STD_HT
            tree%forska_ht = STD_HT + delta_ht*(1.0 - exp(-(s*dbh/delta_ht)))
        else
            ! Shrub form
            delta_ht  = hmax
            tree%forska_ht = delta_ht*(1.0 - exp(-(s*dbh/delta_ht)))
        end if

    end subroutine forska_height

    !:.........................................................................:

    subroutine biomass_c(tree)
        !
        !  Calculates aboveground woody biomass (tC)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    10/10/20     A. C. Foster        Update for shrub forms
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        ! Data dictionary: local varibles
        type (TreeData) :: tree_0    ! Local copy of tree object
        real            :: dbh       ! Diameter at breask height (cm)
        real            :: ht        ! Tree height (m)
        real            :: hcbb      ! Clear branch bole height (m)
        real            :: hrt       ! Rooting depth (m)
        real            :: bulk      ! Wood bulk density (t/m3)
        real            :: stembc    ! Stem biomass (tC)
        real            :: branchbc  ! Branch biomass (tC)
        real            :: abovegr_c ! Aboveground woody biomass (tC)
        real            :: root_c    ! Root biomass (tC)

        ht = tree%forska_ht
        hrt = tree%spec_ptr%rootdepth

        if (tree%form <= 2) then

            ! Tree or tree-like shrub form

            ! Copy tree object
            call copy_tree(tree_0, tree)

            ! Calculate diameter at cbb
            call stem_shape(tree)

            ! Calculate stem biomass (tC)
            stembc   = stem_biomass_c(tree)

            ! Calculate branch biomass (tC)
            branchbc = branch_biomass_c(tree)

            ! Aboveground woody biomass is branch + stem
            abovegr_c = stembc + branchbc

            ! Calculate root biomass (tC)
            root_c    = stembc*hrt/ht + branchbc

            ! Set values
            tree%biomC = abovegr_c + root_c
            tree%stemC = stembc
            tree%branchC = branchbc
            tree%rootC = root_c

        else
            ! Shrub form

            ! Calculate woody biomass (tC)
            abovegr_c = woody_biomass(tree)

            ! Calculate root biomass (tC)
            root_c = 0.3127*abovegr_c*(hrt/ht)

            ! Set values
            tree%biomC = abovegr_c
            tree%stemC = abovegr_c*0.5
            tree%branchC = abovegr_c*0.5
            tree%rootC = root_c

        end if

    end subroutine biomass_c

    !:.........................................................................:

    subroutine biomass_n(tree)
        !
        !  Calculates stem and branch N biomass (tN)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        ! Use stem C:N ratio
        tree%biomN = tree%biomC/STEM_C_N

    end subroutine biomass_n

    !:.........................................................................:

    subroutine leaf_biomass_c(tree)
        !
        !  Calculates leaf biomass (tC)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        tree%leaf_bm = leaf_area(tree)*tree%spec_ptr%leafarea_c*2.0

    end subroutine leaf_biomass_c

    !:.........................................................................:

    real function leaf_area(tree)
        !
        !  Calculates leaf area of a tree (m2)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: tree ! Tree object

        ! Data dictionary: local variables
        real :: dcbb    ! Diameter at clear branch bole height
        real :: leafd_a ! Scalar parameter for DCBB:leaf area relationship

        dcbb = tree%diam_canht
        leafd_a = tree%spec_ptr%leafdiam_a

        leaf_area = (dcbb*dcbb)*leafd_a

    end function leaf_area

    !:.........................................................................:

    subroutine age_survival(tree, a_survive)
        !
        !  Random death check. Trees can reach maximum age as a probability of
        !  1% ('often seen'), 0.1% ('some survive'), or 0.01% ('few seen'), as
        !  set by input age_tol parameter (1-3)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure

        ! Data dictionary: constants
        ! Age-tolerance related parameters
        ! 1% ---> 4.605  0.1% ---> 6.908  0.01% ---> 11.5129
        real, dimension(3), parameter :: CHECK = [4.605, 6.908, 11.51]

        ! Data dictionary: calling arguments
        class(TreeData), intent(in)  :: tree      ! Tree object
        logical,         intent(out) :: a_survive ! Does tree survive?

        ! Data dictionary: local variables
        integer :: k        ! Probability of surviving to max age (1-3; 3 = least likely)
        real    :: agemax   ! Average maximum age (years)
        real    :: rand_val ! Random number (uniform)

        k = tree%spec_ptr%age_tol
        agemax = tree%spec_ptr%max_age

        ! Get random uniform between 0 and 1
        rand_val = urand()

        ! Check against parameter
        if (rand_val < CHECK(k)/agemax) then
            a_survive = .false.
        else
            a_survive = .true.
        end if

    end subroutine age_survival

    !:.........................................................................:

    subroutine growth_survival(tree, g_survive)
        !
        !  Death check by low growth stress.
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure

        ! Data dictionary: constants
        ! Values to check against (depends on stress_tol parameter)
        real, dimension(5), parameter :: CHECK = [0.21, 0.24, 0.27, 0.30, 0.33]

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: tree      ! Tree object
        logical, intent(out)        :: g_survive ! Does tree survive?

        ! Data dictionary: local variables
        integer :: k        ! Relative stress tolerance (1-5; 5 = least tolerant)
        real    :: rand_val ! Random value (uniform)

        k = tree%spec_ptr%stress_tol

        ! Get random uniform between 0 and 1.0
        rand_val = urand()

        ! Check against parameter
        if (tree%mort_marker .and. rand_val < CHECK(k)) then
            g_survive = .false.
        else
            g_survive = .true.
        end if

    end subroutine growth_survival

    !:.........................................................................:

    subroutine newplant_survival(tree, k, p_survive)
      class(TreeData), intent(in) :: tree
      real, intent(in)         :: k    !viability
      logical, intent(out)        :: p_survive
      real    :: rand_val ! Random value (uniform)


      rand_val = urand()
      if(rand_val > k) then
         p_survive = .false.
       else
         p_survive = .true.
       end if

    end subroutine newplant_survival

    !:.........................................................................:

    subroutine fire_survival(tree, I_surf, tau_l, active_crowning, CFB,        &
        f_survive)
    !
    !  Calculates survival of a tree from a fire event.
    !  Adapted from Thonicke et al. 2010 Biogeosciences 7:1991-2010
    !
    !  Record of revisions:
    !      Date       Programmer          Description of change
    !      ====       ==========          =====================
    !    05/01/19     A. C. Foster        Original Code

    ! Data dictionary: constants
    real, parameter :: I_EXP = 0.667 ! Exponent for intensity - scorch height relationship
    real, parameter :: R_CK = 1.0    ! Resistance factor for crown scorch survival
    real, parameter :: P = 3.0       ! Parameter based on defoliation from crown scorch

    ! Data dictionary: calling arguments
    class(TreeData), intent(inout) :: tree            ! Tree object
    real,            intent(in)    :: I_surf          ! Surface fire intensity (kW/m)
    real,            intent(in)    :: tau_l           ! Residence time of fire (min)
    logical,         intent(in)    :: active_crowning ! Active crown fire?
    real,            intent(in)    :: CFB             ! Crown fraction burnt
    logical,         intent(out)   :: f_survive       ! Does tree survive?

    ! Data dictionary: local variables
    real :: cl       ! Crown length (m)
    real :: ck       ! Proportion of crown scorched (0-1)
    real :: bt       ! Bark thickness (cm)
    real :: tau_c    ! Critical time for cambial damage (min)
    real :: sh       ! Scorch height (m)
    real :: tau_frac ! Temporary varible for calculating pm_tau
    real :: pm_tau   ! Probability of mortality from cambial damage
    real :: pm_CK    ! Probability of mortality from crown scorch
    real :: pm_fire  ! Probability of mortality from fire
    real :: rand_val ! Random value (uniform)

    if (.not. active_crowning) then

        ! Tree *may* not be totally scorched

        ! Calculate crown length (m)
        cl = tree%forska_ht - tree%canopy_ht

        ! Calculate scorch height (m)
        sh = tree%spec_ptr%F_i*(I_surf**I_EXP)

        ! Calculate proportion crown scorched (0 to 1)
        ck = (sh - tree%forska_ht + cl)/cl
        ck = max(0.0, min(1.0, ck))
        tree%CK = ck

    else if (active_crowning) then

        ! Scorching from active crown fire
        ck = (sh - tree%forska_ht + cl)/cl
        ck = max(0.0, min(1.0, ck))
        ck = max(ck, CFB)
        tree%CK = ck
    end if

    ! Calculate bark thickness (cm)
    bt = tree%diam_bht*tree%spec_ptr%bark_thick

    ! Calculate critical time for cambial damage (min)
    tau_c = 2.9*(bt**2)

    ! Calculate probability of mortality due to cambial damage
    tau_frac = tau_l/tau_c
    if (tau_frac <= 0.22) then
        pm_tau = 0.0
    else if (tau_frac < 0.22 .and. tau_frac < 2.0) then
        pm_tau = 0.563*tau_frac - 0.125
    else if (tau_frac >= 2.0) then
        pm_tau = 1.0
    end if

    ! Calculate probability of mortality due crown damage
    pm_CK = R_CK*(ck**P)

    ! Calculate overall probability of mortality
    pm_fire = pm_tau + pm_CK - pm_tau*pm_CK

    ! Determine if tree killed by fire
    rand_val = urand()
    if (rand_val < pm_fire) then
        f_survive = .false.
    else
        f_survive = .true.
    endif

    ! Set values
    tree%pm_tau = pm_tau
    tree%pm_CK = pm_CK

end subroutine fire_survival
    !:.........................................................................:

    subroutine wind_survival(tree, w_survive)
        !
        !  Calculates survival from a windthrow event
        !  Adapted from Rich et al. 2007 Journal of Ecology 95:1261-1273
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    08/01/15     A. C. Foster        Original Code

        ! Data dictionary: calling arguments
        class(TreeData), intent(in)  :: tree      ! Tree object
        logical,         intent(out) :: w_survive ! Does tree survive?

        ! Data dictionary: local variables
        real :: logit    ! Temporary variable
        real :: pwind    ! Probability of mortality from windthrow
        real :: rand_val ! Random number (uniform)

        ! Probability of mortality from windthrow
        ! Equation from Rich et al. (2007)
        logit = 0.75*log(tree%diam_bht)
        pwind = 1.0/(1.0+exp(-logit))

        ! Check for survival
        rand_val=urand()
        if (rand_val < pwind) then
            w_survive = .false.
        else
            w_survive = .true.
        endif


    end subroutine wind_survival

    !:.........................................................................:

    subroutine selective_cut(tree, dbh_thresh, cutting_survive)
        !
        !  Calculates whether or not a tree is selectively cut, based on
        !  input dbh threshold
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    02/01/21     A. C. Foster        Original Code

        ! Data dictionary: calling arguments
        class(TreeData), intent(in)  :: tree            ! Tree object
        real,            intent(in)  :: dbh_thresh      ! DBH threshold for cutting
        logical,         intent(out) :: cutting_survive ! Does tree survive?

        if (tree%diam_bht >= dbh_thresh) then
            cutting_survive = .false.
        else
            cutting_survive = .true.
        end if

    end subroutine selective_cut

    !:.........................................................................:

    subroutine max_growth(tree)
        !
        !  Calculates maximum DBH increment growth given "optimal" environmental
        !   conditions
        !  Based on equation from Botkin et al. 1972 Journal of Ecology 60(3):849
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    10/10/20     A. C. Foster        Update for shrub forms

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: tree ! Tree object

        real :: dbh  ! Diameter at breast height (cm)
        real :: ht   ! Tree height (m)
        real :: dmax ! Average maximum dbh (cm)
        real :: hmax ! Average maximum height (m)
        real :: g    ! Growth parameter
        real :: s    ! Growth parameter
        real :: temp ! Temporary variable

        dbh = tree%diam_bht
        ht = tree%forska_ht
        dmax = tree%spec_ptr%max_diam
        hmax = tree%spec_ptr%max_ht
        g = tree%spec_ptr%g
        s = tree%spec_ptr%s

        if (tree%form <= 2) then
            ! Tree or tree-like shrub
            temp = s*exp(-s*dbh/(hmax - STD_HT))*dbh
            tree%diam_opt = g*dbh*(1.0 - dbh*ht/dmax/hmax)/(2.0*ht + temp)
        else
            ! Shrub
            temp = s*exp(-s*dbh/(hmax))*dbh
            tree%diam_opt = g*dbh*(1.0 - dbh*ht/dmax/hmax)/(2.0*ht + temp)
        end if

    end subroutine max_growth

    !:.........................................................................:

    real function stem_biomass_c(tree)
        !
        !  Calculates stem biomass (tC)
        !   Notes from previous model version:
        !       bstemC = bc(bd, 1.0, bulk_density)
        !       Here bd was "basal diameter" - diameter at clear branch bole
        !       height, but for a clear branch bole height of 0.0, so it was
        !       equal to dbh
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !    01/28/21     A. C. Foster        Update to use dbh of tree, without
        !                                       setting cbb to 0.0
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: tree ! Tree object

        real :: dbh          ! Diameter at breast height (cm)
        real :: ht           ! Tree height (m)
        real :: bulk_density ! Wood bulk density (t/m3)
        real :: beta         ! Stem shape parameter
        real :: temp         ! Temporary variable

        dbh = tree%diam_bht
        ht = tree%forska_ht
        bulk_density = tree%spec_ptr%wood_bulk_dens
        beta = tree%spec_ptr%beta

        temp = 1E-4*PI*0.25*B_TO_C*bulk_density*beta/(beta + 2.0)
        stem_biomass_c = temp*dbh*dbh*ht

    end function stem_biomass_c

    !:.........................................................................:

    real function woody_biomass(tree)
        !
        !  Calculates aboveground woody biomass (tC) of a shrub form
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    10/10/20     A. C. Foster        Original code

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: tree ! Tree object

        ! Data dictionary: local variables
        real :: dbh          ! Shrub diameter (cm)
        real :: ht           ! Shrub height (m)
        real :: bulk_density ! Wood bulk density (t/m3)
        real :: beta         ! Stem shape parameter
        real :: temp         ! Temporary variable

        dbh = tree%diam_bht
        ht = tree%forska_ht
        bulk_density = tree%spec_ptr%wood_bulk_dens
        beta = tree%spec_ptr%beta

        temp = 1E-4*PI*0.25*B_TO_C*bulk_density*beta/(beta + 2.0)
        woody_biomass = temp*dbh*dbh*ht*1.8

    end function woody_biomass

    !:.........................................................................:

    real function branch_biomass_c(tree)
        !
        !  Calculates branch biomass of a tree (tC)
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/05     Y. Xiaodong         Original Code
        !    07/27/12     K. Holcomb          Update to OOP structure
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: tree ! Tree object

        ! Data dictionary: local variables
        real :: dcbb         ! Diameter at clear branch bole height (cm)
        real :: hcbb         ! Clear branch bole height (m)
        real :: ht           ! Height (m)
        real :: bulk_density ! Wood bulk density (t/m3)
        real :: beta         ! Stem shape parameter
        real :: temp         ! Temporary variable

        dcbb = tree%diam_canht
        hcbb = tree%canopy_ht
        ht = tree%forska_ht
        bulk_density = tree%spec_ptr%wood_bulk_dens
        beta = tree%spec_ptr%beta

        temp = 1E-4*PI*0.25*B_TO_C*bulk_density*(2.0/(beta + 2.0)-0.33)
        branch_biomass_c = temp*dcbb*dcbb*(ht-hcbb)

    end function branch_biomass_c

    !:.........................................................................:

    subroutine get_diam_category(self, diam_category)
        !
        !  Gets the diameter class of a tree
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    07/27/12     K. Holcomb          Original code
        !    01/28/21     A. C. Foster        Updated to use parameters for DBH
        !                                       bins
        !

        ! Data dictionary: calling arguments
        class(TreeData),         intent(in)  :: self          ! Tree object
        integer, dimension(NHC), intent(out) :: diam_category ! List with 1 in correct bin

        ! Data dictionary: local variables
        real    :: dbh ! Diameter at breast height (cm)
        integer :: i   ! Looping index

        dbh = self%diam_bht

        ! Reset to 0
        diam_category = 0

        ! Find correct bin
        do i = 1, NHC - 1
            if (dbh > DBC_MIN(i) .and. dbh <= DBC_MAX(i)) then
                diam_category(i) = 1
            end if
        end do
        if (dbh > DBC_MIN(NHC)) then
            diam_category(NHC) = 1
        end if

    end subroutine get_diam_category

    !:.........................................................................:

    subroutine timber_get_diam_category(self, diam_category, tree_volume, comp)
        !
        !  Gets the diameter class of a tree
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    07/27/12     K. Holcomb          Original code
        !    01/28/21     A. C. Foster        Updated to use parameters for DBH
        !                                       bins
        !    05/01/22     S. Sundquist        Added timber volume calculation

        ! Data dictionary: calling arguments
        class(TreeData),         intent(in)  :: self          ! Tree object
        integer, dimension(NHC), intent(out) :: diam_category ! List with 1 in correct bin
        real, dimension(2),      intent(out) :: tree_volume   ! CF and BF calculation for individual tree
        character(len = MAX_NLEN)            :: comp          ! Species flag for volume calculation

        ! Data dictionary: local variables
        real    :: dbh ! Diameter at breast height (cm)
        real    :: dbh_in ! Diameter at breast height (in)
        real    :: ht_ft  ! Tree total height (ft)
        integer :: i   ! Looping index

        dbh = self%diam_bht

        ! Reset to 0
        diam_category = 0
        tree_volume = 0.0

        ! Find correct bin
        do i = 1, NSC - 1
            if (dbh > timber_DBC_MIN(i) .and. dbh <= timber_DBC_MAX(i)) then
                diam_category(i) = 1
            end if
        end do
        if (dbh > timber_DBC_MIN(NSC)) then
            diam_category(NSC) = 1
        end if

        ! Volume formulas from Timber Inventory of State Forest Lands in the
        ! Tanana Valley 2013; Table 6, page 18

        dbh_in = dbh*0.3937 !cm to in (aren't americans wonderful!)
        ht_ft = self%forska_ht*3.2808 !and m to ft
        if(diam_category(2) > 0 .or. diam_category(3) > 0) then ! Tree gets added to cubic volume
          !NOR-6 formulas
          if(comp == 'PICEglau' .or. comp == 'PICEmari') then
            tree_volume(1) = -2.055 + 0.2982*dbh_in + 0.00181*(dbh_in**2)*ht_ft
          else if(comp == 'POPUtrem') then
            tree_volume(1) = -0.5553 - 0.02216*(dbh_in**2) + 0.00246*(dbh_in**2)*ht_ft
          else if(comp == 'POPUbals') then
            tree_volume(1) = -3.2187 + 0.8281*dbh_in - 0.05908*(dbh_in**2) - 0.01985*ht_ft + 0.00199*(dbh_in**2)*ht_ft
          else if(comp == 'BETUneoa') then
            tree_volume(1) = -2.5767 + 0.9524*dbh_in - 0.10446*(dbh_in**2) - 0.03303*ht_ft + 0.00282*(dbh_in**2)*ht_ft
          else
            tree_volume(1) = 0.0
          end if

        end if
        if(diam_category(3) > 0) then ! Tree gets added to boardfeet volume
          ! Haack 1963 says deciduous equations are meant for dbh > 11in but TV document says they applied to dbh > 9in
          if(comp == 'PICEglau' .or. comp == 'PICEmari') then
            ! PNW-59 formula
            tree_volume(2) = 39.71 + 4.2659*dbh_in - 0.55865*(dbh_in**2) - 1.1184*ht_ft + 0.016113*(dbh_in**2)*ht_ft - 437.92/(dbh_in**2)
          else if(comp == 'POPUtrem' .or. comp == 'BETUneoa') then
            ! NOR-5 formula
            tree_volume(2) = -27.263 + 0.00995*(dbh_in**2)*ht_ft
          else if(comp == 'POPUbals') then
            ! NOR-5 formula
            tree_volume(2) = -46.7415 + 0.00956*(dbh_in**2)*ht_ft
          else
            tree_volume(2) = 0.0
          end if

        end if

    end subroutine timber_get_diam_category

    !:.........................................................................:

    subroutine write_tree_csv(self, tree_unit)
        !
        !  Helper subroutine for writing tree-level data to the output Tree file
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    05/01/16     A. C. Foster        Original code
        !

        use csv_file

        ! Data dictionary: calling arguments
        class(TreeData), intent(in) :: self      ! Tree object
        integer,         intent(in) :: tree_unit ! Unit number for output tree file

        ! Data dictionary: local variables
        real    :: tla ! Leaf area (m2)
        integer :: i   ! Looping index

        ! Get leaf area
        tla = leaf_area(self)

        ! Write out attributes
!        call csv_write(tree_unit, trim(adjustl(self%spec_ptr%genus_name)),     &
!            .false.)
        call csv_write(tree_unit, self%spec_ptr%unique_id, .false.)
!        call csv_write(tree_unit, self%row, .false.)
!        call csv_write(tree_unit, self%col, .false.)
        call csv_write(tree_unit, self%id, .false.)
        call csv_write(tree_unit, self%tree_age, .false.)
        call csv_write(tree_unit, self%planted, .false.)
        call csv_write(tree_unit, self%diam_bht, .false.)
!        call csv_write(tree_unit, self%diam_canht, .false.)
        call csv_write(tree_unit, self%forska_ht, .true.)
!        call csv_write(tree_unit, self%canopy_ht, .false.)
!        call csv_write(tree_unit, self%leaf_bm, .false.)
!        call csv_write(tree_unit, tla, .false.)
!        call csv_write(tree_unit, self%biomC, .false.)
!        do i = 1, FC_NUM-1
!            call csv_write(tree_unit, self%env_resp(i), .false.)
!        end do
!        call csv_write(tree_unit, self%env_resp(FC_NUM), .true.)

    end subroutine write_tree_csv

    !:.........................................................................:

    subroutine delete_tree(self)
        !
        !  Deletes a tree object, nullifying its pointers
        !
        !  Record of revisions:
        !      Date       Programmer          Description of change
        !      ====       ==========          =====================
        !    02/02/21     A. C. Foster        Original Code
        !

        ! Data dictionary: calling arguments
        class(TreeData), intent(inout) :: self

        if (associated(self%spec_ptr)) then
            nullify(self%spec_ptr)
        end if

    end subroutine delete_tree

    !:.........................................................................:

    subroutine init_age(self)
      class(TreeData), intent(inout) :: self

    ! This creates infinite loop when tree is too small in specieslist params 
      real :: dmax ! Average maximum dbh (cm)
      real :: hmax ! Average maximum height (m)
      real :: delta_ht
      real :: g    ! Growth parameter
      real :: s    ! Growth parameter
      real :: tg ! temporary growth parameter
      real :: tdbh ! temporary dbh
      real :: tht !temporary height
      integer :: age ! Age counter assuming optimal growth

      dmax = self%spec_ptr%max_diam
      hmax = self%spec_ptr%max_ht
      g = self%spec_ptr%g
      s = self%spec_ptr%s

      tdbh = 1.0
      tg = 0.0
      tht = STD_HT
      age = 0

      do while(tdbh < self%diam_bht)
        tg = s*exp(-s*tdbh/(hmax - STD_HT))*tdbh
        delta_ht = hmax - STD_HT
        tht = STD_HT + delta_ht*(1.0 - exp(-(s*tdbh/delta_ht)))
        tdbh = tdbh + g*tdbh*(1.0 - tdbh*tht/dmax/hmax)/(2.0*tht + tg)
        age = age + 1
      end do

      self%tree_age = age
    end subroutine init_age

end module Tree
