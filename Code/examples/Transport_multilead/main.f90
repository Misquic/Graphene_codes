! This program calculates transport for "simply" twisted graphene flake. It means
! that twisted boundary is not taken into account, there is just flip of magnetic
! field. System looks like this (x and y used are for 1 view reference, not physical coordinates)
!    (unfolded view)
! Y   ________________________
! ^  |l|  top              |l|
! |  |e|  B = (0,0,Bz)     |e|
! |  |a|___________________|a|
! |  |d|  bottom           |d|
! |  |1|  B = (0,0,-Bz)    |2|
! |  |_|___________________|_|
! *--------------> X
!
! Coupling between layers is taken into account with Bilayer class, transport is
! calculated using Bubel

program main
  use modscatter
  use FortUtils
  use modshape
  use TransportUtils
  use, intrinsic :: iso_c_binding
  implicit none

  ! Bubel objects
  type(qscatter) :: qt

  ! Shared variables used by internal procedures
  doubleprecision :: middle_y             ! Y coordinate of fold
  doubleprecision :: middle_x             !
  doubleprecision :: yBoundLower          ! Y coordinate where Vg goes from Vgb to linear region
  doubleprecision :: yBoundUpper          ! Y coordinate where Vg goes from linear region to Vgt

  doubleprecision :: Bau                        ! in au
  doubleprecision :: Vgt, Vgb, E0t, E0b, nt, nb ! result from Bilayer
  doubleprecision :: Ef                         ! Fermi energy for calculations
  integer, parameter :: numLeads = 4
  doubleprecision, dimension(numLeads, numLeads) :: T ! Transmission matrix
  integer         :: nx = 90                    ! numbers of atoms / 2 in x direction
  ! integer         :: nx = 20                    ! numbers of atoms / 2 in x direction
                                                ! results in about 196 nm
  integer         :: ny = 122                    ! ~numbers of atoms / 2 in y direction (keep even)
  ! integer         :: ny = 40                    ! ~numbers of atoms / 2 in y direction (keep even)
                                                ! results in about 170 nm
  integer         :: from = 1, to = 1           ! lead indexes for writing transmissions to a file
  integer :: n
  doubleprecision :: r = 0
  ! lead storage
  type(qshape) :: leadShapes(numLeads)
  doubleprecision :: leadTrans(3,numLeads)

!!!!!!!!!!!!!!!!!!!!!!!! main function !!!!!!!!!!!!!!!!!!!!!!!

  call random_seed(size = n)
  print*, "seed_size: ", n
  allocate(seed(n))

  call parseArguments()
  nx = nx * (16 / sf)
  ny = ny * (16 / sf) + 1
  call random_seed(put = seed)

  do n = 1,100
    call random_number(r)
  enddo

  call createSystem()
  if (save_system) then
    call qt%save_system(trim(results_dir)//"/system.xml")
  endif

  if (run_transport) then
    ! Calculate at specific Fermi energy
    Ef = 0.000001D0 ! eV
    ! Ef = 0.00000D0 ! eV
    print*,"========================================"
    print*,"Calculating transport at Ef = ",Ef," eV"
    print*,"========================================"
    T = solveTransportMultilead(qt, Ef, numLeads)
    print *, T
    if (save_densities) then
      call calculateElectronDensity(qt)
    endif

    ! Write single T (lead 1) to file
    open(unit=101, file=trim(results_dir)//"/single_T.dat")
    write(101,"(A)") "Ef[au],T[-],E0t[au],E0b[au],Vgt[au],Vgb[au],nt[au],nb[au]"
    write(101,"(g0,',',g0,',',g0,',',g0,',',g0,',',g0,',',g0,',',g0)") &
      Ef * eV2au, sum(qt%Tn(:)), E0t, E0b, Vgt, Vgb, nt, nb
    close(101)

    ! Write T matrix to file
    open(unit=101, file=trim(results_dir)//"/Transmissions.csv")
    write(101, "(A, I)") "to, from, Transmission, numLeads= ", numLeads
    do from = 1, numLeads
      do to = 1, numLeads
        write(101, "(I,',',I,',',g0)") to, from, T(to, from)
      enddo
    enddo
    close(101)

    if (run_energyScan) then
      print*,""
      print*,"========================================"
      print*,"Performing energy scan..."
      print*,"========================================"
      call performEnergyScan(qt, connect)
    endif


    print*,"Calculation complete!"
  endif

contains

! --------------------------------------------------------------------------------------------------
! Create System
! --------------------------------------------------------------------------------------------------
  subroutine createSystem()
    use Bilayer_interface
    use modscatter
    use modsys
    use modshape
    use modunits
    use, intrinsic :: iso_c_binding
    implicit none

    type(qatom) :: qa
    ! graphene parameters
    doubleprecision,parameter :: alpha30                = 30.0/180.0*M_PI
    integer,parameter         :: atomA = 1, atomB = 2 ! sublattices flags
    doubleprecision,parameter :: carbon_carbon_dist = 0.142 ! nm
    doubleprecision,parameter :: geometric_unit = carbon_carbon_dist * sqrt(3.0)
    doubleprecision,parameter :: geometric_unit2au = geometric_unit * nm2au

    doubleprecision           :: vecs_armchair(2,2)     = (/ (/ 1.0D0,0.0D0 /), (/ sin(alpha30), cos(alpha30) /) /)
    doubleprecision           :: atoms_armchair(2,2)    = (/ (/ 0.0D0,0.0D0 /), (/ 0.0D0, one_over_sqrt_3 /) /)
    doubleprecision           :: pos_offset_armchair(2) = (/ -sin(alpha30), -cos(alpha30) /)

    ! local variables
    integer         :: i, j, atom           ! loop variables
    doubleprecision :: atom_pos(3)
    doubleprecision :: pos_max(2) = (/ 0.0D0, 0.0D0/)
    doubleprecision :: pos_min(2) = (/ 0.0D0, 0.0D0/)
    doubleprecision :: x_min      = 0.0D0, x_max = 0.0D0
    doubleprecision :: y_min      = 0.0D0, y_max = 0.0D0
    type(c_ptr)     :: bilayer
    integer         :: currentLead = 1

    doubleprecision :: yLeadWidth = 0! width of lead
    doubleprecision :: cellSize = 0! width of lead
    doubleprecision :: xLenCut = 0! how much does lead stand out

! --------------------------------------------------------------------------------------------------
    vecs_armchair = vecs_armchair * sf * geometric_unit2au
    atoms_armchair = atoms_armchair * sf * geometric_unit2au
    pos_offset_armchair = pos_offset_armchair * sf * geometric_unit2au

    ! some magic to have nice edges
    pos_max = atoms_armchair(:, 2) + nx * vecs_armchair(:,1) + ny * vecs_armchair(:,2) - 0.001 ! 0.001 is to ommit numerical errors
    pos_max(1) = pos_max(1) - 2 * (ny / 2) * vecs_armchair(1,2)
    pos_max(2) = pos_max(2) + 2 * pos_offset_armchair(2)
    pos_min = pos_offset_armchair * 0.5 + 0.001 ! 0.001 is to ommit numerical errors

    ! settup of Lead width
    yLeadWidth = (pos_max(2) - pos_min(2)) * 0.5 * 0.3 ! 0.5 -> 2 leads on one side, each takes 0.4 of half width
    cellSize = 3 * carbon_carbon_dist * nm2au * sf ! lead must be multiple of size to make good edges of cutss
    yLeadWidth = cellSize * (int(yLeadWidth / cellSize) + 1) ! round up
    xLenCut = (pos_max(1) - pos_min(1)) * 0.05

    x_min = atoms_armchair(1,1)
    x_max = atoms_armchair(1,1)
    y_min = atoms_armchair(1,2)
    y_max = atoms_armchair(1,2)

    print*, "pos_max ", pos_max / nm2au
    print*, "pos_min ", pos_min / nm2au

    !------------------------------------------ Bilayer --------------------------------------------
    bilayer = Bilayer_constructor_default()

    call Bilayer_countAll_B(bilayer, Vt, Vb, Bz, Vgt, Vgb, E0t, E0b, nt, nb)

    ! convert back to au
    Vgt = Vgt * eV2au
    Vgb = Vgb * eV2au
    E0t = E0t * eV2au
    E0b = E0b * eV2au
    nt = nt * inv_cmsq2au
    nb = nb * inv_cmsq2au

    !---------------------------------------- Lattice ----------------------------------------------
    call qt%init_system()
    QSYS_DEBUG_LEVEL = 0
    QSYS_FORCE_SCHUR_DECOMPOSITION  = .false. ! don't use schur method so its quicker
    ! Try with .true.

    ! Generate atoms positions
    do i = 0, nx
      do j = 0, ny
        do atom = atomA, atomB
          atom_pos(1:2) = atoms_armchair(:,atom) + & ! base position, atom choses sublattice
            i * vecs_armchair(:,1) + & ! offset in x direction for i-th "column"
            j * vecs_armchair(:,2) + & ! offset in y dir for j-th "row"
            pos_offset_armchair ! base offset

          ! works only with armchair
          atom_pos(1) = atom_pos(1) - 2 * (j / 2) * vecs_armchair(1,2) ! shift "rows" to make flake rectangular

          if (checkShape(atom_pos, pos_min, pos_max, yLeadWidth, xLenCut)) then

            x_max = max(x_max, atom_pos(1))
            x_min = min(x_min, atom_pos(1))
            y_max = max(y_max, atom_pos(2))
            y_min = min(y_min, atom_pos(2))

            call qa%init( (/atom_pos(1), atom_pos(2), 0.0D0 /), flag=atom) ! create qatom
            call qt%qsystem%add_atom(qa) ! add qatom to system
          endif
        enddo
      enddo
    enddo

    middle_x = 0.5 * (x_min + x_max)
    middle_y = 0.5 * (y_min + y_max)

    yBoundLower = middle_y - (y_max - y_min) * 0.03D0 * 0.5
    yBoundUpper = middle_y + (y_max - y_min) * 0.03D0 * 0.5

    ! yBoundLower = middle_y
    ! yBoundUpper = middle_y

    write(*, "(A,f10.5,A)"), "middle_x    ", middle_x / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "middle_y    ", middle_y / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "yBoundLower ", yBoundLower / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "yBoundUpper ", yBoundUpper / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "x_min       ", x_min / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "x_max       ", x_max / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "y_min       ", y_min / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "y_max       ", y_max / nm2au, " nm"

    !----------------------------------------- Leads -----------------------------------------------

    ! S
    ! call addXInvLeads(x_min, x_max, middle_y, yLeadWidth * 1.05 , vecs_armchair)
    ! call addXInvDiagLeads(x_min, x_max, y_min, y_max, yLeadWidth * 1.05 , vecs_armchair)

    ! H
    ! call addXInvLeads(x_min, x_max, y_max - yLeadWidth/2, yLeadWidth * 1.05 , vecs_armchair)
    ! call addXInvLeads(x_min, x_max, y_min + yLeadWidth/2, yLeadWidth * 1.05, vecs_armchair)

    ! I
    ! lead 1 and 2 are added (voltage)
    call addXInvLeads(currentLead, pos_min(1) + xLenCut, pos_max(1) - xLenCut, middle_y, yLeadWidth * 1.05, vecs_armchair)
    ! lead 3 and 4 are added (current)
    call addXInvUpDownLeads(currentLead, x_min, x_max, y_min, y_max, vecs_armchair)

    print*, "currentLead", currentLead
    !---------------------------------------- Coupling ---------------------------------------------

    ! Coupling between atoms, onsite energies
    qt%qnnbparam%distance = 0.6 * sf * geometric_unit2au
    qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE

    Bau = Bz * T2au

    call qt%qsystem%make_lattice(qt%qnnbparam, c_simple=connect)

    !------------------------------------------ Leads pt 2 -----------------------------------------
    do i = 1, numLeads
      call qt%add_lead(leadShapes(i), leadTrans(:, i))
      if (save_bands) then
        call qt%leads(i)%bands(trim(results_dir)//"/bands" // str(i) // ".dat", &
                              -M_PI / one_over_sqrt_3, +M_PI/ one_over_sqrt_3, M_PI/ one_over_sqrt_3/160.0, & !k_min, k_max, dk
                              -3.0D0 * eV2au, 3.0D0 * eV2au) !E_min, E_max
      endif
    enddo

    write(*,"(A,f8.5,A)") "nt  ", nt / inv_cmsq2au / 1e11, " 10^11 m^-2"
    write(*,"(A,f8.5,A)") "nb  ", nb / inv_cmsq2au / 1e11, " 10^11 m^-2"
    write(*,"(A,f8.5,A)") "Vgt ", Vgt / eV2au, " eV"
    write(*,"(A,f8.5,A)") "Vgb ", Vgb / eV2au, " eV"
    write(*,"(A,f8.5,A)") "ot  ", (- E0t - Vgt) / eV2au, " eV"
    write(*,"(A,f8.5,A)") "ob  ", (- E0b - Vgb) / eV2au, " eV"
    write(*,"(A,f8.5,A)") "E0t ", E0t / eV2au, " eV"
    write(*,"(A,f8.5,A)") "E0b ", E0b / eV2au, " eV"

  end subroutine
! --------------------------------------------------------------------------------------------------

  subroutine leads_init()
    use modshape
    implicit none
    integer :: i
    do i = 1, 4
      call leadShapes(i)%init(SHAPE_NONE)
      leadTrans(:,i) = 0.0D0
    end do
  end subroutine leads_init


  subroutine addLeadRect(idx, xmin, xmax, ymin, ymax, translation)
    use modshape
    implicit none
    integer, intent(in) :: idx
    doubleprecision, intent(in) :: xmin, xmax, ymin, ymax
    doubleprecision, intent(in) :: translation(3)
    if (idx < 1 .or. idx > 4) then
      print *, "addLeadRect: idx out of range", idx
      return
    end if
    call leadShapes(idx)%init_rect(SHAPE_RECTANGLE_XY, xmin, xmax, ymin, ymax)
    leadTrans(:,idx) = translation
  end subroutine addLeadRect


  integer function isInLeads(x, y)
    implicit none
    doubleprecision, intent(in) :: x, y
    doubleprecision :: vec(3)
    integer :: i
    isInLeads = 0
    vec = (/ x, y, 0.0D0 /)
    do i = 1, 4
      if (leadShapes(i)%is_inside(vec)) then
        isInLeads = i
        return
      end if
    end do
  end function isInLeads


! --------------------------------------------------------------------------------------------------
! Add 2 leads on Y sides, Y invariant?
! --------------------------------------------------------------------------------------------------
  subroutine addYInvLeads(idx, y_min, y_max, leadLength, vecs_armchair)
    use modscatter
    use modunits
    use modshape
    implicit none

    integer :: idx
    doubleprecision :: lead_translation(2) = (/ 0.0D0, 0.0D0 /)!
    doubleprecision, intent(in) :: y_min, y_max, leadLength
    doubleprecision, dimension(2,2), intent(in) :: vecs_armchair

! --------------------------------------------------------------------------------------------------

    lead_translation = (/ 0.0D0, vecs_armchair(2,2) /) * 2
    print*, "lead_translation: ", lead_translation

    ! First lead (lower Y side) -> store as lead 1
    call addLeadRect(idx, middle_x - 0.1 - leadLength / 2, &
                     middle_x + 0.1 + leadLength / 2, &
                     y_min - lead_translation(2) * 0.25, &
                     y_min + lead_translation(2) * 0.75, &
                     (/lead_translation(1), lead_translation(2), 0.0D0/))
    idx = idx + 1

    ! Second lead (upper Y side) -> store as lead 2
    call addLeadRect(idx, middle_x - 0.1 - leadLength / 2, &
                     middle_x + 0.1 + leadLength / 2, &
                     y_max - lead_translation(2) * 0.75, &
                     y_max + lead_translation(2) * 0.25, &
                     (/-lead_translation(1), -lead_translation(2), 0.0D0 /))
    idx = idx + 1

  end subroutine
! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------
! Add 2 leads on X sides, x invariant?
! --------------------------------------------------------------------------------------------------
  subroutine addXInvLeads(idx, x_min, x_max, yLeadMiddle, yLeadWidth, vecs_armchair)
    use modscatter
    use modunits
    use modshape
    implicit none

    integer :: idx
    doubleprecision :: lead_translation(2) !
    doubleprecision, intent(in) :: x_min, x_max
    doubleprecision, intent(in) :: yLeadMiddle, yLeadWidth
    doubleprecision, dimension(2,2), intent(in) :: vecs_armchair

! --------------------------------------------------------------------------------------------------

    lead_translation = (/vecs_armchair(1,1), 0.0D0/)
    print*, "lead_translation: ", lead_translation
    ! First lead (lower X)
    call addLeadRect(idx, x_min - 0.1, &
             x_min + lead_translation(1) - 0.1, &
             yLeadMiddle - 0.1 - yLeadWidth / 2, &
             yLeadMiddle + 0.1 + yLeadWidth / 2, &
             (/lead_translation(1), lead_translation(2), 0.0D0/))
    idx = idx + 1

    ! Second lead (higher X)
    call addLeadRect(idx, x_max - lead_translation(1) + 0.1, &
             x_max + 0.1, &
             yLeadMiddle - 0.1 - yLeadWidth / 2, &
             yLeadMiddle + 0.1 + yLeadWidth / 2, &
             (/-lead_translation(1), -lead_translation(2), 0.0D0 /))
    idx = idx + 1

  end subroutine
! --------------------------------------------------------------------------------------------------


  subroutine addXInvDiagLeads(idx, x_min, x_max, y_min, y_max, yLeadWidth, vecs_armchair)
    use modscatter
    use modunits
    use modshape
    implicit none

    integer :: idx
    doubleprecision :: lead_translation(2) !
    doubleprecision, intent(in) :: x_min, x_max, y_min, y_max
    doubleprecision, intent(in) :: yLeadWidth
    doubleprecision, dimension(2,2), intent(in) :: vecs_armchair

! --------------------------------------------------------------------------------------------------

    lead_translation = (/vecs_armchair(1,1), 0.0D0/)
    print*, "lead_translation: ", lead_translation
    ! First lead (lower X) lower Y
    call addLeadRect(idx, &
                     x_min - 0.1, &
                     x_min + lead_translation(1) - 0.1, &
                     y_min - 0.1 , &
                     y_min + 0.1 + yLeadWidth, &
                     (/lead_translation(1), lead_translation(2), 0.0D0 /))
    idx = idx + 1

    ! Second lead (higher X) higher Y
    call addLeadRect(idx, &
                     x_max - lead_translation(1) + 0.1, &
                     x_max + 0.1, &
                     y_max - 0.1 - yLeadWidth, &
                     y_max + 0.1, &
                     (/-lead_translation(1), -lead_translation(2), 0.0D0 /))
    idx = idx + 1

  end subroutine
! --------------------------------------------------------------------------------------------------


  ! add leads that span whole Y side to simulate fold with leads on two layers
  subroutine addXInvUpDownLeads(idx, x_min, x_max, y_min, y_max, vecs_armchair)
    use modscatter
    use modunits
    use modshape
    implicit none

    integer :: idx
    doubleprecision :: lead_translation(2) !
    doubleprecision, intent(in) :: x_min, x_max, y_min, y_max
    doubleprecision, dimension(2,2), intent(in) :: vecs_armchair

! --------------------------------------------------------------------------------------------------

    lead_translation = (/vecs_armchair(1,1), 0.0D0/)
    print*, "lead_translation: ", lead_translation
    ! First lead (lower X spanning Y)
    call addLeadRect(idx, x_min - 0.1, &
                     x_min + lead_translation(1) - 0.1, &
                     y_min - 0.1 , &
                     y_max + 0.1, &
                     (/lead_translation(1), lead_translation(2), 0.0D0/))
    idx = idx + 1

    ! Second lead (higher X)
    call addLeadRect(idx, x_max - lead_translation(1) + 0.1, &
                     x_max + 0.1, &
                     y_min - 0.1, &
                     y_max + 0.1, &
                     (/-lead_translation(1), -lead_translation(2), 0.0D0 /))
    idx = idx + 1

  end subroutine
! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------
! Check if atom position is within bounds
! --------------------------------------------------------------------------------------------------
  logical function checkShape(atom_pos, pos_min, pos_max, yLeadWidth, xLenCut)
    implicit none

    doubleprecision, intent(in) :: atom_pos(3)
    doubleprecision, intent(in) :: pos_min(2)
    doubleprecision, intent(in) :: pos_max(2)
    doubleprecision, intent(in) :: yLeadWidth
    doubleprecision, intent(in) :: xLenCut

    doubleprecision :: middle_y
    doubleprecision :: range(2)

    range = pos_max - pos_min
    middle_y = (pos_max(2) + pos_min(2)) / 2

    ! contacts width is lover then whole width on Y

    !    ______
    !  _|      |_
    ! |          |
    ! |_        _|
    !   |______|

    ! checkShape = &
    !   (atom_pos(2) > pos_min(2) + range(2) * 0.1 .and. &
    !    atom_pos(2) < pos_max(2) - range(2) * 0.1 .and. &
    !    atom_pos(1) > pos_min(1) .and. &
    !    atom_pos(1) < pos_max(1)) &
    !   .or. &
    !   (atom_pos(1) > pos_min(1) + range(1) * 0.05 .and. &
    !    atom_pos(1) < pos_max(1) - range(1) * 0.05 .and. &
    !    atom_pos(2) > pos_min(2) .and. &
    !    atom_pos(2) < pos_max(2))

    !  ________________
    ! |_              _|
    !  _|            |_
    ! |________________|

    ! checkShape = &
    !   ((atom_pos(1) > pos_min(1)) .and. (atom_pos(1) < pos_max(1)) .and. & ! outer x layer
    !    (atom_pos(2) > pos_min(2)) .and. (atom_pos(2) < pos_max(2))) &      ! outer y layer
    !   .and. &
    !   (.not.(((atom_pos(1) < pos_min(1) + xLenCut) .or. &    ! left or cut x
    !           (atom_pos(1) > pos_max(1) - xLenCut)) .and. &  ! right cut x
    !          (atom_pos(2) > pos_min(2) + yLeadWidth) .and. & ! down side of middle section y
    !          (atom_pos(2) < pos_max(2) - yLeadWidth)))       ! uper side of middle section y



    !   |     1      |
    !    ______________
    !   |              | 3
    !   |             _|
    !  _|            |_
    ! |_ ---fold----  _| 2
    !  _|            |
    ! |              |
    ! |______________|   4

    ! checkShape = &
    !   ((atom_pos(1) > pos_min(1) + xLenCut) .and. (atom_pos(1) < pos_max(1) - xLenCut) .and. &
    !    (atom_pos(2) > pos_min(2)) .and. (atom_pos(2) < pos_max(2))) .or. & ! 1
    !   ((atom_pos(1) > pos_min(1)) .and. (atom_pos(1) < pos_max(1)) .and. &
    !    (atom_pos(2) > middle_y - yLeadWidth / 4) .and. (atom_pos(2) < middle_y + yLeadWidth / 4)) .or. & ! 2
    !   ((atom_pos(1) > pos_min(1)) .and. (atom_pos(1) < pos_max(1) - xLenCut) .and. &
    !    (atom_pos(2) > pos_min(2)) .and. (atom_pos(2) < pos_min(2) + yLeadWidth)) .or. & ! 4
    !   ((atom_pos(1) > pos_min(1) + xLenCut) .and. (atom_pos(1) < pos_max(1)) .and. &
    !    (atom_pos(2) > pos_max(2) - yLeadWidth) .and. (atom_pos(2) < pos_max(2))) ! 3



    !     |     1     |
    !  ___________________
    ! |                   | 3
    ! |___             ___|
    !    _|           |_
    !   |_ ---fold---- _| 2
    !  ___|           |___
    ! |                   |
    ! |___________________|   4

    checkShape = &
      ((atom_pos(1) > pos_min(1) + 2 * xLenCut) .and. (atom_pos(1) < pos_max(1) - 2 * xLenCut) .and. &
       (atom_pos(2) > pos_min(2)) .and. (atom_pos(2) < pos_max(2))) .or. & ! 1
      ((atom_pos(1) > pos_min(1) + xLenCut) .and. (atom_pos(1) < pos_max(1) - xLenCut) .and. &
       (atom_pos(2) > middle_y - yLeadWidth / 4) .and. (atom_pos(2) < middle_y + yLeadWidth / 4)) .or. & ! 2
      ((atom_pos(1) > pos_min(1)) .and. (atom_pos(1) < pos_max(1)) .and. &
       (atom_pos(2) > pos_min(2)) .and. (atom_pos(2) < pos_min(2) + yLeadWidth)) .or. & ! 4
      ((atom_pos(1) > pos_min(1)) .and. (atom_pos(1) < pos_max(1)) .and. &
       (atom_pos(2) > pos_max(2) - yLeadWidth) .and. (atom_pos(2) < pos_max(2))) ! 3


  end function


! --------------------------------------------------------------------------------------------------
! Calculate hoping between atoms, here we use Peierls phase
! to simulate magnetic field with gauge: A = (-Bz * y,0,0)
! i.e. B = (0,0,Bz)
! --------------------------------------------------------------------------------------------------
  logical function connect(atomA, atomB, coupling_val, atoms)
    use modcommons
    implicit none

    type(qatom) :: atomA, atomB ! flags of atoms
    type(qatom) :: atoms(:)     ! unused
    complex*16  :: coupling_val ! result coupling

    doubleprecision :: xA, yA, xB, yB
    doubleprecision :: phi
    doubleprecision :: B
    doubleprecision :: t0
    doubleprecision :: Vg
    doubleprecision :: E0
    doubleprecision :: y, x
    doubleprecision :: r
! --------------------------------------------------------------------------------------------------
    if (.false.) atoms(0)%flag = atoms(0)%flag ! supress unused variable warning
    if (.not. (atomA%flag == atomB%flag)) then
      ! hopping
      connect = .true.
      t0 = (- 3.0D0 * eV2au) / sf
      xA = atomA%atom_pos(1)
      yA = atomA%atom_pos(2)
      xB = atomB%atom_pos(1)
      yB = atomB%atom_pos(2)
      B = cosineGradient(y, Bau, yBoundUpper, yBoundLower)

      ! Peierls phase
      phi = - 0.5 * B * (yB + yA - 2 * middle_y) * (xB - xA) ! y x already in au
      coupling_val = t0 * exp(II*phi)
    else
      ! onsite
      connect = .true.
      x = atomA%atom_pos(1)
      y = atomA%atom_pos(2)
      Vg = linear(y, Vgt, Vgb, yBoundUpper, yBoundLower)
      E0 = linear(y, E0t, E0b, yBoundUpper, yBoundLower)
      coupling_val = - E0 - Vg
      if (.not. (isInLeads(x, y))) then
        ! not lead
        call random_number(r)
        ! potencjał andersona, +-5/100 eV dla sf8 -> 1/10 2/10 -----> sf4 - pomnożyć przez 2 sf ||||| rośnie AA maleje
        ! add random number +- 0.05 eV
        coupling_val = coupling_val + (2 * r - 1) * 0.1 * eV2au
        ! coupling_val = coupling_val * (1D0 + (r - 0.5D0) / 10D0)
      endif
    endif
  end function
! --------------------------------------------------------------------------------------------------

end program main

! TODO get "andersonAnmplitude" -> maks wartość zaburzenia to be a parameter
! TODO set andersonAmplitude in leads to 0
