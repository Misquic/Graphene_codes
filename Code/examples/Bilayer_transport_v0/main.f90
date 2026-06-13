! This program calculates transport for "simply" twisted graphene flake. It means
! that twisted boundary is not taken into account, there is just flip of magnetic
! field. System looks like this (x and y used are for 1 view reference, not physical coordinates)
!    (unfolded view)            !   (side view left)             (side view right)        (side view front)
! Y   ________________________  !   ____lead1___________        ____lead2___________     ___                   ___
! ^  |l|  top              |l|  !   |  ___top__gate___ |        |  ___top__gate___ |     |l|  ___top_gate___  |l|
! |  |e|  B = (0,0,Bz)     |e|  !   |  _______top_____ |     ^  | _______top_____  |     |e|_______top________|e|
! |  |a|___________________|a|  !   | /                | B = |  |                \ |     |a|                  |a|
! |  |d|  bottom           |d|  !   | \_____bottom____ |     |  | _____bottom____/ |     |d|_____bottom_______|d|
! |  |1|  B = (0,0,-Bz)    |2|  !   |  __bottom_gate__ |        |  __bottom_gate__ |     |1|  _bottom__gate_  |2|
! |  |_|___________________|_|  !   |__________________|        |__________________|     |_|                  |_|
! *--------------> X            !
!
! Coupling between layers is taken into account with Bilayer class, transport is
! calculated using Bubel

program main
  use modscatter
  use FortUtils
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
  doubleprecision :: T                          ! Transmission
  integer         :: nx = 90                    ! numbers of atoms / 2 in x direction
                                                ! results in about 196 nm
  integer         :: ny = 120                   ! ~numbers of atoms / 2 in y direction (keep even)
                                                ! results in about 170 nm

!!!!!!!!!!!!!!!!!!!!!!!! main function !!!!!!!!!!!!!!!!!!!!!!!
  call parseArguments()
  nx = nx * (16 / sf)
  ny = ny * (16 / sf) + 1

  call createSystem()
  if (save_system) then
    call qt%save_system(trim(results_dir)//"/system.xml")
  endif

  if (run_transport) then
    ! Calculate at specific Fermi energy
    Ef = 0.000001D0 ! eV
    print*,"========================================"
    print*,"Calculating transport at Ef = ",Ef," eV"
    print*,"========================================"
    T =  solveTransport(qt, Ef)
    if (save_densities) then
      call calculateElectronDensity(qt)
      call saveDensities(qt)
    endif

    ! Write to file
    open(unit=101, file=trim(results_dir)//"/single_T.dat")
    write(101,"(A)") "Ef[au],T[-],E0t[au],E0b[au],Vgt[au],Vgb[au],nt[au],nb[au]"
    write(101,"(g0,',',g0,',',g0,',',g0,',',g0,',',g0,',',g0,',',g0)") &
      Ef * eV2au, T, E0t, E0b, Vgt, Vgb, nt, nb

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

  if (plot_results) then
    print*,""
    print*,"========================================"
    print*,"Generating plots..."
    print*,"========================================"
    call generatePlots()
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

! --------------------------------------------------------------------------------------------------
    vecs_armchair = vecs_armchair * sf * geometric_unit2au
    atoms_armchair = atoms_armchair * sf * geometric_unit2au
    pos_offset_armchair = pos_offset_armchair * sf * geometric_unit2au

    ! some magic to have nice edges
    pos_max = atoms_armchair(:, 2) + nx * vecs_armchair(:,1) + ny * vecs_armchair(:,2) - 0.001 ! 0.001 is to ommit numerical errors
    pos_max(1) = pos_max(1) - 2 * (ny / 2) * vecs_armchair(1,2)
    pos_max(2) = pos_max(2) + 2 * pos_offset_armchair(2)
    pos_min = pos_offset_armchair * 0.5 + 0.001 ! 0.001 is to ommit numerical errors
    x_min = atoms_armchair(1,1)
    x_max = atoms_armchair(1,1)
    y_min = atoms_armchair(1,2)
    y_max = atoms_armchair(1,2)

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
    QSYS_DEBUG_LEVEL = 1
    QSYS_FORCE_SCHUR_DECOMPOSITION  = .false. ! don't use schur method so its quicker

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

          if (checkShape(atom_pos, pos_min, pos_max)) then
          ! if (atom_pos(1) > pos_min(1) .and. &
          !     atom_pos(2) > pos_min(2) .and. &
          !     atom_pos(1) < pos_max(1) .and. &
          !     atom_pos(2) < pos_max(2) &
          !     )then

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

    ! yBoundLower = middle_y - (y_max - y_min) * 0.15D0 * 0.5
    ! yBoundUpper = middle_y + (y_max - y_min) * 0.15D0 * 0.5

    yBoundLower = middle_y
    yBoundUpper = middle_y

    write(*, "(A,f10.5,A)"), "middle_x    ", middle_x / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "middle_y    ", middle_y / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "yBoundLower ", yBoundLower / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "yBoundUpper ", yBoundUpper / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "x_min       ", x_min / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "x_max       ", x_max / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "y_min       ", y_min / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "y_max       ", y_max / nm2au, " nm"

    !---------------------------------------- Coupling ---------------------------------------------

    ! Coupling between atoms, onsite energies
    qt%qnnbparam%distance = 0.6 * sf * geometric_unit2au
    qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE

    Bau = Bz * T2au

    call qt%qsystem%make_lattice(qt%qnnbparam, c_simple=connect)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! call addYInvLeads(y_min, y_max, (x_max - x_min) * 1.1 * 0.5, vecs_armchair)

    call addXInvLeads(x_min, x_max, (y_max - y_min) * 1.1, vecs_armchair)
    ! call addYInvLeads(y_min, y_max, (x_max - x_min) * 1.1, vecs_armchair)

    write(*,"(A,f8.5,A)") "nt  ", nt / inv_cmsq2au / 1e11, " 10^11 1/m^2"
    write(*,"(A,f8.5,A)") "nb  ", nb / inv_cmsq2au / 1e11, " 10^11/m^2"
    write(*,"(A,f8.5,A)") "Vgt ", Vgt / eV2au, " eV"
    write(*,"(A,f8.5,A)") "Vgb ", Vgb / eV2au, " eV"
    write(*,"(A,f8.5,A)") "ot  ", (- E0t - Vgt) / eV2au, " eV"
    write(*,"(A,f8.5,A)") "ob  ", (- E0b - Vgb) / eV2au, " eV"
    write(*,"(A,f8.5,A)") "E0t ", E0t / eV2au, " eV"
    write(*,"(A,f8.5,A)") "E0b ", E0b / eV2au, " eV"
    ! print*, " nt ", nt, " 1/nm^2"
    ! print*, " nb ", nb, " 1/nm^2"

  end subroutine
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Add 2 leads on Y sides, Y invariant?
! --------------------------------------------------------------------------------------------------
  subroutine addYInvLeads(y_min, y_max, leadLength, vecs_armchair)
    use modscatter
    use modunits
    use modshape
    implicit none

    type(qshape) :: rect_shape
    doubleprecision :: lead_translation(2) = (/ 0.0D0, 0.0D0 /)!
    doubleprecision, intent(in) :: y_min, y_max, leadLength
    doubleprecision, dimension(2,2), intent(in) :: vecs_armchair

! --------------------------------------------------------------------------------------------------

    lead_translation = (/ 0.0D0, vecs_armchair(2,2) /) * 2
    print*, "lead_translation: ", lead_translation

    ! First lead (lower Y side)
    call rect_shape%init_rect(SHAPE_RECTANGLE_XY, &
                              middle_x - 0.1 - leadLength / 2, &
                              middle_x + 0.1 + leadLength / 2, &
                              y_min - lead_translation(2) * 0.25, &
                              y_min + lead_translation(2) * 0.75)
    call qt%add_lead(rect_shape, (/lead_translation(1), lead_translation(2), 0.0D0 /))

    if (save_bands) then
      call qt%leads(1)%bands(trim(results_dir)//"/bands.dat", &
                            -M_PI / one_over_sqrt_3, +M_PI/ one_over_sqrt_3, M_PI/ one_over_sqrt_3/160.0, & !k_min, k_max, dk
                            -3.0D0 * eV2au, 3.0D0 * eV2au) !E_min, E_max
    endif

    ! Second lead (upper Y side)
    call rect_shape%init_rect(SHAPE_RECTANGLE_XY, &
                              middle_x - 0.1 - leadLength / 2, &
                              middle_x + 0.1 + leadLength / 2, &
                              y_max - lead_translation(2) * 0.75, &
                              y_max + lead_translation(2) * 0.25)

    call qt%add_lead(rect_shape, (/-lead_translation(1), -lead_translation(2), 0.0D0 /))

  end subroutine
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Add 2 leads on X sides, x invariant?
! --------------------------------------------------------------------------------------------------
  subroutine addXInvLeads(x_min, x_max, leadLength, vecs_armchair)
    use modscatter
    use modunits
    use modshape
    implicit none

    type(qshape) :: rect_shape
    doubleprecision :: lead_translation(2) !
    doubleprecision, intent(in) :: x_min, x_max, leadLength
    doubleprecision, dimension(2,2), intent(in) :: vecs_armchair

! --------------------------------------------------------------------------------------------------

    lead_translation = (/vecs_armchair(1,1), 0.0D0/)
    print*, "lead_translation: ", lead_translation
    ! First lead (lower X)
    call rect_shape%init_rect(SHAPE_RECTANGLE_XY, &
                              x_min - 0.1, &
                              x_min + lead_translation(1) - 0.1, &
                              middle_y - 0.1 - leadLength / 2, &
                              middle_y + 0.1 + leadLength / 2)

    call qt%add_lead(rect_shape, (/lead_translation(1), lead_translation(2), 0.0D0 /))

    if (save_bands) then
      call qt%leads(1)%bands(trim(results_dir)//"/bands.dat", &
                            -M_PI / one_over_sqrt_3, M_PI/ one_over_sqrt_3, M_PI/ one_over_sqrt_3/160.0, & !k_min, k_max, dk
                            -3.0D0 * eV2au, 3.0D0 * eV2au) !E_min, E_max
    endif

    ! Second lead (higher X)
    call rect_shape%init_rect(SHAPE_RECTANGLE_XY, &
                              x_max - lead_translation(1) + 0.1, &
                              x_max + 0.1, &
                              middle_y - 0.1 - leadLength / 2, &
                              middle_y + 0.1 + leadLength / 2)

    call qt%add_lead(rect_shape, (/-lead_translation(1), -lead_translation(2), 0.0D0 /))

  end subroutine
! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------
! Check if atom position is within bounds
! --------------------------------------------------------------------------------------------------
  logical function checkShape(atom_pos, pos_min, pos_max)
    implicit none

    doubleprecision, intent(in) :: atom_pos(3)
    doubleprecision, intent(in) :: pos_min(2)
    doubleprecision, intent(in) :: pos_max(2)

    doubleprecision :: range(2)

    range = pos_max - pos_min
    ! contacts width is lover then whole width on Y

    checkShape = &
      (atom_pos(2) > pos_min(2) + range(2) * 0.1 .and. &
       atom_pos(2) < pos_max(2) - range(2) * 0.1 .and. &
       atom_pos(1) > pos_min(1) .and. &
       atom_pos(1) < pos_max(1)) &
      .or. &
      (atom_pos(1) > pos_min(1) + range(1) * 0.05 .and. &
       atom_pos(1) < pos_max(1) - range(1) * 0.05 .and. &
       atom_pos(2) > pos_min(1) .and. &
       atom_pos(2) < pos_max(1))

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
    doubleprecision :: y
! --------------------------------------------------------------------------------------------------
    if (.not. (atomA%flag == atomB%flag)) then
      connect = .true.
      t0 = (3.0D0 * eV2au) / sf
      xA = atomA%atom_pos(1)
      yA = atomA%atom_pos(2)
      xB = atomB%atom_pos(1)
      yB = atomB%atom_pos(2)
      B = Bau
      y = (yB + yA) * 0.5
      if (y < middle_y) B = -Bau ! bottom

      ! Peierls phase
      phi = 0.5 * B * (yB + yA) * (xB - xA) ! y x already in au
      coupling_val = t0 * exp(II*phi)
    else
      connect = .true.
      xA = atomA%atom_pos(1)
      yA = atomA%atom_pos(2)
      xB = atomB%atom_pos(1)
      yB = atomB%atom_pos(2)
      y = (yB + yA) * 0.5
      Vg = linear(y, Vgt, Vgb, yBoundUpper, yBoundLower)
      E0 = linear(y, E0t, E0b, yBoundUpper, yBoundLower)
      coupling_val = - E0 - Vg
    endif
  end function
! --------------------------------------------------------------------------------------------------
end program main
