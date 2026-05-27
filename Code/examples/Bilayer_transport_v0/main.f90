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
  use, intrinsic :: iso_c_binding
  implicit none

  ! Bubel objects
  type(qscatter) :: qt

  ! Shared variables used by internal procedures
  doubleprecision :: middle_y             ! Y coordinate of fold
  doubleprecision :: middle_x             !
  doubleprecision :: yBoundLower          ! Y coordinate where Vg goes from Vgb to linear region
  doubleprecision :: yBoundUpper          ! Y coordinate where Vg goes from linear region to Vgt
  character(len=512) :: results_dir = "./results" ! Directory for output files

  doubleprecision :: Bz = 8                     ! B = (0, 0, Bz) !T
  doubleprecision :: Bau                        ! in au
  doubleprecision :: Vt = 0, Vb = 5             ! eV
  doubleprecision :: Vgt, Vgb, E0t, E0b, nt, nb ! result from Bilayer
  doubleprecision :: Ef                         ! Fermi energy for calculations
  integer         :: sf = 8                     ! scaling factor
  integer         :: nx = 25                    ! numbers of atoms / 2 in x direction
                                                ! results in about 196 nm
  integer         :: ny = 60                    ! ~numbers of atoms / 2 in y direction (keep even)
                                                ! results in about 408 nm

  doubleprecision,parameter :: T2au        = 4.254382E-6          ! B(au) = B(T)*T2au
  doubleprecision,parameter :: eV2au       = 0.03674932587122423  ! V(au)  = V(eV)*eV2au
  doubleprecision,parameter :: nm2au       = 1.0 / 0.0529           ! d(au)  = d(nm)*nm2au
  doubleprecision,parameter :: cm2au       = 1e-2 * 1e9 * nm2au
  doubleprecision,parameter :: inv_cmsq2au = 1. / cm2au / cm2au
  doubleprecision,parameter :: one_over_sqrt_3 = 1.0D0 / sqrt(3.0)

  logical :: run_transport = .true.
  logical :: run_energyScan = .true.
  logical :: plot_results = .true.
  logical :: save_system = .true.
  logical :: save_densities = .true.
  logical :: save_bands = .true.

  integer :: currentArgIndex = 1

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
    Ef = 0.0001D0 ! eV
    print*,"========================================"
    print*,"Calculating transport at Ef = ",Ef," eV"
    print*,"========================================"
    call solveTransport(qt, Ef)
    if (save_densities) then
      call calculateElectronDensity(qt)
      call saveResults(qt)
    endif

    ! Perform energy scan
    if (run_energyScan) then
      print*,""
      print*,"========================================"
      print*,"Performing energy scan..."
      print*,"========================================"
      call performEnergyScan(qt)
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
! Get next command line argument
! --------------------------------------------------------------------------------------------------
  logical function getNextArgument(arg_buffer)
    implicit none
    character(len=512) :: arg_buffer
    integer :: argc

! --------------------------------------------------------------------------------------------------

    argc = command_argument_count()
    getNextArgument = .false.
    if (argc >= currentArgIndex) then
      call get_command_argument(currentArgIndex, arg_buffer)
      currentArgIndex = currentArgIndex + 1
      getNextArgument = .true.
    endif
  end function getNextArgument
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Parse bool from next command line argument
! --------------------------------------------------------------------------------------------------
  logical function parseBoolArg(defaultValue)
    implicit none
    logical defaultValue
    character(len=512) :: arg_buffer

! --------------------------------------------------------------------------------------------------
    if (getNextArgument(arg_buffer)) then
      ! Default to true, set to false if starts with 'f', 'F', or '0'
      if (arg_buffer(1:1) == 't' .or. arg_buffer(1:1) == 'T' .or. arg_buffer(1:1) == '1') then
        parseBoolArg = .true.
      else
        parseBoolArg = .false.
      endif
    else
      parseBoolArg = defaultValue
    endif
  end function parseBoolArg
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Parse integer from next command line argument
! --------------------------------------------------------------------------------------------------
  integer function parseIntArg(defaultValue)
    implicit none
    integer :: defaultValue
    character(len=512) :: arg_buffer

! --------------------------------------------------------------------------------------------------
    if (getNextArgument(arg_buffer)) then
      read(arg_buffer, *) parseIntArg
    else
      parseIntArg = defaultValue
    endif
  end function parseIntArg
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Parse Double Precision from next command line argument
! --------------------------------------------------------------------------------------------------
  doubleprecision function parseDoubleArg(defaultValue)
    implicit none
    double precision :: defaultValue
    character(len=512) :: arg_buffer

! --------------------------------------------------------------------------------------------------
    if (getNextArgument(arg_buffer)) then
      read(arg_buffer, *) parseDoubleArg
    else
      parseDoubleArg = defaultValue
    endif
  end function parseDoubleArg
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Parse Arguments into global variables
! --------------------------------------------------------------------------------------------------
  subroutine parseArguments()
    character(len=512) :: arg_buffer
    character(len=512) :: help_buffer

! --------------------------------------------------------------------------------------------------

    ! check first argument if it's "help" then print help and exit, else its results dir
    if (getNextArgument(arg_buffer)) then
      help_buffer = trim(arg_buffer)
      if (help_buffer == "help") then
        print*, "usage: ./Transport2D <resultsDir> <B in T> <Vb> <Vt> &
                 <save_system> <run_transport> <run_energyScan> &
                 <plot_results> <save_densities> <save_bands> <sf>"
        call exit(0)
      else
        results_dir = trim(arg_buffer)
      endif
    endif

    Bz = parseDoubleArg(Bz)
    Vb = parseDoubleArg(Vb)
    Vt = parseDoubleArg(Vt)
    save_system = parseBoolArg(.true.)
    run_transport = parseBoolArg(.false.)
    run_energyScan = parseBoolArg(run_energyScan)
    plot_results = parseBoolArg(.false.)
    save_densities = parseBoolArg(.false.)
    save_bands = parseBoolArg(.false.)
    sf = parseIntArg(sf)

    print*, "usage: ./Transport2D <resultsDir> <B in T> <Vb> <Vt> &
             <save_system> <run_transport> <run_energyScan> &
             <plot_results> <save_densities> <save_bands> <sf>"
    print*, ""
    print*, "Parsed Arguments"
    print*, ""
    print*, "results_dir: ", trim(results_dir)
    print*, "Bz: ", Bz, " T"
    print*, "Vb: ", Vb, " eV"
    print*, "Vt: ", Vt, " eV"
    print*, "save_system: ", save_system
    print*, "run_transport: ", run_transport
    print*, "run_energyScan: ", run_energyScan
    print*, "plot_results: ", plot_results
    print*, "save_densities: ", save_densities
    print*, "save_bands: ", save_bands
    print*, "sf: ", sf
    print*, ""
    print*, ""
    print*, ""
  end subroutine
! --------------------------------------------------------------------------------------------------



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

    call qt%init_system()
    QSYS_DEBUG_LEVEL = 0
    QSYS_FORCE_SCHUR_DECOMPOSITION  = .true. ! use schur method to calculate modes which is more stable
    ! QSYS_SCATTERING_METHOD = QSYS_SCATTERING_QTBM

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

          if (atom_pos(1) > pos_min(1) .and. &
              atom_pos(2) > pos_min(2) .and. &
              atom_pos(1) < pos_max(1) .and. &
              atom_pos(2) < pos_max(2) &
              )then
              ! atom_pos(1) < x_max + pos_offset_armchair(1) / 2) then

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

    yBoundLower = middle_y - (y_max - y_min) * 0.15D0 * 0.5
    yBoundUpper = middle_y + (y_max - y_min) * 0.15D0 * 0.5

    print*, "middle_x", middle_x / nm2au, " nm"
    print*, "middle_y", middle_y / nm2au, " nm"
    print*, "x_min", x_min / nm2au, " nm"
    print*, "x_max", x_max / nm2au, " nm"
    print*, "y_min", y_min / nm2au, " nm"
    print*, "y_max", y_max / nm2au, " nm"

    ! Coupling between atoms, onsite energies
    qt%qnnbparam%distance = 0.6 * sf * geometric_unit2au
    qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE

    Bau = Bz * T2au

    bilayer = Bilayer_constructor_default()

    call Bilayer_countAll_B(bilayer, Vt, Vb, Bz, Vgt, Vgb, E0t, E0b, nt, nb)

    print*, "Vgt: ", Vgt, " eV"
    print*, " Vgb ", Vgb, " eV"
    print*, "E0t ", E0t, " eV"
    print*, " E0b ", E0b, " eV"
    print*, " nt ", nt, " 1/nm^2"
    print*, " nb ", nb, " 1/nm^2"
    ! convert back to au
    Vgt = Vgt * eV2au
    Vgb = Vgb * eV2au
    E0t = E0t * eV2au
    E0b = E0b * eV2au
    nt = nt * inv_cmsq2au
    nb = nb * inv_cmsq2au

    call qt%qsystem%make_lattice(qt%qnnbparam, c_simple=connect)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! call addYInvLeads(y_min, y_max, (x_max - x_min) * 1.1 * 0.5, vecs_armchair)

    call addXInvLeads(x_min, x_max, (y_max - y_min) * 1.1, vecs_armchair)
    ! call addYInvLeads(y_min, y_max, (x_max - x_min) * 1.1, vecs_armchair)

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
! Calculate linear gradiend between upper and down side
! --------------------------------------------------------------------------------------------------

  doubleprecision function linear(y, bottomValue, topValue)
    implicit none

    doubleprecision, intent(in) :: y
    doubleprecision, intent(in) :: bottomValue
    doubleprecision, intent(in) :: topValue
    doubleprecision :: yRange, dy, VRange

    if (y < yBoundLower) then
      linear = bottomValue
    else if (y > yBoundUpper) then
      linear = topValue
    else
      yRange = yBoundUpper - yBoundLower
      dy = y - yBoundLower
      VRange = topValue - bottomValue
      linear = dy / yRange * VRange + bottomValue
    endif

  end function



! --------------------------------------------------------------------------------------------------
! Calculate hoping between atoms, here we use Peierls phase
! to simulate magnetic field with gauge: A = (-Bz * y,0,0)
! i.e. B = (0,0,Bz)
! --------------------------------------------------------------------------------------------------
  logical function connect(atomA, atomB, coupling_val, atoms)
    use modcommons
    implicit none

    type(qatom) :: atomA, atomB !
    type(qatom) :: atoms(:)     !
    complex*16  :: coupling_val !

    doubleprecision :: xA, yA, xB, yB
    doubleprecision :: phi ! Peirles phase
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
      Vg = linear(y, Vgb, Vgt)
      E0 = linear(y, E0b, E0t)
      coupling_val = + E0 - Vg
    endif
  end function
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Solve transport problem
! --------------------------------------------------------------------------------------------------
  subroutine solveTransport(qt, Ef)
    use modscatter
    implicit none
    type(qscatter) :: qt
    doubleprecision :: Ef
    doubleprecision :: T_total
    integer, parameter :: leadsIds(1) = (/ 1 /)

! --------------------------------------------------------------------------------------------------
    print*,"  Solving transport..."
    call qt%calculate_modes(Ef * eV2au)
    call qt%solve(1, Ef * eV2au) ! TEST if it is faster
    ! call qt%solve_leads(leadsIds, Ef * eV2au) ! TEST if it is faster
    T_total = sum(qt%Tn(:))
    print*,"  Total transmission: ", T_total

    ! Write to file
    open(unit=101, file=trim(results_dir)//"/single_T.dat")
    write(101,"(A)") "Ef(au),T_total,E0t(au),E0b(au),Vgt(au),Vgb(au),nt(au),nb(au)"
    write(101,"(g0,',',g0,',',g0,',',g0,',',g0,',',g0,',',g0,',',g0)") &
      Ef * eV2au, T_total, E0t, E0b, Vgt, Vgb, nt, nb

    close(101)
  end subroutine solveTransport
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Calculate electron density
! --------------------------------------------------------------------------------------------------
  subroutine calculateElectronDensity(qt)
    use modscatter
    implicit none
    type(qscatter) :: qt
    integer :: i

! --------------------------------------------------------------------------------------------------
    print*,"  Calculating electron density..."
    do i = 1, size(qt%qsystem%qauxvec)
      qt%qsystem%qauxvec(i) = sum(qt%qsystem%densities(:,i))
    enddo
  end subroutine calculateElectronDensity
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Save calculated results
! --------------------------------------------------------------------------------------------------
  subroutine saveResults(qt)
    use modscatter
    implicit none
    type(qscatter) :: qt

! --------------------------------------------------------------------------------------------------
    print*,"  Saving results..."
    call qt%qsystem%save_data(trim(results_dir)//"/densities.xml", &
                              array2d=qt%qsystem%densities, &
                              array1d=qt%qsystem%qauxvec)
  end subroutine saveResults
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Perform energy scan
! --------------------------------------------------------------------------------------------------
  subroutine performEnergyScan(qt)
    use modscatter
    implicit none
    type(qscatter) :: qt
    double precision :: E_scan, T_total
    double precision, parameter :: deltaE = 0.001D0
! --------------------------------------------------------------------------------------------------
    ! Open output file for energy scan
    open(unit=100, file=trim(results_dir)//"/T.dat")

    ! Energy scan parameters
    E_scan = (-0.1D0 + 0.0001D0) * eV2au
    do while (E_scan <= 0.1D0 * eV2au)
      ! Update hamiltonian elements
      call qt%qsystem%update_lattice(c_simple=connect)

      ! Calculate modes and solve
      call qt%calculate_modes(E_scan)
      call qt%solve(1, E_scan)

      ! Get total transmission
      T_total = sum(qt%Tn(:))

      ! Write to file
      write(100,"(g0,',',g0)") E_scan, T_total

      E_scan = E_scan + deltaE * eV2au
    enddo

    close(100)
  end subroutine performEnergyScan
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Generate plots using existing Python scripts
! --------------------------------------------------------------------------------------------------
  subroutine generatePlots()
    implicit none

! --------------------------------------------------------------------------------------------------
    if (save_bands) then
      print*,"  Plotting band structure..."
      call execute_command_line("python plot_bands.py "//trim(results_dir)//"/")
    endif
    print*,"  Plotting Transmission..."
    call execute_command_line("python plot_T.py "//trim(results_dir)//"/")
  end subroutine generatePlots
! --------------------------------------------------------------------------------------------------

end program main

! TODO
! - skalowanie np 8 (czy przy skalowaniu wystarczy mniejszy np 48x96)
! - powiększyć układ
! - poszerzyć leady
! - G = 2e^2/h*T

! cyfronet:
!   - module load intel/2023b lub intel/2025b?

! 3D układ (całki przeskoku?)
! więcej leadów (4pkt?)

! Zapytać jeszcze o
! potencjał/relaksacja - najpierw gładkie przejście może wystarczy
! r

! Dlaczego na wykresach V jest skok względem B?