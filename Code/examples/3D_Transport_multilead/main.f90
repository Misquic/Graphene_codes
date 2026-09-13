! This program calculates transport for twisted graphene flake in 3D.
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
  doubleprecision :: z_max      = 0.0D0
  doubleprecision :: foldCenter(3)

  doubleprecision :: Bau                        ! in au
  doubleprecision :: Vgt, Vgb, E0t, E0b, nt, nb ! result from Bilayer
  doubleprecision :: Ef                         ! Fermi energy for calculations
  integer, parameter :: numLeads = 4
  doubleprecision, dimension(numLeads, numLeads) :: T ! Transmission matrix
  integer         :: nx = 200 / 16             ! numbers of atoms / 2 in x direction
  ! integer         :: nx = 90                  ! numbers of atoms / 2 in x direction
  ! results in about 196 nm
  ! integer         :: nx = 180                 ! numbers of atoms / 2 in x direction
  ! integer         :: ny = 122                 ! ~numbers of atoms / 2 in y direction (keep even)
  integer         :: ny = 200 / 16             ! ~numbers of atoms / 2 in y direction (keep even)
  ! integer         :: ny = 60                  ! ~numbers of atoms / 2 in y direction (keep even)
  ! integer         :: ny = 244                 ! ~numbers of atoms / 2 in y direction (keep even)
                                                ! results in about 170 nm
  integer         :: from = 1, to = 1           ! lead indexes for writing transmissions to a file
  integer :: n
  doubleprecision :: r = 0
  doubleprecision :: AA = 0.075D0               ! Anderson potential
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
  write(*, "(A,i,A,i)"), "Nx: ", nx, " Ny: ", ny
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

    if (plot_results) then
      call generatePlots()
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

    ! doubleprecision :: foldRadius = 0.5D0/2
    doubleprecision :: foldRadius = 5.0D0/2
    doubleprecision :: position_scale
    integer         :: i, atom
    integer(c_int)  :: position_count
    type(c_ptr)     :: position_data
    real(c_double), pointer :: positions(:)
    doubleprecision :: atom_pos(3)
    doubleprecision :: x_min      = 0.0D0, x_max = 0.0D0
    doubleprecision :: y_min      = 0.0D0, y_max = 0.0D0
    type(c_ptr)     :: bilayer
    integer(c_int)  :: cutLead = 15, leadWidth = 25
    doubleprecision           :: vecs_armchair(2,2)     = (/ (/ 1.0D0,0.0D0 /), (/ sin(alpha30), cos(alpha30) /) /)
    integer         :: currentLead = 1
    print*, "Creating System"
    ! foldRadius = scaleBasedOnSf(foldRadius / geometric_unit, sf)! in nm. Positions are later scaled by sf
    ! foldRadius = foldRadius! in nm. Positions are later scaled by sf
    ! cutLead = scaleBasedOnSfInt(cutLead, sf)
    cutLead = ny * 0.35
    ! leadWidth = scaleBasedOnSfInt(leadWidth, sf)
    leadWidth = ny * 0.55

    write(*, "(A,f5.3,A,i,A,i)") "R=", foldRadius, " cutLead=", cutLead, " leadWidth=", leadWidth

    vecs_armchair = vecs_armchair * sf * geometric_unit2au

! --------------------------------------------------------------------------------------------------

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

    position_scale = sf * geometric_unit2au
    call Bilayer_generatePositions3D(nx, ny, foldRadius  / geometric_unit / sf, &
                                     cutLead, leadWidth, &
                                     position_data, position_count)
    call c_f_pointer(position_data, positions, [3 * position_count])

    do i = 1, position_count
      atom_pos = position_scale * positions((3 * i - 2) : (3 * i))
      x_max = max(x_max, atom_pos(1))
      x_min = min(x_min, atom_pos(1))
      y_max = max(y_max, atom_pos(2))
      y_min = min(y_min, atom_pos(2))
      z_max = max(z_max, atom_pos(3))
      atom = 1 + merge(0, 1, i <= position_count / 2)
      call qa%init(atom_pos, flag=atom)
      call qt%qsystem%add_atom(qa)
    enddo

    foldCenter = (/ 0.0D0, y_max - foldRadius * nm2au, foldRadius * nm2au/)
#ifdef DEBUG
    print*, "Saving additional Atom at the centre of a fold"
    call qa%init(foldCenter, flag=atom)
    call qt%qsystem%add_atom(qa)
#endif

    middle_x = 0.5 * (x_min + x_max)
    middle_y = 0.5 * (y_min + y_max)

    write(*, "(A,f10.5,A)"), "middle_x    ", middle_x / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "middle_y    ", middle_y / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "x_min       ", x_min / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "x_max       ", x_max / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "y_min       ", y_min / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "y_max       ", y_max / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "z_max       ", z_max / nm2au, " nm"
    write(*, "(A,f10.5,A)"), "AA          ", AA, " eV"

    !----------------------------------------- Leads -----------------------------------------------

    ! I
    ! lead 1 and 2 are added (voltage)
    call leads_init()
    call addXInvFoldLeads(currentLead, x_min, x_max, y_max, foldRadius * nm2au, vecs_armchair, cutLead)
    call addXInvUpDownLeads(currentLead, x_min, x_max, y_min, vecs_armchair, leadWidth)
    !---------------------------------------- Coupling ---------------------------------------------

    ! Coupling between atoms, onsite energies
    qt%qnnbparam%distance = 0.6 * sf * geometric_unit2au
    qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE

    Bau = Bz * T2au

    call qt%qsystem%make_lattice(qt%qnnbparam, c_simple=connect)
    ! ! Remove single bonds if necessary, loop through all atoms and
    ! ! check if the number of bonds is 2. (1st - for on site bond, 2nd for neightbour)
    ! do atom = 1 ,qt%qsystem%no_atoms
    !     if(qt%qsystem%atoms(atom)%no_bonds == 2)   qt%qsystem%atoms(atom)%bActive  = .false.
    ! enddo
    ! call qt%qsystem%make_lattice(qt%qnnbparam,c_simple=connect)

    !------------------------------------------ Leads pt 2 -----------------------------------------
    do i = 1, numLeads
      call qt%add_lead(leadShapes(i), leadTrans(:, i))
      if (save_bands) then
        call qt%leads(i)%bands(trim(results_dir)//"/bands" // trim(str(i)) // ".dat", &
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
    do i = 1, numLeads
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

    write(*, "(A,f10.2,A,f10.2,A,f10.2,A,f10.2)"), "Xmin=", xmin/nm2au, " Xmax=", xmax/nm2au, " Ymin=", ymin/nm2au, " Ymax=", ymax/nm2au
    call leadShapes(idx)%init_rect(SHAPE_RECTANGLE_XY, xmin, xmax, ymin, ymax)
    leadTrans(:,idx) = translation
  end subroutine addLeadRect


  integer function isInLeads(x, y, z)
    implicit none
    doubleprecision, intent(in) :: x, y, z
    doubleprecision :: vec(3)
    integer :: i
    isInLeads = 0
    vec = (/ x, y, z /)
    do i = 1, numLeads
      if (leadShapes(i)%is_inside(vec)) then
        isInLeads = i
        return
      end if
    end do
  end function isInLeads


! --------------------------------------------------------------------------------------------------
! Add 2 leads on X sides on fold?
! --------------------------------------------------------------------------------------------------
  subroutine addXInvFoldLeads(idx, x_min, x_max, y_max, foldRadius, vecs_armchair, cutLead)
    use modscatter
    use modunits
    use modshape
    implicit none

    integer :: idx
    doubleprecision :: lead_translation(2) !
    doubleprecision, intent(in) :: x_min, x_max
    doubleprecision, intent(in) :: y_max
    doubleprecision, intent(in) :: foldRadius
    integer, intent(in) :: cutLead
    doubleprecision, dimension(2,2), intent(in) :: vecs_armchair

! --------------------------------------------------------------------------------------------------

    lead_translation = (/vecs_armchair(1,1), 0.0D0/)
    print*, "lead_translation: ", lead_translation
    ! First lead (lower X)
    call addLeadRect(idx, &
                     x_min - 0.1, &
                     x_min + lead_translation(1) - 0.099, &
                     y_max - 0.1 - (cutLead + 2) * vecs_armchair(2,2) - foldRadius, &
                     y_max + 0.1, &
                     (/lead_translation(1), lead_translation(2), 0.0D0/))
    idx = idx + 1

    ! Second lead (higher X)
    call addLeadRect(idx, &
                     x_max - lead_translation(1) + 0.099, &
                     x_max + 0.1, &
                     y_max - 0.1 - (cutLead + 2) * vecs_armchair(2,2) - foldRadius, &
                     y_max + 0.1, &
                     (/-lead_translation(1), -lead_translation(2), 0.0D0/))
    idx = idx + 1

  end subroutine
! --------------------------------------------------------------------------------------------------


  ! add leads that span two layers
  subroutine addXInvUpDownLeads(idx, x_min, x_max, y_min, vecs_armchair, leadWidth)
    use modscatter
    use modunits
    use modshape
    implicit none

    integer :: idx
    doubleprecision :: lead_translation(2) !
    doubleprecision, intent(in) :: x_min, x_max
    doubleprecision, intent(in) :: y_min
    integer, intent(in) :: leadWidth
    doubleprecision, dimension(2,2), intent(in) :: vecs_armchair

    lead_translation = (/vecs_armchair(1,1), 0.0D0/)
    print*, "lead_translation: ", lead_translation
    ! First lead (lower X)
    call addLeadRect(idx, &
                     x_min - 0.1, &
                     x_min + lead_translation(1) - 0.099, &
                     y_min - 0.1, &
                     y_min + 0.1 + (leadWidth + 3) * vecs_armchair(2,2), &
                     (/lead_translation(1), lead_translation(2), 0.0D0/))
    idx = idx + 1

    ! Second lead (higher X)
    call addLeadRect(idx, &
                     x_max - lead_translation(1) + 0.099, &
                     x_max + 0.1, &
                     y_min - 0.1, &
                     y_min + 0.1 + (leadWidth + 3) * vecs_armchair(2,2), &
                     (/-lead_translation(1), -lead_translation(2), 0.0D0/))
    idx = idx + 1
  end subroutine
! --------------------------------------------------------------------------------------------------


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
    doubleprecision :: t
    doubleprecision, parameter :: beta = 3.37D0, t0 = 3.0D0
    doubleprecision :: Vg
    doubleprecision :: E0
    doubleprecision :: x, y, z
    doubleprecision :: r
    doubleprecision :: dist

! --------------------------------------------------------------------------------------------------
    if (.false.) atoms(0)%flag = atoms(0)%flag ! supress unused variable warning
    if (.not. (atomA%flag == atomB%flag)) then
      ! hopping
      connect = .true.
      ! t = ( - t0 * eV2au ) / sf
      t = ( slaterCoster(atomA%atom_pos, atomB%atom_pos, foldCenter) * eV2au ) / sf
      xA = atomA%atom_pos(1)
      yA = atomA%atom_pos(2)
      xB = atomB%atom_pos(1)
      yB = atomB%atom_pos(2)

      B = Bau
      dist = len(atomA%atom_pos - atomB%atom_pos)
      t = t * exp(- beta * ((dist / (carbonCarbonDistAu * sf)) - 1))
      ! write(*, "(A,f10.5)"), "exp=", exp(- beta * ((dist / (carbonCarbonDistAu * sf)) - 1))

      ! Peierls phase
      phi = - 0.5 * B * (yB + yA - 2 * middle_y) * (xB - xA) ! y x already in au
      coupling_val = t * exp(II*phi)
    else
      ! onsite
      connect = .true.
      x = atomA%atom_pos(1)
      y = atomA%atom_pos(2)
      z = atomA%atom_pos(3)
      Vg = linear(z, Vgt, Vgb, z_max, 0.0D0)
      E0 = linear(z, E0t, E0b, z_max, 0.0D0)
      coupling_val = - E0 - Vg
      if (isInLeads(x, y, z) == 0) then
        ! not lead
        call random_number(r)
        ! potencjał andersona, +-5/100 eV dla sf8 -> 1/10 2/10 -----> sf4 - pomnożyć przez 2 sf ||||| sf rośnie AA maleje
        ! add random number +- 0.05 eV
        coupling_val = coupling_val + (2 * r - 1) * AA * eV2au
      endif
    endif
  end function
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Calculate hoping between atoms for bent connections
! --------------------------------------------------------------------------------------------------
  doubleprecision function slaterCoster(posA, posB, foldCenter) result(t)
    implicit none
    ! args
    doubleprecision :: posA(3), posB(3), foldCenter(3)

    ! local
    doubleprecision, parameter :: Vpp_pi    = - 3.0D0 ! standard hopping
    doubleprecision, parameter :: Vpp_sigma = 1.7 * Vpp_pi !

    doubleprecision :: d(3), nA(3), nB(3)
! --------------------------------------------------------------------------------------------------

    ! vector from first atom to second
    d = posB - posA

    ! versors perpendicular to positions
    call getVersor(posA, foldCenter, nA)
    call getVersor(posB, foldCenter, nB)

    t = Vpp_pi * dot(nA, nB) + &
        (Vpp_sigma - Vpp_pi) * dot(nA, d) * dot(nB, d)

  end function
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Get versor perpendicular to surface at the position pos
! --------------------------------------------------------------------------------------------------
  subroutine getVersor(pos, foldCenter, versor)
    implicit none
    doubleprecision :: pos(3), foldCenter(3), versor(3)
! --------------------------------------------------------------------------------------------------
    ! if not fold
    if (pos(2) < foldCenter(2)) then
      if (pos(3) < foldCenter(3)) then
        versor = (/ 0.0D0, 0.0D0, 1.0D0 /)
      else
        versor = (/ 0.0D0, 0.0D0, -1.0D0 /)
      endif
    else ! if fold
      versor = foldCenter - pos
      versor(1) = 0.0D0 ! get vector to axis of cyllinder
      versor = versor / len(versor)
    endif
  end subroutine
! --------------------------------------------------------------------------------------------------

end program main

! TODO get "andersonAnmplitude" -> maks wartość zaburzenia to be a parameter
! TODO set andersonAmplitude in leads to 0
