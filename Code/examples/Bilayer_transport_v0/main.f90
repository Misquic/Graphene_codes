! This program calculates transport for "simply" twisted graphene flake. It means
! that twisted boundary is not taken into account, there is just flip of magnetic 
! field. System looks like this (x and y used are for 1 view reference, not physical coordinates)
!    (unfolded view)            !   (side view left)             (side view right)        (side view front)                            
! y   ________________________  !   ____lead1___________        ____lead2___________     ___                   ___
! ^  |l|  top              |l|  !   |  ____top_gate___ |        |  ____top_gate___ |     |l|  ___top_gate____  |l|     
! |  |e|  B = (0,0,Bz)     |e|  !   |  _______top_____ |     ^  | _______top_____  |     |e|_______top_________|e| 
! |  |a|___________________|a|  !   | /                | B = |  |                \ |     |a|                   |a|
! |  |d|  bottom           |d|  !   | \_____bottom____ |     |  | _____bottom____/ |     |d|_____bottom________|d|
! |  |1|  B = (0,0,-Bz)    |2|  !   |  _bottom_gate___ |        |  _bottom_gate___ |     |1|  __bottom_gate__  |2|
! |  |_|___________________|_|  !   |__________________|        |__________________|     |_|                   |_|
! *--------------> x            !                                                                                        
!  
! Coupling between layers is taken into account with Bilayer class, transport is
! calculated using Bubel

program main
  use Bilayer_interface
  use modscatter
  use modsys
  use modshape
  use modunits

  implicit none

! Bubel objects
  type(qscatter) :: qt
  type(qatom) :: qa
  type(qshape) :: rect_shape

! graphene parameters
  doubleprecision,parameter :: alpha30 =  30.0/180.0*M_PI

  doubleprecision,parameter :: vecs_armchair(2,2) =  (/  (/ 1.0D0,0.0D0 /) , (/ sin(alpha30) , cos(alpha30) /) /)
  doubleprecision,parameter :: atoms_armchair(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 0.0D0 , 1.0D0/sqrt(3.0) /) /)
  doubleprecision,parameter :: pos_offset(2) =  (/ -5.0D0,0.0D0 /)
 
  ! doubleprecision,parameter :: vecs_zigzag(2,2) =   (/  (/ (3.0/2.0)/sqrt(3.0D0)  ,0.5D0 /) , (/ -(3.0/2.0)/sqrt(3.0D0) , 0.5D0 /) /)
  ! doubleprecision,parameter :: atoms_zigzag(2,2) =  (/  (/ 0.0D0,0.0D0 /) , (/ 1.0D0/sqrt(3.0)  ,0.0D0  /) /)
  ! doubleprecision,parameter :: pos_offset_zigzag(2) =  (/ -5.0D0,-10.0D0 /)

! local variables
  integer,parameter :: nx = 300, ny = 200     ! numbers of atoms in x and y direction (system coordinates)
  integer,parameter :: atomA = 1, atomB = 2 ! sublattices flags

  integer ::         i, j, atom  ! loop variables
  doubleprecision :: atom_pos(3) !
  doubleprecision :: half_y      ! Y coordinate of fold
  doubleprecision :: Bz          ! B = (0, 0, Bz)

!!!!!!!!!!!!!!!!!!!!!!!! main function !!!!!!!!!!!!!!!!!!!!!!!
  call qt%init_system()
  !QSYS_FORCE_SCHUR_DECOMPOSITION  = .true. ! use schur method to calculate modes which is more stable
  
  ! Generate atoms positions
  do i = 0, nx-1
    do j = 0, ny-1
      do atom = atomA, atomB
        atom_pos(1:2) = atoms_armchair(:,atom) + & ! base position, atom choses sublattice
          i * vecs_armchair(:,1) + & ! offset in x direction for i-th "column"
          j * vecs_armchair(:,2) + & ! offset in y dir for j-th "row"
          pos_offset ! base offset
        call qa%init( (/atom_pos(1), atom_pos(2), 0.0D0 /), flag=atom) ! create qatom
        call qt%qsystem%add_atom(qa) ! add qatom to system
      enddo
    enddo
  enddo
  half_y = atoms_armchair(1,2) + (ny-1)*0.5*vecs_armchair(1,2) + pos_offset(2)

  ! Coupling between atoms, onsite energies
  qt%qnnbparam%distance = 0.6
  qt%qnnbparam%NNB_FILTER = QSYS_NNB_FILTER_DISTANCE

  Bz = 0.03 ! atomic units?
  call qt%qsystem%make_lattice(qt%qnnbparam, c_simple=connect)

  call qt%save_system("results/system.xml")

contains
! --------------------------------------------------------
! Calculate hoping between atoms, here we use Peierls phase
! to simulate magnetic field with gauge: A = (-Bz * y,0,0)
! i.e. B = (0,0,Bz)
! --------------------------------------------------------
  logical function connect(atomA, atomB, coupling_val, atoms)
    use modcommons
    
    implicit none

    type(qatom) :: atomA, atomB ! 
    type(qatom) :: atoms(:)     ! 
    complex*16  :: coupling_val ! 

    doubleprecision :: xA, yA, xB, yB
    doubleprecision :: phi ! Peirles phase
    doubleprecision :: B

    connect = .not. (atomA%flag == atomB%flag)
    if (connect) then
      xA = atomA%atom_pos(1)
      yA = atomA%atom_pos(2)
      xB = atomB%atom_pos(1)
      yB = atomB%atom_pos(2)
      B = Bz
      if ((yB+yA)*0.5 < half_y) B = -Bz 
      phi = 0.5*B*(yB+yA)*(xB-xA)
      coupling_val = 1.0D0 * exp(II*phi)
    endif

  end function
end program main