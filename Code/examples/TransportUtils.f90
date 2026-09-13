module TransportUtils
  use FortUtils
  implicit none

  ! unit conversion
  doubleprecision,parameter :: T2au            = 4.254382E-6          ! B(au) = B(T)*T2au
  doubleprecision,parameter :: eV2au           = 0.03674932587122423  ! V(au)  = V(eV)*eV2au
  doubleprecision,parameter :: nm2au           = 1.0 / 0.0529           ! d(au)  = d(nm)*nm2au
  doubleprecision,parameter :: cm2au           = 1e-2 * 1e9 * nm2au
  doubleprecision,parameter :: inv_cmsq2au     = 1. / cm2au / cm2au
  doubleprecision,parameter :: one_over_sqrt_3 = 1.0D0 / sqrt(3.0)
  doubleprecision,parameter :: M_PI            = 3.1415926535898

  ! other
  doubleprecision, parameter :: alpha30            = 30.0/180.0*M_PI
  doubleprecision, parameter :: carbonCarbonDist   = 0.142D0
  doubleprecision, parameter :: carbonCarbonDistAu = carbonCarbonDist * nm2au
  doubleprecision, parameter :: geometric_unit     = carbonCarbonDist * sqrt(3.0D0)
  doubleprecision, parameter :: geometric_unit2au  = geometric_unit * nm2au

abstract interface
  logical function func_simple(atomA,atomB,coupling_val,atoms)
      use modcommons
      implicit none
      type(qatom) :: atomA,atomB
      type(qatom) :: atoms(:)
      complex*16  :: coupling_val
  end function func_simple
end interface

contains


! --------------------------------------------------------------------------------------------------
! Calculate linear gradient between upper and down side
! --------------------------------------------------------------------------------------------------
  doubleprecision function linear(y, topValue, bottomValue, yBoundUpper, yBoundLower)
    implicit none

    doubleprecision, intent(in) :: y
    doubleprecision, intent(in) :: topValue
    doubleprecision, intent(in) :: bottomValue
    doubleprecision, intent(in) :: yBoundUpper
    doubleprecision, intent(in) :: yBoundLower
    doubleprecision :: yRange, dy, VRange

! --------------------------------------------------------------------------------------------------
    if (yBoundLower > yBoundUpper) then
      stop 1
    endif

    if (y <= yBoundLower) then
      linear = bottomValue
    else if (y >= yBoundUpper) then
      linear = topValue
    else
      yRange = yBoundUpper - yBoundLower
      dy = y - yBoundLower
      VRange = topValue - bottomValue
      linear = dy / yRange * VRange + bottomValue
    endif

  end function
! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------
! Calculate cosine gradient between upper and down side
! --------------------------------------------------------------------------------------------------
  doubleprecision function cosineGradient(y, val, yBoundUpper, yBoundLower)
    implicit none

    doubleprecision, intent(in) :: y
    doubleprecision, intent(in) :: val
    doubleprecision, intent(in) :: yBoundUpper
    doubleprecision, intent(in) :: yBoundLower
    doubleprecision :: yRange, dy

! --------------------------------------------------------------------------------------------------
    if (yBoundLower > yBoundUpper) then
      stop 1
    endif

    if (y <= yBoundLower) then
      cosineGradient = val
    else if (y > yBoundUpper) then
      cosineGradient = -val
    else
      yRange = yBoundUpper - yBoundLower
      dy = y - yBoundLower
      cosineGradient = cos(dy / yRange * M_PI) * val
    endif

  end function cosineGradient
! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------
! Solve transport problem
! --------------------------------------------------------------------------------------------------
  doubleprecision function solveTransport(qt, Ef)
    use modscatter
    implicit none
    type(qscatter) :: qt
    doubleprecision :: Ef ! pass in eV
    doubleprecision :: T_total
    integer, parameter :: leadsIds(1) = (/ 1 /)

! --------------------------------------------------------------------------------------------------
    print*,"  Solving transport..."
    call qt%calculate_modes(Ef * eV2au)
    call qt%solve(1, Ef * eV2au)
    ! call qt%solve_leads(leadsIds, Ef * eV2au) ! TEST if it is faster
    T_total = sum(qt%Tn(:))
    print*,"  Total transmission: ", T_total
    solveTransport = T_total

  end function solveTransport
! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------
! Solve transport problem for multiple leads
! --------------------------------------------------------------------------------------------------
  function solveTransportMultilead(qt, Ef, numLeads) result(transmissions)
    use modscatter
    implicit none

    type(qscatter)                :: qt
    doubleprecision               :: Ef ! pass in eV
    integer                       :: numLeads

    doubleprecision, dimension(numLeads, numLeads) :: transmissions

    ! doubleprecision :: T_total
    integer         :: leadId
    integer         :: i, j

    doubleprecision, parameter :: cartesian_vecs(3,3) = (/ (/ 1.0D0, 0.0D0, 0.0D0 /), &
                                                           (/ 0.0D0, 1.0D0, 0.0D0 /), &
                                                           (/ 0.0D0, 0.0D0, 1.0D0 /) /)

! --------------------------------------------------------------------------------------------------
    print*,"  Solving transport..."
    call qt%calculate_modes(Ef * eV2au)

    ! if more then 1 leads, last to calculate is lead 1 so qt%Tn is saved for that one
    do leadId = numLeads, 1, -1
      call qt%solve(leadId, Ef * eV2au)

      if (save_currents) then
        call qt%calculate_currents(leadId)
        call qt%qsystem%save_currents(trim(results_dir)//"/current_"//trim(str(leadId))//".txt", cartesian_vecs)
      endif

    enddo

    if (save_currents) then
      print*, "plotting"
      call execute_command_line("python plot_currents.py "//trim(results_dir))
    endif

    do j = 1 , numLeads ! from
      do i = 1 , numLeads ! to
        ! smatrix(z, do)
        transmissions(i, j) = sum(abs(qt%smatrix(i,j)%Tnm)**2)
        print"(A,i3,A,i3,A,f10.6)","T(",i,",",j,")=", transmissions(i, j)
      enddo
    enddo

    ! Or you can print total tranmission and reflection using auxiliary matrices
    ! but this remember only the last call of qt%solve(lead,Ef) and it contains
    ! summed tranmission probabilies for all leads
    print*,"T     =",sum(qt%Tn(:))
    print*,"R     =",sum(qt%Rn(:))

  end function solveTransportMultilead


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
    call qt%qsystem%save_data(trim(results_dir)//"/densities.xml", &
                              array2d=qt%qsystem%densities, &
                              array1d=qt%qsystem%qauxvec)

  end subroutine calculateElectronDensity
! --------------------------------------------------------------------------------------------------



! --------------------------------------------------------------------------------------------------
! Perform energy scan
! --------------------------------------------------------------------------------------------------
  subroutine performEnergyScan(qt, connect)
    use modscatter
    implicit none

    procedure(func_simple)      :: connect
    type(qscatter)              :: qt
    double precision            :: E_scan, T_total
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
    else
      print*,"  Plotting Transmission..."
      call execute_command_line("python plot_T.py "//trim(results_dir)//"/")
    endif
  end subroutine generatePlots
! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------
! Squares a number
! --------------------------------------------------------------------------------------------------
doubleprecision function pow2(x)
  implicit none
  double precision :: x

  pow2 = x * x

end function pow2

! --------------------------------------------------------------------------------------------------
! Len of vector
! --------------------------------------------------------------------------------------------------
doubleprecision function len(vec)
  implicit none
  double precision :: vec(3)

  len = sqrt(pow2(vec(1)) + pow2(vec(2)) + pow2(vec(3)))

end function len


! --------------------------------------------------------------------------------------------------
! Dot product of vectors
! --------------------------------------------------------------------------------------------------
doubleprecision function dot(vec1, vec2)
  implicit none
  double precision :: vec1(3), vec2(3)

  dot = vec1(1) * vec2(1) + vec1(2) * vec2(2) + vec1(3) * vec2(3)

end function dot


doubleprecision function scaleBasedOnSf(value, sf, maxSfOpt) result(scaled)
  implicit none
  doubleprecision :: value
  integer :: sf
  integer, optional :: maxSfOpt
  integer :: maxSf = 16

  if (present(maxSfOpt)) maxSf = maxSfOpt
  scaled = value * maxSf / sf

end function scaleBasedOnSf


integer function scaleBasedOnSfInt(value, sf, maxSfOpt) result(scaled)
  implicit none
  integer :: value
  integer :: sf
  integer, optional :: maxSfOpt
  integer :: maxSf = 16

  if (present(maxSfOpt)) maxSf = maxSfOpt
  scaled = value * maxSf / sf

end function scaleBasedOnSfInt



endmodule