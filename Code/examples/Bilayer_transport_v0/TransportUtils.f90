module TransportUtils
  use FortUtils
  implicit none

  ! unit conversion
  doubleprecision,parameter :: T2au        = 4.254382E-6          ! B(au) = B(T)*T2au
  doubleprecision,parameter :: eV2au       = 0.03674932587122423  ! V(au)  = V(eV)*eV2au
  doubleprecision,parameter :: nm2au       = 1.0 / 0.0529           ! d(au)  = d(nm)*nm2au
  doubleprecision,parameter :: cm2au       = 1e-2 * 1e9 * nm2au
  doubleprecision,parameter :: inv_cmsq2au = 1. / cm2au / cm2au
  doubleprecision,parameter :: one_over_sqrt_3 = 1.0D0 / sqrt(3.0)

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
! Calculate linear gradiend between upper and down side
! --------------------------------------------------------------------------------------------------
  doubleprecision function linear(y, topValue, bottomValue, yBoundUpper, yBoundLower)
    implicit none

    doubleprecision, intent(in) :: y
    doubleprecision, intent(in) :: topValue
    doubleprecision, intent(in) :: bottomValue
    doubleprecision, intent(in) :: yBoundUpper
    doubleprecision, intent(in) :: yBoundLower
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
! Save calculated resdensitiesults
! --------------------------------------------------------------------------------------------------
  subroutine saveDensities(qt)
    use modscatter
    implicit none
    type(qscatter) :: qt

! --------------------------------------------------------------------------------------------------
    print*,"  Saving results Densities..."
    call qt%qsystem%save_data(trim(results_dir)//"/densities.xml", &
                              array2d=qt%qsystem%densities, &
                              array1d=qt%qsystem%qauxvec)
  end subroutine saveDensities
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
    endif
    print*,"  Plotting Transmission..."
    call execute_command_line("python plot_T.py "//trim(results_dir)//"/")
  end subroutine generatePlots
! --------------------------------------------------------------------------------------------------

endmodule