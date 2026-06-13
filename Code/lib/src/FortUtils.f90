module FortUtils
  implicit none

  integer :: currentArgIndex = 1  ! Global counter for command line arguments

  ! Variables used by parseArguments - declared here so they can be accessed from main program
  character(len=512) :: results_dir = "./results"
  doubleprecision :: Bz = 8.0D0
  doubleprecision :: Vb = -20.0D0
  doubleprecision :: Vt = 0.0D0
  logical :: save_system = .false.
  logical :: run_transport = .false.
  logical :: run_energyScan = .false.
  logical :: plot_results = .false.
  logical :: save_densities = .false.
  logical :: save_bands = .false.
  integer :: sf = 8

contains

  ! Get next command line argument
  ! --------------------------------------------------------------------------------------------------
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
    save_system = parseBoolArg(save_system)
    run_transport = parseBoolArg(run_transport)
    run_energyScan = parseBoolArg(run_energyScan)
    plot_results = parseBoolArg(plot_results)
    save_densities = parseBoolArg(save_densities)
    save_bands = parseBoolArg(save_bands)
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



end module
