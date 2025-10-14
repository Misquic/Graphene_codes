module Bilayer_interface
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    function Bilayer_constructor(nit, nib, &
                                 dt, dg, db, &
                                 et, eg, eb) & 
                                 bind(C, name="Bilayer_constructor") &
                                 result(this_ptr)
      import :: c_double, c_ptr
      real(c_double), intent(in), value :: nit, nib
      real(c_double), intent(in), value :: dt, dg, db
      real(c_double), intent(in), value :: et, eg, eb
      type(c_ptr) :: this_ptr
    end function Bilayer_constructor

    subroutine Bilayer_destructor(this_ptr) bind(C, name="Bilayer_destructor")
      import :: c_ptr
      type(c_ptr), intent(in), value :: this_ptr
    end subroutine Bilayer_destructor

    subroutine Bilayer_countDensities(this_ptr, &
                                      Vt, Vb, &
                                      nt, nb) &
                                      bind(C, name="Bilayer_countDensities")
      import :: c_double, c_ptr
      type(c_ptr), intent(in), value :: this_ptr
      real(c_double), intent(in), value :: Vt, Vb
      real(c_double), intent(out) :: nt, nb
    end subroutine Bilayer_countDensities

    subroutine Bilayer_countDensities_B(this_ptr, &
                                        Vt, Vb, B, &
                                        nt, nb) &
                                        bind(C, name="Bilayer_countDensities")
      import :: c_double, c_ptr
      type(c_ptr), intent(in), value :: this_ptr
      real(c_double), intent(in), value :: Vt, Vb, B
      real(c_double), intent(out) :: nt, nb
    end subroutine Bilayer_countDensities_B

    function count_E_0(n0, B) bind(C, name="count_E_0") result(E_0)
      import :: c_double
      real(c_double), intent(in), value :: n0, B
      real(c_double) :: E_0
    end function
  end interface

end module Bilayer_interface