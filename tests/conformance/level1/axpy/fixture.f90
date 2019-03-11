program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real128
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   double precision x(10), y(10)
   call json%initialize()
   call json%create_array(root, '')

   x = (/1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134,&
         &1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,&
         &-0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,&
         &2.4046533888579508798/)
   y = (/-0.560475646552212603, -0.230177489483279957, 1.558708314149124030,&
         &0.070508391424576003, 0.129287735160946243, 1.715064986883281017,&
         &0.460916205989202299, -1.265061234606533969, -0.686852851893526073,&
         &-0.445661970099958060/)
   call addcase(root, '0', 10, 23.d0, x, 1, y, 1)

   x = (/1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134, .0, .0, .0, .0, .0, .0, .0/)
   y = (/-0.560475646552212603, -0.230177489483279957, 1.558708314149124030, .0, .0, .0, .0, .0, .0, .0/)
   call addcase(root, '0a', 3, 23.d0, x, 1, y, 1)

   x = (/1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134,&
         &1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,&
         &-0.9285670347135380753, -0.2947204467905601977, .0, .0/)
   y = (/-0.560475646552212603, -0.230177489483279957, 1.558708314149124030,&
         &0.070508391424576003, 0.129287735160946243, 1.715064986883281017,&
         &0.460916205989202299, -1.265061234606533969, .0, .0/)
   call addcase(root, '1', 8, 23.d0, x, 1, y, 1)

   x = (/1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134,&
         &1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,&
         &-0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,&
         &2.4046533888579508798/)
   y = (/-0.560475646552212603, -0.230177489483279957, 1.558708314149124030,&
         &0.070508391424576003, 0.129287735160946243, 1.715064986883281017,&
         &0.460916205989202299, -1.265061234606533969, -0.686852851893526073,&
         &-0.445661970099958060/)
   call addcase(root, '2', 5, 23.d0, x, 2, y, 2)

   x = (/1.2629542848807933098, -0.3262333607056494000, 1.3297992629225006134,&
         &1.2724293214294046805, 0.4146414344564082199, -1.5399500419037095433,&
         &-0.9285670347135380753, -0.2947204467905601977, -0.0057671727475369552,&
         &2.4046533888579508798/)
   y = (/-0.560475646552212603, -0.230177489483279957, 1.558708314149124030,&
         &0.070508391424576003, 0.129287735160946243, 1.715064986883281017,&
         &0.460916205989202299, -1.265061234606533969, -0.686852851893526073,&
         &-0.445661970099958060/)
   call addcase(root, '5', 5, 23.d0, x, -2, y, -2)

   call print(root)
contains
   subroutine addcase(root, ncase, n, a, x, incx, y, incy)
      use, intrinsic :: iso_fortran_env, only: wp => real128
      type(json_core) :: json
      type(json_value), pointer :: root, case
      external daxpy
      character(len=*) ncase
      integer n, incx, incy
      double precision x(10), y(10), a
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%add(case, 'a', a)
      call json%add(case, 'x', x)
      call json%add(case, 'incx', incx)
      call json%add(case, 'y', y)
      call json%add(case, 'incy', incy)
      call daxpy(n, a, x, incx, y, incy)
      call json%add(case, 'expect', y)
      nullify (case)
   end subroutine addcase

   subroutine print(root)
      use, intrinsic :: iso_fortran_env, only: real64
      use json_module, CK => json_CK, CK => json_CK
      implicit none
      type(json_core) :: json
      type(json_value), pointer :: root
      logical :: status_ok
      character(kind=CK, len=:), allocatable :: error_msg
      call json%print(root, './tests/fixtures/level1/axpy.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
