program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   double precision x(7), a
   call json%initialize()
   call json%create_array(root, '')

   x = (/1, 2, 3, 4, 0, 0, 0/)
   a = 1
   call addcase(root, '0', 4, a, x, 1)

   a = 2
   x = (/1, 2, 3, 4, 0, 0, 0/)
   call addcase(root, '1', 4, a, x, 1)

   x = (/3, 1, 2, 3, 4, 0, 0/)
   call addcase(root, '2', 5, a, x, 1)

   x = (/-1, 2, 3, 1, 2, 3, 4/)
   call addcase(root, '4', 7, a, x, 1)

   x = (/-1, 2, 3, 1, 2, 3, 4/)
   call addcase(root, '5', 3, a, x, 2)

   x = (/-1, 2, 3, 1, 2, 3, 4/)
   call addcase(root, '6', 6, a, x, 0)

   x = (/-1, 2, 3, 1, 2, 3, 4/)
   call addcase(root, '7', 0, a, x, 2)

   call print(root)
contains
   subroutine addcase(root, ncase, n, a, x, incx)
      type(json_core) :: json
      type(json_value), pointer :: root, case
      external dscal
      character(len=*) ncase
      integer n, incx
      double precision x(7), a
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%add(case, 'a', a)
      call json%add(case, 'x', x)
      call json%add(case, 'incx', incx)
      call dscal(n, a, x, incx)
      call json%add(case, 'expect', x)
      nullify (case)
   end subroutine addcase

   subroutine print(root)
      use, intrinsic :: iso_fortran_env, only: real64
      use json_module, CK => json_CK
      implicit none

      type(json_core) :: json
      type(json_value), pointer :: root
      logical :: status_ok
      character(kind=CK, len=:), allocatable :: error_msg
      call json%print(root, './tests/fixtures/level1/scal.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
