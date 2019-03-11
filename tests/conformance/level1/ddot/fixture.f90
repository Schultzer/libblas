program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real32
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   real x(4), y(4), b
   call json%initialize()
   call json%create_array(root, '')

   x = (/1, 2, 3, 4/)
   y = (/5, 6, 7, 8/)
   b = 2
   call addcase(root, '0', 4, b, x, 1, y, 1)

   x = (/1, 2, 3, 4/)
   y = (/5, 6, 7, 8/)
   call addcase(root, '1', 0, b, x, 1, y, 1)

   x = (/1, 2, 3, 4/)
   y = (/5, 6, 7, 8/)
   b = 3
   call addcase(root, '2', 2, b, x, 2, y, 2)

   x = (/1, 2, 3, 4/)
   y = (/5, 6, 7, 8/)
   call addcase(root, '3', 2, b, x, -1, y, -1)

   x = (/1, 2, 3, 4/)
   y = (/5, 6, 7, 8/)
   call addcase(root, '4', 2, b, x, -1, y, 1)

   call print(root)
contains
   subroutine addcase(root, ncase, n, b, x, incx, y, incy)
      use, intrinsic :: iso_fortran_env, only: wp => real32
      type(json_core) :: json
      type(json_value), pointer :: root, case, array
      external sdsdot
      character(len=*) ncase
      integer n, incx, incy
      real x(4), y(4), b
      double precision dot
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%add(case, 'b', dble(b))
      call json%create_array(array, 'x')
      do i = 1, 4
         call json%add(array, '', dble(x(i)))
      enddo
      call json%add(case, array)
      nullify (array)
      call json%add(case, 'incx', incx)
      call json%create_array(array, 'y')
      do i = 1, 4
         call json%add(array, '', dble(y(i)))
      enddo
      call json%add(case, array)
      nullify (array)
      call json%add(case, 'incy', incy)
      dot = sdsdot(n, b, x, incx, y, incy)
      call json%add(case, 'expect', dot)
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
      call json%print(root, './tests/fixtures/level1/ddot.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
