program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   complex*16 x(6)

   call json%initialize()
   call json%create_array(root, '')
   x = (/(1, 7), (0, 0), (3, 9), (4, 10), (5, 11), (6, 12)/)
   call addcase(root, '0', 6, x, 1)

   x = (/(1, 7), (2, 8), (3, 9), (4, 10), (5, 11), (6, 12)/)
   call addcase(root, '1', 3, x, 2)

   x = (/(1, 7), (2, 8), (3, 9), (4, 10), (5, 11), (6, 12)/)
   call addcase(root, '2', 6, x, 0)

   x = (/(1, 7), (2, 8), (3, 9), (4, 10), (5, 11), (6, 12)/)
   call addcase(root, '3', 0, x, 1)

   x = (/(1, 7), (2, 8), (3, 9), (4, 10), (5, 11), (6, 12)/)
   call addcase(root, '4', 1, x, 1)

   call print(root)
contains
   subroutine addcase(root, ncase, n, x, incx)
      type(json_core) :: json
      type(json_value), pointer :: root, case, array, data
      external izamax
      character(len=*) ncase
      integer n, incx, max
      complex*16 x(6)
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%create_array(array, 'x')
      do i = 1, 6
         call json%create_array(data, '')
         call json%add(data, '', dreal(x(i)))
         call json%add(data, '', dimag(x(i)))
         call json%add(array, data)
         nullify (data)
      enddo
      call json%add(case, array)
      nullify (array)
      call json%add(case, 'incx', incx)
      max = izamax(n, x, incx)
      call json%add(case, 'expect', max)
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
      call json%print(root, './tests/fixtures/level1/complex/iamax.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
