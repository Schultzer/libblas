program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   complex*16 x(6), a
   call json%initialize()
   call json%create_array(root, '')

   x = (/(1.26295428488079331, -0.9285670347135380753),&
      &(-0.32623336070564940, -0.2947204467905601977),&
      &(1.32979926292250061, -0.0057671727475369552),&
      &(1.27242932142940468, 2.4046533888579508798),&
      &(0.41464143445640822, 0.7635934611404595618),&
      &(-1.53995004190370954, -0.7990092489893682037)/)
   a = (2, 4)
   call addcase(root, '0', 6, a, x, 1)

   x = (/(1.26295428488079331, -0.9285670347135380753),&
      &(-0.32623336070564940, -0.2947204467905601977),&
      &(1.32979926292250061, -0.0057671727475369552),&
      &(1.27242932142940468, 2.4046533888579508798),&
      &(0.41464143445640822, 0.7635934611404595618),&
      &(-1.53995004190370954, -0.7990092489893682037)/)
   a = (2, 4)
   call addcase(root, '1', 3, a, x, 2)

   x = (/(1.26295428488079331, -0.9285670347135380753),&
      &(-0.32623336070564940, -0.2947204467905601977),&
      &(1.32979926292250061, -0.0057671727475369552),&
      &(1.27242932142940468, 2.4046533888579508798),&
      &(0.41464143445640822, 0.7635934611404595618),&
      &(-1.53995004190370954, -0.7990092489893682037)/)
   a = (2, 4)
   call addcase(root, '2', 0, a, x, 2)

   x = (/(1.26295428488079331, -0.9285670347135380753),&
      &(-0.32623336070564940, -0.2947204467905601977),&
      &(1.32979926292250061, -0.0057671727475369552),&
      &(1.27242932142940468, 2.4046533888579508798),&
      &(0.41464143445640822, 0.7635934611404595618),&
      &(-1.53995004190370954, -0.7990092489893682037)/)
   a = (2, 4)
   call addcase(root, '3', 3, a, x, 0)

   call print(root)
contains
   subroutine addcase(root, ncase, n, a, x, incx)
      type(json_core) :: json
      type(json_value), pointer :: root, case, array, data
      external zscal
      character(len=*) ncase
      integer n, incx
      complex*16 x(6), a
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%create_array(data, 'a')
      call json%add(data, '', dreal(a))
      call json%add(data, '', dimag(a))
      call json%add(case, data)
      nullify (data)
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
      call zscal(n, a, x, incx)
      call json%create_array(array, 'expect')
      do i = 1, 6
         call json%create_array(data, '')
         call json%add(data, '', dreal(x(i)))
         call json%add(data, '', dimag(x(i)))
         call json%add(array, data)
         nullify (data)
      enddo
      call json%add(case, array)
      nullify (array)
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
      call json%print(root, './tests/fixtures/level1/complex/scal.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
