program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   complex*16 a, b, s
   double precision c
   call json%initialize()
   call json%create_array(root, '')

   a = (11, 19)
   b = (34, 23)
   c = 0
   s = (0, 0)
   call addcase(root, '0', a, b, c, s)

   a = (0, 0)
   b = (-1, -1)
   c = 0
   s = (0, 0)
   call addcase(root, '1', a, b, c, s)
   call print(root)
contains
   subroutine addcase(root, ncase, a, b, c, s)
      type(json_core) :: json
      type(json_value), pointer :: root, case, data, d
      external zrotg
      character(len=*) ncase
      complex*16 a, b, s
      double precision c
      call json%create_object(case, '')
      call json%add(root, case)
      call json%create_array(data, 'a')
      call json%add(data, '', dreal(a))
      call json%add(data, '', dimag(a))
      call json%add(case, data)
      nullify (data)
      call json%create_array(data, 'b')
      call json%add(data, '', dreal(b))
      call json%add(data, '', dimag(b))
      call json%add(case, data)
      nullify (data)
      call json%add(case, 'c', c)
      call json%create_array(data, 's')
      call json%add(data, '', dreal(s))
      call json%add(data, '', dimag(s))
      call json%add(case, data)
      nullify (data)
      call zrotg(a, b, c, s)
      call json%create_object(data, 'expect')
      call json%create_array(d, 'a')
      call json%add(d, '', dreal(a))
      call json%add(d, '', dimag(a))
      call json%add(data, d)
      nullify (d)
      call json%create_array(d, 'b')
      call json%add(d, '', dreal(b))
      call json%add(d, '', dimag(b))
      call json%add(data, d)
      nullify (d)
      call json%add(data, 'c', c)
      call json%create_array(d, 's')
      call json%add(d, '', dreal(s))
      call json%add(d, '', dimag(s))
      call json%add(data, d)
      nullify (d)
      call json%add(case, data)
      nullify (data)
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
      call json%print(root, './tests/fixtures/level1/complex/rotg.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
