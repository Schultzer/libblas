program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   complex*16 x(10)
   call json%initialize()
   call json%create_array(root, '')

   x = (/(0, 0.76359346114045956),&
         &(0, 0),&
         &(1.3297992629225006134, -1.14765700923635139),&
         &(1.2724293214294046805, -0.28946157368822334),&
         &(0.4146414344564082199, -0.29921511789731614),&
         &(-1.5399500419037095433, -0.41151083279506701),&
         &(-0.9285670347135380753, 0.25222344815613229),&
         &(-0.2947204467905601977, -0.89192112728456863),&
         &(-0.0057671727475369552, 0.43568329935571865),&
         &(2.4046533888579508798, -1.23753842192995811)/)
   call addcase(root, '0', 10, x, 1)

   call addcase(root, '1', 0, x, 1)

   call addcase(root, '2', 10, x, 0)

   call print(root)
contains
   subroutine addcase(root, ncase, n, x, incx)
      type(json_core) :: json
      type(json_value), pointer :: root, case, array, data
      external dznrm2
      character(len=*) ncase
      integer n, incx
      complex*16 x(10)
      double precision dznrm2
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%create_array(array, 'x')
      do i = 1, 10
         call json%create_array(data, '')
         call json%add(data, '', dreal(x(i)))
         call json%add(data, '', dimag(x(i)))
         call json%add(array, data)
         nullify (data)
      enddo
      call json%add(case, array)
      nullify (array)
      call json%add(case, 'incx', incx)
      call json%add(case, 'expect', dznrm2(n, x, incx))
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
      call json%print(root, './tests/fixtures/level1/complex/nrm2.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
