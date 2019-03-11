program fixture
   use, intrinsic :: iso_fortran_env, only: wp => real64
   use json_module
   type(json_core) :: json
   type(json_value), pointer :: root
   complex*16 x(10), y(10), x2(1), y2(1), a
   call json%initialize()
   call json%create_array(root, '')

   a = (23, 32)
   x = (/(1.2629542848807933098, 0.76359346114045956),&
         &(-0.3262333607056494000, -0.79900924898936820),&
         &(1.3297992629225006134, -1.14765700923635139),&
         &(1.2724293214294046805, -0.28946157368822334),&
         &(0.4146414344564082199, -0.29921511789731614),&
         &(-1.5399500419037095433, -0.41151083279506701),&
         &(-0.9285670347135380753, 0.25222344815613229),&
         &(-0.2947204467905601977, -0.89192112728456863),&
         &(-0.0057671727475369552, 0.43568329935571865),&
         &(2.4046533888579508798, -1.23753842192995811)/)

   y = (/(-0.560475646552212603, 1.22408179743946155),&
         &(-0.230177489483279957, 0.35981382705736381),&
         &(1.558708314149124030, 0.40077145059405217),&
         &(0.070508391424576003, 0.11068271594511971),&
         &(0.129287735160946243, -0.55584113475407493),&
         &(1.715064986883281017, 1.78691313680307817),&
         &(0.460916205989202299, 0.49785047822923939),&
         &(-1.265061234606533969, -1.96661715662963821),&
         &(-0.686852851893526073, 0.70135590156368555),&
         &(-0.445661970099958060, -0.47279140772793404)/)
   call addcase(root, '0', 10, a, x, 1, y, 1)

   y = (/(-0.560475646552212603, 1.22408179743946155),&
         &(-0.230177489483279957, 0.35981382705736381),&
         &(1.558708314149124030, 0.40077145059405217),&
         &(0.070508391424576003, 0.11068271594511971),&
         &(0.129287735160946243, -0.55584113475407493),&
         &(1.715064986883281017, 1.78691313680307817),&
         &(0.460916205989202299, 0.49785047822923939),&
         &(-1.265061234606533969, -1.96661715662963821),&
         &(-0.686852851893526073, 0.70135590156368555),&
         &(-0.445661970099958060, -0.47279140772793404)/)
   call addcase(root, '1', 5, (23.0_8, 32.0_8), x, 2, y, 2)

   a = (0, 0)
   x2 = (/(3, 4)/)
   y2 = (/(1, 2)/)
   call addcase(root, '2', 1, a, x2, 1, y2, 1)

   a = (1, 2)
   x2 = (/(3, 4)/)
   y2 = (/(1, 2)/)
   call addcase(root, '3', 0, a, x2, 1, y2, 1)

   a = (23, 32)
   x = (/(1.2629542848807933098, 0.76359346114045956),&
      &(-0.3262333607056494000, -0.79900924898936820),&
      &(1.3297992629225006134, -1.14765700923635139),&
      &(1.2724293214294046805, -0.28946157368822334),&
      &(0.4146414344564082199, -0.29921511789731614),&
      &(-1.5399500419037095433, -0.41151083279506701),&
      &(-0.9285670347135380753, 0.25222344815613229),&
      &(-0.2947204467905601977, -0.89192112728456863),&
      &(-0.0057671727475369552, 0.43568329935571865),&
      &(2.4046533888579508798, -1.23753842192995811)/)

   y = (/(-0.560475646552212603, 1.22408179743946155),&
   &(-0.230177489483279957, 0.35981382705736381),&
   &(1.558708314149124030, 0.40077145059405217),&
   &(0.070508391424576003, 0.11068271594511971),&
   &(0.129287735160946243, -0.55584113475407493),&
   &(1.715064986883281017, 1.78691313680307817),&
   &(0.460916205989202299, 0.49785047822923939),&
   &(-1.265061234606533969, -1.96661715662963821),&
   &(-0.686852851893526073, 0.70135590156368555),&
   &(-0.445661970099958060, -0.47279140772793404)/)
   call addcase(root, '4', 5, a, x, -2, y, -2)

   call print(root)
contains
   subroutine addcase(root, ncase, n, a, x, incx, y, incy)
      type(json_core) :: json
      type(json_value), pointer :: root, case, array, data
      external zaxpy
      character(len=*) ncase
      integer n, incx, incy
      complex*16 x(*), y(*), a
      call json%create_object(case, '')
      call json%add(root, case)
      call json%add(case, 'n', n)
      call json%create_array(data, 'a')
      call json%add(data, '', dreal(a))
      call json%add(data, '', dimag(a))
      call json%add(case, data)
      nullify (data)
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
      call json%create_array(array, 'y')
      do i = 1, 10
         call json%create_array(data, '')
         call json%add(data, '', dreal(y(i)))
         call json%add(data, '', dimag(y(i)))
         call json%add(array, data)
         nullify (data)
      enddo
      call json%add(case, array)
      nullify (array)
      call json%add(case, 'incy', incy)
      call zaxpy(n, a, x, incx, y, incy)
      call json%create_array(array, 'expect')
      do i = 1, 10
         call json%create_array(data, '')
         call json%add(data, '', dreal(y(i)))
         call json%add(data, '', dimag(y(i)))
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
      call json%print(root, './tests/fixtures/level1/complex/axpy.json')
      if (json%failed()) then
         call json%check_for_errors(status_ok, error_msg)
         write (*, *) 'Error: '//error_msg
         call json%clear_exceptions()
         call json%destroy(root)
      end if
   end subroutine print
end
