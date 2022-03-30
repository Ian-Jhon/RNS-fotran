! converts EOS files from MeV/fm^3 to CGS units

      program unitconvert 

          implicit none
          integer :: i, k
          double precision, allocatable :: e(:), p(:)
          double precision, allocatable :: enew(:), pnew(:) 
          character :: dummy

          ! open file to read - update with file name
          open(unit=10, file='unconverted.dat', status='unknown')

          ! open file to output - update with file name
          open(unit=20, file='converted.dat', status='unknown')

          read(10,*) dummy
          read(10,*) k 

          allocate(e(k)); allocate(p(k))
          allocate(pnew(k)); allocate(enew(k))

          do i = 3, k+2
            read(10,*) e(i), p(i), dummy
          end do

          do i = 1, k
            enew(i) = e(i) * 1.7827E12
            pnew(i) = p(i) * 1.6022E33
            write(20,*) enew(i), pnew(i)
          end do

          close(unit=10)
          close(unit=20)

      end program
