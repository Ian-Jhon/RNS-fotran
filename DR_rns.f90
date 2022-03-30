!--------------------------------------------------------------------------
! delaney farrell 9/1/21
! differential rotation code for neutron stars -  version 1.0
! needs to be compiled with: iteration.f90, TOV_guess.f90, comp.f90, and
! comp_f_p.f90, and spin_up_freq.f90 (see makefile)  
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! global variables 
!--------------------------------------------------------------------------

      module vars

          implicit none 

          integer :: s_div = 329, m_div = 125, r_div = 900
          integer :: l_max = 10 
          integer :: n_nearest, a_check = 0
          integer :: print_dif = 0, n_tab = 0, print_option 

          double precision :: s_max = 0.9999d0, gamma_p, cf
          double precision :: e_center, p_center 
          double precision :: r_final = 0.d0, m_final = 0.d0
          double precision :: pi = dacos(-1.d0) 
          double precision :: r_min = 1.d-15
          double precision :: mb = 1.66d-24
          double precision :: r_is_final = 0.d0
          double precision :: e_surface, p_surface, k_rescale 
          double precision :: c = 2.9979d10
          double precision :: G = 6.6832d-8
          double precision :: kappa = 1.346790806509621d+13
          double precision :: k_scale = 1.112668301525780d-36
          double precision :: s_e = 0.5d0
          double precision :: m_sun = 1.987e33
          double precision :: mass, r_e, r_e_guess, r_ec, z_p 
          double precision :: gamma_center_old, rho_center_old
          double precision :: gamma_pole_h, gamma_center_h
          double precision :: gamma_equator_h
          double precision :: rho_pole_h, rho_center_h, rho_equator_h 
          double precision :: r_p, s_p, h_center
          double precision :: omega_equator_h, omega_h
          double precision :: enthalpy_min
          double precision :: accu, ds, dm 
          double precision :: r_ratio, velocity_equator 
          double precision :: mass_p, mass_0, jam, n_p, omega_f
          double precision :: W, T
          double precision :: z_b, z_f, im, omega_k, kappa_pole
          double precision :: omega_error, a, m0_error, m0, j_const

          double precision, allocatable :: s_gp(:), s1_s(:), one_s(:) 
          double precision, allocatable :: mu(:), one_m2(:), theta(:) 
          double precision, allocatable :: sin_theta(:) 
          double precision, allocatable :: f_rho(:,:,:) 
          double precision, allocatable :: f_gamma(:,:,:)
          double precision, allocatable :: f_omega(:,:,:)
          double precision, allocatable :: p_2n(:,:), p1_2n_1(:,:) 
          double precision, allocatable :: r_gp(:), r_is_gp(:), m_gp(:) 
          double precision, allocatable :: lambda_gp(:), e_d_gp(:)
          double precision, allocatable :: nu_gp(:) 
          double precision, allocatable :: gamma_guess(:,:)
          double precision, allocatable :: rho_guess(:,:)
          double precision, allocatable :: omega_guess(:,:) 
          double precision, allocatable :: alpha_guess(:,:) 
          double precision, allocatable :: alpha(:,:), gam(:,:) 
          double precision, allocatable :: rho(:,:), omega(:,:)
          double precision, allocatable :: gamma_mu_0(:), rho_mu_0(:) 
          double precision, allocatable :: omega_mu_0(:)
          double precision, allocatable :: omega_h_mu_0(:) 
          double precision, allocatable :: gamma_mu_1(:), rho_mu_1(:) 
          double precision, allocatable :: velocity_sq(:,:)
          double precision, allocatable :: enthalpy(:,:)
          double precision, allocatable :: energy(:,:) 
          double precision, allocatable :: pressure(:,:) 
          double precision, allocatable :: s_rho(:,:), s_gamma(:,:) 
          double precision, allocatable :: s_omega(:,:) 
          double precision, allocatable :: d1_rho(:,:), d1_gamma(:,:) 
          double precision, allocatable :: d1_omega(:,:) 
          double precision, allocatable :: d2_rho(:,:), d2_gamma(:,:)
          double precision, allocatable :: d2_omega(:,:)
          double precision, allocatable :: dgds(:,:), dgdm(:,:)
          double precision, allocatable :: dadm(:,:)
          double precision, allocatable :: log_e_tab(:), log_p_tab(:) 
          double precision, allocatable :: log_h_tab(:), log_n0_tab(:) 

          character :: dummy 

      end module 

!--------------------------------------------------------------------------
! main program 
!--------------------------------------------------------------------------

      program dr_rns

          use vars

          implicit none 
          integer :: ialloc, n_of_models, option, j, I 
          double precision :: err_tol
          double precision :: e_min, e_max, m_fix
          double precision :: m0_const, omega_const, jconst 
          double precision :: ratio_max, ratio_min, ratio_loop
          double precision :: dratio 

          n_nearest = 1 

          ! allocate space for variables 
          allocate(dadm(0:s_div, 0:m_div), stat=ialloc) 
          allocate(dgds(0:s_div, 0:m_div), stat=ialloc)
          allocate(dgdm(0:s_div, 0:m_div), stat=ialloc)
          allocate(d2_rho(0:s_div, 0:l_max), stat=ialloc)
          allocate(d2_gamma(0:s_div, 0:l_max), stat=ialloc)
          allocate(d2_omega(0:s_div, 0:l_max), stat=ialloc)
          allocate(d1_rho(0:l_max, 0:s_div), stat=ialloc)
          allocate(d1_gamma(0:l_max, 0:s_div), stat=ialloc)
          allocate(d1_omega(0:l_max, 0:s_div), stat=ialloc)
          allocate(s_rho(0:s_div, 0:m_div), stat=ialloc)
          allocate(s_omega(0:s_div, 0:m_div), stat=ialloc)
          allocate(s_gamma(0:s_div, 0:m_div), stat=ialloc)
          allocate(pressure(0:s_div, 0:m_div), stat=ialloc)
          allocate(energy(0:s_div, 0:m_div), stat=ialloc)
          allocate(enthalpy(0:s_div, 0:m_div), stat=ialloc)
          allocate(velocity_sq(0:s_div, 0:m_div), stat=ialloc)
          allocate(omega_guess(0:s_div, 0:m_div), stat=ialloc)
          allocate(alpha_guess(0:s_div, 0:m_div), stat=ialloc)
          allocate(gamma_guess(0:s_div, 0:m_div), stat=ialloc)
          allocate(rho_guess(0:s_div, 0:m_div), stat=ialloc)
          allocate(gam(0:s_div, 0:m_div), stat=ialloc)
          allocate(rho(0:s_div, 0:m_div), stat=ialloc)
          allocate(alpha(0:s_div, 0:m_div), stat=ialloc)
          allocate(omega(0:s_div, 0:m_div), stat=ialloc)
          allocate(rho_mu_0(0:s_div), stat=ialloc)
          allocate(rho_mu_1(0:s_div), stat=ialloc)
          allocate(gamma_mu_0(0:s_div), stat=ialloc)
          allocate(gamma_mu_1(0:s_div), stat=ialloc)
          allocate(omega_mu_0(0:s_div), stat=ialloc)
          allocate(omega_h_mu_0(0:s_div), stat=ialloc) 
          allocate(r_is_gp(0:r_div), stat=ialloc)
          allocate(r_gp(0:r_div), stat=ialloc)
          allocate(m_gp(0:r_div), stat=ialloc)
          allocate(lambda_gp(0:r_div), stat=ialloc)
          allocate(e_d_gp(0:r_div), stat=ialloc)
          allocate(nu_gp(0:r_div), stat=ialloc)
          allocate(s_gp(0:s_div), stat=ialloc)
          allocate(s1_s(0:s_div), stat=ialloc)
          allocate(one_s(0:s_div), stat=ialloc)
          allocate(mu(0:m_div), stat=ialloc)
          allocate(one_m2(0:m_div), stat=ialloc)
          allocate(theta(0:m_div), stat=ialloc)
          allocate(sin_theta(0:m_div), stat=ialloc)
          allocate(f_rho(0:s_div, 0:l_max, 0:s_div), stat=ialloc)
          allocate(f_gamma(0:s_div, 0:l_max, 0:s_div), stat=ialloc)
          allocate(f_omega(0:s_div, 0:l_max, 0:s_div), stat=ialloc)
          allocate(p_2n(0:m_div, 0:l_max), stat=ialloc)
          allocate(p1_2n_1(0:m_div, 0:l_max), stat=ialloc)

          ds = s_max/(s_div-1.d0)
          dm = 1.d0/(m_div-1.d0) 

          ! open parameter file 
          open(unit=21, file='parameters.dat', status='unknown') 

          ! parameters for tabulated EOS
          e_surface = 7.8*c*c*k_scale 
          p_surface = 1.01e8*k_scale 
          enthalpy_min = 1.0/(c*c) 

          ! call functions
          call make_grid
          call load_eos
          call comp_f_p

          ! polytropic stars default 
          cf = 1.d0                     ! relaxation constant
          accu = 1.d-4                  ! accuracy 
          err_tol = 1.d-4               ! error tolerance 
          n_of_models = 1               ! # of models to compute
          e_min = 2.66e15*c*c*k_scale   ! min energy density 
          e_max = 1.e16*c*c*k_scale     ! max energy density 
          r_ratio = 0.75                ! ratio r_p / r_e 
          m_fix = 1.4*m_sun             ! gravitational mass
          m0_const = 1.4*m_sun          ! baryonic mass
          omega_const = 1.e-5           ! omega
          j_const = 0.                  ! angular momentum 

          ! reading options (dummy = skipping a line) 
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy     
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) option

          !-----------------------------------------------------------------
          ! program for option 1 (single star) 
          if (option .eq. 1) then 

              ! read suboptions 
              read(21,*) dummy
              read(21,*) e_center, r_ratio, print_option 

              ! screen output 
              print*, '---------------------------------------------'
              print*, 'begin program for tabulated EOS'
              print*, 'option: single star'
              print*, 'central density = ', e_center, ' MeV/fm^3'
              print*, 'ratio r_p/r_e = ', r_ratio 
              print*, '---------------------------------------------'

              ! convert units 
              e_center = e_center*1.7828d12       ! MeV/fm^3 to g/cm^3
              e_center = e_center*c*c*k_scale     ! g/cm^3 to dimensionless

              ! call routines for single star
              call center(e_center)
              call TOV_guess
              call iterate
              call comp_phys
              call print_str

          !-----------------------------------------------------------------    
          ! program for option 2 (kepler frequency sequence)     
          else if (option .eq. 2) then

              ! read suboptions 
              read(21,*) dummy
              read(21,*) e_min, e_max, n_of_models

              ! frequency tolerance 
              omega_error = err_tol 

              ! screen output 
              print*, '---------------------------------------------'
              print*, 'begin program for solving sequence'
              print*, 'option: kepler frequency stars'
              print*, 'number of stars = ', n_of_models
              print*, '---------------------------------------------'

              ! convert units 
              e_min = e_min*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_min = e_min*c*c*k_scale        ! g/cm^3 to dimensionless
              e_max = e_max*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_max = e_max*c*c*k_scale        ! g/cm^3 to dimensionless

            write(*,"(a25)",advance='no') 'central density (MeV/fm^3),'
            write(*,"(a25)",advance='no') 'stellar mass (m_sun),'
            write(*,"(a25)",advance='no') 'baryonic mass (m_sun),'
            write(*,"(a25)",advance='no') 'equatorial radius (km),'
            write(*,"(a25)",advance='no') 'polar radius (km),'
            write(*,"(a25)",advance='no') 'frequency (hz),'
            write(*,"(a25)") 'kepler frequency (hz)'


              ! density increments 
              a = (e_max/e_min)**(1./(n_of_models-1.))

              do j = 1, n_of_models

                e_center = (a**(1.*j-1.))*e_min 

                ! call model for kepler frequency
                call kepler_models

              end do 

          !-----------------------------------------------------------------
          ! program for option 3 (constant baryon mass)
          else if (option .eq. 3) then

              ! read suboptions
              read(21,*) dummy 
              read(21,*) m0, e_min, e_max, n_of_models, print_option

              ! frequency tolerance 
              m0_error = err_tol 

              ! convert mass
              m0 = m0*m_sun

              ! screen output 
              print*, '---------------------------------------------'
              print*, 'begin program for solving sequence'
              print*, 'option: constant baryon mass'
              print*, 'central density range = ', e_min, ' to ', e_max
              print*, 'number of stars = ', n_of_models
              print*, '---------------------------------------------'

              ! convert units 
              e_min = e_min*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_min = e_min*c*c*k_scale        ! g/cm^3 to dimensionless
              e_max = e_max*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_max = e_max*c*c*k_scale        ! g/cm^3 to dimensionless

            write(*,"(a25)",advance='no') 'central density (MeV/fm^3),'
            write(*,"(a25)",advance='no') 'stellar mass (m_sun),'
            write(*,"(a25)",advance='no') 'baryonic mass (m_sun),'
            write(*,"(a25)",advance='no') 'equatorial radius (km),'
            write(*,"(a25)",advance='no') 'polar radius (km),'
            write(*,"(a25)",advance='no') 'frequency (hz),'
            write(*,"(a25)") 'kepler frequency (hz)'

              ! density increments
              a = (e_max/e_min)**(1./(n_of_models-1.)) 

              do j = 1, n_of_models

                e_center = (a**(1.*j-1.))*e_min 

                ! call model for baryon mass
                call m0_models

              end do

          !-----------------------------------------------------------------
          ! program for option 4 (spherical sequence)
          else if (option .eq. 4) then 

              ! read suboptions
              read(21,*) dummy
              read(21,*) e_min, e_max, n_of_models 

              ! screen output
              print*, '---------------------------------------------'
              print*, 'begin program for tabulated EOS'
              print*, 'option: spherical stars'
              print*, 'central density range = ', e_min, ' to ', e_max
              print*, 'number of stars = ', n_of_models
              print*, '---------------------------------------------'

              ! convert units
              e_min = e_min*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_min = e_min*c*c*k_scale        ! g/cm^3 to dimensionless
              e_max = e_max*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_max = e_max*c*c*k_scale        ! g/cm^3 to dimensionless

            write(*,"(a25)",advance='no') 'central density (MeV/fm^3),'
            write(*,"(a25)",advance='no') 'stellar mass (m_sun),'
            write(*,"(a25)",advance='no') 'baryonic mass (m_sun),'
            write(*,"(a25)",advance='no') 'equatorial radius (km),'
            write(*,"(a25)",advance='no') 'polar radius (km),'
            write(*,"(a25)",advance='no') 'frequency (hz),'
            write(*,"(a25)") 'kepler frequency (hz)'

              ! density increments
              a = (e_max/e_min)**(1./(n_of_models-1.))

              do j = 1, n_of_models

                e_center = (a**(1.*j-1.))*e_min 
                r_ratio = 1.e0

                ! spherical functions
                call center(e_center) 
                call TOV_guess
                call iterate
                call comp_phys
                call print_str_sph

              end do

          !-----------------------------------------------------------------
          ! program for option 5 (single star, kepler frequency sequence)
          else if (option .eq. 5) then

              ! read suboptions
              read(21,*) dummy 
              read(21,*) e_min, print_option 

              ! frequency tolerance 
              omega_error = err_tol 

              ! screen output 
              print*, '---------------------------------------------'
              print*, 'begin program for solving sequence'
              print*, 'option: single kepler frequency star'
              print*, 'central density = ', e_min, ' MeV/fm^3'
              print*, '---------------------------------------------'

              ! convert units
              e_min = e_min*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_min = e_min*c*c*k_scale        ! g/cm^3 to dimensionless

              e_center = e_min 

              ! call model for kepler frequency
              call kepler_models_single

          !-----------------------------------------------------------------
          ! program for option 6 (constant central density, variable
          ! frequency)
          else if (option .eq. 6) then

              ! read suboptions
              read(21,*) dummy 
              read(21,*) e_center, n_of_models, print_option 

              ! screen output 
              print*, '---------------------------------------------'
              print*, 'begin program for tabulated EOS'
              print*, 'option: constant central density'
              print*, 'central density = ', e_center, ' MeV/fm^3'
              print*, '---------------------------------------------'

              ! convert units
              e_center = e_center*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_center = e_center*c*c*k_scale        ! g/cm^3 to dimensionless

              ! calculate ratios 
              ratio_max = 1.
              ratio_min = 0.56
              dratio = (ratio_max-ratio_min)/(n_of_models-1.)
              ratio_loop = ratio_min 
              r_ratio = ratio_loop 

              do j = 1, n_of_models

                ! call single star functions
                call center(e_center)
                call TOV_guess
                call iterate
                call comp_phys
                call print_str_const_ec

                ratio_loop = ratio_loop + dratio
                r_ratio = ratio_loop 

              end do 

          !-----------------------------------------------------------------
          ! program for option 7 (spin up)
          else if (option .eq. 7) then

              ! read suboptions
              read(21,*) dummy
              read(21,*) n_of_models

              ! frequency tolerance
              m0_error = err_tol
              call spin_up(n_of_models)

          !-----------------------------------------------------------------
          ! program for option 8 (constant central density range)
          else if (option .eq. 8) then

              ! read suboptions
              read(21,*) dummy
              read(21,*) e_min, e_max, n_of_models, print_option

              ! screen output
              print*, '---------------------------------------------'
              print*, 'begin program for tabulated EOS'
              print*, 'option: constant central density range'
              print*, 'central density range = ', e_min, ' to ', e_max
              print*, 'number of stars = ', n_of_models
              print*, '---------------------------------------------'

              ! convert units
              e_min = e_min*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_min = e_min*c*c*k_scale        ! g/cm^3 to dimensionless
              e_max = e_max*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_max = e_max*c*c*k_scale        ! g/cm^3 to dimensionless

              ! density increments 
              a = (e_max/e_min)**(1./(n_of_models-1.))

              ! for all stars
              do i = 1, n_of_models

                e_center = (a**(1.*j-1.))*e_min

                ! calculate ratios
                ratio_max = 1.
                ratio_min = 0.56
                dratio = (ratio_max-ratio_min)/(n_of_models-1.)
                ratio_loop = ratio_min
                r_ratio = ratio_loop

                do j = 1, n_of_models

                    ! call single star functions
                    call center(e_center)
                    call TOV_guess
                    call iterate
                    call comp_phys
                    call print_str_const_ec

                    ratio_loop = ratio_loop + dratio
                    r_ratio = ratio_loop

                end do

              end do

          !-----------------------------------------------------------------
          ! program for option 9 (constant ratio)
          else if (option .eq. 9) then

              ! read suboptions 
              read(21,*) dummy 
              read(21,*)e_min, e_max, r_ratio, n_of_models, print_option

              ! screen output
              print*, '---------------------------------------------'
              print*, 'begin program for tabulated EOS'
              print*, 'option: constant ratio'
              print*, 'central density range = ', e_min, ' to ', e_max
              print*, 'number of stars = ', n_of_models
              print*, 'ratio r_p/r_e = ', r_ratio 
              print*, '---------------------------------------------'

              ! convert units 
              e_min = e_min*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_min = e_min*c*c*k_scale        ! g/cm^3 to dimensionless
              e_max = e_max*1.7827d12          ! MeV/fm^3 to g/cm^3
              e_max = e_max*c*c*k_scale        ! g/cm^3 to dimensionless

              ! density increments
              a = (e_max/e_min)**(1./(n_of_models-1.)) 

              ! for all stars
              do j = 1, n_of_models

                e_center = (a**(1.*j-1.))*e_min

                ! call single star functions
                call center(e_center)
                call TOV_guess
                call iterate
                call comp_phys
                call print_str_rat

              end do

          end if 

      end program

!--------------------------------------------------------------------------
! grid subroutine
!--------------------------------------------------------------------------

      subroutine make_grid

          use vars

          implicit none
          integer :: i

          do i = 1, s_div

            s_gp(i) = s_max*(1.*real(i)-1.d0)/(s_div-1.d0)
            s1_s(i) = s_gp(i)*(1.d0-s_gp(i))
            one_s(i) = 1.d0-s_gp(i)

          end do

          do i = 1, m_div

            mu(i) = (1.*real(i)-1.d0)/(m_div-1.d0) 
            one_m2(i) = 1.d0-mu(i)**2.d0
            theta(i) = acos(mu(i)) 
            sin_theta(i) = sqrt(one_m2(i)) 

          end do

      end subroutine

!--------------------------------------------------------------------------
! center subroutine
!--------------------------------------------------------------------------

      subroutine center(e_c) 

          use vars 

          implicit none
          double precision :: p_at_e, h_at_p, e_c 

          n_nearest = n_tab/2
          p_center = p_at_e(e_c) 
          h_center = h_at_p(p_center)

      end subroutine 

!--------------------------------------------------------------------------
! functions for center (p_at_e, h_at_p) 
!--------------------------------------------------------------------------

      double precision function p_at_e(ec)

          use vars

          implicit none
          double precision :: ec, argtemp, interp

          argtemp = interp(log_e_tab, log_p_tab, n_tab, log10(ec))
          p_at_e = 10.d0**(argtemp) 

      end function 

      double precision function h_at_p(pc) 

          use vars 

          implicit none
          double precision :: pc, argtemp, interp

          argtemp = interp(log_p_tab, log_h_tab, n_tab, log10(pc))
          h_at_p = 10.d0**(argtemp) 

      end function 

!--------------------------------------------------------------------------
! function rtsec_c
!--------------------------------------------------------------------------

      double precision function rtsec_c(x1, x2, xacc, e_c) 

          implicit none
          double precision :: x1, x2, xacc, e_c
          double precision :: fl, f, dx, swap, xl, rts
          double precision :: e_of_rho0 
          integer :: j, max_iter 

          max_iter = 100
          fl = e_of_rho0(x1)-e_c
          f = e_of_rho0(x2)-e_c 

          if (dabs(fl) .lt. dabs(f)) then 

              rts = x1
              xl = x2
              swap = fl
              fl = f
              f = swap 

          else 

              xl = x1
              rts = x2

          end if 

          do j = 1, max_iter

            dx = (xl-rts)*f/(f-fl) 
            xl = rts
            fl = f
            rts = rts+dx
            f = e_of_rho0(rts)-e_c

            if (dabs(dx) .lt. xacc .or. f .eq. 0.d0) then 
                
                rtsec_c = rts
                goto 10 

            end if 

          end do 

          print*, 'max # of iterations exceeded in rtsec_c'
          rtsec_c = 0.d0

10        continue 

      end function

!--------------------------------------------------------------------------
! function e_of_rho0 
!--------------------------------------------------------------------------

      double precision function e_of_rho0(x) 

          use vars 

          implicit none
          double precision :: x

          e_of_rho0 = (x**gamma_p)/(gamma_p-1.d0)+x

      end function

!--------------------------------------------------------------------------
! function e_at_p
!--------------------------------------------------------------------------

      double precision function e_at_p(pc) 

          use vars

          implicit none
          double precision :: pc, argtemp, interp

          argtemp = interp(log_p_tab, log_e_tab, n_tab, log10(pc))
          e_at_p = 10.d0**argtemp 

      end function

!--------------------------------------------------------------------------
! function n0_at_e
!--------------------------------------------------------------------------

      double precision function n0_at_e(ec)

          use vars

          implicit none
          double precision :: ec, argtemp, interp

          argtemp = interp(log_e_tab, log_n0_tab, n_tab,log10(ec))
          n0_at_e = 10.d0**argtemp

      end function

!--------------------------------------------------------------------------
! function p_at_h
!--------------------------------------------------------------------------

      double precision function p_at_h(h_c)

          use vars

          implicit none
          double precision :: h_c, argtemp, interp

          argtemp = interp(log_h_tab, log_p_tab, n_tab, log10(h_c))
          p_at_h = 10.d0**argtemp

      end function

!--------------------------------------------------------------------------
! function for interpolation (interp) 
!--------------------------------------------------------------------------

      double precision function interp(xp, yp, np, xb) 

          use vars

          implicit none
          integer :: k, m, np
          double precision :: y, xp(0:np), yp(0:np), xb

          m = 4

          call hunt(xp, np, xb, n_nearest) 

          k = min(max(n_nearest-(m-1)/2, 1), np+1-m) 

          if (xb .eq. xp(k) .or. xb .eq. xp(k+1) .or. xb .eq. xp(k+2) & 
              .or. xb .eq. xp(k+3)) then 

              xb = xb+1.d-12 

          end if

          y = (xb-xp(k+1))*(xb-xp(k+2))*(xb-xp(k+3))*yp(k)/ &
              ((xp(k)-xp(k+1))*(xp(k)-xp(k+2))*(xp(k)-xp(k+3))) & 
              +(xb-xp(k))*(xb-xp(k+2))*(xb-xp(k+3))*yp(k+1)/ &
              ((xp(k+1)-xp(k))*(xp(k+1)-xp(k+2))*(xp(k+1)-xp(k+3))) &
              +(xb-xp(k))*(xb-xp(k+1))*(xb-xp(k+3))*yp(k+2)/ & 
              ((xp(k+2)-xp(k))*(xp(k+2)-xp(k+1))*(xp(k+2)-xp(k+3))) &
              +(xb-xp(k))*(xb-xp(k+1))*(xb-xp(k+2))*yp(k+3)/ & 
              ((xp(k+3)-xp(k))*(xp(k+3)-xp(k+1))*(xp(k+3)-xp(k+2)))

          interp = y

      end function

!--------------------------------------------------------------------------
! hunt (nearest gridpoint) subroutine
!--------------------------------------------------------------------------

      subroutine hunt(xx, n, x, j_lo) 

          use vars

          implicit none 
          integer :: n, j_lo, j_hi, j_m, iter
          double precision :: x, xx(0:n) 
          logical :: ascnd 

          ascnd = xx(n) .ge. xx(1)

          if (j_lo .le. 0 .or. j_lo .gt. n) then 

              j_lo = 0
              j_hi = n+1
              goto 3

          end if

          iter = 1

          if (x .ge. xx(j_lo) .eqv. ascnd) then 

1             j_hi = j_lo+iter

              if (j_hi .gt. n) then 

                  j_hi = n+1

              elseif (x .ge. xx(j_hi) .eqv. ascnd) then 

                  j_lo = j_hi
                  iter = iter + iter
                  goto 1

              end if 

          else 

              j_hi = j_lo
2             j_lo = j_hi-iter

              if (j_lo .lt. 1) then 

                  j_lo = 0

              elseif (x .lt. xx(j_lo) .eqv. ascnd) then 

                  j_hi = j_lo
                  iter = iter + iter 
                  goto 2
                  
              end if 

          end if 

3         if (j_hi-j_lo .eq. 1) then 

              if (x .eq. xx(n)) then

                  j_lo = n-1

              end if

              if (x .eq. xx(1)) then

                  j_lo = 1

              end if

              return 

          end if

          j_m = (j_hi+j_lo)/2

          if (x .ge. xx(j_m) .eqv. ascnd) then 

              j_lo = j_m

          else 

              j_hi = j_m 

          end if 
          goto 3

      end subroutine 
      
!--------------------------------------------------------------------------
! print subroutine for spherically sym.
!--------------------------------------------------------------------------

      subroutine print_str_sph 

          use vars

          implicit none 
          integer :: i, j, ji 
          double precision :: r_dim, e_dim, mass0_print 
          double precision :: r_eprint, e_center_print, mass_print
          double precision :: omegaf_print, omegak_print 

          mass0_print = mass_0/m_sun
          r_eprint = r_ec/1.d5
          e_center_print = (e_center/(c*c*k_scale))/(1.7828d12)
          mass_print = mass/m_sun
          omegaf_print = (c/sqrt(kappa))*omega_f/(2.*pi)
          omegak_print = omega_k/(2.*pi) 

          write(*,1010) e_center_print, mass_print, mass0_print, &
              r_eprint, r_ratio*r_eprint, omegaf_print, omegak_print 

1010      format(1F15.5, xx, 1F15.5, xx, 1F15.5, xx, 1F15.5, xx, & 
              1F15.5, xx, 1F15.5, xx, 1F15.5)

      end subroutine

!--------------------------------------------------------------------------
! load EOS subroutine 
!--------------------------------------------------------------------------

      subroutine load_eos 

          use vars

          implicit none
          integer :: i, ialloc
          double precision :: pr, rhor, hr, n0r

          ! open file
          open(unit=91, file='eos.dat', status='unknown') 

          ! read total number of points
          read(91,*) n_tab

          ! allocate variables
          allocate(log_e_tab(1:n_tab+1), stat=ialloc)
          allocate(log_p_tab(1:n_tab+1), stat=ialloc)
          allocate(log_h_tab(1:n_tab+1), stat=ialloc)
          allocate(log_n0_tab(1:n_tab+1), stat=ialloc)

          ! read EOS - CGS system
          do i = 1, n_tab

            read(91,*) rhor, pr, hr, n0r

            log_e_tab(i) = log10(rhor*c*c*k_scale)    ! dimensionless
            log_p_tab(i) = log10(pr*k_scale)          ! dimensionless
            log_h_tab(i) = log10(hr/(c*c))           ! dimensionless
            log_n0_tab(i) = log10(n0r)               ! still w/ dims

          end do 

      end subroutine

!--------------------------------------------------------------------------
! print subroutine
!--------------------------------------------------------------------------

      subroutine print_str

          use vars

          implicit none
          integer :: i, j, ji 
          double precision :: r_dim, e_dim, mass0_print, r_eprint
          double precision :: e_center_print, mass_print, omegaf_print
          double precision :: omegak_print, mass_0_print, p_dim 

          open(unit = 11, file = 'output_struct.dat') 

          mass0_print = mass_0/m_sun
          r_eprint = r_ec/1.d5
          e_center_print = (e_center/(c*c*k_scale))/(1.7827d12) 
          mass_print = mass/m_sun
          omegaf_print = omega_f/(2.*pi) 
          omegak_print = omega_k/(2.*pi) 

          print*, 'central density (MeV/fm^3) = ', e_center_print
          print*, 'stellar mass (m_sun) = ', mass_print
          print*, 'baryonic mass (m_sun) = ', mass0_print
          print*, 'equatorial radius (km) = ', r_eprint
          print*, 'polar radius (km) = ', r_ratio*r_eprint
          print*, 'frequency (hz) = ', omegaf_print
          print*, 'kepler frequency (km) = ', omegak_print

          if (print_option .eq. 1) then 

              do j = 1, m_div, 1

                do i = 1, s_div/2+2

                    r_dim = (((r_ec*(s_gp(i)/(1.d0-s_gp(i))))))/1.d5
                    e_dim = ((energy(i,j))/(c*c*k_scale))/(1.7827d12)
                    p_dim = ((pressure(i,j))/(c*c*k_scale))/(1.7827d12)

                    write(11,1010) theta(j), r_dim, p_dim, alpha(i,j), &
                        gam(i,j), rho(i,j), omega(i,j), velocity_sq(i,j)

                end do

                write(11,*) ' '

              end do

          end if

1010      format(1F15.8,xxx,1E15.8,xxx,1F15.8,xxx,1F15.8,xxx,1F15.8, &
              xxx,1F15.8,xxx,1F15.8,xxx,1F15.8)

      end subroutine

!--------------------------------------------------------------------------
! print for constant ec
!--------------------------------------------------------------------------

      subroutine print_str_const_ec

          use vars

          implicit none
          integer :: i, j, ji
          double precision :: r_dim, e_dim, mass0_print, r_eprint
          double precision :: e_center_print, mass_print, omegaf_print
          double precision :: omegak_print, mass_0_print, p_dim

          open(unit=10, file='output_struct.dat') 

          mass0_print = mass_0/m_sun
          r_eprint = r_ec/1.d5
          e_center_print = (e_center/(c*c*k_scale))/(1.7827d12)
          mass_print = mass/m_sun 
          omegaf_print = omega_f/(2.*pi)
          omegak_print = omega_k/(2.*pi) 

          print*, e_center_print, mass_print, mass0_print, r_eprint, &
              r_ratio*r_eprint, omegaf_print, omegak_print 

      end subroutine

!--------------------------------------------------------------------------
! print for ratio
!--------------------------------------------------------------------------

      subroutine print_str_rat

          use vars

          implicit none
          integer :: i, j, ji
          double precision :: r_dim, e_dim, mass0_print, r_eprint
          double precision :: e_center_print, mass_print, omegaf_print
          double precision :: omegak_print, mass_0_print, p_dim

          mass0_print = mass_0/m_sun
          r_eprint = r_ec/1.d5
          e_center_print = (e_center/(c*c*k_scale))/(1.7827d12)
          mass_print = mass/m_sun
          omegaf_print = omega_f/(2.*pi)
          omegak_print = omega_k/(2.*pi)

          write(*,1013) e_center_print, mass_print, r_eprint, & 
              Im/(1.d45), r_p/r_e, omegaf_print

1013      format(1F25.15,xx,1F25.15,xx,1F25.15,xx,1F25.15,xx,1F25.15, &
              xx,1F25.15) 

      end subroutine

