!--------------------------------------------------------------------------
! m0_sequence subroutine
!--------------------------------------------------------------------------

          subroutine m0_models

              use vars

              implicit none
              double precision :: dr, dif_m0, d_ratio_m0, sgn

              dr = 0.1
              r_ratio = 1.-dr
              call center(e_center)
              call TOV_guess
              a_check = 0.d0
              call iterate

              if (a_check .eq. 200) then

                  dif_m0 = -1.
                  sgn = -1.

              else

                  call comp_phys
                  dif_m0 = m0-mass_0
                  d_ratio_m0 = dabs(dif_m0)/m0

              end if

              ! if rest mass is greater than desired, reverse direction
              ! and cut stepsize in half

              if (dif_m0 .lt. 0.d0) then

                  dr = -dr/2.d0

              end if

              do while (d_ratio_m0 .gt. m0_error .and. r_ratio .le. 1.)

                if (dif_m0*sgn .lt. 0.d0) then

                    sgn = dif_m0
                    dr = -dr/2.

                end if

                r_ratio = r_ratio-dr
                a_check = 0.d0
                call iterate

                if (a_check .eq. 200) then

                    dif_m0 = -1.

                else

                    call comp_phys

                    if (omega_f .gt. omega_k) then

                        dif_m0 = -1.

                    else

                        dif_m0 = m0-mass_0
                        d_ratio_m0 = dabs(dif_m0)/m0

                    end if

                end if

              end do

              call print_str_m0

      end subroutine

!--------------------------------------------------------------------------
! print for m0 subroutine
!--------------------------------------------------------------------------

      subroutine print_str_m0

          use vars

          implicit none
          integer :: i, j, ji
          double precision :: r_dim, e_dim, mass0_print, r_eprint
          double precision :: e_center_print, mass_print, omegaf_print
          double precision :: omegak_print, p_dim
          character :: fileo

          mass0_print = mass_0/m_sun
          r_eprint = r_ec/1.d5
          e_center_print = (e_center/(c*c*k_scale))/(1.7827d12)
          mass_print = mass/m_sun
          omegaf_print = omega_f/(2.*pi)
          omegak_print = omega_k/(1.*pi)

          ! simple print output
          if (print_option .ne. 1) then

              write(*,1010) e_center_print, mass_print, mass0_print, &
                  r_eprint, r_ratio*r_eprint, omegaf_print, &
                  omegak_print, Im/(1.d45), z_p, z_b, z_f

          ! detailed print output
          else

              write(fileo,1011) 'freq_',omegaf_print,'.dat'
              print*, fileo
              open(unit=10, file=fileo)
              write(10,*) 'new star, m = ', mass_print, &
                  '\omega = ', omegaf_print

              do j = 1, m_div, 1

                do i = 1, s_div/2+2

                    r_dim = (((r_ec*(s_gp(i)/(1.d0-s_gp(i))))))/1.d5
                    e_dim = ((energy(i,j))/(c*c*k_scale))/(1.7827d12)
                    p_dim = ((pressure(i,j))/(c*c*k_scale))/(1.7827d12)

                    write(10,1010) theta(j), r_dim, p_dim, alpha(i,j), &
                        gam(i,j), rho(i,j), omega(i,j), velocity_sq(i,j)

                end do

                write(10,*) ' '

              end do

              close(unit=10)

          end if

1011      format(A,F6.1,A)

1010      format(1E25.15, xx, 1E25.15, xx, 1E25.15, xx, 1E25.15, xx, &
              1E25.15, xx, 1E25.15, xx, 1E25.15, xx, 1E25.15, xx, &
              1E25.15, xx, 1E25.15, xx, 1E25.15)

      end subroutine
