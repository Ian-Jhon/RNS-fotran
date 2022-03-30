!--------------------------------------------------------------------------
! iteration subroutine
!--------------------------------------------------------------------------

      subroutine iterate() 

          use vars

          implicit none

          integer :: m, s, n, k, n_of_it, s_temp, kk

          double precision :: s_term, sum_rho, sum_gamma, sum_omega
          double precision :: r_e_old, dif, d_gamma_s, d_gamma_m
          double precision :: d_rho_s, d_rho_m, d_omega_s, d_omega_m
          double precision :: d_gamma_ss, d_gamma_mm, d_gamma_sm
          double precision :: temp1, temp2, temp3, temp4, temp5
          double precision :: temp6, temp7, temp8, m1, s1, s2
          double precision :: ea, r_emb_eq, r_emb_pole, rsm, gsm, esm
          double precision :: psm, v2sm, mum, omsm, sgp, s_1, e_gsm
          double precision :: e_rsm, rho0sm, sq, interp, deriv_s
          double precision :: deriv_m, deriv_sm, p_at_h, e_at_p

          ! initialize 
          n_of_it = 0
          s_term = 0.d0
          sum_rho = 0.d0
          sum_gamma = 0.d0
          sum_omega = 0.d0
          dif = 1.d0

          if (s_max .eq. 1.) then

              s_temp = s_div-1

          else 

              s_temp = s_div

          end if 

          ! setting variables to their guess value obtained from
          ! spherical solution
          do s = 1, s_div

            do m = 1, m_div

                gam(s,m) = gamma_guess(s,m)
                rho(s,m) = rho_guess(s,m)
                alpha(s,m) = alpha_guess(s,m) 
                omega(s,m) = omega_guess(s,m) 

            end do

          end do

          r_e = r_e_guess

          gamma_center_old = 0.d0
          rho_center_old = 0.d0
          print_dif = 0

          do while (dif .gt. accu .or. n_of_it .lt. 2)

            if (print_dif .ne. 0) print*, dif

            ! rescale potentials & make arrays along the polar and equ.
            ! directions
            do s = 1, s_div

                do m = 1, m_div

                    rho(s,m) = rho(s,m)/sq(r_e) 
                    gam(s,m) = gam(s,m)/sq(r_e) 
                    alpha(s,m) = alpha(s,m)/sq(r_e)
                    omega(s,m) = omega(s,m)*r_e

                end do

                ! at equator
                rho_mu_0(s) = rho(s,1)
                gamma_mu_0(s) = gam(s,1)
                omega_mu_0(s) = omega(s,1)

                ! at polar axis
                rho_mu_1(s) = rho(s,m_div) 
                gamma_mu_1(s) = gam(s,m_div) 

            end do

            ! update r_e
            r_e_old = r_e
            r_p = r_ratio*r_e
            s_p = r_p/(r_p+r_e)

            n_nearest = s_div/2
            gamma_pole_h = interp(s_gp, gamma_mu_1, s_div, s_p)
            gamma_equator_h = interp(s_gp, gamma_mu_0, s_div, s_e)
            gamma_center_h = gam(1,1) 
            rho_pole_h=interp(s_gp,rho_mu_1,s_div,s_p)
            rho_equator_h=interp(s_gp,rho_mu_0,s_div,s_e)
            rho_center_h=rho(1,1)

            r_e = sqrt(2*h_center/(gamma_pole_h+rho_pole_h- &
                gamma_center_h-rho_center_h))

            ! compute angular velocity (omega)
            if (r_ratio .eq. 1.) then

                omega_h = 0.d0
                omega_equator_h = 0.d0

            else 

                omega_equator_h = interp(s_gp, omega_mu_0, s_div, s_e)
                omega_h = omega_equator_h + exp(sq(r_e)*rho_equator_h)*&
                    sqrt(1.-exp(sq(r_e)*(gamma_pole_h+rho_pole_h- &
                    gamma_equator_h-rho_equator_h)))

                ! omega_h becomes a matrix for DR

            end if

            n_nearest = n_tab/2

            ! compute velocity, ed, and p
            do s = 1, s_div

                sgp = s_gp(s) 

                do m = 1, m_div

                    rsm = rho(s,m) 

                    if (r_ratio .eq. 1. .or. s .gt. s_div/2+2) then

                        velocity_sq(s,m) = 0.d0

                    else 

                        velocity_sq(s,m) = sq((omega_h-omega(s,m))* &
                            (sgp/(1.-sgp))*sin_theta(m)*exp(-rsm* & 
                            sq(r_e)))

                    end if

                    if (velocity_sq(s,m) .ge. 1.) then

                        velocity_sq(s,m) = 0.d0

                    end if

                    enthalpy(s,m) = enthalpy_min + 0.5*(sq(r_e)* &
                        (gamma_pole_h+rho_pole_h-gam(s,m)-rsm)- &
                        log(1.-velocity_sq(s,m)))

                    ! enthalphy adds a term for DR

                    if (enthalpy(s,m) .le. enthalpy_min .or. sgp & 
                        .gt. s_e) then

                        pressure(s,m) = 0.d0
                        energy(s,m) = 0.d0

                    else 

                        pressure(s,m) = p_at_h(enthalpy(s,m))
                        energy(s,m) = e_at_p(pressure(s,m))

                    end if

                    ! rescale back potentials - except omega
                    rho(s,m) = rho(s,m)*sq(r_e)
                    gam(s,m) = gam(s,m)*sq(r_e)
                    alpha(s,m) = alpha(s,m)*sq(r_e) 

                end do

            end do

            ! evaluation of source terms
            if (s_max .eq. 1) then

                s = s_div

                do m = 1, m_div

                    s_rho(s,m) = 0.d0
                    s_gamma(s,m) = 0.d0
                    s_omega(s,m) = 0.d0

                end do

            end if

            do s = 1, s_temp

                do m = 1, m_div

                    rsm = rho(s,m) 
                    gsm = gam(s,m) 
                    omsm = omega(s,m) 
                    esm = energy(s,m)
                    psm = pressure(s,m) 
                    e_gsm = exp(0.5*gsm) 
                    e_rsm = exp(-rsm) 
                    v2sm = velocity_sq(s,m)
                    mum = mu(m) 
                    m1 = 1.-sq(mum) 
                    sgp = s_gp(s) 
                    s_1 = 1.-sgp
                    s1 = sgp*s_1
                    s2 = sq(sgp/s_1) 
                    ea = 16.*pi*exp(2.*alpha(s,m))*sq(r_e) 

                    if (s .eq. 1) then

                        d_gamma_s = 0.
                        d_gamma_m = 0.
                        d_rho_s = 0.
                        d_rho_m = 0.
                        d_omega_s = 0.
                        d_omega_m = 0.

                    else 

                        d_gamma_s = deriv_s(gam,s,m)
                        d_gamma_m = deriv_m(gam,s,m) 
                        d_rho_s = deriv_s(rho,s,m)
                        d_rho_m = deriv_m(rho,s,m) 
                        d_omega_s = deriv_s(omega,s,m)
                        d_omega_m = deriv_m(omega,s,m) 

                    end if

                    s_rho(s,m) = e_gsm*(0.5*ea*(esm+psm)*s2*(1.+ & 
                        v2sm)/(1.-v2sm)+s2*m1*sq(e_rsm)*(sq(s1* & 
                        d_omega_s)+m1*sq(d_omega_m))+s1*d_gamma_s- &
                        mum*d_gamma_m+0.5*rsm*(ea*psm*s2-s1* &
                        d_gamma_s*(0.5*s1*d_gamma_s+1.)-d_gamma_m* &
                        (0.5*m1*d_gamma_m-mum))) 

                    s_gamma(s,m) = e_gsm*(ea*psm*s2+0.5*gsm*(ea* &
                        psm*s2-0.5*sq(s1*d_gamma_s)-0.5*m1* &
                        sq(d_gamma_m)))

                    s_omega(s,m) = e_gsm*e_rsm*(-ea*(omega_h-omsm)* &
                        (esm+psm)*s2/(1.-v2sm)+omsm*(-0.5*ea*(((1.+ & 
                        v2sm)*esm+2.*v2sm*psm)/(1.-v2sm))*s2-s1*(2.* &
                        d_rho_s+0.5*d_gamma_s)+mum*(2*d_rho_m+0.5* &
                        d_gamma_m)+0.25*sq(s1)*(4*sq(d_rho_s)- &
                        sq(d_gamma_s))+0.25*m1*(4*sq(d_rho_m)- &
                        sq(d_gamma_m))-m1*sq(e_rsm)*(sq(sq(sgp)* & 
                        d_omega_s)+s2*m1*sq(d_omega_m)))) 

                end do

            end do

            ! angular integration
            n = 0
            
            do k = 1, s_div

                do m = 1, m_div-2, 2

                    sum_rho = sum_rho + (dm/3.)*(p_2n(m,n)* &
                        s_rho(k,m)+4.*p_2n(m+1,n)*s_rho(k,m+1)+ &
                        p_2n(m+2,n)*s_rho(k,m+2))

                end do

                d1_rho(n,k) = sum_rho
                d1_gamma(n,k) = 0.
                d1_omega(n,k) = 0.
                sum_rho = 0. 

            end do

            do n = 1, l_max

                do k = 1, s_div

                    do m = 1, m_div-2, 2

                        sum_rho = sum_rho + (dm/3.)*(p_2n(m,n)* &
                            s_rho(k,m)+4.*p_2n(m+1,n)*s_rho(k,m+1)+ &
                            p_2n(m+2,n)*s_rho(k,m+2)) 

                        sum_gamma = sum_gamma + (dm/3.)*(sin((2.*n &
                            -1.)*theta(m))*s_gamma(k,m)+4.*sin((2.*n &
                            -1.)*theta(m+1))*s_gamma(k,m+1)+sin((2.*n &
                            -1.)*theta(m+2))*s_gamma(k,m+2)) 

                        sum_omega = sum_omega + (dm/3.)*(sin_theta(m)* &
                            p1_2n_1(m,n)*s_omega(k,m)+4.* &
                            sin_theta(m+1)*p1_2n_1(m+1,n)* &
                            s_omega(k,m+1)+sin_theta(m+2)* & 
                            p1_2n_1(m+2,n)*s_omega(k,m+2)) 

                    end do

                    d1_rho(n,k) = sum_rho
                    d1_gamma(n,k) = sum_gamma
                    d1_omega(n,k) = sum_omega

                    sum_rho = 0.
                    sum_gamma = 0.
                    sum_omega = 0. 

                end do

            end do

            ! radial integration 
            n = 0

            do s = 1, s_div

                do k = 1, s_div-2, 2

                    sum_rho = sum_rho + (ds/3.)*(f_rho(s,n,k)* &
                        d1_rho(n,k)+4.*f_rho(s,n,k+1)*d1_rho(n,k+1)+ &
                        f_rho(s,n,k+2)*d1_rho(n,k+2)) 

                end do

                d2_rho(s,n) = sum_rho
                d2_gamma(s,n) = 0.
                d2_omega(s,n) = 0.
                sum_rho = 0.

            end do

            do s = 1, s_div

                do n = 1, l_max

                    do k = 1, s_div-2, 2

                        sum_rho = sum_rho + (ds/3.)*(f_rho(s,n,k)* &
                            d1_rho(n,k)+4.*f_rho(s,n,k+1)* & 
                            d1_rho(n,k+1)+f_rho(s,n,k+2)*d1_rho(n,k+2)) 

                        sum_gamma = sum_gamma + (ds/3.)* &
                            (f_gamma(s,n,k)*d1_gamma(n,k)+4.* & 
                            f_gamma(s,n,k+1)*d1_gamma(n,k+1)+ &
                            f_gamma(s,n,k+2)*d1_gamma(n,k+2))

                        sum_omega = sum_omega + (ds/3.)* & 
                            (f_omega(s,n,k)*d1_omega(n,k)+4.* &
                            f_omega(s,n,k+1)*d1_omega(n,k+1)+ & 
                            f_omega(s,n,k+2)*d1_omega(n,k+2)) 

                    end do

                    d2_rho(s,n) = sum_rho
                    d2_gamma(s,n) = sum_gamma
                    d2_omega(s,n) = sum_omega 

                    sum_rho = 0.
                    sum_gamma = 0.
                    sum_omega = 0. 

                end do

            end do

            ! summing coefficients 
            do s = 1, s_div

                do m =1, m_div

                    gsm = gam(s,m) 
                    rsm = rho(s,m)
                    omsm = omega(s,m) 
                    e_gsm = exp(-0.5*gsm)
                    e_rsm = exp(rsm)
                    temp1 = sin_theta(m) 

                    sum_rho = sum_rho - e_gsm*p_2n(m,0)*d2_rho(s,0) 

                    do n = 1, l_max

                        sum_rho = sum_rho - e_gsm*p_2n(m,n)* &
                            d2_rho(s,n) 

                        if (m .eq. m_div) then

                            sum_omega = sum_omega + 0.5*e_rsm*e_gsm* &
                                d2_omega(s,n)

                            sum_gamma = sum_gamma - (2./pi)*e_gsm* &
                                d2_gamma(s,n) 

                        else 

                            sum_omega = sum_omega - e_rsm*e_gsm* &
                                (p1_2n_1(m,n)/(2.*n*(2.*n-1.)*temp1))* &
                                d2_omega(s,n)

                            sum_gamma = sum_gamma - (2./pi)*e_gsm* &
                                (sin((2.*n-1.)*theta(m))/((2.*n-1.)* &
                                temp1))*d2_gamma(s,n) 

                        end if

                    end do

                    rho(s,m) = rsm+cf*(sum_rho-rsm) 
                    gam(s,m) = gsm+cf*(sum_gamma-gsm) 
                    omega(s,m) = omsm+cf*(sum_omega-omsm) 

                    sum_rho = 0.
                    sum_gamma = 0.
                    sum_omega = 0.

                end do

            end do

            ! checking for divergence 
            if (dabs(omega(2,1)) .gt. 100.d0 .or. dabs(rho(2,1)) &
                .gt. 100.d0 .or. dabs(gam(2,1)) .gt. 300.d0) then 

                a_check = 200
                print*, 'b' 
                goto 119

            end if

            ! for spherical case
            if (r_ratio .eq. 1.) then 
                
                do s = 1, s_div

                    do m = 1, m_div

                        rho(s,m) = rho(s,1)
                        gam(s,m) = gam(s,1)
                        omega(s,m) = 0.d0

                    end do

                end do

            end if

            ! treating infinity when s_max = 1.
            if (s_max .eq. 1) then

                do m = 1, m_div

                    rho(s_div,m) = 0.d0
                    gam(s_div,m) = 0.d0
                    omega(s_div,m) = 0.d0

                end do

            end if

            ! computing first derivatives of gamma 
            do s = 1, s_div

                do m = 1, m_div

                    dgds(s,m) = deriv_s(gam,s,m)
                    dgdm(s,m) = deriv_m(gam,s,m) 

                end do

            end do

            ! alpha integration 
            if (r_ratio .eq. 1.) then

                do s = 1, s_div

                    do m = 1, m_div

                        dadm(s,m) = 0.d0

                    end do

                end do

            else 

                do s = 2, s_temp

                    do m = 1, m_div

                        dadm(1,m) = 0.d0
                        sgp = s_gp(s) 
                        s1 = sgp*(1.-sgp)
                        mum = mu(m)
                        m1 = 1.-sq(mum) 

                        d_gamma_s = dgds(s,m) 
                        d_gamma_m = dgdm(s,m)
                        d_rho_s = deriv_s(rho,s,m)
                        d_rho_m = deriv_m(rho,s,m) 
                        d_omega_s = deriv_s(omega,s,m)
                        d_omega_m = deriv_m(omega,s,m) 
                        d_gamma_ss = s1*deriv_s(dgds,s,m)+(1.-2.* &
                            sgp)*d_gamma_s
                        d_gamma_mm = m1*deriv_m(dgdm,s,m)-2.*mum* &
                            d_gamma_m
                        d_gamma_sm = deriv_sm(gam,s,m) 

                        temp1 = 2.*sq(sgp)*(sgp/(1.-sgp))*m1*d_omega_s &
                            *d_omega_m*(1.+s1*d_gamma_s)-(sq(sq(sgp)* &
                            d_omega_s)-sq(sgp*d_omega_m/(1.-sgp))*m1)* &
                            (-mum+m1*d_gamma_m) 

                        temp2 = 1./(m1*sq(1.+s1*d_gamma_s)+sq(-mum+m1 &
                            *d_gamma_m)) 

                        temp3 = s1*d_gamma_ss+sq(s1*d_gamma_s) 

                        temp4 = d_gamma_m*(-mum+m1*d_gamma_m) 

                        temp5 = (sq(s1*(d_rho_s+d_gamma_s))-m1* &
                            sq(d_rho_m+d_gamma_m))*(-mum+m1*d_gamma_m) 

                        temp6 = s1*m1*(0.5*(d_rho_s+d_gamma_s)* & 
                            (d_rho_m+d_gamma_m)+d_gamma_sm+d_gamma_s &
                            *d_gamma_m)*(1.+s1*d_gamma_s) 

                        temp7 = s1*mum*d_gamma_s*(1.+s1*d_gamma_s) 

                        temp8 = m1*exp(-2.*rho(s,m))

                        dadm(s,m) = -0.5*(d_rho_m+d_gamma_m)-temp2* &
                            (0.5*(temp3-d_gamma_mm-temp4)*(-mum+m1* & 
                            d_gamma_m)+0.25*temp5-temp6+temp7+0.25* &
                            temp8*temp1) 

                    end do

                end do

            end if

            do s = 1, s_temp

                alpha(s,1) = 0.d0

                do m = 1, m_div-1

                    alpha(s,m+1) = alpha(s,m) + 0.5*dm*(dadm(s,m+1)+ &
                        dadm(s,m))

                end do

            end do

            do s = 1, s_temp

                do m = 1, m_div

                    alpha(s,m) = alpha(s,m) - alpha(s,m_div)+0.5* &
                        (gam(s,m_div)-rho(s,m_div))

                    if (alpha(s,m) .ge. 300.) then

                        a_check = 200
                        print*, 'a'
                        goto 119

                    end if

                    omega(s,m) = omega(s,m)/r_e 

                end do

            end do

            if (s_max .eq. 1.) then

                do m = 1, m_div

                    alpha(s_div,m) = 0.d0

                end do

            end if

            if (a_check .eq. 200) goto 119

            dif = dabs(r_e_old-r_e)/r_e 
            n_of_it = n_of_it+1

          end do

          do s = 1, s_div

            do m = 1, m_div

                gamma_guess(s,m) = gam(s,m)
                rho_guess(s,m) = rho(s,m)
                omega_guess(s,m) = omega(s,m)
                alpha_guess(s,m) = alpha(s,m) 

            end do

          end do

          r_e_guess = r_e 

119       continue

      end subroutine

!--------------------------------------------------------------------------
! square function
!--------------------------------------------------------------------------

      double precision function sq(x) 

          implicit none
          double precision :: x

          sq = x*x

      end function

!--------------------------------------------------------------------------
! derivative function - s
!--------------------------------------------------------------------------

      double precision function deriv_s(f,s,m) 

          use vars 

          implicit none 
          integer :: s, m
          double precision :: f(0:s_div,0:m_div) 

          if (s .eq. 1) then 

              deriv_s = (f(s+1,m)-f(s,m))/ds

          else if (s .eq. s_div) then 

              deriv_s = (f(s,m)-f(s-1,m))/ds

          else 

              deriv_s = (f(s+1,m)-f(s-1,m))/(2.d0*ds)

          end if

      end function 

!--------------------------------------------------------------------------
! derivative function - m
!--------------------------------------------------------------------------

      double precision function deriv_m(f,s,m) 

          use vars

          implicit none
          integer :: s, m
          double precision :: f(0:s_div,0:m_div)

          if (m .eq. 1) then

              deriv_m = (f(s,m+1)-f(s,m))/dm

          else if (m .eq. m_div) then

              deriv_m = (f(s,m)-f(s,m-1))/dm

          else

              deriv_m = (f(s,m+1)-f(s,m-1))/(2.d0*dm)

          end if

      end function

!--------------------------------------------------------------------------
! derivative function - sm
!--------------------------------------------------------------------------

      double precision function deriv_sm(f,s,m)

          use vars

          implicit none 
          integer :: s, m
          double precision :: f(0:s_div,0:m_div) 

          if (s .eq. 1) then

              if (m .eq. 1) then 

                  deriv_sm =(f(s+1,m+1)-f(s,m+1)-f(s+1,m)+f(s,m))/ & 
                      (dm*ds)

              else if (m .eq. m_div) then

                  deriv_sm = (f(s+1,m)-f(s,m)-f(s+1,m-1)+f(s,m-1))/ &
                      (dm*ds)

              else 

                  deriv_sm = (f(s+1,m+1)-f(s+1,m-1)-f(s,m+1)+f(s,m-1)) &
                      /(2.*dm*ds) 

              end if

          else if (s .eq. s_div) then

              if (m .eq. 1) then 

                  deriv_sm = (f(s,m+1)-f(s,m)-f(s-1,m+1)+f(s-1,m))/ &
                      (dm*ds) 

              else if (m .eq. m_div) then 

                  deriv_sm = (f(s,m)-f(s-1,m)-f(s,m-1)+f(s-1,m-1))/ &
                      (dm*ds) 

              else 

                  deriv_sm = (f(s,m+1)-f(s,m-1)-f(s-1,m+1)+f(s-1,m-1)) &
                      /(2.*dm*ds) 

              end if

          else 

              if (m .eq. 1) then 

                  deriv_sm = (f(s+1,m+1)-f(s-1,m+1)-f(s+1,m)+f(s-1,m)) &
                      /(2.*dm*ds) 

              else if (m .eq. m_div) then 

                  deriv_sm = (f(s+1,m)-f(s-1,m)-f(s+1,m-1)+f(s-1,m-1)) &
                      /(2.*dm*ds) 

              else 

                  deriv_sm = (f(s+1,m+1)-f(s-1,m+1)-f(s+1,m-1)+ &
                      f(s-1,m-1))/(4.*dm*ds)

              end if

          end if

      end function


