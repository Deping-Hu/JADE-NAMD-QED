    subroutine sub_tully_hopping_ana_nac(n_atom, &
                                         mass, &
                                         old_vel_x, &
                                         old_vel_y, &
                                         old_vel_z, &
                                         vel_x, &
                                         vel_y, &
                                         vel_z, &
                                         n_state, &
                                         old_pes_all, &
                                         pes_all, &
                                         pes_ref, &
                                         old_nac_x, &
                                         old_nac_y, &
                                         old_nac_z, &
                                         nac_x, &
                                         nac_y, &
                                         nac_z, &
                                         rho, &
                                         time, &
                                         it, &
                                         dtime, &
                                         ntime_ele, &
                                         index_state, &
                                         cor_dec, &
                                         file_save_ele, &
                                         file_md_out, &
                                         label_hop_prob)

       implicit none

!
!     ----- Parameters -----
!
       include 'param.def'
!
!     ----- Argument -----

       integer, intent(in) :: label_hop_prob, file_save_ele, file_md_out, &

                              n_atom, n_state, ntime_ele, it

       integer, intent(inout) :: index_state

       double precision, intent(in) :: dtime, pes_ref, time, cor_dec

       double precision, intent(in), dimension(n_state) ::  pes_all, old_pes_all

       double precision, intent(in), dimension(n_atom) :: mass, vel_x, vel_y, vel_z, &
                                                          old_vel_x, old_vel_y, old_vel_z

       double precision, intent(in), dimension(n_state, n_state, n_atom) :: nac_x, nac_y, nac_z, &
                                                                            old_nac_x, old_nac_y, old_nac_z

       complex(kind=8), intent(inout), dimension(n_state, n_state) :: rho

!
!  Local

       integer :: i, j, k, l, it_ele, lwork, info, label_rk4

       double precision :: delt, num_random, rho_decress_all

       double precision, dimension(n_state) :: pes_all_eff, ev, pro_hop, pro_sum

       double precision, dimension(n_atom) :: vel_x_eff, vel_y_eff, vel_z_eff

       double precision, dimension(n_state, n_state) :: Hadia_eff

       double precision, dimension(n_state, n_state, n_atom) :: nac_x_eff, nac_y_eff, nac_z_eff

       double precision, allocatable, dimension(:) :: work

       complex(kind=8), dimension(n_state, n_state) :: h_vv, h_vv_conjg, rho_tmp, &
                                                       old_rho, tmp1, tmp2, tmp3, tmp4, &
                                                       ematrix, k1, k2, k3, k4, rho_k_tmp

       complex(kind=8), allocatable, dimension(:, :) :: rwork

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       lwork = n_state*(n_state + 1)

       delt = dtime/ntime_ele

       label_rk4 = 1

       allocate (work(lwork))
       allocate (rwork(3*n_state - 2, 3*n_state - 2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

       old_rho = rho

       pro_hop = 0.d0

       do it_ele = 1, ntime_ele

          pes_all_eff = old_pes_all + &
                        (it_ele - 1)*(pes_all - old_pes_all)/ntime_ele
          nac_x_eff = old_nac_x + &
                      (it_ele - 1)*(nac_x - old_nac_x)/ntime_ele
          nac_y_eff = old_nac_y + &
                      (it_ele - 1)*(nac_y - old_nac_y)/ntime_ele
          nac_z_eff = old_nac_z + &
                      (it_ele - 1)*(nac_z - old_nac_z)/ntime_ele
          vel_x_eff = old_vel_x + &
                      (it_ele - 1)*(vel_x - old_vel_x)/ntime_ele
          vel_y_eff = old_vel_y + &
                      (it_ele - 1)*(vel_y - old_vel_y)/ntime_ele
          vel_z_eff = old_vel_z + &
                      (it_ele - 1)*(vel_z - old_vel_z)/ntime_ele

!!        Construct the effective Hamiltonian
          if (label_rk4 .ne. 1) then

             h_vv(:, :) = cmplx(0.d0, 0.d0)
             do i = 1, n_state
             do j = 1, n_state
                if (i .eq. j) then
                   h_vv(i, j) = pes_all_eff(i)
                   h_vv(i, j) = h_vv(i, j) - pes_ref
                end if
             end do
             end do

             do i = 1, n_state
             do j = 1, n_state
                if (i .ne. j) then
                   do k = 1, n_atom
                      h_vv(i, j) = h_vv(i, j) &
                                   + nac_x_eff(i, j, k)*vel_x_eff(k) &
                                   + nac_y_eff(i, j, k)*vel_y_eff(k) &
                                   + nac_z_eff(i, j, k)*vel_z_eff(k)
                   end do
                   h_vv(i, j) = -c1*h_vv(i, j)
                end if
             end do
             end do

             call ZHEEV('V', 'U', n_state, h_vv, n_state, &
                        ev, work, lwork, rwork, info)

             do k = 1, n_state
             do l = 1, n_state
                ematrix(k, l) = 0.d0
             end do
             end do

             do k = 1, n_state
                ematrix(k, k) = cdexp(-c1*ev(k)*delt)
             end do

             tmp1 = MATMUL(h_vv, ematrix)

             do i = 1, n_state
             do j = 1, n_state
                h_vv_conjg(i, j) = conjg(h_vv(j, i))
             end do
             end do

             tmp2 = MATMUL(tmp1, h_vv_conjg)

             rho_tmp = MATMUL(tmp2, rho)

             do k = 1, n_state
                ematrix(k, k) = cdexp(c1*ev(k)*delt)
             end do

             tmp1 = MATMUL(h_vv, ematrix)
             tmp2 = MATMUL(tmp1, h_vv_conjg)

             rho = MATMUL(rho_tmp, tmp2)

          else if (label_rk4 .eq. 1) then

             Hadia_eff = 0.d0
             do i = 1, n_state
             do j = 1, n_state
                if (i .eq. j) then
                   Hadia_eff(i, j) = pes_all_eff(i)
                end if
             end do
             end do

             h_vv(:, :) = cmplx(0.d0, 0.d0)

             do i = 1, n_state
             do j = 1, n_state
                do k = 1, n_atom
                   h_vv(i, j) = h_vv(i, j) &
                                + nac_x_eff(i, j, k)*vel_x_eff(k) &
                                + nac_y_eff(i, j, k)*vel_y_eff(k) &
                                + nac_z_eff(i, j, k)*vel_z_eff(k)
                end do
             end do
             end do

!!!!!!! 4th Runge-Kutta to propagate the density matrix

!!!!!!!!!!!!! k1
             rho_k_tmp = rho
             tmp1 = MATMUL(Hadia_eff, rho_k_tmp)
             tmp2 = MATMUL(rho_k_tmp, Hadia_eff)

             tmp3 = MATMUL(h_vv, rho_k_tmp)
             tmp4 = MATMUL(rho_k_tmp, h_vv)

             k1 = -c1*(tmp1-tmp2)-(tmp3-tmp4) 

             rho_k_tmp = rho + k1*delt/2.d0

!!!!!!!!!!!!! k2
             tmp1 = MATMUL(Hadia_eff, rho_k_tmp)
             tmp2 = MATMUL(rho_k_tmp, Hadia_eff)

             tmp3 = MATMUL(h_vv, rho_k_tmp)
             tmp4 = MATMUL(rho_k_tmp, h_vv)

             k2 = -c1*(tmp1-tmp2)-(tmp3-tmp4) 

             rho_k_tmp = rho + k2*delt/2.d0

!!!!!!!!!!!!! k3
             tmp1 = MATMUL(Hadia_eff, rho_k_tmp)
             tmp2 = MATMUL(rho_k_tmp, Hadia_eff)

             tmp3 = MATMUL(h_vv, rho_k_tmp)
             tmp4 = MATMUL(rho_k_tmp, h_vv)

             k3 = -c1*(tmp1-tmp2)-(tmp3-tmp4) 

             rho_k_tmp = rho + k3*delt

!!!!!!!!!!!!! k4
             tmp1 = MATMUL(Hadia_eff, rho_k_tmp)
             tmp2 = MATMUL(rho_k_tmp, Hadia_eff)

             tmp3 = MATMUL(h_vv, rho_k_tmp)
             tmp4 = MATMUL(rho_k_tmp, h_vv)

             k4 = -c1*(tmp1-tmp2)-(tmp3-tmp4) 

!!!!! Sum all
             rho = rho +(1/6.d0)*delt*(k1+2*k2+2*k3+k4)
!!!!! End RK

          endif

!            Decoherence corrections during the eletronic propagation
!          if (cor_dec .gt. 0) then
!             call sub_decoherent_corr(n_atom, &
!                                      mass, &
!                                      vel_x, &
!                                      vel_y, &
!                                      vel_z, &
!                                      n_state, &
!                                      pes_all, &
!                                      rho, &
!                                      delt, &
!                                      index_state, &
!                                      cor_dec, &
!                                      file_md_out)
!
!          end if

!!!! Convential Surface Hopping

          if (label_hop_prob .eq. 1) then

             do i = 1, n_state
             do j = 1, n_state
                if (i .eq. j) then
                   h_vv(i, j) = pes_all_eff(i)
                   h_vv(i, i) = h_vv(i, j) - pes_ref
                end if
             end do
             end do

             do i = 1, n_state
             do j = 1, n_state
                h_vv(i, j) = 0.d0
                do k = 1, n_atom
                   h_vv(i, j) = h_vv(i, j) &
                                + nac_x_eff(i, j, k)*vel_x(k) &
                                + nac_y_eff(i, j, k)*vel_y(k) &
                                + nac_z_eff(i, j, k)*vel_z(k)
                end do
                h_vv(i, j) = -c1*h_vv(i, j)
             end do
             end do

             do i = 1, n_state
                if (index_state .ne. i) then
                   h_vv(i, index_state) = 0.d0
                   do k = 1, n_atom
                      h_vv(i, index_state) = h_vv(i, index_state) &
                                             + nac_x_eff(i, index_state, k)*vel_x(k) &
                                             + nac_y_eff(i, index_state, k)*vel_y(k) &
                                             + nac_z_eff(i, index_state, k)*vel_z(k)
                   end do
                   pro_hop(i) = pro_hop(i) &
                                - (2*delt* &
                                   real(conjg(rho(i, index_state))) &
                                   *h_vv(i, index_state) &
                                   )/real(rho(index_state, index_state))
                end if
             end do

          end if

       end do

!!!! Global Flux Surface Hopping (GFSH)

       if (label_hop_prob .eq. 3) then

          if (real(rho(index_state, index_state)) .lt. real(old_rho(index_state, index_state))) then

             rho_decress_all = 0.d0

             do i = 1, n_state
                if (real(rho(i, i)) .lt. real(old_rho(i, i))) then
                   rho_decress_all = rho_decress_all + (real(old_rho(i, i)) - real(rho(i, i)))
                end if
             end do

             do i = 1, n_state

                pro_hop(i) = (real(rho(i, i)) - real(old_rho(i, i))) &
                             /real(old_rho(index_state, index_state)) &
                             *(real(old_rho(index_state, index_state)) - real(rho(index_state, index_state))) &
                             /rho_decress_all

             end do

          end if

       end if


!      Decoherence corrections after the eletronic propagation
       if (cor_dec .gt. 0) then
          call sub_decoherent_corr(n_atom, &
                                  mass, &
                                  vel_x, &
                                  vel_y, &
                                  vel_z, &
                                  n_state, &
                                  pes_all, &
                                  rho, &
                                  dtime, &
                                  index_state, &
                                  cor_dec, &
                                  file_md_out)
       endif
       

       do i = 1, n_state
          if (pro_hop(i) .lt. 0.d0) then
             pro_hop(i) = 0.d0
          end if
       end do

       pro_hop(index_state) = 1
       do i = 1, n_state
          if (index_state .ne. i) then
             pro_hop(index_state) = pro_hop(index_state) - pro_hop(i)
          end if
       end do

       pro_sum = 0.d0
       do i = 1, n_state
          do j = 1, i
             pro_sum(i) = pro_sum(i) + pro_hop(j)
          end do
       end do

       write (file_save_ele, *) "------------------------------------"

       do i = 1, n_state
       do j = 1, n_state
          write (file_save_ele, 9999) it, time*TOFS, &
             "rho", i, j, &
             REAL(rho(i, j)), '+i', AIMAG(rho(i, j))
       end do
       end do

       write (file_save_ele, 9996) it, time*TOFS, &
          "The current state", index_state
       write (file_save_ele, 9997) it, time*TOFS, &
          "Hopping probability", pro_hop(:)
       write (file_save_ele, 9997) it, time*TOFS, &
          "Area for hopping", pro_sum(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       generate random number
!       the init_random_seed may be set before the dynamic cycle.
       call get_random_number(num_random)
       call get_random_number(num_random)
       call get_random_number(num_random)
       call get_random_number(num_random)
       call get_random_number(num_random)

!       Surface-hopping check

       if ((0 .le. num_random) &
           .and. &
           (num_random .lt. pro_sum(1)) &
           ) then
          index_state = 1
          goto 7000
       else
          do i = 2, n_state
             if ((pro_sum(i - 1) .le. num_random) &
                 .and. &
                 (num_random .lt. pro_sum(i)) &
                 ) then
                index_state = i
                goto 7000
             end if
          end do
       end if

7000   write (file_save_ele, 9997) it, time*TOFS, &
          "Random number", num_random
       write (file_save_ele, 9996) it, time*TOFS, &
          "The new state", index_state
       write (file_save_ele, *) "------------------------------------"

9996   format(i10, 1x, f15.8, 3x, a, 1x, i3)
9997   format(i10, 1x, f15.8, 3x, a, 1x, 10(f15.8, 1x))
9999   format(i10, 1x, f15.8, 3x, a, 1x, i3, 1x, i3, 1x, f10.7, a, f10.7)

       deallocate (work)
       deallocate (rwork)

    end
