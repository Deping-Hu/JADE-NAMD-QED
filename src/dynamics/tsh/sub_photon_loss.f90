      subroutine sub_photon_loss(n_state, &
                              n_state_ele, &
                              rho_af, &
                              photon_loss_coef)

      implicit none

!
!     ---- Parameters -----
!
      include 'param.def'
!
!     ---- Argument -----

      integer, intent(in) :: n_state, n_state_ele

      complex(kind=8), intent(inout), dimension(n_state, n_state) :: rho_af

!     local variables

      integer :: i, j, k, ie, ip, je, jp, i_state, j_state

      double precision, dimension(n_state, n_state) :: Jump_operator, Jump_operator_trans

      complex(kind=8), dimension(n_state, n_state) :: rho_tmp1, rho_tmp2, rho_tmp3

      complex(kind=8), dimension(n_state, n_state) :: delta_rho, photon_loss_coef

      photon_loss_coef = 0.0

      do i_state =1, n_state

         do j_state =1, n_state

            call index_P2AF(i_state,n_state_ele,ie,ip)
            call index_P2AF(j_state,n_state_ele,je,jp)

            Jump_operator = 0.d0

            if ((ie .eq. je) .and. ((jp-ip) .eq. 1)) then
         
               Jump_operator(i_state,j_state) = 1

               Jump_operator_trans = transpose(Jump_operator)

               rho_tmp1 = matmul(Jump_operator,rho_af)
               rho_tmp1 = matmul(rho_tmp1, Jump_operator_trans)

               rho_tmp2 = matmul(Jump_operator_trans,Jump_operator)
               rho_tmp2 = matmul(rho_tmp2, rho_af)

               rho_tmp3 = matmul(Jump_operator_trans,Jump_operator)
               rho_tmp3 = matmul(rho_af,rho_tmp3)

               delta_rho =  rho_tmp1 - 0.5*(rho_tmp2+rho_tmp3)

               delta_rho = delta_rho*ip

               photon_loss_coef = photon_loss_coef + delta_rho

            end if

         enddo

      enddo

      return

      end subroutine sub_photon_loss
