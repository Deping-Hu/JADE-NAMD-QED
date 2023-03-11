! Convert electronic Hamiltion from adiabatic states to diabatic states 
!
!...cavity related subroutines
!
subroutine construct_Hcav(Hcav,n_state,n_atom,n_state_ele,pes_all_ele,mu_all,gc,wc,Hamitonian_Compoment,self_int_all)
implicit none
integer, intent(in)::n_state,n_state_ele,n_atom
double precision, intent(in) :: gc, wc
double precision, intent(in), dimension(n_state_ele) :: pes_all_ele
double precision, intent(in), dimension(n_state_ele,n_state_ele) :: mu_all
double precision, intent(in), dimension(n_state_ele,n_state_ele) :: self_int_all

double precision, dimension(n_state, n_state) :: Hcav
character(len=20)::Hamitonian_Compoment
double precision, external:: Kdelta

!!!! local varible
integer::i,j,ie,je,ip,jp,inp,jnp
real*8::Hp,Henp,Henp_JC

Hcav = 0.d0

do i=1,n_state
   do j=1,n_state
      call index_P2AF(i,n_state_ele,ie,ip)
      call index_P2AF(j,n_state_ele,je,jp)

      Hp=Kdelta(ie,je)*Kdelta(ip,jp)*wc*((ip-1)*1.d0) !remove the ZPE
      Henp=gc*mu_all(ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1)+sqrt(jp+0.d0)*Kdelta(ip,jp+1))

      if (ie>je) then
         Henp_JC=gc*mu_all(ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1))
      else if (ie<je) then
         Henp_JC=gc*mu_all(ie,je)*(sqrt(jp+0.d0)*Kdelta(ip,jp+1))
      else if (ie==je) then
         Henp_JC=0.d0
!         Henp_JC=gc*mu_all(ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1)+sqrt(jp+0.d0)*Kdelta(ip,jp+1))
      endif


      if(trim(adjustl(Hamitonian_Compoment))=="e")then
         Hcav(i,j)=pes_all_ele(ie)*Kdelta(ie,je)*Kdelta(ip,jp)
      else if(trim(adjustl(Hamitonian_Compoment))=="e+p")then
         Hcav(i,j)=pes_all_ele(ie)*Kdelta(ie,je)*Kdelta(ip,jp)+Hp
      else if(trim(adjustl(Hamitonian_Compoment))=="JC")then
         Hcav(i,j)=pes_all_ele(ie)*Kdelta(ie,je)*Kdelta(ip,jp)+Hp+Henp_JC
      else if(trim(adjustl(Hamitonian_Compoment))=="Rabi")then
         Hcav(i,j)=pes_all_ele(ie)*Kdelta(ie,je)*Kdelta(ip,jp)+Hp+Henp
      else if(trim(adjustl(Hamitonian_Compoment))=="PF")then
         Hcav(i,j)=pes_all_ele(ie)*Kdelta(ie,je)*Kdelta(ip,jp)+Hp+Henp+&
                Kdelta(ip,jp)*self_int_all(ie,je)*gc*gc/wc
      endif
  enddo
enddo

end subroutine construct_Hcav

!...calculate the gradient of Hcav: gHcav
subroutine  calc_gHcav (gHcav_x, &
                        gHcav_y, &
                        gHcav_z, &
                        pes_all_ele, &
                        Hcav, &
                        n_state, &
                        n_atom, &
                        n_state_ele, &
                        gra_all_ele_x, &
                        gra_all_ele_y, &
                        gra_all_ele_z, &
                        nac_ele_x, &
                        nac_ele_y, &
                        nac_ele_z, &
                        NACv_cav_x, &
                        NACv_cav_y, &
                        NACv_cav_z, &
                        g_mu_x, &
                        g_mu_y, &
                        g_mu_z, &
                        gc,   &
                        wc,   &
                        af_cut, &
                        Hamitonian_Compoment, &
                        g_self_int_x, &
                        g_self_int_y, &
                        g_self_int_z)

integer, intent(in)::n_state,n_state_ele,n_atom
double precision, intent(in) :: gc, wc, af_cut
double precision, intent(in), dimension(n_atom,n_state_ele,n_state_ele) :: g_mu_x
double precision, intent(in), dimension(n_atom,n_state_ele,n_state_ele) :: g_mu_y
double precision, intent(in), dimension(n_atom,n_state_ele,n_state_ele) :: g_mu_z
double precision, intent(in), dimension(n_state_ele,n_atom) :: gra_all_ele_x
double precision, intent(in), dimension(n_state_ele,n_atom) :: gra_all_ele_y
double precision, intent(in), dimension(n_state_ele,n_atom) :: gra_all_ele_z

double precision, intent(in), dimension(n_state_ele,n_state_ele,n_atom) :: nac_ele_x
double precision, intent(in), dimension(n_state_ele,n_state_ele,n_atom) :: nac_ele_y
double precision, intent(in), dimension(n_state_ele,n_state_ele,n_atom) :: nac_ele_z

double precision, intent(in), dimension(n_atom,n_state_ele,n_state_ele) :: g_self_int_x
double precision, intent(in), dimension(n_atom,n_state_ele,n_state_ele) :: g_self_int_y
double precision, intent(in), dimension(n_atom,n_state_ele,n_state_ele) :: g_self_int_z

double precision, intent(in), dimension(n_state_ele) :: pes_all_ele
double precision, intent(in), dimension(n_state, n_state) :: Hcav

character(len=20), intent(in)::Hamitonian_Compoment

double precision, dimension(n_state,n_state,n_atom) :: gHcav_x
double precision, dimension(n_state,n_state,n_atom) :: gHcav_y
double precision, dimension(n_state,n_state,n_atom) :: gHcav_z

double precision, dimension(n_state,n_state,n_atom) :: NACv_cav_x
double precision, dimension(n_state,n_state,n_atom) :: NACv_cav_y
double precision, dimension(n_state,n_state,n_atom) :: NACv_cav_z

double precision, external:: Kdelta

!!!! local varible
integer::i,j,ie,je,ip,jp,idp,i_atom
real*8::dHenp_x,dHenp_y,dHenp_z, dHenp_x_JC, dHenp_y_JC, dHenp_z_JC
real*8::XH(n_state,n_state),HX(n_state,n_state), new_force(n_state,n_state), old_force(n_state,n_state)

gHcav_x=0.d0
gHcav_y=0.d0
gHcav_z=0.d0
NACv_cav_x=0.d0
NACv_cav_y=0.d0
NACv_cav_z=0.d0

do i=1,n_state
   do j=1,n_state
      call index_P2AF(i,n_state_ele,ie,ip)
      call index_P2AF(j,n_state_ele,je,jp)

     do i_atom =1,n_atom

      dHenp_x=gc*g_mu_x(i_atom,ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1)+sqrt(jp+0.d0)*Kdelta(ip,jp+1))
      dHenp_y=gc*g_mu_y(i_atom,ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1)+sqrt(jp+0.d0)*Kdelta(ip,jp+1))
      dHenp_z=gc*g_mu_z(i_atom,ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1)+sqrt(jp+0.d0)*Kdelta(ip,jp+1))


!!!! To be used in JC model
      if (ie>je) then
         dHenp_x_JC=gc*g_mu_x(i_atom,ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1))
         dHenp_y_JC=gc*g_mu_y(i_atom,ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1))
         dHenp_z_JC=gc*g_mu_z(i_atom,ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1))
      else if (ie<je) then
         dHenp_x_JC=gc*g_mu_x(i_atom,ie,je)*(sqrt(jp+0.d0)*Kdelta(ip,jp+1))
         dHenp_y_JC=gc*g_mu_y(i_atom,ie,je)*(sqrt(jp+0.d0)*Kdelta(ip,jp+1))
         dHenp_z_JC=gc*g_mu_z(i_atom,ie,je)*(sqrt(jp+0.d0)*Kdelta(ip,jp+1))
      else if (ie==je) then
!         dHenp_x_JC=gc*g_mu_x(i_atom,ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1)+sqrt(jp+0.d0)*Kdelta(ip,jp+1))
!         dHenp_y_JC=gc*g_mu_y(i_atom,ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1)+sqrt(jp+0.d0)*Kdelta(ip,jp+1))
!         dHenp_z_JC=gc*g_mu_z(i_atom,ie,je)*(sqrt(jp-1.d0)*Kdelta(ip,jp-1)+sqrt(jp+0.d0)*Kdelta(ip,jp+1))

         dHenp_x_JC=0.d0
         dHenp_y_JC=0.d0
         dHenp_z_JC=0.d0
      endif


      if(trim(adjustl(Hamitonian_Compoment))=="e")then
         gHcav_x(i,j,i_atom)=gra_all_ele_x(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)
         gHcav_y(i,j,i_atom)=gra_all_ele_y(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)
         gHcav_z(i,j,i_atom)=gra_all_ele_z(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)
      else if(trim(adjustl(Hamitonian_Compoment))=="e+p")then
         gHcav_x(i,j,i_atom)=gra_all_ele_x(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)
         gHcav_y(i,j,i_atom)=gra_all_ele_y(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)
         gHcav_z(i,j,i_atom)=gra_all_ele_z(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)

      else if(trim(adjustl(Hamitonian_Compoment))=="JC")then
         gHcav_x(i,j,i_atom)=gra_all_ele_x(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)+dHenp_x_JC
         gHcav_y(i,j,i_atom)=gra_all_ele_y(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)+dHenp_y_JC
         gHcav_z(i,j,i_atom)=gra_all_ele_z(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)+dHenp_z_JC

      else if(trim(adjustl(Hamitonian_Compoment))=="Rabi")then
         gHcav_x(i,j,i_atom)=gra_all_ele_x(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)+dHenp_x
         gHcav_y(i,j,i_atom)=gra_all_ele_y(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)+dHenp_y
         gHcav_z(i,j,i_atom)=gra_all_ele_z(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)+dHenp_z
      else if(trim(adjustl(Hamitonian_Compoment))=="PF")then
         gHcav_x(i,j,i_atom)=gra_all_ele_x(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)+dHenp_x+g_self_int_x(i_atom,ie,je)*Kdelta(ip,jp)*gc*gc/wc
         gHcav_y(i,j,i_atom)=gra_all_ele_y(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)+dHenp_y+g_self_int_y(i_atom,ie,je)*Kdelta(ip,jp)*gc*gc/wc
         gHcav_z(i,j,i_atom)=gra_all_ele_z(ie,i_atom)*Kdelta(ie,je)*Kdelta(ip,jp)+dHenp_z+g_self_int_z(i_atom,ie,je)*Kdelta(ip,jp)*gc*gc/wc
      endif
      !derivative coupling matrix
      NACv_cav_x(i,j,i_atom)=nac_ele_x(ie,je,i_atom)*Kdelta(ip,jp)
      NACv_cav_y(i,j,i_atom)=nac_ele_y(ie,je,i_atom)*Kdelta(ip,jp)
      NACv_cav_z(i,j,i_atom)=nac_ele_z(ie,je,i_atom)*Kdelta(ip,jp)

     enddo
   enddo
enddo

!if ((pes_all_ele(2)-pes_all_ele(1)) .lt. af_cut) then

!   write(*,*) "calculate the NAC in PL basis"

   do i_atom=1,n_atom
   
      XH=matmul(NACv_cav_x(:,:,i_atom),Hcav)   !XH
      HX=matmul(Hcav,NACv_cav_x(:,:,i_atom))   !HX
      new_force=XH-HX
      gHcav_x(:,:,i_atom)=gHcav_x(:,:,i_atom)+new_force
   
      XH=matmul(NACv_cav_y(:,:,i_atom),Hcav)   !XH
      HX=matmul(Hcav,NACv_cav_y(:,:,i_atom))   !HX
      new_force=XH-HX
      gHcav_y(:,:,i_atom)=gHcav_y(:,:,i_atom)+new_force
   
      XH=matmul(NACv_cav_z(:,:,i_atom),Hcav)   !XH
      HX=matmul(Hcav,NACv_cav_z(:,:,i_atom))   !HX
      new_force=XH-HX
      gHcav_z(:,:,i_atom)=gHcav_z(:,:,i_atom)+new_force
   
   enddo

!endif


end subroutine calc_gHcav
!
!!
!...change the index from dressed basis to Adiabatic-Fock basis
subroutine index_P2AF(i,ns,is,ip)
implicit none
integer,intent(in)::i   !dressed state index
integer,intent(in)::ns  !# of electronic states
integer,intent(out)::is,ip  !corresponding electronic and photon states

is=mod(i,ns)
if(is==0)is=ns
if(mod(i,ns)/=0)then
  ip=i/ns+1
else
  ip=i/ns
endif

end subroutine index_P2AF
!
!...calculate the non-adiabatic coupling between polariton
!...rotate the gHcav and divide by energy difference

subroutine  calc_pola_NACv(n_state, &
                           n_atom, &
                           pes_all, &
                           gHcav_x,&
                           gHcav_y,&
                           gHcav_z,&
                           gra_all_x,&
                           gra_all_y,&
                           gra_all_z,&
                           nac_x, &
                           nac_y, &
                           nac_z, &
                           U_ss)
implicit none
integer,intent(in)::n_state,n_atom
double precision,  dimension(n_state, n_state) :: U_ss
double precision, intent(in), dimension(n_state) :: pes_all
double precision, intent(in), dimension(n_state, n_state,n_atom) :: gHcav_x
double precision, intent(in), dimension(n_state, n_state,n_atom) :: gHcav_y
double precision, intent(in), dimension(n_state, n_state,n_atom) :: gHcav_z

double precision, dimension(n_state,n_atom) :: gra_all_x
double precision, dimension(n_state,n_atom) :: gra_all_y
double precision, dimension(n_state,n_atom) :: gra_all_z


double precision, dimension(n_state, n_state,n_atom) :: g_pola_H_x
double precision, dimension(n_state, n_state,n_atom) :: g_pola_H_y
double precision, dimension(n_state, n_state,n_atom) :: g_pola_H_z

double precision, dimension(n_state, n_state,n_atom) :: nac_x
double precision, dimension(n_state, n_state,n_atom) :: nac_y
double precision, dimension(n_state, n_state,n_atom) :: nac_z

integer::is,js,i_atom

nac_x=0.d0
nac_y=0.d0
nac_z=0.d0

do i_atom=1,n_atom
g_pola_H_x(:,:,i_atom)=gHcav_x(:,:,i_atom)
g_pola_H_y(:,:,i_atom)=gHcav_y(:,:,i_atom)
g_pola_H_z(:,:,i_atom)=gHcav_z(:,:,i_atom)

call Simil_Trans_real("UtAU",n_state,U_ss,g_pola_H_x(:,:,i_atom))
call Simil_Trans_real("UtAU",n_state,U_ss,g_pola_H_y(:,:,i_atom))
call Simil_Trans_real("UtAU",n_state,U_ss,g_pola_H_z(:,:,i_atom))

!calculate the non-adiabatic coupling
do is=1,n_state-1
   do js=is+1,n_state
      nac_x(is,js,i_atom)=-g_pola_H_x(is,js,i_atom)/(pes_all(is)-pes_all(js))
      nac_x(js,is,i_atom)=-nac_x(is,js,i_atom)

      nac_y(is,js,i_atom)=-g_pola_H_y(is,js,i_atom)/(pes_all(is)-pes_all(js))
      nac_y(js,is,i_atom)=-nac_y(is,js,i_atom)

      nac_z(is,js,i_atom)=-g_pola_H_z(is,js,i_atom)/(pes_all(is)-pes_all(js))
      nac_z(js,is,i_atom)=-nac_z(is,js,i_atom)
   enddo
enddo

do is=1,n_state
  gra_all_x(is,i_atom) = g_pola_H_x(is,is,i_atom)
  gra_all_y(is,i_atom) = g_pola_H_y(is,is,i_atom)
  gra_all_z(is,i_atom) = g_pola_H_z(is,is,i_atom)
enddo

enddo

end subroutine calc_pola_NACv

