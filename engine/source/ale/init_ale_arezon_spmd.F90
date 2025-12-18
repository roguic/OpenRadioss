!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2025 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        the Free Software Foundation, either version 3 of the License, or
!Copyright>        (at your option) any later version.
!Copyright>
!Copyright>        This program is distributed in the hope that it will be useful,
!Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!Copyright>        GNU Affero General Public License for more details.
!Copyright>
!Copyright>        You should have received a copy of the GNU Affero General Public License
!Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!Copyright>
!Copyright>
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      module init_ale_arezon_spmd_mod
        implicit none
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief
!! \details
        subroutine init_ale_arezon_spmd(n2d,numels,numelq,trimat,nmult,ngroup,nparg, &
                                        iparg,need_to_compute,idx_list,elbuf_tab,ale_connect)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use elbufdef_mod
          use spmd_mod
          use initbuf_mod
          use ale_mod , only : ale
          use ale_connectivity_mod , only : t_ale_connectivity










          use debug_mod
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          integer, intent(in) :: n2d !< 0 : 3D, 1 : 2D
          integer, intent(in) :: numels !< number of solid elements
          integer, intent(in) :: numelq !< number of fluid elements
          integer, intent(in) :: trimat !< number of material (law51)
          integer, intent(in) :: nmult !< number of material
          integer, intent(in) :: ngroup !< number of group of element
          integer, intent(in) :: nparg !< first dimension of iparg array
          integer, dimension(nparg,ngroup), intent(in) :: iparg !< group element data
          integer, dimension(4,0:trimat,1:max(1,nmult)), intent(inout) :: need_to_compute !< output array to know if we need to compute the rezoning variable
          integer, dimension(4,0:trimat), intent(inout) :: idx_list !< list of index for each rezoning variable1
          type(elbuf_struct_), dimension(ngroup), intent(in) :: elbuf_tab !< element buffer structure
          type(t_ale_connectivity), intent(inout) :: ale_connect !< ALE data structure for connectivity          
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          logical :: done          
          integer :: itrimat,mat_number,my_check,flag_mat_eos
          integer :: nm,ijk,nvar,idx_nb,ng,idx
          integer :: mtn,llt,nft,iad,ity,npt,jale,ismstr,jeul,jtur
          integer :: jthe,jlag,jmult,jhbe,jivf,nvaux,jpor,jcvt,jclose,jplasol
          integer :: irep,iint,igtp,israt,isrot,icsen,isorth,isorthg,ifailure,jsms
          integer :: nuvar_mat,nuvar_eos,ire,irs
          integer :: ale_group_neigh_nb,ale_group_nb
          integer :: i,k,elem_id,my_address,nb_connected_elm,n_entity
          integer, dimension(4), parameter :: nvar_list = (/2,10,11,12/)
          integer, dimension(4,0:trimat,1:max(1,nmult)) :: need_to_compute_l
          integer, dimension(1:ngroup) :: group_tag
! ----------------------------------------------------------------------------------------------------------------------
!                                                   External functions
! ----------------------------------------------------------------------------------------------------------------------
                                                            
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
          ! itrimat = 0
          ! nvar = 2      / 10 / 11                               / 12
          ! idx =  1-->6  / 1  /NUM_NUVAR_MAT + NUM_NUVAR_EOS     / 1

          ! itrimat =1--> trimat
          ! nvar = 2 itrimat/=4     / 10 / 11               / 12
          ! idx =  1-->6            / 1  /NUM_NUVAR_EOS     / 1 

          if(n2d==0) then
            n_entity = numels
          else
            n_entity = numelq
          end if
          idx_list(1:4,0:trimat) = 0
          idx_list(1,0) = 6
          idx_list(2,0) = 1
          idx_list(3,0) = ale%rezon%num_nuvar_mat + ale%rezon%num_nuvar_eos
          idx_list(4,0) = 1         
          do itrimat=1,trimat
            if(itrimat/=4) idx_list(1,itrimat) = 6
            idx_list(2,itrimat) = 1
            idx_list(3,itrimat) = ale%rezon%num_nuvar_eos
            idx_list(4,itrimat) = 1         
          end do
          
          mat_number = max(1,nmult)
          need_to_compute_l(1:4,0:trimat,1:mat_number) = 0
          need_to_compute(1:4,0:trimat,1:mat_number) = 0          
          do nm=1,mat_number
            do itrimat=0,max(1,trimat)
              do ijk=1,4
                my_check = 0             
                nvar = nvar_list(ijk)
                idx_nb = idx_list(ijk,itrimat)
                do idx=1,idx_nb
                  flag_mat_eos = 0
                  if(nvar==11) then
                    if(itrimat==0.and.idx<=ale%rezon%num_nuvar_mat) flag_mat_eos = 1
                    if(itrimat==0.and.idx>ale%rezon%num_nuvar_mat) flag_mat_eos = 2
                    if(itrimat>0) flag_mat_eos = 2
                  end if
                  do ng=1,ngroup
                    call initbuf(iparg,ng,mtn,llt,nft,iad,ity,npt,jale,ismstr,jeul,jtur,   &
                                jthe,jlag,jmult,jhbe,jivf,nvaux,jpor,jcvt,jclose,jplasol, &
                                irep,iint,igtp,israt,isrot,icsen,isorth,isorthg,ifailure,jsms)
                    nuvar_mat = iparg(81,ng)
                    nuvar_eos = iparg(82,ng)
                    if(itrimat>0.and.mtn/=51) cycle
                    if(jale+jeul==0) cycle
                    if(iparg(8,ng)==1) cycle
                    if(max(1,jmult)<nm) cycle
                    ! pressurer rezoning for outlets (continuity)
                    if(jmult/=0) mtn=iparg(24+nm,ng)
                    !if(nvar==10.and.(mtn==37)) cycle
                    if(nvar==12.and.elbuf_tab(ng)%gbuf%g_temp==0) cycle
                    if(nvar==11)then
                      if(flag_mat_eos==0.or.idx==0) cycle ! %var (mat and eos)
                      if(flag_mat_eos==1)then
                        if(idx>nuvar_mat) cycle ! rezon uvar( i,  1:nuvar_mat) only
                      elseif (flag_mat_eos==2)then
                      if(idx>nuvar_eos) cycle ! rezon uvar( i,  1:nuvar_mat) only
                      endif
                      if(mtn==51 .and. (itrimat==0.or.itrimat==4))cycle
                    endif
                    irs = iparg(15,ng) !rezoning sigma enabled
                    if(nvar==2) then
                      if(irs/=1) cycle
                      if(itrimat==4) cycle
                    endif

                    ire = iparg(16,ng) !rezoning plasticity or burning time
                    if(nvar==10) then
                      if(ire/=1) cycle
                      if(mtn==41) cycle
                      if(mtn==37) cycle
                      if(itrimat/=0.and.itrimat/=4) then
                        ! ok
                      elseif(mtn==5.or.mtn==97.or.mtn==105.or.itrimat==4) then
                        ! ok
                      elseif(mtn==6) then
                        ! ok
                      elseif(mtn>=28.and.mtn/=67.and.mtn/=49) then
                        ! ok
                      else
                        if(mtn==51 .and. itrimat==0)cycle ! nok
                        ! ok
                      endif
                    endif         
                    my_check = 1
                  end do ! loop over the group of element
                end do ! loop over idx
                need_to_compute_l(ijk,itrimat,nm) = my_check
              end do ! loop over nvar
            end do ! loop over itrimat
          end do ! loop over material number
          print*," coucou ",ispmd_debug
          call spmd_allreduce(need_to_compute_l(:,1,1),need_to_compute(:,1,1),4*(trimat+1)*mat_number,spmd_max)
          print*," end coucou ",ispmd_debug
          group_tag(1:ngroup) = 0
          ale_group_nb = 0
          ale_group_neigh_nb = 0
          do ng=1,ngroup
            call initbuf(iparg,ng,mtn,llt,nft,iad,ity,npt,jale,ismstr,jeul,jtur,   &
                         jthe,jlag,jmult,jhbe,jivf,nvaux,jpor,jcvt,jclose,jplasol, &
                         irep,iint,igtp,israt,isrot,icsen,isorth,isorthg,ifailure,jsms)
            if(jale+jeul==0) cycle
            if(iparg(8,ng)==1) cycle
            group_tag(ng) = 1 ! ale group
            ale_group_nb = ale_group_nb + 1
            done = .false.
            do i=1,llt
              elem_id = nft + i
              my_address = ale_connect%ee_connect%iad_connect(elem_id)
              nb_connected_elm = ale_connect%ee_connect%iad_connect(elem_id+1) - ale_connect%ee_connect%iad_connect(elem_id) ! get the number of connected ALE element
              do k=1,nb_connected_elm
                if(.not.done.and.ale_connect%ee_connect%connected(my_address+k-1) > n_entity) then
                  group_tag(ng) = 2
                  ale_group_neigh_nb = ale_group_neigh_nb + 1
                  done = .true.
                end if
              end do
            end do ! loop over the element
          enddo ! loop over the group of element

          ale%ale_group%w_neigh_nb = ale_group_nb - ale_group_neigh_nb
          ale%ale_group%wo_neigh_nb = ale_group_neigh_nb
          allocate(ale%ale_group%list_w_neigh_nb(ale_group_neigh_nb))
          allocate(ale%ale_group%list_wo_neigh_nb(ale_group_nb-ale_group_neigh_nb))

          ale_group_nb = 0
          ale_group_neigh_nb = 0
          do ng=1,ngroup
            if(group_tag(ng)==1) then
              ale_group_nb = ale_group_nb + 1
              ale%ale_group%list_wo_neigh_nb(ale_group_nb) = ng
            endif
            if(group_tag(ng)==2) then
              ale_group_neigh_nb = ale_group_neigh_nb + 1
              ale%ale_group%list_w_neigh_nb(ale_group_neigh_nb) = ng
            end if            
          enddo
          print*,ispmd_debug,"Group : ",ale%ale_group%wo_neigh_nb,ale%ale_group%w_neigh_nb,ngroup
          !
          !
          !  phi(1:numels+nsvois) --> nsvois nombre de voisin element qu'on recoit
          !  
          !  phi(1:numels)
          ! phi qui contient element voisin
          ! com mpi par thread 0
          !
          !  com mpi phi(1:numels) --> s_buff , envoie s_buff, recoit dans r_buff, r_buff --> phi(numels+1:numels+nsvpos)
          !  ! mpi_wait + omp_barrier
          !  calcul avec phi(1:numels+nsvois)
          !
          !
          !
          !

          return
! ----------------------------------------------------------------------------------------------------------------------
        end subroutine init_ale_arezon_spmd
      end module init_ale_arezon_spmd_mod
