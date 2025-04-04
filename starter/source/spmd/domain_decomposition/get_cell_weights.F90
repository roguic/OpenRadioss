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
      module get_cell_weights_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Here is a small description of the routine, [after the header]
!! \details if needed, more details can be added here
        subroutine get_cell_weights( mode,nelem,ncond,iwd,domdec_struct )
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use domain_decomposition_mod , only : domdec_
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Included files
! ----------------------------------------------------------------------------------------------------------------------
#include "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          integer, intent(in) :: mode !< 0 : initialization of the global voxel (empty+non-empty cells), 1 : initialization of the non-empty cells
          integer, intent(in) :: nelem !< total number of elements
          integer, intent(in) :: ncond !< total number of Metis weights
          integer, dimension(ncond*nelem), intent(in) :: iwd !< Metis weights
          type(domdec_), intent(inout) :: domdec_struct !< data structure for the domain decomposition
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: i,j
          integer :: my_position,my_old_position
          integer :: my_address,real_cell_nb
          integer, dimension(3) :: my_cell_position
          integer, dimension(:), allocatable :: tmp_weight
! ------------------------------------------- ---------------------------------------------------------------------------
!                                                   External functions
! ----------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------------------------------------
        ! ---------------------------------
          print*," taille iwd",ncond*nelem,ncond,nelem
          print*," taille %weight",size(domdec_struct%weight)
          print*," taille boite :",domdec_struct%cell_x_nb(1:3)
          print*," mode:", mode  !< Added print statement for mode
        
          if(mode==0) then
            ! ---------------------------------            
            my_address = 0
            do i=1,nelem
              my_cell_position(1:3) = domdec_struct%elm_2_cell_index(1:3,i)
              !print*,"Cell pos :",i,my_cell_position(1),my_cell_position(2),my_cell_position(3)
              my_position = my_cell_position(1)+(my_cell_position(2)-1)*domdec_struct%cell_x_nb(1) + &
              (my_cell_position(3)-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
              !print*,"my pos :",my_position,i
              domdec_struct%active_cell(my_position) = .true.
              do j=1,ncond
                if(my_position-1+j>size(domdec_struct%weight).or.my_position-1+j<0) then
                  print*,"error weight",my_position-1+j,my_position,j,i
                  print*,"error weight 2",my_cell_position(1:3),i
                  stop
                end if
                if(my_address+j>ncond*nelem.or.my_address+j<0) then
                  print*,"error iwd",my_address+j,my_address,j,i
                  stop
                end if
                domdec_struct%weight((my_position-1)*ncond+j) = domdec_struct%weight((my_position-1)*ncond+j) + iwd((i-1)*ncond+j)
              end do
              my_address = my_address + ncond
            end do
            ! ---------------------------------
            do i=1,domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)*domdec_struct%cell_x_nb(3)
              !          print*," --- "
              !          print*," Cell weight",i
              do j=1,ncond
              !            print*,"weight",j,domdec_struct%weight((i-1)*ncond+j)
              enddo
            end do
          else
            ! ---------------------------------
            real_cell_nb = domdec_struct%non_empty_cell_nb
            allocate( tmp_weight(ncond*real_cell_nb) )
            do i=1,real_cell_nb
              my_old_position = domdec_struct%new_cell_2_old_cell(i)
              !print*," cell weight 111 ; ",i,my_old_position
              do j=1,ncond
                tmp_weight((i-1)*ncond+j) = domdec_struct%weight((my_old_position-1)*ncond+j)
              end do
            end do
            call move_alloc(tmp_weight,domdec_struct%weight)
            ! ---------------------------------
          endif
        return 
        end subroutine get_cell_weights
      end module get_cell_weights_mod