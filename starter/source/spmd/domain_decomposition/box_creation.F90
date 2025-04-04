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
      module box_creation_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Here is a small description of the routine, [after the header]
!! \details if needed, more details can be added here
        subroutine box_creation(numnod,ncond,nelem, &
              numels,numelq,numelc,numelt,numelp,numelr,numeltg,  &
              nixs,nixq,nixc,nixt,nixp,nixr,nixtg,                &
              ixs,ixq,ixc,ixt,ixp,ixr,ixtg,                       &
              isolnod,x,domdec_struct)
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use domain_decomposition_mod , only : domdec_,cell_number
          use constant_mod, only : em10
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
          integer, intent(in) :: numnod !< total number of nodes
          integer, intent(in) :: ncond !< total number of Metis weights
          integer, intent(in) :: nelem !< total number of elements)
          integer, intent(in) :: numels !< number of solid elements   
          integer, intent(in) :: numelq !< number of quad elements    
          integer, intent(in) :: numelc !< number of shell elements   
          integer, intent(in) :: numelt !< number of truss elements
          integer, intent(in) :: numelp !< number of beam elements
          integer, intent(in) :: numelr !< number of spring elements
          integer, intent(in) :: numeltg !< number of shell 3n elements 
          integer, intent(in) :: nixs !< first dimension of ixs array  
          integer, intent(in) :: nixq !< first dimension of ixq array
          integer, intent(in) :: nixc !< first dimension of ixc array
          integer, intent(in) :: nixt !< first dimension of ixt array
          integer, intent(in) :: nixp !< first dimension of ixp array
          integer, intent(in) :: nixr !< first dimension of ixr array
          integer, intent(in) :: nixtg !< first dimension of ixtg array
          integer, dimension(nelem) :: isolnod !< solid element kind
          integer, dimension(nixs,numels) :: ixs !< solid element data
          integer, dimension(nixq,numelq) :: ixq !< quad element data
          integer, dimension(nixc,numelc) :: ixc !< shell element data
          integer, dimension(nixt,numelt) :: ixt !< truss element data
          integer, dimension(nixp,numelp) :: ixp !< beam element data
          integer, dimension(nixr,numelr) :: ixr !< spring element data
          integer, dimension(nixtg,numeltg) :: ixtg !< shell 3n element data
          my_real, dimension(3,numnod), intent(in) :: x !< nodal position
          type(domdec_), intent(inout) :: domdec_struct
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: i,j
          integer :: elem_nb_1_dim
          integer :: my_dim,total_cell_nb,local_upper_bound
          integer :: node_nb,elm_type
          integer :: my_offset
          my_real :: my_max,my_size
          my_real, dimension(3) :: barycenter
          my_real, dimension(3) :: min_x,max_x
! ----------------------------------------------------------------------------------------------------------------------
!                                                   External functions
! ----------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------------------------------------
          max_x(1:3) = -huge(max_x(1:3))
          min_x(1:3) = huge(min_x(1:3))
          ! ---------------------------------
          my_offset = 0 ! offset for solid elements
          do i=1,numels
            elm_type = isolnod(i)
            if(elm_type==8.or.elm_type==16.or.elm_type==20) then
              node_nb = 8
            elseif(elm_type==10.or.elm_type==4) then
              node_nb = 4
            else
              node_nb = 8
            endif
            barycenter(1:3) = 0.
            do j=1,node_nb
              barycenter(1:3) = barycenter(1:3) + x(1:3,ixs(1+j,i))
            enddo
            barycenter(1:3) = barycenter(1:3) / node_nb
            do j=1,3
              min_x(j) = min(min_x(j),barycenter(j))
              max_x(j) = max(max_x(j),barycenter(j))
            enddo
          enddo
          ! ---------------------------------
          my_offset = my_offset + numels ! offset for quad elements
          node_nb = 4
          do i=1,numelq
            barycenter(1:3) = 0.
            do j=1,node_nb
              barycenter(1:3) = barycenter(1:3) + x(1:3,ixq(1+j,i))
            enddo
            barycenter(1:3) = barycenter(1:3) / node_nb
            do j=1,3
              min_x(j) = min(min_x(j),barycenter(j))
              max_x(j) = max(max_x(j),barycenter(j))
            enddo
          enddo
          ! ---------------------------------
          my_offset = my_offset + numelq ! offset for shell element
          node_nb = 4
          do i=1,numelc
            barycenter(1:3) = 0.
            do j=1,node_nb
              barycenter(1:3) = barycenter(1:3) + x(1:3,ixc(1+j,i))
            enddo
            barycenter(1:3) = barycenter(1:3) / node_nb
            do j=1,3
              min_x(j) = min(min_x(j),barycenter(j))
              max_x(j) = max(max_x(j),barycenter(j))
            enddo
          enddo
          ! ---------------------------------
          my_offset = my_offset + numelc ! offset for truss element
          node_nb = 3
          do i=1,numelt
            barycenter(1:3) = 0.
            do j=1,node_nb
              barycenter(1:3) = barycenter(1:3) + x(1:3,ixt(1+j,i))
            enddo
            barycenter(1:3) = barycenter(1:3) / node_nb
            do j=1,3
              min_x(j) = min(min_x(j),barycenter(j))
              max_x(j) = max(max_x(j),barycenter(j))
            enddo
          enddo
          ! ---------------------------------
          my_offset = my_offset + numelt ! offset for beam element
          node_nb = 2
          do i=1,numelp
            barycenter(1:3) = 0.
            do j=1,node_nb
              barycenter(1:3) = barycenter(1:3) + x(1:3,ixp(1+j,i))
            enddo
            barycenter(1:3) = barycenter(1:3) / node_nb
            do j=1,3
              min_x(j) = min(min_x(j),barycenter(j))
              max_x(j) = max(max_x(j),barycenter(j))
            enddo
          enddo
          ! ---------------------------------
          my_offset = my_offset + numelp ! offset for spring element
          node_nb = 2
          do i=1,numelr
            barycenter(1:3) = 0.
            do j=1,node_nb
              barycenter(1:3) = barycenter(1:3) + x(1:3,ixr(1+j,i))
            enddo
            barycenter(1:3) = barycenter(1:3) / node_nb
            do j=1,3
              min_x(j) = min(min_x(j),barycenter(j))
              max_x(j) = max(max_x(j),barycenter(j))
            enddo
          enddo
          ! ---------------------------------
          my_offset = my_offset + numelr ! offset for shell 3n elements
  
          node_nb = 3
          do i=1,numeltg
            barycenter(1:3) = 0.
            !!print*," Triangle ",i
            do j=1,node_nb
              barycenter(1:3) = barycenter(1:3) + x(1:3,ixtg(1+j,i))            
            enddo
            
            barycenter(1:3) = barycenter(1:3) / node_nb
            do j=1,3
              min_x(j) = min(min_x(j),barycenter(j))
              max_x(j) = max(max_x(j),barycenter(j))
            enddo
          enddo
          ! ---------------------------------


          print*," nelem :",nelem,nelem**(1./3.)
          elem_nb_1_dim = nelem**(1./3.)
          print*," elem_nb_1_dim : ",elem_nb_1_dim
          ! at least 1 element per cell
          local_upper_bound = min(elem_nb_1_dim,cell_number)

!          do i=1,numnod
!            do j=1,3
!              min_x(j) = min(min_x(j),x(j,i))
!              max_x(j) = max(max_x(j),x(j,i))
!            enddo
            !print*," node pos: ",i,x(1,i),x(2,i),x(3,i)
          !enddo
          print*," min pos: ",min_x(1),min_x(2),min_x(3)
          print*," max pos: ",max_x(1),max_x(2),max_x(3)

          do i=1,3
            min_x(i) = min_x(i) - em10
            max_x(i) = max_x(i) + em10
          enddo
          domdec_struct%min_x(1:3) = min_x(1:3)
          domdec_struct%max_x(1:3) = max_x(1:3)          
          my_max = -huge(my_max)
          do i=1,3
            print*," Dx :",i,domdec_struct%max_x(i)-domdec_struct%min_x(i),domdec_struct%max_x(i),domdec_struct%min_x(i)
            if(domdec_struct%max_x(i)-domdec_struct%min_x(i)>my_max) then
              my_max = domdec_struct%max_x(i)-domdec_struct%min_x(i)
              domdec_struct%dim_order(1) = i
            endif
          enddo
          my_max = -huge(my_max)
          do i=1,3
            if(i/= domdec_struct%dim_order(1).and.domdec_struct%max_x(i)-domdec_struct%min_x(i)>my_max) then
              my_max = domdec_struct%max_x(i)-domdec_struct%min_x(i)
              domdec_struct%dim_order(2) = i
            endif
          enddo

          do i=1,3
            if(i/= domdec_struct%dim_order(1).and.i/= domdec_struct%dim_order(2)) then
              domdec_struct%dim_order(3) = i
            endif
          enddo

          my_dim = domdec_struct%dim_order(1)
          my_size = domdec_struct%max_x(my_dim)-domdec_struct%min_x(my_dim)
          domdec_struct%cell_x_nb(my_dim) = local_upper_bound
          domdec_struct%cell_size = my_size/local_upper_bound

          do i=2,3
            my_dim = domdec_struct%dim_order(i)
            my_size = domdec_struct%max_x(my_dim)-domdec_struct%min_x(my_dim)
            domdec_struct%cell_x_nb(my_dim) = my_size/domdec_struct%cell_size
            domdec_struct%cell_x_nb(my_dim) = max(1,domdec_struct%cell_x_nb(my_dim))
            domdec_struct%cell_x_nb(my_dim) = min(local_upper_bound,domdec_struct%cell_x_nb(my_dim))
          enddo

          if(allocated(domdec_struct%weight)) deallocate(domdec_struct%weight)
          total_cell_nb = domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)*domdec_struct%cell_x_nb(3)
          domdec_struct%total_cell_nb = total_cell_nb
          allocate(domdec_struct%weight(ncond*total_cell_nb))
          domdec_struct%weight(1:ncond*total_cell_nb) = 0
          if(allocated(domdec_struct%cell_proc_connect)) deallocate(domdec_struct%cell_proc_connect)
          allocate(domdec_struct%cell_proc_connect(total_cell_nb))
          domdec_struct%cell_proc_connect(1:total_cell_nb) = 0

          allocate(domdec_struct%active_cell(total_cell_nb))
          domdec_struct%active_cell(1:total_cell_nb) = .false.

          return
        end subroutine box_creation
      end module box_creation_mod
