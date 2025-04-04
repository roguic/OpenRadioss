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
      module get_elm_2_cell_connectivity_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Here is a small description of the routine, [after the header]
!! \details if needed, more details can be added here
        subroutine get_elm_2_cell_connectivity( numnod,nelem,                                       &
                                                numels,numelq,numelc,numelt,numelp,numelr,numeltg,  &
                                                nixs,nixq,nixc,nixt,nixp,nixr,nixtg,                &
                                                ixs,ixq,ixc,ixt,ixp,ixr,ixtg,                       &
                                                isolnod,x,domdec_struct )
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
          integer, intent(in) :: numnod !< total number of nodes
          integer, intent(in) :: nelem !< total number of elements  
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
          integer :: node_nb,elm_type
          integer :: my_offset
          integer, dimension(3) :: my_cell_position
          my_real, dimension(3) :: barycenter
! ------------------------------------------- ---------------------------------------------------------------------------
!                                                   External functions
! ----------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------------------------------------
        if(allocated(domdec_struct%elm_2_cell_index)) deallocate(domdec_struct%elm_2_cell_index)
        allocate(domdec_struct%elm_2_cell_index(3,nelem))
        domdec_struct%elm_2_cell_index(1:3,1:nelem) = 0
        !print*," Booite :",domdec_struct%min_x(1:3),domdec_struct%max_x(1:3)
        !print*," Cell nb :",domdec_struct%cell_x_nb(1:3)
        !print*," Cell size :",domdec_struct%cell_size,nelem

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

          my_cell_position(1:3) = int(  domdec_struct%cell_x_nb(1:3)* &
                                       (barycenter(1:3) - domdec_struct%min_x(1:3)) / &
                                       (domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3)) )
                                       if(my_offset+i>size(domdec_struct%elm_2_cell_index,2).or.my_offset+i<0) then
                                        !print*,"error 1 ",my_offset+i,i
                                        stop
                                      end if
          domdec_struct%elm_2_cell_index(1:3,my_offset+i) = my_cell_position(1:3) + 1 
          if(domdec_struct%elm_2_cell_index(1,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(1).or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)>domdec_struct%cell_x_nb(2).or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(3)) then
            !print*,"error Cell pos 1 ",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
            stop
          end if 
          !print*,"Elm/Cell S",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
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

          my_cell_position(1:3) = int(  domdec_struct%cell_x_nb(1:3)* &
                                       (barycenter(1:3) - domdec_struct%min_x(1:3)) / &
                                       (domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3)) )
                                       if(my_offset+i>size(domdec_struct%elm_2_cell_index,2).or.my_offset+i<0) then
                                        !print*,"error 2 ",my_offset+i,i
                                        stop
                                      end if
          domdec_struct%elm_2_cell_index(1:3,my_offset+i) = my_cell_position(1:3) + 1  
          if(domdec_struct%elm_2_cell_index(1,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(1).or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)>domdec_struct%cell_x_nb(2).or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(3)) then
            !print*,"error Cell pos 2 ",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
            stop
          end if 
          !print*,"Elm/Cell Q ",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)  
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

          my_cell_position(1:3) = int(  domdec_struct%cell_x_nb(1:3)* &
          (barycenter(1:3) - domdec_struct%min_x(1:3)) / &
          (domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3)) )
          if(my_offset+i>size(domdec_struct%elm_2_cell_index,2).or.my_offset+i<0) then
            !print*,"error 3 ",my_offset+i,i
            stop
          end if
          domdec_struct%elm_2_cell_index(1:3,my_offset+i) = my_cell_position(1:3) + 1
          if(domdec_struct%elm_2_cell_index(1,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(1,my_offset+i)>domdec_struct%cell_x_nb(1).or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)>domdec_struct%cell_x_nb(2).or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(3)) then
            !print*,"error Cell pos 3 ",my_offset+i,i
            !print*," Computed cell pos :",domdec_struct%elm_2_cell_index(1:3,my_offset+i)
            !print*," Max cell pos :",domdec_struct%cell_x_nb(1:3)
            do j=1,node_nb
              !print*," Node pos :",ixc(1+j,i),x(1:3,ixc(1+j,i))
            enddo
            !print*," Bary :",barycenter(1:3)
            !print*,"Cell pos comput 1",(barycenter(1:3) - domdec_struct%min_x(1:3)) / &
            !(domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3))
            !print*,"Cell pos comput 2 ",domdec_struct%cell_x_nb(1:3)* &
            !(barycenter(1:3) - domdec_struct%min_x(1:3)) / &
            !(domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3))
            !print*,"Cell pos comput 3 ",int( domdec_struct%cell_x_nb(1:3)* &
            !(barycenter(1:3) - domdec_struct%min_x(1:3)) / &
            !(domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3)) )

            stop
          end if
          !print*,"Elm/Cell C",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
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

          my_cell_position(1:3) = int(  domdec_struct%cell_x_nb(1:3)* &
                                       (barycenter(1:3) - domdec_struct%min_x(1:3)) / &
                                       (domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3)) )
                                       if(my_offset+i>size(domdec_struct%elm_2_cell_index,2).or.my_offset+i<0) then
                                        !print*,"error 4 ",my_offset+i,i
                                        stop
                                      end if
          domdec_struct%elm_2_cell_index(1:3,my_offset+i) = my_cell_position(1:3) + 1    
          if(domdec_struct%elm_2_cell_index(1,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(1).or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)>domdec_struct%cell_x_nb(2).or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(3)) then
            !print*,"error Cell pos 4 ",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
            stop
          end if
          !print*,"Elm/Cell T",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
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

          my_cell_position(1:3) = int(  domdec_struct%cell_x_nb(1:3)* &
                                       (barycenter(1:3) - domdec_struct%min_x(1:3)) / &
                                       (domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3)) )
                                       if(my_offset+i>size(domdec_struct%elm_2_cell_index,2).or.my_offset+i<0) then
                                        !print*,"error 5 ",my_offset+i,i
                                        stop
                                      end if
            domdec_struct%elm_2_cell_index(1:3,my_offset+i) = my_cell_position(1:3) + 1
            if(domdec_struct%elm_2_cell_index(1,my_offset+i)<=0.or. &
            domdec_struct%elm_2_cell_index(1,my_offset+i)<=0.or. &
            domdec_struct%elm_2_cell_index(2,my_offset+i)<0.or. &
            domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(1).or. &
            domdec_struct%elm_2_cell_index(2,my_offset+i)>domdec_struct%cell_x_nb(2).or. &
            domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(3)) then
              !print*,"error Cell pos 5 ",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
              stop
            end if    
            !print*,"Elm/Cell P",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
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

          my_cell_position(1:3) = int(  domdec_struct%cell_x_nb(1:3)* &
                                       (barycenter(1:3) - domdec_struct%min_x(1:3)) / &
                                       (domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3)) )
          if(my_offset+i>size(domdec_struct%elm_2_cell_index,2).or.my_offset+i<0) then
            !print*,"error 6 ",my_offset+i,i,numelr
            stop
          end if
          domdec_struct%elm_2_cell_index(1:3,my_offset+i) = my_cell_position(1:3) + 1
          if(domdec_struct%elm_2_cell_index(1,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(1).or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)>domdec_struct%cell_x_nb(2).or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(3)) then
            !print*,"error Cell pos 6 ",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
            stop
          end if
          !print*,"Elm/Cell R",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
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
          !!print*," barycenter ",barycenter(1),barycenter(2),barycenter(3)
          my_cell_position(1:3) = int(  domdec_struct%cell_x_nb(1:3)* &
                                       (barycenter(1:3) - domdec_struct%min_x(1:3)) / &
                                       (domdec_struct%max_x(1:3) - domdec_struct%min_x(1:3)) )
                                       if(my_offset+i>size(domdec_struct%elm_2_cell_index,2).or.my_offset+i<0) then
                                        !print*,"error 7 ",my_offset+i,i
                                        stop
                                      end if
          domdec_struct%elm_2_cell_index(1:3,my_offset+i) = my_cell_position(1:3) + 1
          if(domdec_struct%elm_2_cell_index(1,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)<=0.or. &
          domdec_struct%elm_2_cell_index(1,my_offset+i)>domdec_struct%cell_x_nb(1).or. &
          domdec_struct%elm_2_cell_index(2,my_offset+i)>domdec_struct%cell_x_nb(2).or. &
          domdec_struct%elm_2_cell_index(3,my_offset+i)>domdec_struct%cell_x_nb(3)) then
            !print*,"error Cell pos 7 ",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
            stop
          end if
          !print*,"Elm/Cell TG",my_offset+i,i,domdec_struct%elm_2_cell_index(1:3,my_offset+i)
          !!print*," pos ",i,my_cell_position(1:3)
        enddo
        ! ---------------------------------

        !!print*," Cell pos 49354 : ",domdec_struct%elm_2_cell_index(1:3,49354)
        return
        end subroutine get_elm_2_cell_connectivity
      end module get_elm_2_cell_connectivity_mod