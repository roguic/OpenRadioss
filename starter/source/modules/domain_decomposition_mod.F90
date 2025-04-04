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
      module domain_decomposition_mod
!=======================================================================================================================
!!\brief 
!=======================================================================================================================
#include "my_real.inc"
        integer, parameter :: cell_number = 128 !< max number of cells
        integer :: newdomdec = 0 
        ! --------------------------------
        ! Structure for domdec_struct 
        type neighbour_cell_
          integer :: cell_nb !< number of cells
          integer, dimension(:), allocatable :: cell_list !< list of cells
        end type neighbour_cell_
        ! -------------------------------- 
        type domdec_
          logical, dimension(:), allocatable :: active_cell !< active cells (non-empty cells)
          integer :: edge_nb
          integer :: total_cell_nb !< total number of cells (empty and non-empty cells)
          integer :: non_empty_cell_nb !< cells with at least 1 element
          integer, dimension(3) :: cell_x_nb !< number of cells in each directions
          integer, dimension(3) :: dim_order !< order of the dimension
          integer, dimension(:), allocatable :: adj_cell_nb !< number of adjacent cells
          integer, dimension(:), allocatable :: adj_cell_connect !< number of adjacent cells
          integer, dimension(:), allocatable :: weight !< Metis weight     
          integer, dimension(:,:), allocatable :: elm_2_cell_index !< connectivity elements --> cells
          integer, dimension(:), allocatable :: cell_proc_connect !< connectivity cells --> processors
          integer, dimension(:), allocatable :: old_cell_2_new_cell !< index to convert old cell id to new cell id (old cell = empty and non-empty cells)
          integer, dimension(:), allocatable :: new_cell_2_old_cell !< index to convert new cell id to old cell id (new cell = non-empty cells)
          my_real :: cell_size !< size of cells
          my_real, dimension(3) :: max_x !< maximum nodal position
          my_real, dimension(3) :: min_x !< minimum nodal position  
          
          type(neighbour_cell_), dimension(:), allocatable :: neighbour_cell !< list of neighbour cells
        end type  domdec_
        ! --------------------------------
! ----------------------------------------------------------------------------------------------------------------------
      end module domain_decomposition_mod
