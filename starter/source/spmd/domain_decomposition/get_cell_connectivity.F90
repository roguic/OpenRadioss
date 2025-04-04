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
      module get_cell_connectivity_mod
      contains
! ======================================================================================================================
!                                                   procedures
! ======================================================================================================================
!! \brief Here is a small description of the routine, [after the header]
!! \details if needed, more details can be added here
        subroutine get_cell_connectivity( domdec_struct )
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
          type(domdec_), intent(inout) :: domdec_struct
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          logical :: break_
          integer :: i,j,k,ijk
          integer :: ix,iy,iz
          integer :: my_size,my_old_size
          integer :: old, next
          integer :: my_offset
          integer :: real_cell_nb
          integer, dimension(3) :: my_cell_position
          integer, dimension(3) :: low_bound,up_bound
          logical, dimension(:),allocatable :: empty_cell
          integer :: my_position,my_address,my_neighbour  
          integer, dimension(:), allocatable :: counter,tmp_adj_cell_nb,tmp_adj_cell_connect
          integer, dimension(:,:), allocatable :: new_cell_2_old_cell
          integer :: level_nb,level,global_level
          integer, dimension(:), allocatable :: level_index,old_list  







          logical :: new_value,break_2
          integer :: my_color,old_next,next_2
          integer, allocatable :: tmp_list(:),color(:),tmp_list_2(:)

! ------------------------------------------- ---------------------------------------------------------------------------
!                                                   External functions
! ----------------------------------------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------------------------------------

        
        ! ---------------------------------
        if(allocated(domdec_struct%adj_cell_nb)) deallocate(domdec_struct%adj_cell_nb)
        my_size = domdec_struct%total_cell_nb
        print*," *** my_size : *** ",my_size
        allocate(domdec_struct%adj_cell_nb(my_size+1))
        domdec_struct%adj_cell_nb(1:my_size+1) = 0
        print*," ***End my_size : *** ",my_size
        ! ---------------------------------

        ! ---------------------------------
        domdec_struct%adj_cell_nb(1) = 1
        ! loop on the cells to compute the number of adjacent cells
        do i=1,domdec_struct%cell_x_nb(1)
          do j=1,domdec_struct%cell_x_nb(2)
            do k=1,domdec_struct%cell_x_nb(3)
              ! get the cell's position
              my_position = i + (j-1)*domdec_struct%cell_x_nb(1) + (k-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
              my_position = my_position + 1 ! offset + 1
              !print*," ** en cours :",my_position
              if(i>1) domdec_struct%adj_cell_nb(my_position) = domdec_struct%adj_cell_nb(my_position) + 1
              if(i<domdec_struct%cell_x_nb(1)) domdec_struct%adj_cell_nb(my_position) = domdec_struct%adj_cell_nb(my_position) + 1
              if(j>1) domdec_struct%adj_cell_nb(my_position) = domdec_struct%adj_cell_nb(my_position) + 1
              if(j<domdec_struct%cell_x_nb(2)) domdec_struct%adj_cell_nb(my_position) = domdec_struct%adj_cell_nb(my_position) + 1
              if(k>1) domdec_struct%adj_cell_nb(my_position) = domdec_struct%adj_cell_nb(my_position) + 1
              if(k<domdec_struct%cell_x_nb(3)) domdec_struct%adj_cell_nb(my_position) = domdec_struct%adj_cell_nb(my_position) + 1
              !print*," ** en cours :",my_position-1, domdec_struct%adj_cell_nb(my_position)
            end do
          end do
        end do
        print*," *** NB cell : *** ",domdec_struct%cell_x_nb(1),domdec_struct%cell_x_nb(2),domdec_struct%cell_x_nb(3)
        ! %adj_cell_nb contains now the address to %adj_cell_connect array
        ! and the number of adjacent cells : number = %adj_cell_nb(n+1) - %adj_cell_nb(n) 
        do i=1,my_size
          domdec_struct%adj_cell_nb(i+1) = domdec_struct%adj_cell_nb(i) + domdec_struct%adj_cell_nb(i+1)
        end do
        do i=1,my_size
          !print*," *** NB cell adjacent *** ",i,  &
          !domdec_struct%adj_cell_nb(i+1)-domdec_struct%adj_cell_nb(i)
        enddo
        ! --------------------------------- 
        domdec_struct%edge_nb = domdec_struct%adj_cell_nb(my_size+1)-domdec_struct%adj_cell_nb(1)
        domdec_struct%edge_nb = domdec_struct%edge_nb / 2
        if(allocated(domdec_struct%adj_cell_connect)) deallocate(domdec_struct%adj_cell_connect)
        allocate(domdec_struct%adj_cell_connect(2*domdec_struct%edge_nb))
        domdec_struct%adj_cell_connect(1:2*domdec_struct%edge_nb) = 0
        ! ---------------------------------
        ! loop on the cells to find the cell's connectivity
        allocate( counter(my_size) )
        counter(1:my_size) = 0
        do i=1,domdec_struct%cell_x_nb(1)
          do j=1,domdec_struct%cell_x_nb(2)
            do k=1,domdec_struct%cell_x_nb(3)
              ! get the cell's position
              my_position = i + (j-1)*domdec_struct%cell_x_nb(1) + (k-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
              if(i>1) then                
                counter(my_position) = counter(my_position) + 1
                my_address = domdec_struct%adj_cell_nb(my_position) - 1 + counter(my_position)
                my_neighbour =(i-1) + (j-1)*domdec_struct%cell_x_nb(1) + (k-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
                domdec_struct%adj_cell_connect(my_address) = my_neighbour
              endif
              if(i<domdec_struct%cell_x_nb(1)) then
                counter(my_position) = counter(my_position) + 1
                my_address = domdec_struct%adj_cell_nb(my_position) - 1 + counter(my_position)
                my_neighbour =(i+1) + (j-1)*domdec_struct%cell_x_nb(1) + (k-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
                domdec_struct%adj_cell_connect(my_address) = my_neighbour
              endif              
              if(j>1) then
                counter(my_position) = counter(my_position) + 1
                my_address = domdec_struct%adj_cell_nb(my_position) - 1 + counter(my_position)
                my_neighbour = i + (j-1-1)*domdec_struct%cell_x_nb(1) + (k-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
                domdec_struct%adj_cell_connect(my_address) = my_neighbour
              endif
                
              if(j<domdec_struct%cell_x_nb(2)) then
                counter(my_position) = counter(my_position) + 1
                my_address = domdec_struct%adj_cell_nb(my_position) - 1 + counter(my_position)
                my_neighbour = i + (j+1-1)*domdec_struct%cell_x_nb(1) + (k-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
                domdec_struct%adj_cell_connect(my_address) = my_neighbour
              endif

              if(k>1) then
                counter(my_position) = counter(my_position) + 1
                my_address = domdec_struct%adj_cell_nb(my_position) - 1 + counter(my_position)
                my_neighbour = i + (j-1)*domdec_struct%cell_x_nb(1) + (k-1-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
                domdec_struct%adj_cell_connect(my_address) = my_neighbour
              endif
              if(k<domdec_struct%cell_x_nb(3)) then
                counter(my_position) = counter(my_position) + 1
                my_address = domdec_struct%adj_cell_nb(my_position) - 1 + counter(my_position)
                my_neighbour = i + (j-1)*domdec_struct%cell_x_nb(1) + (k+1-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
                domdec_struct%adj_cell_connect(my_address) = my_neighbour
              endif
            enddo
          enddo
        enddo

        do i=1,domdec_struct%cell_x_nb(1)
          do j=1,domdec_struct%cell_x_nb(2)
            do k=1,domdec_struct%cell_x_nb(3)
              ! get the cell's position
              my_position = i + (j-1)*domdec_struct%cell_x_nb(1) + (k-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
              !print*," ** Cell pos :",my_position,domdec_struct%adj_cell_nb(my_position+1)-domdec_struct%adj_cell_nb(my_position)
              my_address = domdec_struct%adj_cell_nb(my_position)
              do ijk=1,domdec_struct%adj_cell_nb(my_position+1)-domdec_struct%adj_cell_nb(my_position)                
                my_neighbour = domdec_struct%adj_cell_connect(my_address+ijk-1)
                !print*," Adj cell pos :",my_neighbour
              end do
            enddo
          enddo 
        enddo
        ! ---------------------------------

        ! ---------------------------------
        allocate( empty_cell(domdec_struct%total_cell_nb) )
        allocate( tmp_adj_cell_nb(domdec_struct%total_cell_nb+1) )
        allocate( tmp_adj_cell_connect(125) )
        allocate( new_cell_2_old_cell(domdec_struct%total_cell_nb,4) )
        new_cell_2_old_cell(1:domdec_struct%total_cell_nb,1:4) = 0
        tmp_adj_cell_nb(1:domdec_struct%total_cell_nb+1) = 0
        tmp_adj_cell_connect(1:125) = 0

        tmp_adj_cell_nb(1) = 1
        ! ---------------------------------

        real_cell_nb= 0
        empty_cell(1:domdec_struct%total_cell_nb) = .false.
        do i=1,domdec_struct%total_cell_nb
          if(.not.domdec_struct%active_cell(i)) then
            empty_cell(i) = .true.
          endif          
        end do


        do i=1,domdec_struct%cell_x_nb(1)
          do j=1,domdec_struct%cell_x_nb(2)
            do k=1,domdec_struct%cell_x_nb(3)
              ! get the cell's position
              my_position = i + (j-1)*domdec_struct%cell_x_nb(1) + (k-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
              my_address = domdec_struct%adj_cell_nb(my_position)
              if(empty_cell(my_position)) then
                domdec_struct%adj_cell_connect(my_address) = -my_position
              else
                if(k==1)print*,"Cell Active",i,j,k,my_position,real_cell_nb+1
                domdec_struct%adj_cell_connect(my_address) = my_position
                real_cell_nb = real_cell_nb + 1
                new_cell_2_old_cell(real_cell_nb,1) = my_position
                new_cell_2_old_cell(real_cell_nb,2) = i
                new_cell_2_old_cell(real_cell_nb,3) = j
                new_cell_2_old_cell(real_cell_nb,4) = k
              endif
            enddo
          enddo 
        enddo   

        allocate(domdec_struct%old_cell_2_new_cell(my_size))
        my_size = domdec_struct%total_cell_nb
        allocate(domdec_struct%new_cell_2_old_cell(real_cell_nb))        
        domdec_struct%old_cell_2_new_cell(1:my_size) = 0
        domdec_struct%new_cell_2_old_cell(1:real_cell_nb) = 0
        do i=1,real_cell_nb
          print*," Vrai Cell ",i,new_cell_2_old_cell(i,2:4)
          old = new_cell_2_old_cell(i,1)
          domdec_struct%new_cell_2_old_cell(i) = new_cell_2_old_cell(i,1)
          domdec_struct%old_cell_2_new_cell(old) = i
        enddo

        allocate( level_index(real_cell_nb))
        level_nb = 0
        allocate(domdec_struct%neighbour_cell(real_cell_nb))
        do i=1,real_cell_nb          
          domdec_struct%neighbour_cell(i)%cell_nb = 0
        enddo
        global_level = 0
        do i=1,real_cell_nb
          my_position = new_cell_2_old_cell(i,1)
          my_cell_position(1:3) = new_cell_2_old_cell(i,2:4)
          next = 0
          print*," PoSiTiOn :",i,new_cell_2_old_cell(i,1)," ++ ",my_cell_position(1:3)
          low_bound(1:3) = my_cell_position(1:3)
          up_bound(1:3) = my_cell_position(1:3)
          if( my_cell_position(1)>1 ) low_bound(1) = my_cell_position(1) - 1
          if( my_cell_position(1)<domdec_struct%cell_x_nb(1) ) up_bound(1) = my_cell_position(1) + 1
          if( my_cell_position(2)>1 ) low_bound(2) = my_cell_position(2) - 1
          if( my_cell_position(2)<domdec_struct%cell_x_nb(2) ) up_bound(2) = my_cell_position(2) + 1
          if( my_cell_position(3)>1 ) low_bound(3) = my_cell_position(3) - 1
          if( my_cell_position(3)<domdec_struct%cell_x_nb(3) ) up_bound(3) = my_cell_position(3) + 1

          do iz = low_bound(3),up_bound(3)
            do iy = low_bound(2),up_bound(2)
              do ix = low_bound(1),up_bound(1)
                my_neighbour = ix + (iy-1)*domdec_struct%cell_x_nb(1) + (iz-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
                if(empty_cell(my_neighbour)) then
                  my_neighbour = -1
                else
                  if(my_position==my_neighbour) then
                    my_neighbour = -1
                  else
                    next = next + 1
                    tmp_adj_cell_connect(next) = domdec_struct%old_cell_2_new_cell(my_neighbour)
                  endif
                endif
              enddo
            enddo
          enddo

          if(next>0) then
            domdec_struct%neighbour_cell(i)%cell_nb = next
            if(allocated(domdec_struct%neighbour_cell(i)%cell_list)) deallocate(domdec_struct%neighbour_cell(i)%cell_list)
            allocate(domdec_struct%neighbour_cell(i)%cell_list(next))
            domdec_struct%neighbour_cell(i)%cell_list(1:next) = tmp_adj_cell_connect(1:next)
            global_level = max(global_level,1)
          else
            level = 1
            break_ = .true.
            next = 0
            do while(break_)
              print*," Level en cours :",level,my_cell_position(1:3)
              level = level + 1
              my_size = (1 + 2*level)**3
              if(my_size>size(tmp_adj_cell_connect)) then
                deallocate(tmp_adj_cell_connect)
                allocate(tmp_adj_cell_connect(my_size))
              endif
              my_cell_position(1:3) = new_cell_2_old_cell(i,2:4)
              low_bound(1:3) = my_cell_position(1:3)
              up_bound(1:3) = my_cell_position(1:3)
              if( my_cell_position(1)>1 ) low_bound(1) = my_cell_position(1) - level
              if( my_cell_position(1)<domdec_struct%cell_x_nb(1) ) up_bound(1) = my_cell_position(1) + level
              if( my_cell_position(2)>1 ) low_bound(2) = my_cell_position(2) - level
              if( my_cell_position(2)<domdec_struct%cell_x_nb(2) ) up_bound(2) = my_cell_position(2) + level
              if( my_cell_position(3)>1 ) low_bound(3) = my_cell_position(3) - level
              if( my_cell_position(3)<domdec_struct%cell_x_nb(3) ) up_bound(3) = my_cell_position(3) + level
              do iz = low_bound(3),up_bound(3)
                do iy = low_bound(2),up_bound(2)
                  do ix = low_bound(1),up_bound(1)
                    my_neighbour = ix + (iy-1)*domdec_struct%cell_x_nb(1) + &
                                  (iz-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
                    if(empty_cell(my_neighbour)) then
                      my_neighbour = -1
                    else
                      if(my_position==my_neighbour) then
                        my_neighbour = -1
                      else
                        break_ = .false.
                        next = next + 1
                        tmp_adj_cell_connect(next) = domdec_struct%old_cell_2_new_cell(my_neighbour)
                      endif
                    endif
                  enddo
                enddo
              enddo
              print*," Level en cours WhIlE :",next,my_cell_position(1:3)
            enddo
            global_level = max(global_level,level)
            level_nb = level_nb + 1
            level_index(level_nb) = i
            level = max(level,2)
            domdec_struct%neighbour_cell(i)%cell_nb = next
            if(allocated(domdec_struct%neighbour_cell(i)%cell_list)) deallocate(domdec_struct%neighbour_cell(i)%cell_list)
            allocate(domdec_struct%neighbour_cell(i)%cell_list(next))
            domdec_struct%neighbour_cell(i)%cell_list(1:next) = tmp_adj_cell_connect(1:next)
          endif
        enddo     

!        do i=1,real_cell_nb
!          tmp_adj_cell_nb(i+1) = tmp_adj_cell_nb(i+1) + tmp_adj_cell_nb(i)
!        enddo


        counter(1:real_cell_nb) = 0
        if(global_level>1) then
          print*,' LeVeL :',level_nb, global_level
          do i=1,level_nb
            my_position = level_index(i)
            print*,' LeVeL CeLl 0:',my_position
            do j=1,domdec_struct%neighbour_cell(my_position)%cell_nb
              my_neighbour = domdec_struct%neighbour_cell(my_position)%cell_list(j)
              print*,' LeVeL CeLl 0 neigh:',my_position,my_neighbour
              counter(my_neighbour) = counter(my_neighbour) + 1
            enddo
          enddo
          do i=1,real_cell_nb
            if(counter(i)>0) then
              my_old_size = domdec_struct%neighbour_cell(i)%cell_nb
              print*," Add CELL :",i,counter(i),my_old_size,counter(i)+my_old_size
              my_size = counter(i) + my_old_size
              allocate( old_list(my_size) )               
              old_list(1:my_old_size) = domdec_struct%neighbour_cell(i)%cell_list(1:my_old_size)
              old_list(my_old_size+1:my_size) = 0
              call move_alloc( old_list,domdec_struct%neighbour_cell(i)%cell_list )              
            endif
          enddo

          do i=1,level_nb
            my_position = level_index(i)
            do j=1,domdec_struct%neighbour_cell(my_position)%cell_nb              
              my_neighbour = domdec_struct%neighbour_cell(my_position)%cell_list(j) ! get the neighbour cell 
              print*,' LeVeL CeLl 1 neigh:',my_neighbour,my_position
              domdec_struct%neighbour_cell(my_neighbour)%cell_nb = domdec_struct%neighbour_cell(my_neighbour)%cell_nb + 1
              domdec_struct%neighbour_cell(my_neighbour)%cell_list(domdec_struct%neighbour_cell(my_neighbour)%cell_nb) = my_position
            enddo
          enddo 
        endif
        
        my_size = 0
        do i=1,real_cell_nb
          my_size = domdec_struct%neighbour_cell(i)%cell_nb + my_size
        enddo
        
        domdec_struct%non_empty_cell_nb = real_cell_nb
        print*," non-empty cell nb :",domdec_struct%non_empty_cell_nb,domdec_struct%total_cell_nb
        deallocate(domdec_struct%adj_cell_nb)
        deallocate(domdec_struct%adj_cell_connect)
        allocate( domdec_struct%adj_cell_nb(real_cell_nb+1) )
        allocate( domdec_struct%adj_cell_connect(my_size) )
        domdec_struct%adj_cell_connect(1:my_size) = 0
        domdec_struct%adj_cell_nb(1) = 1
        do i=1,real_cell_nb
          domdec_struct%adj_cell_nb(i+1) = domdec_struct%adj_cell_nb(i) + domdec_struct%neighbour_cell(i)%cell_nb
        enddo

        if(domdec_struct%adj_cell_nb(real_cell_nb+1)-domdec_struct%adj_cell_nb(1)/=my_size) then
          print*,'error adj_cell_nb',domdec_struct%adj_cell_nb(real_cell_nb+1),my_size
          stop
        endif

        do i=1,real_cell_nb
          my_size = domdec_struct%neighbour_cell(i)%cell_nb
          do j=1,my_size
            my_address = domdec_struct%adj_cell_nb(i) - 1 + j
            domdec_struct%adj_cell_connect(my_address) = domdec_struct%neighbour_cell(i)%cell_list(j)
          enddo
        enddo

        do i=1,real_cell_nb+1
          if(domdec_struct%adj_cell_nb(i)<=0) then
            print*,'error reduction cell empty 1',i,domdec_struct%adj_cell_nb(i)
            stop
          endif
        enddo
        do i=1,next
        !    print*,'reduction cell empty 2',i,domdec_struct%adj_cell_connect(i)
        enddo
        do i=1,next
          if(domdec_struct%adj_cell_connect(i)<=0) then
            print*,'error reduction cell empty 2',i,domdec_struct%adj_cell_connect(i)
            stop
          endif
        enddo

        !do i=1,real_cell_nb
        !  print*," ++ Cell :",i," Nb adj cell :",domdec_struct%adj_cell_nb(i+1)-domdec_struct%adj_cell_nb(i)
        !  do j=domdec_struct%adj_cell_nb(i),domdec_struct%adj_cell_nb(i+1)-1
        !    print*,"++ Cell :",i," Adj cell :",j,domdec_struct%adj_cell_connect(j) 
        !    if(domdec_struct%adj_cell_connect(j)<= 0) then
        !      print*,'error reduction cell empty 3',i,j,domdec_struct%adj_cell_connect(j)
        !      stop
        !    endif           
        !  enddo
        !enddo    

        my_color = 1
        break_ = .true.
        allocate( color(real_cell_nb) )
        color(1:real_cell_nb) = 0
        allocate( tmp_list(real_cell_nb) )
        allocate( tmp_list_2(real_cell_nb) )
        tmp_list(1:real_cell_nb) = 0
        tmp_list(1) = 1
        tmp_list_2(1) = 1
        next = 1
        next_2 = 1
        old_next = next
        do while( break_)
          old_next = next
          new_value = .false.
          next = next_2
          print*," NB tmp_list(i)" , next
          tmp_list(1:next_2) = tmp_list_2(1:next_2)
          next_2 = 0 
          do i =1,next
            print*," tmp_list(i) =",tmp_list(i),next
            if(color(tmp_list(i))==0.or.color(tmp_list(i))==my_color) then 
              color(tmp_list(i)) = my_color 
              old_next = next          
              do j=domdec_struct%adj_cell_nb(tmp_list(i)),domdec_struct%adj_cell_nb(tmp_list(i)+1)-1
                if(color(domdec_struct%adj_cell_connect(j))==0) then
                  print*," tmp_list(i) Voisin :", &
                  domdec_struct%adj_cell_connect(j),"color:",color(domdec_struct%adj_cell_connect(j)),i,next
                  color(domdec_struct%adj_cell_connect(j)) = my_color
                  next_2 = next_2 + 1
                  tmp_list_2(next_2) = domdec_struct%adj_cell_connect(j)
                  new_value = .true.
                elseif(color(domdec_struct%adj_cell_connect(j))/=my_color) then
                  print*," erreur de coloration ",my_color, &
                  color(domdec_struct%adj_cell_connect(j)),tmp_list(i),domdec_struct%adj_cell_connect(j)
                  stop
                endif  
              enddo           
            endif
          enddo

          if(.not.new_value) then
            do i=1,next_2
              color(tmp_list_2(i)) = my_color
            enddo
            next = 1
            old_next = 1
            break_ = .false.
            break_2 = .true.
            i = 1
            do while(break_2)
              if(color(i)==0) then
                my_color = my_color + 1
                tmp_list_2(1) = i
                next = 1
                next_2 = 1
                print*," je change :",i,my_color
                if(my_color>real_cell_nb) then
                  print*," Error ColoR :",my_color,real_cell_nb
                  stop
                endif
                break_2 = .false.
                break_ = .true.
              endif
              i = i + 1
              if(i>real_cell_nb) then
                break_2 = .false.
                print*," J ai termine !!! ",my_color
              endif
            enddo
          endif
        enddo
        
        next = 0
        counter(1:real_cell_nb) = 0
        do i=1,real_cell_nb
          if(color(i)/=1) then
            next = next + 1
            tmp_list(next) = i
            counter(i) = 1
            print*," Colour differente",i,color(i)
          endif! 22542       22023
        enddo     
        print*," Nb couleur differente :",next

        do ijk=1,next
          i = tmp_list(ijk)
          my_position = new_cell_2_old_cell(i,1)
          my_cell_position(1:3) = new_cell_2_old_cell(i,2:4)
          low_bound(1:3) = my_cell_position(1:3)
          up_bound(1:3) = my_cell_position(1:3)

          level = 1
          break_ = .true.
          next_2 = 0
          do while(break_)
            print*," -------------------------------- "
            print*,"**** Color en cours :",i,my_position,my_cell_position(1:3)
            next_2 = 0
            level = level + 1
            my_size = (1 + 2*level)**3
            if(my_size>size(tmp_adj_cell_connect)) then
              deallocate(tmp_adj_cell_connect)
              allocate(tmp_adj_cell_connect(my_size))
            endif
            if( my_cell_position(1)>1 ) low_bound(1) = my_cell_position(1) - level
            if( my_cell_position(1)<domdec_struct%cell_x_nb(1) ) up_bound(1) = my_cell_position(1) + level
            if( my_cell_position(2)>1 ) low_bound(2) = my_cell_position(2) - level
            if( my_cell_position(2)<domdec_struct%cell_x_nb(2) ) up_bound(2) = my_cell_position(2) + level
            if( my_cell_position(3)>1 ) low_bound(3) = my_cell_position(3) - level
            if( my_cell_position(3)<domdec_struct%cell_x_nb(3) ) up_bound(3) = my_cell_position(3) + level
            do iz = low_bound(3),up_bound(3)
              do iy = low_bound(2),up_bound(2)
                do ix = low_bound(1),up_bound(1)
                  my_neighbour = ix + (iy-1)*domdec_struct%cell_x_nb(1) + &
                                (iz-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
                  if(empty_cell(my_neighbour)) then
                    my_neighbour = -1
                  else
                    if(my_position==my_neighbour) then
                      my_neighbour = -1
                    else
                      next_2 = next_2 + 1
                      tmp_adj_cell_connect(next_2) = domdec_struct%old_cell_2_new_cell(my_neighbour)
                      print*,"**** VoiSin color :",level,i,domdec_struct%old_cell_2_new_cell(my_neighbour)
                      if(counter(domdec_struct%old_cell_2_new_cell(my_neighbour))==0) then
                        break_ = .false.
                      endif
                    endif
                  endif
                enddo
              enddo
            enddo
            print*," Level en cours WhIlE :",next,my_cell_position(1:3)
          enddo
          domdec_struct%neighbour_cell(i)%cell_nb = next_2
          if(allocated(domdec_struct%neighbour_cell(i)%cell_list)) deallocate(domdec_struct%neighbour_cell(i)%cell_list)
          allocate(domdec_struct%neighbour_cell(i)%cell_list(next_2))
          domdec_struct%neighbour_cell(i)%cell_list(1:next_2) = tmp_adj_cell_connect(1:next_2)
        enddo

        print*," Nb color :",my_color 
        
        counter(1:real_cell_nb) = 0
        do i=1,next
          my_position = tmp_list(i)
          print*,'CoLoR CeLl 0:',my_position
          do j=1,domdec_struct%neighbour_cell(my_position)%cell_nb
            my_neighbour = domdec_struct%neighbour_cell(my_position)%cell_list(j)
            print*,' CoLoR CeLl 0 neigh:',my_position,my_neighbour
            counter(my_neighbour) = counter(my_neighbour) + 1
          enddo
        enddo
        do i=1,real_cell_nb
          if(counter(i)>0) then
            my_old_size = domdec_struct%neighbour_cell(i)%cell_nb
            print*,"CoLoR Add CELL :",i,counter(i),my_old_size,counter(i)+my_old_size
            my_size = counter(i) + my_old_size
            allocate( old_list(my_size) )               
            old_list(1:my_old_size) = domdec_struct%neighbour_cell(i)%cell_list(1:my_old_size)
            old_list(my_old_size+1:my_size) = 0
            call move_alloc( old_list,domdec_struct%neighbour_cell(i)%cell_list )              
          endif
        enddo

        do i=1,next
          my_position = tmp_list(i)
          do j=1,domdec_struct%neighbour_cell(my_position)%cell_nb              
            my_neighbour = domdec_struct%neighbour_cell(my_position)%cell_list(j) ! get the neighbour cell 
            print*,'CoLoR LeVeL CeLl 1 neigh:',my_neighbour,my_position
            domdec_struct%neighbour_cell(my_neighbour)%cell_nb = domdec_struct%neighbour_cell(my_neighbour)%cell_nb + 1
            domdec_struct%neighbour_cell(my_neighbour)%cell_list(domdec_struct%neighbour_cell(my_neighbour)%cell_nb) = my_position
          enddo
        enddo 
        
        my_size = 0
        do i=1,real_cell_nb
          my_size = domdec_struct%neighbour_cell(i)%cell_nb + my_size
        enddo
        
        domdec_struct%non_empty_cell_nb = real_cell_nb
        print*,"22 non-empty cell nb :",domdec_struct%non_empty_cell_nb,domdec_struct%total_cell_nb
        deallocate(domdec_struct%adj_cell_nb)
        deallocate(domdec_struct%adj_cell_connect)
        allocate( domdec_struct%adj_cell_nb(real_cell_nb+1) )
        allocate( domdec_struct%adj_cell_connect(my_size) )
        domdec_struct%adj_cell_connect(1:my_size) = 0
        domdec_struct%adj_cell_nb(1) = 1
        do i=1,real_cell_nb
          domdec_struct%adj_cell_nb(i+1) = domdec_struct%adj_cell_nb(i) + domdec_struct%neighbour_cell(i)%cell_nb
        enddo

        if(domdec_struct%adj_cell_nb(real_cell_nb+1)-domdec_struct%adj_cell_nb(1)/=my_size) then
          print*,'22 error adj_cell_nb',domdec_struct%adj_cell_nb(real_cell_nb+1),my_size
          stop
        endif

        do i=1,real_cell_nb
          my_size = domdec_struct%neighbour_cell(i)%cell_nb
          do j=1,my_size
            my_address = domdec_struct%adj_cell_nb(i) - 1 + j
            domdec_struct%adj_cell_connect(my_address) = domdec_struct%neighbour_cell(i)%cell_list(j)
          enddo
        enddo

        do i=1,real_cell_nb
          print*," ++ TROIZE Cell :",i," Nb adj cell :",domdec_struct%adj_cell_nb(i+1)-domdec_struct%adj_cell_nb(i)
          do j=domdec_struct%adj_cell_nb(i),domdec_struct%adj_cell_nb(i+1)-1
            print*,"++TROIZE Cell :",i," Adj cell :",j,domdec_struct%adj_cell_connect(j) 
            if(domdec_struct%adj_cell_connect(j)<= 0) then
              print*,'TROIZE error reduction cell empty 3',i,j,domdec_struct%adj_cell_connect(j)
              stop
            endif           
          enddo
        enddo    

        do i=1,real_cell_nb+1
          if(domdec_struct%adj_cell_nb(i)<=0) then
            print*,'22 error reduction cell empty 1',i,domdec_struct%adj_cell_nb(i)
            stop
          endif
        enddo
        do i=1,next
        !    print*,'reduction cell empty 2',i,domdec_struct%adj_cell_connect(i)
        enddo
        do i=1,next
          if(domdec_struct%adj_cell_connect(i)<=0) then
            print*,'22 error reduction cell empty 2',i,domdec_struct%adj_cell_connect(i)
            stop
          endif
        enddo















        my_color = 1
        break_ = .true.
        color(1:real_cell_nb) = 0
        tmp_list(1:real_cell_nb) = 0
        tmp_list(1) = 1
        tmp_list_2(1) = 1
        next = 1
        next_2 = 1
        old_next = next
        do while( break_)
          old_next = next
          new_value = .false.
          next = next_2
          print*,"DEUZE NB tmp_list(i)" , next
          tmp_list(1:next_2) = tmp_list_2(1:next_2)
          next_2 = 0 
          do i =1,next
            print*,"DEUZE tmp_list(i) =",tmp_list(i),next
            if(color(tmp_list(i))==0.or.color(tmp_list(i))==my_color) then 
              color(tmp_list(i)) = my_color 
              old_next = next          
              do j=domdec_struct%adj_cell_nb(tmp_list(i)),domdec_struct%adj_cell_nb(tmp_list(i)+1)-1
                if(color(domdec_struct%adj_cell_connect(j))==0) then
                  print*," tmp_list(i) Voisin :", &
                  domdec_struct%adj_cell_connect(j),"color:",color(domdec_struct%adj_cell_connect(j)),i,next
                  color(domdec_struct%adj_cell_connect(j)) = my_color
                  next_2 = next_2 + 1
                  tmp_list_2(next_2) = domdec_struct%adj_cell_connect(j)
                  new_value = .true.
                elseif(color(domdec_struct%adj_cell_connect(j))/=my_color) then
                  print*," erreur de coloration ",my_color, &
                  color(domdec_struct%adj_cell_connect(j)),tmp_list(i),domdec_struct%adj_cell_connect(j)
                  stop
                endif  
              enddo           
            endif
          enddo

          if(.not.new_value) then
            do i=1,next_2
              color(tmp_list_2(i)) = my_color
            enddo
            next = 1
            old_next = 1
            break_ = .false.
            break_2 = .true.
            i = 1
            do while(break_2)
              if(color(i)==0) then
                my_color = my_color + 1
                tmp_list_2(1) = i
                next = 1
                next_2 = 1
                print*,"DEUZE je change :",i,my_color
                if(my_color>real_cell_nb) then
                  print*,"DEUZE Error ColoR :",my_color,real_cell_nb
                  stop
                endif
                break_2 = .false.
                break_ = .true.
              endif
              i = i + 1
              if(i>real_cell_nb) then
                break_2 = .false.
                print*,"DEUZE J ai termine !!! ",my_color
              endif
            enddo
          endif
        enddo

        print*,"DEUZE Nb color :",my_color 

        deallocate( counter )
        deallocate(new_cell_2_old_cell)
!22542       22023
        return
        end subroutine get_cell_connectivity
      end module get_cell_connectivity_mod































      
            !if(.false.) then
            !! ---------------------------------
            !if( my_cell_position(1)>1) then 
            !  break_ = .true.
            !  my_offset = (my_cell_position(2)-1)*domdec_struct%cell_x_nb(1) +  &
            !              (my_cell_position(3)-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
            !  ix = (my_cell_position(1)-1)
            !  my_neighbour = ix + my_offset
            !  if(empty_cell(my_neighbour)) then
            !    my_neighbour = -1
            !    break_ = .false.
            !    do while(break_)            
            !      my_neighbour = ix + my_offset
            !      print*," dow 1",ix
            !      if(empty_cell(my_neighbour)) then
            !        my_neighbour = -1
            !        ix = ix - 1
            !      else
            !        break_ = .false.
            !        my_neighbour = -1    
            !      endif
            !      if(ix<1) then
            !        break_ = .false.     
            !        my_neighbour = -1             
            !      endif
            !    enddo
            !  endif
            !  if(my_neighbour>0) then
            !    if(abs(my_position-my_neighbour)>domdec_struct%cell_x_nb(1)) then
            !     ! print*,"ErRoR 0",my_position,my_neighbour
            !     ! stop
            !    endif
            !    next = next + 1
            !    tmp_adj_cell_connect(next) = domdec_struct%old_cell_2_new_cell(my_neighbour)
            !    tmp_adj_cell_nb(i+1) = tmp_adj_cell_nb(i+1) + 1
            !  endif
            !endif
            !! i j k
            !! i-1 / i+1
            !! i j+1 / j-1 k
            !! i j k-1 / k+1
            !! i+1 j+1/j-1 k
            !! i-1 j+1 / j-1 k
            !! i+1 j k-1 / k+1
            !! i-1 j k-1 / k+1 
            !! ---------------------------------

            !! ---------------------------------
            !if( my_cell_position(1)<domdec_struct%cell_x_nb(1)) then 
            !  break_ = .true.
            !  my_offset = (my_cell_position(2)-1)*domdec_struct%cell_x_nb(1) + &
            !              (my_cell_position(3)-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
            !  ix = (my_cell_position(1)+1)
            !  my_neighbour = ix + my_offset
            !  if(empty_cell(my_neighbour)) then
            !    my_neighbour = -1
            !    break_ = .false.
            !    do while(break_)            
            !      my_neighbour = ix + my_offset
            !      print*," dow 2",ix
            !      if(empty_cell(my_neighbour)) then
            !        my_neighbour = -1
            !        ix = ix + 1
            !      else
            !        break_ = .false.
            !      endif
            !      if(ix>domdec_struct%cell_x_nb(1)) then
            !        break_ = .false.     
            !        my_neighbour = -1             
            !      endif
            !    enddo
            !  endif
            !  if(my_neighbour>0) then
            !    if(abs(my_position-my_neighbour)>domdec_struct%cell_x_nb(1)) then
            !     ! print*,"ErRoR 1",my_position,my_neighbour
            !      !stop
            !    endif
            !    next = next + 1
            !    tmp_adj_cell_connect(next) = domdec_struct%old_cell_2_new_cell(my_neighbour)
            !    tmp_adj_cell_nb(i+1) = tmp_adj_cell_nb(i+1) + 1
            !  endif
            !endif
            !! ---------------------------------

            !! ---------------------------------          
            !if( my_cell_position(2)>1) then 
            !  break_ = .true.
            !  my_offset = my_cell_position(1) + (my_cell_position(3)-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
            !  ijk = 1
            !  iy = (my_cell_position(2)-1-ijk)*domdec_struct%cell_x_nb(1)
            !  print*," Neighbour Y- 0",my_cell_position(1),my_cell_position(2)-ijk,my_cell_position(3)
            !  my_neighbour = iy + my_offset
            !  if(empty_cell(my_neighbour)) then
            !    my_neighbour = -1
            !    break_ = .false.
            !    print*," Neighbour Y- 1",my_cell_position(1),my_cell_position(2)-ijk,my_cell_position(3)
            !    do while(break_)            
            !      my_neighbour = iy + my_offset
            !      print*," dow 3",iy
            !      if(empty_cell(my_neighbour)) then
            !        my_neighbour = -1
            !        ijk = ijk + 1
            !        iy = (my_cell_position(2)-1-ijk)*domdec_struct%cell_x_nb(1) 
            !      else
            !        break_ = .false.
            !      endif
            !      if((my_cell_position(2)-ijk)<1) then
            !        break_ = .false.     
            !        my_neighbour = -1             
            !      endif
            !    enddo
            !  endif
            !  if(my_neighbour>0) then
            !    if(abs(my_position-my_neighbour)>domdec_struct%cell_x_nb(1)) then
            !      !print*,"ErRoR 3",my_position,my_neighbour
            !      !print*," ErRoR 3 cell pos ",my_cell_position(1:3)
            !      !stop
            !    endif
            !    next = next + 1
            !    tmp_adj_cell_connect(next) = domdec_struct%old_cell_2_new_cell(my_neighbour)
            !    tmp_adj_cell_nb(i+1) = tmp_adj_cell_nb(i+1) + 1
            !  endif
            !endif
            !! ---------------------------------

            !! ---------------------------------
            !if( my_cell_position(2)<domdec_struct%cell_x_nb(2)) then 
            !  break_ = .true.
            !  my_offset = my_cell_position(1) + (my_cell_position(3)-1)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
            !  ijk = 1
            !  iy = (my_cell_position(2)-1+ijk)*domdec_struct%cell_x_nb(1) 
            !  print*," Neighbour Y+ 0",my_cell_position(1),my_cell_position(2)+ijk,my_cell_position(3)
            !  my_neighbour = iy + my_offset
            !  if(empty_cell(my_neighbour)) then
            !    my_neighbour = -1
            !    break_ = .false.
            !    !print*," Neighbour Y+ 1",my_cell_position(1),my_cell_position(2)+ijk,my_cell_position(3)
            !    do while(break_)            
            !      my_neighbour = iy + my_offset
            !      print*," dow 4",iy
            !      if(empty_cell(my_neighbour)) then
            !        my_neighbour = -1
            !        ijk = ijk + 1
            !        iy = (my_cell_position(2)-1+ijk)*domdec_struct%cell_x_nb(1) 
            !      else
            !        break_ = .false.
            !      endif
            !      if((my_cell_position(2)+ijk)>domdec_struct%cell_x_nb(2)) then
            !        break_ = .false.     
            !        my_neighbour = -1             
            !      endif
            !    enddo
            !  endif
            !  if(my_neighbour>0) then
            !    if(abs(my_position-my_neighbour)>domdec_struct%cell_x_nb(1)) then                
            !      !print*,"ErRoR 4",my_position,my_neighbour
            !      !print*," bis ",my_cell_position(1),my_cell_position(2)+ijk,my_cell_position(3)
            !      !stop
            !    endif
            !    next = next + 1
            !    tmp_adj_cell_connect(next) = domdec_struct%old_cell_2_new_cell(my_neighbour)
            !    tmp_adj_cell_nb(i+1) = tmp_adj_cell_nb(i+1) + 1
            !  endif
            !endif
            !! ---------------------------------

            !! ---------------------------------
            !if( my_cell_position(3)>1) then 
            !  break_ = .true.
            !  my_offset = my_cell_position(1) + (my_cell_position(2)-1)*domdec_struct%cell_x_nb(1)
            !  ijk = 1
            !  iz = (my_cell_position(3)-1-ijk)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
            !  my_neighbour = iz + my_offset
            !  if(empty_cell(my_neighbour)) then
            !    my_neighbour = -1
            !    break_ = .false.
            !    do while(break_)            
            !      my_neighbour = iz + my_offset
            !      print*," dow 5",iz,my_neighbour,empty_cell(my_neighbour)
            !      if(empty_cell(my_neighbour)) then
            !        my_neighbour = -1
            !        ijk = ijk + 1
            !        iz = (my_cell_position(3)-1-ijk)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
            !      else
            !        break_ = .false.
            !      endif
            !      print*," dow 5 1",iz,(my_cell_position(3)-1-ijk)
            !      if((my_cell_position(3)-ijk)<1) then
            !        break_ = .false.     
            !        my_neighbour = -1             
            !      endif
            !    enddo
            !  endif
            !  if(my_neighbour>0) then
            !    if(abs(my_position-my_neighbour)>domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)) then
            ! !     print*,"ErRoR 5",my_position,my_neighbour
            !  !    stop
            !    endif
            !    next = next + 1              
            !    tmp_adj_cell_connect(next) = domdec_struct%old_cell_2_new_cell(my_neighbour)
            !    print*," ZZ - ",my_cell_position(1:3),my_neighbour,empty_cell(my_neighbour),next, &
            !      domdec_struct%old_cell_2_new_cell(my_neighbour)
            !    tmp_adj_cell_nb(i+1) = tmp_adj_cell_nb(i+1) + 1
            !  endif
            !endif
            !! ---------------------------------

            !! ---------------------------------
            !if( my_cell_position(3)<domdec_struct%cell_x_nb(3)) then 
            !  break_ = .true.
            !  my_offset = my_cell_position(1) + (my_cell_position(2)-1)*domdec_struct%cell_x_nb(1)
            !  ijk = 1
            !  iz = (my_cell_position(3)-1+ijk)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
            !  my_neighbour = iz + my_offset
            !  if(empty_cell(my_neighbour)) then
            !    my_neighbour = -1
            !    break_ = .false.
            !    do while(break_)            
            !      my_neighbour = iz + my_offset
            !      print*," dow 6",my_cell_position(1:3),my_neighbour,empty_cell(my_neighbour)
            !      if(empty_cell(my_neighbour)) then
            !        my_neighbour = -1
            !        ijk = ijk + 1
            !        iz = (my_cell_position(3)-1+ijk)*domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)
            !      else
            !        break_ = .false.
            !      endif
            !      print*," dow 6 1",iz + my_offset,my_cell_position(3)-1+ijk
            !      if((my_cell_position(3)+ijk)>domdec_struct%cell_x_nb(3)) then
            !        break_ = .false.     
            !        my_neighbour = -1             
            !      endif
            !    enddo
            !  endif
            !  if(my_neighbour>0) then
            !    if(abs(my_position-my_neighbour)>domdec_struct%cell_x_nb(1)*domdec_struct%cell_x_nb(2)) then
            ! !       print*,"ErRoR 6",my_position,my_neighbour
            ! !       stop
            !    endif
            !    next = next + 1
            !    tmp_adj_cell_connect(next) = domdec_struct%old_cell_2_new_cell(my_neighbour)
            !    tmp_adj_cell_nb(i+1) = tmp_adj_cell_nb(i+1) + 1
            !    print*," ZZ + ",my_cell_position(1:3),my_neighbour,empty_cell(my_neighbour), & 
            !    next,domdec_struct%old_cell_2_new_cell(my_neighbour)
            !  endif
            !endif  
            !! ---------------------------------      