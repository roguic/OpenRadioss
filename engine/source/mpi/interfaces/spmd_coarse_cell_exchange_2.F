Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2025 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      !||====================================================================
      !||    spmd_coarse_cell_exchange   ../engine/source/mpi/interfaces/spmd_coarse_cell_exchange.F
      !||--- called by ------------------------------------------------------
      !||    inter_prepare_sort          ../engine/source/interfaces/generic/inter_prepare_sort.F
      !||--- calls      -----------------------------------------------------
      !||    check_coarse_grid           ../engine/source/interfaces/generic/check_coarse_grid.F
      !||    spmd_ialltoallv_int         ../engine/source/mpi/generic/spmd_ialltoallv_int.F
      !||--- uses       -----------------------------------------------------
      !||    inter_sorting_mod           ../engine/share/modules/inter_sorting_mod.F
      !||====================================================================
        SUBROUTINE SPMD_COARSE_CELL_EXCHANGE_2(NB_INTER_SORTED,LIST_INTER_SORTED,IRECVFROM,ISENDTO,MODE,
     .              IPARI,SORT_COMM)
!$COMMENT
!       SPMD_COARSE_CELL_EXCHANGE description :
!       for large interfaces : communication of coarse cells with alltoll mpi comm 
!       and check if 2 procs need to echange data
!
!       SPMD_COARSE_CELL_EXCHANGE organization :
!           First part MODE=1 : alltoall comm --> exchange of coarse cell
!           Second part MODE=2 : wait & check if 2 processors for a given interface need to echange data
!$ENDCOMMENT
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
        USE INTER_SORTING_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   M e s s a g e   P a s s i n g
C-----------------------------------------------
#include "spmd.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "param_c.inc"
#include      "task_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
        INTEGER, INTENT(in) :: NB_INTER_SORTED        !   number of interfaces that need to be sorted
        INTEGER, DIMENSION(NB_INTER_SORTED), INTENT(in) :: LIST_INTER_SORTED   !   list of interfaces that need to be sorted
        INTEGER, INTENT(in) :: MODE ! mode : 1 --> end/rcv / 2 --> wait + computation
        INTEGER, DIMENSION(NINTER+1,NSPMD+1), INTENT(in) :: ISENDTO,IRECVFROM ! array for S and R : isendto = nsn ; IRECVFROM = nmn
        INTEGER, DIMENSION(NPARI,NINTER), INTENT(in) ::  IPARI !   interface data
        TYPE(sorting_comm_type), DIMENSION(NINTER), INTENT(inout) :: SORT_COMM   ! structure for interface sorting comm
        !INTEGER, INTENT(inout) :: NB_REQUEST_COARSE_CELL ! number of request
        !INTEGER, DIMENSION(NB_INTER_SORTED), INTENT(inout) :: ARRAY_REQUEST_COARSE_CELL ! array of request
        !INTEGER, DIMENSION(NB_INTER_SORTED), INTENT(inout) :: LIST_INTER_COARSE_CELL    ! list of interface
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
#ifdef MPI        
        INTEGER :: KK,NIN,I,J
        INTEGER :: P,P_LOC,LOCAL_RANK
        INTEGER :: MY_SIZE,OLD_POINTER
        INTEGER :: ADRESS,SHIFT_

        INTEGER IERROR1,STATUS(MPI_STATUS_SIZE),IERROR
        INTEGER :: SIZE_CELL_LIST,TOTAL_RCV_SIZE,TOTAL_SEND_SIZE
        INTEGER :: LOC_PROC,ID_PROC
        INTEGER :: COUNT_COMM_SIZE_CELL,ID_COMM
        INTEGER :: ITIED
        integer, dimension(:), allocatable :: s_buffer
        integer :: my_address,my_value,next,ix,iy,iz,index_r
        integer, parameter :: tag = 13012
!   ----------------------------------------
        LOC_PROC = ISPMD + 1
        !   -------------------------
        !   MODE=1 : alltoall comm --> exchange of coarse cell
        IF(MODE==1) THEN
          DO KK=1,NB_INTER_SORTED
            NIN = LIST_INTER_SORTED(KK)
            IF(SORT_COMM(NIN)%PROC_NUMBER>NSPMD/2) THEN

              IF(.NOT.ALLOCATED(SORT_COMM(NIN)%KEEP_PROC)) THEN
                ALLOCATE(SORT_COMM(NIN)%KEEP_PROC(SORT_COMM(NIN)%PROC_NUMBER) )
              ENDIF
              SORT_COMM(NIN)%KEEP_PROC(1:SORT_COMM(NIN)%PROC_NUMBER) = .FALSE.
              IF(IRECVFROM(NIN,LOC_PROC)==0.AND.ISENDTO(NIN,LOC_PROC)==0) CYCLE

              SORT_COMM(NIN)%nb_coarse_cell_r_request = 0
              IF(.NOT.ALLOCATED(SORT_COMM(NIN)%coarse_cell_r_request)) THEN
                MY_SIZE = SORT_COMM(NIN)%PROC_NUMBER
                ALLOCATE(SORT_COMM(NIN)%coarse_cell_r_request(MY_SIZE))
              ENDIF

              IF(.NOT.ALLOCATED(SORT_COMM(NIN)%coarse_cell_r_index)) THEN
                  MY_SIZE = SORT_COMM(NIN)%PROC_NUMBER
                  ALLOCATE(SORT_COMM(NIN)%coarse_cell_r_index(MY_SIZE))
              ENDIF

              SORT_COMM(NIN)%nb_coarse_cell_s_request = 0
              IF(.NOT.ALLOCATED(SORT_COMM(NIN)%coarse_cell_s_request)) THEN
                    MY_SIZE = SORT_COMM(NIN)%PROC_NUMBER
                    ALLOCATE(SORT_COMM(NIN)%coarse_cell_s_request(MY_SIZE))
              ENDIF
              IF(.NOT.ALLOCATED(SORT_COMM(NIN)%coarse_cell_s_index)) THEN
                  MY_SIZE = SORT_COMM(NIN)%PROC_NUMBER
                 ALLOCATE(SORT_COMM(NIN)%coarse_cell_s_index(MY_SIZE))
              ENDIF

              TOTAL_SEND_SIZE = 2*(NB_BOX_COARSE_GRID**3 + 1)
              allocate(s_buffer(TOTAL_SEND_SIZE))
              s_buffer(:) = 0
              my_address = 0
              ! Initialization of s_buffer with S voxel
              IF(ISENDTO(NIN,LOC_PROC)>0) then
                  s_buffer(1:NB_BOX_COARSE_GRID**3 + 1) = SORT_COMM(NIN)%COARSE_GRID(1:NB_BOX_COARSE_GRID**3 + 1)
                  my_address = my_address + NB_BOX_COARSE_GRID**3 + 1

                  do iz=1,s_buffer(NB_BOX_COARSE_GRID**3 + 1)
                    if(s_buffer(iz)<=0) then
                      print*," Pb EnVoI s_BuFfEr S ",ispmd,iz,s_buffer(my_address+iz)
                    endif
                  enddo
                  !print*," Size S :",ispmd,nin,s_buffer(NB_BOX_COARSE_GRID**3 + 1)
              endif
              next = 0               
                
              ! Initialization of s_buffer with M voxel
              IF(IRECVFROM(NIN,LOC_PROC)>0) then
                  do iz=1,NB_BOX_COARSE_GRID                    
                    do iy=1,NB_BOX_COARSE_GRID
                      do ix=1,NB_BOX_COARSE_GRID
                        if(sort_comm(nin)%main_coarse_grid(ix,iy,iz)/=0) then
                          my_value = IX+(IY-1)*NB_BOX_COARSE_GRID+(IZ-1)*NB_BOX_COARSE_GRID**2
                          next = next + 1
                          s_buffer(my_address+next) = my_value
                        endif
                      enddo
                    enddo
                  enddo
                  s_buffer(my_address+NB_BOX_COARSE_GRID**3 + 1) = next
                  !print*," Size M :",ispmd,nin,s_buffer(my_address+NB_BOX_COARSE_GRID**3 + 1)
                  do iz=1,s_buffer(my_address+NB_BOX_COARSE_GRID**3 + 1)
                    if(s_buffer(my_address+iz)<=0) then
                      print*," Pb EnVoI s_BuFfEr M ",ispmd,my_address+iz,s_buffer(my_address+iz)
                    endif
                  enddo
              endif
!isendto = nsn
!IRECVFROM = nmn
              my_size = 2*(NB_BOX_COARSE_GRID**3 + 1)
              if(.not.allocated(sort_comm(nin)%data_proc)) then
                allocate(sort_comm(nin)%data_proc(nspmd))
              endif
              DO I=1,SORT_COMM(NIN)%PROC_NUMBER
                  ID_PROC = SORT_COMM(NIN)%PROC_LIST(I)
                  if(.not.allocated(sort_comm(nin)%data_proc(id_proc)%coarse_cell)) then
                    allocate( sort_comm(nin)%data_proc(id_proc)%coarse_cell(my_size) )
                  endif
                  sort_comm(nin)%nb_coarse_cell_r_request = sort_comm(nin)%nb_coarse_cell_r_request + 1
                  sort_comm(nin)%coarse_cell_r_index(sort_comm(nin)%nb_coarse_cell_r_request) = i
                  call mpi_irecv(sort_comm(nin)%data_proc(id_proc)%coarse_cell,my_size,MPI_INTEGER,it_spmd(id_proc), 
     .          tag,SPMD_COMM_WORLD,sort_comm(nin)%coarse_cell_r_request(sort_comm(nin)%nb_coarse_cell_r_request),ierror)
              ENDDO

              DO I=1,SORT_COMM(NIN)%PROC_NUMBER
                  ID_PROC = SORT_COMM(NIN)%PROC_LIST(I)
                  sort_comm(nin)%nb_coarse_cell_s_request = sort_comm(nin)%nb_coarse_cell_s_request + 1
                  sort_comm(nin)%coarse_cell_s_index(sort_comm(nin)%nb_coarse_cell_s_request) = i
                  call mpi_isend(s_buffer,my_size,MPI_INTEGER,it_spmd(id_proc), 
     .          tag,SPMD_COMM_WORLD,sort_comm(nin)%coarse_cell_s_request(sort_comm(nin)%nb_coarse_cell_s_request),ierror)
              ENDDO
              deallocate(s_buffer)
            ENDIF                
          ENDDO
        ENDIF
        
        !   -------------------------
        !   MODE=2 : - wait the previous comm
        !            - check if current proc and remote proc need to communicate
        IF(MODE==2) THEN
          DO KK=1,NB_INTER_SORTED
            NIN = LIST_INTER_SORTED(KK)
            ITIED = IPARI(85,NIN)
            IF(SORT_COMM(NIN)%PROC_NUMBER>NSPMD/2) THEN
              DO i=1,sort_comm(nin)%nb_coarse_cell_r_request
                  call mpi_waitany( sort_comm(nin)%nb_coarse_cell_r_request,
     .    sort_comm(nin)%coarse_cell_r_request,index_r,status,ierror)
                  index_r = sort_comm(nin)%coarse_cell_r_index(index_r)
                  ID_PROC =  SORT_COMM(NIN)%PROC_LIST(index_r)
                  
                  !CALL MPI_WAIT(ARRAY_REQUEST_COARSE_CELL(KK),STATUS,IERROR)  
                  ! M local / S remote
                  IF(IRECVFROM(NIN,LOC_PROC)>0.and.ISENDTO(NIN,ID_PROC)>0)
     .            CALL CHECK_COARSE_GRID_2(NIN,SORT_COMM(NIN)%MAIN_COARSE_GRID,
     .    sort_comm(nin)%data_proc(id_proc)%coarse_cell,SORT_COMM,ITIED,index_r)
                  ! S local / M remote                
                  my_address = NB_BOX_COARSE_GRID**3 + 1
                  IF(ISENDTO(NIN,LOC_PROC)>0.and.IRECVFROM(NIN,ID_PROC)>0) 
     .              CALL CHECK_COARSE_GRID_2(NIN,SORT_COMM(NIN)%SECOND_COARSE_GRID,
     .                    sort_comm(nin)%data_proc(id_proc)%coarse_cell(my_address+1),SORT_COMM,ITIED,index_r)
                if(itied/=0) SORT_COMM(NIN)%KEEP_PROC(index_r)=.true.
              ENDDO
              DO i=1,sort_comm(nin)%nb_coarse_cell_s_request
                call mpi_waitany( sort_comm(nin)%nb_coarse_cell_s_request,
     .          sort_comm(nin)%coarse_cell_s_request,index_r,status,ierror)
              ENDDO 
              !SORT_COMM(NIN)%KEEP_PROC(1:SORT_COMM(NIN)%PROC_NUMBER) = .true.
            ENDIF
          ENDDO        
        ENDIF
!   -------------------------

#endif
        END SUBROUTINE SPMD_COARSE_CELL_EXCHANGE_2
