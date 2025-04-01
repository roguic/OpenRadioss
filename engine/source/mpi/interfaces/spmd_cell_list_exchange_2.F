!Copyright>        OpenRadioss
!Copyright>        !Copyright (C) 1986-2025 Altair Engineering Inc.
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
      !||====================================================================
      !||    spmd_cell_list_exchange        ../engine/source/mpi/interfaces/spmd_cell_list_exchange.F
      !||--- called by ------------------------------------------------------
      !||    inter_prepare_sort             ../engine/source/interfaces/generic/inter_prepare_sort.F
      !||--- calls      -----------------------------------------------------
      !||    spmd_cell_size_exchange        ../engine/source/mpi/interfaces/spmd_cell_size_exchange.F
      !||    spmd_cell_size_exchange_init   ../engine/source/mpi/interfaces/spmd_cell_size_exchange_init.F
      !||    spmd_cell_size_post_rcv        ../engine/source/mpi/interfaces/spmd_cell_size_post_rcv.F
      !||    spmd_ialltoallv_int            ../engine/source/mpi/generic/spmd_ialltoallv_int.F
      !||--- uses       -----------------------------------------------------
      !||    inter_sorting_mod              ../engine/share/modules/inter_sorting_mod.F
      !||    inter_struct_mod               ../engine/share/modules/inter_struct_mod.F
      !||    tri7box                        ../engine/share/modules/tri7box.F
      !||====================================================================
        SUBROUTINE SPMD_CELL_LIST_EXCHANGE_2(IRCVFROM,ISENDTO,MODE,WEIGHT,IAD_ELEM, 
     .             FR_ELEM,X,V,MS,TEMP, 
     .             KINET,NODNX_SMS,ITAB,INTBUF_TAB,IPARI, 
     .             NIN,INTER_STRUCT,SORT_COMM,NODNX_SMS_SIZ,TEMP_SIZE,TIMERS)!, GOT_PREVIEW)
!$COMMENT
!       SPMD_CELL_LIST_EXCHANGE description :
!       exchange of remote cell and preparation of comm "send/rcv the data of remote nodes"
!       SPMD_CELL_LIST_EXCHANGE organization :
!       * mode = 1 : send & rcv the number of cell & prepare the array for the next comm
!       * mode = 2 : - wait the previous comm "mode 1"
!                    - send & rcv the cell 
!       * mode = 3 : - wait the previous comm "send/rcv the cell"
!                    - send/rcv the number of secondary nodes   
!                    - prepare the S buffer of secondary nodes (x, v, temp ...)
!       * mode = 4 : - wait the S comm "send the cells"
!                    - deallocation              
!$ENDCOMMENT
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
        USE INTBUFDEF_MOD  
        USE INTER_SORTING_MOD
        USE INTER_STRUCT_MOD
        USE TRI7BOX
        USE TIMER_MOD
!-----------------------------------------------
!   I m p l i c i t   T y p e s
!-----------------------------------------------
#include      "implicit_f.inc"
!-----------------------------------------------
!   M e s s a g e   P a s s i n g
!-----------------------------------------------
#include "spmd.inc"
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "param_c.inc"
#include      "task_c.inc"
#include      "tabsiz_c.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
          type(timer_) :: TIMERS
        INTEGER, INTENT(in) :: NIN
        INTEGER, INTENT(in) :: MODE
        INTEGER, INTENT(in) :: NODNX_SMS_SIZ ! size of NODNX_SMS
        INTEGER, INTENT(in) :: TEMP_SIZE     ! size of TEMP
        INTEGER, DIMENSION(NINTER+1,NSPMD+1), INTENT(in) :: ISENDTO,IRCVFROM
        INTEGER, DIMENSION(NPARI,NINTER), INTENT(in) ::  IPARI
        INTEGER, DIMENSION(NUMNOD), INTENT(inout) :: WEIGHT
        INTEGER, DIMENSION(2,NSPMD+1), INTENT(in) :: IAD_ELEM ! connectivity array iad(P+1)-iad(P) = nb of frontier node
        INTEGER, DIMENSION(SFR_ELEM), INTENT(in) :: FR_ELEM ! frontier node ID
        my_real, DIMENSION(3,NUMNOD), INTENT(in) :: X,V !   position & velocity
        my_real, DIMENSION(NUMNOD), INTENT(in) :: MS   !   mass 
        my_real, DIMENSION(TEMP_SIZE), INTENT(in) :: TEMP !   temperature
        INTEGER, DIMENSION(NUMNOD), INTENT(in) :: ITAB ! global node ID
        INTEGER, DIMENSION(NUMNOD), INTENT(in) :: KINET ! k energy 
        INTEGER, DIMENSION(NODNX_SMS_SIZ), INTENT(in) :: NODNX_SMS ! SMS array
        TYPE(INTBUF_STRUCT_), DIMENSION(NINTER), INTENT(in) :: INTBUF_TAB
        TYPE(inter_struct_type), DIMENSION(NINTER), INTENT(inout) :: INTER_STRUCT   !   structure for interface
        TYPE(sorting_comm_type), DIMENSION(NINTER), INTENT(inout) :: SORT_COMM   ! structure for interface sorting comm
        !INTEGER, INTENT(IN) :: GOT_PREVIEW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------- 
#ifdef MPI        
        LOGICAL :: IS_EXCHANGE_NEEDED
        INTEGER :: KK,I,J,k
        INTEGER :: P,P_LOC,LOCAL_RANK
        INTEGER :: SIZE_,OLD_POINTER
        INTEGER :: ADRESS,SHIFT_
        INTEGER :: INDEX_S,INDEX_R
        INTEGER :: ITIED

        INTEGER IERROR1,STATUS(MPI_STATUS_SIZE),IERROR
        INTEGER :: TOTAL_RCV_SIZE,TOTAL_SEND_SIZE
        INTEGER :: LOC_PROC,ID_PROC
        INTEGER :: COUNT_COMM_SIZE_CELL,ID_COMM
        INTEGER :: DIPLS_
        INTEGER, DIMENSION(:), ALLOCATABLE :: ARRAY_REQUEST
        LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: ALREADY_SEND
        INTEGER, DIMENSION(:), ALLOCATABLE :: INDEX_ALREADY_SEND

        INTEGER :: MSGTYP
        INTEGER :: MSGOFF
        DATA MSGOFF/13000/
        integer, parameter :: tag_1 = 13001
        integer, parameter :: tag_2 = 13002
        integer :: my_size
        integer :: nb_remote_m_cell, nb_remote_s_cell
        integer :: remote_s_node_nb
        integer :: my_shift,my_address,my_address_2,my_value
        integer :: cell_x_id,cell_y_id,cell_z_id,index_proc
        integer :: local_s_buffer(2)
        logical :: keep_comm
!   --------------------------------------------------------------------

        loc_proc = ispmd + 1
        if(mode==1) then
          do p=1,nspmd
            sort_comm(nin)%data_proc(p)%comm_needed = .false.
          enddo

          if(ircvfrom(nin,loc_proc)==0.and.isendto(nin,loc_proc)==0) return
          if(.not.allocated(sort_comm(nin)%nb_cell_proc)) then
            my_size = 2*sort_comm(nin)%proc_number
            allocate(sort_comm(nin)%nb_cell_proc(my_size))
          endif

          sort_comm(nin)%nb_request_r_nb_cell_proc = 0
          if(.not.allocated(sort_comm(nin)%request_r_nb_cell_proc)) then
            my_size = sort_comm(nin)%proc_number
            allocate(sort_comm(nin)%request_r_nb_cell_proc(my_size))
          endif

          sort_comm(nin)%nb_request_s_nb_cell_proc = 0
          if(.not.allocated(sort_comm(nin)%request_s_nb_cell_proc)) then
            my_size = sort_comm(nin)%proc_number
            allocate(sort_comm(nin)%request_s_nb_cell_proc(my_size))
          endif

          if(.not.allocated(sort_comm(nin)%index_rcv)) then
            my_size = sort_comm(nin)%proc_number
            allocate( sort_comm(nin)%index_rcv(my_size) )            
          endif

          if(.not.allocated(sort_comm(nin)%index_send)) then
            my_size = sort_comm(nin)%proc_number
            allocate( sort_comm(nin)%index_send(my_size) )            
          endif

          if(.not.allocated(sort_comm(nin)%index_rcv_2)) then
            my_size = sort_comm(nin)%proc_number
            allocate( sort_comm(nin)%index_rcv_2(my_size) )            
          endif

          if(.not.allocated(sort_comm(nin)%index_send_2)) then
            my_size = sort_comm(nin)%proc_number
            allocate( sort_comm(nin)%index_send_2(my_size) )            
          endif

          do i=1,sort_comm(nin)%proc_number
            id_proc = sort_comm(nin)%proc_list(i)
            my_shift = 2*sort_comm(nin)%proc_number
            my_address = (i-1)*2+1
            if(id_proc/=loc_proc) then
              keep_comm = .true.
              IF(SORT_COMM(NIN)%PROC_NUMBER>NSPMD/2) THEN
                IF(.NOT.SORT_COMM(NIN)%KEEP_PROC(i)) THEN
                  local_s_buffer(1:2) = 0
                  keep_comm = .false.
                ENDIF
              ENDIF
              if(keep_comm) then
                sort_comm(nin)%nb_request_r_nb_cell_proc = sort_comm(nin)%nb_request_r_nb_cell_proc + 1
                sort_comm(nin)%index_rcv(sort_comm(nin)%nb_request_r_nb_cell_proc) = i
                call mpi_irecv(sort_comm(nin)%nb_cell_proc(my_address),2,MPI_INTEGER,it_spmd(id_proc), 
     .           tag_1,SPMD_COMM_WORLD,sort_comm(nin)%request_r_nb_cell_proc(sort_comm(nin)%nb_request_r_nb_cell_proc),ierror)
              else
                sort_comm(nin)%nb_cell_proc(my_address:my_address+1) = 0
              endif
            endif
          enddo

          do i=1,sort_comm(nin)%proc_number
            id_proc = sort_comm(nin)%proc_list(i)
            if(id_proc/=loc_proc) then
              !MPI_ISEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, REQUEST, IERROR)
              local_s_buffer(1:2) = sort_comm(nin)%size_cell_list(1:2)
              keep_comm = .true.
              IF(SORT_COMM(NIN)%PROC_NUMBER>NSPMD/2) THEN
                IF(.NOT.SORT_COMM(NIN)%KEEP_PROC(i)) THEN
                  local_s_buffer(1:2) = 0
                  keep_comm = .false.
                ENDIF
              ENDIF
              if(keep_comm) then
                sort_comm(nin)%nb_request_s_nb_cell_proc = sort_comm(nin)%nb_request_s_nb_cell_proc + 1
                sort_comm(nin)%index_send(sort_comm(nin)%nb_request_s_nb_cell_proc) = i
                call mpi_isend(local_s_buffer,2,MPI_INTEGER,it_spmd(id_proc), 
     .              tag_1,SPMD_COMM_WORLD,sort_comm(nin)%request_s_nb_cell_proc(sort_comm(nin)%nb_request_s_nb_cell_proc),ierror)
              endif
            endif
          enddo
        endif

        if(mode==2) then
          sort_comm(nin)%nb_request_r_cell_proc = 0
          if(.not.allocated(sort_comm(nin)%request_r_cell_proc)) then
            my_size = sort_comm(nin)%proc_number
            allocate(sort_comm(nin)%request_r_cell_proc(my_size))
          endif
        
          sort_comm(nin)%nb_request_s_cell_proc = 0
          if(.not.allocated(sort_comm(nin)%request_s_cell_proc)) then
            my_size = sort_comm(nin)%proc_number
            allocate(sort_comm(nin)%request_s_cell_proc(my_size))
          endif

          do i=1,sort_comm(nin)%nb_request_r_nb_cell_proc
            call mpi_waitany( sort_comm(nin)%nb_request_r_nb_cell_proc,sort_comm(nin)%request_r_nb_cell_proc,index_r,status,ierror)
            j = sort_comm(nin)%index_rcv(index_r)
            id_proc = sort_comm(nin)%proc_list(j)
            my_address = (j-1)*2+1
            if(my_address+1>size(sort_comm(nin)%nb_cell_proc)) then
                print*,"erreur 1 :",ispmd,my_address,size(sort_comm(nin)%nb_cell_proc)
                stop
            endif
            nb_remote_m_cell = sort_comm(nin)%nb_cell_proc(my_address)
            nb_remote_s_cell = sort_comm(nin)%nb_cell_proc(my_address+1)
            my_size = nb_remote_m_cell + 2 * nb_remote_s_cell
            
            allocate( sort_comm(nin)%data_proc(id_proc)%remote_cell(my_size) )
            if(my_size>0) then
              sort_comm(nin)%nb_request_r_cell_proc = sort_comm(nin)%nb_request_r_cell_proc + 1
              sort_comm(nin)%index_rcv_2(sort_comm(nin)%nb_request_r_cell_proc) = j
              call mpi_irecv(sort_comm(nin)%data_proc(id_proc)%remote_cell,my_size,MPI_INTEGER,it_spmd(id_proc), 
     .          tag_2,SPMD_COMM_WORLD,sort_comm(nin)%request_r_cell_proc(sort_comm(nin)%nb_request_r_cell_proc),ierror)
            endif
          enddo

          do i=1,sort_comm(nin)%nb_request_s_nb_cell_proc
            call mpi_waitany( sort_comm(nin)%nb_request_s_nb_cell_proc,sort_comm(nin)%request_s_nb_cell_proc,index_s,status,ierror)
            j = sort_comm(nin)%index_send(index_s)
            id_proc = sort_comm(nin)%proc_list(j)

            my_size = sort_comm(nin)%size_cell_list(1) + 2*sort_comm(nin)%size_cell_list(2)
            IF(SORT_COMM(NIN)%PROC_NUMBER>NSPMD/2) THEN
              IF(.NOT.SORT_COMM(NIN)%KEEP_PROC(j)) THEN
                my_size = 0
              ENDIF
            ENDIF
            if(my_size>0) then
              sort_comm(nin)%nb_request_s_cell_proc = sort_comm(nin)%nb_request_s_cell_proc + 1
              sort_comm(nin)%index_send_2(sort_comm(nin)%nb_request_s_cell_proc) = j
              call mpi_isend(sort_comm(nin)%cell_list,my_size,MPI_INTEGER,it_spmd(id_proc), 
     .       tag_2,SPMD_COMM_WORLD,sort_comm(nin)%request_s_cell_proc(sort_comm(nin)%nb_request_s_cell_proc),ierror)
            endif
          enddo
        endif

        if(mode==3) then
          CALL SPMD_CELL_SIZE_EXCHANGE_INIT_2(1,IRCVFROM,ISENDTO,IPARI,NIN,INTER_STRUCT,SORT_COMM)  
          NSNFI(NIN)%P(1:NSPMD) = 0
          ALLOCATE( ALREADY_SEND(NB_CELL_X,NB_CELL_Y,NB_CELL_Z) )
          IF(ISENDTO(NIN,LOC_PROC)/=0) ALREADY_SEND(1:NB_CELL_X,1:NB_CELL_Y,1:NB_CELL_Z) = .FALSE.
          ALLOCATE( INDEX_ALREADY_SEND(NB_CELL_X*NB_CELL_Y*NB_CELL_Z) )
          do i=1,sort_comm(nin)%nb_request_r_cell_proc
            call mpi_waitany( sort_comm(nin)%nb_request_r_cell_proc,sort_comm(nin)%request_r_cell_proc,index_r,status,ierror)
            CALL STARTIME(TIMERS,157)
            index_proc = sort_comm(nin)%index_rcv_2(index_r)
            id_proc = sort_comm(nin)%proc_list(index_proc)
            my_address = (index_proc-1)*2+1

            nb_remote_m_cell = sort_comm(nin)%nb_cell_proc(my_address)
            nb_remote_s_cell = sort_comm(nin)%nb_cell_proc(my_address+1)
            my_address = nb_remote_m_cell
            my_address_2 = nb_remote_m_cell + nb_remote_s_cell
            if(ircvfrom(nin,loc_proc)>0) then ! check if I have M nodes
              do k=1,nb_remote_s_cell
                my_value = sort_comm(nin)%data_proc(id_proc)%remote_cell(k+my_address)
                if(my_value==0) then
                  print*," value ?",ispmd,k,my_value
                  stop
                endif
                if(my_value>0) then
                  cell_z_id = (my_value - mod(my_value,1000000) ) / 1000000
                  my_value = my_value - cell_z_id * 1000000
                  cell_y_id = (my_value - mod(my_value,1000) ) / 1000
                  my_value = my_value - cell_y_id * 1000
                  cell_x_id = my_value
                  if(cell_x_id>nb_cell_x) then
                    print*," erreur 2 !",cell_x_id,shape(sort_comm(nin)%m_cell_mark,1)
                    stop
                  endif
                  if(sort_comm(nin)%m_cell_mark(cell_x_id,cell_y_id,cell_z_id)/=0) then
                    remote_s_node_nb = sort_comm(nin)%data_proc(id_proc)%remote_cell(k+my_address_2)
                    nsnfi(nin)%p(id_proc) = nsnfi(nin)%p(id_proc) + remote_s_node_nb
                  endif
                else
                  remote_s_node_nb = sort_comm(nin)%data_proc(id_proc)%remote_cell(k+my_address_2)
                  nsnfi(nin)%p(id_proc) = nsnfi(nin)%p(id_proc) + abs(remote_s_node_nb)
                endif
              enddo
            endif
            CALL STOPTIME(TIMERS,157)
            if(isendto(nin,loc_proc)>0) then ! check if I have S nodes            
                is_exchange_needed=.true.
                if(id_proc==loc_proc) is_exchange_needed=.false.
                if(nb_remote_m_cell==0) is_exchange_needed=.false.
                ITIED = IPARI(85,NIN)
                IF(ITIED/=0.AND.IRCVFROM(NIN,id_proc)/=0)  IS_EXCHANGE_NEEDED = .TRUE.   ! itied option : force the exchange    
                if(is_exchange_needed) then
                  
                  CALL STARTIME(TIMERS,156)
                  !   prepare the S buffer of secondary nodes (x, v, temp ...)
                  CALL SPMD_CELL_SIZE_EXCHANGE_2(IRCVFROM,ISENDTO,WEIGHT,
     .              IAD_ELEM,FR_ELEM,X,V,MS,TEMP,KINET,NODNX_SMS,ITAB,INTBUF_TAB,IPARI,NIN,index_proc,
     .              ALREADY_SEND,INDEX_ALREADY_SEND,SORT_COMM,NODNX_SMS_SIZ,TEMP_SIZE)
                  CALL STOPTIME(TIMERS,156)
                endif
          
            endif
          enddo
          do i=1,sort_comm(nin)%nb_request_s_cell_proc
            call mpi_waitany( sort_comm(nin)%nb_request_s_cell_proc,sort_comm(nin)%request_s_cell_proc,index_r,status,ierror)
          enddo
          CALL SPMD_CELL_SIZE_EXCHANGE_INIT_2(2,IRCVFROM,ISENDTO,IPARI,NIN,INTER_STRUCT,SORT_COMM)  
          deallocate( sort_comm(nin)%m_cell_mark )
          DEALLOCATE( ALREADY_SEND )
          DEALLOCATE( INDEX_ALREADY_SEND )
          do i=1,nspmd
            if(allocated( sort_comm(nin)%data_proc(i)%remote_cell )) deallocate(sort_comm(nin)%data_proc(i)%remote_cell)
          enddo
        endif      

        !   ----------------------------
!IF(ISENDTO(NIN,ID_PROC)>0) THEN    !   nsn of proc ID_PROC >0
!IF(IRCVFROM(NIN,ID_PROC)>0) THEN    !   nmn of proc ID_PROC >0        
!   --------------------------------------------------------------------
#endif
        RETURN
        END SUBROUTINE SPMD_CELL_LIST_EXCHANGE_2
