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
!!Copyright>
!!!Copyright>
!!!Copyright>        Commercial Alternative: Altair Radioss Software
!!!Copyright>
!!!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!!!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!!!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      !||====================================================================
      !||    inter_color_voxel    ../engine/source/interfaces/generic/inter_color_voxel.F
      !||--- called by ------------------------------------------------------
      !||    inter_prepare_sort   ../engine/source/interfaces/generic/inter_prepare_sort.F
      !||--- calls      -----------------------------------------------------
      !||    inter_cell_color     ../engine/source/interfaces/generic/inter_cell_color.F
      !||    my_barrier           ../engine/source/system/machine.F
      !||--- uses       -----------------------------------------------------
      !||    inter_sorting_mod    ../engine/share/modules/inter_sorting_mod.F
      !||    inter_struct_mod     ../engine/share/modules/inter_struct_mod.F
      !||    multi_fvm_mod        ../common_source/modules/ale/multi_fvm_mod.F90
      !||====================================================================
        SUBROUTINE INTER_COLOR_VOXEL_2(ITASK,NB_INTER_SORTED,LIST_INTER_SORTED,IPARI,INTBUF_TAB, 
     .                                  X,INTER_STRUCT,SORT_COMM)
!$COMMENT
!       INTER_COLOR_VOXEL description :
!       color the fine cell & coarse cell with main nodes
!
!       INTER_COLOR_VOXEL organization :
!           loop over the interface and call of INTER_CELL_COLOR
!$ENDCOMMENT
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
        USE INTBUFDEF_MOD  
        USE MULTI_FVM_MOD
        USE INTER_SORTING_MOD
        USE INTER_STRUCT_MOD
!-----------------------------------------------
!   I m p l i c i t   T y p e s
!-----------------------------------------------
#include      "implicit_f.inc"
#include      "comlock.inc"
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
#include      "com01_c.inc"
#include      "com04_c.inc"
#include      "param_c.inc"
#include      "task_c.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        INTEGER, INTENT(in) :: ITASK ! thread ID
        INTEGER, INTENT(in) :: NB_INTER_SORTED        !   number of interfaces that need to be sorted
        INTEGER, DIMENSION(NB_INTER_SORTED), INTENT(in) :: LIST_INTER_SORTED   !   list of interfaces that need to be sorted
        INTEGER, DIMENSION(NPARI,NINTER), INTENT(in) :: IPARI    !   interface data
        TYPE(INTBUF_STRUCT_),DIMENSION(NINTER), INTENT(in) :: INTBUF_TAB    ! interface data
        my_real, DIMENSION(3*NUMNOD), INTENT(in) :: X            !   position
        TYPE(inter_struct_type), DIMENSION(NINTER), INTENT(inout) :: INTER_STRUCT   !   structure for interface
        TYPE(sorting_comm_type), DIMENSION(NINTER), INTENT(inout) :: SORT_COMM   ! structure for interface sorting comm
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: KK,I,J
        INTEGER :: NIN
        INTEGER :: IX,IY,IZ
        INTEGER :: SIZE_INDEX_CELL,TOTAL_NB_CELL
        INTEGER :: NRTM,NRTM_T
        INTEGER :: ADRESS,ESHIFT,SHIFT
        INTEGER, DIMENSION(:), ALLOCATABLE :: INDEX_CELL
        LOGICAL :: TYPE18
        INTEGER :: NTY,INACTI
        INTEGER :: VALUE
        my_real :: DGAPLOAD

        my_real :: TZINF,GAPMIN,GAPMAX,GAP,DRAD
        integer, dimension(:,:,:), allocatable :: cell_mark
        integer :: my_size
!   ----------------------------------------
        ! --------------------
        ! allocation of global omp array
        IF(ITASK==0) THEN
            ALLOCATE( NB_LOCAL_CELL(NTHREAD) ) 
            allocate(cell_mark(NB_CELL_X,NB_CELL_Y,NB_CELL_Z)) 
            cell_mark(1:NB_CELL_X,1:NB_CELL_Y,1:NB_CELL_Z) = 0
        ENDIF
        CALL MY_BARRIER()
        ! --------------------
        ! allocation of local omp array
        SIZE_INDEX_CELL = NB_CELL_X*NB_CELL_Y*NB_CELL_Z
        ALLOCATE(INDEX_CELL(SIZE_INDEX_CELL))
        NB_LOCAL_CELL(ITASK+1) = 0

        ! --------------------
        ! loop over the interface
        DO KK=1,NB_INTER_SORTED
            NIN = LIST_INTER_SORTED(KK)
            NRTM = IPARI(4,NIN)
            NRTM_T = NRTM/NTHREAD
            ESHIFT = ITASK*NRTM_T
            ADRESS = 1 + ITASK*(IPARI(4,NIN)/NTHREAD)
            IF(ITASK==NTHREAD-1) NRTM_T= NRTM - ADRESS + 1
            NB_LOCAL_CELL(ITASK+1) = 0
            
            TZINF = INTBUF_TAB(NIN)%VARIABLES(tzinf_index)
            GAP   =INTBUF_TAB(NIN)%VARIABLES(gap_index)
            GAPMIN=INTBUF_TAB(NIN)%VARIABLES(gapmin_index)
            GAPMAX=INTBUF_TAB(NIN)%VARIABLES(gapmax_index)
            DRAD = ZERO
            IF(IPARI(7,NIN)==7) DRAD =INTBUF_TAB(NIN)%VARIABLES(drad_index)
            DGAPLOAD = INTBUF_TAB(NIN)%VARIABLES(bgapemx_index)

            NTY = IPARI(7,NIN)
            INACTI = IPARI(22,NIN)
            TYPE18=.FALSE.
            IF(NTY==7 .AND. INACTI==7)TYPE18=.TRUE.

            IF(ITASK==0) THEN
                COARSE_GRID = .FALSE.
                IF(SORT_COMM(NIN)%PROC_NUMBER>NSPMD/2) THEN
                    COARSE_GRID = .TRUE.
                    IF(.NOT.ALLOCATED(SORT_COMM(NIN)%MAIN_COARSE_GRID) ) THEN
                        ALLOCATE(SORT_COMM(NIN)%MAIN_COARSE_GRID( 
     .                        NB_BOX_COARSE_GRID,NB_BOX_COARSE_GRID,NB_BOX_COARSE_GRID ) )
                   ENDIF
                   SORT_COMM(NIN)%MAIN_COARSE_GRID(:,:,:) = 0
               ENDIF
                if(allocated(sort_comm(nin)%m_cell_mark)) then
                  deallocate( sort_comm(nin)%m_cell_mark )
                endif
                allocate(sort_comm(nin)%m_cell_mark(NB_CELL_X,NB_CELL_Y,NB_CELL_Z))
                sort_comm(nin)%m_cell_mark(1:NB_CELL_X,1:NB_CELL_Y,1:NB_CELL_Z) = 0
            ENDIF

            CALL MY_BARRIER()
            ! --------------------  
            ! cell coloration by omp threads  
            CALL INTER_CELL_COLOR_2( X,IPARI(21,NIN)   ,NRTM_T  ,INTBUF_TAB(NIN)%STFM(1+ESHIFT)   , 
     .                 TZINF,INTER_STRUCT(NIN)%CURV_MAX(ADRESS), 
     .                 GAPMIN ,GAPMAX,INTBUF_TAB(NIN)%GAP_M(1+ESHIFT) , 
     .                 INTBUF_TAB(NIN)%IRECTM(1+4*ESHIFT),GAP,INTBUF_TAB(NIN)%VARIABLES(bgapsmx_index),DRAD, 
     .                 NB_LOCAL_CELL(ITASK+1),SIZE_INDEX_CELL,INDEX_CELL, 
     .                 COARSE_GRID,SORT_COMM(NIN)%MAIN_COARSE_GRID,DGAPLOAD,cell_mark)

            CALL MY_BARRIER()
            ! --------------------   
            ! reduction of number of colored cells
            IF(ITASK==0) THEN
                TOTAL_NB_CELL = 0
                DO I=1,NTHREAD
                    TOTAL_NB_CELL = TOTAL_NB_CELL + NB_LOCAL_CELL(I)
                ENDDO
                SORT_COMM(NIN)%SIZE_CELL_LIST(1) = TOTAL_NB_CELL
                SORT_COMM(NIN)%SIZE_CELL_LIST(2) = sort_comm(nin)%nb_index_to_s_node_nb_cell
                my_size = SORT_COMM(NIN)%SIZE_CELL_LIST(1)+2*SORT_COMM(NIN)%SIZE_CELL_LIST(2)
                ALLOCATE( SORT_COMM(NIN)%CELL_LIST(my_size) )
            ENDIF
            ! --------------------   
            ! flush global array cell_mark to true for the next interface
            DO I=1,NB_LOCAL_CELL(ITASK+1)

                VALUE = INDEX_CELL(I)
                IZ = ( VALUE - MOD(VALUE,1000000) ) / 1000000
                VALUE = VALUE - IZ * 1000000
                IY = ( VALUE - MOD(VALUE,1000) ) / 1000
                VALUE = VALUE - IY * 1000
                IX = VALUE
                sort_comm(nin)%m_cell_mark(ix,iy,iz) = 1
                cell_mark(IX,IY,IZ) = 0
            ENDDO
            CALL MY_BARRIER()
            ! --------------------   

            ! --------------------   
            ! reduction of colored cells
            SHIFT = 0
            IF(ITASK>0) THEN
                DO J=1,ITASK
                    SHIFT = SHIFT + NB_LOCAL_CELL(J)
                ENDDO
            ENDIF
            SORT_COMM(NIN)%CELL_LIST(1+SHIFT:NB_LOCAL_CELL(ITASK+1)+SHIFT) = INDEX_CELL(1:NB_LOCAL_CELL(ITASK+1))
            if(itask==0) THEN
              shift = SORT_COMM(NIN)%SIZE_CELL_LIST(1)
              SORT_COMM(NIN)%CELL_LIST(shift+1:shift+SORT_COMM(NIN)%SIZE_CELL_LIST(2)) =  
     .         sort_comm(nin)%index_to_s_node_nb_cell(1:SORT_COMM(NIN)%SIZE_CELL_LIST(2))
              shift = SORT_COMM(NIN)%SIZE_CELL_LIST(1) + SORT_COMM(NIN)%SIZE_CELL_LIST(2)
              do j=1,SORT_COMM(NIN)%SIZE_CELL_LIST(2)
                VALUE = abs(sort_comm(nin)%index_to_s_node_nb_cell(j))
                IZ = ( VALUE - MOD(VALUE,1000000) ) / 1000000
                VALUE = VALUE - IZ * 1000000
                IY = ( VALUE - MOD(VALUE,1000) ) / 1000
                VALUE = VALUE - IY * 1000
                IX = VALUE
                SORT_COMM(NIN)%CELL_LIST(shift+j) = sort_comm(nin)%s_node_nb_cell(ix,iy,iz)
              enddo
            endif
            CALL MY_BARRIER()
            ! --------------------   
        ENDDO

        CALL MY_BARRIER()
        ! --------------------
        ! deallocation
        IF(ITASK==0) THEN
            DEALLOCATE(cell_mark)
            DEALLOCATE( NB_LOCAL_CELL ) 
        ENDIF
        DEALLOCATE(INDEX_CELL)
        ! --------------------
        RETURN
        END SUBROUTINE INTER_COLOR_VOXEL_2
