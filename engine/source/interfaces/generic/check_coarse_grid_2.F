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
      !||    check_coarse_grid           ../engine/source/interfaces/generic/check_coarse_grid.F
      !||--- called by ------------------------------------------------------
      !||    spmd_coarse_cell_exchange   ../engine/source/mpi/interfaces/spmd_coarse_cell_exchange.F
      !||--- uses       -----------------------------------------------------
      !||    inter_sorting_mod           ../engine/share/modules/inter_sorting_mod.F
      !||====================================================================
        SUBROUTINE CHECK_COARSE_GRID_2(mode,NIN,M_S_COARSE_GRID,SORT_COMM,ITIED)
!$COMMENT
!       CHECK_COARSE_GRID description :
!       for a given interface, check if the correct processor needs to communicate with a remote processor
!
!       CHECK_COARSE_GRID organization :
!           loop over the processor P and :
!               * if P and the current processor ISPMD colored a common cell --> ISPMD and P need to communicate
!               * if P and the current processor ISPMD did not colored a common cell --> ISPMD and P don't need to communicate
!$ENDCOMMENT
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
        USE INTER_SORTING_MOD
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
#include      "comlock.inc"
#include "task_c.inc"
#include      "com01_c.inc"
C-----------------------------------------------
C   C o m m o n   B l o c k s
C-----------------------------------------------
#include      "com04_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
        INTEGER, INTENT(in) :: NIN,mode
        INTEGER, INTENT(in) :: ITIED
        INTEGER, DIMENSION(NB_BOX_COARSE_GRID**3), INTENT(in) :: M_S_COARSE_GRID
        TYPE(sorting_comm_type), DIMENSION(NINTER), INTENT(inout) :: SORT_COMM   ! structure for interface sorting comm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------       
        LOGICAL :: KEEP_THIS_PROC,BREAK_
        INTEGER :: I,J,IJK
        INTEGER :: ID_PROC,DISPL,my_size  
!   ----------------------------------------

        !   allocatio of KEEP_PROC : true if need to send/rcv to proc P
        IF(.NOT.ALLOCATED(SORT_COMM(NIN)%KEEP_PROC)) THEN
            ALLOCATE(SORT_COMM(NIN)%KEEP_PROC(SORT_COMM(NIN)%PROC_NUMBER) )
        ENDIF
        if(mode==1) then
          ! M local / S remote
          !   loop over the processor of the interface : keep P if at least 1 common cell between P and current proc
          DO I=1,SORT_COMM(NIN)%PROC_NUMBER
              ID_PROC = SORT_COMM(NIN)%PROC_LIST(I)
              DISPL = SORT_COMM(NIN)%address_s_coarse_cell(I)
              SORT_COMM(NIN)%KEEP_PROC(I) = .FALSE.
              KEEP_THIS_PROC = .FALSE.
              BREAK_ = .FALSE.
              J = 1
              !   no cell from P --> keep = false
              IF(SORT_COMM(NIN)%nb_s_coarse_cell(I)==0) BREAK_ = .TRUE.

              IF(.NOT.BREAK_) THEN
                  IF(SORT_COMM(NIN)%GLOBAL_COARSE_CELL(DISPL+NB_BOX_COARSE_GRID**3+1)==0) BREAK_ = .TRUE.
              ENDIF

              IF(.NOT.BREAK_) THEN
                  ! loop over the cell of P and check if 1 common cell
                  my_size = SORT_COMM(NIN)%GLOBAL_COARSE_CELL(DISPL+NB_BOX_COARSE_GRID**3+1)
                  if(my_size>4**3.or.my_size<0) then
                    print*,ncycle," PB -1 ",ispmd,id_proc,my_size
                  endif
                  DO J=1,my_size
                    if(DISPL+J>size(SORT_COMM(NIN)%GLOBAL_COARSE_CELL)) then
                      print*,ncycle," PB 0:'",ispmd,id_proc,DISPL,J,size(SORT_COMM(NIN)%GLOBAL_COARSE_CELL)
                    endif
                      IJK = SORT_COMM(NIN)%GLOBAL_COARSE_CELL(DISPL+J)
                      if(ijk>size(M_S_COARSE_GRID).or.ijk<=0) then
                        print*,ncycle," PB 1:'",ispmd,id_proc,ijk,size(M_S_COARSE_GRID)
                      endif
                      IF(M_S_COARSE_GRID(IJK)==1) KEEP_THIS_PROC = .TRUE.
                  ENDDO
              ENDIF
              IF(ITIED/=0) KEEP_THIS_PROC=.TRUE.
              IF(KEEP_THIS_PROC) SORT_COMM(NIN)%KEEP_PROC(I) = .TRUE.
          ENDDO
        else
          ! S local / M remote
          !   loop over the processor of the interface : keep P if at least 1 common cell between P and current proc
          DO I=1,SORT_COMM(NIN)%PROC_NUMBER
            ID_PROC = SORT_COMM(NIN)%PROC_LIST(I)
            DISPL = SORT_COMM(NIN)%address_m_coarse_cell(I)
            KEEP_THIS_PROC = .FALSE.
            BREAK_ = .FALSE.
            
            IF(SORT_COMM(NIN)%KEEP_PROC(I)) then
              KEEP_THIS_PROC = .TRUE.
              BREAK_ = .TRUE.
            ENDIF

            J = 1
            !   no cell from P --> keep = false
            IF(SORT_COMM(NIN)%nb_m_coarse_cell(I)==0) BREAK_ = .TRUE.

            IF(.NOT.BREAK_) THEN
              my_size = SORT_COMM(NIN)%GLOBAL_COARSE_CELL(DISPL+NB_BOX_COARSE_GRID**3+1)
              IF(my_size==0) BREAK_ = .TRUE.
            ENDIF

            IF(.NOT.BREAK_) THEN
                ! loop over the cell of P and check if 1 common cell
                my_size = SORT_COMM(NIN)%GLOBAL_COARSE_CELL(DISPL+NB_BOX_COARSE_GRID**3+1)
                if(my_size>4**3.or.my_size<0) then
                  print*,ncycle," PB -1 S ",ispmd,id_proc,my_size
                endif
                DO J=1,my_size
                  if(DISPL+J>size(SORT_COMM(NIN)%GLOBAL_COARSE_CELL)) then
                    print*,ncycle," PB 0 S'",ispmd,id_proc,DISPL,J,size(SORT_COMM(NIN)%GLOBAL_COARSE_CELL)
                  endif
                    IJK = SORT_COMM(NIN)%GLOBAL_COARSE_CELL(DISPL+J)
                    if(ijk>size(M_S_COARSE_GRID).or.ijk<=0) then
                      print*,ncycle," PB 1 S'",ispmd,id_proc,ijk,size(M_S_COARSE_GRID),my_size
                    endif
                    IF(M_S_COARSE_GRID(IJK)==1) KEEP_THIS_PROC = .TRUE.
                ENDDO
            ENDIF
            IF(ITIED/=0) KEEP_THIS_PROC=.TRUE.
            IF(KEEP_THIS_PROC) SORT_COMM(NIN)%KEEP_PROC(I) = .TRUE.
        ENDDO
        endif
        !   ------------------
        RETURN
        END SUBROUTINE CHECK_COARSE_GRID_2
