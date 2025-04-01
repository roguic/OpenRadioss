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
        SUBROUTINE CHECK_COARSE_GRID_2(NIN,M_S_COARSE_GRID,remote_M_S_COARSE_GRID,SORT_COMM,ITIED,index_proc)
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
        INTEGER, INTENT(in) :: NIN
        INTEGER, INTENT(in) :: ITIED,index_proc
        INTEGER, DIMENSION(NB_BOX_COARSE_GRID**3), INTENT(in) :: M_S_COARSE_GRID
        INTEGER, DIMENSION(NB_BOX_COARSE_GRID**3+1), INTENT(in) :: remote_M_S_COARSE_GRID
        TYPE(sorting_comm_type), DIMENSION(NINTER), INTENT(inout) :: SORT_COMM   ! structure for interface sorting comm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------       
        LOGICAL :: KEEP_THIS_PROC,BREAK_
        INTEGER :: J,IJK
        INTEGER :: ID_PROC,DISPL,my_size  
!   ----------------------------------------

        ! M local / S remote
        !   loop over the processor of the interface : keep P if at least 1 common cell between P and current proc
        KEEP_THIS_PROC = .FALSE.
        BREAK_ = .FALSE.       

        IF(.NOT.BREAK_) THEN
          IF(remote_M_S_COARSE_GRID(NB_BOX_COARSE_GRID**3+1)==0) BREAK_ = .TRUE.
        ENDIF

        IF(.NOT.BREAK_) THEN
          ! loop over the cell of P and check if 1 common cell
          my_size = remote_M_S_COARSE_GRID(NB_BOX_COARSE_GRID**3+1)
          if(my_size>4**3.or.my_size<0) then
                    print*,ncycle," PB -1 ",ispmd,id_proc,my_size
          endif
          DO J=1,my_size
            IJK = remote_M_S_COARSE_GRID(J)
            if(ijk<=0.or.ijk>size(M_S_COARSE_GRID)) then
              print*," pb 1!",ispmd,ijk,size(M_S_COARSE_GRID)
            endif
            IF(M_S_COARSE_GRID(IJK)/=0) KEEP_THIS_PROC = .TRUE.
          ENDDO
        ENDIF
        IF(ITIED/=0) KEEP_THIS_PROC=.TRUE.
        IF(KEEP_THIS_PROC) SORT_COMM(NIN)%KEEP_PROC(index_proc) = .TRUE.
        !   ------------------
        RETURN
        END SUBROUTINE CHECK_COARSE_GRID_2
