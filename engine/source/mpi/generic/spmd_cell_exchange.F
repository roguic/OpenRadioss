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
      !||    spmd_cell_exchange         ../engine/source/mpi/generic/spmd_cell_exchange.F
      !||--- called by ------------------------------------------------------
      !||    inter_sort_07              ../engine/source/interfaces/int07/inter_sort_07.F
      !||--- calls      -----------------------------------------------------
      !||    ancmsg                     ../engine/source/output/message/message.F
      !||    arret                      ../engine/source/system/arret.F
      !||    fill_voxel_local_partial   ../engine/source/interfaces/intsort/fill_voxel.F90
      !||    fill_voxel_remote          ../engine/source/interfaces/intsort/fill_voxel.F90
      !||    startime                   ../engine/source/system/timer_mod.F90
      !||    stoptime                   ../engine/source/system/timer_mod.F90
      !||--- uses       -----------------------------------------------------
      !||    fill_voxel_mod             ../engine/source/interfaces/intsort/fill_voxel.F90
      !||    inter_sorting_mod          ../engine/share/modules/inter_sorting_mod.F
      !||    inter_struct_mod           ../engine/share/modules/inter_struct_mod.F
      !||    message_mod                ../engine/share/message_module/message_mod.F
      !||    multi_fvm_mod              ../common_source/modules/ale/multi_fvm_mod.F90
      !||    timer_mod                  ../engine/source/system/timer_mod.F90
      !||    tri7box                    ../engine/share/modules/tri7box.F
      !||====================================================================
      SUBROUTINE SPMD_CELL_EXCHANGE(TIMERS, NIN,ISENDTO,IRCVFROM,NSN,NSNR,IGAP,
     1                               IFQ,INACTI,NSNFIOLD,INTTH,ITYP,STFNS, NSV,
     2                               NRTM, X,
     2                               ITIED,NMN,INTER_STRUCT,SORT_COMM, GOT_PREVIEW)
!$COMMENT
!       SPMD_CELL_EXCHANGE description :
!       exchange of secondary node data (x, v, temp...)
!       SPMD_CELL_EXCHANGE organization :
!       proc P needs to :
!                  * send data if local NSN > 0 & remote NMN > 0 (--> SORT_COMM(NIN)%NB(P)>0)
!                  * rcv data if local NMN > 0 & remote NSN > 0 (--> given by SORT_COMM(NIN)%NBIRECV)
!$ENDCOMMENT
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
        USE TIMER_MOD
        USE FILL_VOXEL_MOD
        USE TRI7BOX
        USE MESSAGE_MOD
        USE MULTI_FVM_MOD
        USE INTER_SORTING_MOD
        USE INTER_STRUCT_MOD
        USE COMMAND_LINE_ARGS_MOD, ONLY : GOT_NEWSORT
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
#include      "task_c.inc"
#include      "timeri_c.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
        TYPE(TIMER_) :: TIMERS
        INTEGER NIN, IFQ, INACTI, IGAP,INTTH,NSN,NSNR,
     .        ITIED,
     .        NSNFIOLD(NSPMD),
     .        ISENDTO(NINTER+1,NSPMD+1), IRCVFROM(NINTER+1,NSPMD+1),
     .        ITYP
        INTEGER :: GOT_PREVIEW !< flag to indicate if -preview is available
        INTEGER :: NMN
        TYPE(inter_struct_type), DIMENSION(NINTER), INTENT(inout) :: INTER_STRUCT   !   structure for interface
        TYPE(sorting_comm_type), DIMENSION(NINTER), INTENT(inout) :: SORT_COMM   ! structure for interface sorting comm
        my_real :: stfns(nsn)
        INTEGER :: NSV(NSN)
        INTEGER :: NRTM
        my_real :: X(3,NUMNOD)
C-----------------------------------------------
C   L o c a l  V a r i a b l e s
C-----------------------------------------------
#ifdef MPI
      INTEGER MSGTYP,INFO,I,NOD, DT_CST, LOC_PROC,P,IDEB,
     .        SIZ,J, L, BUFSIZ, LEN, NB, IERROR1, IAD,
     .        STATUS(MPI_STATUS_SIZE),IERROR,REQ_SB(NSPMD),
     .        REQ_RB(NSPMD),KK,NBIRECV,IRINDEXI(NSPMD),
     .        REQ_RD(NSPMD),REQ_SD(NSPMD),
     .        REQ_RC(NSPMD),REQ_SC(NSPMD),
     .        INDEXI,ISINDEXI(NSPMD),
     .        MSGOFF, MSGOFF2, MSGOFF3, MSGOFF4, MSGOFF5,
     .        RSIZ, ISIZ, L2, REQ_RD2(NSPMD),
     .        LEN2, RSHIFT, ISHIFT, ND, JDEB, Q, NBB
      integer :: k
      INTEGER :: P_LOC
      INTEGER :: KEY,CODE
     
      DATA MSGOFF4/6026/ 
      DATA MSGOFF5/6027/ 
        
      my_real
     .        XMAXB,YMAXB,ZMAXB,XMINB,YMINB,ZMINB
      
      INTEGER, DIMENSION(:), ALLOCATABLE :: ITAGNSNFI  
      my_real, DIMENSION(:,:), ALLOCATABLE :: XTMP
      INTEGER :: ADRESS, LOCAL_RANK
      INTEGER :: SIZE_S,my_address,next
      INTEGER :: OFFSET(NSPMD)
      INTEGER :: DUMMY(1)
C-----------------------------------------------
C   S o u r c e  L i n e s
C-----------------------------------------------
C
C================================================================
C     tag des boites contenant des facettes
C     et creation des candidats
C================================================================
      LOC_PROC = ISPMD + 1
      ! save the old NSNFI values
      IF(INACTI==5.OR.INACTI==6.OR.INACTI==7.OR.IFQ>0
     .   .OR.ITIED/=0.OR.ITYP==23.OR.ITYP==24   
     .   .OR.ITYP==25) THEN
         DO P = 1, NSPMD
           NSNFIOLD(P) = INTER_STRUCT(NIN)%NSNFIOLD(P)
         END DO
      END IF
      NSNR = SORT_COMM(NIN)%NSNR
      OFFSET(1:NSPMD) = NSNR+1 


      IF(.NOT. (IRCVFROM(NIN,LOC_PROC)==0.AND.ISENDTO(NIN,LOC_PROC)==0))  THEN

      IF (IMONM > 0) CALL STARTIME(TIMERS,25)

      RSIZ = SORT_COMM(NIN)%RSIZ
      ISIZ = SORT_COMM(NIN)%ISIZ  

      ! ---------------------------------------
      ! prepare the rcv --> local number of main node > 0
      ! allocation of buffer + reception
      IF(IRCVFROM(NIN,LOC_PROC)/=0) THEN  !   local nmn>0

       IF(NSNR>0) THEN   !   nsn remote > 0 --> only on proc with nmn>0
        if(GOT_NEWSORT==0) then
          ALLOCATE(XREM(RSIZ,NSNR),STAT=IERROR)
          ALLOCATE(IREM(ISIZ,NSNR),STAT=IERROR)          
        else
          IERROR = 0
          do i=1,nspmd
            if(ispmd/=i-1) then
              allocate(sort_comm(nin)%data_proc(i)%irem(NSNFI(NIN)%P(i)*ISIZ+1))
              allocate(sort_comm(nin)%data_proc(i)%xrem(NSNFI(NIN)%P(i)*RSIZ+1))
              sort_comm(nin)%data_proc(i)%irem(:) = -2345678
              sort_comm(nin)%data_proc(i)%xrem(:) = -2345678
            endif
          enddo
        endif
        
          IF(IERROR/=0) THEN
            CALL ANCMSG(MSGID=20,ANMODE=ANINFO)
            CALL ARRET(2)
          ENDIF
          IDEB = 1
          DO L = 1, SORT_COMM(NIN)%NBIRECV
            P = SORT_COMM(NIN)%ISINDEXI(L)            
            OFFSET(L) = IDEB
            if(GOT_NEWSORT==0) then
              MSGTYP = MSGOFF4
              LEN = NSNFI(NIN)%P(P)*RSIZ
              CALL MPI_IRECV(
     1        XREM(1,IDEB),LEN,REAL,IT_SPMD(P),
     2        MSGTYP,SPMD_COMM_WORLD,REQ_RD(L),IERROR)
     
              MSGTYP = MSGOFF5
              LEN2 = NSNFI(NIN)%P(P)*ISIZ
              CALL MPI_IRECV(
     1        IREM(1,IDEB),LEN2,MPI_INTEGER,IT_SPMD(P),
     2        MSGTYP,SPMD_COMM_WORLD,REQ_RD2(L),IERROR)
              IDEB = IDEB + NSNFI(NIN)%P(P)     
            else
              !print*," R :",ncycle,nin,ispmd,P-1,NSNFI(NIN)%P(P)
              MSGTYP = MSGOFF4
              LEN = NSNFI(NIN)%P(P)*RSIZ + 1
              CALL MPI_IRECV(
     1        sort_comm(nin)%data_proc(p)%xrem,LEN,REAL,IT_SPMD(P),
     2        MSGTYP,SPMD_COMM_WORLD,REQ_RD(L),IERROR)
     
              MSGTYP = MSGOFF5      
              LEN2 = NSNFI(NIN)%P(P)*ISIZ + 1
              CALL MPI_IRECV(
     1        sort_comm(nin)%data_proc(p)%irem,LEN2,MPI_INTEGER,IT_SPMD(P),
     2        MSGTYP,SPMD_COMM_WORLD,REQ_RD2(L),IERROR)
              IDEB = IDEB + NSNFI(NIN)%P(P)   
            endif       
          ENDDO
        ENDIF
      ENDIF
      ! ---------------------------------------

      ! ---------------------------------------
      ! prepare the send --> local number of secondary node > 0 & remote number of main node > 0
      DO P=1,NSPMD
        IF(P/=LOC_PROC) THEN
            if(got_newsort/=0) then
              if(sort_comm(nin)%data_proc(p)%comm_needed) then
                MSGTYP = MSGOFF4
                SIZE_S = SORT_COMM(NIN)%NB(P) * SORT_COMM(NIN)%RSIZ + 1                
                !print*," S :",ncycle,nin,ispmd,P-1,SORT_COMM(NIN)%NB(P)
                CALL MPI_ISEND(
     1            SORT_COMM(NIN)%DATA_PROC(P)%RBUF,SIZE_S,REAL,IT_SPMD(P),MSGTYP,
     2            SPMD_COMM_WORLD,SORT_COMM(NIN)%REQ_SD2(P),ierror)
                MSGTYP = MSGOFF5              
                SIZE_S = SORT_COMM(NIN)%NB(P) * SORT_COMM(NIN)%ISIZ + 1              

                CALL MPI_ISEND(
     1              SORT_COMM(NIN)%DATA_PROC(P)%IBUF,SIZE_S,MPI_INTEGER,
     2              IT_SPMD(P),MSGTYP,
     3              SPMD_COMM_WORLD,SORT_COMM(NIN)%REQ_SD3(P),ierror)
              endif
            else
              IF(SORT_COMM(NIN)%NB(P)/=0 ) THEN
                MSGTYP = MSGOFF4

                SIZE_S = SORT_COMM(NIN)%NB(P) * SORT_COMM(NIN)%RSIZ
                CALL MPI_ISEND(
     1            SORT_COMM(NIN)%DATA_PROC(P)%RBUF,SIZE_S,REAL,IT_SPMD(P),MSGTYP,
     2            SPMD_COMM_WORLD,SORT_COMM(NIN)%REQ_SD2(P),ierror)
                MSGTYP = MSGOFF5       
                SIZE_S = SORT_COMM(NIN)%NB(P) * SORT_COMM(NIN)%ISIZ
                CALL MPI_ISEND(
     1            SORT_COMM(NIN)%DATA_PROC(P)%IBUF,SIZE_S,MPI_INTEGER,
     2            IT_SPMD(P),MSGTYP,
     3            SPMD_COMM_WORLD,SORT_COMM(NIN)%REQ_SD3(P),ierror)
              ENDIF
            endif
        ENDIF
      ENDDO


        ! -----------------------------
        if(GOT_PREVIEW == 1) THEN
        ! finish the last groups of secondary nodes, if any
            CALL FILL_VOXEL_LOCAL_PARTIAL(nsn,nsv,nsnr,nrtm,numnod,x,stfns,INTER_STRUCT(NIN),DUMMY,0) 
        ENDIF
 

      ! ---------------------------------------
      ! wait the rcv comm
      IF(IRCVFROM(NIN,LOC_PROC)/=0) THEN  !   nmn > 0
       IF(NSNR>0) THEN   !   nsnr>0 only on proc with nmn>0
        !print*," Wait R :",ncycle,nin,ispmd,SORT_COMM(NIN)%NBIRECV
         DO L = 1, SORT_COMM(NIN)%NBIRECV
          CALL MPI_WAITANY(SORT_COMM(NIN)%NBIRECV,REQ_RD,INDEXI,STATUS,IERROR) !XREM

         ! CALL MPI_WAITANY(SORT_COMM(NIN)%NBIRECV,REQ_RD2,INDEXI,STATUS,IERROR)

          if(GOT_PREVIEW==1) then
          call fill_voxel_remote(
     .                    OFFSET(INDEXI), 
     .                    OFFSET(INDEXI+1)-1, 
     .                    nsn,
     .                    nsnr,
     .                    inter_struct(nin)%nbx,
     .                    inter_struct(nin)%nby,
     .                    inter_struct(nin)%nbz,
     .                    size(XREM,1),
     .                    inter_struct(nin)%voxel,
     .                    inter_struct(nin)%next_nod,
     .                    inter_struct(nin)%size_node,
     .                    inter_struct(nin)%nb_voxel_on,
     .                    inter_struct(nin)%list_nb_voxel_on,
     .                    inter_struct(nin)%last_nod,
     .                    xrem,
     .                    inter_struct(nin)%box_limit_main)
          endif
          CALL MPI_WAIT(REQ_RD2(INDEXI),STATUS,IERROR) !IREM

    
         ENDDO
         !set specifics IREM and XREM indexes for INT24 sorting
         !IGAPXREMP = IREM(4,1)
         !I24XREMP  = IREM(5,1)
         !I24IREMP  = IREM(6,1)
         if(GOT_NEWSORT/=0) then
          SORT_COMM(NIN)%NSNR = 0
          do l=1,sort_comm(nin)%nbirecv
            P = SORT_COMM(NIN)%ISINDEXI(L)          
            nsnfi(nin)%p(p) = SORT_COMM(NIN)%DATA_PROC(P)%irem(1)
            SORT_COMM(NIN)%NSNR = SORT_COMM(NIN)%NSNR + nsnfi(nin)%p(p)          
          enddo
          nsnr = SORT_COMM(NIN)%NSNR
          ALLOCATE(XREM(RSIZ,NSNR),STAT=IERROR)
          ALLOCATE(IREM(ISIZ,NSNR),STAT=IERROR)
          !xrem(:,:) = -12345678
          !irem(:,:) = -12345678
          my_address = 0
          do l=1,sort_comm(nin)%nbirecv
            P = SORT_COMM(NIN)%ISINDEXI(L)
            next = 0
            do j=1,nsnfi(nin)%p(p)
              do k=1,RSIZ
                next = next + 1
          !      if(SORT_COMM(NIN)%DATA_PROC(P)%xrem(1+next)==-2345678) then
          !        !  pb -2   2   1  0         110          12           1
          !        print*," pb -2",nin,ispmd,p-1,1+next,nsnfi(nin)%p(p),l
          !        stop
          !      endif
                xrem(k,my_address+j) = SORT_COMM(NIN)%DATA_PROC(P)%xrem(1+next)
              enddo
            enddo
            next = 0
            do j=1,nsnfi(nin)%p(p)
              do k=1,ISIZ
                next = next + 1
          !      if(SORT_COMM(NIN)%DATA_PROC(P)%irem(1+next)==-2345678) then
          !        print*," pb -1",nin,ispmd,p-1,1+next,nsnfi(nin)%p(p),l
          !        stop
          !      endif
                irem(k,my_address+j) = SORT_COMM(NIN)%DATA_PROC(P)%irem(1+next)
              enddo
            enddo
            my_address = my_address + nsnfi(nin)%p(p)
          enddo  
!          ideb = 0
!          DO l = 1, SORT_COMM(NIN)%NBIRECV
!            P = SORT_COMM(NIN)%ISINDEXI(L)      
!
!          do j=1,nsnfi(nin)%p(p)
!            do k=1,RSIZ              
!              if(xrem(k,j+ideb)==-12345678) then
!                print*," pb 1",k,j,xrem(k,j)
!                stop
!              endif
!            enddo
!          enddo
!          do j=1,nsnfi(nin)%p(p)
!            do k=1,isiz
!              if(irem(k,j+ideb)==-12345678) then
!                print*," pb 2",k,j,irem(k,j)
!                stop
!              endif
!            enddo
!          enddo
!          ideb = ideb + nsnfi(nin)%p(p)
!          enddo
         endif           

        ENDIF
      ENDIF

      ! ---------------------------------------
      ! wait the send comm : only for proc with nsn>0 & nmn>0
      ! for proc with nsn>0 & nmn=0, wait is done after the sort
      IF(ISENDTO(NIN,LOC_PROC)/=0) THEN   !   nsn >0
        DO P = 1, NSPMD
          IF(IRCVFROM(NIN,P)/=0) THEN     !   nmn >0
            IF(P/=LOC_PROC) THEN
              IF(SORT_COMM(NIN)%NB(P)/=0) THEN
                IF(NMN/=0) THEN
                    CALL MPI_WAIT(SORT_COMM(NIN)%REQ_SD2(P),STATUS,IERROR) 
                    DEALLOCATE(SORT_COMM(NIN)%DATA_PROC(P)%RBUF)
                    !   can be moved if local nsn>0 & local nmn=0  
                    CALL MPI_WAIT(SORT_COMM(NIN)%REQ_SD3(P),STATUS,IERROR) 
                    DEALLOCATE(SORT_COMM(NIN)%DATA_PROC(P)%IBUF)         
                    SORT_COMM(NIN)%NB(P) = 0             
                ENDIF
              END IF
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      ! ---------------------------------------
      ELSE 

        if(GOT_PREVIEW == 1) THEN
        ! finish the last groups of secondary nodes, if any
            CALL FILL_VOXEL_LOCAL_PARTIAL(nsn,nsv,nsnr,nrtm,numnod,x,stfns,INTER_STRUCT(NIN),DUMMY,0) 
        ENDIF

       ENDIF
      IF (IMONM > 0) CALL STOPTIME(TIMERS,25)
#endif
      RETURN

      END SUBROUTINE SPMD_CELL_EXCHANGE
