!     Force calculation
!     Leonid Zhigilei, 2001
      SUBROUTINE Forces()
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'

        FD(1:3,1:NNF)=0.0d0
        POT(1:NNF)=0.0d0
        STEN(1:NNF,1:3,1:3)=0.0d0

        IF(KEYBS.EQ.4) THEN
           IF(KBOUND.EQ.2) THEN
              CALL F_EAM_br()
           ELSEIF(KBOUND.EQ.11) THEN
              CALL F_EAM_br_test()
           ELSE
              CALL F_EAM()
           ENDIF
           
        ELSEIF (KEYBS.EQ.5) THEN
           IF(LIDZ.EQ.0.AND.KBOUND.EQ.2) THEN
              IF(NTYPE.GE.2) THEN
                 CALL F_SW()
                 CALL F_EAM_br()
                 CALL F_pair()
              ELSE
                 CALL F_SW()
              ENDIF
           ELSE
              IF(NTYPE.GE.2) THEN
                 CALL F_SW()
                 CALL F_EAM()
                 CALL F_pair()
              ELSE
                 CALL F_SW()
              ENDIF
           ENDIF
           
        ELSEIF(KEYBS.EQ.6) THEN
           IF (mypid.EQ.0) print *,"Not implemented"
           CALL MPI_FINALIZE(info)
           STOP
           
        ELSEIF(KEYBS.EQ.0) THEN
           CALL F_pair()
           
        ELSEIF(KEYBS.EQ.7) THEN
           IF(NTYPE.EQ.1) THEN
              IF(LIDZ.EQ.0.AND.KBOUND.EQ.2) THEN
                 IF(mypid.eq.0) THEN
                    print *,"Not implemented in Forces"
                    STOP
                 ENDIF
              ELSE
                 CALL F_TF()
              ENDIF
           ELSEIF(NTYPE.EQ.2) THEN
              CALL F_TF()
              IF(LIDZ.EQ.0.AND.KBOUND.EQ.2) THEN
                 CALL F_EAM_br()
              ELSE
                 CALL F_EAM()
              ENDIF
              CALL F_pair()
           ELSE
              IF(mypid.eq.0) THEN
                 print *,"Not implemente in Forces"
                 STOP
              ENDIF
           ENDIF

        ELSEIF (KEYBS.EQ.8) THEN
           IF(NTYPE.EQ.1) THEN
              IF(LIDZ.EQ.0.AND.KBOUND.EQ.2) THEN
                 IF(mypid.eq.0) THEN
                    print *,"Not implemented in Forces"
                    STOP
                 ENDIF
              ELSE
                 CALL F_TFM()
              ENDIF
           ELSEIF(NTYPE.EQ.2) THEN
              CALL F_TFM()
              IF(LIDZ.EQ.0.AND.KBOUND.EQ.2) THEN
                 CALL F_EAM_br()
              ELSE
                 CALL F_EAM()
              ENDIF
              CALL F_pair()
           ENDIF
        ELSE
           IF(mypid.eq.0) THEN
              print *,"The potential is not implemented"
!              call Flushout(6)
           ENDIF
           CALL MPI_FINALIZE(info)
           STOP
        ENDIF

        
!       Additional forces from TTM
        IF (LGOT.EQ.2.AND.KEYDEP.NE.2) CALL F_TTM()
        
!       Additional forces from matter-substrate interactions, used here
!       specifically for the project with Chichkov
        IF(KEYDEP.NE.0) CALL F_SUB()

!       STEN is the static portion of stress tensor x atomic volume
        STEN(1:NN1,1:3,1:3)=STEN(1:NN1,1:3,1:3)*0.5d0

        RETURN
      END SUBROUTINE Forces
