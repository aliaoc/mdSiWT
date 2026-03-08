!     Evaluation of forces for the Stillinger Weber potential for Silicon
!     [F. H. Stillinger and T. A. Weber, Phys. Rev. B 31, 5262-5271 (1985)]
!     Avinash Dongare and Leonid Zhigilei, 2003
      SUBROUTINE F_SW()
        INCLUDE 'common.h' 
        REAL*8 LAM

!       Clear 3-body part of the forces
        FDSW(1:3,1:NNF) = 0.0d0
        STENSW(1:NNF,1:3,1:3) = 0.0d0

        outer_loop: DO I=1,NN1        !loop over ALL particles
           IF (NNG(I).gt.MAXNNB) THEN
              print *, ' TOO MANY NEIGHBOURS: ',NNG(I),' > ', maxnnb
              call flush(6)
              STOP
           ENDIF
           ktypei=KTYPE(I)

           IF(ktypei.EQ.KWT) CYCLE outer_loop

           inner_loop: DO KJ=1,NNG(I)  !loop over known possible neighbors of i
              J=NNNG(KJ,i)              !extract the atom number
              KTYPEJ = KTYPE(J)

              IF(ktypej.eq.KWT) CYCLE inner_loop
              
              DXIJ=XD(1,J)-XD(1,I)      ! distances Xij, Yij and Zij
              DYIJ=XD(2,J)-XD(2,I)
              DZIJ=XD(3,J)-XD(3,I)
              
              IJ = ijindex(KTYPEI, KTYPE(J))     ! potential type for the pair
              RRIJ=DXIJ*DXIJ+DYIJ*DYIJ+DZIJ*DZIJ
              RDIJ=SQRT(RRIJ)                    ! Distance Rij
                
              IF (RDIJ.GE.RM(IJ)) CYCLE inner_loop
             
!             We are using full neighbor list. 
!             Only do pair potential if j>i

              JGTI: IF (J.GT.i) THEN
                 P=RDIJ*DXR(IJ)
                 KF=P
                 P=P-KF

!                Potential energy (give half the potential to each particle)
                 PENER=0.5d0*(UT(IJ,KF)+(UT(IJ,KF+1)-UT(IJ,KF))*P)
                 POT(I)=POT(I)+PENER
                 POT(J)=POT(J)+PENER 
!                End of potential energy
       
                 P=FT(IJ,KF)+(FT(IJ,KF+1)-FT(IJ,KF))*P

                 P=P/RDIJ

                 DX=P*DXIJ                        ! x component of force
                 FD(1,I)=FD(1,I)-DX
                 FD(1,J)=FD(1,J)+DX
                 DZ=P*DZIJ                        ! z component of force
                 FD(3,I)=FD(3,I)-DZ
                 FD(3,J)=FD(3,J)+DZ
                 DY=P*DYIJ                        ! y component of force
                 FD(2,I)=FD(2,I)-DY
                 FD(2,J)=FD(2,J)+DY

!                Static Portion of Stress Tensor
                 STEN(I,1,1)=STEN(I,1,1)+P*DXIJ*DXIJ
                 STEN(I,2,2)=STEN(I,2,2)+P*DYIJ*DYIJ
                 STEN(I,3,3)=STEN(I,3,3)+P*DZIJ*DZIJ
                 
                 STEN(J,1,1)=STEN(J,1,1)+P*DXIJ*DXIJ
                 STEN(J,2,2)=STEN(J,2,2)+P*DYIJ*DYIJ
                 STEN(J,3,3)=STEN(J,3,3)+P*DZIJ*DZIJ
              END IF JGTI

!             Now we add the three body part due to the neighbors of i
!             Since this information will be used in MPI, we collect it
!             in an additional force FDSW array
              inner_loop2: DO KK=KJ+1,NNG(I) ! Loop over all neighbors of atom I
                 K=NNNG(KK,I)   
                 KTYPEK = KTYPE(K)

                 IF(KTYPEK.EQ.KWT) CYCLE inner_loop2
                 
                 DXIK=XD(1,K)-XD(1,I)
                 DYIK=XD(2,K)-XD(2,I)
                 DZIK=XD(3,K)-XD(3,I)
                 IK = ijindex(KTYPEI, KTYPE(K))
                 RRIK=DXIK*DXIK+DYIK*DYIK+DZIK*DZIK
                 RDIK=SQRT(RRIK)
                 IF (RDIK.GE.RM(IK)) CYCLE inner_loop2
   
!                Calculations of the distance between J and K atoms
                 DXJK=XD(1,K)-XD(1,J) ! with atom I. Hence K>J required
                 DYJK=XD(2,K)-XD(2,J)
                 DZJK=XD(3,K)-XD(3,J) ! distances Xik, Yik and Zik

                 RRJK = DXJK*DXJK + DYJK*DYJK + DZJK*DZJK
                 RDJK = SQRT(RRJK)
                 JK = ijindex(KTYPE(J),KTYPE(K))

!                This now forms the triplet of I, J and K atoms
                 eij = SW_SIG(IJ)/(RDIJ - RM(IJ))
                 eik = SW_SIG(IK)/(RDIK - RM(IK))

                 deij = -SW_SIG(IJ)/(RDIJ-RM(IJ))**2
                 deik = -SW_SIG(IK)/(RDIK-RM(IK))**2 

!                Use the cosine rule to calculate cos(theta)
                 costh = (RRIK + RRIJ - RRJK)/(2.0D+00*RDIJ*RDIK)
                
                 dcosdrij = 1.0D+00/RDIK - costh/RDIJ
                 dcosdrik = 1.0D+00/RDIJ - costh/RDIK
                 dcosdrjk = -RDJK/(RDIJ * RDIK)

                 cosTHIRD = costh + 1.0D+00/3.0D+00

                 IF(cosTHIRD.EQ.0.0D+00) THEN
                    TWcTH = 0.0D+00
                 ELSE
                    TWcTH = 2.0D+00/cosTHIRD
                 ENDIF
                 
                 EPS = SQRT(SW_EPS(IJ) * SW_EPS(IK))
                 LAM = (SW_LAM(KTYPEJ)*((SW_LAM(KTYPEI))**2.0d0)*SW_LAM(KTYPEK))**0.25d0
                 potegy= EPS*LAM*exp(SW_GAM(IJ)*eij+SW_GAM(IJ)*eik)*cosTHIRD**2

                 PEGY   = potegy/ENUNIT
!                Arbitarily giving all the three body energy to the central atom
                 POT(I)= POT(I)+PEGY

!                Calculate the Forces
!                Rij derivatives:
                 FFIJ = -potegy * (SW_GAM(IJ)*deij + TWcTH * dcosdrij)/RDIJ
                 FFIJ =FFIJ/ENUNIT               ! Program units

                 FDSW(1,I) = FDSW(1,I) - FFIJ*DXIJ
                 FDSW(1,J) = FDSW(1,J) + FFIJ*DXIJ

                 FDSW(3,I) = FDSW(3,I) - FFIJ*DZIJ
                 FDSW(3,J) = FDSW(3,J) + FFIJ*DZIJ

                 FDSW(2,I) = FDSW(2,I) - FFIJ*DYIJ
                 FDSW(2,J) = FDSW(2,J) + FFIJ*DYIJ

!                Rik derivatives:
                 FFIK = -potegy * (SW_GAM(IK)*deik + TWcTH * dcosdrik)/RDIK
                 FFIK =FFIK/ENUNIT               ! Program units
   
                 FDSW(1,I) = FDSW(1,I) - FFIK*DXIK
                 FDSW(1,K) = FDSW(1,K) + FFIK*DXIK

                 FDSW(3,I) = FDSW(3,I) - FFIK*DZIK
                 FDSW(3,K) = FDSW(3,K) + FFIK*DZIK

                 FDSW(2,I) = FDSW(2,I) - FFIK*DYIK
                 FDSW(2,K) = FDSW(2,K) + FFIK*DYIK

!                RJk derivatives:
                 FFJK =  -potegy * (TWcTH * dcosdrjk)/RDJK
                 FFJK =FFJK/ENUNIT               ! Program units
  
                 FDSW(1,J) = FDSW(1,J) - FFJK*DXJK
                 FDSW(1,K) = FDSW(1,K) + FFJK*DXJK

                 FDSW(3,J) = FDSW(3,J) - FFJK*DZJK
                 FDSW(3,K) = FDSW(3,K) + FFJK*DZJK

                 FDSW(2,J) = FDSW(2,J) - FFJK*DYJK
                 FDSW(2,K) = FDSW(2,K) + FFJK*DYJK

!                Static portion of stress tensor - 3 body 
                 STENSW(I,1,1)=STENSW(I,1,1)+FFIJ*DXIJ*DXIJ+FFIK*DXIK*DXIK
                 STENSW(I,2,2)=STENSW(I,2,2)+FFIJ*DYIJ*DYIJ+FFIK*DYIK*DYIK
                 STENSW(I,3,3)=STENSW(I,3,3)+FFIJ*DZIJ*DZIJ+FFIK*DZIK*DZIK
              
                 STENSW(J,1,1)=STENSW(J,1,1)+FFIJ*DXIJ*DXIJ+FFJK*DXJK*DXJK
                 STENSW(J,2,2)=STENSW(J,2,2)+FFIJ*DYIJ*DYIJ+FFJK*DYJK*DYJK
                 STENSW(J,3,3)=STENSW(J,3,3)+FFIJ*DZIJ*DZIJ+FFJK*DZJK*DZJK

                 STENSW(K,1,1)=STENSW(K,1,1)+FFIK*DXIK*DXIK+FFJK*DXJK*DXJK
                 STENSW(K,2,2)=STENSW(K,2,2)+FFIK*DYIK*DYIK+FFJK*DYJK*DYJK
                 STENSW(K,3,3)=STENSW(K,3,3)+FFIK*DZIK*DZIK+FFJK*DZJK*DZJK

              END DO inner_loop2
           END DO inner_loop

        END DO outer_loop

!       We shall now share the force contribution into the border layer of
!       each node because of the skin layer part taken from the related node
        CALL SHARE_SW()
        
!       Finally, add all forces up within the processor's border
        FD(1:3,1:NN1) = FD(1:3,1:NN1) + FDSW(1:3,1:NN1)
        STEN(1:NN1,1:3,1:3) = STEN(1:NN1,1:3,1:3) + STENSW(1:NN1,1:3,1:3)

        RETURN
      END SUBROUTINE F_SW
