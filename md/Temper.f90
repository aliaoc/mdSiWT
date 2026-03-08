!     Calculation of temperature and averaged energies
!     Energies in eV/particle, Temperature in K.
!     You can modify this to calculate temperature 
!     for atoms with KHIST(J)=2 only.
!     Pressure calculation added, April 2002
!     Leonid Zhigilei, 2000
      SUBROUTINE TEMPER(QINT,POTT,TEMPTR,TEMPTR1,TEMPTR2)
        INCLUDE 'common.h'
        
        POTT=0.0d0
        QINT=0.0d0
        QINT1 = 0.0d0
        QINT2 = 0.0d0
        TEMPTR = 0.0d0
        TEMPTR1 = 0.0d0
        TEMPTR2 = 0.0d0
        
        ent_loop: DO I=1,NN1
           POTT=POTT+POT(I)
           QINT = QINT + QIN(I)
           IF(KTYPE(I).EQ.1) THEN
              QINT1=QINT1+QIN(I)
           ELSEIF(KTYPE(I).EQ.2) THEN
              QINT2=QINT2+QIN(I)
           ELSE
              print *,"Temper.f90 is out of order",mypid
              STOP
           ENDIF
        END DO ent_loop
        
        QINT=QINT*ENUNIT
        QINT1 = QINT1*ENUNIT
        QINT2 = QINT2*ENUNIT
        POTT=POTT*ENUNIT

        TEMPTR=QINT*2.0d0/BK/DBLE(NDIM)/DBLE(NAN)
        TEMPTR1=QINT1*2.0d0/BK/DBLE(NDIM)/DBLE(NMETG)
        TEMPTR2=QINT2*2.0d0/BK/DBLE(NDIM)/DBLE(NWATG)
        
        RETURN
      END SUBROUTINE TEMPER
