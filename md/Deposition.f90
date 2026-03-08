!     Adding a cluster with velocity Vcl [m/s], z coordinate of Zshift. 
!     Cluster configuration is in the file connected to UNIT 8. 
!     Random number generator is used to choose X and Y positions of 
!     the cluster as well as the orientation of the cluster relatively to
!     the substrate.
      SUBROUTINE Deposition
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'

        IDDEP = 1
        
        IF(mypid.EQ.0) THEN
           NN1 = NN1 + 1
           KTYPE(NN1) = 2
           XD(1,NN1) = -130.0d0
           x_dep = -127.0d0
           y_dep = -127.0d0
           XD(2,NN1) = -127.0d0
           XD(3,NN1) = 5.0d0
           Q1D(1:2,NN1) = 0.0d0
           Q1D(3,NN1) = -5.0d0*DELTA
           KHIST(NN1) = 9
           Q2(NN1:NN1+2) = 0.0d0
           Q3(NN1:NN1+2) = 0.0d0
           Q4(NN1:NN1+2) = 0.0d0
           
           CALL LINKLIST()
           CALL NBLIST_MPI()
        ENDIF
        CALL MPI_REDUCE(NN1,NAN,1,MPI_INTEGER,MPI_SUM, &
             0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(NAN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        CALL MPI_BCAST(x_dep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(y_dep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        RETURN
      END SUBROUTINE Deposition
          
