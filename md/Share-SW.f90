!     This subroutine is responsible for sharing of the forces information
!     for particles that are in the skin of the current preocessor in the
!     case of SW potential
      SUBROUTINE SHARE_SW()
        INCLUDE 'common.h'
        INCLUDE 'mpif.h'
        INTEGER status(MPI_STATUS_SIZE)
        REAL*8, DIMENSION(:), ALLOCATABLE :: smsgbuf,rmsgbuf

        ALLOCATE(smsgbuf(MSGSWP),rmsgbuf(MSGSWP))
        
        ixoffset = mynodex*nlcx
        iyoffset = mynodey*nlcy
        izoffset = mynodez*nlcz

        neigh6: DO nbb = 1,26

           IF(nbb.le.18) THEN
              nb = longarray(nbb)
              nodecheck = longnode
           ELSEIF(nbb.gt.24)THEN
              nb = korarray(nbb-24)
              nodecheck = kornode
           ELSE
              nb = midlarray(nbb-18)
              nodecheck = midlnode
           ENDIF

           l = 1

!          Prevent communication based on PBC
           IF(ncom(nb).EQ.1) THEN

!             Loop through skin cells
              DO icn = 1,nscell(nb)
                 ic = kscell(nb,icn)
                 i = ltop(ic)
              
                 IF(i.gt.0) THEN

!                   The following information will be used to map the particle
                    smsgbuf(l)     = XD(1,i) + ximage(nb)
                    smsgbuf(l + 1) = XD(2,i) + yimage(nb)
                    smsgbuf(l + 2) = XD(3,i) + zimage(nb)

                    l = l + 3

10                  CONTINUE
                    smsgbuf(l)     = FDSW(1,i)
                    smsgbuf(l + 1) = FDSW(2,i)
                    smsgbuf(l + 2) = FDSW(3,i)
                    l = l + 3

!                   With STENSW passing the message is about twice longer,
!                   which slows down the run on at least 5% in total,
!                   so we will not pass STENSW until we do need it
                    IF(MOD(ISTEP,NEPRT).EQ.0.OR.ISTEP.EQ.1) THEN
                       smsgbuf(l)     = STENSW(i,1,1)
                       smsgbuf(l + 1) = STENSW(i,2,2)
                       smsgbuf(l + 2) = STENSW(i,3,3)
                       l = l + 3
                    ENDIF
     
                    i = link(i)
                    
                    IF(i.gt.0) GOTO 10
                    
                 ENDIF
                 
              ENDDO

           ENDIF

           msgsend = l - 1

!          Use blocking MPI to essure for a reliable passing of the information
!          Use index "nb" of the loop as a tag fo send/receive pair
           IF(mypid.NE.nbpid(nb)) THEN
              IF(MOD(nodecheck,2).EQ.0) THEN
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),52+nb,MPI_COMM_WORLD,info)
                 msgrecv = MSGSWP
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),52+nb,MPI_COMM_WORLD,status,info)
!                Need to count the lenght since it vvalue is not known
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ibelong = status(MPI_SOURCE)
              ELSE
                 msgrecv = MSGSWP
                 CALL MPI_RECV(rmsgbuf,msgrecv,MPI_DOUBLE_PRECISION, &
                      nbres(nb),52+nb,MPI_COMM_WORLD,status,info)
!                Need to count the lenght since it vvalue is not known 
                 CALL MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,icount,info)
                 msgrecv = icount
                 ibelong = status(MPI_SOURCE)
                 CALL MPI_SEND(smsgbuf,msgsend,MPI_DOUBLE_PRECISION, &
                      nbpid(nb),52+nb,MPI_COMM_WORLD,info)
              ENDIF
           ELSE
              msgrecv = msgsend
              rmsgbuf(1:msgrecv) = smsgbuf(1:msgsend)
              ibelong = mypid
           ENDIF

           l = 1

           IF(msgrecv.ne.0) THEN  

12            CONTINUE

!             Map the incoming particle to identify the related
!             link cell in order to distribute forces
              ixg = INT((rmsgbuf(l) - xlstart)*nxsign*nlcxg/XL + 1)
              iyg = INT((rmsgbuf(l + 1) - ylstart)*nysign*nlcyg/YL + 1)
              izg = INT((rmsgbuf(l + 2) - zlstart)*nzsign*nlczg/ZL + 1)
              
              ix = ixg - ixoffset
              iy = iyg - iyoffset
              iz = izg - izoffset

              ip = ix + nlcx2*(iy + nlcy2*iz)

              idx = ltop(ip)

              l = l + 3

11            CONTINUE

!             Having a prticle mapped once, we continue with the
!             current link cell util it is completed
              FDSW(1,idx) = FDSW(1,idx) + rmsgbuf(l)
              FDSW(2,idx) = FDSW(2,idx) + rmsgbuf(l + 1)
              FDSW(3,idx) = FDSW(3,idx) + rmsgbuf(l + 2)
              l = l + 3

              IF(MOD(ISTEP,NEPRT).EQ.0.OR.ISTEP.EQ.1) THEN
                 STENSW(idx,1,1) = STENSW(idx,1,1) + rmsgbuf(l)
                 STENSW(idx,2,2) = STENSW(idx,2,2) + rmsgbuf(l + 1)
                 STENSW(idx,3,3) = STENSW(idx,3,3) + rmsgbuf(l + 2)
                 l= l + 3
              ENDIF

              idx = link(idx)

              IF(idx.ne.0) GOTO 11

              IF((l-1).lt.msgrecv) GOTO 12

           ENDIF

        ENDDO neigh6

        DEALLOCATE(smsgbuf,rmsgbuf)

        RETURN
      END SUBROUTINE SHARE_SW
