      program gdr_real 

      implicit double precision (a-h,o-z)

      parameter(nax=512,nsx=3,ncorx=nsx*(nsx+1)/2)
      parameter(nofgx=600,mic=2,mic3=mic*mic*mic)

      integer alpha,beta

      common/metr/ ht(3,3),htm1(3,3)

      dimension tau(3,nax,nsx),stau(3,nax,nsx)
      dimension staur(3,mic3*nax,nsx),staug(3,mic3*nax,nsx)
      dimension sdist(3),rdist(3)

      dimension gdr(ncorx), aneigh(ncorx)
      dimension gdr_tmp(2*ncorx+1,nofgx)
      dimension sgdr_tmp(2*ncorx+1,nofgx)
      dimension ngdr(nofgx,ncorx), ngdrt(nofgx,ncorx)
      integer na(nsx)
      integer is1(ncorx), is2(ncorx)
      integer ierr

      character*25 fileato,filecel,dummy

      read(5,*,IOSTAT=ierr) fileato, filecel
      if(ierr.ne.0) then
        write(6,*) 'usage: gdr.x < input'
        write(6,*) '  input layout'
        write(6,*) '  file.pos file.cel'
        write(6,*) '  nsp nout'
        write(6,*) '  na(1) .. na(nsp)'
        write(6,*) '  is1(1)    is2(1)'
        write(6,*) '  ...'
        write(6,*) '  is1(nout) is2(nout)'
        write(6,*) '  npt  istart istride iostride'
        stop
      endif

      read(5,*) nsp, nout
      if(nsp .GT. nsx) then
        call error(' gdr ',' nsp > nsx ',nsx)
      end if
      if(nout .GT. ncorx) then
        call error(' gdr ',' nout > ncorx ',ncorx)
      end if
      read(5,*) (na(i),i=1,nsp)
      do i=1,nsp
        if(na(i) .GT. nax) then
          call error(' gdr ',' na(is) > nax ',i)
        end if
      end do
      do i=1,nout
        read(5,*) is1(i), is2(i)
        if(is1(i) .GT. nsp) then
          call error(' gdr ',' is1(i) > nsp ',i)
        end if
        if(is2(i) .GT. nsp) then
          call error(' gdr ',' is2(i) > nsp ',i)
        end if
      end do

      read(5,*) npt, istart, istride, iostride


      nat = 0
      do i = 1, nsp
        nat = nat + na(i)
      end do

      nofg   = 300
      pi     = 4.d0*datan(1.d0)
      box    = mic
      boxby2 = box/2.0d0
      res    = nofg/boxby2
     
      ICOUN=0
      IINF=0
      ISUP=0
      IF(MIC.EQ.1) THEN
        IINF=0
        ISUP=0
      ELSE IF(MIC.EQ.2) THEN
        IINF=0
        ISUP=1
      ELSE IF(MIC.EQ.3) THEN
        IINF=-1
        ISUP=1
      ELSE
        STOP 'MIC NOT PROGRAMMED'
      END IF

      DO iout = 1, nout
        GDR(iout)   =0.D0
        ANEIGH(iout)=0.D0
        DO ir = 1, nofgx
          NGDR(ir, iout) = 0
          NGDRT(ir,iout) = 0
        END DO
      END DO

C=======================================================================

      open(unit=1,file=fileato,status='old')
      open(unit=2,file=filecel,status='old')
      open(unit=3,file="gofr.out",status='unknown')

      omegat = 0.0d0
      nnpt = 0
      dmic = dble(mic)

      do ipt = 1, npt

        nfi = ipt * istride

        read(1,'(a)') dummy
        do is = 1, nsp
          do ia = 1, na(is)
            read(1,*) (tau(j,ia,is),j=1,3)
          end do
        end do

        read(2,'(a)') dummy
        do i=1,3
          read(2,*) (ht(i,j),j=1,3)
        end do

        if(ipt.ge.istart .and. mod(ipt,iostride).eq.0) then

          call inv3(ht,htm1,omega)
          omegat = omegat + omega
          nnpt = nnpt + 1

          box    = omega**(1.0d0/3.0d0) * mic
          boxby2 = box/2.0d0
          res    = nofg/boxby2
 
          do is = 1,nsp
            do ia=1,na(is)
               call r_to_s(tau(1,ia,is),stau(1,ia,is))
            end do  
          end do

C======================================================================= 
C==                             G(R)                                  ==
C==  ATT!  SCALED VARIABLES ARE USED UNTIL REAL DISTANCES ARE NEEDED  ==
C======================================================================= 

          DO ic = 1, nout
            DO ir = 1, nofg
              ngdrt(ir,ic) = 0
            END DO
          END DO

          DO IS=1,NSP
            DO IA=1,na(is)
              CALL PBCS(STAU(1,IA,is),STAU(2,IA,is),STAU(3,IA,is),
     C             STAUR(1,IA,is),STAUR(2,IA,is),STAUR(3,IA,is),1)
            END DO
          END DO

          DO IS=1,NSP
            IC = 0
            DO I=IINF,ISUP
              DO J=IINF,ISUP
                DO K=IINF,ISUP
                  DO IA=1,NA(IS)
                    IC = IC + 1
                    STAUG(1,IC,IS) = STAUR(1,IA,IS) + I
                    STAUG(2,IC,IS) = STAUR(2,IA,IS) + J
                    STAUG(3,IC,IS) = STAUR(3,IA,IS) + K
                  END DO
                END DO
              END DO
            END DO
          END DO

          DO iout = 1, nout

            K = is1(iout)
            J = is2(iout)

            IF(K .EQ. J) THEN
              LAX = NA(K) * MIC3 - 1
              IAD = 2
            ELSE 
              LAX = NA(K) * MIC3
              IAD = 1
            END IF

            DO L = 1, LAX

              IF(K .EQ. J) THEN
                INF = L + 1
              ELSE
                INF = 1
              END IF

              DO M = INF, NA(J) * MIC3

                XIJ = STAUG(1,L,K) - STAUG(1,M,J)
                YIJ = STAUG(2,L,K) - STAUG(2,M,J)
                ZIJ = STAUG(3,L,K) - STAUG(3,M,J)

                SDIST(1) = XIJ - DNINT(XIJ/DMIC)*DMIC
                SDIST(2) = YIJ - DNINT(YIJ/DMIC)*DMIC
                SDIST(3) = ZIJ - DNINT(ZIJ/DMIC)*DMIC

                CALL S_TO_R(SDIST,RDIST)
                ERRE2 = RDIST(1)**2 + RDIST(2)**2 + RDIST(3)**2

                NWH   = NINT(DSQRT(ERRE2)*RES)
                IF(NWH.GT.NOFG) NWH = NOFG

                NGDRT(NWH,iout) = NGDRT(NWH,iout) + IAD

              END DO
            END DO
          END DO

          DO ic = 1, nout
            DO ir = 1, nofg
              ngdr(ir,ic) = ngdr(ir,ic) + ngdrt(ir,ic)
            END DO
          END DO

        end if

      end do

      close(unit=1)
      close(unit=2)

C=======================================================================
C==                         PRINT G(R)                                ==
C=======================================================================

c     WRITE(6,1747)
      omega  = omegat / dble(nnpt)
      alat   = omega**(1.0d0/3.0d0)
      dens   = NAT / omega
      box    = alat * mic
      boxby2 = box / 2.0d0
      res    = nofg / boxby2

      DO iout = 1, nout
        DO ir = 1, nofg

          erre   = DBLE(ir) / res
          erre2  = erre * erre
          vshell = 4. * pi * erre2 / res
          naij   = na(is1(iout))
          IF (is1(iout).NE.is2(iout)) THEN
            dens   = DBLE(naij)  / omega
          ELSE
            dens   = DBLE(naij-1)/ omega
          END IF 

          gdr(iout)    = ngdr(ir,iout) / 
     c      DBLE(nnpt*mic3*naij*dens*vshell)
          aneigh(iout) = aneigh(iout) + ngdr(ir,iout) / 
     c      DBLE(nnpt*naij*mic3)
          gdr_tmp(1       ,ir) = erre * 0.52917
          gdr_tmp(2*iout  ,ir) = GDR(iout)
          gdr_tmp(2*iout+1,ir) = ANEIGH(iout)

        END DO
      END DO

      do ir = 1, nofg-1
         WRITE(3,1748) (gdr_tmp(iout,ir), iout = 1, 2*nout+1)
      end do

c      call smussa_gdr(gdr_tmp,sgdr_tmp,nv)
c      do i=1,nv
c         WRITE(4,1748) (SGDR_TMP(J,i),J=1,3)
c      end do


 1747 FORMAT(//4X,'R',8X,'G(R)',8X,'N(R)')
 1748 FORMAT(1X,13(F8.4,2X))

      stop 'ok!'
      end

C=======================================================================
C==            COMPUTES THE INVERSE OF A MATRIX 3*3                   ==
C=======================================================================

      SUBROUTINE INV3(A,C,DEN)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      REAL*8 A(3,3),B(3,3),C(3,3)

      DEN = 0.D0
      S   = 1.D0
      I   = 1
      J   = 2
      K   = 3

    1 DO IPERM=1,3
        DEN = DEN + S*A(1,I)*A(2,J)*A(3,K)
        L   = I
        I   = J
        J   = K
        K   = L
      END DO

      I = 2
      J = 1
      K = 3
      S = - S
      IF(S.LT.0.D0) GO TO 1

      IF(DABS(DEN).LT.1.D-20) CALL ERROR('INV3','SING. MATRIX: ',0)

      I = 1
      J = 2
      K = 3

      DO IR=1,3
        B(IR,1) = (A(2,J)*A(3,K) - A(2,K)*A(3,J)) / DEN
        B(IR,2) = (A(3,J)*A(1,K) - A(3,K)*A(1,J)) / DEN
        B(IR,3) = (A(1,J)*A(2,K) - A(1,K)*A(2,J)) / DEN
        L = I
        I = J
        J = K
        K = L
      END DO

      DO L=1,3
        DO K=1,3
          C(K,L) = B(K,L)
        END DO
      END DO

      RETURN
      END

C=======================================================================
C==  TRANSFORMS SCALED VECTOR S IN REAL VECTOR R (THROUGH MATRIX HT)  ==
C=======================================================================

      SUBROUTINE S_TO_R (S,R)

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION     S(*),R(*)

      COMMON/METR/ HT(3,3),HTM1(3,3)

c      DO I=1,3
c        R(I) = 0.D0
c        DO J=1,3
c          R(I) = R(I) + S(J)*HT(J,I)
c        END DO
c      END DO

      R(1) = S(1)*HT(1,1) + S(2)*HT(2,1) + S(3)*HT(3,1)
      R(2) = S(1)*HT(1,2) + S(2)*HT(2,2) + S(3)*HT(3,2)
      R(3) = S(1)*HT(1,3) + S(2)*HT(2,3) + S(3)*HT(3,3)

      RETURN
      END

C=======================================================================
C== TRANSFORMS REAL VECTOR R IN SCALED VECTOR S (THROUGH MATRIX HTM1) ==
C=======================================================================

      SUBROUTINE R_TO_S (R,S)

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION     S(*),R(*)

      COMMON/METR/ HT(3,3),HTM1(3,3)

      DO I=1,3
        S(I) = 0.D0

        DO J=1,3
          S(I) = S(I) + R(J)*HTM1(J,I)
        END DO

      END DO

      RETURN
      END

C=======================================================================
C==    COMPUTES THE PERIODIC BOUNDARY CONDITIONS IN THE SCALED        ==
C==                    VARIABLES SYSTEM                               ==
C=======================================================================

      SUBROUTINE PBCS(X1,Y1,Z1, X2,Y2,Z2, M)

      IMPLICIT REAL*8 (A-H,O-Z)

      REAL*8 MIC

      MIC = DBLE(M)

      X2 = X1 - DNINT(X1/MIC)*MIC
      Y2 = Y1 - DNINT(Y1/MIC)*MIC
      Z2 = Z1 - DNINT(Z1/MIC)*MIC
c      X2 = X1 - INT(X1/MIC)*MIC
c      Y2 = Y1 - INT(Y1/MIC)*MIC
c      Z2 = Z1 - INT(Z1/MIC)*MIC

      RETURN
      END

C===============================================================

      SUBROUTINE ERROR(A,B,N)

      CHARACTER*(*) A,B
      WRITE(6,1) A,B,N
    1 FORMAT(//' PROGRAM ',A,':',A,'.',8X,I8,8X,'STOP')
      STOP


      end


      SUBROUTINE MSD

!=======================================================================
!===    Calcola il MSD e il displacement per particella.
!===    Ha bisogno in input del file di traiettorie continue degli atomi  
!===    in coordinate reali.
!===    Output files: fort.2   disp. max per particella specie a 
!===                  fort.12  disp. max per particella specie b 
!===                  fort.3   disp. per particella specie a
!===                  fort.13  disp. per particella specie b
!===                  fort.4   MSD   specie a , b
!===                  
!===    Ultima modifica Carlo Cavazzoni. 
!=======================================================================

      IMPLICIT NONE

      real*8, parameter :: picosecond = 0.2418901D-4


      real*8, allocatable :: taui(:,:,:)
      real*8, allocatable :: tau(:,:,:)
      real*8, allocatable :: sqd(:,:)
      real*8, allocatable :: msqd(:)
      real*8              :: cdm(3),cdmi(3), r(3)
      real*8, allocatable :: mass(:)
      real*8              :: dt, masst, toffset, time
      integer, allocatable :: na(:)
      integer              :: nsp, nstep, nskip, nat, nax
      integer              :: is,ia,k,t,ierr
      character*80         :: atofile, dummy

!=======================================================================

      read(5,*,IOSTAT=ierr) atofile,nstep,nskip,dt,nsp
      if(ierr.ne.0) then
        write(6,*) 'usage: msd.x < input'
        write(6,*) '  input layout'
        write(6,*) '  file.pos nstep nskip timestep nsp'
        write(6,*) '  na(1) mass(1)'
        write(6,*) '  .'
        write(6,*) '  .'
        write(6,*) '  na(nsp) mass(nsp)'
        stop
      endif
 
      toffset = 0.0d0
      allocate(na(nsp))
      allocate(mass(nsp))

      do is = 1, nsp
        read(5,*) na(is),mass(is)
      end do

!=======================================================================

      nat   = sum(na)
      masst = sum(dble(na)*mass)
      nax   = maxval(na)
      allocate(tau(3,nax,nsp))
      allocate(taui(3,nax,nsp))
      allocate(sqd(nax,nsp))
      allocate(msqd(nsp))

      open(unit=34,file=atofile,status='old')

      tau = 0.0d0
      taui = 0.0d0

      do t = 1, nskip
        read(34,'(a)') dummy
        do is = 1, nsp
          do ia = 1, na(is)
            read(34,*) (tau(k,ia,is),k=1,3) 
          end do
        end do
      end do

      do t = nskip+1, nstep

        read(34,'(a)') dummy
        do is = 1, nsp
          do ia = 1, na(is)
            read(34,*) (tau(k,ia,is),k=1,3) 
          end do
        end do

        cdm = 0.d0
        do is = 1, nsp
          do ia = 1, na(is)
            cdm(:) = cdm(:) + tau(:,ia,is)*mass(is)
          end do
        end do
        cdm = cdm / masst

        if(t.eq.(nskip+1)) then
          cdmi = cdm
          do is = 1, nsp
            do ia = 1, na(is)
              taui(:,ia,is) = tau(:,ia,is) - cdmi(:)
            end do
          end do
        end if

        time = dble(t) * dt + toffset
        sqd  = 0.D0
        msqd = 0.D0
        do is = 1, nsp
          do ia = 1, na(is)
            r(:) = tau(:,ia,is) - cdm(:)
            sqd(ia,is)  = sum((r(:)-taui(:,ia,is))**2)
          end do
          msqd(is) = sum(sqd(:,is))/dble(na(is)) 
        end do
        write(4,100) time,time*picosecond,(msqd(is),is=1,nsp)
      end do
      close(unit=34)
      deallocate(na)
      deallocate(mass)
      deallocate(tau)
      deallocate(taui)
      deallocate(sqd)
      deallocate(msqd)

 100  format(f9.1,2x,f9.7,8f12.6)

      END SUBROUTINE
