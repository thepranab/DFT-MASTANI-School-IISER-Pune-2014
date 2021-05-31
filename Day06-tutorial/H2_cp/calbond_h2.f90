PROGRAM calculate

IMPLICIT NONE
INTEGER :: ns, nax, nat, is, ia, i, j, k, nsteps, nbonds, nangs, ib, iang, nstart
REAL(KIND=8) :: rmag, rang, r(3), r1(3), r2(3), dt, rmag2
REAL(KIND=8), ALLOCATABLE :: pos(:,:,:)
INTEGER, ALLOCATABLE :: na(:) !, bndx_is(:,:), bndx_ia(:,:), angx_is(:,:), ang_ia(:,:)
CHARACTER(LEN=50) :: fname

READ(5,*) fname
READ(5,*) ns
ALLOCATE(na(ns))
READ(5,*) na(:)
nax = MAXVAL(na)
WRITE(6,*) nax
READ(5,*) nsteps, nstart, dt

!READ(5,*) nbonds
!ALLOCATE(bndx_ia(nbonds,2), bndx_is(nbonds,2))
!DO i=1,nbonds
! READ(5,*) bndx_is(i,1), bndx_ia(i,1), bndx_is(i,2), bndx_ia(i,2)
!END DO

!READ(5,*) nangs
!ALLOCATE(angx_ia(nangs,3), angx_is(nangs,3))
!DO i=1,nbonds
! READ(5,*) angx_is(i,1), angx_ia(i,1), angx_is(i,2), angx_ia(i,2), angx_is(i,3), angx_ia(i,3)
!END DO

ALLOCATE(pos(ns,nax,3))

j = 0


OPEN(UNIT=1,FILE=ADJUSTL(fname))
OPEN(UNIT=2,FILE="bonds.out")
DO i=1, nsteps
 
 READ(1,*) 
 DO is = 1,ns
  DO ia=1,na(is)
   READ(1,*) pos(is,ia,:)
  END DO
 END DO

 IF (i < nstart) CYCLE

 j = j + 1

!**** Specific to H2O molecule
 WRITE(2,"(I8,2X,F10.6,2X)", ADVANCE="NO") j, (j-1)*dt
 r1 = pos(1,2,:)-pos(1,1,:)

 rmag = DSQRT(DOT_PRODUCT(r1,r1))
 WRITE(2,"(F10.6,2X)") rmag

END DO
CLOSE(UNIT=2)
CLOSE(UNIT=1)

END PROGRAM calculate
