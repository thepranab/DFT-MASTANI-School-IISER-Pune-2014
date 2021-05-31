!---------------------------------------------------------------------------
      program resp_mat
!---------------------------------------------------------------------------
!
! This program construct the response matrices to compute the Hubbard U.
! The responses dnI/dalphaJ are written in the files dn_(1:nat)_da_(1:ntyp) 
! for the scf response and dn0_(1:nat)_da_(1:ntyp) for the bare one.
!
      implicit none

      real(kind=8)   primv(3,3), atp_cry(3,200), atp_car(3,200),        &
                     d_cry(3), d_car(3), sum, sum0, atp_car_s(3,5000),   &
                     atp_cry_s(3,5000), daux
      real(kind=8), allocatable ::  chi(:,:), chi0(:,:), chim1(:,:),    &
     &                              chi0m1(:,:), alpha(:),              &
     &                              nn(:), dist(:,:), d0(:,:,:),        &
     &                              d0_s(:,:,:), dist_s(:,:),           &
     &                              chi_s(:,:), chi0_s(:,:),            &
     &                              chim1_s(:,:), chi0m1_s(:,:),        &
     &                              kp(:,:,:), kp0(:,:,:)

      integer  ntyp, na(10), nalfa, nat, nat0, nat_s, n1, n2, n3,       &
     &         na_s(10)
      integer  i, j, k, l, m, l1, l2, l3, i0
      integer, allocatable :: typ(:), typ_s(:), mol(:,:,:),             &
                              mol_s(:,:,:), i00(:,:), i000(:,:),        &
                              spin(:), spin_s(:), spinp(:,:,:),         &
                              spinp_s(:,:,:)
      integer ind, idx(10,10), done(10,10,10,2)

      character*20 filepos, back, filednda, filedndab, filednbdab
      character*20, allocatable :: filename(:,:), filename0(:,:)
      character*1 id1, id2, id3
      logical :: magn

      namelist /input_mat/ ntyp, na, nalfa, filepos, back, n1, n2, n3,  &
                           filednda, magn

      read(5,input_mat)

      nat = 0
      do i = 1,ntyp
         nat = nat + na(i)
      end do
      nat_s = nat*n1*n2*n3

      allocate(chi(nat+1,nat+1))
      allocate(chi0(nat+1,nat+1))
      allocate(chim1(nat+1,nat+1))
      allocate(chi0m1(nat+1,nat+1))
      allocate(alpha(nalfa))
      allocate(nn(nalfa))
      allocate(dist(nat,nat))
      allocate(typ(nat))
      allocate(mol(ntyp,ntyp,100))
      allocate(d0(ntyp,ntyp,100))
      allocate(typ_s(nat_s))
      allocate(mol_s(ntyp,ntyp,100))
      allocate(dist_s(nat_s,nat_s))
      allocate(d0_s(ntyp,ntyp,100))
      allocate(chi_s(nat_s+1,nat_s+1))
      allocate(chi0_s(nat_s+1,nat_s+1))
      allocate(chim1_s(nat_s+1,nat_s+1))
      allocate(chi0m1_s(nat_s+1,nat_s+1))
      allocate(kp(ntyp,ntyp,100))
      allocate(kp0(ntyp,ntyp,100))
      allocate(filename(ntyp,nat))
      allocate(filename0(ntyp,nat))
      allocate(i00(ntyp,ntyp))
      allocate(i000(ntyp,ntyp))
      allocate(spinp(ntyp,ntyp,100))
      allocate(spinp_s(ntyp,ntyp,100))
      allocate(spin(nat))
      allocate(spin_s(nat_s))

      open(11,file=filepos,status='unknown')
      
      spin(:) = 0
      spin_s(:) = 0
      do i = 1,3
         read(11,*) (primv(j,i),j=1,3)
      end do
      nat0 = 0
      if (magn) then
         do i = 1,ntyp
            do j = 1,na(i)
               nat0 = nat0 + 1
               read(11,*) (atp_cry(k,nat0),k=1,3), spin(nat0)
               typ(nat0) = i
            end do
         end do
      else
         do i = 1,ntyp
            do j = 1,na(i)
               nat0 = nat0 + 1
               read(11,*) (atp_cry(k,nat0),k=1,3)
               typ(nat0) = i
            end do
         end do
      end if
      close(11)

      open(11,file=filednda,status='unknown')

      do i = 1,ntyp
         do j = 1,nat
            read(11,*) filename(i,j), filename0(i,j)
         end do
      end do

      close(11)

      chi = 0.d0
      chi0 = 0.d0

! let's construct the first column of the response matrices

      nat0 = 0
      do i = 1,ntyp
         do j = 1,nat
            open(11,file=filename(i,j),status='unknown')
            do k = 1,nalfa
               read(11,*) alpha(k), nn(k)
            end do
            chi(j,nat0+1) = (nn(nalfa) - nn(1))/(alpha(nalfa)-alpha(1))
!            if (abs(nn(nalfa) - nn(1)).le.4.d-6) chi(j,nat0+1)=0.0d0
!            chi(j,nat0+1) = (nn(nalfa-1) - nn(2))/(alpha(nalfa-1)-alpha(2))
!            if (abs(nn(nalfa-1) - nn(2)).le.1.8d-4) chi(j,nat0+1)=0.0d0
!            chi(j,nat0+1) = (nn(nalfa) - nn(2))/(alpha(nalfa)-alpha(2))
!            if (abs(nn(nalfa-1) - nn(2)).le.3.d-6) chi(j,nat0+1)=0.0d0
            close(11)
            open(11,file=filename0(i,j),status='unknown')
            do k = 1,nalfa
               read(11,*) alpha(k), nn(k)
            end do
            chi0(j,nat0+1) = (nn(nalfa) - nn(1))/(alpha(nalfa)-alpha(1))
!            if (abs(nn(nalfa) - nn(1)).le.8.d-6) chi0(j,nat0+1)=0.0d0
!            chi0(j,nat0+1) = (nn(nalfa-1) - nn(2))/(alpha(nalfa-1)-alpha(2))
!            if (abs(nn(nalfa-1) - nn(2)).le.2.d-4) chi0(j,nat0+1)=0.0d0
!            chi0(j,nat0+1) = (nn(nalfa) - nn(2))/(alpha(nalfa)-alpha(2))
!            if (abs(nn(nalfa-1) - nn(2)).le.6.d-6) chi0(j,nat0+1)=0.0d0
            close(11)
         end do
         nat0 = nat0 + na(i)
      end do

! let's construct the other elements of the matrices; we need interatomic 
! distances

      atp_car = 0.d0
!      do i = 1,nat
!         do k = 1,3
!            if (atp_cry(k,i).ge.1.d0) atp_cry(k,i) = atp_cry(k,i) - 1.d0
!            if (atp_cry(k,i).lt.0.d0) atp_cry(k,i) = atp_cry(k,i) + 1.d0
!         end do
!      end do
      do i = 1,nat
         do j = 1,3
            do k = 1,3
               atp_car(j,i) = atp_car(j,i)+primv(j,k)*atp_cry(k,i)
            end do
         end do
      end do

      do i = 1,nat
         dist(i,i) = 0.d0
         do j = i+1,nat
            dist(i,j) = 1000.d0
            do l1 = -1,1
               do l2 = -1,1
                  do l3 = -1,1
                     daux = 0.d0
                     do k = 1,3
                        daux = daux + (atp_car(k,i) - atp_car(k,j) -    &
     &                         dfloat(l1)*primv(k,1) -                  &
     &                         dfloat(l2)*primv(k,2) -                  &
     &                         dfloat(l3)*primv(k,3) )**2
                     end do
                     daux = dsqrt(daux)
                     dist(i,j) = min(dist(i,j),daux)
                  end do
               end do
            end do
            dist(j,i) = dist(i,j)
         end do
      end do

      do j = 1,nat
         do l = 1,nat
            if (chi(j,l).eq.0.d0) then
               nat0 = 0
!               do k = 1,typ(l)-1
!                  nat0 = nat0+na(k)
!               end do
!               nat0 = nat0 + 1
               do k = 1,nat
                  if (typ(k).eq.typ(l)) then
                     nat0 = k
                     go to 35
                  end if
               end do
 35            continue
               do k = 1,nat
                  if (abs(dist(k,nat0)-dist(j,l)).le.1.d-5 .and.        &
     &                typ(k).eq.typ(j).and.typ(nat0).eq.typ(l) .and.    &
     &                spin(j)*spin(l).eq.spin(nat0)*spin(k)) then
                     chi(j,l)=chi(k,nat0)
                     chi0(j,l)=chi0(k,nat0)
                  end if
               end do
            end if
         end do
      end do

! symmetrization

      do i = 1,nat
         do j = 1,nat
            chi(i,j) = 0.5d0*(chi(i,j)+chi(j,i))
            chi(j,i) = chi(i,j)
            chi0(i,j) = 0.5d0*(chi0(i,j)+chi0(j,i))
            chi0(j,i) = chi0(i,j)
         end do
      end do

! molteplicities

      mol = 0
      do i = 1,ntyp
         do l = 1,nat
            if (typ(l).eq.i) then 
               i0 = l
               go to 9
            end if
         end do
 9       continue
         do j = i,ntyp
            i00(i,j) = 0
            do nat0 = 1,nat
               if (typ(nat0).eq.j) then
                  if (i00(i,j).eq.0) then
                     d0(i,j,1) = dist(i0,nat0)
                     mol(i,j,1) = 1
                     i00(i,j) = 1
                     kp(i,j,1) = chi(nat0,i0)
                     kp0(i,j,1) = chi0(nat0,i0)
                     spinp(i,j,1) = spin(nat0)*spin(i0)
                  else
                     do k = i00(i,j),1,-1
                        if (abs(dist(i0,nat0)-d0(i,j,k)).le.1.d-5 .and. &
     &                      spin(i0)*spin(nat0).eq.spinp(i,j,k)) then
                           mol(i,j,k) = mol(i,j,k) + 1
                           go to 10
                        end if
                     end do
                     i00(i,j) = i00(i,j) + 1
                     d0(i,j,i00(i,j)) = dist(i0,nat0)
                     mol(i,j,i00(i,j)) = 1
                     kp(i,j,i00(i,j)) = chi(nat0,i0)
                     kp0(i,j,i00(i,j)) = chi0(nat0,i0)
                     spinp(i,j,i00(i,j)) = spin(nat0)*spin(i0)
 10                  continue
                  end if
               end if
            end do
            i00(j,i) = i00(i,j)
            d0(j,i,:) = d0(i,j,:)
            mol(j,i,:) = mol(i,j,:)
            spinp(j,i,:) = spinp(i,j,:)
         end do
      end do
               
! background term
     
      if (back.eq.'neutral') then
         sum = 0.d0
         sum0 = 0.d0
         do i = 1,nat
            chi(i,nat+1) = 0.d0
            chi0(i,nat+1) = 0.d0
            do j = 1,nat
               chi(i,nat+1) = chi(i,nat+1) - chi(i,j)
               chi0(i,nat+1) = chi0(i,nat+1) - chi0(i,j)
            end do
            chi(nat+1,i) = chi(i,nat+1)
            chi0(nat+1,i) = chi0(i,nat+1)
            sum = sum - chi(i,nat+1)
            sum0 = sum0 - chi0(i,nat+1)
         end do
         chi(nat+1,nat+1) = sum
         chi0(nat+1,nat+1) = sum0
         do i = 1,nat+1
            do j = 1,nat+1
               chi(i,j) = chi(i,j) + 0.01
               chi0(i,j) = chi0(i,j) + 0.01
            end do
         end do
      else if (back.eq.'computed') then
         open(11,file=filedndab,status='unknown')
         do j = 1,nat
            read(11,*) filename(1,j), filename0(1,j)
         end do
         close(11)
         do i = 1,1
            do j = 1,nat
               open(11,file=filename(i,j),status='unknown')
               do k = 1,nalfa
                  read(11,*) alpha(k), nn(k)
               end do
               chi(j,nat+1) = (nn(nalfa) - nn(1))/(alpha(nalfa)-alpha(1))
               chi(nat+1,j) = chi(j,nat+1)
               close(11)
               open(11,file=filename0(i,j),status='unknown')
               do k = 1,nalfa
                  read(11,*) alpha(k), nn(k)
               end do
               chi0(j,nat+1) = (nn(nalfa) - nn(1))/(alpha(nalfa)-alpha(1))
               chi0(nat+1,j) = chi0(j,nat+1)
               close(11)
            end do
         end do
         open(11,file=filednbdab,status='unknown')
         do j = 1,nat
            read(11,*) filename(1,j), filename0(1,j)
         end do
         close(11)
         do i = 1,1
            do j = 1,nat
               open(11,file=filename(i,j),status='unknown')
               do k = 1,nalfa
                  read(11,*) alpha(k), nn(k)
               end do
               chi(nat+1,nat+1) = chi(nat+1,nat+1) +                    &
     &                       (nn(nalfa) - nn(1))/(alpha(nalfa)-alpha(1))
               close(11)
               open(11,file=filename0(i,j),status='unknown')
               do k = 1,nalfa
                  read(11,*) alpha(k), nn(k)
               end do
               chi0(nat+1,nat+1) = chi0(nat+1,nat+1) +                  &
     &                       (nn(nalfa) - nn(1))/(alpha(nalfa)-alpha(1))
               close(11)
            end do
         end do
      else if (back.eq.'no'.or.back.eq.' ') then
         chi(nat+1,nat+1) = 1.d0
         chi0(nat+1,nat+1) = 1.d0
      end if 

      call invmat(chi,chim1,nat+1)
      call invmat(chi0,chi0m1,nat+1)

      if (back.eq.'neutral') then
         do i = 1,nat+1
            do j = 1,nat+1
               daux = 0.d0
               do l = 1,nat+1
                  daux = daux + chi0(i,l)*chi0m1(l,j)
               end do
               chi(i,j) = chi(i,j) - 0.01
               chi0(i,j) = chi0(i,j) - 0.01
            end do
         end do
      end if

! Let's write the result for U in the primitive cell.

      open(12,file='Umat.out',status='unknown')

      write(12,*)
      write(12,*) ' number of atoms in the primitive cell:', nat
      write(12,*) 
      if (nat.le.32) then
         write(12,*)
         write(12,*) ' CHI_0 Matrix '
         write(12,*)
         do i = 1,nat+1
            write(12,'(8(x,f8.4))') (chi0(i,j),j=1,nat+1)
         end do
         write(12,*)
         
         write(12,*)
         write(12,*) ' CHI Matrix '
         write(12,*)
         do i = 1,nat+1
            write(12,'(8(x,f8.4))') (chi(i,j),j=1,nat+1)
         end do
         write(12,*)
   
         write(12,*)
         write(12,*) ' CHI0^-1 - CHI^-1 Matrix '
         write(12,*)
         do i = 1,nat+1
            write(12,'(8(x,f9.4))') (chi0m1(i,j)-chim1(i,j),j=1,nat+1)
         end do
      else

         write(12,*)
         write(12,*) ' CHI0 Matrix '
         write(12,*)
         do i = 1,ntyp
            do nat0 = 1, nat
               if (typ(nat0).eq.i) then
                  do j = 1,nat
                     write(12,'(i3,x,i3,x,i5,x,i5,2(x,f12.6))') i, typ(j), &
     &                  nat0, j, dist(nat0,j), chi0(nat0,j)
                  end do
                  go to 19
               end if
            end do
 19         continue
         end do

         write(12,*)
         write(12,*) ' CHI Matrix '
         write(12,*)
         do i = 1,ntyp
            do nat0 = 1, nat
               if (typ(nat0).eq.i) then
                  do j = 1,nat
                     write(12,'(i3,x,i3,x,i5,x,i5,2(x,f12.6))') i, typ(j), &
     &                  nat0, j, dist(nat0,j), chi(nat0,j)
                  end do
                  go to 29
               end if
            end do
 29         continue
         end do

         write(12,*)
         write(12,*) ' CHI0^-1 - CHI^-1 Matrix '
         write(12,*)
         do i = 1,ntyp
            do nat0 = 1, nat
               if (typ(nat0).eq.i) then
                  do j = 1,nat
                     write(12,'(4(i3,x),i5,x,i5,2(x,f12.6))') i, typ(j), &
     &                 spin(nat0), spin(j), nat0, j,                &
     &                 dist(nat0,j), chi0m1(nat0,j)-chim1(nat0,j)
                  end do
                  go to 39
               end if
            end do
 39         continue
         end do
      end if

      nat0 = 0
      do i = 1,ntyp
         do j = 1,nat
            if (typ(j).eq.i) then
               write(12,*) ' type: ', i, ' U0 = ', chi0m1(j,j)-          &
     &                     chim1(j,j)
               go to 15
            end if
         end do
 15      continue
      end do

      write(12,*)
      write(12,*)
      write(12,*)
      write(12,*) ' Extrapolation to a supercell ', n1, n2, n3
      
! extrapolation to supercells

      na_s = 0
      nat_s = 0
      do i = 1,n1
         do j = 1,n2
            do k = 1,n3
               do l = 1,nat
                  nat_s = nat_s + 1
                  atp_cry_s(1,nat_s) = atp_cry(1,l) + dfloat(i-1)
                  atp_cry_s(2,nat_s) = atp_cry(2,l) + dfloat(j-1)
                  atp_cry_s(3,nat_s) = atp_cry(3,l) + dfloat(k-1)
                  typ_s(nat_s) = typ(l)
                  spin_s(nat_s) = spin(l)
                  na_s(typ(l)) = na_s(typ(l)) + 1
               end do
            end do
         end do
      end do

      atp_car_s = 0.d0
      do i = 1,nat_s
         do j = 1,3
            do k = 1,3
               atp_car_s(j,i) = atp_car_s(j,i)+primv(j,k)*atp_cry_s(k,i)
            end do
         end do
      end do

      do i = 1,nat_s
         dist_s(i,i) = 0.d0
         do j = i+1,nat_s
            dist_s(i,j) = 1000.d0
            do l1 = -n1,n1,n1
               do l2 = -n2,n2,n2
                  do l3 = -n3,n3,n3
                     daux = 0.d0
                     do k =1,3
                        daux = daux + ( atp_car_s(k,i) - atp_car_s(k,j) &
     &                                - dfloat(l1)*primv(k,1)           &
     &                                - dfloat(l2)*primv(k,2)           &
     &                                - dfloat(l3)*primv(k,3) )**2
                     end do
                     daux = dsqrt(daux)
                     dist_s(i,j) = min(dist_s(i,j),daux)
                  end do
               end do
            end do
            dist_s(j,i) = dist_s(i,j)
         end do
      end do

! molteplicities

      write(12,*) 
      write(12,*) ' number of atoms in the supercell:', nat_s
      write(12,*) 

      mol_s = 0
      do i = 1,ntyp
         do l = 1,nat_s
            if (typ_s(l).eq.i) then 
               i0 = l
               go to 11
            end if
         end do
 11      continue
         do j = i,ntyp
            i000(i,j) = 0
            do nat0 = 1,nat_s
               if (typ_s(nat0).eq.j) then
                  if (i000(i,j).eq.0) then
                     d0_s(i,j,1) = dist_s(i0,nat0)
                     mol_s(i,j,1) = 1
                     i000(i,j) = 1
                     spinp_s(i,j,1) = spin_s(nat0)*spin_s(i0)
                  else
                     do k = i000(i,j),1,-1
                        if (abs(dist_s(i0,nat0)-d0_s(i,j,k)).le.1.d-5   &
     &                      .and. spin_s(i0)*spin_s(nat0).eq.spinp_s(i,j,k)) then
                           mol_s(i,j,k) = mol_s(i,j,k) + 1
                           go to 12
                        end if
                     end do
                     i000(i,j) = i000(i,j) + 1
                     d0_s(i,j,i000(i,j)) = dist_s(i0,nat0)
                     mol_s(i,j,i000(i,j)) = 1
                     spinp_s(i,j,i000(i,j)) = spin_s(nat0)*spin_s(i0)
 12                  continue
                  end if
               end if
            end do
            i000(j,i) = i000(i,j)
            d0_s(j,i,:) = d0_s(i,j,:)
            mol_s(j,i,:) = mol_s(i,j,:)
            spinp_s(j,i,:) = spinp_s(i,j,:)
            if (j.ge.i) then
               write(12,*)
               write(12,*) ' type1: ', i, ' type2: ', j
               write(12,*)
               write(12,*) ' number of shells: ', i000(i,j)
               write(12,*)
               write(12,*) ' shell  #n: pc, sc, spinp   dist       chi_pc     chi_sc '
               write(12,*)
               do k = 1,i000(i,j)
                  do m = 1,i00(i,j)
                     if (abs(d0_s(i,j,k)-d0(i,j,m)).le.1.d-5.and.       &
     &                   spinp(i,j,m).eq.spinp_s(i,j,k)) then
                        write(12,'(2x,i3,5x,3(2x,i2),3x,3(3x,f8.4))')   &
     &                         k, mol(i,j,m),mol_s(i,j,k), spinp(i,j,m),    &
     &                         d0_s(i,j,k), &
     &                         kp0(i,j,m), kp0(i,j,m)*dfloat(mol(i,j,m))&
     &                         /dfloat(mol_s(i,j,k))
                        go to 89
                     end if
                  end do
                  write(12,'(2x,i3,5x,3(2x,i2),3x,3(3x,f8.4))')         &
     &                  k, mol(i,j,m+1),mol_s(i,j,k), spinp_s(i,j,k),       &
     &                  d0_s(i,j,k),      &
     &                  kp0(i,j,m+1), kp0(i,j,m+1)*dfloat(mol(i,j,m+1)) &
     &                  /dfloat(mol_s(i,j,k))
 89               continue
               end do
            end if
         end do
      end do
               
               
      do i = 1,ntyp
         do j = i,ntyp
            do k = 1,i00(i,j)
                mol(j,i,k) = mol(i,j,k)
                d0(j,i,k) = d0(i,j,k)
                kp(j,i,k) = kp(i,j,k)
                kp0(j,i,k) = kp0(i,j,k)
                i00(j,i) = i00(i,j)
            end do
            do k = 1,i000(i,j)
                mol_s(j,i,k) = mol_s(i,j,k)
                d0_s(j,i,k) = d0_s(i,j,k)
                i000(j,i) = i000(i,j)
            end do
         end do
      end do

      do i = 1,nat_s
         do k = i, nat_s
            do l = 1,i00(typ_s(i),typ_s(k))
               do j = 1,i000(typ_s(i),typ_s(k))
                  if (abs(d0_s(typ_s(i),typ_s(k),j)-                    &
     &                 d0(typ_s(i),typ_s(k),l)).le.1.d-5 .and.          &
     &                   spinp(typ_s(i),typ_s(k),l).eq.spinp_s(typ_s(i),typ_s(k),j)) then
                     m = j
                     go to 333
                  end if
               end do
 333           continue
               if (abs(dist_s(i,k)-d0(typ_s(i),typ_s(k),l)).le.1.d-5     &
     &             .and. spin_s(i)*spin_s(k).eq.spinp(typ_s(i),typ_s(k),l)) then
                  chi_s(i,k) = kp(typ_s(i),typ_s(k),l)*                  &
                                  dfloat(mol(typ_s(i),typ_s(k),l))/         &
                                  dfloat(mol_s(typ_s(i),typ_s(k),m))
                  chi0_s(i,k) = kp0(typ_s(i),typ_s(k),l)*                &
                                   dfloat(mol(typ_s(i),typ_s(k),l))/        &
                                   dfloat(mol_s(typ_s(i),typ_s(k),m))
                  go to 14
               end if
            end do
            chi_s(i,k) = 0.d0
            chi0_s(i,k) = 0.d0
 14         continue
            chi_s(k,i) = chi_s(i,k)
            chi0_s(k,i) = chi0_s(i,k)
         end do
      end do

! background term
     
      if (back.eq.'neutral') then
         sum = 0.d0
         sum0 = 0.d0
         do i = 1,nat_s
            chi_s(i,nat_s+1) = 0.d0
            chi0_s(i,nat_s+1) = 0.d0
            do j = 1,nat_s
               chi_s(i,nat_s+1) = chi_s(i,nat_s+1) - chi_s(i,j)
               chi0_s(i,nat_s+1) = chi0_s(i,nat_s+1) - chi0_s(i,j)
            end do
            chi_s(nat_s+1,i) = chi_s(i,nat_s+1)
            chi0_s(nat_s+1,i) = chi0_s(i,nat_s+1)
            sum = sum - chi_s(i,nat_s+1)
            sum0 = sum0 - chi0_s(i,nat_s+1)
         end do
         chi_s(nat_s+1,nat_s+1) = sum
         chi0_s(nat_s+1,nat_s+1) = sum0
         do i = 1,nat_s+1
            do j = 1,nat_s+1
               chi_s(i,j) = chi_s(i,j) + 0.01
               chi0_s(i,j) = chi0_s(i,j) + 0.01
            end do
         end do
      else if (back.eq.'computed') then
         do i = 1,nat_s
            do k = 1,ntyp
               if (typ_s(i).eq.k) then
                  nat0 = 0
                  do l = 1,k-1
                     nat0 = nat0+na(l)
                  end do
                  nat0 = nat0+1
               end if
            end do
            chi_s(nat_s+1,i) = chi(nat+1,nat0+1)
            chi_s(i,nat_s+1) = chi_s(nat_s+1,i)
            chi0_s(nat_s+1,i) = chi0(nat+1,nat0+1)
            chi0_s(i,nat_s+1) = chi0_s(nat_s+1,i)
         end do
         chi_s(nat_s+1,nat_s+1) = chi(nat+1,nat+1)
         chi0_s(nat_s+1,nat_s+1) = chi0(nat+1,nat+1)
      else if (back.eq.'no'.or.back.eq.' ') then
         chi_s(nat_s+1,nat_s+1) = 1.d0
         chi0_s(nat_s+1,nat_s+1) = 1.d0
      end if 
     
      call invmat(chi_s,chim1_s,nat_s+1)
      call invmat(chi0_s,chi0m1_s,nat_s+1)

      if (back.eq.'neutral') then
         do i = 1,nat_s+1
            do j = 1,nat_s+1
               daux = 0.d0
               do l = 1,nat_s+1
                  daux = daux + chi0_s(i,l)*chi0m1_s(l,j)
               end do
               chi_s(i,j) = chi_s(i,j) - 0.01
               chi0_s(i,j) = chi0_s(i,j) - 0.01
            end do
         end do
      end if

! Let's write the result for U in the supercell.

      if (nat_s.le.32) then
         write(12,*)
         write(12,*) ' CHI_0 Matrix '
         write(12,*)
         do i = 1,nat_s+1
            write(12,'(8(x,f8.4))') (chi0_s(i,j),j=1,nat_s+1)
         end do
         write(12,*)
         
         write(12,*)
         write(12,*) ' CHI Matrix '
         write(12,*)
         do i = 1,nat_s+1
            write(12,'(8(x,f8.4))') (chi_s(i,j),j=1,nat_s+1)
         end do
         write(12,*)
   
         write(12,*)
         write(12,*) ' CHI0^-1 - CHI^-1 Matrix '
         write(12,*)
         do i = 1,nat_s+1
!            write(12,*) (chi0m1(i,j),j=1,nat+1)
            write(12,'(8(x,f9.4))') (chi0m1_s(i,j)-chim1_s(i,j),j=1,nat_s+1)
         end do
      else
         write(12,*)
         write(12,*) ' CHI0 Matrix '
         write(12,*)
         do i = 1,ntyp
            do nat0 = 1, nat_s
               if (typ_s(nat0).eq.i) then
                  do j = 1,nat_s
                     write(12,'(4(i3,x),i5,x,i5,2(x,f12.6))') i, typ_s(j), &
     &                 spin_s(nat0), spin_s(j), nat0, j, dist_s(nat0,j), chi0_s(nat0,j)
                  end do
                  go to 13
               end if
            end do
 13         continue
         end do

         write(12,*)
         write(12,*) ' CHI Matrix '
         write(12,*)
         do i = 1,ntyp
            do nat0 = 1, nat_s
               if (typ_s(nat0).eq.i) then
                  do j = 1,nat_s
                     write(12,'(4(i3,x),i5,x,i5,2(x,f12.6))') i, typ_s(j), &
     &                 spin_s(nat0), spin_s(j), nat0, j, dist_s(nat0,j), chi_s(nat0,j)
                  end do
                  go to 23
               end if
            end do
 23         continue
         end do

         write(12,*)
         write(12,*) ' CHI0^-1 - CHI^-1 Matrix '
         write(12,*)
         do i = 1,ntyp
            do nat0 = 1, nat_s
               if (typ_s(nat0).eq.i) then
                  do j = 1,nat_s
                     write(12,'(4(i3,x),i5,x,i5,2(x,f12.6))') i, typ_s(j), &
     &                 spin_s(nat0), spin_s(j), nat0, j,                &
     &                 dist_s(nat0,j), chi0m1_s(nat0,j)-chim1_s(nat0,j)
                  end do
                  go to 33
               end if
            end do
 33         continue
         end do
      end if   

      do i = 1,ntyp
         do j = 1,nat_s
            if (typ_s(j).eq.i) then
               write(12,*) ' type: ', i, ' U1 = ', chi0m1_s(j,j)-          &
     &                     chim1_s(j,j)
               go to 355
            end if
         end do
 355     continue
      end do

      close(12)

!      open(11,file='Chimat',status='unknown')
!      do i = 1,nat_s
!         do j = i,nat_s
!            write(11,'(2(x,i6),3(x,f12.8))') i, j, chi0_s(i,j), &
!                 chi_s(i,j)!, chi0m1_s(i,j)-chim1_s(i,j)
!         end do
!      end do  
!      close(11)

      deallocate(chi, chi0, chim1, chi0m1, alpha, nn, dist, typ,        &
     &           mol, d0, typ_s, mol_s, dist_s, d0_s, chi_s, chi0_s,    &
     &           chim1_s, chi0m1_s, kp, kp0, filename, filename0,       &
     &           i00, i000, spin, spin_s )

      stop
 
      end 

!----------------------------------------------------------------------
      subroutine invmat (a, a_inv, n)
!----------------------------------------------------------------------
! computes the inverse "a_inv" of matrix "a", both dimensioned (n,n)
! matrix "a" is unchanged on output - LAPACK
!
!  use parameters
      implicit none
      integer :: n, i, j
      real(kind=8) :: a(n, n), a_inv(n, n)
                                                                                
      integer, parameter :: lwork = 500000
! lwork is the dimension of work and ipiv. Must be >= n
      integer :: info, lda, ipiv (lwork)
! info=0: inversion was successful
! lda   : leading dimension (the same as n)
! ipiv  : work space for pivoting
      real(kind=8) :: work (lwork)
                                                                                
      if (n.gt.lwork) write(6,*) 'work space is too small ', n
      lda = n
      call DCOPY (n * n, a, 1, a_inv, 1)
      call DGETRF (n, n, a_inv, lda, ipiv, info)
      call DGETRI (n, a_inv, lda, ipiv, work, lwork, info)

      do i = 1,n
         do j = 1,n
!            if (i.le.4.and.j.le.4) write(6,*) 'cazzo: ', i, j, a(i,j)
         end do
      end do
                                                                                
      return
                                                                                
      end subroutine invmat
