
program molcasto47

  ! convert calculated data from an OpenMolcas run and the resulting
  ! H5 file into a generaic 'FILE47' input for the NBO program. A
  ! separate Molcas-formatted orbital data file (ascii) can be
  ! specified, in which case the orbitals and density matrix in FILE47
  ! will be generated from the orbital coefficients and the
  ! occupations in the orbital file.

  ! (c) 2024, Aleksandr Zaichenko and Jochen Autschbach

  use definitions
  use hdf5
  use h5extractor
  use orbitals
!  use functions

  implicit none

  ! ============================================================================

  real(KREAL), dimension(:), allocatable :: overlap_out, charges, fock_out, &
    mo_out, occ_out, contraction

  real(KREAL), dimension(:), allocatable :: prim, dens_out, coord, DS, &
    EXPONENTS, CS, CP, CD, CF, CG, CH, CI

  integer(HID_T), dimension(:), allocatable :: basis, atoms

  integer, dimension(:), allocatable :: IDSP, IDSB, label, n0nums

  character*(LCHARS) :: h5file, orbfile, f_flag

  real(KREAL), dimension(:,:), allocatable :: C, n, F, S, A, P, D

  integer(KINT) :: Nbas, np, ns, spn, Nshell, lbas, mbas

  integer(KINT) :: i, j, k, m, l, Lmax, irrep, ios

  integer, dimension(:,:), allocatable :: basis_ident
  
  real(KREAL) :: rtemp

  logical :: SYM, have_orb, save_fock

  integer, dimension(:), allocatable :: nprim, npntr, ncomps

  real(KREAL), parameter :: zero = 0.0d0
 
  ! ============================================================================

  ! check command line options, check if input files can be opened:

  if (iargc().lt.1) stop 'need at least one file name argument. aborting'
  if (iargc().gt.3) stop 'looks like there are too many arguments. aborting'

  call getarg(1,h5file)
  ios = 0
  open(iuh, file=trim(h5file), iostat=ios)
  if (ios /= 0) then
    stop 'file '//trim(h5file)//' cannot be opened. aborting'
  else
    write(out,'(/1x,a)') 'will use Molcas H5 file '//trim(h5file)
    close(iuh)
  endif

  save_fock = .false.
  have_orb = .false. ! flag for whether to use orbital file data
  if (iargc().ge.2) then
    call getarg(2,orbfile)
    if (trim(orbfile) == '-fock') then
       save_fock = .true.
       write(out,*) 'Fock matrix will be written!'
    else
    ios = 0
    open(iuo, file=trim(orbfile), iostat=ios)
    if (ios /= 0) then
      orbfile = 'filenotfound'
      write(out,*) 'orbitals file cannot be opened. Will proceed suspiciously'
      have_orb = .false.
    else
      read(iuo, *, iostat=ios)
      if (ios /= 0) then
        stop 'cannot read from orbitals file. aborting'
      end if
      ! note: orbital file unit iuo remains open
      have_orb = .true.
    endif ! ios
    endif
  end if ! iargc()

  if (iargc().eq.3) then
    call getarg(3,f_flag)
    if (trim(f_flag) == '-fock') then
       save_fock = .true.
       write(out,*) 'Fock matrix will be written!'
    end if
  endif 

  ! Extraction and ordering of the data from h5 file

  call h5molcas(h5file, atoms, charges, coord, overlap_out, fock_out,&
    & mo_out, occ_out, basis, prim, IDSP, IDSB, dens_out, DS, SYM)

  Nbas = int(sum(basis),kind(Nbas))

  allocate(C(Nbas,Nbas))
  allocate(n(Nbas,Nbas))
  allocate(A(Nbas,Nbas))
  allocate(P(Nbas,Nbas))
  allocate(S(Nbas,Nbas))
  allocate(F(Nbas,Nbas))
  allocate(n0nums(size(prim)/2))
  allocate(label(Nbas))

  ! Initilize variables 
  C = 0 
  n = 0 
  A = 0 
  S = 0
  F = 0

  l = 0 
  m = 0
  if (SYM) then
    k = 0
    do irrep = 1, size(basis)
      do i = 1,int(basis(irrep),kind(i))
        if (i<0) stop ('Error: possible integer overflow in i-loop!')
        do j = 1,int(basis(irrep),kind(j))
          if (j<0) stop ('Error: possible integer overflow in j-loop!')
          k = k + 1
          C(m+j,m+i) = mo_out(k)
          S(m+j,m+i) = overlap_out(k)
          F(m+j,m+i) = fock_out(k)
          if (i == j) then
            l = l+1
            n(m+i,m+j) = occ_out(l)
          endif
        enddo
      enddo
      m = m + int(basis(irrep),kind(m))
    enddo
  else
    k = 0
    do i = 1,Nbas
      do j = 1,Nbas
        k = k + 1
        C(j,i) = mo_out(k)
        S(j,i) = overlap_out(k)
        F(j,i) = fock_out(k)
        if (i == j) then
          n(i,j) = occ_out(i)
        endif
      enddo
    enddo
  endif

  ! If an orbital file was provided, then use the data.
  ! Arrays C and n will be overwritten. Afterwards, we
  ! close the orbital data file.

  if (have_orb) then
    write(out,*) 'will use orbital file ', trim(orbfile)
    call getorbitals(nbas,C,n,SYM,basis)
    close (iuo)
  end if

 ! write(*,*) 'MO:'
 ! write(*,'(10F8.4)') C
 ! write(*,*) 'Occ:'
 !  write(*,'(10F8.4)') (n(j,j), j = 1,Nbas)


  ! Aleksandr, please comment what is done next. Looks like a density
  ! matrix is calculated. Also, we appear to be placing the
  ! desymmetrization matrix DS into array D but all you would need to
  ! do is an array reshape to go from 1D arrays to 2D arrays instead
  ! of duplicating the arrays. We can do that later, however, after
  ! everything is confirmed to work as intended.

  ! Yes. Here and in some other cases can be used reshape or pointers
  ! but I am not very familiar with the reshape function. I will try 
  ! to do it later

  if (SYM) then
    !write the desymmetrization matrix
    
    allocate(D(Nbas,Nbas))
    k = 0
    D = 0
    do i = 1,Nbas
      do j = 1,Nbas
        k = k + 1
        D(j,i) = DS(k)
      enddo
    enddo
   
   !calculate density matrix  
    call DGEMM('N','T',Nbas,Nbas,Nbas,1.0D0,n,Nbas,C,Nbas,0.0D0,A,Nbas)
    call DGEMM('N','N',Nbas,Nbas,Nbas,1.0D0,C,Nbas,A,Nbas,0.0D0,P,Nbas)
   !dessymetrize all
    call DGEMM('N','T',Nbas,Nbas,Nbas,1.0D0,P,Nbas,D,Nbas,0.0D0,A,Nbas)
    call DGEMM('N','N',Nbas,Nbas,Nbas,1.0D0,D,Nbas,A,Nbas,0.0D0,P,Nbas)

    call DGEMM('N','T',Nbas,Nbas,Nbas,1.0D0,S,Nbas,D,Nbas,0.0D0,A,Nbas)
    call DGEMM('N','N',Nbas,Nbas,Nbas,1.0D0,D,Nbas,A,Nbas,0.0D0,S,Nbas)
    
    call DGEMM('N','T',Nbas,Nbas,Nbas,1.0D0,F,Nbas,D,Nbas,0.0D0,A,Nbas)
    call DGEMM('N','N',Nbas,Nbas,Nbas,1.0D0,D,Nbas,A,Nbas,0.0D0,F,Nbas)

  else
    !calculate density matrix in symmetry free case 

    call DGEMM('N','T',Nbas,Nbas,Nbas,1.0D0,n,Nbas,C,Nbas,0.0D0,A,Nbas)
    call DGEMM('N','N',Nbas,Nbas,Nbas,1.0D0,C,Nbas,A,Nbas,0.0D0,P,Nbas)

  endif ! sym

  ! assemble basis set information:

  np = 0
  spn = 0
  ns = 1
  m = 0
  k = 0

  
  allocate(EXPONENTS(size(prim)/2))
  allocate(CONTRACTION(size(prim)/2))
  

  do i = 1,size(prim),2
    np = np+1
    EXPONENTS(np) = prim(i)
    CONTRACTION(np) = prim(i+1)
    n0nums(np) = (i+1)/2
  enddo
  

  allocate(CS(np))
  allocate(CP(np))
  allocate(CD(np))
  allocate(CF(np))
  allocate(CG(np))
  allocate(CH(np))
  allocate(CI(np))
  

  CS = 0.0D0
  CP = 0.0D0
  CD = 0.0D0
  CF = 0.0D0
  CG = 0.0D0
  CH = 0.0D0
  CI = 0.0D0

  ! Construct arrays with contraction coefficients, shells and pointers
  ! with degeneracy info
  ! this part of code was changed and optimized
  
  ns = 0
  k = 0
  Nshell = 0 
  
  
  do i = 1,SIZE(IDSP),3
     if (i+4 < SIZE(IDSP)) then
        if ((IDSP(i) .ne. IDSP(i+3)) .or. (IDSP(i+1) .ne. IDSP(i+4)) .or. & 
           (IDSP(i+2) .ne. IDSP(i+5))) then 
           Nshell = Nshell + 1
         endif 
      else
           Nshell = Nshell + 1 
      exit 
      endif 
  enddo
  
  allocate(basis_ident(5,Nshell))
  
  basis_ident = 0
  basis_ident(5,1)=1
  
  do i = 1,SIZE(IDSP),3
     ns = ns+1
     if (i+4 < SIZE(IDSP)) then
        if ((IDSP(i) .ne. IDSP(i+3)) .or. (IDSP(i+1) .ne. IDSP(i+4)) .or. & 
           (IDSP(i+2) .ne. IDSP(i+5))) then 
           k = k + 1
           basis_ident(1,k) = IDSP(i)
           basis_ident(2,k) = IDSP(i+1)
           basis_ident(3,k) = IDSP(i+2)
           basis_ident(4,k) = ns 
           basis_ident(5,k+1) = basis_ident(5,k)+ns
           ns = 0
         endif 
      else
           k = k + 1 
           basis_ident(1,k) = IDSP(i)
           basis_ident(2,k) = IDSP(i+1)
           basis_ident(3,k) = IDSP(i+2)
           basis_ident(4,k) = ns
      exit 
      endif 
  enddo

  do i = 1,np

    lbas = IDSP(3*n0nums(i)-1) ! basis function ang. mom.

    if (lbas == 0) then
      CS(i) = CONTRACTION(i)
    else if (lbas == 1) then
      CP(i) = CONTRACTION(i)
    else if (lbas == 2) then
      CD(i) = CONTRACTION(i)
    else if (lbas == 3) then
      CF(i) = CONTRACTION(i)
    else if (lbas == 4) then
      CG(i) = CONTRACTION(i)
    else if (lbas == 5) then
      CH(i) = CONTRACTION(i)
    else if (lbas == 6) then
      CI(i) = CONTRACTION(i)
    else
      write(err,*) 'basis exceeds max. angular momentum of 6: ',lbas
      stop 'error termination'
    endif
  enddo

  allocate(nprim(Nbas))
  allocate(npntr(Nbas))
  allocate(ncomps(Nbas))
  
  j = 0 
  do i = 1,SIZE(IDSB),4
     j = j + 1
     do k = 1,Nshell
          if (basis_ident(1,k) == IDSB(i)   .and. &
              basis_ident(2,k) == IDSB(i+2) .and. &
              basis_ident(3,k) == IDSB(i+1)) then 
              nprim(j) = basis_ident(4,k)
              npntr(j) = basis_ident(5,k) 
           endif
     enddo
  enddo

  ncomps = 1

  ! assign NBO FILE-47 input labels to the basis functions:

  Lmax = 0
  label = 0
  k = 0
  do i=1,size(IDSB),4
    k = k +1

    lbas = idsb(i+2) ! angular momentum quantum number of basis fct.
    mbas = idsb(i+3) ! |m| of basis fct., sign indicates cos vs. sin combos

    ! In NBO, for f and higher angular momenta, the function
    ! assignment (spherical) is systematic in the sense that label l51
    ! is for m=0 where l = 3,4,5,6 [e.g., 351 is f0], and labels l52,
    ! l53, l54, l55, etc, are for the cos and sin linear combinations
    ! with |m|=1, |m|=2, etc. [e.g. 351 and 353 is the c1 and s1
    ! linear combination, respectively]. The p- and d-functions follow
    ! a different order. There appears to be on page B-75 of the NBO-7
    ! manual a minor misprint. Function 252=xz corresponds to c1 (not
    ! s1, as written), and function 253=yz corresponds to s1 (not c1,
    ! as written). This has been checked against analytic formulas of
    ! the functions by implementing the recursion formulas from IOData
    ! (https://iodata.readthedocs.io/en/latest/basis.html) in
    ! Mathematica, and by calculating the overlap matrix from the
    ! function assignments below and confirming it corresponds to what
    ! Molcas generates internally.

    if (lbas==0.and.mbas==0) then ! s
      label(k) = 51

    else if (lbas==1.and.mbas==1) then ! px
      label(k) = 151

    else if (lbas==1.and.mbas==-1) then ! py
      label(k) = 152

    else if (lbas==1.and.mbas==0) then ! pz
      label(k) = 153

    else if (lbas==2.and.mbas==-2) then ! dxy = d(-2)
      label(k) = 251

    else if (lbas==2.and.mbas==-1) then ! dyz = d(-1)
      label(k) = 253

    else if (lbas==2.and.mbas==0) then ! dz2 = 2zz-xx-yy = d(0)
      label(k) = 255

    else if (lbas==2.and.mbas==1) then ! dxz = d(+1)
      label(k) = 252

    else if (lbas==2.and.mbas==2) then ! dx2-y2 = xx-yy = d(+2)
      label(k) = 254

    else if (lbas==3.and.mbas==-3) then
      label(k) = 357

    else if (lbas==3.and.mbas==-2) then
      label(k) = 355

    else if (lbas==3.and.mbas==-1) then
      label(k) = 353

    else if (lbas==3.and.mbas==0) then
      label(k) = 351

    else if (lbas==3.and.mbas==1) then
      label(k) = 352

    else if (lbas==3.and.mbas==2) then
      label(k) = 354

    else if (lbas==3.and.mbas==3) then
      label(k) = 356

    else if (lbas==4.and.mbas==-4) then
      label(k) = 459

    else if (lbas==4.and.mbas==-3) then
      label(k) = 457

    else if (lbas==4.and.mbas==-2) then
      label(k) = 455

    else if (lbas==4.and.mbas==-1) then
      label(k) = 453

    else if (lbas==4.and.mbas==0) then
      label(k) = 451

    else if (lbas==4.and.mbas==1) then
      label(k) = 452

    else if (lbas==4.and.mbas==2) then
      label(k) = 454

    else if (lbas==4.and.mbas==3) then
      label(k) = 456

    else if (lbas==4.and.mbas==4) then
      label(k) = 458

    else if (lbas==5.and.mbas==-5) then
      label(k) = 561

    else if (lbas==5.and.mbas==-4) then
      label(k) = 559

    else if (lbas==5.and.mbas==-3) then
      label(k) = 557

    else if (lbas==5.and.mbas==-2) then
      label(k) = 555

    else if (lbas==5.and.mbas==-1) then
      label(k) = 553

    else if (lbas==5.and.mbas==0) then
      label(k) = 551

    else if (lbas==5.and.mbas==1) then
      label(k) = 552

    else if (lbas==5.and.mbas==2) then
      label(k) = 554

    else if (lbas==5.and.mbas==3) then
      label(k) = 556

    else if (lbas==5.and.mbas==4) then
      label(k) = 558

    else if (lbas==5.and.mbas==5) then
      label(k) = 560

    else if (lbas==6.and.mbas==-6) then
      label(k) = 663

    else if (lbas==6.and.mbas==-5) then
      label(k) = 661

    else if (lbas==6.and.mbas==-4) then
      label(k) = 659

    else if (lbas==6.and.mbas==-3) then
      label(k) = 657

    else if (lbas==6.and.mbas==-2) then
      label(k) = 655

    else if (lbas==6.and.mbas==-1) then
      label(k) = 653

    else if (lbas==6.and.mbas==0) then
      label(k) = 651

    else if (lbas==6.and.mbas==1) then
      label(k) = 652

    else if (lbas==6.and.mbas==2) then
      label(k) = 654

    else if (lbas==6.and.mbas==3) then
      label(k) = 656

    else if (lbas==6.and.mbas==4) then
      label(k) = 658

    else if (lbas==6.and.mbas==5) then
      label(k) = 660

    else if (lbas==6.and.mbas==6) then
      label(k) = 662

    else
      ! we ran out of options for the basis. exit with an error
      write(err,*) 'i, k =',i,k
      stop 'cannot assign basis function'
    endif

    if (lbas > Lmax) then
      Lmax = lbas
   endif
 enddo

 ! calculate trace(density matrix * overlap matrix) in desymmetrized (?)
 ! form; must give the number of electrons


 call DGEMM('N','N',Nbas,Nbas,Nbas,1.0D0,P,Nbas,S,Nbas,0.0D0,A,Nbas)

 rtemp = zero
 do i = 1, Nbas
   rtemp = rtemp + A(i,i)
 enddo
 write(*,'(1x,a,f15.6/1x,a/)') 'Total number of electrons from tr[P S]:', &
   & rtemp, 'Do not use File.47 if this result does not look right!'


 ! write the data to an NBO FILE47 formatted file

 open(unit=i47,file='File.47',status="unknown", action="write")


 write(i47,'(A,I8,A,I8,A)') '$GENNBO  NATOMS=',size(atoms),' NBAS='&
   &,Nbas,' BODM BOHR $END'
 write(i47,'(A)') '$NBO AOINFO MULORB BNDIDX NBCP NLMO NPA FILE=NBO&
   &-molcas $END'
 write(i47,'(A)') '$COORD'
 write(i47,'(A,A,F12.8)') trim(h5file), ' Tr[PS]=', rtemp

 k = 0
 do i = 1, size(atoms)
   write(i47,'(I10,1x,I10,3(1x,F15.9))') atoms(i), int(charges(i)), coord(k+1),&
     & coord(k+2), coord(k+3)
   k = k + 3
 enddo
 write(i47,'(A)') '$END'

 write(i47,'(A)') '$BASIS'
 write(i47,'(1x,''CENTER='',10(1x,I7))') (IDSB(i), i = 1,size(IDSB),4)
 write(i47,'(1x,''LABEL='',10(1x,I4))') (label(i), i = 1, size(label))
 write(i47,'(A)') '$END'

 write(i47,'(A)') '$CONTRACT'
 write(i47,*) 'NSHELL = ', Nbas
 write(i47,*) 'NEXP = ', np
 write(i47,'(1x,''NCOMP='',10(1x,I5))') (ncomps(i), i = 1, Nbas)
 write(i47,'(1x,''NPRIM='',10(1x,I5))') (nprim(i), i = 1, Nbas)
 write(i47,'(1x,''NPTR='',10(1x,I5))') (npntr(i), i = 1, Nbas)
 write(i47,'(1x,''EXP='',3(1x,E24.15))') (EXPONENTS(i), i = 1, np)

 ! write contraction coefficients. We assume that the basis has at the
 ! very least a bunch of S functions (l=0), and we assume that if the
 ! basis Lmax has a certain value that all angular momenta up to that
 ! value are represented in the basis. A case with only P or D or F
 ! etc., or S and D but no P and similar situations, would have to be
 ! treated differently but is unlikely to occur.

 write(i47,'(1x,''CS='',3(1x,E24.15))') (CS(i), i = 1, np)

 ! check for angular momenta > 0:

 if (Lmax > 0) then
   write(i47,'(1x,''CP='',3(1x,E24.15))') (CP(i), i = 1, np)
 endif

 if (Lmax > 1) then
   write(i47,'(1x,''CD='',3(1x,E24.15))') (CD(i), i = 1, np)
 endif

 if (Lmax > 2) then
   write(i47,'(1x,''CF='',3(1x,E24.15))') (CF(i), i = 1, np)
 endif

 if (Lmax > 3) then
   write(i47,'(1x,''CG='',3(1x,E24.15))') (CG(i), i = 1, np)
 endif

 if (Lmax > 4) then
   write(i47,'(1x,''CH='',3(1x,E24.15))') (CH(i), i = 1, np)
 endif

 if (Lmax > 5) then
   write(i47,'(1x,''CI='',3(1x,E24.15))') (CI(i), i = 1, np)
 endif

 write(i47,'(A)') '$END'

 write(i47,'(A)') '$OVERLAP'
 write(i47,'(4E24.16)') S
 write(i47,'(A)') '$END'

 write(i47,'(A)') '$DENSITY'
 write(i47,'(4E24.16)') P
 write(i47,'(A)') '$END'
 
 if (save_fock) then 
   write(i47,'(A)') '$FOCK'
   write(i47,'(4E24.16)') F
   write(i47,'(A)') '$END'
 endif  

 ! For LCAMO implementation the C matrices should be 
 ! additionally desymmetrized!
 !write(i47,'(A)') '$LCAOMO'
 !write(i47,'(4E24.16)') C
 !write(i47,'(A)') '$END'


 close(i47)

 deallocate(EXPONENTS)
 deallocate(CONTRACTION)
 deallocate(CS)
 deallocate(CP)
 deallocate(CD)
 deallocate(CF)
 deallocate(CG)
 deallocate(CH)
 deallocate(CI)
 deallocate(overlap_out)
 deallocate(fock_out)
 deallocate(mo_out)
 deallocate(basis)
 deallocate(nprim)
 deallocate(npntr)
 deallocate(ncomps)
 deallocate(basis_ident)

 stop 'normal termination of molcasto47'

end program molcasto47
