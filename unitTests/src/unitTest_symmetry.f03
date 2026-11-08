      program unitTest_symmetry
!
!     unitTest_symmetry: Tests for the symmetry handling of MQC matrices and
!     MQC SCF integrals:
!       T1  MQC_Matrix_Test_Symmetric (full and packed storage, non-square)
!       T2  getSymm (symmetry implied by the storage format)
!       T3  MQC_Matrix_Test_Pair_Symmetric (off-diagonal block pairs)
!       T4  put with the symmetry options (full, packed upper/lower/diagonal)
!       T5  HermMatMul for matrices (value and output storage, forceSym),
!           MQC_Matrix_Packed_Block_Split against MQC_Matrix_Packed_Block_Join
!       T6  mqc_integral_test_symmetric (space/spin/general, missing blocks,
!           non-square dimensions)
!       T7  HermMatMul for integrals
!       T8  mqc_integral_allocate with blocks omitted by keyword
!       T9  integral symmetry flag: setSymm/getSymm, compression keeps every
!           element, copies keep the flag, 'unknown' expands, 'asymmetric'
!       T10 operations on a compressed integral give the same result as on the
!           same integral uncompressed
!       T11 put rules on flagged integrals (diagonal, symm keyword, mismatch),
!           checked against the same writes on the full matrix
!       T12 'asymmetric'/'fullsymmetric'/'diagonal' flags set by allocate,
!           zero blocks not stored, put rules on zero and diagonal integrals,
!           integral%imp (split by the integral dimensions, type promotion,
!           symmetry from the matrix storage) and mqc_matrix_to_integral,
!           allocate with symmetry= and one off-diagonal block (implied block)
!       T13 flag propagation through dagger, transpose, conjg and the scalar
!           products (values, expected flag, flag consistent with the full test)
!       T14 flag propagation through A+B, A-B and A .ewp. B for every pair of
!           flags, and through HermMatMul (values against the full matrices,
!           expected flag, flag consistent with the full test)
!       T15 identity, identity_func, initialize and eigenvalues%diagf with
!           integer, real, complex and mqc_scalar values, spin and rectangular
!           blocks (values, flag, flag consistent with the full test)
!       T16 inverse, power and pseudo-inverse of matrices (value from A^-1 A = I,
!           (A^1/2)^2 = A, repeated products and the Moore-Penrose conditions;
!           output storage) and of flagged and unflagged integrals (values,
!           flag, flag consistent with the full test)
!       T17 swap_offDiagBlocks, pad (matrix and integral), output_orbitals and
!           spinBlockGHF on flagged integrals (values against the unflagged
!           integral or the full matrix, flag, storage)
!     Every test is checked against an independent reference: the
!     full-matrix test, the plain MatMul(Dagger(B),MatMul(A,B)), or the
!     expected value written out explicitly.
!
!     -M. M. F. Moraes, 2026.
!
!     USE Connections
!
      use MQC_EST
      use iso_fortran_env, only: int64, real64
!
      implicit none
      integer::nfail
      integer(kind=int64)::n,m
      real(kind=real64),parameter::tol=1.0e-12,tol2=1.0e-10
      character(len=16),dimension(5),parameter::opts=['symmetric    ', &
        'antisymmetric','hermitian    ','antihermitian','any          ']
!
      type(mqc_matrix)::Z0,E0,X,Y,P,Q,H,AH,S,AS,RS,RA,R,Rect,B,BR,M1,M2
      type(mqc_matrix)::Ha,Hb,Hab,Za,Zb,Tn,Ga,Gb,Kab,Kba,Ba,Bb,Bab,Bba,Hp,Tt,Tz
      type(mqc_scf_integral)::A,BI,AG,AZ,AD,ADr,ADi,AI
      type(mqc_matrix)::G2,Da,Db,Dr,Di,MF2,HH,D2,X7,Hc,Dc
      integer(kind=int64)::k8
      type(mqc_scalar)::sc8
      complex(kind=real64)::z
      integer,parameter::nOps=13
      type(mqc_scf_integral),dimension(nOps)::OPS
      character(len=16),dimension(nOps)::FLG
      character(len=24),dimension(nOps)::NAM
      logical,dimension(nOps)::REALD
      integer::i14,j14
      type(mqc_scf_eigenvalues)::EV
      type(mqc_vector)::Va,Vb
      character(len=16)::evName,evType
      type(mqc_matrix)::PD,PD2,RPD,Hs,Dz,XA
      integer(kind=int64)::sd8
!
 9000 format(/,1x,'Tests completed: ',i0,' failure(s)')
!
      nfail = 0
      n = 4
      m = 3
      call Z0%init(n,n,(0.0,0.0))
      X = cfill(n,n,1)
      Y = cfill(n,n,2)
      H = X + dagger(X)
      AH = X - dagger(X)
      S = X + transpose(X)
      AS = X - transpose(X)
      R = rfill(n,n,3)
      RS = R + transpose(R)
      RA = R - transpose(R)
      Rect = cfill(n,m,4)
!
!--------------------------------------------------------------------
!     T1  MQC_Matrix_Test_Symmetric: expected truth table for
!     (symmetric, antisymmetric, hermitian, antihermitian, any), on
!     full and on packed storage.
!--------------------------------------------------------------------
!
      call t1('T1 complex hermitian',H,[.false.,.false.,.true.,.false.,.true.])
      call t1('T1 complex antihermitian',AH,[.false.,.false.,.false.,.true.,.true.])
      call t1('T1 complex symmetric',S,[.true.,.false.,.false.,.false.,.true.])
      call t1('T1 complex antisymmetric',AS,[.false.,.true.,.false.,.false.,.true.])
      call t1('T1 real symmetric',RS,[.true.,.false.,.true.,.false.,.true.])
      call t1('T1 real antisymmetric',RA,[.false.,.true.,.false.,.true.,.true.])
      call t1('T1 general',X,[.false.,.false.,.false.,.false.,.false.])
      call t1('T1 non-square',Rect,[.false.,.false.,.false.,.false.,.false.])
!
!--------------------------------------------------------------------
!     T2  getSymm
!--------------------------------------------------------------------
!
      call t2('T2 full',X,'unknown')
      call t2('T2 hermitian packed',packed(H),'hermitian')
      call t2('T2 antihermitian packed',packed(AH),'antihermitian')
      call t2('T2 complex symmetric packed',packed(S),'symmetric')
      call t2('T2 complex antisymmetric packed',packed(AS),'antisymmetric')
      call t2('T2 real symmetric packed',packed(RS),'hermitian')
      call t2('T2 real antisymmetric packed',packed(RA),'antihermitian')
      call M1%init(n,n,(0.0,0.0))
      call M1%put((1.0,0.5),1,1)
      call M1%put((2.0,0.0),3,3)
      call mqc_matrix_full2diag(M1)
      call t2('T2 diagonal',M1,'diagonal')
      call t2('T2 unallocated',E0,'unallocated')
!
!--------------------------------------------------------------------
!     T3  MQC_Matrix_Test_Pair_Symmetric(A,B): A = op(B)
!--------------------------------------------------------------------
!
      P = cfill(n,m,5)
      Q = cfill(m,n,6)
      call t3('T3 A = transpose(B)',P,transpose(P),[.true.,.false.,.false.,.false.,.true.])
      call t3('T3 A = -transpose(B)',P,neg(transpose(P)),[.false.,.true.,.false.,.false.,.true.])
      call t3('T3 A = dagger(B)',P,dagger(P),[.false.,.false.,.true.,.false.,.true.])
      call t3('T3 A = -dagger(B)',P,neg(dagger(P)),[.false.,.false.,.false.,.true.,.true.])
      call t3('T3 unrelated',P,Q,[.false.,.false.,.false.,.false.,.false.])
      call t3('T3 B unallocated, A non-zero',P,E0,[.false.,.false.,.false.,.false.,.false.])
      call t3('T3 B unallocated, A zero',Z0,E0,[.true.,.true.,.true.,.true.,.true.])
      call t3('T3 both unallocated',E0,E0,[.true.,.true.,.true.,.true.,.true.])
      call t3('T3 shape mismatch',P,P,[.false.,.false.,.false.,.false.,.false.])
!
!--------------------------------------------------------------------
!     T4  put with the symmetry options
!--------------------------------------------------------------------
!
      z = (1.0,2.0)
      call t4off('T4 packed hermitian upper',H,'hermitian',1,3,'hermitian')
      call t4off('T4 packed hermitian lower',H,'hermitian',3,1,'hermitian')
      call t4off('T4 packed antihermitian upper',AH,'antihermitian',1,3,'antihermitian')
      call t4off('T4 packed antihermitian lower',AH,'antihermitian',3,1,'antihermitian')
      call t4off('T4 packed symmetric upper',S,'symmetric',1,3,'symmetric')
      call t4off('T4 packed symmetric lower',S,'symmetric',3,1,'symmetric')
      call t4off('T4 packed antisymmetric upper',AS,'antisymmetric',1,3,'antisymmetric')
      call t4off('T4 packed antisymmetric lower',AS,'antisymmetric',3,1,'antisymmetric')
      call t4off('T4 packed hermitian, symmetric option',H,'symmetric',1,3,'unknown')
      call t4off('T4 packed symmetric, element option',S,'element',1,3,'unknown')
      call t4diag('T4 packed hermitian diag real',H,(3.0,0.0),'hermitian')
      call t4diag('T4 packed hermitian diag complex',H,z,'unknown')
      call t4diag('T4 packed antihermitian diag imaginary',AH,(0.0,4.0),'antihermitian')
      call t4diag('T4 packed antihermitian diag real',AH,(3.0,0.0),'unknown')
      call t4diag('T4 packed symmetric diag complex',S,z,'symmetric')
      call t4diag('T4 packed antisymmetric diag zero',AS,(0.0,0.0),'antisymmetric')
      call t4diag('T4 packed antisymmetric diag non-zero',AS,z,'unknown')
      call t4full('T4 full element off-diag','element',1,3)
      call t4full('T4 full hermitian off-diag','hermitian',1,3)
      call t4full('T4 full antihermitian off-diag','antihermitian',1,3)
      call t4full('T4 full symmetric off-diag','symmetric',1,3)
      call t4full('T4 full antisymmetric off-diag','antisymmetric',1,3)
      call t4full('T4 full hermitian diag','hermitian',2,2)
      call t4full('T4 full antisymmetric diag','antisymmetric',2,2)
!
!--------------------------------------------------------------------
!     T5  HermMatMul for matrices: value vs MatMul and output storage
!--------------------------------------------------------------------
!
      B = cfill(n,m,7)
      BR = rfill(n,m,8)
      call t5('T5 hermitian full',H,B,'hermitian')
      call t5('T5 hermitian packed',packed(H),B,'hermitian')
      call t5('T5 antihermitian packed',packed(AH),B,'antihermitian')
      call t5('T5 real symmetric packed, real B',packed(RS),BR,'hermitian')
      call t5('T5 real antisymmetric packed, real B',packed(RA),BR,'antihermitian')
      call t5('T5 real symmetric packed, complex B',packed(RS),B,'hermitian')
      call t5('T5 general',X,B,'unknown')
      Hc = H
      do k8 = 1, n
        sc8 = Hc%at(k8,k8)
        call Hc%put(sc8%cval()+(0.0,1.0),k8,k8)
      end do
      call t5('T5 hermitian off-diagonal, complex diagonal, full',Hc,B,'unknown')
      M1 = Hc
      call mqc_matrix_pack(M1,'hermitian')
      call t5('T5 StorHerm with complex diagonal',M1,B,'unknown')
      call Da%init(n,n,(0.0,0.0))
      do k8 = 1, n
        call Da%put(cmplx(real(k8,kind=real64),0.5*k8,kind=real64),k8,k8)
      end do
      M1 = Da
      call mqc_matrix_full2diag(M1)
      call t5('T5 StorDiag with complex diagonal',M1,B,'unknown')
      call t5f('T5 forceSym hermitian, full',H,B,'hermitian')
      call t5f('T5 forceSym antihermitian, full',AH,B,'antihermitian')
      call t5f('T5 forceSym antihermitian, packed',packed(AH),B,'antihermitian')
      call t5f('T5 forceSym antihermitian, real full',RA,BR,'antihermitian')
      call tsplit('T5 split hermitian 4+3',packed(cfill(7,7,20)+dagger(cfill(7,7,20))),4,3)
      call tsplit('T5 split antihermitian 3+4',packed(cfill(7,7,21)-dagger(cfill(7,7,21))),3,4)
      call tsplit('T5 split real symmetric 4+3',packed(rfill(7,7,22)+transpose(rfill(7,7,22))),4,3)
      call tsplit('T5 split complex symmetric 2+5',packed(cfill(7,7,23)+transpose(cfill(7,7,23))),2,5)
      call tsplit('T5 split hermitian 7+0',packed(cfill(7,7,20)+dagger(cfill(7,7,20))),7,0)
!
!--------------------------------------------------------------------
!     T6  mqc_integral_test_symmetric vs the full-matrix test
!--------------------------------------------------------------------
!
      Ha = H
      Hb = Y + dagger(Y)
      Hab = cfill(n,n,9)
      Za = AH
      Zb = Y - dagger(Y)
      Ga = X
      Gb = Y
      Kab = cfill(n,n,10)
      Kba = cfill(n,n,11)
      Tn = neg(dagger(Hab))
      call mqc_integral_allocate(A,'a','space',Ha)
      call t6('T6 space hermitian',A)
      call mqc_integral_allocate(A,'a','space',Za)
      call t6('T6 space antihermitian',A)
      call mqc_integral_allocate(A,'a','space',S)
      call t6('T6 space complex symmetric',A)
      call mqc_integral_allocate(A,'a','space',Ga)
      call t6('T6 space general',A)
      call mqc_integral_allocate(A,'a','spin',Ha,Hb)
      call t6('T6 spin hermitian',A)
      call mqc_integral_allocate(A,'a','spin',Ha,Gb)
      call t6('T6 spin hermitian/general',A)
      call mqc_integral_allocate(A,'a','spin',Za,Zb)
      call t6('T6 spin antihermitian',A)
      call mqc_integral_allocate(A,'a','spin',RS,RS)
      call t6('T6 spin real symmetric',A)
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,dagger(Hab))
      call t6('T6 general hermitian',A)
      call mqc_integral_allocate(A,'a','general',Za,Zb,Hab,Tn)
      call t6('T6 general antihermitian',A)
      call mqc_integral_allocate(A,'a','general',S,S,Hab,transpose(Hab))
      call t6('T6 general complex symmetric',A)
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,Kab)
      call t6('T6 general, bad off-diagonal',A)
      call mqc_integral_allocate(A,'a','general',Ga,Gb,Kab,Kba)
      call t6('T6 general, no symmetry',A)
      call mqc_integral_allocate(A,'a','general',Ha,Hb,alphaBeta=Hab)
      call t6('T6 general, only alphaBeta (non-zero)',A)
      call mqc_integral_allocate(A,'a','general',Ha,Hb,alphaBeta=Z0)
      call t6('T6 general, only alphaBeta (zero)',A)
      call mqc_integral_allocate(A,'a','general',Ha,Hb,betaAlpha=Z0)
      call t6('T6 general, only betaAlpha (zero)',A)
      call mqc_integral_allocate(A,'a','general',Ha,Hb)
      call t6('T6 general, no off-diagonal blocks',A)
      Hp = packed(Ha)
      Tz = Hab + dagger(Hab)
      Tt = packed(Tz)
      call mqc_integral_allocate(A,'a','general',Hp,Hb,Tt,Tz)
      call t6('T6 general, packed blocks',A)
      call mqc_integral_allocate(A,'a','spin',nDimA1=n,nDimA2=m,nDimB1=n,nDimB2=m)
      call t6dims('T6 spin, no blocks, non-square dims',A,.false.)
      call mqc_integral_allocate(A,'a','spin',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      call t6dims('T6 spin, no blocks, square dims',A,.true.)
!
!--------------------------------------------------------------------
!     T7  HermMatMul for integrals vs MatMul(Dagger(B),MatMul(A,B))
!--------------------------------------------------------------------
!
      Ba = cfill(n,m,12)
      Bb = cfill(n,m,13)
      Bab = cfill(n,m,14)
      Bba = cfill(n,m,15)
      call mqc_integral_allocate(A,'a','space',Ha)
      call mqc_integral_allocate(BI,'b','space',Ba)
      call t7('T7 space x space',A,BI)
      call mqc_integral_allocate(A,'a','spin',Ha,Hb)
      call mqc_integral_allocate(BI,'b','spin',Ba,Bb)
      call t7('T7 spin x spin',A,BI)
      call mqc_integral_allocate(A,'a','space',Ha)
      call t7('T7 space x spin',A,BI)
      call mqc_integral_allocate(A,'a','spin',Ha,Hb)
      call mqc_integral_allocate(BI,'b','space',Ba)
      call t7('T7 spin x space',A,BI)
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,dagger(Hab))
      call mqc_integral_allocate(BI,'b','general',Ba,Bb,Bab,Bba)
      call t7('T7 general x general',A,BI)
      call mqc_integral_allocate(A,'a','spin',Ha,Hb)
      call t7('T7 spin x general',A,BI)
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,dagger(Hab))
      call mqc_integral_allocate(BI,'b','spin',Ba,Bb)
      call t7('T7 general x spin',A,BI)
      call mqc_integral_allocate(A,'a','general',Ga,Gb,Kab,Kba)
      call mqc_integral_allocate(BI,'b','general',Ba,Bb,Bab,Bba)
      call t7('T7 non-hermitian general',A,BI)
      call mqc_integral_allocate(A,'a','spin',Za,Zb)
      call mqc_integral_allocate(BI,'b','spin',Ba,Bb)
      call t7('T7 antihermitian spin',A,BI)
      call mqc_integral_allocate(A,'a','general',Za,Zb,Hab,Tn)
      call mqc_integral_allocate(BI,'b','general',Ba,Bb,Bab,Bba)
      call t7('T7 antihermitian general',A,BI)
      call mqc_integral_allocate(A,'a','space',RS)
      call mqc_integral_allocate(BI,'b','space',rfill(n,m,16))
      call t7('T7 real space x space',A,BI)
!
!--------------------------------------------------------------------
!     T8  mqc_integral_allocate with blocks omitted by keyword
!--------------------------------------------------------------------
!
      call mqc_integral_allocate(A,'a','general',Ha,Hb,betaAlpha=Hab)
      call check('T8 general, betaAlpha only by keyword', &
        A%blockSize('alpha',1).eq.n.and.A%blockSize('beta',2).eq.n)
      call mqc_integral_allocate(A,'a','spin',beta=Hb,nDimA1=n,nDimA2=n)
      call check('T8 spin, beta only by keyword', &
        A%blockSize('alpha',1).eq.n.and.A%blockSize('beta',1).eq.n)
      call mqc_integral_allocate(A,'a','space',nDimA1=n,nDimA2=m)
      call check('T8 space, dimensions only', &
        A%blockSize('alpha',1).eq.n.and.A%blockSize('alpha',2).eq.m)
!
!
!--------------------------------------------------------------------
!     T9  integral symmetry flag
!--------------------------------------------------------------------
!
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,dagger(Hab))
      call check('T9 allocate -> unknown',A%getSymm().eq.'unknown')
      call t9('T9 general hermitian',A,'hermitian')
      call mqc_integral_allocate(A,'a','general',Ha,Hb,betaAlpha=Z0)
      call t9('T9 general hermitian, only a zero betaAlpha stored',A,'hermitian')
      call mqc_integral_allocate(A,'a','general',Za,Zb,Hab,Tn)
      call t9('T9 general antihermitian',A,'antihermitian')
      call mqc_integral_allocate(A,'a','general',S,S,Hab,transpose(Hab))
      call t9('T9 general symmetric',A,'symmetric')
      call mqc_integral_allocate(A,'a','general',AS,AS,Hab,neg(transpose(Hab)))
      call t9('T9 general antisymmetric',A,'antisymmetric')
      call mqc_integral_allocate(A,'a','spin',Ha,Hb)
      call t9('T9 spin hermitian',A,'hermitian')
      call mqc_integral_allocate(A,'a','spin',RS,RS)
      call t9('T9 spin real symmetric as hermitian',A,'hermitian')
      call t9('T9 spin real symmetric as symmetric',A,'symmetric')
      call mqc_integral_allocate(A,'a','space',Ha)
      call t9('T9 space hermitian',A,'hermitian')
      call mqc_integral_allocate(A,'a','general',Ga,Gb,Kab,Kba)
      call A%setSymm('asymmetric')
      call check('T9 general asymmetric -> asymmetric',A%getSymm().eq.'asymmetric')
      call check('T9 general asymmetric: test any is false', &
        .not.mqc_integral_test_symmetric(A,'any'))
!
!--------------------------------------------------------------------
!     T10 operations on compressed integrals
!--------------------------------------------------------------------
!
      call mqc_integral_allocate(BI,'b','general',Ga,Gb,Kab,Kba)
      G2 = cfill(2*n,2*n,17)
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,dagger(Hab))
      call t10('T10 general hermitian',A,'hermitian')
      call mqc_integral_allocate(A,'a','general',Za,Zb,Hab,Tn)
      call t10('T10 general antihermitian',A,'antihermitian')
      call mqc_integral_allocate(A,'a','general',S,S,Hab,transpose(Hab))
      call t10('T10 general symmetric',A,'symmetric')
      call mqc_integral_allocate(A,'a','general',AS,AS,Hab,neg(transpose(Hab)))
      call t10('T10 general antisymmetric',A,'antisymmetric')
      call mqc_integral_allocate(A,'a','spin',Ha,Hb)
      call t10('T10 spin hermitian',A,'hermitian')
      call mqc_integral_allocate(A,'a','space',Ha)
      call t10('T10 space hermitian',A,'hermitian')
!
!--------------------------------------------------------------------
!     T11 put rules on flagged integrals (full indices; alpha = 1..4,
!     beta = 5..8)
!--------------------------------------------------------------------
!
      call mqc_integral_allocate(AG,'a','general',Ha,Hb,Hab,dagger(Hab))
      call t11('T11 herm: real alpha diagonal',AG,'hermitian',(3.0,0.0),2,2,'','hermitian')
      call t11('T11 herm: complex alpha diagonal',AG,'hermitian',z,2,2,'','unknown')
      call t11('T11 herm: real beta diagonal',AG,'hermitian',(3.0,0.0),6,6,'','hermitian')
      call t11('T11 herm: alpha off-diagonal, no symm',AG,'hermitian',z,1,3,'','unknown')
      call t11('T11 herm: alpha off-diagonal, herm',AG,'hermitian',z,1,3,'herm','hermitian')
      call t11('T11 herm: alpha lower, herm',AG,'hermitian',z,3,1,'hermitian','hermitian')
      call t11('T11 herm: beta off-diagonal, herm',AG,'hermitian',z,5,7,'herm','hermitian')
      call t11('T11 herm: alphaBeta element, herm',AG,'hermitian',z,6,2,'herm','hermitian')
      call t11('T11 herm: betaAlpha element, herm',AG,'hermitian',z,2,6,'herm','hermitian')
      call t11('T11 herm: alphaBeta (1,1), herm',AG,'hermitian',z,5,1,'herm','hermitian')
      call t11('T11 herm: alphaBeta element, no symm',AG,'hermitian',z,6,2,'','unknown')
      call t11('T11 herm: off-diagonal, sym (mismatch)',AG,'hermitian',z,1,3,'sym','unknown')
      call t11('T11 herm: complex diagonal, herm',AG,'hermitian',z,2,2,'herm','unknown')
      call mqc_integral_allocate(AG,'a','general',Za,Zb,Hab,Tn)
      call t11('T11 antiherm: imaginary diagonal',AG,'antihermitian',(0.0,4.0),2,2,'','antihermitian')
      call t11('T11 antiherm: real diagonal',AG,'antihermitian',(3.0,0.0),2,2,'','unknown')
      call t11('T11 antiherm: off-diagonal, antiherm',AG,'antihermitian',z,1,3,'antiherm','antihermitian')
      call t11('T11 antiherm: alphaBeta element, antiherm',AG,'antihermitian',z,6,2,'antiherm', &
        'antihermitian')
      call mqc_integral_allocate(AG,'a','general',S,S,Hab,transpose(Hab))
      call t11('T11 sym: complex diagonal',AG,'symmetric',z,2,2,'','symmetric')
      call t11('T11 sym: off-diagonal, sym',AG,'symmetric',z,1,3,'sym','symmetric')
      call t11('T11 sym: betaAlpha element, sym',AG,'symmetric',z,2,6,'sym','symmetric')
      call mqc_integral_allocate(AG,'a','general',AS,AS,Hab,neg(transpose(Hab)))
      call t11('T11 antisym: zero diagonal',AG,'antisymmetric',(0.0,0.0),2,2,'','antisymmetric')
      call t11('T11 antisym: non-zero diagonal',AG,'antisymmetric',z,2,2,'','unknown')
      call t11('T11 antisym: off-diagonal, antisym',AG,'antisymmetric',z,1,3,'antisym','antisymmetric')
      call mqc_integral_allocate(AG,'a','spin',Ha,Hb)
      call t11('T11 spin herm: real diagonal',AG,'hermitian',(3.0,0.0),6,6,'','hermitian')
      call t11('T11 spin herm: alphaBeta element, herm',AG,'hermitian',z,6,2,'herm','hermitian')
      call mqc_integral_allocate(AG,'a','space',Ha)
      call t11('T11 space herm: off-diagonal, herm',AG,'hermitian',z,1,3,'herm','hermitian')
      call t11('T11 space herm: beta diagonal',AG,'hermitian',(3.0,0.0),6,6,'','unknown')
      call mqc_integral_allocate(AG,'a','general',Ga,Gb,Kab,Kba)
      call t11('T11 asymmetric: diagonal',AG,'asymmetric',(3.0,0.0),2,2,'','unknown')
      call mqc_integral_allocate(AG,'a','general',Ha,Hb,Hab,dagger(Hab))
      call t11('T11 unknown: off-diagonal, herm',AG,'unknown',z,1,3,'herm','unknown')
!
!--------------------------------------------------------------------
!     T12 case 5 flags, zero blocks, zero/diagonal put rules and imp
!--------------------------------------------------------------------
!
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,dagger(Hab))
      call check('T12 allocate, square with blocks -> unknown',A%getSymm().eq.'unknown')
      call mqc_integral_allocate(A,'a','spin',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      call check('T12 allocate, square, no blocks -> fullsymmetric',A%getSymm().eq.'fullsymmetric')
      call mqc_integral_allocate(A,'a','general',Z0,Z0,Z0,Z0)
      call check('T12 allocate, square, zero blocks -> fullsymmetric',A%getSymm().eq.'fullsymmetric')
      call mqc_integral_allocate(A,'a','spin',Ba,Bb)
      call check('T12 allocate, non-square -> asymmetric',A%getSymm().eq.'asymmetric')
      call A%setSymm('unknown')
      call check('T12 non-square: unknown keeps asymmetric',A%getSymm().eq.'asymmetric')
      call mqc_integral_allocate(A,'a','spin',nDimA1=n,nDimA2=m,nDimB1=n,nDimB2=m)
      call check('T12 allocate, non-square, no blocks -> asymmetric',A%getSymm().eq.'asymmetric')
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Z0,Z0)
      call mqc_integral_allocate(BI,'b','general',Ha,Hb)
      call check('T12 zero off-diagonal blocks dropped, values kept', &
        A%getSymm().eq.'unknown'.and.maxdiff(A%getBlock('full'),BI%getBlock('full')).lt.tol)
      call mqc_set_very_small(-1.0)
      call mqc_integral_allocate(A,'a','general',Z0,Z0,Z0,Z0)
      call check('T12 negative MQC_very_small: zero blocks kept -> unknown',A%getSymm().eq.'unknown')
      call mqc_set_very_small(1.0e-12)
!
      call timpl('T12 implied: herm, only alphaBeta',Ha,Hb,'hermitian',Hab,E0,Hab,dagger(Hab), &
        'hermitian')
      call timpl('T12 implied: herm, only betaAlpha',Ha,Hb,'hermitian',E0,dagger(Hab),Hab,dagger(Hab), &
        'hermitian')
      call timpl('T12 implied: antiherm, only alphaBeta',Za,Zb,'antihermitian',Hab,E0,Hab,Tn, &
        'antihermitian')
      call timpl('T12 implied: antiherm, only betaAlpha',Za,Zb,'antihermitian',E0,Tn,Hab,Tn, &
        'antihermitian')
      call timpl('T12 implied: sym, only betaAlpha',S,S,'symmetric',E0,transpose(Hab),Hab, &
        transpose(Hab),'symmetric')
      call timpl('T12 implied: antisym, only alphaBeta',AS,AS,'antisymmetric',Hab,E0,Hab, &
        neg(transpose(Hab)),'antisymmetric')
      call timpl('T12 implied: real herm, only alphaBeta',RS,RS,'hermitian',R,E0,R,transpose(R), &
        'hermitian')
      call timpl('T12 implied: herm fails, only alphaBeta',Ga,Gb,'hermitian',Hab,E0,Hab,dagger(Hab), &
        'unknown')
      call timpl('T12 implied: herm fails, only betaAlpha',Ga,Gb,'hermitian',E0,dagger(Hab),Hab, &
        dagger(Hab),'unknown')
      call timpl('T12 implied: herm, zero alphaBeta given',Ha,Hb,'hermitian',Z0,E0,Z0,Z0,'hermitian')
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,Kab,symmetry='hermitian')
      call check('T12 both off-diagonal blocks given, inconsistent -> unknown', &
        A%getSymm().eq.'unknown'.and.maxdiff(A%getBlock('beta-alpha'),Kab).lt.tol)
      call mqc_integral_allocate(A,'a','general',Ba,Bb,alphaBeta=Rect,symmetry='hermitian')
      call check('T12 implied: non-square, missing block zero -> asymmetric', &
        A%getSymm().eq.'asymmetric'.and.isZeroMat(A%getBlock('beta-alpha')))
!
      call Da%init(n,n,(0.0,0.0))
      call Db%init(n,n,(0.0,0.0))
      call Dr%init(n,n,0.0)
      do k8 = 1, n
        call Da%put(cmplx(real(k8,kind=real64),0.5*k8,kind=real64),k8,k8)
        call Db%put(cmplx(2.0*k8,-1.0*k8,kind=real64),k8,k8)
        call Dr%put(real(3*k8,kind=real64),k8,k8)
      end do
      call mqc_integral_allocate(AD,'d','general',Da,Db)
      call t9('T12 general complex diagonal',AD,'complexdiag')
      call mqc_integral_allocate(ADr,'d','general',Dr,Dr)
      call t9('T12 general real diagonal',ADr,'realdiag')
      call Di%init(n,n,(0.0,0.0))
      do k8 = 1, n
        call Di%put(cmplx(0.0,real(k8,kind=real64),kind=real64),k8,k8)
      end do
      call mqc_integral_allocate(ADi,'d','general',Di,Di)
      call t9('T12 general imaginary diagonal',ADi,'imagdiag')
      call mqc_integral_allocate(AZ,'z','general',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      call t9('T12 general zero',AZ,'fullsymmetric')
      call t1diag('T12 matrix test, complex diagonal',Da)
      call t1diag('T12 matrix test, real diagonal',Dr)
      call t1diag('T12 matrix test, imaginary diagonal',Di)
!
      call t11('T12 zero: zero element',AZ,'fullsymmetric',(0.0,0.0),2,2,'','fullsymmetric')
      call t11('T12 zero: plain alpha diagonal',AZ,'fullsymmetric',z,2,2,'','complexdiag')
      call t11('T12 zero: plain beta diagonal',AZ,'fullsymmetric',z,6,6,'','complexdiag')
      call t11('T12 zero: plain real diagonal',AZ,'fullsymmetric',(3.0,0.0),2,2,'','realdiag')
      call t11('T12 zero: plain imaginary diagonal',AZ,'fullsymmetric',(0.0,3.0),2,2,'','imagdiag')
      call t11('T12 zero: plain off-diagonal',AZ,'fullsymmetric',z,1,3,'','unknown')
      call t11('T12 zero: off-diagonal, herm',AZ,'fullsymmetric',z,1,3,'herm','hermitian')
      call t11('T12 zero: alphaBeta element, herm',AZ,'fullsymmetric',z,6,2,'herm','hermitian')
      call t11('T12 zero: complex diagonal, herm',AZ,'fullsymmetric',z,2,2,'herm','unknown')
      call t11('T12 zero: complex diagonal, sym',AZ,'fullsymmetric',z,2,2,'sym','symmetric')
      call t11('T12 diagonal: plain diagonal',AD,'complexdiag',z,2,2,'','complexdiag')
      call t11('T12 real diagonal: plain real diagonal',ADr,'realdiag',(5.0,0.0),2,2,'','realdiag')
      call t11('T12 real diagonal: plain complex diagonal',ADr,'realdiag',z,2,2,'','complexdiag')
      call t11('T12 diagonal: plain off-diagonal',AD,'complexdiag',z,1,3,'','unknown')
      call t11('T12 diagonal: off-diagonal herm, complex diagonal',AD,'complexdiag',z,1,3,'herm', &
        'unknown')
      call t11('T12 diagonal: off-diagonal sym, complex diagonal',AD,'complexdiag',z,1,3,'sym', &
        'symmetric')
      call t11('T12 real diagonal: off-diagonal herm',ADr,'realdiag',z,1,3,'herm','hermitian')
      call t11('T12 real diagonal: alphaBeta element, herm',ADr,'realdiag',z,6,2,'herm','hermitian')
      call t11('T12 imaginary diagonal: off-diagonal antiherm',ADi,'imagdiag',z,1,3,'antiherm', &
        'antihermitian')
      call t11('T12 imaginary diagonal: off-diagonal herm',ADi,'imagdiag',z,1,3,'herm','unknown')
!
      MF2 = cfill(2*n,2*n,18)
      call mqc_integral_allocate(AI,'i','general',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      call timp('T12 imp general',AI,MF2,'general','unknown')
      call MF2%classify()
      call check('T12 imp general source classified asymmetric',&
        MF2%getSymm().eq.'asymmetric')
      call timp('T12 imp classified asymmetric',AI,MF2,'general','unknown')
      call mqc_integral_allocate(AI,'i','spin',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      call timp('T12 imp spin promoted to general',AI,MF2,'general','unknown')
      call mqc_integral_allocate(A,'a','spin',Ha,Hb)
      call mqc_integral_allocate(AI,'i','space',nDimA1=n,nDimA2=n)
      call timp('T12 imp space promoted to spin',AI,A%getBlock('full'),'spin','')
      call mqc_integral_allocate(A,'a','spin',Ha,Ha)
      call timp('T12 imp space, beta = alpha',AI,A%getBlock('full'),'space','')
      HH = MF2 + dagger(MF2)
      call mqc_integral_allocate(AI,'i','general',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      call timp('T12 imp packed hermitian matrix',AI,packed(HH),'general','hermitian')
      call D2%init(2*n,2*n,(0.0,0.0))
      do k8 = 1, 2*n
        call D2%put(cmplx(real(k8,kind=real64),1.0,kind=real64),k8,k8)
      end do
      call mqc_matrix_full2diag(D2)
      call timp('T12 imp diagonal matrix',AI,D2,'general','complexdiag')
      X7 = cfill(7,7,19)
      call mqc_integral_allocate(AI,'i','general',nDimA1=4,nDimA2=4, &
        nDimB1=3,nDimB2=3)
      call timp('T12 imp uneven alpha/beta',AI,X7,'general','unknown')
      call mqc_matrix_to_integral(MF2,AI,'i','spin',n)
      call check('T12 matrix_to_integral spin promoted to general', &
        AI%type().eq.'general'.and.maxdiff(AI%getBlock('full'),MF2).lt.tol)
      call tblock('T12 getBlock storage, general hermitian',Ha,Hb,Hab,dagger(Hab),'hermitian','hermitian')
      call tblock('T12 getBlock storage, general antihermitian',Za,Zb,Hab,Tn,'antihermitian', &
        'antihermitian')
      call tblock('T12 getBlock storage, general symmetric',S,S,Hab,transpose(Hab),'symmetric','symmetric')
      call mqc_integral_allocate(A,'d','general',Da,Db)
      call A%setSymm('complexdiag')
      call check('T12 getBlock storage, complexdiag full -> diagonal', &
        mqc_matrix_get_symmetry(A%getBlock('full')).eq.'diagonal')
      call check('T12 getBlock storage, complexdiag alpha -> diagonal', &
        mqc_matrix_get_symmetry(A%getBlock('alpha')).eq.'diagonal')
      call check('T12 getBlock storage, complexdiag full nosym -> unknown', &
        mqc_matrix_get_symmetry(A%getBlock('full',nosym=.true.)).eq.'unknown')
      call mqc_integral_allocate(BI,'b','general',Da,Db)
      call check('T12 getBlock complexdiag full values', &
        maxdiff(A%getBlock('full'),BI%getBlock('full')).lt.tol)
      call mqc_integral_allocate(BI,'b','general',Ba,Bb,Bab,Bba)
      call t7('T7b HermMatMul complexdiag x general',A,BI)
      call mqc_integral_allocate(A,'d','general',Dr,Dr)
      call A%setSymm('realdiag')
      call t7('T7b HermMatMul realdiag x general',A,BI)
      call mqc_integral_allocate(A,'d','general',Di,Di)
      call A%setSymm('imagdiag')
      call t7('T7b HermMatMul imagdiag x general',A,BI)
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,dagger(Hab))
      call A%setSymm('hermitian')
      call t7('T7b HermMatMul compressed hermitian x general',A,BI)
      call mqc_integral_allocate(A,'a','spin',Ha,Hb)
      call A%setSymm('hermitian')
      call t7('T7b HermMatMul compressed spin hermitian x general',A,BI)
      call mqc_matrix_to_integral(packed(H),AI,'i','space')
      call check('T12 matrix_to_integral space, packed hermitian', &
        AI%type().eq.'space'.and.AI%getSymm().eq.'hermitian'.and. &
        maxdiff(AI%getBlock('alpha'),H).lt.tol)
!
!--------------------------------------------------------------------
!     T13 flag propagation (dagger, transpose, conjg, scalar products)
!--------------------------------------------------------------------
!
      call mqc_integral_allocate(A,'a','general',Ha,Hb,Hab,dagger(Hab))
      call t13('T13 hermitian complex',A,'hermitian')
      call mqc_integral_allocate(A,'a','spin',RS,RS)
      call t13('T13 hermitian real',A,'hermitian')
      call mqc_integral_allocate(A,'a','general',Za,Zb,Hab,Tn)
      call t13('T13 antihermitian complex',A,'antihermitian')
      call mqc_integral_allocate(A,'a','spin',RA,RA)
      call t13('T13 antihermitian real',A,'antihermitian')
      call mqc_integral_allocate(A,'a','general',S,S,Hab,transpose(Hab))
      call t13('T13 symmetric',A,'symmetric')
      call mqc_integral_allocate(A,'a','general',AS,AS,Hab,neg(transpose(Hab)))
      call t13('T13 antisymmetric',A,'antisymmetric')
      call mqc_integral_allocate(A,'a','general',Dr,Dr)
      call t13('T13 realdiag',A,'realdiag')
      call mqc_integral_allocate(A,'a','general',Di,Di)
      call t13('T13 imagdiag',A,'imagdiag')
      call mqc_integral_allocate(A,'a','general',Da,Db)
      call t13('T13 complexdiag',A,'complexdiag')
      call mqc_integral_allocate(A,'a','general',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      call t13('T13 fullsymmetric',A,'fullsymmetric')
      call mqc_integral_allocate(A,'a','general',Ga,Gb,Kab,Kba)
      call t13('T13 asymmetric',A,'asymmetric')
      call mqc_integral_allocate(A,'a','space',Ha)
      call t13('T13 space hermitian',A,'hermitian')
      call Dc%init(n,n,(0.0,0.0))
      do k8 = 1, n
        call Dc%put(cmplx(real(k8,kind=real64),real(k8,kind=real64),kind=real64),k8,k8)
      end do
      call mqc_integral_allocate(A,'a','general',Dc,Dc)
      call A%setSymm('complexdiag')
      AG = (1.0,-1.0)*A
      call check('T13 complexdiag (1+i)k times (1-i) -> realdiag',AG%getSymm().eq.'realdiag')
!
!--------------------------------------------------------------------
!     T14 flag propagation (A+B, A-B, A .ewp. B, HermMatMul)
!--------------------------------------------------------------------
!
      call mqc_integral_allocate(OPS(1),'o','general',Ha,Hb,Hab,dagger(Hab))
      FLG(1) = 'hermitian'
      NAM(1) = 'herm'
      call mqc_integral_allocate(OPS(2),'o','spin',RS,RS)
      FLG(2) = 'hermitian'
      NAM(2) = 'real herm (spin)'
      call mqc_integral_allocate(OPS(3),'o','general',Za,Zb,Hab,Tn)
      FLG(3) = 'antihermitian'
      NAM(3) = 'antiherm'
      call mqc_integral_allocate(OPS(4),'o','spin',RA,RA)
      FLG(4) = 'antihermitian'
      NAM(4) = 'real antiherm (spin)'
      call mqc_integral_allocate(OPS(5),'o','general',S,S,Hab,transpose(Hab))
      FLG(5) = 'symmetric'
      NAM(5) = 'sym'
      call mqc_integral_allocate(OPS(6),'o','general',AS,AS,Hab,neg(transpose(Hab)))
      FLG(6) = 'antisymmetric'
      NAM(6) = 'antisym'
      call mqc_integral_allocate(OPS(7),'o','general',Dr,Dr)
      FLG(7) = 'realdiag'
      NAM(7) = 'realdiag'
      call mqc_integral_allocate(OPS(8),'o','general',Di,Di)
      FLG(8) = 'imagdiag'
      NAM(8) = 'imagdiag'
      call mqc_integral_allocate(OPS(9),'o','general',Da,Db)
      FLG(9) = 'complexdiag'
      NAM(9) = 'complexdiag'
      call mqc_integral_allocate(OPS(10),'o','general',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      FLG(10) = 'fullsymmetric'
      NAM(10) = 'zero'
      call mqc_integral_allocate(OPS(11),'o','general',Ga,Gb,Kab,Kba)
      FLG(11) = 'asymmetric'
      NAM(11) = 'asymmetric'
      call mqc_integral_allocate(OPS(12),'o','general',Ha,Hb,Hab,dagger(Hab))
      FLG(12) = 'unknown'
      NAM(12) = 'unflagged herm'
      call mqc_integral_allocate(OPS(13),'o','space',Ha)
      FLG(13) = 'hermitian'
      NAM(13) = 'herm (space)'
      do i14 = 1, nOps
        if(FLG(i14).ne.'unknown') call OPS(i14)%setSymm(FLG(i14))
        REALD(i14) = .not.mqc_matrix_havecomplex(OPS(i14)%getBlock('full',nosym=.true.))
      end do
      do i14 = 1, nOps
        do j14 = 1, nOps
          call t14(i14,j14)
        end do
      end do
      AG = mqc_integral_multiply(OPS(1),OPS(5))
      AI = OPS(1) .ewp. OPS(5)
      call check('T14 mqc_integral_multiply = .ewp.', &
        maxdiff(AG%getBlock('full',nosym=.true.),AI%getBlock('full',nosym=.true.)).lt.tol.and. &
        AG%getSymm().eq.AI%getSymm())
!
      call mqc_integral_allocate(BI,'b','general',Ba,Bb,Bab,Bba)
      call t14h('T14 HermMatMul herm',OPS(1),BI,'hermitian')
      call t14h('T14 HermMatMul unflagged herm',OPS(12),BI,'hermitian')
      call t14h('T14 HermMatMul real herm (spin)',OPS(2),BI,'hermitian')
      call t14h('T14 HermMatMul antiherm',OPS(3),BI,'antihermitian')
      call mqc_integral_allocate(A,'a','general',Za,Zb,Hab,Tn)
!     Without a flag only the Hermitian structure of a full A is detected.
      call t14h('T14 HermMatMul unflagged antiherm',A,BI,'unknown')
      call A%setSymm('antihermitian')
      call t14h('T14 HermMatMul antiherm flagged after setSymm',A,BI,'antihermitian')
      call t14h('T14 HermMatMul sym',OPS(5),BI,'unknown')
      call t14h('T14 HermMatMul realdiag',OPS(7),BI,'hermitian')
      call t14h('T14 HermMatMul imagdiag',OPS(8),BI,'antihermitian')
      call t14h('T14 HermMatMul complexdiag',OPS(9),BI,'unknown')
      call t14h('T14 HermMatMul zero',OPS(10),BI,'fullsymmetric')
      call t14h('T14 HermMatMul asymmetric',OPS(11),BI,'unknown')
      call mqc_integral_allocate(BI,'b','spin',Ba,Bb)
      call t14h('T14 HermMatMul herm x spin',OPS(1),BI,'hermitian')
      call t14h('T14 HermMatMul real herm (spin) x spin',OPS(2),BI,'hermitian')
      call mqc_integral_allocate(BI,'b','space',Ba)
      call t14h('T14 HermMatMul herm (space) x space',OPS(13),BI,'hermitian')
!
!--------------------------------------------------------------------
!     T15 identity, identity_func and initialize
!--------------------------------------------------------------------
!
      call AI%identity(4,4)
      call t15('T15 identity default',AI,ref15(4,4,4,4,(1.0,0.0),.true.),'realdiag')
      call AI%identity(4,4,scalar=2.0)
      call t15('T15 identity real 2',AI,ref15(4,4,4,4,(2.0,0.0),.true.),'realdiag')
      call AI%identity(4,4,scalar=3)
      call t15('T15 identity integer 3',AI,ref15(4,4,4,4,(3.0,0.0),.true.),'realdiag')
      call AI%identity(4,4,scalar=(2.0,0.0))
      call t15('T15 identity complex (2,0)',AI,ref15(4,4,4,4,(2.0,0.0),.true.),'realdiag')
      call AI%identity(4,4,scalar=(0.0,2.0))
      call t15('T15 identity complex (0,2)',AI,ref15(4,4,4,4,(0.0,2.0),.true.),'imagdiag')
      call AI%identity(4,4,scalar=(1.0,1.0))
      call t15('T15 identity complex (1,1)',AI,ref15(4,4,4,4,(1.0,1.0),.true.),'complexdiag')
      sc8 = (0.0,2.0)
      call AI%identity(4,4,scalar=sc8)
      call t15('T15 identity mqc_scalar (0,2)',AI,ref15(4,4,4,4,(0.0,2.0),.true.),'imagdiag')
      sc8 = 2.0
      call AI%identity(4,4,scalar=sc8)
      call t15('T15 identity mqc_scalar real 2',AI,ref15(4,4,4,4,(2.0,0.0),.true.),'realdiag')
      call AI%identity(4,3,scalar=(0.0,2.0))
      call t15('T15 identity spin 4/3 (0,2)',AI,ref15(4,4,3,3,(0.0,2.0),.true.),'imagdiag')
      call AI%identity(4,3)
      call t15('T15 identity spin 4/3 default',AI,ref15(4,4,3,3,(1.0,0.0),.true.),'realdiag')
      call AI%identity(4,4,scalar=0.0)
      call t15('T15 identity scalar 0',AI,ref15(4,4,4,4,(0.0,0.0),.true.),'fullsymmetric')
      call AI%identity(4,4,nAlpha2=3,nBeta2=3)
      call t15('T15 identity rectangular',AI,ref15(4,3,4,3,(1.0,0.0),.true.),'asymmetric')
      AI = mqc_integral_identity_func(4,4)
      call t15('T15 identity_func square',AI,ref15(4,4,4,4,(1.0,0.0),.true.),'realdiag')
      AI = mqc_integral_identity_func(4,3,scalar=(0.0,2.0))
      call t15('T15 identity_func spin (0,2)',AI,ref15(4,4,3,3,(0.0,2.0),.true.),'imagdiag')
      AI = mqc_integral_identity_func(4,3,nAlpha2=4,nBeta2=2)
      call t15('T15 identity_func rectangular',AI,ref15(4,4,3,2,(1.0,0.0),.true.),'asymmetric')
      call AI%init(4,4)
      call t15('T15 init default',AI,ref15(4,4,4,4,(0.0,0.0),.false.),'fullsymmetric')
      call AI%init(4,4,scalar=3.0)
      call t15('T15 init real 3',AI,ref15(4,4,4,4,(3.0,0.0),.false.),'hermitian')
      call AI%init(4,4,scalar=2)
      call t15('T15 init integer 2',AI,ref15(4,4,4,4,(2.0,0.0),.false.),'hermitian')
      call AI%init(4,4,scalar=(3.0,0.0))
      call t15('T15 init complex (3,0)',AI,ref15(4,4,4,4,(3.0,0.0),.false.),'hermitian')
      call AI%init(4,4,scalar=(0.0,3.0))
      call t15('T15 init complex (0,3)',AI,ref15(4,4,4,4,(0.0,3.0),.false.),'antihermitian')
      call AI%init(4,4,scalar=(1.0,2.0))
      call t15('T15 init complex (1,2)',AI,ref15(4,4,4,4,(1.0,2.0),.false.),'symmetric')
      sc8 = (0.0,3.0)
      call AI%init(4,4,scalar=sc8)
      call t15('T15 init mqc_scalar (0,3)',AI,ref15(4,4,4,4,(0.0,3.0),.false.),'antihermitian')
      call AI%init(4,3,scalar=(0.0,3.0))
      call t15('T15 init spin 4/3 (0,3)',AI,ref15(4,4,3,3,(0.0,3.0),.false.),'antihermitian')
      call AI%init(4,4,scalar=0.0)
      call t15('T15 init scalar 0',AI,ref15(4,4,4,4,(0.0,0.0),.false.),'fullsymmetric')
      call AI%init(4,4,scalar=1.0,nAlpha2=3,nBeta2=3)
      call t15('T15 init rectangular',AI,ref15(4,3,4,3,(1.0,0.0),.false.),'asymmetric')
      call AI%init(4,3,scalar=(1.0,2.0),nAlpha2=4,nBeta2=2)
      call t15('T15 init spin, rectangular beta',AI,ref15(4,4,3,2,(1.0,2.0),.false.),'asymmetric')
!
!     eigenvalues%diagf (mqc_eigenvalues_to_integral); the names are variables since
!     mqc_eigenvalues_allocate changes the case of its arguments in place.
      evName = ''
      evType = 'space'
      Va = evec([(1.0,0.0),(-2.0,0.0),(3.0,0.0),(0.5,0.0)],.true.)
      call mqc_eigenvalues_allocate(EV,evName,evType,Va)
      AI = EV%diagf()
      M1 = MQC_Matrix_DiagMatrix_Put_Vector_Func(Va)
      call mqc_integral_allocate(BI,'b','space',M1)
      call t15('T15 diagf real space',AI,BI%getBlock('full',nosym=.true.),'realdiag')
      call check('T15 diagf real space: alpha StorDiag', &
        mqc_matrix_get_symmetry(AI%getBlock('alpha')).eq.'diagonal')
      Va = evec([(1.0,1.0),(-2.0,0.5),(3.0,0.0),(0.5,-1.0)],.false.)
      call mqc_eigenvalues_allocate(EV,evName,evType,Va)
      AI = EV%diagf()
      M1 = MQC_Matrix_DiagMatrix_Put_Vector_Func(Va)
      call mqc_integral_allocate(BI,'b','space',M1)
      call t15('T15 diagf complex space',AI,BI%getBlock('full',nosym=.true.),'complexdiag')
      Va = evec([(0.0,1.0),(0.0,-2.0),(0.0,3.0),(0.0,0.5)],.false.)
      call mqc_eigenvalues_allocate(EV,evName,evType,Va)
      AI = EV%diagf()
      M1 = MQC_Matrix_DiagMatrix_Put_Vector_Func(Va)
      call mqc_integral_allocate(BI,'b','space',M1)
      call t15('T15 diagf imaginary space',AI,BI%getBlock('full',nosym=.true.),'imagdiag')
      evType = 'spin'
      Va = evec([(1.0,0.0),(-2.0,0.0),(3.0,0.0),(0.5,0.0)],.true.)
      Vb = evec([(4.0,0.0),(5.0,0.0),(-6.0,0.0)],.true.)
      call mqc_eigenvalues_allocate(EV,evName,evType,Va,Vb)
      AI = EV%diagf()
      M1 = MQC_Matrix_DiagMatrix_Put_Vector_Func(Va)
      M2 = MQC_Matrix_DiagMatrix_Put_Vector_Func(Vb)
      call mqc_integral_allocate(BI,'b','spin',M1,M2)
      call t15('T15 diagf real spin 4/3',AI,BI%getBlock('full',nosym=.true.),'realdiag')
      Vb = evec([(4.0,1.0),(5.0,0.0),(-6.0,0.0)],.false.)
      call mqc_eigenvalues_allocate(EV,evName,evType,Va,Vb)
      AI = EV%diagf()
      M2 = MQC_Matrix_DiagMatrix_Put_Vector_Func(Vb)
      call mqc_integral_allocate(BI,'b','spin',M1,M2)
      call t15('T15 diagf spin, complex beta',AI,BI%getBlock('full',nosym=.true.),'complexdiag')
      evType = 'space'
      Va = evec([(0.0,0.0),(0.0,0.0),(0.0,0.0),(0.0,0.0)],.true.)
      call mqc_eigenvalues_allocate(EV,evName,evType,Va)
      AI = EV%diagf()
      M1 = MQC_Matrix_DiagMatrix_Put_Vector_Func(Va)
      call mqc_integral_allocate(BI,'b','space',M1)
      call t15('T15 diagf zero',AI,BI%getBlock('full',nosym=.true.),'fullsymmetric')
!
!--------------------------------------------------------------------
!     T16 inverse, power and pseudo-inverse
!--------------------------------------------------------------------
!
      PD = hermPD(n,24)
      PD2 = hermPD(2*n,25)
      RPD = matmul(transpose(R),R) + eyeR(n)
      Hs = matmul(cfill(n,1,26),dagger(cfill(n,1,26)))
      call Dz%init(n,n,(0.0,0.0))
      call Dz%put((2.0,1.0),1,1)
      call Dz%put((0.0,-3.0),3,3)
      call mqc_matrix_full2diag(Dz)
      XA = AH + eyeC(n)
!     Matrix inverse.
      call tinv('T16 inv hermitian packed',packed(PD),'hermitian')
      call tinv('T16 inv real symmetric packed',packed(RPD),'hermitian')
      call tinv('T16 inv complex symmetric packed',packed(S),'symmetric')
      call tinv('T16 inv antihermitian packed',packed(AH),'antihermitian')
      call tinv('T16 inv real antisymmetric packed',packed(RA),'antihermitian')
      call tinv('T16 inv complex antisymmetric packed',packed(AS),'antisymmetric')
      M1 = Da
      call mqc_matrix_full2diag(M1)
      call tinv('T16 inv complex diagonal',M1,'diagonal')
      M1 = Dr
      call mqc_matrix_full2diag(M1)
      call tinv('T16 inv real diagonal',M1,'diagonal')
!     The array assignment after LU packs a result that tests symmetric.
      call tinv('T16 inv hermitian full',PD,'hermitian')
      call tinv('T16 inv general full',XA,'asymmetric')
!     Matrix power.
      call tpowh('T16 power 1/2 hermitian packed',packed(PD),'hermitian')
      call tpowh('T16 power 1/2 hermitian full',PD,'hermitian')
      call tpowh('T16 power 1/2 real symmetric packed',packed(RPD),'hermitian')
      call tpowh('T16 power 1/2 hermitian, negative eigenvalues',packed(H),'unknown')
      call tpowh('T16 power 1/2 real symmetric, negative eigenvalues',packed(RS),'unknown')
      call tpowh('T16 power 1/2 complex diagonal',Dz,'diagonal')
      call tpowk('T16 power 2 general',X,2,'unknown')
      call tpowk('T16 power 3 general',XA,3,'unknown')
      call tpowk('T16 power 2 complex symmetric',S,2,'unknown')
      call tpowk('T16 power 2 hermitian packed',packed(H),2,'hermitian')
      call tpowk('T16 power 2 antihermitian packed',packed(AH),2,'hermitian')
      call tpowk('T16 power 3 antihermitian packed',packed(AH),3,'antihermitian')
      call tpowk('T16 power 3 antihermitian full',AH,3,'antihermitian')
      call tpowk('T16 power -1 hermitian packed',packed(PD),-1,'hermitian')
      call tpowk('T16 power 2 complex diagonal',Dz,2,'diagonal')
!     Matrix pseudo-inverse.
      call tpinv('T16 pinv hermitian packed, full rank',packed(PD),n,'hermitian')
      call tpinv('T16 pinv hermitian packed, rank 1',packed(Hs),1,'hermitian')
      call tpinv('T16 pinv hermitian full, rank 1',Hs,1,'hermitian')
      call tpinv('T16 pinv antihermitian packed',packed(AH),n,'antihermitian')
      call tpinv('T16 pinv real symmetric packed',packed(RPD),n,'hermitian')
      call tpinv('T16 pinv diagonal with zeros',Dz,2,'diagonal')
      call tpinv('T16 pinv general (SVD)',X,n,'unknown')
      call tpinv('T16 pinv complex symmetric (SVD)',S,n,'unknown')
!     Integrals.
      call mqc_integral_allocate(AI,'i','general',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      call AI%imp(PD2)
      call AI%setSymm('hermitian')
      call tint9('T16 general hermitian (flagged)',AI,'hermitian','hermitian','hermitian')
      call mqc_integral_allocate(AG,'i','general',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      call AG%imp(PD2)
      call tint9('T16 general hermitian (unflagged)',AG,'hermitian','hermitian','hermitian')
      call mqc_integral_allocate(AI,'i','spin',PD,hermPD(n,27))
      call AI%setSymm('hermitian')
      call tint9('T16 spin hermitian',AI,'hermitian','hermitian','hermitian')
      call mqc_integral_allocate(AI,'i','space',packed(PD))
      call AI%setSymm('hermitian')
      call tint9('T16 space hermitian',AI,'hermitian','hermitian','hermitian')
      call mqc_integral_allocate(AI,'i','general',Dr,Dr)
      call AI%setSymm('realdiag')
      call tint9('T16 realdiag',AI,'realdiag','realdiag','realdiag')
      call mqc_integral_allocate(AI,'i','general',Da,Db)
      call AI%setSymm('complexdiag')
      call tint9('T16 complexdiag',AI,'complexdiag','complexdiag','complexdiag')
      call mqc_integral_allocate(AI,'i','general',Di,Di)
      call AI%setSymm('imagdiag')
      call tint9('T16 imagdiag',AI,'imagdiag','complexdiag','imagdiag')
      call mqc_integral_allocate(AI,'i','general',Za,Zb,Hab,Tn)
      call AI%setSymm('antihermitian')
      call tint9('T16 general antihermitian',AI,'antihermitian','unknown','antihermitian')
      AG = AI
      call AG%power(2)
      AZ = MatMul(AI,AI)
      call t15('T16 general antihermitian: power 2',AG,AZ%getBlock('full',nosym=.true.),'hermitian')
      AG = AI
      call AG%power(3)
      AZ = MatMul(AZ,AI)
      call check('T16 general antihermitian: power 3: values', &
        maxdiff(AG%getBlock('full',nosym=.true.),AZ%getBlock('full',nosym=.true.)).lt.tol2)
      call t15('T16 general antihermitian: power 3',AG,AG%getBlock('full',nosym=.true.),'antihermitian')
      call mqc_integral_allocate(AI,'i','general',nDimA1=n,nDimA2=n,nDimB1=n,nDimB2=n)
      AG = AI%pinv(spaceDim=sd8)
      call check('T16 zero: pinv -> fullsymmetric, spaceDim 0',AG%getSymm().eq.'fullsymmetric'.and.sd8.eq.0)
!
!--------------------------------------------------------------------
!     T17 swap_offDiagBlocks, pad, output_orbitals, spinBlockGHF
!--------------------------------------------------------------------
!
      do i14 = 1, nOps
        call t17swap(i14)
      end do
      Bab = cfill(n,m,14)
      Bba = cfill(n,m,15)
      call mqc_integral_allocate(AI,'i','general',Ba,Bb,Bab,Bba)
      AG = mqc_integral_swap_offDiagBlocks(AI)
      call check('T17 swap rectangular same-shape blocks: values and flag', &
        maxdiff(AG%getBlock('alpha-beta'),Bba).lt.tol.and.maxdiff(AG%getBlock('beta-alpha'),Bab).lt.tol.and. &
        maxdiff(AG%getBlock('alpha'),Ba).lt.tol.and.AG%getSymm().eq.'asymmetric')
      call t17pad('T17 pad full 6, herm',OPS(1),6,6,'full','hermitian')
      call t17pad('T17 pad alpha 6, herm',OPS(1),6,6,'alpha','hermitian')
      call t17pad('T17 pad beta 6, antiherm',OPS(3),6,6,'beta','antihermitian')
      call t17pad('T17 pad full 6, sym',OPS(5),6,6,'full','symmetric')
      call t17pad('T17 pad full 6, realdiag',OPS(7),6,6,'full','realdiag')
      call t17pad('T17 pad full 6, complexdiag',OPS(9),6,6,'full','complexdiag')
      call t17pad('T17 pad full 6, real herm (spin)',OPS(2),6,6,'full','hermitian')
      call t17pad('T17 pad alpha 6, herm (space)',OPS(13),6,6,'alpha','hermitian')
      call t17pad('T17 pad full 6x5, herm',OPS(1),6,5,'full','asymmetric')
      call t17pad('T17 pad alpha 6x5, herm',OPS(1),6,5,'alpha','asymmetric')
      call t17pad('T17 pad full 6, unflagged',OPS(12),6,6,'full','unknown')
      call t17pad('T17 pad full 6, zero',OPS(10),6,6,'full','fullsymmetric')
      M1 = packed(H)
      M2 = M1%mpad(6_int64,6_int64)
      call check('T17 matrix pad, packed hermitian 6x6: storage', &
        M2%getSymm().eq.'hermitian'.and.MQC_Matrix_HaveComplex(M2))
      call check('T17 matrix pad, packed hermitian 6x6: values', &
        maxdiff(M2%mat([1,4],[1,4]),H).lt.tol.and.isZeroMat(M2%mat([5,6],[1,6])))
      M2 = H%mpad(6_int64,5_int64)
      call check('T17 matrix pad, complex full 6x5: type and values', &
        MQC_Matrix_HaveComplex(M2).and.maxdiff(M2%mat([1,4],[1,4]),H).lt.tol)
      M1 = Dz
      M2 = M1%mpad(6_int64,6_int64)
      call check('T17 matrix pad, diagonal 6x6: storage',M2%getSymm().eq.'diagonal')
      do i14 = 1, nOps
        AI = OPS(i14)
        AG = OPS(i14)
        call AG%setSymm('unknown')
        AZ = mqc_integral_output_orbitals(AI,'complex')
        BI = mqc_integral_output_orbitals(AG,'complex')
        call check('T17 output_orbitals complex, '//trim(NAM(i14)), &
          maxdiff(AZ%getBlock('full',nosym=.true.),BI%getBlock('full',nosym=.true.)).lt.tol)
      end do
      call mqc_integral_allocate(AI,'i','general',Ha,Hb,Hab,dagger(Hab))
      AG = AI
      call AG%setSymm('hermitian')
      call mqc_integral_spinBlockGHF(AI,4_int64,1_int64)
      call mqc_integral_spinBlockGHF(AG,4_int64,1_int64)
      call check('T17 spinBlockGHF flagged = unflagged', &
        maxdiff(AG%getBlock('full',nosym=.true.),AI%getBlock('full',nosym=.true.)).lt.tol)
!
!
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1
!
      contains
!
!     Deterministic complex and real fills (no random numbers, so the
!     output is reproducible).
!
      function cfill(r,c,seed) result(Mat)
      integer(kind=int64),intent(in)::r,c
      integer,intent(in)::seed
      type(mqc_matrix)::Mat
      integer(kind=int64)::i,j
      call Mat%init(r,c,(0.0,0.0))
      do j = 1, c
        do i = 1, r
          call Mat%put(cmplx(sin(real(i+2*j+seed,kind=real64)), &
            cos(real(3*i-j+2*seed,kind=real64)),kind=real64),i,j)
        end do
      end do
      end function cfill
!
      function rfill(r,c,seed) result(Mat)
      integer(kind=int64),intent(in)::r,c
      integer,intent(in)::seed
      type(mqc_matrix)::Mat
      integer(kind=int64)::i,j
      call Mat%init(r,c,0.0)
      do j = 1, c
        do i = 1, r
          call Mat%put(sin(real(2*i+j+seed,kind=real64)),i,j)
        end do
      end do
      end function rfill
!
      function neg(Mat) result(MatOut)
      type(mqc_matrix),intent(in)::Mat
      type(mqc_matrix)::MatOut,Zero
      call Zero%init(size(Mat,1),size(Mat,2),(0.0,0.0))
      MatOut = Zero - Mat
      end function neg
!
      function packed(Mat) result(MatOut)
      type(mqc_matrix),intent(in)::Mat
      type(mqc_matrix)::MatOut
      MatOut = Mat
      call mqc_matrix_full2symm(MatOut)
      end function packed
!
      function val(Mat,i,j) result(c)
      type(mqc_matrix),intent(in)::Mat
      integer,intent(in)::i,j
      complex(kind=real64)::c
      type(mqc_scalar)::sc
      sc = Mat%at(i,j)
      c = sc%cval()
      end function val
!
      function maxdiff(M1,M2) result(d)
      type(mqc_matrix),intent(in)::M1,M2
      real(kind=real64)::d
      type(mqc_matrix)::Dm
      type(mqc_scalar)::sc
      Dm = M1 - M2
      sc = Dm%norm('M')
      d = abs(sc%cval())
      end function maxdiff
!
      function mirror(op,c) result(mc)
      character(len=*),intent(in)::op
      complex(kind=real64),intent(in)::c
      complex(kind=real64)::mc
      select case (op)
      case ('hermitian')
        mc = conjg(c)
      case ('antihermitian')
        mc = -conjg(c)
      case ('symmetric')
        mc = c
      case ('antisymmetric')
        mc = -c
      case default
        mc = (0.0,0.0)
      end select
      end function mirror
!
      subroutine check(label,ok)
      character(len=*),intent(in)::label
      logical,intent(in)::ok
      if(ok) then
        write(*,'(1x,a,a)') 'PASS  ',trim(label)
      else
        write(*,'(1x,a,a)') 'FAIL  ',trim(label)
        nfail = nfail + 1
      endIf
      end subroutine check
!
      subroutine t1(label,Mat,expected)
      character(len=*),intent(in)::label
      type(mqc_matrix),intent(in)::Mat
      logical,dimension(5),intent(in)::expected
      type(mqc_matrix)::Mp
      integer::k
      do k = 1, 5
        call check(label//' full: '//trim(opts(k)), &
          mqc_matrix_test_symmetric(Mat,trim(opts(k))).eqv.expected(k))
      end do
      if(size(Mat,1).ne.size(Mat,2)) return
      Mp = packed(Mat)
      do k = 1, 5
        call check(label//' packed: '//trim(opts(k)), &
          mqc_matrix_test_symmetric(Mp,trim(opts(k))).eqv.expected(k))
      end do
      end subroutine t1
!
      subroutine t2(label,Mat,expected)
      character(len=*),intent(in)::label,expected
      type(mqc_matrix),intent(in)::Mat
      call check(label//' -> '//expected,Mat%getSymm().eq.expected)
      end subroutine t2
!
      subroutine t3(label,MatA,MatB,expected)
      character(len=*),intent(in)::label
      type(mqc_matrix),intent(in)::MatA,MatB
      logical,dimension(5),intent(in)::expected
      integer::k
      do k = 1, 5
        call check(label//': '//trim(opts(k)), &
          mqc_matrix_test_pair_symmetric(MatA,MatB,trim(opts(k))).eqv.expected(k))
      end do
      end subroutine t3
!
      subroutine t4off(label,Mat,op,i,j,storeAfter)
      character(len=*),intent(in)::label,op,storeAfter
      type(mqc_matrix),intent(in)::Mat
      integer,intent(in)::i,j
      type(mqc_matrix)::Mp
      logical::ok
      Mp = packed(Mat)
      call Mp%put(z,i,j,op)
      ok = abs(val(Mp,i,j)-z).lt.tol.and.Mp%getSymm().eq.storeAfter
      if(op.ne.'element') ok = ok.and.abs(val(Mp,j,i)-mirror(op,z)).lt.tol
      call check(label,ok)
      end subroutine t4off
!
      subroutine t4diag(label,Mat,c,storeAfter)
      character(len=*),intent(in)::label,storeAfter
      type(mqc_matrix),intent(in)::Mat
      complex(kind=real64),intent(in)::c
      type(mqc_matrix)::Mp
      Mp = packed(Mat)
      call Mp%put(c,2,2)
      call check(label,abs(val(Mp,2,2)-c).lt.tol.and.Mp%getSymm().eq.storeAfter)
      end subroutine t4diag
!
      subroutine t4full(label,op,i,j)
      character(len=*),intent(in)::label,op
      integer,intent(in)::i,j
      type(mqc_matrix)::Mf
      logical::ok
      call Mf%init(n,n,(0.0,0.0))
      call Mf%put(z,i,j,op)
!     a zero matrix (StorZero) takes the symmetry of its first write with a flag
      ok = abs(val(Mf,i,j)-z).lt.tol.and.(Mf%getSymm().eq.'unknown'.or.(op.ne.'element'.and. &
        (Mf%getSymm().eq.op.or.(op.eq.'symmetric'.and.Mf%getSymm().eq.'hermitian'))))
      if(i.ne.j) ok = ok.and.abs(val(Mf,j,i)-mirror(op,z)).lt.tol
      call check(label,ok)
      end subroutine t4full
!
      subroutine t5(label,MatA,MatB,storeAfter)
      character(len=*),intent(in)::label,storeAfter
      type(mqc_matrix),intent(in)::MatA,MatB
      type(mqc_matrix)::C1,C2
      C1 = HermMatMul(MatA,MatB)
      C2 = MatMul(Dagger(MatB),MatMul(MatA,MatB))
      call check(label//' value',maxdiff(C1,C2).lt.tol)
      call check(label//' storage -> '//storeAfter,C1%getSymm().eq.storeAfter)
      end subroutine t5
!
      subroutine t5f(label,MatA,MatB,storeAfter)
      character(len=*),intent(in)::label,storeAfter
      type(mqc_matrix),intent(in)::MatA,MatB
      type(mqc_matrix)::C1,C2
      C1 = HermMatMul(MatA,MatB,forceSym=storeAfter)
      C2 = MatMul(Dagger(MatB),MatMul(MatA,MatB))
      call check(label//' value',maxdiff(C1,C2).lt.tol)
      call check(label//' storage -> '//storeAfter,C1%getSymm().eq.storeAfter)
      end subroutine t5f
!
      subroutine tsplit(label,Mat,na,nb)
      character(len=*),intent(in)::label
      type(mqc_matrix),intent(in)::Mat
      integer,intent(in)::na,nb
      type(mqc_matrix)::Ma,Mab,Mb,Mj,Mf
      integer(kind=int64)::na8,nb8
      logical::ok
      na8 = na
      nb8 = nb
      call MQC_Matrix_Packed_Block_Split(Mat,na8,nb8,Ma,Mab,Mb)
      Mf = Mat
      call mqc_matrix_symm2full(Mf)
      ok = .true.
      if(na.gt.0) ok = ok.and.Ma%getSymm().eq.Mat%getSymm().and. &
        maxdiff(Ma,Mf%mat([1,na8],[1,na8])).lt.tol
      if(nb.gt.0) ok = ok.and.Mb%getSymm().eq.Mat%getSymm().and. &
        maxdiff(Mb,Mf%mat([na8+1,na8+nb8],[na8+1,na8+nb8])).lt.tol
!     the off-diagonal block is non-square: StorNoSym ('asymmetric')
      if(na.gt.0.and.nb.gt.0) ok = ok.and.Mab%getSymm().eq.'asymmetric'.and. &
        maxdiff(Mab,Mf%mat([na8+1,na8+nb8],[1,na8])).lt.tol
      call check(label//': blocks',ok)
      Mj = MQC_Matrix_Packed_Block_Join(Ma,Mab,Mb,na8,nb8,Mat%getSymm())
      call check(label//': join round trip',Mj%getSymm().eq.Mat%getSymm().and.maxdiff(Mj,Mat).lt.tol)
      end subroutine tsplit
!
      subroutine t6(label,Int)
      character(len=*),intent(in)::label
      type(mqc_scf_integral),intent(in)::Int
      type(mqc_matrix)::F
      integer::k
      F = Int%getBlock('full')
      do k = 1, 5
        call check(label//': '//trim(opts(k)), &
          mqc_integral_test_symmetric(Int,trim(opts(k))).eqv. &
          mqc_matrix_test_symmetric(F,trim(opts(k))))
      end do
      end subroutine t6
!
      subroutine t6dims(label,Int,expected)
      character(len=*),intent(in)::label
      type(mqc_scf_integral),intent(in)::Int
      logical,intent(in)::expected
      integer::k
      do k = 1, 4
        call check(label//': '//trim(opts(k)), &
          mqc_integral_test_symmetric(Int,trim(opts(k))).eqv.expected)
      end do
      end subroutine t6dims
!
      subroutine t7(label,IntA,IntB)
      character(len=*),intent(in)::label
      type(mqc_scf_integral),intent(in)::IntA,IntB
      type(mqc_scf_integral)::C1,C2
      type(mqc_matrix)::F1,F2
      C1 = HermMatMul(IntA,IntB)
      C2 = MatMul(Dagger(IntB),MatMul(IntA,IntB))
      F1 = C1%getBlock('full')
      F2 = C2%getBlock('full')
      call check(label//' value',maxdiff(F1,F2).lt.tol)
      call check(label//' type',C1%type().eq.C2%type())
      end subroutine t7
!
      subroutine t9(label,Int,sym)
      character(len=*),intent(in)::label,sym
      type(mqc_scf_integral),intent(in)::Int
      type(mqc_scf_integral)::IC,ICopy
      type(mqc_matrix)::F0,F1,B0
      type(mqc_scalar)::s0,s1
      integer(kind=int64)::i,j
      integer::k
      logical::ok
      character(len=16),dimension(4),parameter::blocks=['alpha    ','beta     ', &
        'alphabeta','betaalpha']
      F0 = Int%getBlock('full')
      IC = Int
      call IC%setSymm(sym)
      call check(label//': flag -> '//sym,IC%getSymm().eq.sym)
      F1 = IC%getBlock('full')
      call check(label//': full matrix kept',maxdiff(F0,F1).lt.tol)
      do k = 1, 4
        if(Int%type().ne.'general'.and.k.gt.2) exit
        B0 = Int%getBlock(trim(blocks(k)))
        ok = .true.
        do j = 1, size(B0,2)
          do i = 1, size(B0,1)
            s0 = B0%at(i,j)
            s1 = IC%at(i,j,trim(blocks(k)))
            if(abs(s0%cval()-s1%cval()).gt.tol) ok = .false.
          end do
        end do
        call check(label//': every '//trim(blocks(k))//' element by at',ok)
      end do
      do k = 1, 5
        call check(label//': test '//trim(opts(k)), &
          mqc_integral_test_symmetric(IC,trim(opts(k))).eqv. &
          mqc_matrix_test_symmetric(F0,trim(opts(k))))
        call check(label//': test math '//trim(opts(k)), &
          mqc_integral_test_symmetric(IC,trim(opts(k)),.true.).eqv. &
          mqc_matrix_test_symmetric(F0,trim(opts(k)),.true.))
      end do
      ICopy = IC
      call check(label//': copy keeps the flag',ICopy%getSymm().eq.sym)
      call ICopy%setSymm('unknown')
      F1 = ICopy%getBlock('full')
      call check(label//': unknown expands, values kept', &
        ICopy%getSymm().eq.'unknown'.and.maxdiff(F0,F1).lt.tol)
      end subroutine t9
!
      subroutine t10(label,Int,sym)
      character(len=*),intent(in)::label,sym
      type(mqc_scf_integral),intent(in)::Int
      type(mqc_scf_integral)::IU,IC,R1,R2
      type(mqc_scalar)::sc,c1,c2
      IU = Int
      IC = Int
      call IC%setSymm(sym)
      sc = (2.0,1.0)
      R1 = dagger(IC)
      R2 = dagger(IU)
      call cmpInt(label//': dagger',R1,R2)
      R1 = transpose(IC)
      R2 = transpose(IU)
      call cmpInt(label//': transpose',R1,R2)
      R1 = conjg(IC)
      R2 = conjg(IU)
      call cmpInt(label//': conjg',R1,R2)
      R1 = sc*IC
      R2 = sc*IU
      call cmpInt(label//': scalar*integral',R1,R2)
      R1 = IC*sc
      R2 = IU*sc
      call cmpInt(label//': integral*scalar',R1,R2)
      R1 = IC + BI
      R2 = IU + BI
      call cmpInt(label//': A+B',R1,R2)
      R1 = BI + IC
      R2 = BI + IU
      call cmpInt(label//': B+A',R1,R2)
      R1 = IC - BI
      R2 = IU - BI
      call cmpInt(label//': A-B',R1,R2)
      R1 = BI - IC
      R2 = BI - IU
      call cmpInt(label//': B-A',R1,R2)
      R1 = MatMul(IC,BI)
      R2 = MatMul(IU,BI)
      call cmpInt(label//': MatMul(A,B)',R1,R2)
      R1 = MatMul(BI,IC)
      R2 = MatMul(BI,IU)
      call cmpInt(label//': MatMul(B,A)',R1,R2)
      R1 = MatMul(IC,G2)
      R2 = MatMul(IU,G2)
      call cmpInt(label//': MatMul(A,matrix)',R1,R2)
      R1 = MatMul(G2,IC)
      R2 = MatMul(G2,IU)
      call cmpInt(label//': MatMul(matrix,A)',R1,R2)
      R1 = HermMatMul(IC,BI)
      R2 = HermMatMul(IU,BI)
      call cmpInt(label//': HermMatMul(A,B)',R1,R2)
      R1 = IC .ewp. BI
      R2 = IU .ewp. BI
      call cmpInt(label//': A .ewp. B',R1,R2)
      c1 = Contraction(IC,BI)
      c2 = Contraction(IU,BI)
      call check(label//': Contraction(A,B)',abs(c1%cval()-c2%cval()).lt.tol2)
      call check(label//': getBlock ab', &
        maxdiff(IC%getBlock('alpha-beta'),IU%getBlock('alpha-beta')).lt.tol)
      call check(label//': getBlock ba', &
        maxdiff(IC%getBlock('beta-alpha'),IU%getBlock('beta-alpha')).lt.tol)
      R1 = IC
      call R1%updateType()
      R2 = IU
      call R2%updateType()
      call cmpInt(label//': updateType',R1,R2)
      R1 = IC%changeType('general')
      R2 = IU%changeType('general')
      call cmpInt(label//': changeType general',R1,R2)
      end subroutine t10
!
      subroutine cmpInt(label,Int1,Int2)
      character(len=*),intent(in)::label
      type(mqc_scf_integral),intent(in)::Int1,Int2
      call check(label,maxdiff(Int1%getBlock('full'),Int2%getBlock('full')).lt.tol2.and. &
        Int1%type().eq.Int2%type())
      end subroutine cmpInt
!
      subroutine t11(label,Int,sym,c,i,j,symmOpt,flagAfter)
      character(len=*),intent(in)::label,sym,symmOpt,flagAfter
      type(mqc_scf_integral),intent(in)::Int
      complex(kind=real64),intent(in)::c
      integer,intent(in)::i,j
      type(mqc_scf_integral)::IC
      type(mqc_matrix)::F
      type(mqc_scalar)::sc
      character(len=16)::matOpt
      IC = Int
      if(sym.ne.'unknown') call IC%setSymm(sym)
      F = IC%getBlock('full')
      sc = c
      if(symmOpt.eq.'') then
        call IC%put(sc,i,j)
        call F%put(sc,i,j)
      else
        select case (symmOpt)
        case ('herm')
          matOpt = 'hermitian'
        case ('antiherm')
          matOpt = 'antihermitian'
        case ('sym')
          matOpt = 'symmetric'
        case ('antisym')
          matOpt = 'antisymmetric'
        case default
          matOpt = symmOpt
        end select
        call IC%put(sc,i,j,symm=symmOpt)
        call F%put(sc,i,j,trim(matOpt))
        if(Int%type().eq.'space'.and.i.le.n.and.j.le.n) &
          call F%put(sc,i+n,j+n,trim(matOpt))
      endIf
      call check(label//': values',maxdiff(IC%getBlock('full'),F).lt.tol)
      call check(label//': flag -> '//flagAfter,IC%getSymm().eq.flagAfter)
      end subroutine t11
!
      subroutine t1diag(label,Mat)
      character(len=*),intent(in)::label
      type(mqc_matrix),intent(in)::Mat
      type(mqc_matrix)::Md
      integer::k
      Md = Mat
      call mqc_matrix_full2diag(Md)
      do k = 1, 5
        call check(label//': StorDiag = full, '//trim(opts(k)), &
          mqc_matrix_test_symmetric(Md,trim(opts(k))).eqv.mqc_matrix_test_symmetric(Mat,trim(opts(k))))
        call check(label//': StorDiag = full, math, '//trim(opts(k)), &
          mqc_matrix_test_symmetric(Md,trim(opts(k)),.true.).eqv. &
          mqc_matrix_test_symmetric(Mat,trim(opts(k)),.true.))
      end do
      end subroutine t1diag
!
      subroutine tblock(label,Ma,Mb,Mab,Mba,sym,storage)
      character(len=*),intent(in)::label,sym,storage
      type(mqc_matrix),intent(in)::Ma,Mb,Mab,Mba
      type(mqc_scf_integral)::IC,IU
      call mqc_integral_allocate(IU,'u','general',Ma,Mb,Mab,Mba)
      IC = IU
      call IC%setSymm(sym)
      call check(label//': full storage -> '//storage, &
        mqc_matrix_get_symmetry(IC%getBlock('full')).eq.storage)
      call check(label//': full values',maxdiff(IC%getBlock('full'),IU%getBlock('full')).lt.tol)
      call check(label//': alpha storage -> '//storage, &
        mqc_matrix_get_symmetry(IC%getBlock('alpha')).eq.storage)
      call check(label//': alpha values',maxdiff(IC%getBlock('alpha'),IU%getBlock('alpha')).lt.tol)
      call check(label//': beta storage -> '//storage, &
        mqc_matrix_get_symmetry(IC%getBlock('beta')).eq.storage)
      call check(label//': alphaBeta storage -> unknown', &
        mqc_matrix_get_symmetry(IC%getBlock('alpha-beta')).eq.'unknown')
      call check(label//': full nosym storage -> unknown', &
        mqc_matrix_get_symmetry(IC%getBlock('full',nosym=.true.)).eq.'unknown')
      call check(label//': full nosym values', &
        maxdiff(IC%getBlock('full',nosym=.true.),IU%getBlock('full')).lt.tol)
      call check(label//': alpha nosym storage -> unknown', &
        mqc_matrix_get_symmetry(IC%getBlock('alpha',nosym=.true.)).eq.'unknown')
      call check(label//': alpha nosym values', &
        maxdiff(IC%getBlock('alpha',nosym=.true.),IU%getBlock('alpha')).lt.tol)
      end subroutine tblock
!
      subroutine t13(label,Int,sym)
      character(len=*),intent(in)::label,sym
      type(mqc_scf_integral),intent(in)::Int
      type(mqc_scf_integral)::IC,IU,R1,R2
      type(mqc_scalar)::sc
      complex(kind=real64),dimension(4)::cs
      character(len=16)::expected
      logical::realData
      integer::k
      cs = [(2.0,0.0),(0.0,1.5),(1.0,2.0),(0.0,0.0)]
      IU = Int
      IC = Int
      call IC%setSymm(sym)
      realData = .not.mqc_matrix_havecomplex(Int%getBlock('full',nosym=.true.))
      R1 = dagger(IC)
      R2 = dagger(IU)
      call t13check(label//': dagger',R1,R2,sym)
      R1 = transpose(IC)
      R2 = transpose(IU)
      call t13check(label//': transpose',R1,R2,sym)
      R1 = conjg(IC)
      R2 = conjg(IU)
      call t13check(label//': conjg',R1,R2,sym)
      do k = 1, 4
        sc = cs(k)
        expected = scaledFlag(sym,cs(k),realData)
        R1 = sc*IC
        R2 = sc*IU
        call t13check(label//': c*A, c='//cstr(cs(k)),R1,R2,expected)
        R1 = IC*sc
        R2 = IU*sc
        call t13check(label//': A*c, c='//cstr(cs(k)),R1,R2,expected)
      end do
      end subroutine t13
!
      subroutine t13check(label,R1,R2,expected)
      character(len=*),intent(in)::label,expected
      type(mqc_scf_integral),intent(in)::R1,R2
      type(mqc_matrix)::F
      integer::k
      logical::ok
      F = R2%getBlock('full',nosym=.true.)
      call check(label//': values',maxdiff(R1%getBlock('full',nosym=.true.),F).lt.tol2)
      if(expected.ne.'') call check(label//': flag -> '//trim(expected),R1%getSymm().eq.expected)
      ok = .true.
      do k = 1, 5
        if(.not.(mqc_integral_test_symmetric(R1,trim(opts(k)),.true.).eqv. &
          mqc_matrix_test_symmetric(F,trim(opts(k)),.true.))) ok = .false.
      end do
      call check(label//': flag consistent with full test',ok)
      end subroutine t13check
!
      function scaledFlag(sym,c,realData) result(flag)
      character(len=*),intent(in)::sym
      complex(kind=real64),intent(in)::c
      logical,intent(in)::realData
      character(len=16)::flag
      logical::isReal,isImag
      isReal = aimag(c).eq.0.0
      isImag = real(c).eq.0.0.and.aimag(c).ne.0.0
      flag = sym
      if(c.eq.(0.0,0.0)) then
        if(sym.ne.'asymmetric') flag = 'fullsymmetric'
        if(sym.eq.'asymmetric') flag = 'fullsymmetric'
        return
      endIf
      if(isReal) return
      select case (sym)
      case('hermitian')
        if(isImag) then
          flag = 'antihermitian'
        elseIf(realData) then
          flag = 'symmetric'
        else
          flag = 'unknown'
        endIf
      case('antihermitian')
        if(isImag) then
          flag = 'hermitian'
        elseIf(realData) then
          flag = 'antisymmetric'
        else
          flag = 'unknown'
        endIf
      case('realdiag')
        if(isImag) then
          flag = 'imagdiag'
        else
          flag = 'complexdiag'
        endIf
      case('imagdiag')
        if(isImag) then
          flag = 'realdiag'
        else
          flag = 'complexdiag'
        endIf
      end select
      end function scaledFlag
!
      function cstr(c) result(str)
      complex(kind=real64),intent(in)::c
      character(len=16)::str
      write(str,'(a,f3.1,a,f3.1,a)') '(',real(c),',',aimag(c),')'
      end function cstr
!
      subroutine t14(i,j)
      integer,intent(in)::i,j
      type(mqc_scf_integral)::R
      type(mqc_matrix)::Fi,Fj,Fr
      character(len=64)::lab
      Fi = OPS(i)%getBlock('full',nosym=.true.)
      Fj = OPS(j)%getBlock('full',nosym=.true.)
      lab = 'T14 '//trim(NAM(i))//' , '//trim(NAM(j))
      R = OPS(i) + OPS(j)
      Fr = Fi + Fj
      call t14check(trim(lab)//': A+B',R,Fr,expect14('sum',i,j,Fr))
      R = OPS(i) - OPS(j)
      Fr = Fi - Fj
      call t14check(trim(lab)//': A-B',R,Fr,expect14('sum',i,j,Fr))
      R = OPS(i) .ewp. OPS(j)
      Fr = MQC_ElementMatrixProduct(Fi,Fj)
      call t14check(trim(lab)//': A.ewp.B',R,Fr,expect14('ewp',i,j,Fr))
      end subroutine t14
!
      subroutine t14h(label,IntA,IntB,expected)
      character(len=*),intent(in)::label,expected
      type(mqc_scf_integral),intent(in)::IntA,IntB
      type(mqc_scf_integral)::R1,R2
      R1 = HermMatMul(IntA,IntB)
      R2 = MatMul(Dagger(IntB),MatMul(IntA,IntB))
      call check(label//': values against MatMul', &
        maxdiff(R1%getBlock('full',nosym=.true.),R2%getBlock('full',nosym=.true.)).lt.tol2)
!     The MatMul reference is (anti-)Hermitian only to round-off, so the flag is
!     checked against the full matrix of the result itself.
      call t14check(label,R1,R1%getBlock('full',nosym=.true.),expected)
      end subroutine t14h
!
      subroutine t14check(label,R,F,expected)
      character(len=*),intent(in)::label,expected
      type(mqc_scf_integral),intent(in)::R
      type(mqc_matrix),intent(in)::F
      integer::k
      logical::ok
      call check(label//': values',maxdiff(R%getBlock('full',nosym=.true.),F).lt.tol2)
      call check(label//': flag -> '//trim(expected),R%getSymm().eq.expected)
      ok = .true.
      do k = 1, 5
        if(.not.(mqc_integral_test_symmetric(R,trim(opts(k)),.true.).eqv. &
          mqc_matrix_test_symmetric(F,trim(opts(k)),.true.))) ok = .false.
      end do
      call check(label//': flag consistent with full test',ok)
      end subroutine t14check
!
!     Expected flag of A+B, A-B ('sum') or A .ewp. B ('ewp') from the symmetries
!     of the operands (symmetric, antisymmetric, hermitian, antihermitian,
!     diagonal); F is the full result.
      function expect14(op,i,j,F) result(flag)
      character(len=*),intent(in)::op
      integer,intent(in)::i,j
      type(mqc_matrix),intent(in)::F
      character(len=16)::flag
      logical,dimension(5)::si,sj,so
      type(mqc_scalar)::sc
      complex(kind=real64)::c
      logical::realOK,imagOK
      integer(kind=int64)::k
      si = symmetriesOf(FLG(i),REALD(i))
      sj = symmetriesOf(FLG(j),REALD(j))
      if(op.eq.'ewp') then
        so(1) = (si(1).and.sj(1)).or.(si(2).and.sj(2))
        so(2) = (si(1).and.sj(2)).or.(si(2).and.sj(1))
        so(3) = (si(3).and.sj(3)).or.(si(4).and.sj(4))
        so(4) = (si(3).and.sj(4)).or.(si(4).and.sj(3))
        so(5) = si(5).or.sj(5)
      else
        so = si.and.sj
      endIf
      sc = F%norm('M')
      if(abs(sc%cval()).lt.tol) then
        flag = 'fullsymmetric'
      elseIf(so(5)) then
        realOK = .true.
        imagOK = .true.
        do k = 1, size(F,1)
          sc = F%at(k,k)
          c = sc%cval()
          if(abs(aimag(c)).gt.tol) realOK = .false.
          if(abs(real(c)).gt.tol) imagOK = .false.
        end do
        if(realOK) then
          flag = 'realdiag'
        elseIf(imagOK) then
          flag = 'imagdiag'
        else
          flag = 'complexdiag'
        endIf
      elseIf(so(3)) then
        flag = 'hermitian'
      elseIf(so(4)) then
        flag = 'antihermitian'
      elseIf(so(1)) then
        flag = 'symmetric'
      elseIf(so(2)) then
        flag = 'antisymmetric'
      elseIf(op.eq.'sum'.and.((FLG(i).eq.'asymmetric'.and.FLG(j).eq.'fullsymmetric').or. &
        (FLG(i).eq.'fullsymmetric'.and.FLG(j).eq.'asymmetric'))) then
        flag = 'asymmetric'
      else
        flag = 'unknown'
      endIf
      end function expect14
!
      function symmetriesOf(sym,realData) result(s)
      character(len=*),intent(in)::sym
      logical,intent(in)::realData
      logical,dimension(5)::s
      select case (sym)
      case('hermitian')
        s = [realData,.false.,.true.,.false.,.false.]
      case('antihermitian')
        s = [.false.,realData,.false.,.true.,.false.]
      case('symmetric')
        s = [.true.,.false.,realData,.false.,.false.]
      case('antisymmetric')
        s = [.false.,.true.,.false.,realData,.false.]
      case('realdiag')
        s = [.true.,.false.,.true.,.false.,.true.]
      case('imagdiag')
        s = [.true.,.false.,.false.,.true.,.true.]
      case('complexdiag')
        s = [.true.,.false.,.false.,.false.,.true.]
      case('fullsymmetric')
        s = .true.
      case default
        s = .false.
      end select
      end function symmetriesOf
!
      function isZeroMat(Mat) result(isZero)
      type(mqc_matrix),intent(in)::Mat
      logical::isZero
      type(mqc_scalar)::sc
      isZero = .true.
      if(MQC_Matrix_isAllocated(Mat)) then
        sc = Mat%norm('M')
        isZero = abs(sc%cval()).lt.tol
      endIf
      end function isZeroMat
!
!     T16 helpers.
      function hermPD(nh,seed) result(Mat)
      integer(kind=int64),intent(in)::nh
      integer,intent(in)::seed
      type(mqc_matrix)::Mat,Xh
      Xh = cfill(nh,nh,seed)
      Mat = matmul(dagger(Xh),Xh) + eyeC(nh)
      end function hermPD
!
      function eyeC(ne) result(Mat)
      integer(kind=int64),intent(in)::ne
      type(mqc_matrix)::Mat
      integer(kind=int64)::i
      call Mat%init(ne,ne,(0.0,0.0))
      do i = 1, ne
        call Mat%put((1.0,0.0),i,i)
      end do
      end function eyeC
!
      function eyeR(ne) result(Mat)
      integer(kind=int64),intent(in)::ne
      type(mqc_matrix)::Mat
      integer(kind=int64)::i
      call Mat%init(ne,ne,0.0)
      do i = 1, ne
        call Mat%put(1.0,i,i)
      end do
      end function eyeR
!
      function fullOf(Mat) result(MatOut)
      type(mqc_matrix),intent(in)::Mat
      type(mqc_matrix)::MatOut
      MatOut = Mat
      if(MatOut%getSymm().eq.'diagonal') then
        call mqc_matrix_diag2full(MatOut)
      elseIf(MatOut%getSymm().ne.'unknown') then
        call mqc_matrix_symm2full(MatOut)
      endIf
      end function fullOf
!
      subroutine tinv(label,Mat,storeAfter)
      character(len=*),intent(in)::label,storeAfter
      type(mqc_matrix),intent(in)::Mat
      type(mqc_matrix)::Mi,Mc
      Mc = Mat
      Mi = Mc%inv()
      call check(label//': A^-1 A = I',maxdiff(matmul(fullOf(Mi),fullOf(Mat)),eyeC(size(Mat,1))).lt.tol2)
      call check(label//': storage -> '//storeAfter,Mi%getSymm().eq.storeAfter)
      end subroutine tinv
!
      subroutine tpowh(label,Mat,storeAfter)
      character(len=*),intent(in)::label,storeAfter
      type(mqc_matrix),intent(in)::Mat
      type(mqc_matrix)::Mh,Mf
      Mh = Mat
      call Mh%power(0.5)
      Mf = fullOf(Mh)
      call check(label//': (A^1/2)^2 = A',maxdiff(matmul(Mf,Mf),fullOf(Mat)).lt.tol2)
      call check(label//': storage -> '//storeAfter,Mh%getSymm().eq.storeAfter)
      end subroutine tpowh
!
      subroutine tpowk(label,Mat,k,storeAfter)
      character(len=*),intent(in)::label,storeAfter
      type(mqc_matrix),intent(in)::Mat
      integer,intent(in)::k
      type(mqc_matrix)::Mk,Mr,Mf
      integer::j
      Mk = Mat
      call Mk%power(k)
      Mf = fullOf(Mat)
      if(k.gt.0) then
        Mr = Mf
        do j = 2, k
          Mr = matmul(Mr,Mf)
        end do
        call check(label//': value against products',maxdiff(fullOf(Mk),Mr).lt.tol2)
      else
        Mr = matmul(fullOf(Mk),Mf)
        call check(label//': A^-1 A = I',maxdiff(Mr,eyeC(size(Mat,1))).lt.tol2)
      endIf
      call check(label//': storage -> '//storeAfter,Mk%getSymm().eq.storeAfter)
      end subroutine tpowk
!
      subroutine tpinv(label,Mat,rank,storeAfter)
      character(len=*),intent(in)::label,storeAfter
      type(mqc_matrix),intent(in)::Mat
      integer,intent(in)::rank
      type(mqc_matrix)::Mp,Mc,Af,Pf,AP,PA
      integer(kind=int64)::sd
      Mc = Mat
      Mp = Mc%pinv(spaceDim=sd)
      Af = fullOf(Mat)
      Pf = fullOf(Mp)
      AP = matmul(Af,Pf)
      PA = matmul(Pf,Af)
      call check(label//': A P A = A',maxdiff(matmul(AP,Af),Af).lt.tol2)
      call check(label//': P A P = P',maxdiff(matmul(PA,Pf),Pf).lt.tol2)
      call check(label//': A P and P A hermitian', &
        maxdiff(AP,dagger(AP)).lt.tol2.and.maxdiff(PA,dagger(PA)).lt.tol2)
      call check(label//': spaceDim',sd.eq.rank)
      call check(label//': storage -> '//storeAfter,Mp%getSymm().eq.storeAfter)
      end subroutine tpinv
!
!     Integral inverse, power 1/2 and pinv against the matrix routines on the full
!     matrix (values); expected flags of the three results, checked against the full
!     matrix of each result (the references are symmetric only to round-off).
      subroutine tint9(label,Int,invFlag,sqrtFlag,pinvFlag)
      character(len=*),intent(in)::label,invFlag,sqrtFlag,pinvFlag
      type(mqc_scf_integral),intent(in)::Int
      type(mqc_scf_integral)::R
      type(mqc_matrix)::F,Fr
      F = Int%getBlock('full',nosym=.true.)
      R = Int%inv()
      Fr = F%inv()
      call check(label//': inv: values',maxdiff(R%getBlock('full',nosym=.true.),fullOf(Fr)).lt.tol2)
      call t15(label//': inv',R,R%getBlock('full',nosym=.true.),invFlag)
      R = Int
      call R%power(0.5)
      Fr = F
      call Fr%power(0.5)
      call check(label//': power 1/2: values',maxdiff(R%getBlock('full',nosym=.true.),fullOf(Fr)).lt.tol2)
      call t15(label//': power 1/2',R,R%getBlock('full',nosym=.true.),sqrtFlag)
      R = Int%pinv()
      Fr = F%pinv()
      call check(label//': pinv: values',maxdiff(R%getBlock('full',nosym=.true.),fullOf(Fr)).lt.tol2)
      call t15(label//': pinv',R,R%getBlock('full',nosym=.true.),pinvFlag)
      end subroutine tint9
!
!     T17: swap of every T14 operand against the explicit block swap.
      subroutine t17swap(i)
      integer,intent(in)::i
      type(mqc_scf_integral)::R,Ref
      type(mqc_matrix)::Ba1,Bb1,Bab1,Bba1
      character(len=16)::expected
      R = mqc_integral_swap_offDiagBlocks(OPS(i))
      if(OPS(i)%type().eq.'general') then
        Ba1 = OPS(i)%getBlock('alpha',nosym=.true.)
        Bb1 = OPS(i)%getBlock('beta',nosym=.true.)
        Bab1 = OPS(i)%getBlock('alpha-beta')
        Bba1 = OPS(i)%getBlock('beta-alpha')
        call mqc_integral_allocate(Ref,'r','general',Ba1,Bb1,Bba1,Bab1)
        expected = FLG(i)
      else
        Ref = OPS(i)
        expected = FLG(i)
      endIf
      call t15('T17 swap '//trim(NAM(i)),R,Ref%getBlock('full',nosym=.true.),expected)
      end subroutine t17swap
!
!     T17: integral pad against the pad of the unflagged integral.
      subroutine t17pad(label,Int,ax1,ax2,block,expected)
      character(len=*),intent(in)::label,block,expected
      type(mqc_scf_integral),intent(in)::Int
      integer,intent(in)::ax1,ax2
      type(mqc_scf_integral)::R,IU,RU
      integer(kind=int64)::a1,a2
      a1 = ax1
      a2 = ax2
      R = mqc_integral_pad(Int,a1,a2,block)
      IU = Int
      call IU%setSymm('unknown')
      RU = mqc_integral_pad(IU,a1,a2,block)
      call t15(label,R,RU%getBlock('full',nosym=.true.),expected)
      end subroutine t17pad
!
!     Eigenvalue vector for T15 (real data when isReal).
      function evec(c,isReal) result(V)
      complex(kind=real64),dimension(:),intent(in)::c
      logical,intent(in)::isReal
      type(mqc_vector)::V
      integer(kind=int64)::i,nv
      nv = size(c)
      if(isReal) then
        call V%init(nv,0.0)
        do i = 1, nv
          call V%put(real(c(i)),i)
        end do
      else
        call V%init(nv,(0.0,0.0))
        do i = 1, nv
          call V%put(c(i),i)
        end do
      endIf
      end function evec
!
!     Block-diagonal reference for T15: alpha nA1 x nA2 and beta nB1 x nB2 filled
!     with c (on the diagonal only when diagOnly).
      function ref15(nA1,nA2,nB1,nB2,c,diagOnly) result(F)
      integer,intent(in)::nA1,nA2,nB1,nB2
      complex(kind=real64),intent(in)::c
      logical,intent(in)::diagOnly
      type(mqc_matrix)::F
      integer(kind=int64)::i,j,r,s
      r = nA1+nB1
      s = nA2+nB2
      call F%init(r,s,(0.0,0.0))
      do i = 1, r
        do j = 1, s
          if((i.le.nA1).neqv.(j.le.nA2)) cycle
          if(diagOnly) then
            if((i.le.nA1.and.i.eq.j).or.(i.gt.nA1.and.i-nA1.eq.j-nA2)) call F%put(c,i,j)
          else
            call F%put(c,i,j)
          endIf
        end do
      end do
      end function ref15
!
      subroutine t15(label,Int,F,flag)
      character(len=*),intent(in)::label,flag
      type(mqc_scf_integral),intent(in)::Int
      type(mqc_matrix),intent(in)::F
      integer::k
      logical::ok
      call check(label//': values',maxdiff(Int%getBlock('full',nosym=.true.),F).lt.tol)
      call check(label//': flag -> '//flag,Int%getSymm().eq.flag)
      ok = .true.
      do k = 1, 5
        if(.not.(mqc_integral_test_symmetric(Int,trim(opts(k)),.true.).eqv. &
          mqc_matrix_test_symmetric(F,trim(opts(k)),.true.))) ok = .false.
      end do
      call check(label//': flag consistent with full test',ok)
      end subroutine t15
!
      subroutine timpl(label,Ma,Mb,sym,Mab,Mba,Rab,Rba,flagAfter)
      character(len=*),intent(in)::label,sym,flagAfter
      type(mqc_matrix),intent(in)::Ma,Mb,Mab,Mba,Rab,Rba
      type(mqc_scf_integral)::IC,IR
      if(MQC_Matrix_isAllocated(Mab)) then
        call mqc_integral_allocate(IC,'c','general',Ma,Mb,alphaBeta=Mab,symmetry=sym)
      else
        call mqc_integral_allocate(IC,'c','general',Ma,Mb,betaAlpha=Mba,symmetry=sym)
      endIf
      call mqc_integral_allocate(IR,'r','general',Ma,Mb,Rab,Rba)
      call check(label//': values (implied block)', &
        maxdiff(IC%getBlock('full',nosym=.true.),IR%getBlock('full')).lt.tol)
      call check(label//': flag -> '//flagAfter,IC%getSymm().eq.flagAfter)
      if(flagAfter.ne.'unknown') call check(label//': alpha packed', &
        mqc_matrix_get_symmetry(IC%getBlock('alpha')).ne.'unknown')
      end subroutine timpl
!
      subroutine timp(label,Int,Mat,typeAfter,flagAfter)
      character(len=*),intent(in)::label,typeAfter,flagAfter
      type(mqc_scf_integral),intent(in)::Int
      type(mqc_matrix),intent(in)::Mat
      type(mqc_scf_integral)::IC
      IC = Int
      call IC%imp(Mat)
      call check(label//': values',maxdiff(IC%getBlock('full'),Mat).lt.tol)
      call check(label//': type -> '//typeAfter,IC%type().eq.typeAfter)
      if(flagAfter.ne.'') call check(label//': flag -> '//flagAfter,IC%getSymm().eq.flagAfter)
      end subroutine timp
!
      end program unitTest_symmetry
