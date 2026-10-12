      program unitTest_NEF_empty_space_transform
!
!     Reproduces a suspected mqc_twoeris_nef_r4_transform bug: builds a
!     NEF_Matrix with the SAME structure as IMat%NEF_natOrb (one row
!     space spanning the full AO basis with alpha/beta in parallel;
!     THREE column spaces -- closed/active/virtual -- where the middle
!     (active) space is genuinely empty, 0 alpha + 0 beta orbitals,
!     while closed and virtual each have 2 alpha + 2 beta, matching a
!     minimal-basis H4-like system with no active space), then calls
!     mqc_twoeris_nef_r4_transform with that single coefficient on all
!     4 legs -- exactly mirroring inter_matrix_convert_NatOrb_2e_eri's
!     NEF-branch call -- to see whether the "Column length badly
!     defined in MQC_Matrix_Matrix_Put" crash (seen inside
!     TwoERI_Trans_Full's getBlock('beta') calls) reproduces here, in
!     isolation from the rest of SNOWMAN.
!
!     -M. M. F. Moraes, 2026.
!
      use MQC_EST
      use MQC_NEF
      use iso_fortran_env, only: int64, real64
!
      implicit none
      type(mqc_scf_eigenvalues),dimension(:),allocatable::rowMap,colMap
      type(mqc_scf_integral)::flatA,flatB,flatCoeff
      type(mqc_twoERIs)::testERI
      type(mqc_nef_matrix)::natOrbLike
      type(mqc_nef_r4)::transformResult
      type(mqc_scalar)::sc
      integer(kind=int64)::i
!
!--------------------------------------------------------------------
!     Row map: ONE space spanning the full basis (nBasis=4), alpha and
!     beta both covering all 4 indices in parallel -- mirrors
!     fullBasisMap in inter_matrix_compute_natural_orbitals.
!--------------------------------------------------------------------
!
      allocate(rowMap(1))
      call rowMap(1)%init(0,0)
      do i = 1, 4
        call rowMap(1)%push(i,'alpha')
        call rowMap(1)%push(i,'beta')
      end do
!
!--------------------------------------------------------------------
!     Column map: THREE spaces -- closed (2+2), active (0+0, EMPTY),
!     virtual (2+2) -- mirrors IMat%nat_orb_space for a minimal-basis
!     system with no active space.
!--------------------------------------------------------------------
!
      allocate(colMap(3))
      call colMap(1)%init(0,0)
      call colMap(1)%push(1,'alpha'); call colMap(1)%push(2,'alpha')
      call colMap(1)%push(1,'beta');  call colMap(1)%push(2,'beta')
      call colMap(2)%init(0,0)
      call colMap(3)%init(0,0)
      call colMap(3)%push(3,'alpha'); call colMap(3)%push(4,'alpha')
      call colMap(3)%push(3,'beta');  call colMap(3)%push(4,'beta')
!
!--------------------------------------------------------------------
!     flatCoeff: a dense 4x4 complex 'spin'-type coefficient (distinct
!     alpha/beta values, like a genuine unrestricted NatOrb coefficient
!     matrix), imported into natOrbLike exactly as
!     inter_matrix_compute_natural_orbitals imports total_mo into
!     IMat%NEF_natOrb.
!--------------------------------------------------------------------
!
      call flatCoeff%init(4,4)
      do i = 1, 4
        sc = cmplx(1.0d0,0.0d0,kind=real64)
        call flatCoeff%put(sc,i,i,'alpha')
        sc = cmplx(0.9d0,0.1d0,kind=real64)
        call flatCoeff%put(sc,i,i,'beta')
      end do
!
      call natOrbLike%init(rowMap,initial_map_matrix_col=colMap)
      call natOrbLike%imp(flatCoeff)
!
!--------------------------------------------------------------------
!     Test ERI: flat, 'space'-type, nBasis=4.
!--------------------------------------------------------------------
!
!     NOTE: puts below deliberately omit the spinBlockIn label -- mqc_integral_put
!     promotes a 'space'-type integral to 'spin' the moment it is called WITH an
!     explicit 'alpha' label (mqc_est.F03:23365-23368), which would make testERI
!     'spin'-typed and fail TwoERI_Trans_Full's 'regular'/'space'-only precondition.
      call flatA%init(4,4)
      sc = cmplx(0.4d0,0.0d0,kind=real64); call flatA%put(sc,1,1)
      sc = cmplx(0.1d0,0.0d0,kind=real64); call flatA%put(sc,2,2)
      sc = cmplx(0.2d0,0.0d0,kind=real64); call flatA%put(sc,3,3)
      sc = cmplx(0.3d0,0.0d0,kind=real64); call flatA%put(sc,4,4)
      call flatB%init(4,4)
      sc = cmplx(-0.1d0,0.0d0,kind=real64); call flatB%put(sc,1,1)
      sc = cmplx(0.2d0,0.0d0,kind=real64); call flatB%put(sc,2,2)
      sc = cmplx(0.3d0,0.0d0,kind=real64); call flatB%put(sc,3,3)
      sc = cmplx(-0.4d0,0.0d0,kind=real64); call flatB%put(sc,4,4)
      testERI = flatA%outer2ERI(flatB)
!
!--------------------------------------------------------------------
!     Reproduce: single coefficient on all 4 legs, exactly matching
!     inter_matrix_convert_NatOrb_2e_eri's NEF-branch call.
!--------------------------------------------------------------------
!
      write(*,*) 'Calling mqc_twoeris_nef_r4_transform with an empty '//&
        'middle natural-orbital space...'
      transformResult = mqc_twoeris_nef_r4_transform(testERI,natOrbLike)
      write(*,*) 'PASS: mqc_twoeris_nef_r4_transform completed without error'
!
      end program unitTest_NEF_empty_space_transform
