      program unitTest_NEF_twoeris_r4_transform
!
!     unitTest_NEF_twoeris_r4_transform: checks that
!     mqc_twoeris_nef_r4_transform (the flat-AO-ERI, four-independent-
!     NEF_Matrix-coefficient transform that returns an MQC_NEF_R4,
!     built directly on TwoERI_Trans_Full) agrees with an INDEPENDENT
!     reference built entirely from already-established NEF primitives:
!     outer2ERI/outer2R4 (to build the test ERI and the per-axis-pair
!     coefficient outer products) and PartialContraction2R4 (to contract
!     the AO legs away -- the ERI-first overload for the flat test ERI,
!     the NEF_R4-vs-NEF_R4 overload for the second stage), with Conjg
!     applied to the "bra" (P,R-leg) coefficients to match
!     mqc_twoeris_nef_r4_transform's own conjugation convention (via
!     TwoERI_Trans_Full's Dagger(MO_Coeff_1/3)) on those axes.
!
!     Unlike unitTest_NEF_r4_full_transform (whose NEF_Matrix coefficients
!     have TWO row-spaces, exercising MQC_NEF_R4%transform's fully general
!     multi-row-space case), mqc_twoeris_nef_r4_transform's own design
!     assumes a SINGLE row-space per coefficient spanning the whole flat
!     AO basis (NumSpace_row=1) -- the same shape IMat%NEF_natOrb always
!     has (InterMatrixMod.f03:6479-6488) -- since it extracts each
!     coefficient's per-column-space block via blockIndex([1,1],[2,p]),
!     hardcoding row-space index 1. So this fixture uses ONE row-space
!     (the full 2-orbital AO basis) and TWO column-spaces (1 orbital
!     each), rather than copying the sibling test's two-row-space
!     fixture verbatim.
!
!     Fixture: flatA/flatB (dense, complex, genuinely 'space'-type --
!     see the NOTE below on avoiding mqc_integral_put's type-promotion
!     quirk) build the test ERI via outer2ERI. The four transform
!     coefficients are a genuine MIX of types -- NEFCoeff_1, NEFCoeff_3
!     'space'-type; NEFCoeff_2, NEFCoeff_4 'spin'-type -- each fully
!     dense (including off-diagonal space blocks, to exercise genuine
!     space-mixing) with distinct complex values, driving the resolved
!     output type to 'spin'.
!
!     -M. M. F. Moraes, 2026.
!
      use MQC_EST
      use MQC_NEF
      use iso_fortran_env, only: int64, real64
!
      implicit none
      integer::nfail
      real(kind=real64),parameter::tol=1.0d-8
!
      type(mqc_scf_eigenvalues),dimension(:),allocatable::rowMap,colMap
      type(mqc_scf_integral)::flatA,flatB,flatC1,flatC2,flatC3,flatC4
      type(mqc_twoERIs)::testERI
      type(mqc_nef_matrix)::NEFCoeff_1,NEFCoeff_2,NEFCoeff_3,NEFCoeff_4
      type(mqc_nef_matrix)::coeff1Conjg,coeff3Conjg
      type(mqc_nef_r4)::transformResult,stage1,refResult,outerPQ,outerRS
      type(mqc_scalar)::sc
      complex(kind=real64)::cGot,cExpected
      integer(kind=int64)::P,Q,R,S
      integer(kind=int64),dimension(4)::set1,set2,set3,set4,spinCombo
      character(len=80)::testLabel
      character(len=5)::spinName
      integer::spinIdx
!
 1000 format(1x,'PASS  ',a)
 1010 format(1x,'FAIL  ',a,'  got=(',f14.8,',',f14.8,')  expected=(',f14.8,',',f14.8,')')
 9000 format(/,1x,'Tests completed: ',i0,' failure(s)')
!
      nfail = 0
!
!--------------------------------------------------------------------
!     Row map: ONE space spanning the full 2-orbital AO basis, alpha
!     and beta both covering both indices in parallel -- mirrors
!     fullBasisMap in inter_matrix_compute_natural_orbitals.
!--------------------------------------------------------------------
!
      allocate(rowMap(1))
      call rowMap(1)%init(0,0)
      call rowMap(1)%push(1,'alpha'); call rowMap(1)%push(2,'alpha')
      call rowMap(1)%push(1,'beta');  call rowMap(1)%push(2,'beta')
!
!--------------------------------------------------------------------
!     Column map: TWO spaces, 1 orbital each -- mirrors a 2-space
!     natural-orbital-style partition (e.g. occupied/virtual) of the
!     same 2-orbital AO basis.
!--------------------------------------------------------------------
!
      allocate(colMap(2))
      call colMap(1)%init(0,0)
      call colMap(1)%push(1,'alpha'); call colMap(1)%push(1,'beta')
      call colMap(2)%init(0,0)
      call colMap(2)%push(2,'alpha'); call colMap(2)%push(2,'beta')
!
!--------------------------------------------------------------------
!     flatA, flatB: dense, complex, used only to build the flat test
!     ERI via outer2ERI.
!
!     NOTE: puts below deliberately omit the spinBlockIn label --
!     mqc_integral_put promotes a 'space'-type integral to 'spin' the
!     moment it is called WITH an explicit 'alpha' label
!     (mqc_est.F03:23365-23368), which would make testERI 'spin'-typed
!     and fail TwoERI_Trans_Full's 'regular'/'space'-only precondition.
!--------------------------------------------------------------------
!
      call flatA%init(2,2)
      sc = cmplx(0.4d0,0.2d0,kind=real64); call flatA%put(sc,1,1)
      sc = cmplx(0.1d0,-0.3d0,kind=real64); call flatA%put(sc,1,2)
      sc = cmplx(-0.2d0,0.5d0,kind=real64); call flatA%put(sc,2,1)
      sc = cmplx(0.3d0,0.1d0,kind=real64); call flatA%put(sc,2,2)
!
      call flatB%init(2,2)
      sc = cmplx(-0.1d0,0.4d0,kind=real64); call flatB%put(sc,1,1)
      sc = cmplx(0.2d0,0.2d0,kind=real64); call flatB%put(sc,1,2)
      sc = cmplx(0.3d0,-0.1d0,kind=real64); call flatB%put(sc,2,1)
      sc = cmplx(-0.4d0,0.3d0,kind=real64); call flatB%put(sc,2,2)
!
      testERI = flatA%outer2ERI(flatB)
!
!--------------------------------------------------------------------
!     Four coefficients, each (full 2-orbital AO row) x (2 one-orbital
!     column spaces): NEFCoeff_1, NEFCoeff_3 'space'-type (plain puts,
!     no label, per the NOTE above); NEFCoeff_2, NEFCoeff_4 'spin'-type
!     (explicit alpha+beta puts are intentional here -- distinct values
!     per spin is exactly what makes them genuinely 'spin'-type), each
!     fully dense and complex.
!--------------------------------------------------------------------
!
      call flatC1%init(2,2)
      sc = cmplx(0.5d0,0.1d0,kind=real64); call flatC1%put(sc,1,1)
      sc = cmplx(-0.2d0,0.3d0,kind=real64); call flatC1%put(sc,1,2)
      sc = cmplx(0.1d0,-0.4d0,kind=real64); call flatC1%put(sc,2,1)
      sc = cmplx(0.3d0,0.2d0,kind=real64); call flatC1%put(sc,2,2)
      call NEFCoeff_1%init(rowMap,initial_map_matrix_col=colMap)
      call NEFCoeff_1%imp(flatC1)
!
      call flatC2%init(2,2)
      sc = cmplx(0.2d0,-0.1d0,kind=real64); call flatC2%put(sc,1,1,'alpha')
      sc = cmplx(0.3d0,0.2d0,kind=real64); call flatC2%put(sc,1,2,'alpha')
      sc = cmplx(-0.1d0,0.1d0,kind=real64); call flatC2%put(sc,2,1,'alpha')
      sc = cmplx(0.4d0,-0.2d0,kind=real64); call flatC2%put(sc,2,2,'alpha')
      sc = cmplx(-0.3d0,0.2d0,kind=real64); call flatC2%put(sc,1,1,'beta')
      sc = cmplx(0.1d0,0.3d0,kind=real64); call flatC2%put(sc,1,2,'beta')
      sc = cmplx(0.2d0,-0.3d0,kind=real64); call flatC2%put(sc,2,1,'beta')
      sc = cmplx(-0.4d0,0.1d0,kind=real64); call flatC2%put(sc,2,2,'beta')
      call NEFCoeff_2%init(rowMap,initial_map_matrix_col=colMap)
      call NEFCoeff_2%imp(flatC2)
!
      call flatC3%init(2,2)
      sc = cmplx(-0.3d0,0.2d0,kind=real64); call flatC3%put(sc,1,1)
      sc = cmplx(0.4d0,-0.1d0,kind=real64); call flatC3%put(sc,1,2)
      sc = cmplx(0.2d0,0.3d0,kind=real64); call flatC3%put(sc,2,1)
      sc = cmplx(-0.1d0,-0.2d0,kind=real64); call flatC3%put(sc,2,2)
      call NEFCoeff_3%init(rowMap,initial_map_matrix_col=colMap)
      call NEFCoeff_3%imp(flatC3)
!
      call flatC4%init(2,2)
      sc = cmplx(0.1d0,0.4d0,kind=real64); call flatC4%put(sc,1,1,'alpha')
      sc = cmplx(-0.2d0,-0.1d0,kind=real64); call flatC4%put(sc,1,2,'alpha')
      sc = cmplx(0.3d0,0.2d0,kind=real64); call flatC4%put(sc,2,1,'alpha')
      sc = cmplx(-0.4d0,0.3d0,kind=real64); call flatC4%put(sc,2,2,'alpha')
      sc = cmplx(0.2d0,-0.2d0,kind=real64); call flatC4%put(sc,1,1,'beta')
      sc = cmplx(0.1d0,0.1d0,kind=real64); call flatC4%put(sc,1,2,'beta')
      sc = cmplx(-0.3d0,0.4d0,kind=real64); call flatC4%put(sc,2,1,'beta')
      sc = cmplx(0.2d0,0.3d0,kind=real64); call flatC4%put(sc,2,2,'beta')
      call NEFCoeff_4%init(rowMap,initial_map_matrix_col=colMap)
      call NEFCoeff_4%imp(flatC4)
!
!--------------------------------------------------------------------
!     Method 1: the function under test.
!--------------------------------------------------------------------
!
      transformResult = mqc_twoeris_nef_r4_transform(testERI,NEFCoeff_1,&
           NEFCoeff_2,NEFCoeff_3,NEFCoeff_4)
!
!--------------------------------------------------------------------
!     Method 2: independent reference built from Conjg + outer2ERI/
!     outer2R4 + PartialContraction2R4. PQ-leg contracted first (ERI-
!     first overload, since testERI is still flat at this stage), then
!     RS-leg (NEF_R4-vs-NEF_R4 overload), landing directly in
!     (P,Q,R,S) axis order.
!--------------------------------------------------------------------
!
      coeff1Conjg = Conjg(NEFCoeff_1)
      coeff3Conjg = Conjg(NEFCoeff_3)
      outerPQ = coeff1Conjg%outer2R4(NEFCoeff_2)
      outerRS = coeff3Conjg%outer2R4(NEFCoeff_4)
      stage1 = PartialContraction2R4(testERI,outerPQ,[1,2,3,4],[1,3,2,4])
      refResult = PartialContraction2R4(stage1,outerRS,[1,2,3,4],[1,3,2,4])
!
!--------------------------------------------------------------------
!     Compare over all space-quadruples (2 column-spaces per leg) and
!     all 4 spin combinations (alpha, beta, alphaBeta, betaAlpha).
!--------------------------------------------------------------------
!
      do P = 1, 2
        do Q = 1, 2
          do R = 1, 2
            do S = 1, 2
              do spinIdx = 1, 4
                select case(spinIdx)
                case(1)
                  spinCombo = [1,1,1,1]; spinName = 'alpha'
                case(2)
                  spinCombo = [-1,-1,-1,-1]; spinName = 'beta'
                case(3)
                  spinCombo = [1,1,-1,-1]; spinName = 'aBeta'
                case default
                  spinCombo = [-1,-1,1,1]; spinName = 'bAlph'
                end select
                set1 = [1,P,1,spinCombo(1)]
                set2 = [2,Q,1,spinCombo(2)]
                set3 = [3,R,1,spinCombo(3)]
                set4 = [4,S,1,spinCombo(4)]
                sc = transformResult%at(set1,set2,set3,set4); cGot = sc
                sc = refResult%at(set1,set2,set3,set4); cExpected = sc
                write(testLabel,'(a,i0,a,i0,a,i0,a,i0,a,a)') &
                  'twoeris_transform vs outer+partialContraction: space(',P,',',Q,',',R,',',S,') ',&
                  trim(spinName)
                if(abs(cGot-cExpected).le.tol) then
                  write(*,1000) trim(testLabel)
                else
                  write(*,1010) trim(testLabel),real(cGot),aimag(cGot),real(cExpected),aimag(cExpected)
                  nfail = nfail+1
                end if
              end do
            end do
          end do
        end do
      end do
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1
!
      end program unitTest_NEF_twoeris_r4_transform
