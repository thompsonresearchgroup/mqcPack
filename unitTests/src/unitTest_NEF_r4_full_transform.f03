      program unitTest_NEF_r4_full_transform
!
!     unitTest_NEF_r4_full_transform: checks that MQC_NEF_R4%transform
!     (the four-independent-NEF_Matrix-coefficient, space-mixing AO->MO
!     two-electron integral transform) agrees with an INDEPENDENT
!     reference built entirely from already-established NEF primitives:
!     outer2R4 (to build both the test ERI itself and the per-axis-pair
!     coefficient outer products) and PartialContraction2R4 (to contract
!     the AO legs away), with Conjg applied to the "bra" (P,R-leg)
!     coefficients to match %transform's own conjugation convention on
!     those axes (PartialContraction2R4 itself does a plain, non-
!     conjugating sum).
!
!     Fixture: 2 spaces, 1 alpha + 1 beta orbital each. The ERI itself is
!     NEFR4A = matA%outer2R4(matB), matA/matB being dense, complex
!     'space'-type NEF_Matrix (a deliberately simple, separable test
!     tensor -- any nonzero 4-index object works as a transform input).
!     The four transform coefficients are a genuine MIX of types --
!     NEFCoeff_1, NEFCoeff_3 'space'-type; NEFCoeff_2, NEFCoeff_4
!     'spin'-type -- each fully dense (including off-diagonal space
!     blocks, to exercise genuine space-mixing) with distinct complex
!     values, driving the resolved output type to 'spin'.
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
      type(mqc_scf_eigenvalues),dimension(:),allocatable::rowMap
      type(mqc_scf_integral)::flatA,flatB,flatC1,flatC2,flatC3,flatC4
      type(mqc_nef_matrix)::matA,matB,NEFCoeff_1,NEFCoeff_2,NEFCoeff_3,NEFCoeff_4
      type(mqc_nef_matrix)::coeff1Conjg,coeff3Conjg
      type(mqc_nef_r4)::NEFR4A,transformResult,stage1,refResult,outerPQ,outerRS
      type(mqc_scalar)::sc
      complex(kind=real64)::cGot,cExpected
      integer(kind=int64)::P,Q,R,S,spaceIdx
      integer(kind=int64),dimension(4)::set1,set2,set3,set4,spinCombo
      character(len=64)::testLabel
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
!     2-space row/col map, 1 alpha + 1 beta orbital per space.
!--------------------------------------------------------------------
!
      allocate(rowMap(2))
      do spaceIdx = 1, 2
        call rowMap(spaceIdx)%init(0,0)
        call rowMap(spaceIdx)%push(spaceIdx,'alpha')
        call rowMap(spaceIdx)%push(spaceIdx,'beta')
      end do
!
!--------------------------------------------------------------------
!     matA, matB: dense, complex, 'space'-type (alpha-only puts), used
!     only to build the test ERI via outer2R4.
!--------------------------------------------------------------------
!
      call flatA%init(2,2)
      sc = cmplx(0.4d0,0.2d0,kind=real64); call flatA%put(sc,1,1,'alpha')
      sc = cmplx(0.1d0,-0.3d0,kind=real64); call flatA%put(sc,1,2,'alpha')
      sc = cmplx(-0.2d0,0.5d0,kind=real64); call flatA%put(sc,2,1,'alpha')
      sc = cmplx(0.3d0,0.1d0,kind=real64); call flatA%put(sc,2,2,'alpha')
      call matA%init(rowMap)
      call matA%imp(flatA)
!
      call flatB%init(2,2)
      sc = cmplx(-0.1d0,0.4d0,kind=real64); call flatB%put(sc,1,1,'alpha')
      sc = cmplx(0.2d0,0.2d0,kind=real64); call flatB%put(sc,1,2,'alpha')
      sc = cmplx(0.3d0,-0.1d0,kind=real64); call flatB%put(sc,2,1,'alpha')
      sc = cmplx(-0.4d0,0.3d0,kind=real64); call flatB%put(sc,2,2,'alpha')
      call matB%init(rowMap)
      call matB%imp(flatB)
!
      NEFR4A = matA%outer2R4(matB)
!
!--------------------------------------------------------------------
!     Four coefficients: NEFCoeff_1, NEFCoeff_3 'space'-type; NEFCoeff_2,
!     NEFCoeff_4 'spin'-type, each dense (including off-diagonal space
!     blocks) and complex.
!--------------------------------------------------------------------
!
      call flatC1%init(2,2)
      sc = cmplx(0.5d0,0.1d0,kind=real64); call flatC1%put(sc,1,1,'alpha')
      sc = cmplx(-0.2d0,0.3d0,kind=real64); call flatC1%put(sc,1,2,'alpha')
      sc = cmplx(0.1d0,-0.4d0,kind=real64); call flatC1%put(sc,2,1,'alpha')
      sc = cmplx(0.3d0,0.2d0,kind=real64); call flatC1%put(sc,2,2,'alpha')
      call NEFCoeff_1%init(rowMap)
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
      call NEFCoeff_2%init(rowMap)
      call NEFCoeff_2%imp(flatC2)
!
      call flatC3%init(2,2)
      sc = cmplx(-0.3d0,0.2d0,kind=real64); call flatC3%put(sc,1,1,'alpha')
      sc = cmplx(0.4d0,-0.1d0,kind=real64); call flatC3%put(sc,1,2,'alpha')
      sc = cmplx(0.2d0,0.3d0,kind=real64); call flatC3%put(sc,2,1,'alpha')
      sc = cmplx(-0.1d0,-0.2d0,kind=real64); call flatC3%put(sc,2,2,'alpha')
      call NEFCoeff_3%init(rowMap)
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
      call NEFCoeff_4%init(rowMap)
      call NEFCoeff_4%imp(flatC4)
!
!--------------------------------------------------------------------
!     Method 1: the function under test.
!--------------------------------------------------------------------
!
      transformResult = NEFR4A%transform(NEFCoeff_1,NEFCoeff_2,NEFCoeff_3,NEFCoeff_4)
!
!--------------------------------------------------------------------
!     Method 2: independent reference built from Conjg + outer2R4 +
!     PartialContraction2R4. PQ-leg contracted first, then RS-leg, so the
!     result lands directly in (P,Q,R,S) axis order.
!--------------------------------------------------------------------
!
      coeff1Conjg = Conjg(NEFCoeff_1)
      coeff3Conjg = Conjg(NEFCoeff_3)
      outerPQ = coeff1Conjg%outer2R4(NEFCoeff_2)
      outerRS = coeff3Conjg%outer2R4(NEFCoeff_4)
      stage1 = PartialContraction2R4(NEFR4A,outerPQ,[1,2,3,4],[1,3,2,4])
      refResult = PartialContraction2R4(stage1,outerRS,[1,2,3,4],[1,3,2,4])
!
!--------------------------------------------------------------------
!     Compare over all space-quadruples and all 4 spin combinations
!     (alpha, beta, alphaBeta, betaAlpha).
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
                  'transform vs outer+partialContraction: space(',P,',',Q,',',R,',',S,') ',&
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
      end program unitTest_NEF_r4_full_transform
