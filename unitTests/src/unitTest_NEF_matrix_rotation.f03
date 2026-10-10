      program unitTest_NEF_matrix_rotation
!
!     unitTest_NEF_matrix_rotation: checks that MQC_NEF_Matrix%getRotation
!     (the cyclic-Jacobi inter-space block-diagonalizer) actually satisfies
!     its defining property: given a Hermitian NEF_Matrix W and
!     R = W%getRotation(), Dagger(R)*W*R has every inter-space off-diagonal
!     block equal to zero, and R itself is unitary (Dagger(R)*R = identity).
!
!     Fixture: 3 spaces, 1 alpha + 1 beta orbital each, Hermitian
!     ('spin'-type, built via independent 'alpha'/'beta' puts), with
!     genuinely complex coupling on BOTH pairs (1,2) and (1,3) (sharing
!     space 1, so the sweep must actually iterate/interact across pairs,
!     not just perform one isolated rotation), pair (2,3) and the diagonal
!     left exactly zero.
!
!     -M. M. F. Moraes, 2026.
!
!     USE Connections
!
      use MQC_EST
      use MQC_NEF
      use iso_fortran_env, only: int64, real64
!
!     Variable Declarations
!
      implicit none
      integer::nfail
      real(kind=real64),parameter::tol=1.0d-8
!
      type(mqc_scf_eigenvalues),dimension(:),allocatable::rowMap
      type(mqc_scf_integral)::fullInt
      type(mqc_nef_matrix)::NEFMat,rotFull,Wrot,Uout
      type(mqc_scalar)::sc
      complex(kind=real64)::cGot,cExpected
      integer::row,col,s,spinIdx
      character(len=5),dimension(2)::spinLabels
      character(len=64)::testLabel
!
!     Format Statements
!
 1000 format(1x,'PASS  ',a)
 1010 format(1x,'FAIL  ',a,'  |got|=',f16.10,'  |expected|=',f16.10)
 9000 format(/,1x,'Tests completed: ',i0,' failure(s)')
!
      nfail = 0
      spinLabels(1) = 'alpha'
      spinLabels(2) = 'beta'
!
!--------------------------------------------------------------------
!     Build a 3-space row/col map, 1 alpha + 1 beta orbital per space.
!--------------------------------------------------------------------
!
      allocate(rowMap(3))
      do s = 1, 3
        call rowMap(s)%init(0,0)
        call rowMap(s)%push(s,'alpha')
        call rowMap(s)%push(s,'beta')
      end do
!
!--------------------------------------------------------------------
!     Build a 3x3(alpha) + 3x3(beta) Hermitian flat integral with
!     genuinely complex, nonzero content on pairs (1,2) and (1,3),
!     zero everywhere else (including (2,3) and the diagonal).
!--------------------------------------------------------------------
!
      call fullInt%init(3,3)
      sc = cmplx(0.3d0,0.4d0,kind=real64)
      call fullInt%put(sc,1,2,'alpha')
      sc = cmplx(0.3d0,-0.4d0,kind=real64)
      call fullInt%put(sc,2,1,'alpha')
      sc = cmplx(-0.5d0,0.2d0,kind=real64)
      call fullInt%put(sc,1,3,'alpha')
      sc = cmplx(-0.5d0,-0.2d0,kind=real64)
      call fullInt%put(sc,3,1,'alpha')
!
      sc = cmplx(-0.2d0,0.1d0,kind=real64)
      call fullInt%put(sc,1,2,'beta')
      sc = cmplx(-0.2d0,-0.1d0,kind=real64)
      call fullInt%put(sc,2,1,'beta')
      sc = cmplx(0.4d0,-0.3d0,kind=real64)
      call fullInt%put(sc,1,3,'beta')
      sc = cmplx(0.4d0,0.3d0,kind=real64)
      call fullInt%put(sc,3,1,'beta')
!
      call NEFMat%init(rowMap)
      call NEFMat%imp(fullInt)
!
      rotFull = NEFMat%getRotation()
      Wrot = matmul(Dagger(rotFull),matmul(NEFMat,rotFull))
      Uout = matmul(Dagger(rotFull),rotFull)
!
!--------------------------------------------------------------------
!     Test 1: every inter-space off-diagonal block of Dagger(R)*W*R is
!     zero.
!--------------------------------------------------------------------
!
      do row = 1, 3
        do col = row+1, 3
          do spinIdx = 1, 2
            sc = Wrot%at(1,1,row,col,trim(spinLabels(spinIdx)))
            cGot = sc
            write(testLabel,'(a,i0,a,i0,a,a)') &
              'Dagger(R)*W*R block-diagonal: space(',row,',',col,&
              ') ',trim(spinLabels(spinIdx))
            if(abs(cGot).le.tol) then
              write(*,1000) trim(testLabel)
            else
              write(*,1010) trim(testLabel),abs(cGot),0.0d0
              nfail = nfail+1
            endIf
          end do
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 2: R is unitary (Dagger(R)*R = identity).
!--------------------------------------------------------------------
!
      do row = 1, 3
        do col = 1, 3
          do spinIdx = 1, 2
            sc = Uout%at(1,1,row,col,trim(spinLabels(spinIdx)))
            cGot = sc
            if(row.eq.col) then
              cExpected = cmplx(1.0d0,0.0d0,kind=real64)
            else
              cExpected = cmplx(0.0d0,0.0d0,kind=real64)
            end if
            write(testLabel,'(a,i0,a,i0,a,a)') &
              'Dagger(R)*R unitary: space(',row,',',col,&
              ') ',trim(spinLabels(spinIdx))
            if(abs(cGot-cExpected).le.tol) then
              write(*,1000) trim(testLabel)
            else
              write(*,1010) trim(testLabel),abs(cGot),abs(cExpected)
              nfail = nfail+1
            endIf
          end do
        end do
      end do
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1
!
      end program unitTest_NEF_matrix_rotation
