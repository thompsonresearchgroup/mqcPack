      program unitTest_NEF_merge_unmerge_roundtrip
!
!     unitTest_NEF_merge_unmerge_roundtrip: checks that
!     MQC_NEF_Set%mergeSpace followed by MQC_NEF_Set%separateSpace
!     (called on the merged result, with the ORIGINAL set passed as the
!     template) exactly reproduces the original (NumI x NumJ) grid of
!     NEF_Vector/NEF_Matrix cells -- a pure data-movement round trip, so
!     every recovered value must match the original exactly (to the
!     tolerance below).
!
!     Both fixtures are deliberately NON-SYMMETRICAL (NumI /= NumJ, and
!     for the matrix fixture, NumSpace_row /= NumSpace_col) with
!     NON-UNIFORM per-block alpha/beta lengths (sizes differ block to
!     block), and every orbital index across every block/spin/row/col is
!     given a GLOBALLY distinct raw index so that any axis-offset or
!     row/col-swap bug would land a value in the wrong slot rather than
!     accidentally landing on an equal one. Every element is a distinct
!     complex number encoding (I,J,block/row,col,el) so a wrong-slot bug
!     is guaranteed to be visible.
!
!     Vector fixture: tests both jointIJ=.true. (grid folds into ONE
!     combined NEF_Vector) and jointIJ=.false. (grid's I-axis folds into
!     row blocks, J-axis becomes the single alpha-labeled matrix column,
!     per the "J always lands on alpha" convention).
!
!     Matrix fixture: tests jointIJ=.false. (both I and J fold, I into
!     rows, J into columns), and jointIJ=.true. with jointAxis=1 (both
!     fold into rows) and jointAxis=2 (both fold into columns).
!
!     -M. M. F. Moraes, 2026.
!
      use MQC_EST
      use MQC_NEF
      use iso_fortran_env, only: int64, real64
!
      implicit none
      integer::nfail
      real(kind=real64),parameter::tol=1.0d-10
!
      type(mqc_scf_eigenvalues),dimension(:),allocatable::rowMap,colMap
      type(mqc_nef_set)::origVecSet,mergedVecSet,roundVecSet
      type(mqc_nef_set)::origMatSet,mergedMatSet,roundMatSet
      type(mqc_nef_vector)::cellVec,roundCellVec
      type(mqc_nef_matrix)::cellMat,roundCellMat
      type(mqc_vector)::flatVec,blockOrig,blockRound
      type(mqc_scf_integral)::flatInt,blockInt1,blockInt2
      type(mqc_scalar)::sc
      integer(kind=int64)::I,J,b,r,c,el
      integer(kind=int64),dimension(:),allocatable::alphaIdx,betaIdx
      complex(kind=real64)::cGot,cExpected
      character(len=120)::testLabel
      integer::nRowBlk,nColBlk,nRowBasis,nColBasis
      logical::myJointIJ
      integer::myJointAxis
!
 1000 format(1x,'PASS  ',a)
 1010 format(1x,'FAIL  ',a,'  got=(',f14.8,',',f14.8,')  expected=(',f14.8,',',f14.8,')')
 9000 format(/,1x,'Tests completed: ',i0,' failure(s)')
!
      nfail = 0
!
!--------------------------------------------------------------------
!     VECTOR fixture: 2 blocks, non-uniform sizes, globally distinct
!     raw orbital indices (1..6) so every (block,spin,el) slot is
!     unambiguous in the flat vector used by impv.
!       block 1: alpha=[1,2] (len 2), beta=[3]   (len 1)
!       block 2: alpha=[4]   (len 1), beta=[5,6] (len 2)
!     NumI=3, NumJ=2 (non-symmetrical grid).
!--------------------------------------------------------------------
!
      allocate(rowMap(2))
      call rowMap(1)%init(0,0)
      call rowMap(1)%push(1_int64,'alpha'); call rowMap(1)%push(2_int64,'alpha')
      call rowMap(1)%push(3_int64,'beta')
      call rowMap(2)%init(0,0)
      call rowMap(2)%push(4_int64,'alpha')
      call rowMap(2)%push(5_int64,'beta'); call rowMap(2)%push(6_int64,'beta')
!
      call origVecSet%init(3,2)
      do I = 1, 3
        do J = 1, 2
          call flatVec%init(6)
          do el = 1, 6
            sc = cmplx(10.0d0*I+1.0d0*J+0.1d0*el,&
                 5.0d0*I-2.0d0*J-0.05d0*el,kind=real64)
            call flatVec%put(sc,el)
          end do
          call cellVec%init(rowMap)
          call cellVec%impv(flatVec)
          call origVecSet%put(cellVec,int(I),int(J))
        end do
      end do
!
!     jointIJ = .true.
!
      mergedVecSet = origVecSet%mergeSpace(jointIJ=.true.)
      roundVecSet = mergedVecSet%separateSpace(origVecSet,jointIJ=.true.)
      do I = 1, 3
        do J = 1, 2
          cellVec = origVecSet%atVector(int(I),int(J))
          roundCellVec = roundVecSet%atVector(int(I),int(J))
          do b = 1, 2
            blockOrig = cellVec%getBlockVec(indexes=[b],blockName='alpha')
            blockRound = roundCellVec%getBlockVec(indexes=[b],blockName='alpha')
            do el = 1, size(blockOrig)
              cGot = blockRound%at(el); cExpected = blockOrig%at(el)
              write(testLabel,'(a,i0,a,i0,a,i0,a,i0)') &
                'vector jointIJ=T roundtrip: I=',I,' J=',J,' block=',b,' alpha el=',el
              if(abs(cGot-cExpected).le.tol) then
                write(*,1000) trim(testLabel)
              else
                write(*,1010) trim(testLabel),real(cGot),aimag(cGot),real(cExpected),aimag(cExpected)
                nfail = nfail+1
              end if
            end do
            blockOrig = cellVec%getBlockVec(indexes=[b],blockName='beta')
            blockRound = roundCellVec%getBlockVec(indexes=[b],blockName='beta')
            do el = 1, size(blockOrig)
              cGot = blockRound%at(el); cExpected = blockOrig%at(el)
              write(testLabel,'(a,i0,a,i0,a,i0,a,i0)') &
                'vector jointIJ=T roundtrip: I=',I,' J=',J,' block=',b,' beta el=',el
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
!
!     jointIJ = .false.
!
      mergedVecSet = origVecSet%mergeSpace(jointIJ=.false.)
      roundVecSet = mergedVecSet%separateSpace(origVecSet,jointIJ=.false.)
      do I = 1, 3
        do J = 1, 2
          cellVec = origVecSet%atVector(int(I),int(J))
          roundCellVec = roundVecSet%atVector(int(I),int(J))
          do b = 1, 2
            blockOrig = cellVec%getBlockVec(indexes=[b],blockName='alpha')
            blockRound = roundCellVec%getBlockVec(indexes=[b],blockName='alpha')
            do el = 1, size(blockOrig)
              cGot = blockRound%at(el); cExpected = blockOrig%at(el)
              write(testLabel,'(a,i0,a,i0,a,i0,a,i0)') &
                'vector jointIJ=F roundtrip: I=',I,' J=',J,' block=',b,' alpha el=',el
              if(abs(cGot-cExpected).le.tol) then
                write(*,1000) trim(testLabel)
              else
                write(*,1010) trim(testLabel),real(cGot),aimag(cGot),real(cExpected),aimag(cExpected)
                nfail = nfail+1
              end if
            end do
            blockOrig = cellVec%getBlockVec(indexes=[b],blockName='beta')
            blockRound = roundCellVec%getBlockVec(indexes=[b],blockName='beta')
            do el = 1, size(blockOrig)
              cGot = blockRound%at(el); cExpected = blockOrig%at(el)
              write(testLabel,'(a,i0,a,i0,a,i0,a,i0)') &
                'vector jointIJ=F roundtrip: I=',I,' J=',J,' block=',b,' beta el=',el
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
!
!--------------------------------------------------------------------
!     MATRIX fixture: non-square block structure (2 row-blocks, 3
!     col-blocks), non-uniform sizes, globally distinct raw indices.
!       row block 1: alpha=[1,2] (len 2), beta=[3]   (len 1)
!       row block 2: alpha=[4]   (len 1), beta=[5,6] (len 2)
!       col block 1: alpha=[1]   (len 1), beta=[2]   (len 1)
!       col block 2: alpha=[3,4] (len 2), beta=[5]   (len 1)
!       col block 3: alpha=[6]   (len 1), beta=[7,8] (len 2)
!     NumI=3, NumJ=2 (non-symmetrical grid), 'spin'-type (distinct
!     alpha/beta values, no alphaBeta/betaAlpha cross terms).
!--------------------------------------------------------------------
!
      allocate(colMap(3))
      call colMap(1)%init(0,0)
      call colMap(1)%push(1_int64,'alpha')
      call colMap(1)%push(2_int64,'beta')
      call colMap(2)%init(0,0)
      call colMap(2)%push(3_int64,'alpha'); call colMap(2)%push(4_int64,'alpha')
      call colMap(2)%push(5_int64,'beta')
      call colMap(3)%init(0,0)
      call colMap(3)%push(6_int64,'alpha')
      call colMap(3)%push(7_int64,'beta'); call colMap(3)%push(8_int64,'beta')
!
      nRowBasis = 6
      nColBasis = 8
      call origMatSet%init(3,2)
      do I = 1, 3
        do J = 1, 2
          call flatInt%init(nRowBasis,nRowBasis,nAlpha2=nColBasis,nBeta2=nColBasis)
          do r = 1, nRowBasis
            do c = 1, nColBasis
              sc = cmplx(10.0d0*I+1.0d0*J+0.1d0*r+0.01d0*c,&
                   5.0d0*I-2.0d0*J-0.05d0*r+0.02d0*c,kind=real64)
              call flatInt%put(sc,r,c,'alpha')
              sc = cmplx(-10.0d0*I-1.0d0*J+0.1d0*r+0.01d0*c,&
                   -5.0d0*I+2.0d0*J-0.05d0*r+0.02d0*c,kind=real64)
              call flatInt%put(sc,r,c,'beta')
            end do
          end do
          call cellMat%init(rowMap,initial_map_matrix_col=colMap)
          call cellMat%imp(flatInt)
          call origMatSet%put(cellMat,int(I),int(J))
        end do
      end do
!
!     jointIJ = .false.
!
      mergedMatSet = origMatSet%mergeSpace(jointIJ=.false.)
      roundMatSet = mergedMatSet%separateSpace(origMatSet,jointIJ=.false.)
      call compareMatrixSets(origMatSet,roundMatSet,'matrix jointIJ=F')
!
!     jointIJ = .true., jointAxis = 1
!
      mergedMatSet = origMatSet%mergeSpace(jointIJ=.true.,jointAxis=1)
      roundMatSet = mergedMatSet%separateSpace(origMatSet,jointIJ=.true.,jointAxis=1)
      call compareMatrixSets(origMatSet,roundMatSet,'matrix jointIJ=T axis=1')
!
!     jointIJ = .true., jointAxis = 2
!
      mergedMatSet = origMatSet%mergeSpace(jointIJ=.true.,jointAxis=2)
      roundMatSet = mergedMatSet%separateSpace(origMatSet,jointIJ=.true.,jointAxis=2)
      call compareMatrixSets(origMatSet,roundMatSet,'matrix jointIJ=T axis=2')
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1
!
      contains
!
      subroutine compareMatrixSets(setA,setB,label)
      implicit none
      type(mqc_nef_set),intent(in)::setA,setB
      character(len=*),intent(in)::label
      type(mqc_nef_matrix)::matA,matB
      type(mqc_scf_integral)::blkA,blkB
      integer(kind=int64)::ii,jj,rr,cc,pp,qq
      complex(kind=real64)::got,expected
      character(len=160)::lbl
!
      do ii = 1, 3
        do jj = 1, 2
          matA = setA%atMatrix(int(ii),int(jj))
          matB = setB%atMatrix(int(ii),int(jj))
          do rr = 1, 2
            do cc = 1, 3
              blkA = matA%iat(int(rr),int(cc))
              blkB = matB%iat(int(rr),int(cc))
              do pp = 1, blkA%blockSize('alpha',1)
                do qq = 1, blkA%blockSize('alpha',2)
                  got = blkB%at(pp,qq,'alpha'); expected = blkA%at(pp,qq,'alpha')
                  write(lbl,'(a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') trim(label),&
                    ' roundtrip: I=',ii,' J=',jj,' row=',rr,' col=',cc,' p=',pp,' q=',qq
                  if(abs(got-expected).le.tol) then
                    write(*,'(1x,a,a)') 'PASS  ',trim(lbl)//' alpha'
                  else
                    write(*,'(1x,a,a,a,f14.8,a,f14.8,a,f14.8,a,f14.8,a)') 'FAIL  ',trim(lbl)//' alpha',&
                      '  got=(',real(got),',',aimag(got),')  expected=(',real(expected),',',aimag(expected),')'
                    nfail = nfail+1
                  end if
                end do
              end do
              do pp = 1, blkA%blockSize('beta',1)
                do qq = 1, blkA%blockSize('beta',2)
                  got = blkB%at(pp,qq,'beta'); expected = blkA%at(pp,qq,'beta')
                  write(lbl,'(a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') trim(label),&
                    ' roundtrip: I=',ii,' J=',jj,' row=',rr,' col=',cc,' p=',pp,' q=',qq
                  if(abs(got-expected).le.tol) then
                    write(*,'(1x,a,a)') 'PASS  ',trim(lbl)//' beta'
                  else
                    write(*,'(1x,a,a,a,f14.8,a,f14.8,a,f14.8,a,f14.8,a)') 'FAIL  ',trim(lbl)//' beta',&
                      '  got=(',real(got),',',aimag(got),')  expected=(',real(expected),',',aimag(expected),')'
                    nfail = nfail+1
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
      end subroutine compareMatrixSets
!
      end program unitTest_NEF_merge_unmerge_roundtrip
