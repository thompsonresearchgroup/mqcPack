      program unitTest_NEF_collapse_expand
!
!     unitTest_NEF_collapse_expand: Round-trip tests for the MQC_NEF
!     collapse/expand primitives --
!       MQC_NEF_Matrix%collapse  / MQC_NEF_Vector%expand  (Matrix<->Vector)
!       MQC_NEF_Matrix%expand    / MQC_NEF_R4%collapse    (Matrix<->R4)
!       MQC_NEF_Set%collapse+%outerCollapse vs %fullCollapse (Set-level)
!     Tests 1-2 build a small MQC_NEF_Matrix with exactly 1 alpha + 1
!     beta orbital per block (every intra-block index is trivially 1).
!     Tests 3-4 repeat the same two round trips on a richer fixture
!     with 2 alpha + 2 beta orbitals per block, so the intra-block
!     collapseIndex arithmetic is exercised with varying indices.
!     Test 5 checks that MQC_NEF_Set%collapse() followed by
!     %outerCollapse() reproduces %fullCollapse()'s single-pass result
!     exactly, on a 2x2 grid of 2x2 MQC_NEF_Matrix cells.
!     Test 6/7 round-trip 'triplet'/'all' spin cases (needs a
!     'general'-type integral for real alphabeta/betaalpha content).
!     Test 8 checks collapse's swapOrder=.true. against the transpose
!     relationship it implies (expand has no swapOrder parameter).
!     Test 9 checks a non-default Order=[2,1,4,3] round trip.
!     Test 10 checks MQC_NEF_Set%outerCollapse's two jointIJ modes
!     against hand-derived expected values.
!     Test 11 round-trips Matrix<->R4 with BOTH combined-axis sizes
!     (map1/map2) simultaneously >1, unlike Tests 2/4/9's one-trivial-
!     axis fixtures.
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
      real(kind=real64),parameter::tol=1.0d-12
!
      type(mqc_scf_eigenvalues),dimension(:),allocatable::rowMap
      type(mqc_scf_eigenvalues),dimension(:),allocatable::map1,map2,map3,map4
      type(mqc_scf_eigenvalues),dimension(:),allocatable::rowMap3
      type(mqc_scf_eigenvalues),dimension(:),allocatable::map1b,map2b,map3b,map4b
      type(mqc_scf_eigenvalues),dimension(:),allocatable::rowMap4x
      type(mqc_scf_eigenvalues),dimension(:),allocatable::map1c,map2c,map3c,map4c
      type(mqc_scf_integral)::fullInt,fullInt3,cellInt,fullInt4
      type(mqc_nef_matrix)::NEFMat,MatBack,MatBack2
      type(mqc_nef_matrix)::NEFMat3,MatBack3,MatBack4,cellMat
      type(mqc_nef_matrix)::NEFMat4,MatBack5,MatBack6,MatBackSwap,MatBack7
      type(mqc_nef_matrix)::NEFMat5,MatBack8,outerMatDefault,outerMatJoint
      type(mqc_nef_vector)::VecOut,VecOut3,vecPathA,vecPathB
      type(mqc_nef_vector)::VecOut4,VecOut5,VecSwap
      type(mqc_nef_r4)::R4Out,R4Out2,R4Out3,R4Out3b
      type(mqc_nef_set)::NEFSet,SetCollapsed,SetOuter,SetFull
      type(mqc_nef_set)::SetOuterDefault,SetOuterJoint
      type(mqc_vector)::flatA,flatB
      type(mqc_scalar)::sc
      real(kind=real64)::got,expected,base,expectedVal
      integer::row,col,ispin,i,j,rr,cc,ci,cj,kk,nTot
      integer::origRow,origCol,newRow,newCol,G
      character(len=5),dimension(2)::spinLabels
      character(len=64)::testLabel
!
!     Format Statements
!
 1000 format(1x,'PASS  ',a)
 1010 format(1x,'FAIL  ',a,'  got=',f16.10,'  expected=',f16.10)
 9000 format(/,1x,'Tests completed: ',i0,' failure(s)')
!
      nfail = 0
      spinLabels(1) = 'alpha'
      spinLabels(2) = 'beta'
!
!--------------------------------------------------------------------
!     Build a 2x2 'space'-type integral with known, distinct values
!     and a matching 2-space NEF row/col map (1 alpha + 1 beta orbital
!     per space), then import it into a MQC_NEF_Matrix.
!--------------------------------------------------------------------
!
      call fullInt%init(2,2)
      sc = 1.0 ; call fullInt%put(sc,1,1)
      sc = 2.0 ; call fullInt%put(sc,1,2)
      sc = 3.0 ; call fullInt%put(sc,2,1)
      sc = 4.0 ; call fullInt%put(sc,2,2)
!
      allocate(rowMap(2))
      call rowMap(1)%init(0,0)
      call rowMap(1)%push(1,'alpha')
      call rowMap(1)%push(1,'beta')
      call rowMap(2)%init(0,0)
      call rowMap(2)%push(2,'alpha')
      call rowMap(2)%push(2,'beta')
!
      call NEFMat%init(rowMap)
      call NEFMat%imp(fullInt)
!
!--------------------------------------------------------------------
!     Test 1: Matrix -(collapse)-> Vector -(expand)-> Matrix round trip
!--------------------------------------------------------------------
!
      VecOut = NEFMat%collapse(spinType='singlet')
      MatBack = VecOut%expand(rowMap,rowMap,spinType='singlet')
!
      do row = 1, 2
        do col = 1, 2
          do ispin = 1, 2
            sc = MatBack%at(1,1,row,col,trim(spinLabels(ispin)))
            got = sc%rval()
            sc = NEFMat%at(1,1,row,col,trim(spinLabels(ispin)))
            expected = sc%rval()
            write(testLabel,'(a,i0,a,i0,a,a)') &
              'Matrix->Vector->Matrix: (',row,',',col,') ',trim(spinLabels(ispin))
            if(abs(got-expected).le.tol) then
              write(*,1000) trim(testLabel)
            else
              write(*,1010) trim(testLabel),got,expected
              nfail = nfail+1
            endIf
          end do
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 2: Matrix -(expand)-> R4 -(collapse)-> Matrix round trip
!     map1/map3 are trivial (1 space, 1 alpha+1 beta orbital) and
!     map2/map4 reproduce rowMap's 2-space structure, so that
!     NumSpace_row = size(map1)*size(map2) = 2 and likewise for col,
!     matching NEFMat's own 2x2 grid exactly.
!--------------------------------------------------------------------
!
      allocate(map1(1))
      call map1(1)%init(0,0)
      call map1(1)%push(1,'alpha')
      call map1(1)%push(1,'beta')
!
      allocate(map2(2))
      call map2(1)%init(0,0)
      call map2(1)%push(1,'alpha')
      call map2(1)%push(1,'beta')
      call map2(2)%init(0,0)
      call map2(2)%push(2,'alpha')
      call map2(2)%push(2,'beta')
!
      map3 = map1
      map4 = map2
!
      R4Out = NEFMat%expand(map1,map2,map3,map4,spinType='singlet')
      MatBack2 = R4Out%collapse(spinType='singlet')
!
      do row = 1, 2
        do col = 1, 2
          do ispin = 1, 2
            sc = MatBack2%at(1,1,row,col,trim(spinLabels(ispin)))
            got = sc%rval()
            sc = NEFMat%at(1,1,row,col,trim(spinLabels(ispin)))
            expected = sc%rval()
            write(testLabel,'(a,i0,a,i0,a,a)') &
              'Matrix->R4->Matrix: (',row,',',col,') ',trim(spinLabels(ispin))
            if(abs(got-expected).le.tol) then
              write(*,1000) trim(testLabel)
            else
              write(*,1010) trim(testLabel),got,expected
              nfail = nfail+1
            endIf
          end do
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 3: multi-orbital-per-block Matrix -(collapse)-> Vector
!     -(expand)-> Matrix round trip. Each of the 2 spaces now holds
!     2 alpha + 2 beta orbitals, so the intra-block collapseIndex
!     arithmetic is exercised with i,j both varying (1,2), unlike
!     Tests 1-2 above where every block held exactly 1 orbital.
!--------------------------------------------------------------------
!
      call fullInt3%init(4,4)
      do rr = 1, 4
        do cc = 1, 4
          sc = real((rr-1)*4+cc)
          call fullInt3%put(sc,rr,cc)
        end do
      end do
!
      allocate(rowMap3(2))
      call rowMap3(1)%init(0,0)
      call rowMap3(1)%push(1,'alpha')
      call rowMap3(1)%push(2,'alpha')
      call rowMap3(1)%push(1,'beta')
      call rowMap3(1)%push(2,'beta')
      call rowMap3(2)%init(0,0)
      call rowMap3(2)%push(3,'alpha')
      call rowMap3(2)%push(4,'alpha')
      call rowMap3(2)%push(3,'beta')
      call rowMap3(2)%push(4,'beta')
!
      call NEFMat3%init(rowMap3)
      call NEFMat3%imp(fullInt3)
!
      VecOut3 = NEFMat3%collapse(spinType='singlet')
      MatBack3 = VecOut3%expand(rowMap3,rowMap3,spinType='singlet')
!
      do row = 1, 2
        do col = 1, 2
          do i = 1, 2
            do j = 1, 2
              do ispin = 1, 2
                sc = MatBack3%at(i,j,row,col,trim(spinLabels(ispin)))
                got = sc%rval()
                sc = NEFMat3%at(i,j,row,col,trim(spinLabels(ispin)))
                expected = sc%rval()
                write(testLabel,'(a,i0,a,i0,a,i0,a,i0,a,a)') &
                  'multi-orbital Matrix->Vector->Matrix: (',row,',',col,&
                  ') [',i,',',j,'] ',trim(spinLabels(ispin))
                if(abs(got-expected).le.tol) then
                  write(*,1000) trim(testLabel)
                else
                  write(*,1010) trim(testLabel),got,expected
                  nfail = nfail+1
                endIf
              end do
            end do
          end do
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 4: same multi-orbital NEFMat3, round-tripped through R4
!     instead (Matrix -(expand)-> R4 -(collapse)-> Matrix). map1b/
!     map3b are trivial (1 space, 1 alpha+1 beta orbital); map2b/
!     map4b carry the real 2-space/2-orbital structure, so that
!     size(map1b)*size(map2b) = 2 = NEFMat3's own row/col count and
!     the per-block orbital counts line up exactly.
!--------------------------------------------------------------------
!
      allocate(map1b(1))
      call map1b(1)%init(0,0)
      call map1b(1)%push(1,'alpha')
      call map1b(1)%push(1,'beta')
!
      allocate(map2b(2))
      call map2b(1)%init(0,0)
      call map2b(1)%push(1,'alpha')
      call map2b(1)%push(2,'alpha')
      call map2b(1)%push(1,'beta')
      call map2b(1)%push(2,'beta')
      call map2b(2)%init(0,0)
      call map2b(2)%push(3,'alpha')
      call map2b(2)%push(4,'alpha')
      call map2b(2)%push(3,'beta')
      call map2b(2)%push(4,'beta')
!
      map3b = map1b
      map4b = map2b
!
      R4Out2 = NEFMat3%expand(map1b,map2b,map3b,map4b,spinType='singlet')
      MatBack4 = R4Out2%collapse(spinType='singlet')
!
      do row = 1, 2
        do col = 1, 2
          do i = 1, 2
            do j = 1, 2
              do ispin = 1, 2
                sc = MatBack4%at(i,j,row,col,trim(spinLabels(ispin)))
                got = sc%rval()
                sc = NEFMat3%at(i,j,row,col,trim(spinLabels(ispin)))
                expected = sc%rval()
                write(testLabel,'(a,i0,a,i0,a,i0,a,i0,a,a)') &
                  'multi-orbital Matrix->R4->Matrix: (',row,',',col,&
                  ') [',i,',',j,'] ',trim(spinLabels(ispin))
                if(abs(got-expected).le.tol) then
                  write(*,1000) trim(testLabel)
                else
                  write(*,1010) trim(testLabel),got,expected
                  nfail = nfail+1
                endIf
              end do
            end do
          end do
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 5: MQC_NEF_Set collapse + outerCollapse vs fullCollapse.
!     Builds a 2x2 NEF_Set of 2x2 NEF_Matrix cells (each cell reusing
!     rowMap from Tests 1-2, with distinct values per cell), then
!     checks that collapse()+outerCollapse() (two steps, default
!     jointIJ) reproduces fullCollapse() (one condensed step) exactly,
!     element-by-element over the flattened combined Vector.
!--------------------------------------------------------------------
!
      call NEFSet%init(2,2)
      do ci = 1, 2
        do cj = 1, 2
          base = real(10*((ci-1)*2+cj))
          call cellInt%init(2,2)
          sc = base+1.0 ; call cellInt%put(sc,1,1)
          sc = base+2.0 ; call cellInt%put(sc,1,2)
          sc = base+3.0 ; call cellInt%put(sc,2,1)
          sc = base+4.0 ; call cellInt%put(sc,2,2)
          call cellMat%init(rowMap)
          call cellMat%imp(cellInt)
          call NEFSet%put(cellMat,ci,cj)
        end do
      end do
!
      SetCollapsed = NEFSet%collapse(spinType='singlet')
      SetOuter = SetCollapsed%outerCollapse()
      SetFull = NEFSet%fullCollapse(spinType='singlet')
!
      vecPathA = SetOuter%atVector(1,1)
      vecPathB = SetFull%atVector(1,1)
!
      flatA = vecPathA%getBlockVec()
      flatB = vecPathB%getBlockVec()
!
      if(size(flatA).ne.size(flatB)) then
        write(*,1010) 'collapse+outerCollapse vs fullCollapse: size mismatch',&
          real(size(flatA)),real(size(flatB))
        nfail = nfail+1
      else
        nTot = size(flatA)
        do kk = 1, nTot
          sc = flatA%at(kk)
          got = sc%rval()
          sc = flatB%at(kk)
          expected = sc%rval()
          write(testLabel,'(a,i0)') 'collapse+outerCollapse vs fullCollapse: elem ',kk
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected
            nfail = nfail+1
          endIf
        end do
      end if
!
!--------------------------------------------------------------------
!     Test 6: 'triplet' Matrix -(collapse)-> Vector -(expand)-> Matrix
!     round trip. Needs real alphabeta/betaalpha content, so NEFMat4
!     is built from a 'general'-type integral (distinct values in all
!     four spin blocks). Triplet only carries alphabeta/betaalpha
!     content; pure alpha/beta blocks are untouched by this path.
!--------------------------------------------------------------------
!
      call fullInt4%init(2,2)
      sc = 1.0  ; call fullInt4%put(sc,1,1,'alpha')
      sc = 2.0  ; call fullInt4%put(sc,1,2,'alpha')
      sc = 3.0  ; call fullInt4%put(sc,2,1,'alpha')
      sc = 4.0  ; call fullInt4%put(sc,2,2,'alpha')
      sc = 5.0  ; call fullInt4%put(sc,1,1,'beta')
      sc = 6.0  ; call fullInt4%put(sc,1,2,'beta')
      sc = 7.0  ; call fullInt4%put(sc,2,1,'beta')
      sc = 8.0  ; call fullInt4%put(sc,2,2,'beta')
      sc = 9.0  ; call fullInt4%put(sc,1,1,'alphabeta')
      sc = 10.0 ; call fullInt4%put(sc,1,2,'alphabeta')
      sc = 11.0 ; call fullInt4%put(sc,2,1,'alphabeta')
      sc = 12.0 ; call fullInt4%put(sc,2,2,'alphabeta')
      sc = 13.0 ; call fullInt4%put(sc,1,1,'betaalpha')
      sc = 14.0 ; call fullInt4%put(sc,1,2,'betaalpha')
      sc = 15.0 ; call fullInt4%put(sc,2,1,'betaalpha')
      sc = 16.0 ; call fullInt4%put(sc,2,2,'betaalpha')
!
      call NEFMat4%init(rowMap)
      call NEFMat4%imp(fullInt4)
!
      VecOut4 = NEFMat4%collapse(spinType='triplet')
      MatBack5 = VecOut4%expand(rowMap,rowMap,spinType='triplet')
!
      do row = 1, 2
        do col = 1, 2
          sc = MatBack5%at(1,1,row,col,'alphabeta')
          got = sc%rval()
          sc = NEFMat4%at(1,1,row,col,'alphabeta')
          expected = sc%rval()
          write(testLabel,'(a,i0,a,i0,a)') &
            'triplet Matrix->Vector->Matrix: (',row,',',col,') alphabeta'
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected
            nfail = nfail+1
          endIf
!
          sc = MatBack5%at(1,1,row,col,'betaalpha')
          got = sc%rval()
          sc = NEFMat4%at(1,1,row,col,'betaalpha')
          expected = sc%rval()
          write(testLabel,'(a,i0,a,i0,a)') &
            'triplet Matrix->Vector->Matrix: (',row,',',col,') betaalpha'
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected
            nfail = nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 7: 'all' Matrix -(collapse)-> Vector -(expand)-> Matrix
!     round trip, same general-type NEFMat4. 'all' reconstructs every
!     spin block (alpha, beta, alphabeta, betaalpha).
!--------------------------------------------------------------------
!
      VecOut5 = NEFMat4%collapse(spinType='all')
      MatBack6 = VecOut5%expand(rowMap,rowMap,spinType='all')
!
      do row = 1, 2
        do col = 1, 2
          sc = MatBack6%at(1,1,row,col,'alpha')
          got = sc%rval() ; sc = NEFMat4%at(1,1,row,col,'alpha') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0,a)') 'all Matrix->Vector->Matrix: (',row,',',col,') alpha'
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
!
          sc = MatBack6%at(1,1,row,col,'beta')
          got = sc%rval() ; sc = NEFMat4%at(1,1,row,col,'beta') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0,a)') 'all Matrix->Vector->Matrix: (',row,',',col,') beta'
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
!
          sc = MatBack6%at(1,1,row,col,'alphabeta')
          got = sc%rval() ; sc = NEFMat4%at(1,1,row,col,'alphabeta') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0,a)') 'all Matrix->Vector->Matrix: (',row,',',col,') alphabeta'
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
!
          sc = MatBack6%at(1,1,row,col,'betaalpha')
          got = sc%rval() ; sc = NEFMat4%at(1,1,row,col,'betaalpha') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0,a)') 'all Matrix->Vector->Matrix: (',row,',',col,') betaalpha'
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 8: swapOrder=.true. on Matrix%collapse. Expand has no
!     swapOrder parameter, so re-expanding a swap-collapsed Vector
!     with plain expand() re-reads it via the UN-swapped formula; for
!     a square grid this yields the TRANSPOSE of the original matrix
!     (block(row,col) ends up holding the original's block(col,row)).
!--------------------------------------------------------------------
!
      VecSwap = NEFMat4%collapse(spinType='singlet',swapOrder=.true.)
      MatBackSwap = VecSwap%expand(rowMap,rowMap,spinType='singlet')
!
      do row = 1, 2
        do col = 1, 2
          sc = MatBackSwap%at(1,1,row,col,'alpha')
          got = sc%rval()
          sc = NEFMat4%at(1,1,col,row,'alpha')
          expected = sc%rval()
          write(testLabel,'(a,i0,a,i0,a)') &
            'swapOrder transpose check: (',row,',',col,') alpha'
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
!
          sc = MatBackSwap%at(1,1,row,col,'beta')
          got = sc%rval()
          sc = NEFMat4%at(1,1,col,row,'beta')
          expected = sc%rval()
          write(testLabel,'(a,i0,a,i0,a)') &
            'swapOrder transpose check: (',row,',',col,') beta'
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 9: non-default Order=[2,1,4,3] on Matrix%expand / R4%collapse.
!     Order only relabels which physical R4 axis each role (ia,ja,ka,la)
!     lands on; using the SAME Order on both legs keeps the round trip
!     self-consistent regardless of permutation. Reuses Test 2's
!     map1-4/NEFMat fixture.
!--------------------------------------------------------------------
!
      R4Out3 = NEFMat%expand(map1,map2,map3,map4,spinType='singlet',&
                              Order=[2_int64,1_int64,4_int64,3_int64])
      MatBack7 = R4Out3%collapse(spinType='singlet',&
                                  Order=[2_int64,1_int64,4_int64,3_int64])
!
      do row = 1, 2
        do col = 1, 2
          do ispin = 1, 2
            sc = MatBack7%at(1,1,row,col,trim(spinLabels(ispin)))
            got = sc%rval()
            sc = NEFMat%at(1,1,row,col,trim(spinLabels(ispin)))
            expected = sc%rval()
            write(testLabel,'(a,i0,a,i0,a,a)') &
              'Order=[2,1,4,3] Matrix->R4->Matrix: (',row,',',col,') ',trim(spinLabels(ispin))
            if(abs(got-expected).le.tol) then
              write(*,1000) trim(testLabel)
            else
              write(*,1010) trim(testLabel),got,expected
              nfail = nfail+1
            endIf
          end do
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 10: MQC_NEF_Set%outerCollapse on the MATRIX-kind NEFSet
!     from Test 5, both jointIJ modes, checked against hand-derived
!     expected values (not just a round trip).
!       jointIJ=.false. (default): newRow=(I-1)*origNumSpace_row+r,
!         newCol=(J-1)*origNumSpace_col+c (I,J independently fold).
!       jointIJ=.true. (jointAxis=1 default): G=(I-1)*NumJ+J folds
!         into the row; col axis is left untouched (=c).
!--------------------------------------------------------------------
!
      SetOuterDefault = NEFSet%outerCollapse()
      outerMatDefault = SetOuterDefault%atMatrix(1,1)
      do ci = 1, 2
        do cj = 1, 2
          do origRow = 1, 2
            do origCol = 1, 2
              base = real(10*((ci-1)*2+cj))
              expectedVal = base + real((origRow-1)*2+origCol)
              newRow = (ci-1)*2+origRow
              newCol = (cj-1)*2+origCol
              sc = outerMatDefault%at(1,1,newRow,newCol,'alpha')
              got = sc%rval()
              write(testLabel,'(a,i0,a,i0,a)') &
                'outerCollapse jointIJ=F: (',newRow,',',newCol,') alpha'
              if(abs(got-expectedVal).le.tol) then
                write(*,1000) trim(testLabel)
              else
                write(*,1010) trim(testLabel),got,expectedVal
                nfail = nfail+1
              endIf
            end do
          end do
        end do
      end do
!
      SetOuterJoint = NEFSet%outerCollapse(jointIJ=.true.)
      outerMatJoint = SetOuterJoint%atMatrix(1,1)
      do ci = 1, 2
        do cj = 1, 2
          do origRow = 1, 2
            do origCol = 1, 2
              base = real(10*((ci-1)*2+cj))
              expectedVal = base + real((origRow-1)*2+origCol)
              G = (ci-1)*2+cj
              newRow = (G-1)*2+origRow
              newCol = origCol
              sc = outerMatJoint%at(1,1,newRow,newCol,'alpha')
              got = sc%rval()
              write(testLabel,'(a,i0,a,i0,a)') &
                'outerCollapse jointIJ=T: (',newRow,',',newCol,') alpha'
              if(abs(got-expectedVal).le.tol) then
                write(*,1000) trim(testLabel)
              else
                write(*,1010) trim(testLabel),got,expectedVal
                nfail = nfail+1
              endIf
            end do
          end do
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 11: Matrix -(expand)-> R4 -(collapse)-> Matrix round trip
!     where BOTH combined-axis sizes are >1 simultaneously (map1c,
!     map2c each 2 spaces), unlike Tests 2/4/9 where one side was
!     always a trivial size-1 dummy axis. NEFMat5 is a 4x4 single-
!     orbital-per-block matrix (NumSpace_row=NumSpace_col=4=
!     size(map1c)*size(map2c)).
!--------------------------------------------------------------------
!
      allocate(rowMap4x(4))
      do kk = 1, 4
        call rowMap4x(kk)%init(0,0)
        call rowMap4x(kk)%push(kk,'alpha')
        call rowMap4x(kk)%push(kk,'beta')
      end do
      call NEFMat5%init(rowMap4x)
      call NEFMat5%imp(fullInt3)
!
      allocate(map1c(2))
      call map1c(1)%init(0,0)
      call map1c(1)%push(1,'alpha')
      call map1c(1)%push(1,'beta')
      call map1c(2)%init(0,0)
      call map1c(2)%push(1,'alpha')
      call map1c(2)%push(1,'beta')
!
      allocate(map2c(2))
      call map2c(1)%init(0,0)
      call map2c(1)%push(1,'alpha')
      call map2c(1)%push(1,'beta')
      call map2c(2)%init(0,0)
      call map2c(2)%push(1,'alpha')
      call map2c(2)%push(1,'beta')
!
      map3c = map1c
      map4c = map2c
!
      R4Out3b = NEFMat5%expand(map1c,map2c,map3c,map4c,spinType='singlet')
      MatBack8 = R4Out3b%collapse(spinType='singlet')
!
      do row = 1, 4
        do col = 1, 4
          do ispin = 1, 2
            sc = MatBack8%at(1,1,row,col,trim(spinLabels(ispin)))
            got = sc%rval()
            sc = NEFMat5%at(1,1,row,col,trim(spinLabels(ispin)))
            expected = sc%rval()
            write(testLabel,'(a,i0,a,i0,a,a)') &
              'R4 both-axes-nontrivial Matrix->R4->Matrix: (',row,',',col,') ',&
              trim(spinLabels(ispin))
            if(abs(got-expected).le.tol) then
              write(*,1000) trim(testLabel)
            else
              write(*,1010) trim(testLabel),got,expected
              nfail = nfail+1
            endIf
          end do
        end do
      end do
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1
!
      end program unitTest_NEF_collapse_expand
