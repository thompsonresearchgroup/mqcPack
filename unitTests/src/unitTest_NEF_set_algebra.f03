      program unitTest_NEF_set_algebra
!
!     unitTest_NEF_set_algebra: Tests for the MQC_NEF_Set algebra built
!     this session on top of MQC_NEF_Vector/Matrix/R4 -- MatMul against
!     plain mqc_vector/mqc_matrix and against a single shared
!     MQC_NEF_Matrix (full grid reduction), MatMul between two NEF_Sets
!     (outer-product inner op, rank-raising), Contraction/
!     FullContraction/PartialContraction/PartialContraction2R4 between
!     two NEF_Sets (grid-level (I,J)x(J,K)=(I,K), only the shared J axis
!     contracted). Grids are non-square (2x3 / 3x2) throughout; inner
!     cells are non-square too wherever that is the actual point of the
!     test (the new mixed Matrix/Vector PartialContraction variants).
!
!     Each grid-level function is checked by manually replicating its
!     "sum over J of <base-level op>" logic in this file using the same
!     (already-established) base-level MatMul/Contraction/
!     PartialContraction/outer2Matrix/outer2R4 calls, then comparing
!     against the new grid-level function's result -- this verifies the
!     new orchestration code directly without re-deriving the underlying
!     physics (e.g. R4,R4 Contraction's coulomb/exchange diagrams) from
!     scratch.
!
!     -M. M. F. Moraes, 2026.
!
!     USE Connections
!
      use MQC_EST
      use MQC_NEF
      use iso_fortran_env, only: int64, real64
!
      implicit none
      integer::nfail
      real(kind=real64),parameter::tol=1.0d-10
!
      type(mqc_scf_eigenvalues),dimension(:),allocatable::map2,map3
      type(mqc_scf_eigenvalues),dimension(:),allocatable::mapTriv1,mapTriv2
      type(mqc_scf_integral)::fullInt3,cellInt
      type(mqc_nef_matrix)::cellMat,gotMat,expMat
      type(mqc_nef_vector)::cellVec,gotVec,expVec
      type(mqc_nef_r4)::cellR4,gotR4,expR4
      type(mqc_nef_set)::NEFSetA,NEFSetB,NEFSetC
      type(mqc_nef_set)::NEFSetAv,NEFSetBv,NEFSetCv
      type(mqc_nef_set)::NEFSetAr4,NEFSetBr4,NEFSetCr4
      type(mqc_nef_set)::NEFSetR,NEFSetS
      type(mqc_nef_set)::NEFSetOut
      type(mqc_vector)::plainVec3,plainVec2
      type(mqc_matrix)::plainMat32,plainMat22,outMat
      type(mqc_scalar)::sc,outScalar,expScalar
      real(kind=real64)::got,expected,base
      integer::I,J,K,ii,jj
      character(len=64)::testLabel
!
 1000 format(1x,'PASS  ',a)
 1010 format(1x,'FAIL  ',a,'  got=',f18.8,'  expected=',f18.8)
 9000 format(/,1x,'Tests completed: ',i0,' failure(s)')
!
      nfail = 0
!
!--------------------------------------------------------------------
!     Shared building blocks.
!--------------------------------------------------------------------
!
      allocate(map2(2))
      call map2(1)%init(0,0) ; call map2(1)%push(1,'alpha') ; call map2(1)%push(1,'beta')
      call map2(2)%init(0,0) ; call map2(2)%push(2,'alpha') ; call map2(2)%push(2,'beta')
!
      allocate(map3(3))
      call map3(1)%init(0,0) ; call map3(1)%push(1,'alpha') ; call map3(1)%push(1,'beta')
      call map3(2)%init(0,0) ; call map3(2)%push(2,'alpha') ; call map3(2)%push(2,'beta')
      call map3(3)%init(0,0) ; call map3(3)%push(3,'alpha') ; call map3(3)%push(3,'beta')
!
      allocate(mapTriv1(1))
      call mapTriv1(1)%init(0,0) ; call mapTriv1(1)%push(1,'alpha') ; call mapTriv1(1)%push(1,'beta')
      allocate(mapTriv2(2))
      call mapTriv2(1)%init(0,0) ; call mapTriv2(1)%push(1,'alpha') ; call mapTriv2(1)%push(1,'beta')
      call mapTriv2(2)%init(0,0) ; call mapTriv2(2)%push(2,'alpha') ; call mapTriv2(2)%push(2,'beta')
!
!     Non-square (2x3) grid of 2x2 square matrix-kind cells, and a
!     compatible (3x2) grid, for the same-shape-required operations
!     (Contraction, FullContraction, MatMul-outer, PartialContraction
!     default index1). Cell(I,J)/(J,K) values are distinct everywhere.
!
      call NEFSetA%init(2,3)
      do I = 1, 2
        do J = 1, 3
          base = real(10*((I-1)*3+J))
          call cellInt%init(2,2)
          sc = base+1.0 ; call cellInt%put(sc,1,1)
          sc = base+2.0 ; call cellInt%put(sc,1,2)
          sc = base+3.0 ; call cellInt%put(sc,2,1)
          sc = base+4.0 ; call cellInt%put(sc,2,2)
          call cellMat%init(map2)
          call cellMat%imp(cellInt)
          call NEFSetA%put(cellMat,I,J)
        end do
      end do
!
      call NEFSetB%init(3,2)
      do J = 1, 3
        do K = 1, 2
          base = real(1000+100*((J-1)*2+K))
          call cellInt%init(2,2)
          sc = base+1.0 ; call cellInt%put(sc,1,1)
          sc = base+2.0 ; call cellInt%put(sc,1,2)
          sc = base+3.0 ; call cellInt%put(sc,2,1)
          sc = base+4.0 ; call cellInt%put(sc,2,2)
          call cellMat%init(map2)
          call cellMat%imp(cellInt)
          call NEFSetB%put(cellMat,J,K)
        end do
      end do
!
      call NEFSetC%init(2,3)
      do I = 1, 2
        do J = 1, 3
          base = real(5000+100*((I-1)*3+J))
          call cellInt%init(2,2)
          sc = base+1.0 ; call cellInt%put(sc,1,1)
          sc = base+2.0 ; call cellInt%put(sc,1,2)
          sc = base+3.0 ; call cellInt%put(sc,2,1)
          sc = base+4.0 ; call cellInt%put(sc,2,2)
          call cellMat%init(map2)
          call cellMat%imp(cellInt)
          call NEFSetC%put(cellMat,I,J)
        end do
      end do
!
!     Vector-kind (length-2) analogs, same grid shapes.
!
      call NEFSetAv%init(2,3)
      do I = 1, 2
        do J = 1, 3
          base = real(10*((I-1)*3+J))
          call plainVec2%init(2)
          sc = base+1.0 ; call plainVec2%put(sc,1)
          sc = base+2.0 ; call plainVec2%put(sc,2)
          call cellVec%init(map2)
          call cellVec%impv(plainVec2)
          call NEFSetAv%put(cellVec,I,J)
        end do
      end do
!
      call NEFSetBv%init(3,2)
      do J = 1, 3
        do K = 1, 2
          base = real(1000+100*((J-1)*2+K))
          call plainVec2%init(2)
          sc = base+1.0 ; call plainVec2%put(sc,1)
          sc = base+2.0 ; call plainVec2%put(sc,2)
          call cellVec%init(map2)
          call cellVec%impv(plainVec2)
          call NEFSetBv%put(cellVec,J,K)
        end do
      end do
!
      call NEFSetCv%init(2,3)
      do I = 1, 2
        do J = 1, 3
          base = real(5000+100*((I-1)*3+J))
          call plainVec2%init(2)
          sc = base+1.0 ; call plainVec2%put(sc,1)
          sc = base+2.0 ; call plainVec2%put(sc,2)
          call cellVec%init(map2)
          call cellVec%impv(plainVec2)
          call NEFSetCv%put(cellVec,I,J)
        end do
      end do
!
!     R4-kind analogs, built by expanding the same 2x2 matrix cells
!     via a trivial(1)+2-space axis split, same grid shapes.
!
      call NEFSetAr4%init(2,3)
      do I = 1, 2
        do J = 1, 3
          cellMat = NEFSetA%atMatrix(I,J)
          cellR4 = cellMat%expand(mapTriv1,mapTriv2,mapTriv1,mapTriv2,spinType='singlet')
          call NEFSetAr4%put(cellR4,I,J)
        end do
      end do
!
      call NEFSetBr4%init(3,2)
      do J = 1, 3
        do K = 1, 2
          cellMat = NEFSetB%atMatrix(J,K)
          cellR4 = cellMat%expand(mapTriv1,mapTriv2,mapTriv1,mapTriv2,spinType='singlet')
          call NEFSetBr4%put(cellR4,J,K)
        end do
      end do
!
      call NEFSetCr4%init(2,3)
      do I = 1, 2
        do J = 1, 3
          cellMat = NEFSetC%atMatrix(I,J)
          cellR4 = cellMat%expand(mapTriv1,mapTriv2,mapTriv1,mapTriv2,spinType='singlet')
          call NEFSetCr4%put(cellR4,I,J)
        end do
      end do
!
!     Non-square-CELL grids (2x3 / 3x2 inner matrices, via NEFMatR/
!     NEFMatC's row/col maps) for the mixed-kind PartialContraction
!     tests, where non-square inner structure is the actual point.
!
      call NEFSetR%init(2,2)
      do I = 1, 2
        do J = 1, 2
          base = real(10*((I-1)*2+J))
          call fullInt3%init(3,3)
          do ii = 1, 3
            do jj = 1, 3
              sc = base+real((ii-1)*3+jj)
              call fullInt3%put(sc,ii,jj)
            end do
          end do
          call cellMat%init(map2,initial_map_matrix_col=map3)
          call cellMat%imp(fullInt3)
          call NEFSetR%put(cellMat,I,J)
        end do
      end do
!
      call NEFSetS%init(2,2)
      do J = 1, 2
        do K = 1, 2
          base = real(500+10*((J-1)*2+K))
          call fullInt3%init(3,3)
          do ii = 1, 3
            do jj = 1, 3
              sc = base+real((ii-1)*3+jj)
              call fullInt3%put(sc,ii,jj)
            end do
          end do
          call cellMat%init(map3,initial_map_matrix_col=map2)
          call cellMat%imp(fullInt3)
          call NEFSetS%put(cellMat,J,K)
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 1: MQC_NEFSet_Vector_Multiply. plainVec3=[1,0,0] picks out
!     J=1 only, so NEFSetOut(I,1) should equal NEFSetA(I,1).
!--------------------------------------------------------------------
!
      call plainVec3%init(3)
      sc = 1.0 ; call plainVec3%put(sc,1)
      sc = 0.0 ; call plainVec3%put(sc,2)
      sc = 0.0 ; call plainVec3%put(sc,3)
!
      NEFSetOut = MatMul(NEFSetA,plainVec3)
      do I = 1, 2
        gotMat = NEFSetOut%atMatrix(I,1)
        expMat = NEFSetA%atMatrix(I,1)
        do ii = 1, 2
          do jj = 1, 2
            sc = gotMat%at(1,1,ii,jj,'alpha') ; got = sc%rval()
            sc = expMat%at(1,1,ii,jj,'alpha') ; expected = sc%rval()
            write(testLabel,'(a,i0,a,i0,a,i0,a)') &
              'T1 NEFSet*vector: I=',I,' (',ii,',',jj,')'
            if(abs(got-expected).le.tol) then
              write(*,1000) trim(testLabel)
            else
              write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
            endIf
          end do
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 2: MQC_Vector_NEFSet_Multiply. plainVec2=[1,0] picks out
!     I=1 only, so NEFSetOut(1,J) should equal NEFSetA(1,J).
!--------------------------------------------------------------------
!
      call plainVec2%init(2)
      sc = 1.0 ; call plainVec2%put(sc,1)
      sc = 0.0 ; call plainVec2%put(sc,2)
!
      NEFSetOut = MatMul(plainVec2,NEFSetA)
      do J = 1, 3
        gotMat = NEFSetOut%atMatrix(1,J)
        expMat = NEFSetA%atMatrix(1,J)
        sc = gotMat%at(1,1,1,1,'alpha') ; got = sc%rval()
        sc = expMat%at(1,1,1,1,'alpha') ; expected = sc%rval()
        write(testLabel,'(a,i0)') 'T2 vector*NEFSet: J=',J
        if(abs(got-expected).le.tol) then
          write(*,1000) trim(testLabel)
        else
          write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
        endIf
      end do
!
!--------------------------------------------------------------------
!     Test 3: MQC_NEFSet_Matrix_Multiply. plainMat32's columns are the
!     unit vectors [1,0,0] and [0,1,0], so NEFSetOut(I,1)=NEFSetA(I,1)
!     and NEFSetOut(I,2)=NEFSetA(I,2).
!--------------------------------------------------------------------
!
      call plainMat32%init(3,2)
      sc=1.0 ; call plainMat32%put(sc,1,1)
      sc=0.0 ; call plainMat32%put(sc,2,1)
      sc=0.0 ; call plainMat32%put(sc,3,1)
      sc=0.0 ; call plainMat32%put(sc,1,2)
      sc=1.0 ; call plainMat32%put(sc,2,2)
      sc=0.0 ; call plainMat32%put(sc,3,2)
!
      NEFSetOut = MatMul(NEFSetA,plainMat32)
      do I = 1, 2
        do K = 1, 2
          gotMat = NEFSetOut%atMatrix(I,K)
          expMat = NEFSetA%atMatrix(I,K)
          sc = gotMat%at(1,1,1,1,'alpha') ; got = sc%rval()
          sc = expMat%at(1,1,1,1,'alpha') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T3 NEFSet*matrix: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 4: MQC_Matrix_NEFSet_Multiply. plainMat22=identity, so
!     NEFSetOut should equal NEFSetA unchanged.
!--------------------------------------------------------------------
!
      call plainMat22%init(2,2)
      sc=1.0 ; call plainMat22%put(sc,1,1)
      sc=0.0 ; call plainMat22%put(sc,1,2)
      sc=0.0 ; call plainMat22%put(sc,2,1)
      sc=1.0 ; call plainMat22%put(sc,2,2)
!
      NEFSetOut = MatMul(plainMat22,NEFSetA)
      do I = 1, 2
        do J = 1, 3
          gotMat = NEFSetOut%atMatrix(I,J)
          expMat = NEFSetA%atMatrix(I,J)
          sc = gotMat%at(1,1,1,1,'alpha') ; got = sc%rval()
          sc = expMat%at(1,1,1,1,'alpha') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T4 matrix*NEFSet: I=',I,' J=',J
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 5/6: MQC_NEFSet_NEFMatrix_Multiply / MQC_NEFMatrix_NEFSet_
!     Multiply (per-cell broadcast of a shared NEFMatrix against every
!     (I,J) cell -- NOT a grid reduction), both against a 2x2 identity
!     NEFMatrix, verified cell-by-cell against the already-trusted
!     base-level MatMul.
!--------------------------------------------------------------------
!
      call cellInt%init(2,2)
      sc=1.0 ; call cellInt%put(sc,1,1)
      sc=0.0 ; call cellInt%put(sc,1,2)
      sc=0.0 ; call cellInt%put(sc,2,1)
      sc=1.0 ; call cellInt%put(sc,2,2)
      call cellMat%init(map2)
      call cellMat%imp(cellInt)
!
      NEFSetOut = MatMul(NEFSetA,cellMat)
      do I = 1, 2
        do J = 1, 3
          expMat = MatMul(NEFSetA%atMatrix(I,J),cellMat)
          gotMat = NEFSetOut%atMatrix(I,J)
          do ii = 1, 2
            do jj = 1, 2
              sc = gotMat%at(1,1,ii,jj,'alpha') ; got = sc%rval()
              sc = expMat%at(1,1,ii,jj,'alpha') ; expected = sc%rval()
              write(testLabel,'(a,i0,a,i0,a,i0,a,i0,a)') &
                   'T5 NEFSet*NEFMatrix: I=',I,' J=',J,' (',ii,',',jj,')'
              if(abs(got-expected).le.tol) then
                write(*,1000) trim(testLabel)
              else
                write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
              endIf
            end do
          end do
        end do
      end do
!
      NEFSetOut = MatMul(cellMat,NEFSetA)
      do I = 1, 2
        do J = 1, 3
          expMat = MatMul(cellMat,NEFSetA%atMatrix(I,J))
          gotMat = NEFSetOut%atMatrix(I,J)
          do ii = 1, 2
            do jj = 1, 2
              sc = gotMat%at(1,1,ii,jj,'alpha') ; got = sc%rval()
              sc = expMat%at(1,1,ii,jj,'alpha') ; expected = sc%rval()
              write(testLabel,'(a,i0,a,i0,a,i0,a,i0,a)') &
                   'T6 NEFMatrix*NEFSet: I=',I,' J=',J,' (',ii,',',jj,')'
              if(abs(got-expected).le.tol) then
                write(*,1000) trim(testLabel)
              else
                write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
              endIf
            end do
          end do
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 7: MQC_NEFSet_NEFSet_Multiply, vector-kind (outer2Matrix
!     inner op). Verified against a manual Sum_J outer2Matrix.
!--------------------------------------------------------------------
!
      NEFSetOut = MatMul(NEFSetAv,NEFSetBv)
      do I = 1, 2
        do K = 1, 2
          expVec = NEFSetAv%atVector(I,1)
          cellVec = NEFSetBv%atVector(1,K)
          expMat = mqc_nefvector_nefvector_outer_product(expVec,cellVec)
          do J = 2, 3
            expVec = NEFSetAv%atVector(I,J)
            cellVec = NEFSetBv%atVector(J,K)
            expMat = expMat + mqc_nefvector_nefvector_outer_product(expVec,cellVec)
          end do
          gotMat = NEFSetOut%atMatrix(I,K)
          sc = gotMat%at(1,1,1,1,'alpha') ; got = sc%rval()
          sc = expMat%at(1,1,1,1,'alpha') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T7 NEFSet*NEFSet vector-outer: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 8: MQC_NEFSet_NEFSet_Multiply, matrix-kind (outer2R4 inner
!     op). Verified against a manual Sum_J outer2R4.
!--------------------------------------------------------------------
!
      NEFSetOut = MatMul(NEFSetA,NEFSetB)
      do I = 1, 2
        do K = 1, 2
          expMat = NEFSetA%atMatrix(I,1)
          cellMat = NEFSetB%atMatrix(1,K)
          expR4 = mqc_nefmatrix_nefmatrix_outer_product(expMat,cellMat)
          do J = 2, 3
            expMat = NEFSetA%atMatrix(I,J)
            cellMat = NEFSetB%atMatrix(J,K)
            expR4 = expR4 + mqc_nefmatrix_nefmatrix_outer_product(expMat,cellMat)
          end do
          gotR4 = NEFSetOut%atR4(I,K)
          sc = gotR4%at([1,1,1,1],[2,1,1,1],&
                        [3,1,1,1],[4,1,1,1])
          got = sc%rval()
          sc = expR4%at([1,1,1,1],[2,1,1,1],&
                        [3,1,1,1],[4,1,1,1])
          expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T8 NEFSet*NEFSet matrix-outer: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 9: MQC_NEFVector_NEFVector_Contraction (base-level dot
!     product). map2 pushes the SAME orbital index for 'alpha' and
!     'beta' in each space, so impv duplicates each input value into
!     both spin slots: [1,2] -> (alpha1=1,beta1=1,alpha2=2,beta2=2),
!     [3,4] -> (3,3,4,4). Hand-computed dot = 1*3+1*3+2*4+2*4 = 22.
!--------------------------------------------------------------------
!
      call plainVec2%init(2)
      sc=1.0 ; call plainVec2%put(sc,1)
      sc=2.0 ; call plainVec2%put(sc,2)
      call cellVec%init(map2)
      call cellVec%impv(plainVec2)
!
      call plainVec2%init(2)
      sc=3.0 ; call plainVec2%put(sc,1)
      sc=4.0 ; call plainVec2%put(sc,2)
      call expVec%init(map2)
      call expVec%impv(plainVec2)
!
      sc = Contraction(cellVec,expVec)
      got = sc%rval() ; expected = 22.0d0
      if(abs(got-expected).le.tol) then
        write(*,1000) 'T9 bare Vector,Vector Contraction'
      else
        write(*,1010) 'T9 bare Vector,Vector Contraction',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 10/11/12: MQC_NEFSet_NEFSet_Contraction, vector/matrix/r4
!     kinds, each verified against a manual Sum_J Contraction(cell,cell).
!--------------------------------------------------------------------
!
      outMat = Contraction(NEFSetAv,NEFSetBv)
      do I = 1, 2
        do K = 1, 2
          outScalar = Contraction(NEFSetAv%atVector(I,1),NEFSetBv%atVector(1,K))
          do J = 2, 3
            outScalar = outScalar + Contraction(NEFSetAv%atVector(I,J),NEFSetBv%atVector(J,K))
          end do
          sc = outMat%at(I,K) ; got = sc%rval()
          expected = outScalar%rval()
          write(testLabel,'(a,i0,a,i0)') 'T10 NEFSet Contraction vector: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
      outMat = Contraction(NEFSetA,NEFSetB)
      do I = 1, 2
        do K = 1, 2
          outScalar = Contraction(NEFSetA%atMatrix(I,1),NEFSetB%atMatrix(1,K))
          do J = 2, 3
            outScalar = outScalar + Contraction(NEFSetA%atMatrix(I,J),NEFSetB%atMatrix(J,K))
          end do
          sc = outMat%at(I,K) ; got = sc%rval()
          expected = outScalar%rval()
          write(testLabel,'(a,i0,a,i0)') 'T11 NEFSet Contraction matrix: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
      outMat = Contraction(NEFSetAr4,NEFSetBr4)
      do I = 1, 2
        do K = 1, 2
          outScalar = Contraction(NEFSetAr4%atR4(I,1),NEFSetBr4%atR4(1,K),'coulomb')
          do J = 2, 3
            outScalar = outScalar + Contraction(NEFSetAr4%atR4(I,J),NEFSetBr4%atR4(J,K),'coulomb')
          end do
          sc = outMat%at(I,K) ; got = sc%rval()
          expected = outScalar%rval()
          write(testLabel,'(a,i0,a,i0)') 'T12 NEFSet Contraction r4: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 13/14/15: MQC_NEFSet_NEFSet_Full_Contraction, vector/
!     matrix/r4 kinds (NEFSetA/Av/Ar4 against NEFSetC/Cv/Cr4, same
!     shape), verified against a manual Sum_I Sum_J Contraction.
!--------------------------------------------------------------------
!
      expScalar = Contraction(NEFSetAv%atVector(1,1),NEFSetCv%atVector(1,1))
      do J = 2, 3
        expScalar = expScalar + Contraction(NEFSetAv%atVector(1,J),NEFSetCv%atVector(1,J))
      end do
      do I = 2, 2
        do J = 1, 3
          expScalar = expScalar + Contraction(NEFSetAv%atVector(I,J),NEFSetCv%atVector(I,J))
        end do
      end do
      outScalar = FullContraction(NEFSetAv,NEFSetCv)
      got = outScalar%rval() ; expected = expScalar%rval()
      if(abs(got-expected).le.tol) then
        write(*,1000) 'T13 NEFSet FullContraction vector'
      else
        write(*,1010) 'T13 NEFSet FullContraction vector',got,expected ; nfail=nfail+1
      endIf
!
      expScalar = Contraction(NEFSetA%atMatrix(1,1),NEFSetC%atMatrix(1,1))
      do J = 2, 3
        expScalar = expScalar + Contraction(NEFSetA%atMatrix(1,J),NEFSetC%atMatrix(1,J))
      end do
      do I = 2, 2
        do J = 1, 3
          expScalar = expScalar + Contraction(NEFSetA%atMatrix(I,J),NEFSetC%atMatrix(I,J))
        end do
      end do
      outScalar = FullContraction(NEFSetA,NEFSetC)
      got = outScalar%rval() ; expected = expScalar%rval()
      if(abs(got-expected).le.tol) then
        write(*,1000) 'T14 NEFSet FullContraction matrix'
      else
        write(*,1010) 'T14 NEFSet FullContraction matrix',got,expected ; nfail=nfail+1
      endIf
!
      expScalar = Contraction(NEFSetAr4%atR4(1,1),NEFSetCr4%atR4(1,1),'coulomb')
      do J = 2, 3
        expScalar = expScalar + Contraction(NEFSetAr4%atR4(1,J),NEFSetCr4%atR4(1,J),'coulomb')
      end do
      do I = 2, 2
        do J = 1, 3
          expScalar = expScalar + Contraction(NEFSetAr4%atR4(I,J),NEFSetCr4%atR4(I,J),'coulomb')
        end do
      end do
      outScalar = FullContraction(NEFSetAr4,NEFSetCr4)
      got = outScalar%rval() ; expected = expScalar%rval()
      if(abs(got-expected).le.tol) then
        write(*,1000) 'T15 NEFSet FullContraction r4'
      else
        write(*,1010) 'T15 NEFSet FullContraction r4',got,expected ; nfail=nfail+1
      endIf
!
!--------------------------------------------------------------------
!     Test 16: PartialContraction matrix,matrix with DEFAULT index1
!     (=[2,3,1,4]) -- must equal ordinary MatMul(cellA,cellB) summed
!     over J, since the default reduces exactly to that.
!--------------------------------------------------------------------
!
      NEFSetOut = PartialContraction(NEFSetA,NEFSetB)
      do I = 1, 2
        do K = 1, 2
          expMat = MatMul(NEFSetA%atMatrix(I,1),NEFSetB%atMatrix(1,K))
          do J = 2, 3
            expMat = expMat + MatMul(NEFSetA%atMatrix(I,J),NEFSetB%atMatrix(J,K))
          end do
          gotMat = NEFSetOut%atMatrix(I,K)
          sc = gotMat%at(1,1,1,1,'alpha') ; got = sc%rval()
          sc = expMat%at(1,1,1,1,'alpha') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T16 PartialContraction matrix,matrix default: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 17: PartialContraction matrix,matrix with NON-DEFAULT,
!     NON-SQUARE-cell index1=[1,4,3,2]. NEFSetR cells are 2x3 (map2
!     row x map3 col); NEFSetS cells are 3x2 (map3 row x map2 col).
!     a=1 (MatA's row, size 2) contracted against b=4 (MatB's col,
!     size 2) -- the only non-default axis pair whose sizes actually
!     conform here. Survivors: MatA's col (size 3, global 2) and
!     MatB's row (size 3, global 3); c=3,d=2 requests the SWAPPED
!     output order (MatB's survivor first). Verified against the
!     equivalent Transpose/MatMul expression directly, since that IS
!     what the combined-index convention reduces to by construction.
!--------------------------------------------------------------------
!
      NEFSetOut = PartialContraction(NEFSetR,NEFSetS,index1=[1,4,3,2])
      do I = 1, 2
        do K = 1, 2
          expMat = MatMul(Transpose(NEFSetR%atMatrix(I,1)),Transpose(NEFSetS%atMatrix(1,K)))
          do J = 2, 2
            gotMat = MatMul(Transpose(NEFSetR%atMatrix(I,J)),Transpose(NEFSetS%atMatrix(J,K)))
            expMat = expMat + gotMat
          end do
          expMat = Transpose(expMat)
          gotMat = NEFSetOut%atMatrix(I,K)
          sc = gotMat%at(1,1,1,1,'alpha') ; got = sc%rval()
          sc = expMat%at(1,1,1,1,'alpha') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T17 PartialContraction matrix,matrix custom: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 18/19: PartialContraction matrix,vector and vector,matrix,
!     default index1, verified against ordinary MatMul summed over J.
!--------------------------------------------------------------------
!
      NEFSetOut = PartialContraction(NEFSetA,NEFSetBv)
      do I = 1, 2
        do K = 1, 2
          expVec = MatMul(NEFSetA%atMatrix(I,1),NEFSetBv%atVector(1,K))
          do J = 2, 3
            expVec = expVec + MatMul(NEFSetA%atMatrix(I,J),NEFSetBv%atVector(J,K))
          end do
          gotVec = NEFSetOut%atVector(I,K)
          plainVec2 = gotVec%getBlockVec() ; sc = plainVec2%at(1) ; got = sc%rval()
          plainVec2 = expVec%getBlockVec() ; sc = plainVec2%at(1) ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T18 PartialContraction matrix,vector: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
      NEFSetOut = PartialContraction(NEFSetAv,NEFSetB)
      do I = 1, 2
        do K = 1, 2
          expVec = MatMul(NEFSetB%atMatrix(1,K),NEFSetAv%atVector(I,1))
          do J = 2, 3
            expVec = expVec + MatMul(NEFSetB%atMatrix(J,K),NEFSetAv%atVector(I,J))
          end do
          gotVec = NEFSetOut%atVector(I,K)
          plainVec2 = gotVec%getBlockVec() ; sc = plainVec2%at(1) ; got = sc%rval()
          plainVec2 = expVec%getBlockVec() ; sc = plainVec2%at(1) ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T19 PartialContraction vector,matrix: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 20: PartialContraction r4,matrix (reuses MQC_NEF_R4_NEF_
!     Matrix_Contraction via index1 as its indexA), verified against a
!     manual Sum_J Contraction(r4cell,matrixcell). NEFSetAr4's cells
!     were built via expand(mapTriv1,mapTriv2,mapTriv1,mapTriv2,...),
!     so axes 1,3 are trivial (size 1) and axes 2,4 carry the real
!     size-2 content -- index1=[2,4,1,3] contracts the real axes
!     (2,4) against NEFSetB's 2x2 row/col, leaving the trivial axes
!     (1,3) as survivors (a degenerate 1x1 output, fine for a plumbing
!     check). The default [1,2,3,4] would contract a trivial size-1
!     axis against NEFSetB's size-2 row and fail to conform. label=
!     'coulomb' avoids the default 'doublebar' label's additional
!     exchange diagram, which internally swaps the ja/la roles and
!     would hit the same trivial-vs-size-2 mismatch the other way.
!--------------------------------------------------------------------
!
      NEFSetOut = PartialContraction(NEFSetAr4,NEFSetB,index1=[2,4,1,3],label='coulomb')
      do I = 1, 2
        do K = 1, 2
          expMat = Contraction(NEFSetAr4%atR4(I,1),NEFSetB%atMatrix(1,K),&
                                label='coulomb',indexA=[2,4,1,3])
          do J = 2, 3
            expMat = expMat + Contraction(NEFSetAr4%atR4(I,J),NEFSetB%atMatrix(J,K),&
                                           label='coulomb',indexA=[2,4,1,3])
          end do
          gotMat = NEFSetOut%atMatrix(I,K)
          sc = gotMat%at(1,1,1,1,'alpha') ; got = sc%rval()
          sc = expMat%at(1,1,1,1,'alpha') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T20 PartialContraction r4,matrix: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 21: PartialContraction r4,r4 -> Matrix-kind NEF_Set,
!     verified against a manual Sum_J PartialContraction(r4,r4).
!--------------------------------------------------------------------
!
      NEFSetOut = PartialContraction(NEFSetAr4,NEFSetBr4)
      do I = 1, 2
        do K = 1, 2
          expMat = PartialContraction(NEFSetAr4%atR4(I,1),NEFSetBr4%atR4(1,K))
          do J = 2, 3
            expMat = expMat + PartialContraction(NEFSetAr4%atR4(I,J),NEFSetBr4%atR4(J,K))
          end do
          gotMat = NEFSetOut%atMatrix(I,K)
          sc = gotMat%at(1,1,1,1,'alpha') ; got = sc%rval()
          sc = expMat%at(1,1,1,1,'alpha') ; expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T21 PartialContraction r4,r4->Matrix: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
!--------------------------------------------------------------------
!     Test 22: PartialContraction2R4 r4,r4 -> R4-kind NEF_Set,
!     verified against a manual Sum_J PartialContraction2R4(r4,r4).
!--------------------------------------------------------------------
!
      NEFSetOut = PartialContraction2R4(NEFSetAr4,NEFSetBr4)
      do I = 1, 2
        do K = 1, 2
          expR4 = PartialContraction2R4(NEFSetAr4%atR4(I,1),NEFSetBr4%atR4(1,K))
          do J = 2, 3
            expR4 = expR4 + PartialContraction2R4(NEFSetAr4%atR4(I,J),NEFSetBr4%atR4(J,K))
          end do
          gotR4 = NEFSetOut%atR4(I,K)
          sc = gotR4%at([1,1,1,1],[2,1,1,1],&
                        [3,1,1,1],[4,1,1,1])
          got = sc%rval()
          sc = expR4%at([1,1,1,1],[2,1,1,1],&
                        [3,1,1,1],[4,1,1,1])
          expected = sc%rval()
          write(testLabel,'(a,i0,a,i0)') 'T22 PartialContraction2R4 r4,r4->R4: I=',I,' K=',K
          if(abs(got-expected).le.tol) then
            write(*,1000) trim(testLabel)
          else
            write(*,1010) trim(testLabel),got,expected ; nfail=nfail+1
          endIf
        end do
      end do
!
      write(*,9000) nfail
      if(nfail.gt.0) stop 1
!
      end program unitTest_NEF_set_algebra
