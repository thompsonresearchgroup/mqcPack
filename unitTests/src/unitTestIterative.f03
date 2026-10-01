Program unitTestIterative
  Use iso_fortran_env,Only:int64,real64
  Use MQC_General
  Use MQC_Algebra
  Implicit None

  Integer(Kind=int64),Parameter::Dimension=24_int64,NumberOfRoots=4_int64
  Integer(Kind=int64),Parameter::MetricRank=20_int64
  Type(MQC_Matrix)::RealMatrix,ComplexMatrix,MetricMatrix,DenseHamiltonian,&
    DenseMetric,DenseVectors,DavidsonVectors,SingularHamiltonian,SingularMetric,&
    ReducedHamiltonian,ReducedMetric
  Type(MQC_Vector)::DenseValues,DavidsonValues,ResidualNorms
  Type(MQC_Explicit_Matrix_Operator)::Operator,MetricOperator
  Real(Kind=real64)::Value,ErrorValue,MaxEigenvalueError,MaxResidual
  Real(Kind=real64)::MaxMetricError
  Complex(Kind=real64)::ComplexValue
  Complex(Kind=real64),Allocatable::RootVectors(:,:),MetricArray(:,:),RootMetric(:,:)
  Real(Kind=real64),Allocatable::Rotation(:,:),HouseholderVector(:),BaseHamiltonian(:,:),&
    BaseMetric(:,:),RotatedHamiltonian(:,:),RotatedMetric(:,:)
  Integer(Kind=int64)::i,j,Iterations
  Logical::Converged

  Call RealMatrix%Init(Dimension,Dimension)
  Do j=1,Dimension
    Do i=j,Dimension
      Value=Sin(Real(3*i+5*j,Kind=real64))+&
        0.1_real64*Cos(Real(7*i-2*j,Kind=real64))
      If(i.eq.j) Value=Value+Real(i,Kind=real64)
      Call RealMatrix%Put(Value,i,j)
      Call RealMatrix%Put(Value,j,i)
    EndDo
  EndDo
  Write(*,'(1x,A)') 'Iterative test: dense real diagonalization'
  Call RealMatrix%Diag(DenseValues,DenseVectors)
  Write(*,'(1x,A)') 'Iterative test: real Davidson'
  Call Operator%Init(RealMatrix)
  Call MQC_Davidson(Operator,NumberOfRoots,DavidsonValues,DavidsonVectors,&
    Tolerance=1.0e-10_real64,MaximumSubspace=12_int64,&
    MaximumIterations=100_int64,Converged=Converged,Iterations=Iterations,&
    ResidualNorms=ResidualNorms)
  Call AssertTrue(Converged,'Real symmetric Davidson calculation did not converge')
  MaxEigenvalueError=0.0_real64
  Do i=1,NumberOfRoots
    ErrorValue=DavidsonValues%At(i)-DenseValues%At(i)
    MaxEigenvalueError=Max(MaxEigenvalueError,Abs(ErrorValue))
  EndDo
  MaxResidual=ResidualNorms%MaxVal()
  Call AssertTrue(MaxEigenvalueError.lt.1.0e-9_real64,&
    'Real symmetric Davidson eigenvalues disagree with dense diagonalization')
  Call AssertTrue(MaxResidual.lt.1.0e-9_real64,&
    'Real symmetric Davidson residual is too large')

  Call ComplexMatrix%Init(Dimension,Dimension,&
    Scalar=(0.0_real64,0.0_real64))
  Do j=1,Dimension
    Do i=j,Dimension
      If(i.eq.j) Then
        ComplexValue=Cmplx(Real(i,Kind=real64)+&
          0.2_real64*Sin(Real(i,Kind=real64)),0.0_real64,Kind=real64)
      Else
        ComplexValue=Cmplx(0.2_real64*Sin(Real(2*i+j,Kind=real64)),&
          0.15_real64*Cos(Real(i-3*j,Kind=real64)),Kind=real64)
      EndIf
      Call ComplexMatrix%Put(ComplexValue,i,j)
      Call ComplexMatrix%Put(Conjg(ComplexValue),j,i)
    EndDo
  EndDo
  Write(*,'(1x,A)') 'Iterative test: dense complex diagonalization'
  Call ComplexMatrix%Diag(DenseValues,DenseVectors)
  Write(*,'(1x,A)') 'Iterative test: complex Davidson'
  Call Operator%Init(ComplexMatrix)
  Call MQC_Davidson(Operator,NumberOfRoots,DavidsonValues,DavidsonVectors,&
    Tolerance=1.0e-10_real64,MaximumSubspace=12_int64,&
    MaximumIterations=100_int64,Converged=Converged,Iterations=Iterations,&
    ResidualNorms=ResidualNorms)
  Call AssertTrue(Converged,'Complex Hermitian Davidson calculation did not converge')
  MaxEigenvalueError=0.0_real64
  Do i=1,NumberOfRoots
    ErrorValue=DavidsonValues%At(i)-DenseValues%At(i)
    MaxEigenvalueError=Max(MaxEigenvalueError,Abs(ErrorValue))
  EndDo
  MaxResidual=ResidualNorms%MaxVal()
  Call AssertTrue(MaxEigenvalueError.lt.1.0e-9_real64,&
    'Complex Hermitian Davidson eigenvalues disagree with dense diagonalization')
  Call AssertTrue(MaxResidual.lt.1.0e-9_real64,&
    'Complex Hermitian Davidson residual is too large')

  Call MetricMatrix%Init(Dimension,Dimension)
  Do j=1,Dimension
    Do i=j,Dimension
      If(i.eq.j) Then
        Value=1.0_real64+0.01_real64*Real(i,Kind=real64)
      Else
        Value=0.08_real64/Real(1+Abs(i-j),Kind=real64)
      EndIf
      Call MetricMatrix%Put(Value,i,j)
      Call MetricMatrix%Put(Value,j,i)
    EndDo
  EndDo
  DenseHamiltonian=RealMatrix
  DenseMetric=MetricMatrix
  Write(*,'(1x,A)') 'Iterative test: dense generalized diagonalization'
  Call DenseHamiltonian%Eigensys(DenseMetric,DenseValues,DenseVectors)
  Call Operator%Init(RealMatrix)
  Call MetricOperator%Init(MetricMatrix)
  Write(*,'(1x,A)') 'Iterative test: generalized Davidson'
  Call MQC_Davidson(Operator,NumberOfRoots,DavidsonValues,DavidsonVectors,&
    MetricOperator=MetricOperator,Tolerance=1.0e-10_real64,&
    MaximumSubspace=12_int64,MaximumIterations=150_int64,&
    Converged=Converged,Iterations=Iterations,ResidualNorms=ResidualNorms)
  Call AssertTrue(Converged,'Generalized Davidson calculation did not converge')
  MaxEigenvalueError=0.0_real64
  Do i=1,NumberOfRoots
    ErrorValue=DavidsonValues%At(i)-DenseValues%At(i)
    MaxEigenvalueError=Max(MaxEigenvalueError,Abs(ErrorValue))
  EndDo
  MaxResidual=ResidualNorms%MaxVal()
  Call AssertTrue(MaxEigenvalueError.lt.1.0e-9_real64,&
    'Generalized Davidson eigenvalues disagree with dense diagonalization')
  Call AssertTrue(MaxResidual.lt.1.0e-9_real64,&
    'Generalized Davidson residual is too large')
  Allocate(RootVectors(Dimension,NumberOfRoots),MetricArray(Dimension,Dimension),&
    RootMetric(NumberOfRoots,NumberOfRoots))
  Do j=1,Dimension
    Do i=1,Dimension
      MetricArray(i,j)=MetricMatrix%At(i,j)
    EndDo
  EndDo
  Do j=1,NumberOfRoots
    Do i=1,Dimension
      RootVectors(i,j)=DavidsonVectors%At(i,j)
    EndDo
  EndDo
  RootMetric=MatMul(Conjg(Transpose(RootVectors)),MatMul(MetricArray,RootVectors))
  MaxMetricError=0.0_real64
  Do j=1,NumberOfRoots
    Do i=1,NumberOfRoots
      If(i.eq.j) Then
        MaxMetricError=Max(MaxMetricError,Abs(RootMetric(i,j)-1.0_real64))
      Else
        MaxMetricError=Max(MaxMetricError,Abs(RootMetric(i,j)))
      EndIf
    EndDo
  EndDo
  Call AssertTrue(MaxMetricError.lt.1.0e-10_real64,&
    'Generalized Davidson eigenvectors are not metric orthonormal')

  Allocate(Rotation(Dimension,Dimension),HouseholderVector(Dimension),&
    BaseHamiltonian(Dimension,Dimension),BaseMetric(Dimension,Dimension),&
    RotatedHamiltonian(Dimension,Dimension),RotatedMetric(Dimension,Dimension))
  Rotation=0.0_real64
  BaseHamiltonian=0.0_real64
  BaseMetric=0.0_real64
  Do i=1,Dimension
    Rotation(i,i)=1.0_real64
    HouseholderVector(i)=Sin(Real(i,Kind=real64))
  EndDo
  Rotation=Rotation-2.0_real64*MatMul(Reshape(HouseholderVector,[Dimension,1]),&
    Reshape(HouseholderVector,[1,Dimension]))/Dot_Product(HouseholderVector,HouseholderVector)
  Call ReducedHamiltonian%Init(MetricRank,MetricRank)
  Call ReducedMetric%Init(MetricRank,MetricRank)
  Do j=1,MetricRank
    Do i=1,MetricRank
      BaseHamiltonian(i,j)=RealMatrix%At(i,j)
      BaseMetric(i,j)=MetricMatrix%At(i,j)
      Call ReducedHamiltonian%Put(BaseHamiltonian(i,j),i,j)
      Call ReducedMetric%Put(BaseMetric(i,j),i,j)
    EndDo
  EndDo
  Do i=MetricRank+1,Dimension
    BaseHamiltonian(i,i)=100.0_real64+Real(i,Kind=real64)
  EndDo
  RotatedHamiltonian=MatMul(Transpose(Rotation),MatMul(BaseHamiltonian,Rotation))
  RotatedMetric=MatMul(Transpose(Rotation),MatMul(BaseMetric,Rotation))
  Call SingularHamiltonian%Init(Dimension,Dimension)
  Call SingularMetric%Init(Dimension,Dimension)
  Do j=1,Dimension
    Do i=1,Dimension
      Call SingularHamiltonian%Put(RotatedHamiltonian(i,j),i,j)
      Call SingularMetric%Put(RotatedMetric(i,j),i,j)
    EndDo
  EndDo
  Call ReducedHamiltonian%Eigensys(ReducedMetric,DenseValues,DenseVectors)
  Call Operator%Init(SingularHamiltonian)
  Call MetricOperator%Init(SingularMetric)
  Write(*,'(1x,A)') 'Iterative test: rank-deficient generalized Davidson'
  Call MQC_Davidson(Operator,NumberOfRoots,DavidsonValues,DavidsonVectors,&
    MetricOperator=MetricOperator,Tolerance=1.0e-10_real64,&
    MaximumSubspace=16_int64,MaximumIterations=200_int64,&
    Converged=Converged,Iterations=Iterations,ResidualNorms=ResidualNorms)
  Call AssertTrue(Converged,'Rank-deficient Davidson calculation did not converge')
  MaxEigenvalueError=0.0_real64
  Do i=1,NumberOfRoots
    ErrorValue=DavidsonValues%At(i)-DenseValues%At(i)
    MaxEigenvalueError=Max(MaxEigenvalueError,Abs(ErrorValue))
  EndDo
  MaxResidual=ResidualNorms%MaxVal()
  Call AssertTrue(MaxEigenvalueError.lt.1.0e-9_real64,&
    'Rank-deficient Davidson eigenvalues disagree with reduced dense problem')
  Call AssertTrue(MaxResidual.lt.1.0e-9_real64,&
    'Rank-deficient Davidson residual is too large')

  Call RealMatrix%Diag(DenseValues,DenseVectors)
  Call Operator%Init(RealMatrix)
  Write(*,'(1x,A)') 'Iterative test: optional Davidson results'
  Call MQC_Davidson(Operator,Eigenvalues=DavidsonValues,&
    Tolerance=1.0e-10_real64,Converged=Converged)
  Call AssertTrue(Converged,'Default one-root Davidson calculation did not converge')
  Call AssertTrue(Size(DavidsonValues).eq.1_int64,&
    'Default Davidson root count is not one')
  Call AssertTrue(Abs(DavidsonValues%At(1)-DenseValues%At(1)).lt.1.0e-9_real64,&
    'Default Davidson root disagrees with dense diagonalization')
  Call MQC_Davidson(Operator,NumberOfRoots=2_int64,Eigenvectors=DavidsonVectors,&
    Tolerance=1.0e-10_real64,Converged=Converged)
  Call AssertTrue(Converged,'Eigenvector-only Davidson calculation did not converge')
  Call AssertTrue(Size(DavidsonVectors,1).eq.Dimension.and.&
    Size(DavidsonVectors,2).eq.2_int64,'Incorrect eigenvector-only result dimensions')
  Call MQC_Davidson(Operator,Tolerance=1.0e-10_real64,Converged=Converged)
  Call AssertTrue(Converged,'Result-free Davidson calculation did not converge')

  Write(*,'(1x,A)') 'Iterative test: matrix-bound Davidson'
  Call RealMatrix%Davidson(DavidsonValues,DavidsonVectors,&
    NumberOfRoots=NumberOfRoots,Tolerance=1.0e-10_real64,Converged=Converged)
  Call AssertTrue(Converged,'Matrix-bound Davidson calculation did not converge')
  MaxEigenvalueError=0.0_real64
  Do i=1,NumberOfRoots
    ErrorValue=DavidsonValues%At(i)-DenseValues%At(i)
    MaxEigenvalueError=Max(MaxEigenvalueError,Abs(ErrorValue))
  EndDo
  Call AssertTrue(MaxEigenvalueError.lt.1.0e-9_real64,&
    'Matrix-bound Davidson eigenvalues disagree with dense diagonalization')

  DenseHamiltonian=RealMatrix
  DenseMetric=MetricMatrix
  Call DenseHamiltonian%Eigensys(DenseMetric,DenseValues,DenseVectors)
  Write(*,'(1x,A)') 'Iterative test: matrix-bound generalized Davidson'
  Call RealMatrix%Davidson(DavidsonValues,DavidsonVectors,&
    NumberOfRoots=NumberOfRoots,MetricMatrix=MetricMatrix,&
    Tolerance=1.0e-10_real64,MaximumSubspace=12_int64,&
    MaximumIterations=150_int64,Converged=Converged)
  Call AssertTrue(Converged,&
    'Matrix-bound generalized Davidson calculation did not converge')
  MaxEigenvalueError=0.0_real64
  Do i=1,NumberOfRoots
    ErrorValue=DavidsonValues%At(i)-DenseValues%At(i)
    MaxEigenvalueError=Max(MaxEigenvalueError,Abs(ErrorValue))
  EndDo
  Call AssertTrue(MaxEigenvalueError.lt.1.0e-9_real64,&
    'Matrix-bound generalized Davidson eigenvalues disagree with dense result')

  Write(*,'(1x,A,I0,A,ES12.5,A,ES12.5)') 'unitTestIterative: PASS (',&
    Iterations,' iterations, eigenvalue error ',MaxEigenvalueError,&
    ', residual ',MaxResidual

Contains

  Subroutine AssertTrue(Condition,Message)
    Logical,Intent(In)::Condition
    Character(Len=*),Intent(In)::Message

    If(.not.Condition) Then
      Write(*,'(1x,A)') 'FAIL: '//Trim(Message)
      Error Stop 1
    EndIf
  End Subroutine AssertTrue

End Program unitTestIterative
