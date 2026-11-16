!> \brief Implementations of the MQC_Algebra matrix-free eigensolver API.
Submodule(MQC_Algebra) MQC_Algebra_Iterative
  Implicit None
Contains

!> Store a dense matrix in the explicit linear-operator adapter.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Explicit_Operator_Init
    If(Size(Matrix,1).ne.Size(Matrix,2)) &
      Call MQC_Error('Explicit linear operator requires a square matrix')
    Operator%Matrix=Matrix
  End Procedure MQC_Explicit_Operator_Init

!> Return the dimension of an explicit linear operator.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Explicit_Operator_Dimension
    Dimension=Size(Operator%Matrix,1)
  End Procedure MQC_Explicit_Operator_Dimension

!> Apply the explicit matrix to a block of trial vectors.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Explicit_Operator_Apply
    Complex(Kind=real64)::Element
    Integer(Kind=int64)::i,j

    If(Size(Vectors,1).ne.Operator%Dimension()) &
      Call MQC_Error('Vector dimension disagrees with explicit operator')
    If(Any(Shape(Products).ne.Shape(Vectors))) &
      Call MQC_Error('Product shape disagrees with input vector block')
    Products=(0.0_real64,0.0_real64)
    Do j=1,Operator%Dimension()
      Do i=1,Operator%Dimension()
        Element=Operator%Matrix%At(i,j)
        Products(i,:)=Products(i,:)+Element*Vectors(j,:)
      EndDo
    EndDo
  End Procedure MQC_Explicit_Operator_Apply

!> Return the diagonal of an explicit matrix operator.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Explicit_Operator_Diagonal
    Type(MQC_Vector)::MatrixDiagonal
    Integer(Kind=int64)::i

    If(Size(Diagonal).ne.Operator%Dimension()) &
      Call MQC_Error('Diagonal dimension disagrees with explicit operator')
    MatrixDiagonal=Operator%Matrix%Diagonal()
    Do i=1,Operator%Dimension()
      Diagonal(i)=MatrixDiagonal%At(i)
    EndDo
  End Procedure MQC_Explicit_Operator_Diagonal

!> Compute the lowest algebraic eigenpairs of a Hermitian linear operator
!> using thick-restarted block Davidson iteration.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Davidson
    Complex(Kind=real64),Allocatable::Basis(:,:),OperatorBasis(:,:),MetricBasis(:,:),&
      ProjectedVectors(:,:),RitzVectors(:,:),OperatorRitz(:,:),Residuals(:,:),&
      Corrections(:,:),Diagonal(:),MetricDiagonal(:),WorkVector(:),MetricRitz(:,:)
    Real(Kind=real64),Allocatable::ProjectedValues(:),RootResiduals(:)
    Type(MQC_Matrix)::ProjectedMatrix,ProjectedMetric,ProjectedEigenvectors
    Type(MQC_Vector)::ProjectedEigenvalues,DiagonalValues,DiagonalOrder
    Complex(Kind=real64)::Denominator
    Real(Kind=real64)::ConvergenceTolerance,MetricCutoffTolerance,NormValue,&
      DenominatorThreshold
    Integer(Kind=int64)::Dimension,MaxSubspace,MaxIterations,SubspaceSize,&
      CorrectionCount,IterationCount,i,j,k,InitialCount,NewCount,ProjectedRank,&
      RequestedRoots
    Logical::AllConverged

    Dimension=Operator%Dimension()
    RequestedRoots=1_int64
    If(Present(NumberOfRoots)) RequestedRoots=NumberOfRoots
    If(Present(MetricOperator)) Then
      If(MetricOperator%Dimension().ne.Dimension) &
        Call MQC_Error('Operator and metric dimensions disagree in MQC_Davidson')
    EndIf
    If(RequestedRoots.lt.1.or.RequestedRoots.gt.Dimension) &
      Call MQC_Error('Invalid number of roots in MQC_Davidson')
    ConvergenceTolerance=1.0e-8_real64
    If(Present(Tolerance)) ConvergenceTolerance=Tolerance
    MetricCutoffTolerance=1.0e-12_real64
    If(Present(MetricTolerance)) MetricCutoffTolerance=MetricTolerance
    MaxSubspace=Min(Dimension,Max(4_int64*RequestedRoots,20_int64))
    If(Present(MaximumSubspace)) MaxSubspace=Min(Dimension,MaximumSubspace)
    MaxSubspace=Min(Dimension,Max(MaxSubspace,RequestedRoots+1_int64))
    MaxIterations=100_int64
    If(Present(MaximumIterations)) MaxIterations=MaximumIterations
    If(ConvergenceTolerance.le.0.0_real64.or.MetricCutoffTolerance.le.0.0_real64.or.&
      MaxIterations.lt.1) &
      Call MQC_Error('Invalid convergence controls in MQC_Davidson')

    Allocate(Basis(Dimension,MaxSubspace),OperatorBasis(Dimension,MaxSubspace),&
      MetricBasis(Dimension,MaxSubspace),MetricRitz(Dimension,RequestedRoots),&
      ProjectedVectors(MaxSubspace,MaxSubspace),RitzVectors(Dimension,RequestedRoots),&
      OperatorRitz(Dimension,RequestedRoots),Residuals(Dimension,RequestedRoots),&
      Corrections(Dimension,RequestedRoots),Diagonal(Dimension),MetricDiagonal(Dimension),&
      WorkVector(Dimension),ProjectedValues(MaxSubspace),&
      RootResiduals(RequestedRoots))
    Basis=(0.0_real64,0.0_real64)
    OperatorBasis=(0.0_real64,0.0_real64)
    MetricBasis=(0.0_real64,0.0_real64)
    Call Operator%Diagonal(Diagonal)
    If(Present(MetricOperator)) Then
      Call MetricOperator%Diagonal(MetricDiagonal)
    Else
      MetricDiagonal=(1.0_real64,0.0_real64)
    EndIf
    Call DiagonalValues%Init(Dimension)
    Do i=1,Dimension
      If(Abs(MetricDiagonal(i)).gt.1.0e-14_real64) Then
        Call DiagonalValues%Put(Real(Diagonal(i)/MetricDiagonal(i),Kind=real64),i)
      Else
        Call DiagonalValues%Put(Huge(1.0_real64),i)
      EndIf
    EndDo
    DiagonalOrder=DiagonalValues%ArgSort()

    InitialCount=RequestedRoots
    If(Present(InitialVectors)) Then
      If(Size(InitialVectors,1).ne.Dimension.or.&
        Size(InitialVectors,2).lt.RequestedRoots) &
        Call MQC_Error('Invalid initial-vector dimensions in MQC_Davidson')
      Do j=1,RequestedRoots
        Do i=1,Dimension
          Basis(i,j)=InitialVectors%At(i,j)
        EndDo
      EndDo
    Else
      Do j=1,RequestedRoots
        i=DiagonalOrder%At(j)
        Basis(i,j)=(1.0_real64,0.0_real64)
      EndDo
    EndIf
    Call MQC_Iterative_Orthonormalize_Block(Basis,InitialCount)
    If(InitialCount.lt.RequestedRoots) &
      Call MQC_Error('Initial Davidson vectors are linearly dependent')
    SubspaceSize=InitialCount
    Call Operator%Apply(Basis(:,1:SubspaceSize),OperatorBasis(:,1:SubspaceSize))
    If(Present(MetricOperator)) Then
      Call MetricOperator%Apply(Basis(:,1:SubspaceSize),MetricBasis(:,1:SubspaceSize))
    Else
      MetricBasis(:,1:SubspaceSize)=Basis(:,1:SubspaceSize)
    EndIf
    AllConverged=.False.

    Do IterationCount=1,MaxIterations
      Call MQC_Iterative_Projected_Matrix(Basis(:,1:SubspaceSize),&
        OperatorBasis(:,1:SubspaceSize),ProjectedMatrix)
      If(Present(MetricOperator)) Then
        Call MQC_Iterative_Projected_Matrix(Basis(:,1:SubspaceSize),&
          MetricBasis(:,1:SubspaceSize),ProjectedMetric)
        Call MQC_Iterative_Generalized_Ritz(ProjectedMatrix,ProjectedMetric,&
          ProjectedEigenvalues,ProjectedEigenvectors,MetricCutoffTolerance,ProjectedRank)
        If(ProjectedRank.lt.RequestedRoots) &
          Call MQC_Error('Davidson metric rank is smaller than the number of requested roots')
      Else
        Call ProjectedMatrix%Diag(ProjectedEigenvalues,ProjectedEigenvectors)
        ProjectedRank=SubspaceSize
      EndIf
      Do i=1,RequestedRoots
        ProjectedValues(i)=ProjectedEigenvalues%At(i)
        Do j=1,SubspaceSize
          ProjectedVectors(j,i)=ProjectedEigenvectors%At(j,i)
        EndDo
      EndDo
      Call MQC_Iterative_Gemm('N',Basis(:,1:SubspaceSize),&
        ProjectedVectors(1:SubspaceSize,1:RequestedRoots),RitzVectors,&
        (1.0_real64,0.0_real64),(0.0_real64,0.0_real64))
      Call MQC_Iterative_Gemm('N',OperatorBasis(:,1:SubspaceSize),&
        ProjectedVectors(1:SubspaceSize,1:RequestedRoots),OperatorRitz,&
        (1.0_real64,0.0_real64),(0.0_real64,0.0_real64))
      Call MQC_Iterative_Gemm('N',MetricBasis(:,1:SubspaceSize),&
        ProjectedVectors(1:SubspaceSize,1:RequestedRoots),MetricRitz,&
        (1.0_real64,0.0_real64),(0.0_real64,0.0_real64))
      CorrectionCount=0
      Do k=1,RequestedRoots
        Residuals(:,k)=OperatorRitz(:,k)-ProjectedValues(k)*MetricRitz(:,k)
        RootResiduals(k)=Sqrt(Real(Dot_Product(Residuals(:,k),Residuals(:,k)),&
          Kind=real64))
        If(RootResiduals(k).le.ConvergenceTolerance) Cycle
        WorkVector=-Residuals(:,k)
        DenominatorThreshold=Max(1.0e-12_real64,&
          1.0e-10_real64*Max(1.0_real64,Abs(ProjectedValues(k))))
        Do i=1,Dimension
          Denominator=Diagonal(i)-ProjectedValues(k)*MetricDiagonal(i)
          If(Abs(Denominator).lt.DenominatorThreshold) Then
            If(Real(Denominator,Kind=real64).lt.0.0_real64) Then
              Denominator=Cmplx(-DenominatorThreshold,0.0_real64,Kind=real64)
            Else
              Denominator=Cmplx(DenominatorThreshold,0.0_real64,Kind=real64)
            EndIf
          EndIf
          WorkVector(i)=WorkVector(i)/Denominator
        EndDo
        CorrectionCount=CorrectionCount+1
        Corrections(:,CorrectionCount)=WorkVector
      EndDo
      If(CorrectionCount.gt.0) Call MQC_Iterative_Block_Orthogonalize(&
        Corrections(:,1:CorrectionCount),CorrectionCount,Basis(:,1:SubspaceSize))
      AllConverged=All(RootResiduals.le.ConvergenceTolerance)
      If(AllConverged.or.CorrectionCount.eq.0) Exit

      If(SubspaceSize+CorrectionCount.gt.MaxSubspace) Then
        Basis=(0.0_real64,0.0_real64)
        Basis(:,1:RequestedRoots)=RitzVectors
        SubspaceSize=RequestedRoots
        Call MQC_Iterative_Orthonormalize_Block(Basis,SubspaceSize)
        NewCount=Min(CorrectionCount,MaxSubspace-SubspaceSize)
        Basis(:,SubspaceSize+1:SubspaceSize+NewCount)=Corrections(:,1:NewCount)
        Call MQC_Iterative_Block_Orthogonalize(&
          Basis(:,SubspaceSize+1:SubspaceSize+NewCount),NewCount,&
          Basis(:,1:SubspaceSize))
        SubspaceSize=SubspaceSize+NewCount
        Call Operator%Apply(Basis(:,1:SubspaceSize),OperatorBasis(:,1:SubspaceSize))
        If(Present(MetricOperator)) Then
          Call MetricOperator%Apply(Basis(:,1:SubspaceSize),MetricBasis(:,1:SubspaceSize))
        Else
          MetricBasis(:,1:SubspaceSize)=Basis(:,1:SubspaceSize)
        EndIf
      Else
        Basis(:,SubspaceSize+1:SubspaceSize+CorrectionCount)=&
          Corrections(:,1:CorrectionCount)
        Call Operator%Apply(Basis(:,SubspaceSize+1:SubspaceSize+CorrectionCount),&
          OperatorBasis(:,SubspaceSize+1:SubspaceSize+CorrectionCount))
        If(Present(MetricOperator)) Then
          Call MetricOperator%Apply(&
            Basis(:,SubspaceSize+1:SubspaceSize+CorrectionCount),&
            MetricBasis(:,SubspaceSize+1:SubspaceSize+CorrectionCount))
        Else
          MetricBasis(:,SubspaceSize+1:SubspaceSize+CorrectionCount)=&
            Basis(:,SubspaceSize+1:SubspaceSize+CorrectionCount)
        EndIf
        SubspaceSize=SubspaceSize+CorrectionCount
      EndIf
    EndDo

    If(Present(Eigenvalues)) Then
      Call Eigenvalues%Init(RequestedRoots)
      Do k=1,RequestedRoots
        Call Eigenvalues%Put(ProjectedValues(k),k)
      EndDo
    EndIf
    If(Present(Eigenvectors)) Then
      Call Eigenvectors%Init(Dimension,RequestedRoots,&
        Scalar=(0.0_real64,0.0_real64))
      Do k=1,RequestedRoots
        Do i=1,Dimension
          Call Eigenvectors%Put(RitzVectors(i,k),i,k)
        EndDo
      EndDo
    EndIf
    If(Present(ResidualNorms)) Then
      Call ResidualNorms%Init(RequestedRoots)
      Do k=1,RequestedRoots
        Call ResidualNorms%Put(RootResiduals(k),k)
      EndDo
    EndIf
    If(Present(Converged)) Converged=AllConverged
    If(Present(Iterations)) Iterations=Min(IterationCount,MaxIterations)
  End Procedure MQC_Davidson

!> Iteratively diagonalize a stored matrix using the explicit linear-operator
!> adapter.  A supplied metric matrix defines the generalized eigenproblem.
!> \author L. M. Thompson
!> \date 2026
  Module Procedure MQC_Matrix_Davidson
    Type(MQC_Explicit_Matrix_Operator)::Operator,MetricOperator

    Call Operator%Init(Matrix)
    If(Present(MetricMatrix)) Then
      Call MetricOperator%Init(MetricMatrix)
      Call MQC_Davidson(Operator,NumberOfRoots,Eigenvalues,Eigenvectors,&
        MetricOperator=MetricOperator,InitialVectors=InitialVectors,Tolerance=Tolerance,&
        MetricTolerance=MetricTolerance,MaximumSubspace=MaximumSubspace,&
        MaximumIterations=MaximumIterations,Converged=Converged,Iterations=Iterations,&
        ResidualNorms=ResidualNorms)
    Else
      Call MQC_Davidson(Operator,NumberOfRoots,Eigenvalues,Eigenvectors,&
        InitialVectors=InitialVectors,Tolerance=Tolerance,MetricTolerance=MetricTolerance,&
        MaximumSubspace=MaximumSubspace,MaximumIterations=MaximumIterations,&
        Converged=Converged,Iterations=Iterations,ResidualNorms=ResidualNorms)
    EndIf
  End Procedure MQC_Matrix_Davidson

!> Solve the Hermitian projected generalized eigenproblem by symmetric metric
!> orthogonalization.  This avoids the non-Hermitian ZGGEV fallback and returns
!> Ritz coefficients that are orthonormal in the projected metric.
!> \author L. M. Thompson
!> \date 2026
  Subroutine MQC_Iterative_Generalized_Ritz(ProjectedMatrix,ProjectedMetric,&
    Eigenvalues,Eigenvectors,MetricTolerance,MetricRank)
    Type(MQC_Matrix),Intent(In)::ProjectedMatrix,ProjectedMetric
    Type(MQC_Vector),Intent(Out)::Eigenvalues
    Type(MQC_Matrix),Intent(Out)::Eigenvectors
    Real(Kind=real64),Intent(In)::MetricTolerance
    Integer(Kind=int64),Intent(Out)::MetricRank
    Type(MQC_Vector)::MetricEigenvalues,OrthogonalEigenvalues
    Type(MQC_Matrix)::MetricEigenvectors,OrthogonalMatrix,OrthogonalEigenvectors
    Complex(Kind=real64),Allocatable::Hamiltonian(:,:),MetricVectors(:,:),&
      Transformation(:,:),OrthogonalHamiltonian(:,:),OrthogonalVectors(:,:),&
      GeneralizedVectors(:,:)
    Real(Kind=real64)::MetricValue,MetricCutoff,MaximumMetricValue
    Integer(Kind=int64)::i,j,Dimension,RankIndex

    Dimension=Size(ProjectedMatrix,1)
    Call ProjectedMetric%Diag(MetricEigenvalues,MetricEigenvectors)
    MaximumMetricValue=0.0_real64
    Do i=1,Dimension
      MetricValue=MetricEigenvalues%At(i)
      MaximumMetricValue=Max(MaximumMetricValue,MetricValue)
    EndDo
    MetricCutoff=Max(1.0e-14_real64,MetricTolerance*MaximumMetricValue)
    If(MaximumMetricValue.le.0.0_real64) &
      Call MQC_Error('Projected Davidson metric has no positive support')
    MetricRank=0
    Do i=1,Dimension
      If(MetricEigenvalues%At(i).lt.-MetricCutoff) &
        Call MQC_Error('Projected Davidson metric is indefinite')
      If(MetricEigenvalues%At(i).gt.MetricCutoff) MetricRank=MetricRank+1
    EndDo

    Allocate(Hamiltonian(Dimension,Dimension),MetricVectors(Dimension,MetricRank),&
      Transformation(Dimension,MetricRank),OrthogonalHamiltonian(MetricRank,MetricRank),&
      OrthogonalVectors(MetricRank,MetricRank),GeneralizedVectors(Dimension,MetricRank))
    RankIndex=0
    Do j=1,Dimension
      MetricValue=MetricEigenvalues%At(j)
      Do i=1,Dimension
        Hamiltonian(i,j)=ProjectedMatrix%At(i,j)
      EndDo
      If(MetricValue.le.MetricCutoff) Cycle
      RankIndex=RankIndex+1
      Do i=1,Dimension
        MetricVectors(i,RankIndex)=MetricEigenvectors%At(i,j)
      EndDo
      Transformation(:,RankIndex)=MetricVectors(:,RankIndex)/Sqrt(MetricValue)
    EndDo
    OrthogonalHamiltonian=MatMul(Conjg(Transpose(Transformation)),&
      MatMul(Hamiltonian,Transformation))
    Call OrthogonalMatrix%Init(MetricRank,MetricRank,&
      Scalar=(0.0_real64,0.0_real64))
    Do j=1,MetricRank
      Do i=1,MetricRank
        Call OrthogonalMatrix%Put(OrthogonalHamiltonian(i,j),i,j)
      EndDo
    EndDo
    Call OrthogonalMatrix%Diag(OrthogonalEigenvalues,OrthogonalEigenvectors)
    Do j=1,MetricRank
      Do i=1,MetricRank
        OrthogonalVectors(i,j)=OrthogonalEigenvectors%At(i,j)
      EndDo
    EndDo
    GeneralizedVectors=MatMul(Transformation,OrthogonalVectors)
    Eigenvalues=OrthogonalEigenvalues
    Call Eigenvectors%Init(Dimension,MetricRank,Scalar=(0.0_real64,0.0_real64))
    Do j=1,MetricRank
      Do i=1,Dimension
        Call Eigenvectors%Put(GeneralizedVectors(i,j),i,j)
      EndDo
    EndDo
  End Subroutine MQC_Iterative_Generalized_Ritz

!> Construct the Hermitian Ritz matrix in the current Davidson subspace.
!> \author L. M. Thompson
!> \date 2026
  Subroutine MQC_Iterative_Projected_Matrix(Basis,OperatorBasis,ProjectedMatrix)
    Complex(Kind=real64),Intent(In)::Basis(:,:),OperatorBasis(:,:)
    Type(MQC_Matrix),Intent(Out)::ProjectedMatrix
    Complex(Kind=real64),Allocatable::Projected(:,:)
    Integer(Kind=int64)::Dimension

    Dimension=Size(Basis,2)
    Allocate(Projected(Dimension,Dimension))
    Call MQC_Iterative_Gemm('C',Basis,OperatorBasis,Projected,&
      (1.0_real64,0.0_real64),(0.0_real64,0.0_real64))
    Projected=0.5_real64*(Projected+Conjg(Transpose(Projected)))
    ProjectedMatrix=Projected
  End Subroutine MQC_Iterative_Projected_Matrix

!> Orthonormalize a block and discard numerically dependent columns.
!> \author L. M. Thompson
!> \date 2026
  Subroutine MQC_Iterative_Orthonormalize_Block(Vectors,NumberOfVectors)
    Complex(Kind=real64),Intent(InOut)::Vectors(:,:)
    Integer(Kind=int64),Intent(InOut)::NumberOfVectors
    Complex(Kind=real64),Allocatable::WorkVector(:)
    Real(Kind=real64)::NormValue
    Integer(Kind=int64)::Source,Target

    Allocate(WorkVector(Size(Vectors,1)))
    Target=0
    Do Source=1,NumberOfVectors
      WorkVector=Vectors(:,Source)
      If(Target.gt.0) Call MQC_Iterative_Orthogonalize(WorkVector,Vectors(:,1:Target))
      NormValue=Sqrt(Real(Dot_Product(WorkVector,WorkVector),Kind=real64))
      If(NormValue.gt.1.0e-12_real64) Then
        Target=Target+1
        Vectors(:,Target)=WorkVector/NormValue
      EndIf
    EndDo
    NumberOfVectors=Target
  End Subroutine MQC_Iterative_Orthonormalize_Block

!> Reorthogonalize one vector against a basis using two MGS passes.
!> \author L. M. Thompson
!> \date 2026
  Subroutine MQC_Iterative_Orthogonalize(Vector,Basis)
    Complex(Kind=real64),Intent(InOut)::Vector(:)
    Complex(Kind=real64),Intent(In)::Basis(:,:)
    Complex(Kind=real64),Allocatable::Overlap(:)
    Integer(Kind=int32)::N,M
    Integer::Pass
    External::zgemv

    N=Int(Size(Basis,1),Kind=int32)
    M=Int(Size(Basis,2),Kind=int32)
    If(M.eq.0) Return
    Allocate(Overlap(M))
    Do Pass=1,2
      Call zgemv('C',N,M,(1.0_real64,0.0_real64),Basis,N,Vector,1_int32,&
        (0.0_real64,0.0_real64),Overlap,1_int32)
      Call zgemv('N',N,M,(-1.0_real64,0.0_real64),Basis,N,Overlap,1_int32,&
        (1.0_real64,0.0_real64),Vector,1_int32)
    EndDo
  End Subroutine MQC_Iterative_Orthogonalize

!> Orthogonalize a block of vectors against an orthonormal basis (two
!> classical Gram-Schmidt passes with ZGEMM), then orthonormalize the block
!> and drop dependent columns.  Kept columns are packed to the front.
  Subroutine MQC_Iterative_Block_Orthogonalize(Vectors,NumberOfVectors,Basis)
    Complex(Kind=real64),Intent(InOut)::Vectors(:,:)
    Integer(Kind=int64),Intent(InOut)::NumberOfVectors
    Complex(Kind=real64),Intent(In)::Basis(:,:)
    Complex(Kind=real64),Allocatable::Overlap(:,:)
    Integer::Pass

    If(Size(Basis,2).gt.0.and.NumberOfVectors.gt.0) Then
      Allocate(Overlap(Size(Basis,2),NumberOfVectors))
      Do Pass=1,2
        Call MQC_Iterative_Gemm('C',Basis,Vectors(:,1:NumberOfVectors),Overlap,&
          (1.0_real64,0.0_real64),(0.0_real64,0.0_real64))
        Call MQC_Iterative_Gemm('N',Basis,Overlap,Vectors(:,1:NumberOfVectors),&
          (-1.0_real64,0.0_real64),(1.0_real64,0.0_real64))
      EndDo
    EndIf
    Call MQC_Iterative_Orthonormalize_Block(Vectors,NumberOfVectors)
  End Subroutine MQC_Iterative_Block_Orthogonalize

!> C = Alpha*op(A)*B + Beta*C with BLAS ZGEMM, op(A) = A ('N') or A^H ('C').
  Subroutine MQC_Iterative_Gemm(TransA,A,B,C,Alpha,Beta)
    Character(Len=1),Intent(In)::TransA
    Complex(Kind=real64),Intent(In)::A(:,:),B(:,:),Alpha,Beta
    Complex(Kind=real64),Intent(InOut)::C(:,:)
    Integer(Kind=int32)::M,N,K,LDA,LDB,LDC
    External::zgemm

    M=Int(Size(C,1),Kind=int32)
    N=Int(Size(C,2),Kind=int32)
    If(TransA.eq.'N') Then
      K=Int(Size(A,2),Kind=int32)
    Else
      K=Int(Size(A,1),Kind=int32)
    EndIf
    If(M.eq.0.or.N.eq.0) Return
    If(K.eq.0) Then
      C=Beta*C
      Return
    EndIf
    LDA=Int(Max(1_int64,Int(Size(A,1),Kind=int64)),Kind=int32)
    LDB=Int(Max(1_int64,Int(Size(B,1),Kind=int64)),Kind=int32)
    LDC=Int(Max(1_int64,Int(Size(C,1),Kind=int64)),Kind=int32)
    Call zgemm(TransA,'N',M,N,K,Alpha,A,LDA,B,LDB,Beta,C,LDC)
  End Subroutine MQC_Iterative_Gemm

End Submodule MQC_Algebra_Iterative
