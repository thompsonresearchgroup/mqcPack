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
      RitzVectors=MatMul(Basis(:,1:SubspaceSize),&
        ProjectedVectors(1:SubspaceSize,1:RequestedRoots))
      OperatorRitz=MatMul(OperatorBasis(:,1:SubspaceSize),&
        ProjectedVectors(1:SubspaceSize,1:RequestedRoots))
      MetricRitz=MatMul(MetricBasis(:,1:SubspaceSize),&
        ProjectedVectors(1:SubspaceSize,1:RequestedRoots))
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
        Call MQC_Iterative_Orthogonalize(WorkVector,Basis(:,1:SubspaceSize))
        If(CorrectionCount.gt.0) Call MQC_Iterative_Orthogonalize(WorkVector,&
          Corrections(:,1:CorrectionCount))
        NormValue=Sqrt(Real(Dot_Product(WorkVector,WorkVector),Kind=real64))
        If(NormValue.gt.1.0e-12_real64) Then
          CorrectionCount=CorrectionCount+1
          Corrections(:,CorrectionCount)=WorkVector/NormValue
        EndIf
      EndDo
      AllConverged=All(RootResiduals.le.ConvergenceTolerance)
      If(AllConverged.or.CorrectionCount.eq.0) Exit

      If(SubspaceSize+CorrectionCount.gt.MaxSubspace) Then
        Basis=(0.0_real64,0.0_real64)
        Basis(:,1:RequestedRoots)=RitzVectors
        SubspaceSize=RequestedRoots
        Call MQC_Iterative_Orthonormalize_Block(Basis,SubspaceSize)
        NewCount=0
        Do k=1,Min(CorrectionCount,MaxSubspace-SubspaceSize)
          WorkVector=Corrections(:,k)
          Call MQC_Iterative_Orthogonalize(WorkVector,Basis(:,1:SubspaceSize))
          If(NewCount.gt.0) Call MQC_Iterative_Orthogonalize(WorkVector,&
            Basis(:,SubspaceSize+1:SubspaceSize+NewCount))
          NormValue=Sqrt(Real(Dot_Product(WorkVector,WorkVector),Kind=real64))
          If(NormValue.gt.1.0e-12_real64) Then
            NewCount=NewCount+1
            Basis(:,SubspaceSize+NewCount)=WorkVector/NormValue
          EndIf
        EndDo
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
    Complex(Kind=real64)::Element
    Integer(Kind=int64)::i,j,Dimension

    Dimension=Size(Basis,2)
    Call ProjectedMatrix%Init(Dimension,Dimension,&
      Scalar=(0.0_real64,0.0_real64))
    Do j=1,Dimension
      Do i=j,Dimension
        Element=Dot_Product(Basis(:,i),OperatorBasis(:,j))
        If(i.eq.j) Element=Cmplx(Real(Element,Kind=real64),0.0_real64,Kind=real64)
        Call ProjectedMatrix%Put(Element,i,j)
        Call ProjectedMatrix%Put(Conjg(Element),j,i)
      EndDo
    EndDo
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
    Complex(Kind=real64)::Overlap
    Integer(Kind=int64)::Pass,i

    Do Pass=1,2
      Do i=1,Size(Basis,2)
        Overlap=Dot_Product(Basis(:,i),Vector)
        Vector=Vector-Basis(:,i)*Overlap
      EndDo
    EndDo
  End Subroutine MQC_Iterative_Orthogonalize

End Submodule MQC_Algebra_Iterative
