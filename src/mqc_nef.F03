

!>    \brief <b> MQC NEF contains objects for the manipulation of intermediates used in
!>    the Natural Excitation Framework </b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC NEF contains objects that speed up and enhance the implementation of 
!>    orbital and excitation subspaces in Fortan. The objects allow for automatic 
!>    determination of orbitais subspace maintaining the spin structure .
!>    The NEF derived types defined in this package are:
!>
!>    1. MQC_NEF_Vector:         NEF
!>    1. MQC_NEF_Matrix:         NEF
!>    1. MQC_NEF_R4:             NEF
!>    1. MQC_NEF_Set:            NEF
!>
!>    This module is level 3 in the MQC hierarchy and so depends on level
!>     0, 1 and 2 modules. 
!>
!>    \endverbatim
!
      Module MQC_NEF
!
!     ***********************************************************************
!     ***********************************************************************
!     **                                                                   **
!     **               The Merced Quantum Chemistry Package                **
!     **                            (MQCPack)                              **
!     **                                                                   **
!     ** Written By:                                                       **
!     **    Lee M. Thompson, Xianghai Sheng, Andrew D. Mahler, Dave        **
!     **    Mullally, and Hrant P. Hratchian                               **
!     **                                                                   **
!     **                          Version 26.8.0                           **
!     **                          August 3, 2026                           **
!     **                                                                   **
!     **                                                                   **
!     **                                                                   **
!     ** Modules beloning to MQCPack:                                      **
!     **    1. MQC_General                                                 **
!     **    2. MQC_Binary                                                  **
!     **    3. MQC_DataStructures                                          **
!     **    4. MQC_Integrals                                               **
!     **    5. MQC_Algebra                                                 **
!     **    6. MQC_Files                                                   **
!     **    7. MQC_Molecule                                                **
!     **    8. MQC_EST                                                     **
!     **    9. MQC_Gaussian                                                **
!     **   10. MQC_NEF                                                     **
!     **                                                                   **
!     ***********************************************************************
!     ***********************************************************************
!
      use mqc_est
      Use iso_fortran_env, only: int32, int64, real64
!
!----------------------------------------------------------------
!                                                               |
!     TYPE AND CLASS DEFINITIONS                                |
!                                                               |
!----------------------------------------------------------------
!
!
!>    \brief <b> </b>
      Type MQC_NEF_Vector
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Main_Vector
        Logical,Dimension(:),Allocatable,Private::Zero_Flag
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Map_Matrix
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha,Size_vector_beta
        Integer,Private::NumSpace=0
      Contains
        Procedure, Public::print => mqc_print_nef_vector
        Procedure, Public::init => mqc_nef_vector_initialize
        Procedure, public::getBlockVec => mqc_nef_vector_output_block_vector
        Procedure, Public::impv => mqc_nef_vector_import_vector
        Procedure, Public::expand => mqc_nef_vector_dimention_expand
        Procedure, Public::ins => mqc_nef_vector_insert
        Procedure, Public::dimen => mqc_nef_vector_dimension
!>    \brief <b> Outer Product of two variables into a MQC NEF Matrix</b>
        Procedure, Public::outer2Matrix => mqc_nefvector_nefvector_outer_product
      End Type
!
!>    \brief <b> </b>
      Type MQC_NEF_Matrix
        Type(MQC_SCF_Integral),Dimension(:),Allocatable,Private::Main_Matrix
        Logical,Dimension(:),Allocatable,Private::Zero_Flag
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Map_Matrix_Col,Map_Matrix_Row
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_col,Size_vector_beta_col
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_row,Size_vector_beta_row
        Integer,Private::NumSpace_col=0,NumSpace_row=0
      Contains
        Procedure, Public::print => mqc_print_nef_matrix
        Procedure, Public::init => mqc_nef_matrix_initialize
        Procedure, public::collapseIndex => mqc_nef_matrix_collapse_index
        Procedure, public::collapseMap => mqc_nef_matrix_collapse_map
        Procedure, public::blockIndex => mqc_nef_matrix_index
        Procedure, public::at => mqc_nef_matrix_at
        Procedure, Public::imp => mqc_nef_matrix_import_integral
        Procedure, Public::export => mqc_nef_matrix_export_integral
        Procedure, Public::slice => mqc_nef_matrix_dimention_slice
        Procedure, Public::collapse => mqc_nef_matrix_dimention_collapse
        Procedure, Public::expand => mqc_nef_matrix_dimention_expand
        Procedure, Public::dimen => mqc_nef_matrix_dimension
!>    \brief <b> Outer Product of two variables into a two NEF R4</b>
        Procedure, Public::outer2R4 => mqc_NEFMatrix_NEFMatrix_outer_product
!>    \brief <b> Outer Product of two variables into a NEF R4</b>
        Procedure, Public::crossouter2R4 => mqc_NEFMatrix_NEFMatrix_crossouter_product
        Procedure, public::getBlockMat => mqc_nef_matrix_output_block_matrix
        Procedure, public::getMap => mqc_nef_matrix_get_map
        Procedure, Public::pinvMat => mqc_nef_matrix_block_diagonal_pinv
      End Type
!
!>    \brief <b> </b>
      Type MQC_NEF_R4
        Type(MQC_TwoERIs),Dimension(:),Allocatable,Private::Main_R4
        Logical,Dimension(:),Allocatable,Private::Zero_Flag
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Map_Matrix_1,Map_Matrix_2,&
                                                                    Map_Matrix_3,Map_Matrix_4
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_1,Size_vector_beta_1
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_2,Size_vector_beta_2
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_3,Size_vector_beta_3
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_4,Size_vector_beta_4
        Integer,Private::NumSpace_1=0,NumSpace_2=0,NumSpace_3=0,NumSpace_4=0
      Contains
        Procedure, Public::print => mqc_print_nef_r4
        Procedure, public::at => mqc_nef_r4_at
        Procedure, public::blockSize => mqc_nef_r4_dimension
        Procedure, public::blockSizeVec => mqc_nef_r4_dimension_vector
        Procedure, public::blockMap => mqc_nef_r4_map
        Procedure, public::blockIndex => mqc_nef_r4_index
        Procedure, Public::collapse => mqc_nef_r4_dimention_collapse
        Procedure, Public::collapseIndex => mqc_nef_r4_collapse_index
        Procedure, public::collapseMap => mqc_nef_r4_collapse_map
        procedure, public::swap => mqc_nef_r4_dimension_swap
        Procedure, Public::pinvR4 => mqc_nef_r4_block_diagonal_pinv
        Procedure, Public::put => mqc_nef_r4_put
      End Type
!
!>    \brief <b> Holds a dynamic (I,J)-indexed grid of exactly one NEF kind, determined lazily on first write </b>
      Type MQC_NEF_Set
        Character(Len=64),Private::Array_Type='open'   ! 'open','vector','matrix','r4'
        Integer,Private::NumI=0,NumJ=0
        Type(MQC_NEF_Vector),Dimension(:,:),Allocatable,Private::Set_Vector
        Type(MQC_NEF_Matrix),Dimension(:,:),Allocatable,Private::Set_Matrix
        Type(MQC_NEF_R4),Dimension(:,:),Allocatable,Private::Set_R4
        Logical,Dimension(:,:),Allocatable,Private::Zero_Flag
        Integer,Dimension(4)::NumSpace = -1
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_1,Size_vector_beta_1
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_2,Size_vector_beta_2
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_3,Size_vector_beta_3
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_4,Size_vector_beta_4
      Contains
        Procedure, Public::init => mqc_nef_set_initialize
        Procedure, Public::isOpen => mqc_nef_set_is_open
        Procedure, Public::isVector => mqc_nef_set_is_vector
        Procedure, Public::isMatrix => mqc_nef_set_is_matrix
        Procedure, Public::isR4 => mqc_nef_set_is_r4
        Procedure, Public::isZero => mqc_nef_set_is_zero
        Procedure, Public::Type => mqc_nef_set_array_type
        Procedure, Public::putVector => mqc_nef_set_put_vector
        Procedure, Public::putMatrix => mqc_nef_set_put_matrix
        Procedure, Public::putR4 => mqc_nef_set_put_r4
        Generic,   Public::put => putVector, putMatrix, putR4
        Procedure, Public::atVector => mqc_nef_set_at_vector
        Procedure, Public::atMatrix => mqc_nef_set_at_matrix
        Procedure, Public::atR4 => mqc_nef_set_at_r4
        Procedure, Public::collapse => mqc_nef_set_collapse
        Procedure, Public::expand => mqc_nef_set_expand
        Procedure, Public::outerCollapse => mqc_nef_set_outer_collapse
        Procedure, Public::outerExpand => mqc_nef_set_outer_expand
        Procedure, Public::fullCollapse => mqc_nef_set_full_collapse
        Procedure, Public::fullExpand => mqc_nef_set_full_expand
        Procedure, Public::getMap => mqc_nef_set_get_map
      End Type
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
!
!>    \brief <b> Multiplies two objects</b>
      Interface MatMul
        Module Procedure MQC_NEFMatrix_Integral_Multiply
        Module Procedure MQC_Integral_NEFMatrix_Multiply
        Module Procedure MQC_NEFMatrix_NEFMatrix_Multiply
        Module Procedure MQC_NEFMatrix_NEFVector_Multiply
        Module Procedure MQC_NEFSet_Vector_Multiply
        Module Procedure MQC_Vector_NEFSet_Multiply
        Module Procedure MQC_NEFSet_Matrix_Multiply
        Module Procedure MQC_Matrix_NEFSet_Multiply
        Module Procedure MQC_NEFSet_NEFMatrix_Multiply
        Module Procedure MQC_NEFMatrix_NEFSet_Multiply
        Module Procedure MQC_NEFSet_NEFSet_Multiply
      End Interface
!
!>    \brief <b> Returns the Hermitian conjugate</b>
      Interface Dagger
        Module Procedure MQC_NEFMatrix_Conjugate_Transpose
        Module Procedure MQC_NEFR4_Conjugate_Transpose
        Module Procedure MQC_NEFSet_Conjugate_Transpose
      End Interface
!
!>    \brief <b> Returns the conjugate</b>
      Interface Conjg
        Module Procedure MQC_NEFMatrix_Conjugate
      End Interface
!
!>    \brief <b> Returns the transpose</b>
      Interface Transpose
        Module Procedure MQC_NEFMatrix_Transpose
        Module Procedure MQC_NEFSet_Transpose
      End Interface
!
!>    \brief <b> Contracts two objects</b>
      Interface Contraction
        Module Procedure mqc_NEFVector_NEFVector_contraction
        Module Procedure mqc_NEFMatrix_NEFMatrix_contraction
        Module Procedure mqc_NEF_R4_NEF_R4_contraction
        Module Procedure MQC_NEF_R4_NEF_Matrix_Contraction
        Module Procedure MQC_NEFSet_NEFSet_Contraction
      End Interface
!
!>    \brief <b> Contracts two objects</b>
      Interface PartialContraction
        Module Procedure mqc_NEF_R4_NEF_R4_partial_contraction
        Module Procedure MQC_NEF_R4_NEF_Matrix_Contraction
        Module Procedure MQC_NEFSet_NEFSet_Partial_Contraction
      End Interface
!
!>    \brief <b> Contracts two objects</b>
      Interface PartialContraction2R4
        Module Procedure mqc_NEF_R4_NEF_R4_partial_contraction_2NEF_R4
        Module Procedure mqc_NEF_R4_ERI_partial_contraction_2NEF_R4
        Module Procedure mqc_ERI_NEF_R4_partial_contraction_2NEF_R4
        Module Procedure MQC_NEFSet_NEFSet_Partial_Contraction_2NEF_R4
      End Interface
!
!>    \brief <b> Fully contracts two objects to a scalar</b>
      Interface FullContraction
        Module Procedure mqc_NEFSet_NEFSet_full_contraction
      End Interface
!
!>    \brief <b> Returns the block length along the specified dimension</b>
      Interface Size
        Module Procedure mqc_nef_r4_dimension
      End Interface
!

!
!
!----------------------------------------------------------------
!                                                               |
!     OPERATOR INTERFACES                                       |
!                                                               |
!----------------------------------------------------------------
!
!
!>    \brief <b> Sums two variables</b>
      Interface Operator (+)
        Module Procedure MQC_NEF_Vector_Sum
        Module Procedure MQC_NEFMatrix_Sum
        Module Procedure MQC_NEF_R4_Sum
        Module Procedure MQC_NEFSet_NEFVector_Sum
        Module Procedure MQC_NEFVector_NEFSet_Sum
        Module Procedure MQC_NEFSet_NEFMatrix_Sum
        Module Procedure MQC_NEFMatrix_NEFSet_Sum
        Module Procedure MQC_NEFSet_NEFR4_Sum
        Module Procedure MQC_NEFR4_NEFSet_Sum
        Module Procedure MQC_NEFSet_NEFSet_Sum
      End Interface
!      
!>    \brief <b> Multiplies two variables</b>
      Interface Operator (*)
        Module Procedure MQC_Scalar_NEF_Vector_Multiply
        Module Procedure MQC_NEF_Vector_Scalar_Multiply
        Module Procedure MQC_Scalar_NEFMatrix_Multiply
        Module Procedure MQC_NEFMatrix_Scalar_Multiply
        Module Procedure MQC_Scalar_NEFR4_Multiply
        Module Procedure MQC_NEFR4_Scalar_Multiply
        Module Procedure MQC_Scalar_NEFSet_Multiply
        Module Procedure MQC_NEFSet_Scalar_Multiply
      End Interface

      !>    \brief <b> Subtracts two variables</b>
      Interface Operator (-)
        Module Procedure MQC_NEF_Vector_Difference
        Module Procedure MQC_NEFMatrix_Difference
        Module Procedure MQC_NEF_R4_Difference
        Module Procedure MQC_NEFSet_NEFVector_Difference
        Module Procedure MQC_NEFVector_NEFSet_Difference
        Module Procedure MQC_NEFSet_NEFMatrix_Difference
        Module Procedure MQC_NEFMatrix_NEFSet_Difference
        Module Procedure MQC_NEFSet_NEFR4_Difference
        Module Procedure MQC_NEFR4_NEFSet_Difference
        Module Procedure MQC_NEFSet_NEFSet_Difference
      End Interface

!    
!
!
!----------------------------------------------------------------
!                                                               |
!     SUBROUTINES AND FUNCTIONS                                 |
!                                                               |
!----------------------------------------------------------------
!
!
      CONTAINS
!
!
!     PROCEDURE MQC_Print_NEF_Vector
!
!>    \brief <b> MQC_Print_NEF_Vector is a subroutine used to print an MQC
!>    NEF Vector object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_NEF_Vector is a subroutine used to print an MQC NEF Vector object. 
!>    Blank_At_Top and Blank_At_Bottom are optional logical arguments to print blank
!>    lines before or after output.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVec
!>    \verbatim
!>        Integral is Class(MQC_NEF_Vector)
!>        The MQC NEF Vector variable to be printed.
!>    \endverbatim
!>
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64) 
!>        The Fortran file number to print to. 
!>    \endverbatim
!>
!>    \param[in] Header
!>    \verbatim
!>        Header is Character(Len=*)
!>        The title to print along with Integral.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Top
!>    \verbatim
!>        Blank_At_Top is Logical,Optional
!>        = .True.:  print blank line above output
!>        = .False.: do not print blank line above output.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Bottom
!>    \verbatim
!>        Blank_At_Bottom is Logical,Optional
!>        = .True.:  print blank line below output
!>        = .False.: do not print blank line below output.
!>    \endverbatim
!>
!>    \param[in] FormatStr
!>    \verbatim
!>        FormatStr is Character(Len=*),Optional
!>        Formatting statement for printing.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M.M.F. Moraes
!>    \date 2026
!
      subroutine mqc_print_nef_vector(NEFVec,iOut,header, &
          blank_at_top,blank_at_bottom,formatStr)
!
        implicit none
        class(mqc_nef_vector)::NEFVec
        integer(kind=int64),intent(in)::iOut
        character(len=*),intent(in)::header
        character(len=64)::arrayType
        logical,intent(In),optional::blank_at_top,blank_at_bottom
        character(len=*),intent(in),optional::formatStr
        integer::i
!
 1000 Format(1x,A)
 1020 Format( " " )
!  
        if(present(blank_at_top)) then
          if(blank_at_top) write(iout,1020)
        endif
        write(iout,1000) trim(header)
        do i = 1, NEFVec%NumSpace
          if (NEFVec%Zero_Flag(i)) then
            write(iout,1000) 'Sub-Space '//trim(num2char(i))//' (all zero)'
          else
            call mqc_print(NEFVec%Main_vector(i),iOut,'Sub-Space '//trim(num2char(i)))
          end if
        end do
        
      end subroutine mqc_print_nef_vector
!
!
!
      subroutine mqc_nef_vector_initialize(NEFVec,initial_map_matrix,scalar)
!
        implicit none
        class(mqc_nef_vector),intent(inOut)::NEFVec
        type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in)::initial_map_matrix
        Class(*),Optional::scalar
        integer::i
        type(mqc_scalar)::scalarCheck
        logical::scalarIsZero
!
        if(Present(Scalar)) then
          select type(scalar)
          type is (integer)
            scalarCheck = scalar
          type is (real)
            scalarCheck = scalar
          type is (complex)
            scalarCheck = scalar
          type is (mqc_scalar)
            scalarCheck = scalar
          class default
            call mqc_error('mqc_nef_vector_initialize: scalar type not defined')
          end select
          scalarIsZero = (abs(scalarCheck).eq.0.0)
        else
          scalarIsZero = .true.
        end if

        if(allocated(NEFVec%Main_vector)) deallocate(NEFVec%Main_vector)
        allocate(NEFVec%Main_vector(size(initial_map_matrix)))
        if(allocated(NEFVec%Zero_Flag)) deallocate(NEFVec%Zero_Flag)
        allocate(NEFVec%Zero_Flag(size(initial_map_matrix)))
        if(allocated(NEFVec%Size_vector_alpha)) deallocate(NEFVec%Size_vector_alpha)
        allocate(NEFVec%Size_vector_alpha(size(initial_map_matrix)))
        if(allocated(NEFVec%Size_vector_beta)) deallocate(NEFVec%Size_vector_beta)
        allocate(NEFVec%Size_vector_beta(size(initial_map_matrix)))
        NEFVec%NumSpace = size(initial_map_matrix)
        NEFVec%Map_Matrix = initial_map_matrix
        do i = 1, NEFVec%NumSpace
          NEFVec%Zero_Flag(i) = scalarIsZero
          select case(initial_map_matrix(i)%type())
          case ('space')
            if (.not.scalarIsZero) &
              call NEFVec%Main_vector(i)%init(&
                     Size(initial_map_matrix(i),'alpha'),&
                     scalar=scalar)
            NEFVec%Size_vector_alpha(i) = Size(initial_map_matrix(i),'alpha')
            NEFVec%Size_vector_beta(i) = Size(initial_map_matrix(i),'alpha')
          case ('spin')
            if (.not.scalarIsZero) &
              call NEFVec%Main_vector(i)%init(&
                     Size(initial_map_matrix(i),'alpha'),&
                     scalar=scalar,&
                     lengthBeta=Size(initial_map_matrix(i),'beta'))
            NEFVec%Size_vector_alpha(i) = Size(initial_map_matrix(i),'alpha')
            NEFVec%Size_vector_beta(i) = Size(initial_map_matrix(i),'beta')
          case default
            call mqc_error_A('type not valid in mqc_nef_vector_initialize', 6, &
              'initial_map_matrix(i)%type()', initial_map_matrix(i)%type() )
          end select
        end do
!
      end subroutine mqc_nef_vector_initialize
!
!     PROCEDURE MQC_Print_NEF_Matrix
!
!>    \brief <b> MQC_Print_NEF_Matrix is a subroutine used to print an MQC
!>    NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_NEF_Matrix is a subroutine used to print an MQC NEF Matrix object. 
!>    Blank_At_Top and Blank_At_Bottom are optional logical arguments to print blank
!>    lines before or after output.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        Integral is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix variable to be printed.
!>    \endverbatim
!>
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64) 
!>        The Fortran file number to print to. 
!>    \endverbatim
!>
!>    \param[in] Header
!>    \verbatim
!>        Header is Character(Len=*)
!>        The title to print along with Integral.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Top
!>    \verbatim
!>        Blank_At_Top is Logical,Optional
!>        = .True.:  print blank line above output
!>        = .False.: do not print blank line above output.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Bottom
!>    \verbatim
!>        Blank_At_Bottom is Logical,Optional
!>        = .True.:  print blank line below output
!>        = .False.: do not print blank line below output.
!>    \endverbatim
!>
!>    \param[in] FormatStr
!>    \verbatim
!>        FormatStr is Character(Len=*),Optional
!>        Formatting statement for printing.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M.M.F. Moraes
!>    \date 2026
!
      subroutine mqc_print_nef_matrix(NEFMat,iOut,header, &
          blank_at_top,blank_at_bottom,formatStr,blockPrint)
!
      implicit none
      integer(kind=int64) :: nef_bidx
      class(mqc_nef_matrix)::NEFMat
      integer(kind=int64),intent(in)::iOut
      character(len=*),intent(in)::header
      logical,optional::blockPrint
      character(len=64)::arrayType
      type(mqc_scf_integral)::tmp_int
      logical,intent(In),optional::blank_at_top,blank_at_bottom
      character(len=*),intent(in),optional::formatStr
      integer::i,j,sub_i,sub_j,line_num,col_num
      logical::my_block_print
      Character(Len=256)::Line
!
      if (present(blockPrint)) then
        my_block_print = blockPrint
      else
        my_block_print = .true.
      end if
!
 1000 Format(1x,A)
 1020 Format( " " )
!  
      if(present(blank_at_top)) then
        if(blank_at_top) write(iout,1020)
      endif
      write(iout,1000) trim(header)
      line_num = 0
      if (my_block_print) then
        write(6,*) 'Alpha-Alpha'
        do j = 1, NEFMat%NumSpace_row
          if (j.eq.1) then 
            Line = ' '
            do i = 1, NEFMat%NumSpace_col
              do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                Line = trim(Line)//'           '//num2char(NEFMat%Map_Matrix_Col(i)%at(sub_i,'alpha'))
              end do
            end do
            write(*,*) Line
          end if
          do sub_j = 1, NEFMat%Size_vector_alpha_row(j)
            line_num = line_num + 1
            col_num = 1 
            Line = num2char(NEFMat%Map_Matrix_Row(j)%at(sub_j,'alpha'))
            do i = 1, NEFMat%NumSpace_col
              nef_bidx = NEFMat%blockIndex([1,j],[2,i])
              if (NEFMat%Zero_Flag(nef_bidx)) then
                do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                  Line = trim(Line)//'     '//num2char(0.00000)
                  if (sub_i.ne.NEFMat%Size_vector_alpha_col(i)) Line = trim(Line)//' '
                end do
              else
              tmp_int = NEFMat%Main_matrix(nef_bidx)
              do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                if (tmp_int%at(sub_j,sub_i,'alpha').eq.0.0) then
                  Line = trim(Line)//'     '//num2char(0.00000)
                else if (tmp_int%at(sub_j,sub_i,'alpha').lt.0.0) then
                  Line = trim(Line)//'    '//num2char(tmp_int%at(sub_j,sub_i,'alpha'))
                else
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'alpha'))
                end if
                if (sub_i.ne.NEFMat%Size_vector_alpha_col(i)) Line = trim(Line)//' '
              end do
              end if
              if(i.ne.NEFMat%NumSpace_col) Line = trim(Line)//'|'
            end do
            write(6,*) Line
          end do
          if(j.ne.NEFMat%NumSpace_row) &
             write(6,*) '----------------------------------------------------------------------&
                    &----------------------------------------------------------------------'    
        end do
        line_num = 0
        write(6,*) 'Beta-Beta'
        do j = 1, NEFMat%NumSpace_row
          do sub_j = 1, NEFMat%Size_vector_beta_row(j)
            line_num = line_num + 1
            col_num = 1
            Line = num2char(NEFMat%Map_Matrix_Row(j)%at(sub_j,'beta'))
            do i = 1, NEFMat%NumSpace_col
              nef_bidx = NEFMat%blockIndex([1,j],[2,i])
              if (NEFMat%Zero_Flag(nef_bidx)) then
                do sub_i = 1, NEFMat%Size_vector_beta_col(i)
                  Line = trim(Line)//'     '//num2char(0.00000)
                  if (sub_i.ne.NEFMat%Size_vector_beta_col(i)) Line = trim(Line)//' '
                end do
              else
              tmp_int = NEFMat%Main_matrix(nef_bidx)
              do sub_i = 1, NEFMat%Size_vector_beta_col(i)
                if (tmp_int%at(sub_j,sub_i,'beta').eq.0.0) then
                  Line = trim(Line)//'     '//num2char(0.00000)
                else if (tmp_int%at(sub_j,sub_i,'beta').lt.0.0) then
                  Line = trim(Line)//'    '//num2char(tmp_int%at(sub_j,sub_i,'beta'))
                else
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'beta'))
                end if
                if (sub_i.ne.NEFMat%Size_vector_beta_col(i)) Line = trim(Line)//' '
              end do
              end if
              if(i.ne.NEFMat%NumSpace_col) Line = trim(Line)//'|'
            end do
            write(6,*) trim(Line)
          end do
          if(j.ne.NEFMat%NumSpace_row) &
             write(6,*) '----------------------------------------------------------------------&
                    &----------------------------------------------------------------------'
        end do
        line_num = 0
        write(6,*) 'Alpha-Beta'
        do j = 1, NEFMat%NumSpace_row
          do sub_j = 1, NEFMat%Size_vector_beta_row(j)
            line_num = line_num + 1
            col_num = 1
            Line = num2char(NEFMat%Map_Matrix_Row(j)%at(sub_j,'beta'))
            do i = 1, NEFMat%NumSpace_col
              nef_bidx = NEFMat%blockIndex([1,j],[2,i])
              if (NEFMat%Zero_Flag(nef_bidx)) then
                do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                  Line = trim(Line)//'     '//num2char(0.00000)
                  if (sub_i.ne.NEFMat%Size_vector_beta_col(i)) Line = trim(Line)//' '
                end do
              else
              tmp_int = NEFMat%Main_matrix(nef_bidx)
              do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                if (tmp_int%at(sub_j,sub_i,'alphabeta').eq.0.0) then
                  Line = trim(Line)//'     '//num2char(0.00000)
                else if (tmp_int%at(sub_j,sub_i,'alphabeta').lt.0.0) then
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'alphabeta'))
                else
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'alphabeta'))
                end if
                if (sub_i.ne.NEFMat%Size_vector_beta_col(i)) Line = trim(Line)//' '
              end do
              end if
              if(i.ne.NEFMat%NumSpace_col) Line = trim(Line)//'|'
            end do
            write(6,*) trim(Line)
          end do
          if(j.ne.NEFMat%NumSpace_row) &
             write(6,*) '----------------------------------------------------------------------&
                    &----------------------------------------------------------------------'    
        end do
        line_num = 0
        write(6,*) 'Beta-Alpha'
        do j = 1, NEFMat%NumSpace_row
          do sub_j = 1, NEFMat%Size_vector_alpha_row(j)
            line_num = line_num + 1
            col_num = 1
            Line = num2char(NEFMat%Map_Matrix_Row(j)%at(sub_j,'alpha'))
            do i = 1, NEFMat%NumSpace_col
              nef_bidx = NEFMat%blockIndex([1,j],[2,i])
              if (NEFMat%Zero_Flag(nef_bidx)) then
                do sub_i = 1, NEFMat%Size_vector_beta_col(i)
                  Line = trim(Line)//'     '//num2char(0.00000)
                  if (sub_i.ne.NEFMat%Size_vector_alpha_col(i)) Line = trim(Line)//' '
                end do
              else
              tmp_int = NEFMat%Main_matrix(nef_bidx)
              do sub_i = 1, NEFMat%Size_vector_beta_col(i)
                if (tmp_int%at(sub_j,sub_i,'betaalpha').eq.0.0) then
                  Line = trim(Line)//'     '//num2char(0.00000)
                else if (tmp_int%at(sub_j,sub_i,'betaalpha').lt.0.0) then
                  Line = trim(Line)//'    '//num2char(tmp_int%at(sub_j,sub_i,'betaalpha'))
                else
                  Line = trim(Line)//'     '//num2char(tmp_int%at(sub_j,sub_i,'betaalpha'))
                end if
                if (sub_i.ne.NEFMat%Size_vector_alpha_col(i)) Line = trim(Line)//' '
              end do
              end if
              if(i.ne.NEFMat%NumSpace_col) Line = trim(Line)//'|'
            end do
            write(6,*) trim(Line)
          end do
          if(j.ne.NEFMat%NumSpace_row) &
             write(6,*) '----------------------------------------------------------------------&
                    &----------------------------------------------------------------------'
        end do
        write(6,*) ' '
      else
        do i = 1, NEFMat%NumSpace_col
          do j = 1, NEFMat%NumSpace_row
            nef_bidx = NEFMat%blockIndex([1,j],[2,i])
            if (NEFMat%Zero_Flag(nef_bidx)) then
              write(iOut,1000) 'Sub-Space '//trim(num2char(j))//','//trim(num2char(i))//' (all zero)'
            else
              call mqc_print(NEFMat%Main_matrix(nef_bidx),iOut,&
                    'Sub-Space '//trim(num2char(j))//','//trim(num2char(i)))
            end if
          end do
        end do
      end if
! 
      end subroutine mqc_print_nef_matrix
!
!
!     PROCEDURE MQC_NEF_Matrix_Initialize
!
!>    \brief <b> MQC_NEF_Matrix_Initialize is a subroutine that initializes a
!>    NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_Matrix_Initialize is a subroutine that initializes a MQC_NEF_Matrix. Default
!>    element values are 0.0, or otherwise vector can be initialized with optional
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        Integral is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to initialize.
!>    \endverbatim
!>
!>    \param[in] initial_map_matrix_row
!>    \verbatim
!>        initial_map_matrix_row is a type(mqc_scf_eigenvalues),dimension(:),allocatable
!>        Maps each element from an original mqc_scf_integral basis into the nef basis.
!>    \endverbatim
!>
!>    \param[in] scalar
!>    \verbatim
!>        Scalar is Class(*),Optional
!>        Value to set each element of Vector. If not present,
!>        the value is set to 0.0. Can be of type integer, real,
!>        complex or MQC_Scalar.
!>    \endverbatim
!>
!>    \param[in] initial_map_matrix_col
!>    \verbatim
!>        initial_map_matrix_row is a type(mqc_scf_eigenvalues),dimension(:),allocatable,optional
!>        Maps each element from an original mqc_scf_integral basis into the nef basis on the 
!>        second axis (columns). If not present, the mapping is the same in both axes.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_matrix_initialize(NEFMat,initial_map_matrix_row,scalar,initial_map_matrix_col)
!
        implicit none
        integer(kind=int64) :: nef_bidx
        class(mqc_nef_matrix),intent(inOut)::NEFMat
        type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in)::initial_map_matrix_row
        type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in),optional::initial_map_matrix_col
        type(mqc_scf_eigenvalues),dimension(:),allocatable::active_map_matrix_col
        Class(*),Optional::scalar
        logical::same_vec
        integer::i,j
        type(mqc_scalar)::scalarCheck
        logical::scalarIsZero

        if(Present(Scalar)) then
          select type(scalar)
          type is (integer)
            scalarCheck = scalar
          type is (real)
            scalarCheck = scalar
          type is (complex)
            scalarCheck = scalar
          type is (mqc_scalar)
            scalarCheck = scalar
          class default
            call mqc_error('mqc_nef_matrix_initialize: scalar type not defined')
          end select
          scalarIsZero = (abs(scalarCheck).eq.0.0)
        else
          scalarIsZero = .true.
        end if

        NEFMat%NumSpace_row = size(initial_map_matrix_row)
        NEFMat%Map_Matrix_row = initial_map_matrix_row
        if(present(initial_map_matrix_col)) then
          NEFMat%NumSpace_col = size(initial_map_matrix_col)
          active_map_matrix_col = initial_map_matrix_col
          NEFMat%Map_Matrix_col = initial_map_matrix_col
          same_vec = .false.
        else
          NEFMat%NumSpace_col = NEFMat%NumSpace_row
          active_map_matrix_col = initial_map_matrix_row
          NEFMat%Map_Matrix_col = initial_map_matrix_row
          same_vec = .true.
        end if
        if(allocated(NEFMat%Size_vector_alpha_row)) deallocate(NEFMat%Size_vector_alpha_row)
        if(allocated(NEFMat%Size_vector_alpha_col)) deallocate(NEFMat%Size_vector_alpha_col)
        if(allocated(NEFMat%Size_vector_beta_row)) deallocate(NEFMat%Size_vector_beta_row)
        if(allocated(NEFMat%Size_vector_beta_col)) deallocate(NEFMat%Size_vector_beta_col)
        if(allocated(NEFMat%Main_matrix)) deallocate(NEFMat%Main_matrix)
        if(allocated(NEFMat%Zero_Flag)) deallocate(NEFMat%Zero_Flag)
        allocate(NEFMat%Size_vector_alpha_row(NEFMat%NumSpace_row))
        allocate(NEFMat%Size_vector_beta_row(NEFMat%NumSpace_row))
        allocate(NEFMat%Size_vector_alpha_col(NEFMat%NumSpace_col))
        allocate(NEFMat%Size_vector_beta_col(NEFMat%NumSpace_col))

        allocate(NEFMat%Main_matrix(NEFMat%NumSpace_row*NEFMat%NumSpace_col))
        allocate(NEFMat%Zero_Flag(NEFMat%NumSpace_row*NEFMat%NumSpace_col))
        do i = 1, NEFMat%NumSpace_row
          NEFMat%Size_vector_alpha_row(i) = Size(initial_map_matrix_row(i),'alpha')
          NEFMat%Size_vector_beta_row(i) = Size(initial_map_matrix_row(i),'beta')
          do j = 1, NEFMat%NumSpace_col
            if(i.eq.1) then
              NEFMat%Size_vector_alpha_col(j) = Size(active_map_matrix_col(j),'alpha')
              NEFMat%Size_vector_beta_col(j) = Size(active_map_matrix_col(j),'beta')
            end if
            if (same_vec.and.i.eq.j) then
              nef_bidx = NEFMat%blockIndex([1,i],[2,i])
              NEFMat%Zero_Flag(nef_bidx) = scalarIsZero
              if (.not.scalarIsZero) then
                call NEFMat%Main_matrix(nef_bidx)%init(&
                       Size(initial_map_matrix_row(i),'alpha'),&
                       Size(initial_map_matrix_row(i),'beta'),&
                       scalar=scalar)
              end if
            else
              nef_bidx = NEFMat%blockIndex([1,i],[2,j])
              NEFMat%Zero_Flag(nef_bidx) = scalarIsZero
              if (.not.scalarIsZero) then
                call NEFMat%Main_matrix(nef_bidx)%init(&
                       Size(initial_map_matrix_row(i),'alpha'),&
                       Size(initial_map_matrix_row(i),'beta'),&
                       scalar=scalar,&
                       nAlpha2=Size(active_map_matrix_col(j),'alpha'),&
                       nBeta2=Size(active_map_matrix_col(j),'beta'))
              end if
            end if
          end do
        end do
      end subroutine mqc_nef_matrix_initialize
!
!
!     PROCEDURE MQC_NEF_Matrix_Import_Integral
!
!>    \brief <b> MQC_NEF_Matrix_Import_Integral is a subroutine that import the values 
!>    from a MQC SCF Integral object to an initiallized NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_Matrix_Import_Integral is a subroutine that import the values 
!>    from a MQC SCF Integral object to an initiallized NEF Matrix object
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to initialize.
!>    \endverbatim
!>
!>    \param[in] integral
!>    \verbatim
!>        integral is a type(mqc_scf_integral)
!>        A MQC SCF Integral object containing the NEF Matrix elements in the 
!>        original order (in general AO or MO)
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_matrix_import_integral(NEFMat,integral)
!
        implicit none
      integer(kind=int64) :: nef_bidx
        class(mqc_nef_matrix),intent(inOut)::NEFMat
        type(mqc_scf_integral)::integral
        type(mqc_scf_integral)::tmp_integral
        type(mqc_vector)::tmp_alpha,tmp_beta
        integer::i,j,row,col,rowSize
!
        do row = 1, NEFMat%NumSpace_row
          tmp_alpha = NEFMat%Map_Matrix_Row(row)%getBlock('alpha')
          tmp_beta = NEFMat%Map_Matrix_Row(row)%getBlock('beta')
          rowSize = size(tmp_alpha)+size(tmp_beta)
          tmp_integral = &
             integral%orbitals('useStrings',&
              alphaOrbsIn = tmp_alpha%getArrayI(),&
               BetaOrbsIn = tmp_beta%getArrayI(),axis=1)
          do col = 1, NEFMat%NumSpace_col
            tmp_alpha = NEFMat%Map_Matrix_Col(col)%getBlock('alpha')
            tmp_beta = NEFMat%Map_Matrix_Col(col)%getBlock('beta')
            nef_bidx = NEFMat%blockIndex([1,row],[2,col])
            NEFMat%Main_matrix(nef_bidx)=&
                tmp_integral%orbitals('useStrings',&
              alphaOrbsIn = tmp_alpha%getArrayI(),&
               BetaOrbsIn = tmp_beta%getArrayI(),axis=2)
            NEFMat%Zero_Flag(nef_bidx) = &
             (rowSize.eq.0.or.(size(tmp_alpha)+size(tmp_beta)).eq.0)
          end do
        end do
!
      end subroutine mqc_nef_matrix_import_integral
!
!
!
!     PROCEDURE MQC_NEF_Matrix_Export_Integral
!
!>    \brief <b> MQC_NEF_Matrix_Export_Integral is a subroutine that export the values 
!>    to a MQC SCF Integral object from a NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_Matrix_Export_Integral is a subroutine that export the values 
!>    to a MQC SCF Integral object from a NEF Matrix object
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix.
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_export_integral(NEFMat) result(integral)
!
        implicit none
        class(mqc_nef_matrix),intent(inOut)::NEFMat
        type(mqc_scf_integral)::integral
        type(mqc_scf_integral)::tmp_integral
        type(mqc_vector)::tmp_alpha,tmp_beta
        integer::id,ia,ja,ib,jb,row,col,maxColAlpha,maxColBeta,maxRowAlpha,maxRowBeta
!
        maxColAlpha = 0
        maxColBeta  = 0
        maxRowAlpha = 0
        maxRowBeta  = 0
        do row = 1, NEFMat%NumSpace_row
          maxRowAlpha = maxRowAlpha + NEFMat%Size_vector_alpha_row(row)
          maxRowBeta = maxRowBeta + NEFMat%Size_vector_beta_row(row)
        end do
        do col = 1, NEFMat%NumSpace_col
          maxColAlpha = maxColAlpha + NEFMat%Size_vector_alpha_col(col)
          maxColBeta = maxColBeta + NEFMat%Size_vector_beta_col(col)
        end do
!
        call integral%init(maxRowAlpha,maxRowBeta,nAlpha2=maxColAlpha,nBeta2=maxColBeta)
        do row = 1, NEFMat%NumSpace_row
          do col = 1, NEFMat%NumSpace_col
            id = NEFMat%blockIndex([1,row],[2,col])
            do ia = 1, NEFMat%Size_vector_alpha_row(row)
              do ja = 1, NEFMat%Size_vector_alpha_col(col)
                call integral%put(NEFMat%Main_matrix(id)%at(ia,ja,'alpha'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'alpha'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'alpha'),'alpha')
              end do
              do ja = 1, NEFMat%Size_vector_beta_col(col)
                call integral%put(NEFMat%Main_matrix(id)%at(ia,ja,'betaalpha'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'alpha'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'beta'),'betaalpha')
              end do
            end do
            do ia = 1, NEFMat%Size_vector_beta_row(row)
              do ja = 1, NEFMat%Size_vector_alpha_col(col)
                call integral%put(NEFMat%Main_matrix(id)%at(ia,ja,'alphabeta'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'beta'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'alpha'),'alphabeta')
              end do
              do ja = 1, NEFMat%Size_vector_beta_col(col)
                call integral%put(NEFMat%Main_matrix(id)%at(ia,ja,'beta'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'beta'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'beta'),'beta')
              end do
            end do
          end do
        end do
!
      end function mqc_nef_matrix_export_integral
!
!
!     PROCEDURE MQC_NEF_Vector_Import_Vector
!
!>    \brief <b> MQC_NEF_Vector_Import_Vector is a subroutine that import the values 
!>    from a MQC Vector object to an initiallized NEF Vector object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    MQC_NEF_Vector_Import_Vector is a subroutine that import the values 
!>    from a MQC Vector object to an initiallized NEF Vector object
!>    \verbatim
!>       
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVec
!>    \verbatim
!>        NEFVec is Class(MQC_NEF_Vector)
!>        The MQC NEF Vector to initialize.
!>    \endverbatim
!>
!>    \param[in] vector
!>    \verbatim
!>        vector is a type(mqc_vector)
!>        A MQC Vector object containing the NEF Vector elements in the 
!>        original order (in general AO or MO)
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_vector_import_vector(NEFVec,vector)
!
        implicit none
        class(mqc_nef_vector),intent(inOut)::NEFVec
        type(mqc_vector)::vector
        type(mqc_vector)::tmp_alpha,tmp_beta
        integer::row,el
        type(mqc_scalar)::zero
!
        zero = 0.0
        do row = 1, NEFVec%NumSpace
          tmp_alpha = NEFVec%Map_Matrix(row)%getBlock('alpha')
          tmp_beta = NEFVec%Map_Matrix(row)%getBlock('beta')
          call NEFVec%Main_Vector(row)%init(size(tmp_alpha),zero,size(tmp_beta))
          do el = 1, size(tmp_alpha)
             call NEFVec%Main_Vector(row)%put(vector%at(tmp_alpha%at(el)),el,'alpha')
          end do
          do el = 1, size(tmp_beta)
             call NEFVec%Main_Vector(row)%put(vector%at(tmp_beta%at(el)),el,'beta')
          end do
          NEFVec%Zero_Flag(row) = (size(tmp_alpha)+size(tmp_beta)).eq.0
        end do
!
      end subroutine mqc_nef_vector_import_vector
!     PROCEDURE MQC_NEFMatrix_NEFMatrix_Multiply
!       
!>    \brief <b> MQC_NEFMatrix_NEFMatrix_Multiply is used to multiply two MQC NEF Matrix
!>    type variables and returns an MQC NEF Matrix</b> 
!
!>    \par Purpose:
!     =============
!>      
!>    \verbatim
!>      
!>    MQC_NEFMatrix_NEFMatrix_Multiply is used to multiply two MQC NEF Matrix type 
!>    variables and returns an MQC NEF Matrix</b> 
!>
!>    \endverbatim
!       
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMAtA is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix type variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix)
!>        The second MQC NEF Matrix variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_nefmatrix_multiply(NEFMatrixA,NEFMatrixB) result(NEFMatrixOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatrixA,NEFMatrixB
      type(mqc_nef_matrix)::NEFMatrixOut
      integer(kind=int64)::rowA,colB,colA,final_index,nef_bidx,nef_bidx2
!
      if(NEFMatrixA%NumSpace_col.ne.NEFMatrixB%NumSpace_row) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_nefmatrix_multiply',6,&
        'NEFMatA%NumSpace_col',NEFMatrixA%NumSpace_col,&
        'NEFMatB%NumSpace_row',NEFMatrixB%NumSpace_row)

      do colA = 1, NEFMatrixA%NumSpace_col
        if(NEFMatrixA%Size_vector_alpha_col(colA)&
           .ne.NEFMatrixB%Size_vector_alpha_row(colA).or. &
           NEFMatrixA%Size_vector_beta_col(colA)&
           .ne.NEFMatrixB%Size_vector_beta_row(colA)) &
          call mqc_error_i('NEF subspace dimensions are &
          &not conformable in mqc_nefmatrix_nefmatrix_multiply',6,&
          'Sub-space index (colA)', colA,&
          'NEFMatA%Size_vector_alpha_col(colA)',NEFMatrixA%Size_vector_alpha_col(colA),&
          'NEFMatB%Size_vector_alpha_row(colA)',NEFMatrixB%Size_vector_alpha_row(colA),&
          'NEFMatA%Size_vector_beta_col(colA)',NEFMatrixA%Size_vector_beta_col(colA),&
          'NEFMatB%Size_vector_beta_row(colA)',NEFMatrixB%Size_vector_beta_row(colA))
      end do

      call NEFMatrixOut%init(NEFMatrixA%Map_Matrix_Row,&
          initial_map_matrix_col=NEFMatrixB%Map_Matrix_Col)
      do rowA = 1, NEFMatrixA%NumSpace_row
        do colB = 1, NEFMatrixB%NumSpace_col
          final_index = NEFMatrixOut%blockIndex([1,rowA],[2,colB])
          do colA = 1, NEFMatrixA%NumSpace_col
            nef_bidx = NEFMatrixA%blockIndex([1,rowA],[2,colA])
            nef_bidx2 = NEFMatrixB%blockIndex([1,colA],[2,colB])
            if (NEFMatrixA%Zero_Flag(nef_bidx).or.NEFMatrixB%Zero_Flag(nef_bidx2)) cycle
            if (NEFMatrixOut%Zero_Flag(final_index)) then
              NEFMatrixOut%Main_Matrix(final_index) =  &
                matmul(NEFMatrixA%Main_Matrix(nef_bidx),&
                       NEFMatrixB%Main_Matrix(nef_bidx2))
              NEFMatrixOut%Zero_Flag(final_index) = .false.
            else
              NEFMatrixOut%Main_Matrix(final_index) = &
                 NEFMatrixOut%Main_Matrix(final_index) + &
                matmul(NEFMatrixA%Main_Matrix(nef_bidx),&
                       NEFMatrixB%Main_Matrix(nef_bidx2))
            end if
          end do
        end do
      end do
!
      end function mqc_nefmatrix_nefmatrix_multiply
!
!     PROCEDURE MQC_NEFMatrix_NEFVector_Multiply
!
!>    \brief <b> MQC_NEFMatrix_NEFVector_Multiply is used to multiply a MQC NEF Matrix
!>    and a MQC NEF Vector type variable and returns an MQC NEF Vector</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_NEFVector_Multiply is used to multiply a MQC NEF Matrix type
!>    variable by a MQC NEF Vector type variable and returns an MQC NEF Vector
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatrixA
!>    \verbatim
!>        NEFMatrixA is type(mqc_nef_matrix)
!>        The MQC NEF Matrix type variable to multiply.
!>    \endverbatim
!>
!>    \param[in] NEFVectorB
!>    \verbatim
!>        NEFVectorB is type(mqc_nef_vector)
!>        The MQC NEF Vector type variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_nefvector_multiply(NEFMatrixA,NEFVectorB) result(NEFVectorOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatrixA
      type(mqc_nef_vector),intent(in)::NEFVectorB
      type(mqc_nef_vector)::NEFVectorOut
      integer(kind=int64)::rowA,colA,nef_bidx
!
      if(NEFMatrixA%NumSpace_col.ne.NEFVectorB%NumSpace) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_nefvector_multiply',6,&
        'NEFMatA%NumSpace_col',NEFMatrixA%NumSpace_col,&
        'NEFVecB%NumSpace',NEFVectorB%NumSpace)

      do colA = 1, NEFMatrixA%NumSpace_col
        if(NEFMatrixA%Size_vector_alpha_col(colA)&
           .ne.NEFVectorB%Size_vector_alpha(colA).or. &
           NEFMatrixA%Size_vector_beta_col(colA)&
           .ne.NEFVectorB%Size_vector_beta(colA)) &
          call mqc_error_i('NEF subspace dimensions are &
          &not conformable in mqc_nefmatrix_nefvector_multiply',6,&
          'Sub-space index (colA)', colA,&
          'NEFMatA%Size_vector_alpha_col(colA)',NEFMatrixA%Size_vector_alpha_col(colA),&
          'NEFVecB%Size_vector_alpha(colA)',NEFVectorB%Size_vector_alpha(colA),&
          'NEFMatA%Size_vector_beta_col(colA)',NEFMatrixA%Size_vector_beta_col(colA),&
          'NEFVecB%Size_vector_beta(colA)',NEFVectorB%Size_vector_beta(colA))
      end do

      call NEFVectorOut%init(NEFMatrixA%Map_Matrix_Row)
      do rowA = 1, NEFMatrixA%NumSpace_row
        do colA = 1, NEFMatrixA%NumSpace_col
          nef_bidx = NEFMatrixA%blockIndex([1,rowA],[2,colA])
          if (NEFMatrixA%Zero_Flag(nef_bidx).or.NEFVectorB%Zero_Flag(colA)) cycle
          if (NEFVectorOut%Zero_Flag(rowA)) then
            NEFVectorOut%Main_Vector(rowA) = &
              NEFMatrixA%Main_Matrix(nef_bidx).dot.NEFVectorB%Main_Vector(colA)
            NEFVectorOut%Zero_Flag(rowA) = .false.
          else
            NEFVectorOut%Main_Vector(rowA) = NEFVectorOut%Main_Vector(rowA) + &
              (NEFMatrixA%Main_Matrix(nef_bidx).dot.NEFVectorB%Main_Vector(colA))
          end if
        end do
      end do
!
      end function mqc_nefmatrix_nefvector_multiply
!
!     PROCEDURE MQC_NEFMatrix_Integral_Multiply
!       
!>    \brief <b> MQC_NEFMatrix_Integral_Multiply is used to multiply a MQC NEF Matrix
!>    and MQC SCF Integral type variables and returns an MQC NEF Matrix</b> 
!
!>    \par Purpose:
!     =============
!>      
!>    \verbatim
!>      
!>    MQC_NEFMatrix_Integral_Multiply is used to multiply a MQC NEF Matrix
!>    and MQC SCF Integral type variables and returns an MQC NEF Matrix
!>
!>    \endverbatim
!       
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMAt is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix type variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] IntegralB
!>    \verbatim
!>        Integral is type(mqc_scf_integral)
!>        The second MQC SCF Integral variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_integral_multiply(NEFMatA,IntegralB) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA
      type(mqc_scf_integral),intent(in)::IntegralB
      type(mqc_nef_matrix)::NEFMatOut
      type(mqc_scf_eigenvalues),dimension(:),allocatable::IntegralMap
      integer(kind=int64)::row,i
!
      if(NEFMatA%NumSpace_col.ne.1) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_integral_multiply',6,&
        'NEFMatA%NumSpace_col',NEFMatA%NumSpace_col)
     
      if(NEFMatA%Size_vector_alpha_col(1)&
         .ne.size(IntegralB,'alpha',1).or. &
         NEFMatA%Size_vector_beta_col(1)&
         .ne.size(IntegralB,'beta',1)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nefmatrix_integral_multiply',6,&
        'NEFMatA%Size_vector_alpha_col',NEFMatA%Size_vector_alpha_col(1),&
        'size(IntegralB,alpha,1)',size(IntegralB,'alpha',1),&
        'NEFMatA%Size_vector_beta_col',NEFMatA%Size_vector_beta_col(1),&
        'size(IntegralB,beta,1)',size(IntegralB,'beta',1))

      allocate(IntegralMap(1))
      call IntegralMap(1)%init(0,0)
      !! HERE WE SHOULD CHECK IF THE INTEGRAL IS SPACE OR NOT 
      do i = 1, size(IntegralB,'alpha',2)
        call IntegralMap(1)%push(i,'alpha')
      end do
      do i = 1, size(IntegralB,'beta',2)
        call IntegralMap(1)%push(i,'beta')
      end do


      call NEFMatOut%init(NEFMatA%Map_Matrix_Row,&
          initial_map_matrix_col=IntegralMap)
      do row = 1, NEFMatA%NumSpace_row
        if (.not.NEFMatA%Zero_Flag(row)) then
          NEFMatOut%Main_Matrix(row) =  matmul(NEFMatA%Main_Matrix(row),IntegralB)
          NEFMatOut%Zero_Flag(row) = .false.
        end if
      end do
!
      end function mqc_nefmatrix_integral_multiply
!       
!     PROCEDURE MQC_Integral_NEFMatrix_Multiply
!       
!>    \brief <b> MQC_Integral_NEFMatrix_Multiply is used to multiply a MQC NEF Matrix
!>    and MQC SCF Integral type variables and returns an MQC NEF Matrix</b> 
!
!>    \par Purpose:
!     =============
!>      
!>    \verbatim
!>      
!>    MQC_Integral_NEFMatrix_Multiply is used to multiply a MQC NEF Matrix
!>    and MQC SCF Integral type variables and returns an MQC NEF Matrix
!>
!>    \endverbatim
!       
!     Arguments:
!     ==========
!>    \param[in] IntegralA
!>    \verbatim
!>        Integral is type(mqc_scf_integral)
!>        The second MQC SCF Integral variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMAt is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix type variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_integral_nefmatrix_multiply(IntegralA,NEFMatB) result(NEFMatOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::IntegralA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_matrix)::NEFMatOut
      type(mqc_scf_eigenvalues),dimension(:),allocatable::IntegralMap
      integer(kind=int64)::col,i
!
      if(NEFMatB%NumSpace_row.ne.1) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_integral_nefmatrix_multiply',6,&
        'NEFMatA%NumSpace_row',NEFMatB%NumSpace_row)
 
      if(NEFMatB%Size_vector_alpha_row(1)&
         .ne.size(IntegralA,'alpha',2).or. & 
         NEFMatB%Size_vector_beta_row(1)&
         .ne.size(IntegralA,'beta',2)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_integral_nefmatrix_multiply',6,&
        'NEFMatB%Size_vector_alpha_row',NEFMatB%Size_vector_alpha_row(1),&
        'size(IntegralA,alpha,2)',size(IntegralA,'alpha',2),&
        'NEFMatB%Size_vector_beta_row',NEFMatB%Size_vector_beta_row(1),&
        'size(IntegralA,beta,2)',size(IntegralA,'beta',2))

      allocate(IntegralMap(1))
      call IntegralMap(1)%init(0,0)
      !! HERE WE SHOULD CHECK IF THE INTEGRAL IS SPACE OR NOT 
      do i = 1, size(IntegralA,'alpha',1)
        call IntegralMap(1)%push(i,'alpha')
      end do
      do i = 1, size(IntegralA,'beta',1)
        call IntegralMap(1)%push(i,'beta')
      end do

      call NEFMatOut%init(IntegralMap,&
          initial_map_matrix_col=NEFMatB%Map_Matrix_Col)
      do col = 1, NEFMatB%NumSpace_col
        if (.not.NEFMatB%Zero_Flag(col)) then
          NEFMatOut%Main_Matrix(col) = matmul(IntegralA,NEFMatB%Main_Matrix(col))
          NEFMatOut%Zero_Flag(col) = .false.
        end if
      end do
!
      end function mqc_integral_nefmatrix_multiply
!
!
!     PROCEDURE MQC_NEF_Vector_Sum
!
!>    \brief <b> MQC_NEF_Vector_Sum is used to sum two MQC NEF Vector type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Sum is used to sum two MQC NEF vector type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVecA
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>        The first MQC NEF Vector variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFVecB
!>    \verbatim               
!>        NEFVecB is type(mqc_nef_vector)
!>        The first MQC NEF Vector variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_sum(NEFVecA,NEFVecB) result(NEFVecOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecA,NEFVecB
      type(mqc_nef_vector)::NEFVecOut
      integer::col
!
      if(NEFVecA%NumSpace.ne.NEFVecB%NumSpace) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_vector_sum',6,&
        'NEFVecA%NumSpace',NEFVecA%NumSpace,&
        'NEFVecB%NumSpace',NEFVecB%NumSpace)
!
      do col = 1, NEFVecA%NumSpace
        if (NEFVecA%Size_vector_alpha(col).ne.NEFVecB%Size_vector_alpha(col) .or. &
                NEFVecA%Size_vector_beta(col) .ne.NEFVecB%Size_vector_beta(col)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_vector_sum',6,&
          'Sub-space index (col)', col,&
          'NEFVecA%Size_vector_alpha(col)',NEFVecA%Size_vector_alpha(col),&
          'NEFVecB%Size_vector_alpha(col)',NEFVecB%Size_vector_alpha(col),&
          'NEFVecA%Size_vector_beta(col)',NEFVecA%Size_vector_beta(col),&
          'NEFVecB%Size_vector_beta(col)',NEFVecB%Size_vector_beta(col))
      end do
!
      call NEFVecOut%init(NEFVecA%Map_Matrix)
      do col = 1, NEFVecB%NumSpace
        if (NEFVecA%Zero_Flag(col).and.NEFVecB%Zero_Flag(col)) then
          continue
        else if (NEFVecA%Zero_Flag(col)) then
          NEFVecOut%Main_Vector(col) = NEFVecB%Main_Vector(col)
          NEFVecOut%Zero_Flag(col) = .false.
        else if (NEFVecB%Zero_Flag(col)) then
          NEFVecOut%Main_Vector(col) = NEFVecA%Main_Vector(col)
          NEFVecOut%Zero_Flag(col) = .false.
        else
          NEFVecOut%Main_Vector(col) = NEFVecA%Main_Vector(col) + NEFVecB%Main_Vector(col)
          NEFVecOut%Zero_Flag(col) = .false.
        end if
      end do
!
      end function mqc_nef_vector_sum
!
!
!     PROCEDURE MQC_NEFMatrix_Sum
!
!>    \brief <b> MQC_NEFMatrix_Sum is used to sum two MQC NEF matrix type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Sum is used to sum two MQC NEF matrix type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim               
!>        NEFMatB is type(mqc_nef_matrix)
!>        The second MQC NEF Matrix variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_sum(NEFMatA,NEFMatB) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA,NEFMatB
      type(mqc_nef_matrix)::NEFMatOut
      integer::col,row,id
!
      if(NEFMatA%NumSpace_col.ne.NEFMatB%NumSpace_col.or.&
         NEFMatA%NumSpace_row.ne.NEFMatB%NumSpace_row) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_sum',6,&
        'NEFMatA%NumSpace_col',NEFMatA%NumSpace_col,&
        'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col,&
        'NEFMatA%NumSpace_row',NEFMatA%NumSpace_row,&
        'NEFMatB%NumSpace_row',NEFMatB%NumSpace_row)

      do col = 1, NEFMatA%NumSpace_col
        if(NEFMatA%Size_vector_alpha_col(col)&
           .ne.NEFMatB%Size_vector_alpha_col(col).or. &
               NEFMatA%Size_vector_beta_col(col)  &
           .ne.NEFMatB%Size_vector_beta_col(col)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_sum',6,&
          'Sub-space index (col)', col,&
          'NEFMatA%Size_vector_alpha_col(col)',NEFMatA%Size_vector_alpha_col(col),&
          'NEFMatB%Size_vector_alpha_col(col)',NEFMatB%Size_vector_alpha_col(col),&
          'NEFMatA%Size_vector_beta_col(col)',NEFMatA%Size_vector_beta_col(col),&
          'NEFMatB%Size_vector_beta_col(col)',NEFMatB%Size_vector_beta_col(col))
      end do
      do row = 1, NEFMatA%NumSpace_row
        if(NEFMatA%Size_vector_alpha_row(row)&
           .ne.NEFMatB%Size_vector_alpha_row(row).or. &
               NEFMatA%Size_vector_beta_row(row)  &
           .ne.NEFMatB%Size_vector_beta_row(row)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_sum',6,&
          'Sub-space index (row)', row,&
          'NEFMatA%Size_vector_alpha_row(row)',NEFMatA%Size_vector_alpha_row(row),&
          'NEFMatB%Size_vector_alpha_row(row)',NEFMatB%Size_vector_alpha_row(row),&
          'NEFMatA%Size_vector_beta_row(row)',NEFMatA%Size_vector_beta_row(row),&
          'NEFMatB%Size_vector_beta_row(row)',NEFMatB%Size_vector_beta_row(row))
      end do

      call NEFMatOut%init(NEFMatA%Map_Matrix_Row,&
          initial_map_matrix_col=NEFMatB%Map_Matrix_Col)
      do row = 1, NEFMatA%NumSpace_row
        do col = 1, NEFMatB%NumSpace_col
          id = NEFMatOut%blockIndex([1,row],[2,col])
          if (NEFMatA%Zero_Flag(id).and.NEFMatB%Zero_Flag(id)) then
            continue
          else if (NEFMatA%Zero_Flag(id)) then
            NEFMatOut%Main_Matrix(id) = NEFMatB%Main_Matrix(id)
            NEFMatOut%Zero_Flag(id) = .false.
          else if (NEFMatB%Zero_Flag(id)) then
            NEFMatOut%Main_Matrix(id) = NEFMatA%Main_Matrix(id)
            NEFMatOut%Zero_Flag(id) = .false.
          else
            NEFMatOut%Main_Matrix(id) = NEFMatA%Main_Matrix(id)+NEFMatB%Main_Matrix(id)
            NEFMatOut%Zero_Flag(id) = .false.
          end if
        end do
      end do
!
      end function mqc_nefmatrix_sum
!
!
!     PROCEDURE MQC_NEF_R4_Sum
!
!>    \brief <b> MQC_NEF_R4_Sum is used to sum two MQC NEF R4 type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Sum is used to sum two MQC NEF R4 type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4)
!>        The first MQC NEF R4 variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4)
!>        The second MQC NEF R4 variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_sum(NEFR4A,NEFR4B) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      type(mqc_nef_r4)::NEFR4Out
      integer::id1,id2,id3,id4,id
!
      if(NEFR4A%NumSpace_1.ne.NEFR4B%NumSpace_1)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_sum',6,&
        'NEFR4A%NumSpace_1',NEFR4A%NumSpace_1,&
        'NEFR4B%NumSpace_1',NEFR4B%NumSpace_1)
      if(NEFR4A%NumSpace_2.ne.NEFR4B%NumSpace_2)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_sum',6,&
        'NEFR4A%NumSpace_2',NEFR4A%NumSpace_2,&
        'NEFR4B%NumSpace_2',NEFR4B%NumSpace_2)
      if(NEFR4A%NumSpace_3.ne.NEFR4B%NumSpace_3)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_sum',6,&
        'NEFR4A%NumSpace_3',NEFR4A%NumSpace_3,&
        'NEFR4B%NumSpace_3',NEFR4B%NumSpace_3)
      if(NEFR4A%NumSpace_4.ne.NEFR4B%NumSpace_4)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_sum',6,&
        'NEFR4A%NumSpace_4',NEFR4A%NumSpace_4,&
        'NEFR4B%NumSpace_4',NEFR4B%NumSpace_4)
     
      do id1 = 1, NEFR4A%NumSpace_1
        if(NEFR4A%Size_vector_alpha_1(id1)&
           .ne.NEFR4B%Size_vector_alpha_1(id1).or.&
           NEFR4A%Size_vector_beta_1(id1)&
           .ne.NEFR4B%Size_vector_beta_1(id1))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_sum',6,&
          'Sub-space index (id1)', id1,&
          'NEFR4A%Size_vector_alpha_1(id1)',NEFR4A%Size_vector_alpha_1(id1),&
          'NEFR4B%Size_vector_alpha_1(id1)',NEFR4B%Size_vector_alpha_1(id1),&
          'NEFR4A%Size_vector_beta_1(id1)',NEFR4A%Size_vector_beta_1(id1),&
          'NEFR4B%Size_vector_beta_1(id1)',NEFR4B%Size_vector_beta_1(id1))
      end do
      do id2 = 1, NEFR4A%NumSpace_2
        if(NEFR4A%Size_vector_alpha_2(id2)&
           .ne.NEFR4B%Size_vector_alpha_2(id2).or.&
           NEFR4A%Size_vector_beta_2(id2)&
           .ne.NEFR4B%Size_vector_beta_2(id2))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_sum',6,&
          'Sub-space index (id2)', id2,&
          'NEFR4A%Size_vector_alpha_2(id2)',NEFR4A%Size_vector_alpha_2(id2),&
          'NEFR4B%Size_vector_alpha_2(id2)',NEFR4B%Size_vector_alpha_2(id2),&
           'NEFR4A%Size_vector_beta_2(id2)', NEFR4A%Size_vector_beta_2(id2),&
           'NEFR4B%Size_vector_beta_2(id2)', NEFR4B%Size_vector_beta_2(id2))
      end do
      do id3 = 1, NEFR4A%NumSpace_3
        if(NEFR4A%Size_vector_alpha_3(id3)&
           .ne.NEFR4B%Size_vector_alpha_3(id3).or.&
           NEFR4A%Size_vector_beta_3(id3)&
           .ne.NEFR4B%Size_vector_beta_3(id3))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_sum',6,&
          'Sub-space index (id3)', id3,&
          'NEFR4A%Size_vector_alpha_3(id3)',NEFR4A%Size_vector_alpha_3(id3),&
          'NEFR4B%Size_vector_alpha_3(id3)',NEFR4B%Size_vector_alpha_3(id3),&
           'NEFR4A%Size_vector_beta_3(id3)', NEFR4A%Size_vector_beta_3(id3),&
           'NEFR4B%Size_vector_beta_3(id3)', NEFR4B%Size_vector_beta_3(id3))
      end do
      do id4 = 1, NEFR4A%NumSpace_4
        if(NEFR4A%Size_vector_alpha_4(id4)&
           .ne.NEFR4B%Size_vector_alpha_4(id4).or.&
           NEFR4A%Size_vector_beta_4(id4)&
           .ne.NEFR4B%Size_vector_beta_4(id4))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_sum',6,&
          'Sub-space index (id4)', id4,&
          'NEFR4A%Size_vector_alpha_4(id4)',NEFR4A%Size_vector_alpha_4(id4),&
          'NEFR4B%Size_vector_alpha_4(id4)',NEFR4B%Size_vector_alpha_4(id4),&
           'NEFR4A%Size_vector_beta_4(id4)', NEFR4A%Size_vector_beta_4(id4),&
           'NEFR4B%Size_vector_beta_4(id4)', NEFR4B%Size_vector_beta_4(id4))
      end do


      NEFR4Out = NEFR4A
      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              id = NEFR4Out%blockIndex([4,id4],[3,id3],[2,id2],[1,id1])
              if (NEFR4A%Zero_Flag(id).and.NEFR4B%Zero_Flag(id)) then
                continue
              else if (NEFR4A%Zero_Flag(id)) then
                NEFR4Out%Main_R4(id) = NEFR4B%Main_R4(id)
                NEFR4Out%Zero_Flag(id) = .false.
              else if (NEFR4B%Zero_Flag(id)) then
                NEFR4Out%Main_R4(id) = NEFR4A%Main_R4(id)
                NEFR4Out%Zero_Flag(id) = .false.
              else
                NEFR4Out%Main_R4(id) = NEFR4A%Main_R4(id) + NEFR4B%Main_R4(id)
                NEFR4Out%Zero_Flag(id) = .false.
              end if
           end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_sum
!
!     PROCEDURE MQC_NEFMatrix_Conjugate
!
!>    \brief <b> MQC_NEFMatrix_Conjugate is used to compute the conjugate 
!>    of a MQC NEF Matrix type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Conjugate is used to compute the conjugate of 
!>    a MQC NEF Matrix type variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix type variable to conjugate.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_NEFMatrix_conjugate(NEFMat) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::row,col,nef_bidx,nef_bidx2
!
      call NEFMatOut%init(NEFMat%Map_Matrix_row,&
           initial_map_matrix_col=NEFMat%Map_Matrix_col)
      do row = 1, NEFMat%NumSpace_row
        do col = 1, NEFMat%NumSpace_col
          nef_bidx = NEFMatOut%blockIndex([1,row],[2,col])
          nef_bidx2 = NEFMat%blockIndex([1,row],[2,col])
          if (.not.NEFMat%Zero_Flag(nef_bidx2)) then
            NEFMatOut%Main_Matrix(nef_bidx) = &
              conjg(NEFMat%Main_Matrix(nef_bidx2))
            NEFMatOut%Zero_Flag(nef_bidx) = .false.
          end if
        end do
      end do
!
      end function mqc_NEFMatrix_conjugate
!
!     PROCEDURE MQC_NEFMatrix_Conjugate_Transpose
!
!>    \brief <b> MQC_NEFMatrix_Conjugate_Transpose is used to compute the conjugate 
!>    transpose of a MQC NEF Matrix type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Conjugate_Transpose is used to compute the conjugate transpose of 
!>    a MQC NEF Matrix type variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix type variable to conjugate transpose.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_NEFMatrix_conjugate_transpose(NEFMat) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::row,col,nef_bidx,nef_bidx2
!
      call NEFMatOut%init(NEFMat%Map_Matrix_col,&
           initial_map_matrix_col=NEFMat%Map_Matrix_row)

      do row = 1, NEFMat%NumSpace_row
        do col = 1, NEFMat%NumSpace_col
          nef_bidx = NEFMatOut%blockIndex([2,row],[1,col])
          nef_bidx2 = NEFMat%blockIndex([1,row],[2,col])
          if (.not.NEFMat%Zero_Flag(nef_bidx2)) then
            NEFMatOut%Main_Matrix(nef_bidx) = &
              dagger(NEFMat%Main_Matrix(nef_bidx2))
            NEFMatOut%Zero_Flag(nef_bidx) = .false.
          end if
        end do
      end do
!
      end function mqc_NEFMatrix_conjugate_transpose
!
!     PROCEDURE MQC_NEFR4_Conjugate_Transpose
!
!>    \brief <b> MQC_NEFR4_Conjugate_Transpose is used to compute the conjugate 
!>    transpose of a MQC NEF R4 type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFR4_Conjugate_Transpose is used to compute the conjugate transpose of 
!>    a MQC NEF R4 type variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is type(mqc_nef_r4)
!>        The MQC NEF R4 type variable to conjugate transpose.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_NEFR4_conjugate_transpose(NEFR4) result(NEFR4Out)
!
      implicit none
      integer(kind=int64) :: nef_bidx, nef_bidx2
      type(mqc_nef_r4),intent(in)::NEFR4
      type(mqc_nef_r4)::NEFR4Out
      integer(kind=int64)::id1,id2,id3,id4,N1234
!
      NEFR4Out%NumSpace_1 = NEFR4%NumSpace_4
      NEFR4Out%NumSpace_2 = NEFR4%NumSpace_3
      NEFR4Out%NumSpace_3 = NEFR4%NumSpace_2
      NEFR4Out%NumSpace_4 = NEFR4%NumSpace_1
      N1234 = NEFR4Out%NumSpace_1*NEFR4Out%NumSpace_2*NEFR4Out%NumSpace_3*NEFR4Out%NumSpace_4
      allocate(NEFR4Out%Main_R4(N1234))
      allocate(NEFR4Out%Zero_Flag(N1234))
      NEFR4Out%Zero_Flag = .true.
      NEFR4Out%Map_Matrix_1 = NEFR4%Map_Matrix_4
      NEFR4Out%Map_Matrix_2 = NEFR4%Map_Matrix_3
      NEFR4Out%Map_Matrix_3 = NEFR4%Map_Matrix_2
      NEFR4Out%Map_Matrix_4 = NEFR4%Map_Matrix_1
      NEFR4Out%Size_vector_alpha_1 = NEFR4%Size_vector_alpha_4
      NEFR4Out%Size_vector_alpha_2 = NEFR4%Size_vector_alpha_3
      NEFR4Out%Size_vector_alpha_3 = NEFR4%Size_vector_alpha_2
      NEFR4Out%Size_vector_alpha_4 = NEFR4%Size_vector_alpha_1
      NEFR4Out%Size_vector_beta_1 = NEFR4%Size_vector_beta_4
      NEFR4Out%Size_vector_beta_2 = NEFR4%Size_vector_beta_3
      NEFR4Out%Size_vector_beta_3 = NEFR4%Size_vector_beta_2
      NEFR4Out%Size_vector_beta_4 = NEFR4%Size_vector_beta_1

      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              nef_bidx = NEFR4Out%blockIndex([4,id4],[3,id3],[2,id2],[1,id1])
              nef_bidx2 = NEFR4%blockIndex([1,id4],[2,id3],[3,id2],[4,id1])
              if (.not.NEFR4%Zero_Flag(nef_bidx2)) then
                NEFR4Out%Main_R4(nef_bidx) = dagger(NEFR4%Main_R4(nef_bidx2))
                NEFR4Out%Zero_Flag(nef_bidx) = .false.
              end if
            end do
          end do
        end do
      end do
!
      end function mqc_NEFR4_conjugate_transpose
!
!
!     PROCEDURE MQC_NEFMatrix_Transpose
!
!>    \brief <b> MQC_NEFMatrix_Transpose is used to compute the 
!>    transpose of a MQC NEF Matrix type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Transpose is used to compute the transpose of 
!>    a MQC NEF Matrix type variable.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix type variable to transpose.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_NEFMatrix_transpose(NEFMat) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::row,col,nef_bidx,nef_bidx2
!
      call NEFMatOut%init(NEFMat%Map_Matrix_col,&
           initial_map_matrix_col=NEFMat%Map_Matrix_row)

      do row = 1, NEFMat%NumSpace_row
        do col = 1, NEFMat%NumSpace_col
          nef_bidx = NEFMatOut%blockIndex([2,row],[1,col])
          nef_bidx2 = NEFMat%blockIndex([1,row],[2,col])
          if (.not.NEFMat%Zero_Flag(nef_bidx2)) then
            NEFMatOut%Main_Matrix(nef_bidx) = &
              transpose(NEFMat%Main_Matrix(nef_bidx2))
            NEFMatOut%Zero_Flag(nef_bidx) = .false.
          end if
        end do
      end do
!
      end function mqc_NEFMatrix_transpose
!
!     
!     PROCEDURE MQC_NEF_Scalar_Vector_Multiply
!       
!>    \brief <b> MQC_NEF_Scalar_Vector_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Vector type variable and returns an MQC NEF Vector </b> 
!       
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_NEF_Scalar_Vector_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Vector type variable and returns an MQC NEF Vector.
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] ScalarA
!>    \verbatim
!>        ScalarB class(*)
!>        The variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFVecB
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>        The MQC NEF Vector variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_scalar_nef_vector_multiply(scalarA,NEFVecB) result(NEFVecOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecB
      class(*),intent(in)::scalarA
      type(mqc_scalar)::scalar
      type(mqc_nef_vector)::NEFVecOut
      integer::row
!
      select type(scalarA)
      type is (integer)
        scalar = scalarA
      type is (real)
        scalar = scalarA
      type is (complex)
        scalar = scalarA
      type is (mqc_scalar)
        scalar = scalarA
      class default
        call mqc_error_I('scalarIn type not defined in mqc_NEF_Vector_scalar_multiply',6)
      end select
!
      NEFVecOut = NEFVecB
      if (abs(scalar).eq.0.0) then
        NEFVecOut%Zero_Flag = .true.
      else
        do row = 1, NEFVecB%NumSpace
          if (.not.NEFVecB%Zero_Flag(row)) &
            NEFVecOut%Main_vector(row) = NEFVecB%Main_vector(row)*scalar
        end do
      end if
!
      end function mqc_scalar_nef_vector_multiply
!
!     
!     PROCEDURE MQC_NEF_Vector_Scalar_Multiply
!       
!>    \brief <b> MQC_NEF_Vector_Scalar_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Vector type variable and returns an MQC NEF Vector </b> 
!       
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_NEF_Vector_Scalar_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Vector type variable and returns an MQC NEF Vector.
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] NEFVecA
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>        The MQC NEF Vector variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] ScalarB
!>    \verbatim
!>        ScalarB class(*)
!>        The variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_scalar_multiply(NEFVecA,scalarB) result(NEFVecOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecA
      class(*),intent(in)::scalarB
      type(mqc_scalar)::scalar
      type(mqc_nef_vector)::NEFVecOut
      integer::row
!
      select type(scalarB)
      type is (integer)
        scalar = scalarB
      type is (real)
        scalar = scalarB
      type is (complex)
        scalar = scalarB
      type is (mqc_scalar)
        scalar = scalarB
      class default
        call mqc_error_I('scalarIn type not defined in mqc_scalar_NEF_Vector_multiply',6)
      end select
!
      NEFVecOut = NEFVecA
      if (abs(scalar).eq.0.0) then
        NEFVecOut%Zero_Flag = .true.
      else
        do row = 1, NEFVecA%NumSpace
          if (.not.NEFVecA%Zero_Flag(row)) &
            NEFVecOut%Main_vector(row) = scalar*NEFVecA%Main_vector(row)
        end do
      end if
!
      end function mqc_nef_vector_scalar_multiply
!
!     
!     PROCEDURE MQC_Scalar_NEFMatrix_Multiply
!       
!>    \brief <b> MQC_Scalar_NEFMatrix_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Matrix type variable and returns an MQC NEF Matrix</b> 
!       
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_Scalar_NEFMatrix_Multiply is used to multiply an MQC scalar with an MQC 
!>    NEF Matrix type variable and returns an MQC NEF Matrix
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] ScalarIn
!>    \verbatim
!>        ScalarIn class(*)
!>        The variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_scalar_nefmatrix_multiply(scalarIn,NEFMat) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      class(*),intent(in)::scalarIn
      type(mqc_scalar)::scalar
      type(mqc_nef_matrix)::NEFMatOut
      integer::col,row,id
!
      select type(scalarIn)
      type is (integer)
        scalar = scalarIn
      type is (real)
        scalar = scalarIn
      type is (complex)
        scalar = scalarIn
      type is (mqc_scalar)
        scalar = scalarIn
      class default
        call mqc_error_I('scalarIn type not defined in mqc_scalar_NEFMat_multiply',6)
      end select
!
      NEFMatOut = NEFMat
      if (abs(scalar).eq.0.0) then
        NEFMatOut%Zero_Flag = .true.
      else
        do row = 1, NEFMat%NumSpace_row
          do col = 1, NEFMat%NumSpace_col
            id=NEFMat%blockIndex([1,row],[2,col])
            if (.not.NEFMat%Zero_Flag(id)) &
              NEFMatOut%Main_matrix(id)=scalar*NEFMat%Main_matrix(id)
          end do
        end do
      end if
!
      end function mqc_scalar_nefmatrix_multiply
!
!     
!     PROCEDURE MQC_NEFMatrix_Scalar_Multiply
!       
!>    \brief <b> MQC_NEFMatrix_Scalar_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF Matrix type variable and returns an MQC NEF Matrix</b> 
!       
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_NEFMatrix_Scalar_Multiply is used to multiply an MQC scalar with an MQC 
!>    NEF Matrix type variable and returns an MQC NEF Matrix
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        The MQC NEF Matrix variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] ScalarB
!>    \verbatim
!>        ScalarB class(*)
!>        The variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_scalar_multiply(NEFMat,scalarB) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      class(*),intent(in)::scalarB
      type(mqc_scalar)::scalar
      type(mqc_nef_matrix)::NEFMatOut
      integer::col,row,id
!
      select type(scalarB)
      type is (integer)
        scalar = scalarB
      type is (real)
        scalar = scalarB
      type is (complex)
        scalar = scalarB
      type is (mqc_scalar)
        scalar = scalarB
      class default
        call mqc_error_I('scalarIn type not defined in mqc_scalar_NEFMat_multiply',6)
      end select
!
      NEFMatOut = NEFMat
      if (abs(scalar).eq.0.0) then
        NEFMatOut%Zero_Flag = .true.
      else
        do row = 1, NEFMat%NumSpace_row
          do col = 1, NEFMat%NumSpace_col
            id = NEFMat%blockIndex([1,row],[2,col])
            if (.not.NEFMat%Zero_Flag(id)) &
              NEFMatOut%Main_matrix(id)=NEFMat%Main_matrix(id)*scalar
          end do
        end do
      end if
!
      end function mqc_nefmatrix_scalar_multiply
!
!     
!     PROCEDURE MQC_Scalar_NEFR4_Multiply
!       
!>    \brief <b> MQC_Scalar_NEFR4_Multiply is used to multiply an MQC scalar with
!>    an MQC NEF R4 type variable and returns an MQC NEF R4</b> 
!     
!>    \par Purpose:
!     =============
!>    
!>    \verbatim
!>    
!>    MQC_Scalar_NEFR4_Multiply is used to multiply an MQC scalar with an MQC 
!>    NEF R4 type variable and returns an MQC NEF R4
!>    
!>    \endverbatim
!     
!     Arguments:
!     ==========
!>    \param[in] ScalarIn
!>    \verbatim
!>        ScalarIn class(*)
!>        The variable to multiply.
!>    \endverbatim
!>    
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is type(mqc_nef_r4)
!>        The MQC NEF R4 variable to multiply.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_scalar_nefr4_multiply(scalarIn,NEFR4) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4
      class(*),intent(in)::scalarIn
      type(mqc_scalar)::scalar
      type(mqc_nef_r4)::NEFR4Out
      integer::id,id1,id2,id3,id4
!
      select type(scalarIn)
      type is (integer)
        scalar = scalarIn
      type is (real)
        scalar = scalarIn
      type is (complex)
        scalar = scalarIn
      type is (mqc_scalar)
        scalar = scalarIn
      class default
        call mqc_error_I('scalarIn type not defined in mqc_scalar_NEFR4_multiply',6)
      end select
!
      NEFR4Out = NEFR4
      if (abs(scalar).eq.0.0) then
        NEFR4Out%Zero_Flag = .true.
      else
        do id1 = 1, NEFR4%NumSpace_1
          do id2 = 1, NEFR4%NumSpace_2
            do id3 = 1, NEFR4%NumSpace_3
              do id4 = 1, NEFR4%NumSpace_4
                id = NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
                if (.not.NEFR4%Zero_Flag(id)) &
                  NEFR4Out%Main_R4(id) = scalar*NEFR4Out%Main_R4(id)
              end do
            end do
          end do
        end do
      end if
!
      end function mqc_scalar_nefr4_multiply
      function mqc_nefr4_scalar_multiply(NEFR4,scalarB) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4
      class(*),intent(in)::scalarB
      type(mqc_scalar)::scalar
      type(mqc_nef_r4)::NEFR4Out
      integer::id,id1,id2,id3,id4
!
      select type(scalarB)
      type is (integer)
        scalar = scalarB
      type is (real)
        scalar = scalarB
      type is (complex)
        scalar = scalarB
      type is (mqc_scalar)
        scalar = scalarB
      class default
        call mqc_error_I('scalarIn type not defined in mqc_NEFR4_scalar_multiply',6)
      end select
!
      NEFR4Out = NEFR4
      if (abs(scalar).eq.0.0) then
        NEFR4Out%Zero_Flag = .true.
      else
        do id1 = 1, NEFR4%NumSpace_1
          do id2 = 1, NEFR4%NumSpace_2
            do id3 = 1, NEFR4%NumSpace_3
              do id4 = 1, NEFR4%NumSpace_4
                id = NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
                if (.not.NEFR4%Zero_Flag(id)) &
                  NEFR4Out%Main_R4(id) = NEFR4Out%Main_R4(id)*scalar
              end do
            end do
          end do
        end do
      end if
!
      end function mqc_nefr4_scalar_multiply
!
!
!
!     PROCEDURE MQC_NEF_Vector_Difference
!
!>    \brief <b> MQC_NEF_Vector_Difference is used to subtract two MQC NEF Vector type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Difference is used to subtract two MQC NEF vector type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVecA
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>        The first MQC NEF Vector variable to subtract.
!>    \endverbatim
!>    
!>    \param[in] NEFVecB
!>    \verbatim               
!>        NEFVecB is type(mqc_nef_vector)
!>        The second MQC NEF Vector variable to subtract.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      function mqc_nef_vector_difference(NEFVecA,NEFVecB) result(NEFVecOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecA,NEFVecB
      type(mqc_nef_vector)::NEFVecOut
      type(mqc_nef_vector)::tmpVectorAlpha,tmpVectorBeta
      integer::col
!
      if(NEFVecA%NumSpace.ne.NEFVecB%NumSpace) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_vector_difference',6,&
        'NEFVecA%NumSpace',NEFVecA%NumSpace,&
        'NEFVecB%NumSpace',NEFVecB%NumSpace)
!
        do col = 1, NEFVecA%NumSpace
          if (NEFVecA%Size_vector_alpha(col).ne.NEFVecB%Size_vector_alpha(col) .or. &
              NEFVecA%Size_vector_beta(col) .ne.NEFVecB%Size_vector_beta(col)) &
                call mqc_error_i('NEF subspace dimensions are not conformable &
                        &in mqc_nef_vector_difference',6,&
                        'Sub-space index (col)', col,&
                        'NEFVecA%Size_vector_alpha(col)',NEFVecA%Size_vector_alpha(col),&
                        'NEFVecB%Size_vector_alpha(col)',NEFVecB%Size_vector_alpha(col),&
                        'NEFVecA%Size_vector_beta(col)',NEFVecA%Size_vector_beta(col),&
                        'NEFVecB%Size_vector_beta(col)',NEFVecB%Size_vector_beta(col))      
        end do
!
      call NEFVecOut%init(NEFVecA%Map_Matrix)
      do col = 1, NEFVecB%NumSpace
        if (NEFVecA%Zero_Flag(col).and.NEFVecB%Zero_Flag(col)) then
          continue
        else if (NEFVecA%Zero_Flag(col)) then
          NEFVecOut%Main_Vector(col) = (-1.0)*NEFVecB%Main_Vector(col)
          NEFVecOut%Zero_Flag(col) = .false.
        else if (NEFVecB%Zero_Flag(col)) then
          NEFVecOut%Main_Vector(col) = NEFVecA%Main_Vector(col)
          NEFVecOut%Zero_Flag(col) = .false.
        else
          NEFVecOut%Main_Vector(col) = NEFVecA%Main_Vector(col) - NEFVecB%Main_Vector(col)
          NEFVecOut%Zero_Flag(col) = .false.
        end if
      end do
!
!
      end function mqc_nef_vector_difference
!

!
!     PROCEDURE MQC_NEFMatrix_Difference
!
!>    \brief <b> MQC_NEFMatrix_Difference is used to subtract two MQC NEF matrix type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Difference is used to subtract two MQC NEF matrix type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix)
!>        The first MQC NEF Matrix variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim               
!>        NEFMatB is type(mqc_nef_matrix)
!>        The second MQC NEF Matrix variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      function mqc_nefmatrix_difference(NEFMatA,NEFMatB) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA,NEFMatB
      type(mqc_nef_matrix)::NEFMatOut
      integer::col,row,id
!
      if(NEFMatA%NumSpace_col.ne.NEFMatB%NumSpace_col.or.&
         NEFMatA%NumSpace_row.ne.NEFMatB%NumSpace_row) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_difference',6,&
        'NEFMatA%NumSpace_col',NEFMatA%NumSpace_col,&
        'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col,&
        'NEFMatA%NumSpace_row',NEFMatA%NumSpace_row,&
        'NEFMatB%NumSpace_row',NEFMatB%NumSpace_row)

      do col = 1, NEFMatA%NumSpace_col
        if(NEFMatA%Size_vector_alpha_col(col)&
           .ne.NEFMatB%Size_vector_alpha_col(col).or. &
               NEFMatA%Size_vector_beta_col(col)  &
           .ne.NEFMatB%Size_vector_beta_col(col)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_difference',6,&
          'Sub-space index (col)', col,&
          'NEFMatA%Size_vector_alpha_col(col)',NEFMatA%Size_vector_alpha_col(col),&
          'NEFMatB%Size_vector_alpha_col(col)',NEFMatB%Size_vector_alpha_col(col),&
          'NEFMatA%Size_vector_beta_col(col)',NEFMatA%Size_vector_beta_col(col),&
          'NEFMatB%Size_vector_beta_col(col)',NEFMatB%Size_vector_beta_col(col))
      end do
      do row = 1, NEFMatA%NumSpace_row
        if(NEFMatA%Size_vector_alpha_row(row)&
           .ne.NEFMatB%Size_vector_alpha_row(row).or. &
               NEFMatA%Size_vector_beta_row(row)  &
           .ne.NEFMatB%Size_vector_beta_row(row)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_difference',6,&
          'Sub-space index (row)', row,&
          'NEFMatA%Size_vector_alpha_row(row)',NEFMatA%Size_vector_alpha_row(row),&
          'NEFMatB%Size_vector_alpha_row(row)',NEFMatB%Size_vector_alpha_row(row),&
          'NEFMatA%Size_vector_beta_row(row)',NEFMatA%Size_vector_beta_row(row),&
          'NEFMatB%Size_vector_beta_row(row)',NEFMatB%Size_vector_beta_row(row))
      end do

      call NEFMatOut%init(NEFMatA%Map_Matrix_Row,&
          initial_map_matrix_col=NEFMatB%Map_Matrix_Col)
      do row = 1, NEFMatA%NumSpace_row
        do col = 1, NEFMatB%NumSpace_col
          id = NEFMatOut%blockIndex([1,row],[2,col])
          if (NEFMatA%Zero_Flag(id).and.NEFMatB%Zero_Flag(id)) then
            continue
          else if (NEFMatA%Zero_Flag(id)) then
            NEFMatOut%Main_Matrix(id) = (-1.0)*NEFMatB%Main_Matrix(id)
            NEFMatOut%Zero_Flag(id) = .false.
          else if (NEFMatB%Zero_Flag(id)) then
            NEFMatOut%Main_Matrix(id) = NEFMatA%Main_Matrix(id)
            NEFMatOut%Zero_Flag(id) = .false.
          else
            NEFMatOut%Main_Matrix(id) = NEFMatA%Main_Matrix(id)-NEFMatB%Main_Matrix(id)
            NEFMatOut%Zero_Flag(id) = .false.
          end if
        end do
      end do
      end function mqc_nefmatrix_difference
!
!
!     PROCEDURE MQC_NEF_R4_Difference
!
!>    \brief <b> MQC_NEF_R4_Difference is used to subtract two MQC NEF R4 type variables</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Difference is used to subtract two MQC NEF R4 type variables
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4)
!>        The first MQC NEF R4 variable to sum.
!>    \endverbatim
!>    
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4)
!>        The second MQC NEF R4 variable to sum.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      function mqc_nef_r4_difference(NEFR4A,NEFR4B) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      type(mqc_nef_r4)::NEFR4Out
      integer::id1,id2,id3,id4,id
!
      if(NEFR4A%NumSpace_1.ne.NEFR4B%NumSpace_1)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_difference',6,&
        'NEFR4A%NumSpace_1',NEFR4A%NumSpace_1,&
        'NEFR4B%NumSpace_1',NEFR4B%NumSpace_1)
      if(NEFR4A%NumSpace_2.ne.NEFR4B%NumSpace_2)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_difference',6,&
        'NEFR4A%NumSpace_2',NEFR4A%NumSpace_2,&
        'NEFR4B%NumSpace_2',NEFR4B%NumSpace_2)
      if(NEFR4A%NumSpace_3.ne.NEFR4B%NumSpace_3)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_difference',6,&
        'NEFR4A%NumSpace_3',NEFR4A%NumSpace_3,&
        'NEFR4B%NumSpace_3',NEFR4B%NumSpace_3)
      if(NEFR4A%NumSpace_4.ne.NEFR4B%NumSpace_4)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_difference',6,&
        'NEFR4A%NumSpace_4',NEFR4A%NumSpace_4,&
        'NEFR4B%NumSpace_4',NEFR4B%NumSpace_4)
     
      do id1 = 1, NEFR4A%NumSpace_1
        if(NEFR4A%Size_vector_alpha_1(id1)&
           .ne.NEFR4B%Size_vector_alpha_1(id1).or.&
           NEFR4A%Size_vector_beta_1(id1)&
           .ne.NEFR4B%Size_vector_beta_1(id1))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_difference',6,&
          'Sub-space index (id1)', id1,&
          'NEFR4A%Size_vector_alpha_1(id1)',NEFR4A%Size_vector_alpha_1(id1),&
          'NEFR4B%Size_vector_alpha_1(id1)',NEFR4B%Size_vector_alpha_1(id1),&
          'NEFR4A%Size_vector_beta_1(id1)',NEFR4A%Size_vector_beta_1(id1),&
          'NEFR4B%Size_vector_beta_1(id1)',NEFR4B%Size_vector_beta_1(id1))
      end do
      do id2 = 1, NEFR4A%NumSpace_2
        if(NEFR4A%Size_vector_alpha_2(id2)&
           .ne.NEFR4B%Size_vector_alpha_2(id2).or.&
           NEFR4A%Size_vector_beta_2(id2)&
           .ne.NEFR4B%Size_vector_beta_2(id2))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_difference',6,&
          'Sub-space index (id2)', id2,&
          'NEFR4A%Size_vector_alpha_2(id2)',NEFR4A%Size_vector_alpha_2(id2),&
          'NEFR4B%Size_vector_alpha_2(id2)',NEFR4B%Size_vector_alpha_2(id2),&
           'NEFR4A%Size_vector_beta_2(id2)', NEFR4A%Size_vector_beta_2(id2),&
           'NEFR4B%Size_vector_beta_2(id2)', NEFR4B%Size_vector_beta_2(id2))
      end do
      do id3 = 1, NEFR4A%NumSpace_3
        if(NEFR4A%Size_vector_alpha_3(id3)&
           .ne.NEFR4B%Size_vector_alpha_3(id3).or.&
           NEFR4A%Size_vector_beta_3(id3)&
           .ne.NEFR4B%Size_vector_beta_3(id3))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_difference',6,&
          'Sub-space index (id3)', id3,&
          'NEFR4A%Size_vector_alpha_3(id3)',NEFR4A%Size_vector_alpha_3(id3),&
          'NEFR4B%Size_vector_alpha_3(id3)',NEFR4B%Size_vector_alpha_3(id3),&
           'NEFR4A%Size_vector_beta_3(id3)', NEFR4A%Size_vector_beta_3(id3),&
           'NEFR4B%Size_vector_beta_3(id3)', NEFR4B%Size_vector_beta_3(id3))
      end do
      do id4 = 1, NEFR4A%NumSpace_4
        if(NEFR4A%Size_vector_alpha_4(id4)&
           .ne.NEFR4B%Size_vector_alpha_4(id4).or.&
           NEFR4A%Size_vector_beta_4(id4)&
           .ne.NEFR4B%Size_vector_beta_4(id4))&
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nef_r4_difference',6,&
          'Sub-space index (id4)', id4,&
          'NEFR4A%Size_vector_alpha_4(id4)',NEFR4A%Size_vector_alpha_4(id4),&
          'NEFR4B%Size_vector_alpha_4(id4)',NEFR4B%Size_vector_alpha_4(id4),&
           'NEFR4A%Size_vector_beta_4(id4)', NEFR4A%Size_vector_beta_4(id4),&
           'NEFR4B%Size_vector_beta_4(id4)', NEFR4B%Size_vector_beta_4(id4))
      end do


      NEFR4Out = NEFR4A
      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              id = NEFR4Out%blockIndex([4,id4],[3,id3],[2,id2],[1,id1])
              if (NEFR4A%Zero_Flag(id).and.NEFR4B%Zero_Flag(id)) then
                continue
              else if (NEFR4A%Zero_Flag(id)) then
                NEFR4Out%Main_R4(id) = (-1.0)*NEFR4B%Main_R4(id)
                NEFR4Out%Zero_Flag(id) = .false.
              else if (NEFR4B%Zero_Flag(id)) then
                NEFR4Out%Main_R4(id) = NEFR4A%Main_R4(id)
                NEFR4Out%Zero_Flag(id) = .false.
              else
                NEFR4Out%Main_R4(id) = NEFR4A%Main_R4(id) - NEFR4B%Main_R4(id)
                NEFR4Out%Zero_Flag(id) = .false.
              end if
           end do
          end do
        end do
      end do
      end function mqc_nef_r4_difference
!
!
!     PROCEDURE MQC_Print_NEF_R4
!
!>    \brief <b> MQC_Print_NEF_R4 is a subroutine used to print an MQC
!>    NEF R4 Tensor object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_Print_NEF_R4 is a subroutine used to print an MQC NEF R4 Tensor object. 
!>    Blank_At_Top and Blank_At_Bottom are optional logical arguments to print blank
!>    lines before or after output.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        Integral is Class(MQC_NEF_R4)
!>        The MQC NEF R4 tensor variable to be printed.
!>    \endverbatim
!>
!>    \param[in] IOut
!>    \verbatim
!>        IOut is Integer(kind=int64) 
!>        The Fortran file number to print to. 
!>    \endverbatim
!>
!>    \param[in] Header
!>    \verbatim
!>        Header is Character(Len=*)
!>        The title to print along with Integral.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Top
!>    \verbatim
!>        Blank_At_Top is Logical,Optional
!>        = .True.:  print blank line above output
!>        = .False.: do not print blank line above output.
!>    \endverbatim
!>
!>    \param[in] Blank_At_Bottom
!>    \verbatim
!>        Blank_At_Bottom is Logical,Optional
!>        = .True.:  print blank line below output
!>        = .False.: do not print blank line below output.
!>    \endverbatim
!>
!>    \param[in] FormatStr
!>    \verbatim
!>        FormatStr is Character(Len=*),Optional
!>        Formatting statement for printing.
!>    \endverbatim
!>
!>    \param[in] Short_Print_Threshold
!>    \verbatim
!>        Short_Print_Threshold is Real,Optional
!>        if present: print only elements with absolute value larger than Short_Print_Threshold
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M.M.F. Moraes
!>    \date 2026
!
      subroutine mqc_print_nef_r4(NEFR4,iOut,header, &
          blank_at_top,blank_at_bottom,formatStr,short_print_threshold)
!
      implicit none
      class(mqc_nef_r4)::NEFR4
      integer(kind=int64),intent(in)::iOut
      character(len=*),intent(in)::header
      character(len=64)::arrayType
      logical,intent(In),optional::blank_at_top,blank_at_bottom
      character(len=*),intent(in),optional::formatStr
      Real(kind=real64),intent(in),optional::short_print_threshold
      Real(kind=real64)::print_threshold
      integer::id,id1,id2,id3,id4
!
 1000 Format(1x,A)
 1020 Format( " " )
!  
!  
      if(present(blank_at_top)) then
        if(blank_at_top) write(iout,1020)
      endif
      If(Present(short_print_threshold)) then
        print_threshold = short_print_threshold
      else
        print_threshold = -1.0
      EndIf

      write(iout,1000) trim(header)
      do id1 = 1, NEFR4%NumSpace_1
        do id2 = 1, NEFR4%NumSpace_2
          do id3 = 1, NEFR4%NumSpace_3
            do id4 = 1, NEFR4%NumSpace_4
              id = NEFR4%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              if (NEFR4%Zero_Flag(id)) then
                write(iOut,1000) &
   'Sub-Space ('//trim(num2char(id1))//trim(num2char(id2))//'|'//trim(num2char(id3))//trim(num2char(id4))//') '&
               &' (all zero)'
              else
              call mqc_print(NEFR4%Main_R4(id),iOut,&
   'Sub-Space ('//trim(num2char(id1))//trim(num2char(id2))//'|'//trim(num2char(id3))//trim(num2char(id4))//') '&
               &' '//trim(NEFR4%Main_R4(id)%storType())//' '//trim(NEFR4%Main_R4(id)%type()),&
              short_print_threshold=print_threshold)
              end if
            end do
          end do
        end do
      end do
        
      end subroutine mqc_print_nef_r4
!
!     
!     
!     PROCEDURE MQC_NEFMatrix_NEFMatrix_Outer_Product
!     
!>    \brief <b> MQC_NEFMatrix_NEFMatrix_Outer_Product is used to calculate the
!>    outer product between two MQC NEF Matrix type variabel and returns an MQC NEF R4</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_NEFMatrix_Outer_Product is used to calculate the outer product between
!>    two MQC NEF Matrix and returns an MQC NEF R4. 
!>
!>    \endverbatim
!     
!     Arguments: 
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix)
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix)
!>    \endverbatim
!
!         
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!     
      function mqc_nefmatrix_nefmatrix_outer_product(NEFMatA,NEFMatB) result(NEFR4Out)
!
      implicit none
      integer(kind=int64) :: nef_bidx
      class(mqc_nef_matrix),intent(in)::NEFMatA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_r4)::NEFR4Out

      integer::idA,idB,id,id1,id2,id3,id4,N1,N12,N123,N1234

      NEFR4Out%NumSpace_2 = NEFMatA%NumSpace_col
      NEFR4Out%NumSpace_1 = NEFMatA%NumSpace_row
      NEFR4Out%NumSpace_4 = NEFMatB%NumSpace_col
      NEFR4Out%NumSpace_3 = NEFMatB%NumSpace_row
      N1 = NEFR4Out%NumSpace_4
      N12 = N1*NEFR4Out%NumSpace_3
      N123 = N12*NEFR4Out%NumSpace_2
      N1234 = N123*NEFR4Out%NumSpace_1

      allocate(NEFR4Out%Main_R4(N1234))
      allocate(NEFR4Out%Zero_Flag(N1234))
      NEFR4Out%Zero_Flag = .true.
      NEFR4Out%Map_Matrix_2 = NEFMatA%Map_Matrix_Col
      NEFR4Out%Map_Matrix_1 = NEFMatA%Map_Matrix_Row
      NEFR4Out%Map_Matrix_4 = NEFMatB%Map_Matrix_Col
      NEFR4Out%Map_Matrix_3 = NEFMatB%Map_Matrix_Row
      NEFR4Out%Size_vector_alpha_2 = NEFMatA%Size_vector_alpha_Col
      NEFR4Out%Size_vector_alpha_1 = NEFMatA%Size_vector_alpha_Row
      NEFR4Out%Size_vector_alpha_4 = NEFMatB%Size_vector_alpha_Col
      NEFR4Out%Size_vector_alpha_3 = NEFMatB%Size_vector_alpha_Row
      NEFR4Out%Size_vector_beta_2 = NEFMatA%Size_vector_beta_Col
      NEFR4Out%Size_vector_beta_1 = NEFMatA%Size_vector_beta_Row
      NEFR4Out%Size_vector_beta_4 = NEFMatB%Size_vector_beta_Col
      NEFR4Out%Size_vector_beta_3 = NEFMatB%Size_vector_beta_Row

      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          idA = NEFMatA%blockIndex([1,id1],[2,id2])
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              idB = NEFMatB%blockIndex([1,id3],[2,id4])
              nef_bidx = NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              if (.not.(NEFMatA%Zero_Flag(idA).or.NEFMatB%Zero_Flag(idB))) then
                NEFR4Out%Main_R4(nef_bidx) = &
                  NEFMatA%Main_Matrix(idA)%outer2ERI(NEFMatB%Main_Matrix(idB))
                NEFR4Out%Zero_Flag(nef_bidx) = .false.
              end if
            end do
          end do
        end do
      end do
!
      end function mqc_nefmatrix_nefmatrix_outer_product
!
!
!     PROCEDURE MQC_NEFVector_NEFVector_Outer_Product
!
!>    \brief <b> MQC_NEFVector_NEFVector_Outer_Product is used to calculate the
!>    outer product between two MQC NEF Vector type variables and returns an
!>    MQC NEF Matrix</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFVector_NEFVector_Outer_Product is used to calculate the outer
!>    product between two MQC NEF Vector and returns an MQC NEF Matrix,
!>    mirroring MQC_NEFMatrix_NEFMatrix_Outer_Product one rank lower
!>    (NEFVecA's own space axis becomes the output's row axis, NEFVecB's
!>    becomes the column axis).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVecA
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>    \endverbatim
!>
!>    \param[in] NEFVecB
!>    \verbatim
!>        NEFVecB is type(mqc_nef_vector)
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefvector_nefvector_outer_product(NEFVecA,NEFVecB) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_vector),intent(in)::NEFVecA
      type(mqc_nef_vector),intent(in)::NEFVecB
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::id,id1,id2
!
      NEFMatOut%NumSpace_row = NEFVecA%NumSpace
      NEFMatOut%NumSpace_col = NEFVecB%NumSpace
      allocate(NEFMatOut%Main_Matrix(NEFMatOut%NumSpace_row*NEFMatOut%NumSpace_col))
      allocate(NEFMatOut%Zero_Flag(size(NEFMatOut%Main_Matrix)))
      NEFMatOut%Zero_Flag = .true.
      NEFMatOut%Map_Matrix_Row = NEFVecA%Map_Matrix
      NEFMatOut%Map_Matrix_Col = NEFVecB%Map_Matrix
      NEFMatOut%Size_vector_alpha_row = NEFVecA%Size_vector_alpha
      NEFMatOut%Size_vector_beta_row = NEFVecA%Size_vector_beta
      NEFMatOut%Size_vector_alpha_col = NEFVecB%Size_vector_alpha
      NEFMatOut%Size_vector_beta_col = NEFVecB%Size_vector_beta
!
      do id1 = 1, NEFMatOut%NumSpace_row
        do id2 = 1, NEFMatOut%NumSpace_col
          id = NEFMatOut%blockIndex([1,id1],[2,id2])
          if (.not.(NEFVecA%Zero_Flag(id1).or.NEFVecB%Zero_Flag(id2))) then
            NEFMatOut%Main_Matrix(id) = NEFVecA%Main_Vector(id1)%outer(NEFVecB%Main_Vector(id2))
            NEFMatOut%Zero_Flag(id) = .false.
          end if
        end do
      end do
!
      end function mqc_nefvector_nefvector_outer_product
!
!
!
!     PROCEDURE MQC_NEFMatrix_NEFMatrix_Crossouter_Product
!     
!>    \brief <b> MQC_NEFMatrix_NEFMatrix_crossouter_Product is used to calculate the
!>    crossouter product between two MQC NEF Matrix type variabel and returns an MQC NEF R4</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_NEFMatrix_Crossouter_Product is used to calculate the crossouter product between
!>    two MQC NEF Matrix and returns an MQC NEF R4. 
!>
!>    \endverbatim
!     
!     Arguments: 
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix)
!>    \endverbatim
!>    
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix)
!>    \endverbatim
!
!         
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!     
      function mqc_nefmatrix_nefmatrix_crossouter_product(NEFMatA,NEFMatB) result(NEFR4Out)
!
      implicit none
      integer(kind=int64) :: nef_bidx
      class(mqc_nef_matrix),intent(in)::NEFMatA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_r4)::NEFR4Out

      integer::idA,idB,id,id1,id2,id3,id4,N1,N12,N123,N1234

      NEFR4Out%NumSpace_1 = NEFMatA%NumSpace_row
      NEFR4Out%NumSpace_4 = NEFMatA%NumSpace_col
      NEFR4Out%NumSpace_3 = NEFMatB%NumSpace_row
      NEFR4Out%NumSpace_2 = NEFMatB%NumSpace_col
      N1 = NEFR4Out%NumSpace_4
      N12 = N1*NEFR4Out%NumSpace_3
      N123 = N12*NEFR4Out%NumSpace_2
      N1234 = N123*NEFR4Out%NumSpace_1

      allocate(NEFR4Out%Main_R4(N1234))
      allocate(NEFR4Out%Zero_Flag(N1234))
      NEFR4Out%Zero_Flag = .true.
      NEFR4Out%Map_Matrix_1 = NEFMatA%Map_Matrix_Row
      NEFR4Out%Map_Matrix_4 = NEFMatA%Map_Matrix_Col
      NEFR4Out%Map_Matrix_3 = NEFMatB%Map_Matrix_Row
      NEFR4Out%Map_Matrix_2 = NEFMatB%Map_Matrix_Col
      NEFR4Out%Size_vector_alpha_1 = NEFMatA%Size_vector_alpha_Row
      NEFR4Out%Size_vector_alpha_4 = NEFMatA%Size_vector_alpha_Col
      NEFR4Out%Size_vector_alpha_3 = NEFMatB%Size_vector_alpha_Row
      NEFR4Out%Size_vector_alpha_2 = NEFMatB%Size_vector_alpha_Col
      NEFR4Out%Size_vector_beta_1 = NEFMatA%Size_vector_beta_Row
      NEFR4Out%Size_vector_beta_4 = NEFMatA%Size_vector_beta_Col
      NEFR4Out%Size_vector_beta_3 = NEFMatB%Size_vector_beta_Row
      NEFR4Out%Size_vector_beta_2 = NEFMatB%Size_vector_beta_Col

      do id1 = 1, NEFR4Out%NumSpace_1
        do id4 = 1, NEFR4Out%NumSpace_4
          idA = NEFMatA%blockIndex([1,id1],[2,id4])
          do id3 = 1, NEFR4Out%NumSpace_3
            do id2 = 1, NEFR4Out%NumSpace_2
              idB = NEFMatB%blockIndex([1,id3],[2,id2])
              nef_bidx = NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              if (.not.(NEFMatA%Zero_Flag(idA).or.NEFMatB%Zero_Flag(idB))) then
                NEFR4Out%Main_R4(nef_bidx) = &
                  NEFMatA%Main_Matrix(idA)%crossouter2ERI(NEFMatB%Main_Matrix(idB))
                NEFR4Out%Zero_Flag(nef_bidx) = .false.
              end if
            end do
          end do
        end do
      end do
!
      end function mqc_nefmatrix_nefmatrix_crossouter_product
!
!
!     PROCEDURE MQC_NEFVector_NEFVector_Contraction
!
!
!>    \brief <b> MQC_NEFVector_NEFVector_Contraction is used to return the
!>    contraction (dot product) of two MQC NEF Vector objects.</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFVector_NEFVector_Contraction is used to return the
!>    contraction of two MQC NEF Vector objects. A vector has only one
!>    axis, so "contract" and "multiply" coincide here -- this sums the
!>    bare-vector dot product (Transpose(block)).dot.block, via the
!>    plain mqc_vector Operator(.dot.), over every block.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVecA
!>    \verbatim
!>        NEFVecA is type(mqc_nef_vector)
!>        The first MQC NEF Vector to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFVecB
!>    \verbatim
!>        NEFVecB is type(mqc_nef_vector)
!>        The second MQC NEF Vector to be contracted.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefvector_nefvector_contraction(NEFVecA,NEFVecB) result(output)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecA,NEFVecB
      type(mqc_scalar)::output
      integer::id
!
      if(NEFVecA%NumSpace.ne.NEFVecB%NumSpace) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefvector_nefvector_contraction',6,&
        'NEFVecA%NumSpace',NEFVecA%NumSpace,&
        'NEFVecB%NumSpace',NEFVecB%NumSpace)
!
      do id = 1, NEFVecA%NumSpace
        if(NEFVecA%Size_vector_alpha(id).ne.NEFVecB%Size_vector_alpha(id).or. &
           NEFVecA%Size_vector_beta(id).ne.NEFVecB%Size_vector_beta(id)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefvector_nefvector_contraction',6,&
          'Sub-space index', id,&
          'NEFVecA%Size_vector_alpha(id)',NEFVecA%Size_vector_alpha(id),&
          'NEFVecB%Size_vector_alpha(id)',NEFVecB%Size_vector_alpha(id),&
          'NEFVecA%Size_vector_beta(id)',NEFVecA%Size_vector_beta(id),&
          'NEFVecB%Size_vector_beta(id)',NEFVecB%Size_vector_beta(id))
      end do
!
      output = 0.0
      do id = 1, NEFVecA%NumSpace
        if(NEFVecA%Zero_Flag(id).or.NEFVecB%Zero_Flag(id)) cycle
        output = output + (Transpose(NEFVecA%Main_Vector(id)%getBlock())&
                           .dot.NEFVecB%Main_Vector(id)%getBlock())
      end do
!
      end function mqc_nefvector_nefvector_contraction
!
!
!     PROCEDURE MQC_NEFMatrix_NEFMatrix_Contraction
!
!
!>    \brief <b> MQC_NEFMatrix_NEFMatrix_Contraction is used to return the contraction of two 
!>    MQC NEF Matrix objects.</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim 
!>
!>    MQC_NEFMatrix_NEFMatrix_Contraction is used to return the contraction of two 
!>    MQC NEF Matrix objects.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is type(mqc_nef_matrix) 
!>        The first MQC NEF Matrix to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix) 
!>        The second MQC NEF Matrix to be contracted.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_nefmatrix_contraction(NEFMatA,NEFMatB) result(output)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA,NEFMatB
      type(mqc_scalar)::output
      integer::row,col,id
!
      if(NEFMatA%NumSpace_col.ne.NEFMatB%NumSpace_col.or.&
         NEFMatA%NumSpace_row.ne.NEFMatB%NumSpace_row) &
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nefmatrix_sum',6,&
        'NEFMatA%NumSpace_col',NEFMatA%NumSpace_col,&
        'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col,&
        'NEFMatA%NumSpace_row',NEFMatA%NumSpace_row,&
        'NEFMatB%NumSpace_row',NEFMatB%NumSpace_row)

      do col = 1, NEFMatA%NumSpace_col
        if(NEFMatA%Size_vector_alpha_col(col)&
           .ne.NEFMatB%Size_vector_alpha_col(col).or. &
               NEFMatA%Size_vector_beta_col(col)  &
           .ne.NEFMatB%Size_vector_beta_col(col)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_sum',6,&
          'Sub-space index (col)', col,&
          'NEFMatA%Size_vector_alpha_col(col)',NEFMatA%Size_vector_alpha_col(col),&
          'NEFMatB%Size_vector_alpha_col(col)',NEFMatB%Size_vector_alpha_col(col),&
          'NEFMatA%Size_vector_beta_col(col)',NEFMatA%Size_vector_beta_col(col),&
          'NEFMatB%Size_vector_beta_col(col)',NEFMatB%Size_vector_beta_col(col))
      end do
      do row = 1, NEFMatA%NumSpace_row
        if(NEFMatA%Size_vector_alpha_row(row)&
           .ne.NEFMatB%Size_vector_alpha_row(row).or. &
               NEFMatA%Size_vector_beta_row(row)  &
           .ne.NEFMatB%Size_vector_beta_row(row)) &
          call mqc_error_i('NEF subspace dimensions are &
           &not conformable in mqc_nefmatrix_sum',6,&
          'Sub-space index (row)', row,&
          'NEFMatA%Size_vector_alpha_row(row)',NEFMatA%Size_vector_alpha_row(row),&
          'NEFMatB%Size_vector_alpha_row(row)',NEFMatB%Size_vector_alpha_row(row),&
          'NEFMatA%Size_vector_beta_row(row)',NEFMatA%Size_vector_beta_row(row),&
          'NEFMatB%Size_vector_beta_row(row)',NEFMatB%Size_vector_beta_row(row))
      end do
      output = 0.0
      do row = 1, NEFMatA%NumSpace_row
        do col = 1, NEFMatA%NumSpace_col
          id = NEFMatA%blockIndex([1,row],[2,col])
          if (NEFMatA%Zero_Flag(id).or.NEFMatB%Zero_Flag(id)) cycle
          output = output + Contraction(NEFMatA%Main_Matrix(id),NEFMatB%Main_Matrix(id))
        end do
      end do
!
      end function mqc_nefmatrix_nefmatrix_contraction
!
!
!     PROCEDURE MQC_NEF_R4_NEF_R4_Contraction
!
!
!>    \brief <b> MQC_NEF_R4_NEF_R4_Contraction is used to return the contraction of two 
!>    MQC NEF R4 objects.</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_NEF_R4_Contraction is used to return the contraction of two 
!>    MQC NEF R4 objects.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4) 
!>        The second MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is character(len=*),optional 
!>        = 'coulomb':   contracts with coulomb integrals.
!>        = 'exchange':  contracts with exchange integrals.
!>        = 'doublebar': contracts with double bar integrals
!>                       (default).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_nef_r4_contraction(NEFR4A,NEFR4B,label) result(output)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      character(len=*),optional,intent(in)::label
      character(len=64)::myLabel
      type(mqc_scalar)::output
      integer::id1,id2,id3,id4,id,idex
      type(mqc_twoeris):: debbug_ERI
!
      if(present(label)) then
        call string_change_case(label,'L',myLabel)
      else
        myLabel = 'doublebar'
      endIf
!
      if(myLabel.ne.'coulomb'.and.myLabel.ne.'exchange'.and.&
         myLabel.ne.'doublebar') &
        call mqc_error_a('Unrecognized 2ERI integral interaction specified&
                        &in mqc_nef_r4_nef_r4_contraction.',&
                          6,'myLabel',myLabel)
!
      if(NEFR4A%NumSpace_1.ne.NEFR4B%NumSpace_1)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
        'NEFR4A%NumSpace_1',NEFR4A%NumSpace_1,&
        'NEFR4B%NumSpace_1',NEFR4B%NumSpace_1)
      if(NEFR4A%NumSpace_2.ne.NEFR4B%NumSpace_2)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
        'NEFR4A%NumSpace_2',NEFR4A%NumSpace_2,&
        'NEFR4B%NumSpace_2',NEFR4B%NumSpace_2)
      if(NEFR4A%NumSpace_3.ne.NEFR4B%NumSpace_3)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
        'NEFR4A%NumSpace_3',NEFR4A%NumSpace_3,&
        'NEFR4B%NumSpace_3',NEFR4B%NumSpace_3)
      if(NEFR4A%NumSpace_4.ne.NEFR4B%NumSpace_4)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
        'NEFR4A%NumSpace_4',NEFR4A%NumSpace_4,&
        'NEFR4B%NumSpace_4',NEFR4B%NumSpace_4)
     
      do id1 = 1, NEFR4A%NumSpace_1
        if(NEFR4A%Size_vector_alpha_1(id1)&
           .ne.NEFR4B%Size_vector_alpha_1(id1).or.&
           NEFR4A%Size_vector_beta_1(id1)&
           .ne.NEFR4B%Size_vector_beta_1(id1))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
          'Sub-space index (id1)', id1,&
          'NEFR4A%Size_vector_alpha_1(id1)',NEFR4A%Size_vector_alpha_1(id1),&
          'NEFR4B%Size_vector_alpha_1(id1)',NEFR4B%Size_vector_alpha_1(id1),&
          'NEFR4A%Size_vector_beta_1(id1)',NEFR4A%Size_vector_beta_1(id1),&
          'NEFR4B%Size_vector_beta_1(id1)',NEFR4B%Size_vector_beta_1(id1))
      end do
      do id2 = 1, NEFR4A%NumSpace_2
        if(NEFR4A%Size_vector_alpha_2(id2)&
           .ne.NEFR4B%Size_vector_alpha_2(id2).or.&
           NEFR4A%Size_vector_beta_2(id2)&
           .ne.NEFR4B%Size_vector_beta_2(id2))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
          'Sub-space index (id2)', id2,&
          'NEFR4A%Size_vector_alpha_2(id2)',NEFR4A%Size_vector_alpha_2(id2),&
          'NEFR4B%Size_vector_alpha_2(id2)',NEFR4B%Size_vector_alpha_2(id2),&
           'NEFR4A%Size_vector_beta_2(id2)', NEFR4A%Size_vector_beta_2(id2),&
           'NEFR4B%Size_vector_beta_2(id2)', NEFR4B%Size_vector_beta_2(id2))
      end do
      do id3 = 1, NEFR4A%NumSpace_3
        if(NEFR4A%Size_vector_alpha_3(id3)&
           .ne.NEFR4B%Size_vector_alpha_3(id3).or.&
           NEFR4A%Size_vector_beta_3(id3)&
           .ne.NEFR4B%Size_vector_beta_3(id3))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
          'Sub-space index (id3)', id3,&
          'NEFR4A%Size_vector_alpha_3(id3)',NEFR4A%Size_vector_alpha_3(id3),&
          'NEFR4B%Size_vector_alpha_3(id3)',NEFR4B%Size_vector_alpha_3(id3),&
           'NEFR4A%Size_vector_beta_3(id3)', NEFR4A%Size_vector_beta_3(id3),&
           'NEFR4B%Size_vector_beta_3(id3)', NEFR4B%Size_vector_beta_3(id3))
      end do
      do id4 = 1, NEFR4A%NumSpace_4
        if(NEFR4A%Size_vector_alpha_4(id4)&
           .ne.NEFR4B%Size_vector_alpha_4(id4).or.&
           NEFR4A%Size_vector_beta_4(id4)&
           .ne.NEFR4B%Size_vector_beta_4(id4))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_r4_contraction',6,&
          'Sub-space index (id4)', id4,&
          'NEFR4A%Size_vector_alpha_4(id4)',NEFR4A%Size_vector_alpha_4(id4),&
          'NEFR4B%Size_vector_alpha_4(id4)',NEFR4B%Size_vector_alpha_4(id4),&
           'NEFR4A%Size_vector_beta_4(id4)', NEFR4A%Size_vector_beta_4(id4),&
           'NEFR4B%Size_vector_beta_4(id4)', NEFR4B%Size_vector_beta_4(id4))
      end do
!
      output = 0.0
      do id1 = 1, NEFR4A%NumSpace_1
        do id2 = 1, NEFR4A%NumSpace_2
          do id3 = 1, NEFR4A%NumSpace_3
            do id4 = 1, NEFR4A%NumSpace_4
              id = NEFR4A%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
                if (.not.(NEFR4A%Zero_Flag(id).or.NEFR4B%Zero_Flag(id))) &
                  output = output + contraction(NEFR4A%Main_R4(id),NEFR4B%Main_R4(id),'coulomb')
              end if
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
                 idex = NEFR4A%blockIndex([1,id1],[4,id2],[3,id3],[2,id4])
                 if (.not.(NEFR4A%Zero_Flag(idex).or.NEFR4B%Zero_Flag(id))) &
                   output = output - &
                     contraction(NEFR4A%Main_R4(idex),NEFR4B%Main_R4(id)%swap(2,4),'coulomb')
              end if
            end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_nef_r4_contraction
!
!
!     PROCEDURE MQC_NEF_R4_NEF_Matrix_Contraction
!
!
!>    \brief <b> MQC_NEF_R4_NEF_Matrix_Contraction is used to return the contraction of a
!>    MQC NEF R4 and MQC NEF Matrix objects.</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_NEF_Matrix_Contraction is used to return the contraction of a
!>    MQC NEF R4 and MQC NEF Matrix objects.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is type(mqc_nef_matrix) 
!>        The MQC NEF Matrix to be contracted.
!>    \endverbatim
!>
!>    \param[in] Label
!>    \verbatim
!>        Label is character(len=*),optional 
!>        = 'coulomb':   contracts with coulomb integrals.
!>        = 'exchange':  contracts with exchange integrals.
!>        = 'doublebar': contracts with double bar integrals
!>                       (default).
!>    \endverbatim
!>
!>    \param[in] IndexA
!>    \verbatim
!>        IndexA is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIsA to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_nef_matrix_contraction(NEFR4A,NEFMatB,label,indexA) result(output)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A
      type(mqc_nef_matrix),intent(in)::NEFMatB
      character(len=*),optional,intent(in)::label
      integer(kind=int64),dimension(4),intent(in),optional::indexA
      integer(kind=int64),dimension(4)::my_indexA
      character(len=64)::myLabel
      type(mqc_nef_matrix)::output
      integer::id1,id2,id3,id4,idA,idB,idex,idOut
      type(mqc_twoeris):: debbug_ERI
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la
!
      if(present(label)) then
        call string_change_case(label,'L',myLabel)
      else
        myLabel = 'doublebar'
      endIf
!
      if(myLabel.ne.'coulomb'.and.myLabel.ne.'exchange'.and.&
         myLabel.ne.'doublebar') &
        call mqc_error_a('Unrecognized 2ERI integral interaction specified&
                        &in mqc_nef_r4_nef_r4_contraction.',&
                          6,'myLabel',myLabel)
!
      if(present(indexA)) then
        ia=indexA(1)
        ja=indexA(2)
        ka=indexA(3)
        la=indexA(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_Matrix_Contraction',&
                                        6,'ka', ka,'la',la)
        
        my_indexA = indexA
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_indexA = [1,2,3,4]
      endif
!
      if(NEFR4A%blockSize(ia).ne.NEFMatB%NumSpace_row)&
        call mqc_error_i('NEF number of subspaces are &
        &not conformable in mqc_nef_r4_nef_Matrix_contraction',6,&
        'NEFR4A%blockSize(ia)',NEFR4A%blockSize(ia),&
        'NEFMatB%NumSpace_row',NEFMatB%NumSpace_row)
      if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
        if(NEFR4A%blockSize(ja).ne.NEFMatB%NumSpace_col)&
          call mqc_error_i('NEF number of subspaces are &
          &not conformable in mqc_nef_r4_nef_Matrix_contraction',6,&
          'NEFR4A%blockSize(ja)',NEFR4A%blockSize(ja),&
          'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col)
      end if
      if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
        if(NEFR4A%blockSize(la).ne.NEFMatB%NumSpace_col)&
          call mqc_error_i('NEF number of subspaces are &
          &not conformable for exchange in mqc_nef_r4_nef_Matrix_contraction',6,&
          'NEFR4A%blockSize(la)',NEFR4A%blockSize(la),&
          'NEFMatB%NumSpace_col',NEFMatB%NumSpace_col)
      end if
     
      do id1 = 1, NEFR4A%blockSize(ia)
        if(NEFR4A%blockSize(ia,id1,'alpha')&
           .ne.NEFMatB%Size_vector_alpha_row(id1).or.&
           NEFR4A%blockSize(ia,id1,'beta')&
           .ne.NEFMatB%Size_vector_beta_row(id1))&
          call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nef_r4_nef_Matrix_contraction',6,&
          'Sub-space index (id1)', id1,&
          'NEFR4A%blockSize(ia,id1,alpha)',NEFR4A%blockSize(ia,id1,'alpha'),&
          'NEFMatB%Size_vector_alpha_row(id1)',NEFMatB%Size_vector_alpha_row(id1),&
          'NEFR4A%blockSize(ia,id1,beta)',NEFR4A%blockSize(ia,id1,'beta'),&
          'NEFMatB%Size_vector_beta_row(id1)',NEFMatB%Size_vector_beta_row(id1))
      end do

      if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
        do id1 = 1, NEFR4A%blockSize(ja)
          if(NEFR4A%blockSize(ja,id1,'alpha')&
             .ne.NEFMatB%Size_vector_alpha_col(id1).or.&
             NEFR4A%blockSize(ja,id1,'beta')&
             .ne.NEFMatB%Size_vector_beta_col(id1))&
            call mqc_error_i('NEF subspace dimensions are &
          &not conformable in mqc_nef_r4_nef_Matrix_contraction',6,&
            'Sub-space index (id1)', id1,&
            'NEFR4A%blockSize(ja,id1,alpha)',NEFR4A%blockSize(ja,id1,'alpha'),&
            'NEFMatB%Size_vector_alpha_col(id1)',NEFMatB%Size_vector_alpha_col(id1),&
            'NEFR4A%blockSize(ja,id1,beta)',NEFR4A%blockSize(ja,id1,'beta'),&
            'NEFMatB%Size_vector_beta_col(id1)',NEFMatB%Size_vector_beta_col(id1))
        end do
      end if
      if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
        do id1 = 1, NEFR4A%blockSize(la)
          if(NEFR4A%blockSize(la,id1,'alpha')&
             .ne.NEFMatB%Size_vector_alpha_col(id1).or.&
             NEFR4A%blockSize(la,id1,'beta')&
             .ne.NEFMatB%Size_vector_beta_col(id1))&
            call mqc_error_i('NEF subspace dimensions are &
          &not conformable for exchange in mqc_nef_r4_nef_Matrix_contraction',6,&
            'Sub-space index (id1)', id1,&
            'NEFR4A%blockSize(la,id1,alpha)',NEFR4A%blockSize(la,id1,'alpha'),&
            'NEFMatB%Size_vector_alpha_col(id1)',NEFMatB%Size_vector_alpha_col(id1),&
            'NEFR4A%blockSize(la,id1,beta)',NEFR4A%blockSize(la,id1,'beta'),&
            'NEFMatB%Size_vector_beta_col(id1)',NEFMatB%Size_vector_beta_col(id1))
        end do
      end if
!
      Map_Matrix_Col = NEFR4A%blockMap(la)
      Map_Matrix_Row = NEFR4A%blockMap(ka)
      call output%init(Map_Matrix_Row,initial_map_matrix_col=Map_Matrix_Col)
      do id4 = 1, NEFR4A%blockSize(la)
        do id3 = 1, NEFR4A%blockSize(ka)
          do id2 = 1, NEFR4A%blockSize(ja)
            do id1 = 1, NEFR4A%blockSize(ia)
              idA = NEFR4A%blockIndex([ia,id1],[ja,id2],[ka,id3],[la,id4])
              idB = NEFMatB%blockIndex([1,id1],[2,id2])
              idOut = NEFMatB%blockIndex([1,id3],[2,id4])
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
                if ((NEFR4A%blockSize(ia,id1,'alpha')+NEFR4A%blockSize(ia,id1,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ja,id2,'alpha')+NEFR4A%blockSize(ja,id2,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ka,id3,'alpha')+NEFR4A%blockSize(ka,id3,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(la,id4,'alpha')+NEFR4A%blockSize(la,id4,'beta')).eq.0) cycle
                if (.not.(NEFR4A%Zero_Flag(idA).or.NEFMatB%Zero_Flag(idB))) then
                  if (ka.gt.la) then
                    if (output%Zero_Flag(idOut)) then
                      output%Main_Matrix(idOut) = &
                          transpose(contraction(NEFR4A%Main_R4(idA),&
                                        NEFMatB%Main_Matrix(idB),[ia,ja]))
                    else
                      output%Main_Matrix(idOut) = output%Main_Matrix(idOut) &
                          + transpose(contraction(NEFR4A%Main_R4(idA),&
                                        NEFMatB%Main_Matrix(idB),[ia,ja]))
                    end if
                  else
                    if (output%Zero_Flag(idOut)) then
                      output%Main_Matrix(idOut) = &
                          contraction(NEFR4A%Main_R4(idA),&
                                        NEFMatB%Main_Matrix(idB),[ia,ja])
                    else
                      output%Main_Matrix(idOut) = output%Main_Matrix(idOut) &
                          + contraction(NEFR4A%Main_R4(idA),&
                                        NEFMatB%Main_Matrix(idB),[ia,ja])
                    end if
                  end if
                  output%Zero_Flag(idOut) = .false.
                end if
              end if
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
                if ((NEFR4A%blockSize(ia,id1,'alpha')+NEFR4A%blockSize(ia,id1,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(la,id2,'alpha')+NEFR4A%blockSize(la,id2,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ka,id3,'alpha')+NEFR4A%blockSize(ka,id3,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ja,id4,'alpha')+NEFR4A%blockSize(ja,id4,'beta')).eq.0) cycle
                idex = NEFR4A%blockIndex([ia,id1],[la,id2],[ka,id3],[ja,id4])
                idB = NEFMatB%blockIndex([1,id1],[2,id2])
                idOut = NEFMatB%blockIndex([1,id3],[2,id4])
                if (.not.(NEFR4A%Zero_Flag(idex).or.NEFMatB%Zero_Flag(idB))) then
                  if (output%Zero_Flag(idOut)) then
                    output%Main_Matrix(idOut) = &
                        (-1.0)*contraction(NEFR4A%Main_R4(idex),&
                                      NEFMatB%Main_Matrix(idB),[ia,la])
                  else
                    output%Main_Matrix(idOut) = output%Main_Matrix(idOut) &
                        - contraction(NEFR4A%Main_R4(idex),&
                                      NEFMatB%Main_Matrix(idB),[ia,la])
                  end if
                  output%Zero_Flag(idOut) = .false.
                end if
              end if
            end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_nef_matrix_contraction
!
!
!     PROCEDURE MQC_NEF_Matrix_Dimension_Slice
!
!
!>    \brief <b> MQC_NEF_Matrix_Dimension_Slice is used to return a slice of
!>    the input NEF Matrix object </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Dimension_Slice is used to return a slice of
!>    the input NEF Matrix object.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix) 
!>        NEF Matrix Object to be Sliced.
!>    \endverbatim
!>
!>    \param[in] rows
!>    \verbatim
!>        rows is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the remaining rows after the slicing.
!>        Default: all rows
!>    \endverbatim
!>
!>    \param[in] columns
!>    \verbatim
!>        columns is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the remaining columns after the slicing.
!>        Default: all columns
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_dimention_slice(NEFMat,rows,columns) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64),dimension(:),optional,intent(in)::rows,columns
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map_Col,New_Map_Row
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,id,idOld
      logical::full_row,full_col
!
      if(present(rows)) then
        my_row = rows
        full_row=.false.
        if(present(columns)) then
          my_col = columns
          full_col=.false.
        else
          full_col=.true.
        endIf
      else
        full_row=.true.
        if(present(columns)) then
          my_col = columns
          full_col=.false.
        else
          NEFMatOut = NEFMat
          return
        endIf
      endIf
!
      if(full_col) then
        New_Map_Col = NEFMat%Map_Matrix_Col
      else
        allocate(New_Map_Col(size(my_col)))
        do id = 1, size(my_col)
          if(my_col(id).gt.NEFMat%NumSpace_col.or.my_col(id).le.0) &
        call mqc_error_i('Column index out-of-bound in mqc_nef_matrix_dimention_slice',6,&
             'my_col(id)',my_col(id),'NEFMat%NumSpace_col',NEFMat%NumSpace_col)
          New_Map_Col(id) = NEFMat%Map_Matrix_Col(my_col(id))
        end do
      end if
!
      if(full_row) then
        New_Map_Row = NEFMat%Map_Matrix_Row
      else
        allocate(New_Map_Row(size(my_row)))
        do id = 1, size(my_row)
          if(my_row(id).gt.NEFMat%NumSpace_row.or.my_row(id).le.0) &
        call mqc_error_i('Row index out-of-bound in mqc_nef_matrix_dimention_slice',6,&
             'my_row(id)',my_row(id),'NEFMat%NumSpace_row',NEFMat%NumSpace_row)
          New_Map_Row(id) = NEFMat%Map_Matrix_Row(my_row(id))
        end do
      end if
!
      call NEFMatOut%init(New_Map_Row,&
          initial_map_matrix_col=New_Map_Col)
      if(full_row) then
        do col = 1, size(my_col)
          do row = 1, size(NEFMat%Map_Matrix_Row)
            id = NEFMatOut%blockIndex([1,row],[2,col])
            idOld = NEFMat%blockIndex([1,row],[2,my_col(col)])
            if (.not.NEFMat%Zero_Flag(idOld)) then
              NEFMatOut%Main_matrix(id) = NEFMat%Main_matrix(idOld)
              NEFMatOut%Zero_Flag(id) = .false.
            end if
          end do
        end do
      else if (full_col) then
        do row = 1, size(my_row)
          do col = 1, size(NEFMat%Map_Matrix_Col)
            id = NEFMatOut%blockIndex([1,row],[2,col])
            idOld = NEFMat%blockIndex([1,my_row(row)],[2,col])
            if (.not.NEFMat%Zero_Flag(idOld)) then
              NEFMatOut%Main_matrix(id) = NEFMat%Main_matrix(idOld)
              NEFMatOut%Zero_Flag(id) = .false.
            end if
          end do
        end do
      else
        do row = 1, size(my_row)
          do col = 1, size(my_col)
            id = NEFMatOut%blockIndex([1,row],[2,col])
            idOld = NEFMat%blockIndex([1,my_row(row)],[2,my_col(col)])
            if (.not.NEFMat%Zero_Flag(idOld)) then
              NEFMatOut%Main_matrix(id) = NEFMat%Main_matrix(idOld)
              NEFMatOut%Zero_Flag(id) = .false.
            end if
          end do
        end do
      end if
!
      end function mqc_nef_matrix_dimention_slice
!
!
!     PROCEDURE MQC_NEF_VECTOR_DIMENSION
!
!>    \brief <b> MQC_NEF_VECTOR_Dimension is a function that returns the size of a dimension
!>    MQC NEF VECTOR object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_VECTOR_Dimension is a function that returns the size of a dimension or subdimension 
!>    of an MQC NEF VECTOR object. The required argument label specifies which block is queried. 
!>    The optional arguments axis, subaxis and spin specifies which dimension is returned. 
!>    The following options are available:
!>
!>    1. .not.present(SubAxis) returns the dimension of sub-spaceS (default).
!>    2. SubAxis = 'n' selects the n-th dimension of sub-spaces.
!>
!>    1. Spin = returns the 'alpha' dimension of the n-th dimension of  sub-spaces in <Axis> (default).
!>    2. Spin = returns the 'beta' dimension of the n-th dimension of  sub-spaces in <Axis>.
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] MatIn
!>    \verbatim
!>        VecIn is Class(MQC_NEF_VECTOR)
!>        The name of the MQC NEF VECTOR variable.
!>    \endverbatim
!>
!>    \param[in] SubAxis
!>    \verbatim
!>        SubAxis is Integer(kind=int64),Optional
!>       .not.present(SubAxis) returns the dimension of sub-spaces(default).
!>       SubAxis = 'n' selects the n-th dimension of sub-spaces.
!>    \endverbatim
!>
!>    \param[in] 
!>    \verbatim Spin
!>        Spin is Character(Len=*),optional
!>        The spin block to query.
!>        = 'alpha','a': alpha axis dimension.
!>        = 'beta','b':  beta axis dimension.
!>        Returns an error if present while subaxis is not present
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      function mqc_nef_vector_dimension(VecIn,subaxis,spin) result(dimBlock)
!
      implicit none
      class(mqc_nef_vector),intent(in)::VecIn
      integer(kind=int64)::dimBlock
      Character(Len=*),intent(in),optional::spin
      integer(kind=int64),intent(in),optional::subaxis
      Character(Len=64)::my_spin
      integer(kind=int64)::my_sub_axis
      logical::do_sub
!
      if(present(subaxis)) then
        my_sub_axis = subaxis
        do_sub=.true.
      else
        do_sub=.false.
      endIf
      if(present(spin)) then
        if (.not.do_sub) &
           call mqc_error('Spin specified without subaxis in MQC_NEF_Vector_Dimension',6)
        call string_change_case(spin,'l',my_spin)
      end if
!     
      if (do_sub) then
        select case(my_spin)
        case('alpha','a')
            if (my_sub_axis.gt.VecIn%NumSpace) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_Vector_Dimension',6,&
                'VecIn%NumSpace',VecIn%NumSpace,'my_sub_axis',my_sub_axis)
            dimBlock = VecIn%Size_vector_alpha(my_sub_axis)
        case('beta','b')
            if (my_sub_axis.gt.VecIn%NumSpace) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_Vector_Dimension',6,&
                'VecIn%NumSpace_col',VecIn%NumSpace,'my_sub_axis',my_sub_axis)
            dimBlock = VecIn%Size_vector_beta(my_sub_axis)
        case default
          call mqc_error_A('Spin label not valid in mqc_nef_vector_dimension', 6, &
               'my_spin', my_spin )
        end select
      else
          dimBlock = VecIn%NumSpace
      end if
      end function mqc_nef_vector_dimension
!
!     PROCEDURE MQC_NEF_MATRIX_DIMENSION
!
!>    \brief <b> MQC_NEF_MATRIX_Dimension is a function that returns the size of a dimension
!>    MQC NEF MATRIX object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_MATRIX_Dimension is a function that returns the size of a dimension or subdimension 
!>    of an MQC NEF MATRIX object. The required argument label specifies which block is queried. 
!>    The optional arguments axis, subaxis and spin specifies which dimension is returned. 
!>    The following options are available:
!>
!>    1. Axis = '1' selects the first dimension (column)(default).
!>    2. Axis = '2' selects the second dimension (row).
!>
!>    1. .not.present(SubAxis) returns the dimension of sub-spaces in <Axis> (default).
!>    2. SubAxis = 'n' selects the n-th dimension of sub-spaces in <Axis>.
!>
!>    1. Spin = returns the 'alpha' dimension of the n-th dimension of  sub-spaces in <Axis> (default).
!>    2. Spin = returns the 'beta' dimension of the n-th dimension of  sub-spaces in <Axis>.
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] MatIn
!>    \verbatim
!>        MatIn is Class(MQC_NEF_MATRIX)
!>        The name of the MQC NEF MATRIX variable.
!>    \endverbatim
!>
!>    \param[in] Axis
!>    \verbatim
!>        Axis is Integer(kind=int64),Optional
!>        The spin block axis to return.
!>        = 1: number of first dimension (column)(default).
!>        = 2: number of second dimension (row).
!>    \endverbatim
!>
!>    \param[in] SubAxis
!>    \verbatim
!>        SubAxis is Integer(kind=int64),Optional
!>       .not.present(SubAxis) returns the dimension of sub-spaces in <Axis> (default).
!>       SubAxis = 'n' selects the n-th dimension of sub-spaces in <Axis>.
!>    \endverbatim
!>
!>    \param[in] 
!>    \verbatim Spin
!>        Spin is Character(Len=*),optional
!>        The spin block to query.
!>        = 'alpha','a': alpha axis dimension.
!>        = 'beta','b':  beta axis dimension.
!>        Returns an error if present while subaxis is not present
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      function mqc_nef_matrix_dimension(MatIn,axis,subaxis,spin) result(dimBlock)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::MatIn
      integer(kind=int64)::dimBlock
      Character(Len=*),intent(in),optional::spin
      integer(kind=int64),intent(in),optional::axis,subaxis
      Character(Len=64)::my_spin
      integer(kind=int64)::my_axis,my_sub_axis
      logical::do_sub
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      endIf
      if(my_axis.le.0.or.my_axis.gt.2) &
        call mqc_error_I('Unknown axis specified in MQC_NEF_Matrix_Dimension',6,'my_axis',my_axis)

      if(present(subaxis)) then
        my_sub_axis = subaxis
        do_sub=.true.
      else
        do_sub=.false.
      endIf
      if(present(spin)) then
        if (.not.do_sub) &
           call mqc_error('Spin specified without subaxis in MQC_NEF_Matrix_Dimension',6)
        call string_change_case(spin,'l',my_spin)
      end if
!     
      if (do_sub) then
        select case(my_spin)
        case('alpha','a')
          if(my_axis.eq.1) then
            if (my_sub_axis.gt.MatIn%NumSpace_col) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_Matrix_Dimension',6,&
                'MatIn%NumSpace_col',MatIn%NumSpace_col,'my_sub_axis',my_sub_axis)
            dimBlock = MatIn%Size_vector_alpha_col(my_sub_axis)
          else if(my_axis.eq.2) then
            if (my_sub_axis.gt.MatIn%NumSpace_row) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_Matrix_Dimension',6,&
                'MatIn%NumSpace_row',MatIn%NumSpace_row,'my_sub_axis',my_sub_axis)
            dimBlock = MatIn%Size_vector_alpha_row(my_sub_axis)
          end if
        case('beta','b')
          if(my_axis.eq.1) then
            if (my_sub_axis.gt.MatIn%NumSpace_col) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_Matrix_Dimension',6,&
                'MatIn%NumSpace_col',MatIn%NumSpace_col,'my_sub_axis',my_sub_axis)
            dimBlock = MatIn%Size_vector_beta_col(my_sub_axis)
          else if(my_axis.eq.2) then
            if (my_sub_axis.gt.MatIn%NumSpace_row) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_Matrix_Dimension',6,&
                'MatIn%NumSpace_row',MatIn%NumSpace_row,'my_sub_axis',my_sub_axis)
            dimBlock = MatIn%Size_vector_beta_row(my_sub_axis)
          end if
        case default
          call mqc_error_A('Spin label not valid in mqc_nef_matrix_dimension', 6, &
               'my_spin', my_spin )
        end select
      else
        if(my_axis.eq.1) then
          dimBlock = MatIn%NumSpace_col
        else if(my_axis.eq.2) then
          dimBlock = MatIn%NumSpace_row
        end if
      end if
    end function mqc_nef_matrix_dimension
!
!
!     PROCEDURE MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4
!
!>    \brief <b> MQC_NEF_R4_NEF_R4_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC NEF R4 </b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_NEF_R4_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC NEF R4 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4) 
!>        The second MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] Index1
!>    \verbatim
!>        Index1 is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIs1 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>        Ex: If Index1=[3,2,4,1], the third and second indexes of ERIs1 are contracted.
!>            The fourth and first are multiplied to form the first and second indices
!>            of the result matrix, respectivelly.
!>    \endverbatim
!>
!>    \param[in] Index2
!>    \verbatim
!>        Index2 is Integer(kind=int64),dimension(4),optional
!>        Contain the indexes of ERIs2 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>
!>        Ex: index1=[3,2,4,1] and index2=[2,1,4,3]
!>        sum_{i,j} R4A_{pjiq} R4B_{jisr} = R4Out_{pqrs}
!>
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_nef_r4_partial_contraction_2nef_r4(NEFR4A,NEFR4B,index1,index2) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_r4)::NEFR4Out
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      Integer(kind=int64)::idiab,idjab,A_index,B_index,out_index,N1234
!
      if(present(index1)) then
        ia=index1(1)
        ja=index1(2)
        ka=index1(3)
        la=index1(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ka', ka,'la',la)
        my_index1 = index1
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_index1 = [1,2,3,4]
      endif
      if(present(index2)) then
        ib=index2(1)
        jb=index2(2)
        kb=index2(3)
        lb=index2(4)
        if (ib.eq.jb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'jb',jb)
        if (ib.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'kb',kb)
        if (ib.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'lb',lb)
        if (jb.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'jb', jb,'kb',kb)
        if (jb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'jb', jb,'lb',lb)
        if (kb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4',&
                                        6,'kb', kb,'lb',lb)
        my_index2 = index2
      else
        ib=1
        jb=2
        kb=3
        lb=4
        my_index2 = [1,2,3,4]
      endif
      if (ia.gt.4.or.ia.le.0.or.ib.gt.4.or.ib.le.0.or.&
          ja.gt.4.or.ja.le.0.or.jb.gt.4.or.jb.le.0.or.&
          ka.gt.4.or.ka.le.0.or.kb.gt.4.or.kb.le.0.or.&
          la.gt.4.or.la.le.0.or.lb.gt.4.or.lb.le.0.) Call MQC_Error(&
           'Index out of range in MQC_NEF_R4_NEF_R4_Partial_Contraction_2R4', 6)

      if (nefr4A%blockSize(ia).ne.nefr4B%blockSize(ib).or.&
          nefr4A%blockSize(ja).ne.nefr4B%blockSize(jb)) &
          call mqc_error_i('NEFR4 not correctly dimensioned in mqc_nef_r4_nef_r4_contraction_2R4',6,&
           'dim1',nefr4A%blockSize(ia)-nefr4B%blockSize(ib),&
           'dim2',nefr4A%blockSize(ja)-nefr4B%blockSize(jb))

      NEFR4Out%NumSpace_1 = NEFR4A%blockSize(ka)
      NEFR4Out%NumSpace_2 = NEFR4A%blockSize(la)
      NEFR4Out%NumSpace_3 = NEFR4B%blockSize(kb)
      NEFR4Out%NumSpace_4 = NEFR4B%blockSize(lb)
      N1234 = NEFR4Out%NumSpace_1*NEFR4Out%NumSpace_2*NEFR4Out%NumSpace_3*NEFR4Out%NumSpace_4

      allocate(NEFR4Out%Main_R4(N1234))
      allocate(NEFR4Out%Zero_Flag(N1234))
      NEFR4Out%Zero_Flag = .true.
      NEFR4Out%Map_Matrix_1 = NEFR4A%blockMap(ka)
      NEFR4Out%Map_Matrix_2 = NEFR4A%blockMap(la)
      NEFR4Out%Map_Matrix_3 = NEFR4B%blockMap(kb)
      NEFR4Out%Map_Matrix_4 = NEFR4B%blockMap(lb)
      NEFR4Out%Size_vector_alpha_1 = NEFR4A%blockSizeVec(ka,'alpha')
      NEFR4Out%Size_vector_alpha_2 = NEFR4A%blockSizeVec(la,'alpha')
      NEFR4Out%Size_vector_alpha_3 = NEFR4B%blockSizeVec(kb,'alpha')
      NEFR4Out%Size_vector_alpha_4 = NEFR4B%blockSizeVec(lb,'alpha')
      NEFR4Out%Size_vector_beta_1 = NEFR4A%blockSizeVec(ka,'beta')
      NEFR4Out%Size_vector_beta_2 = NEFR4A%blockSizeVec(la,'beta')
      NEFR4Out%Size_vector_beta_3 = NEFR4B%blockSizeVec(kb,'beta')
      NEFR4Out%Size_vector_beta_4 = NEFR4B%blockSizeVec(lb,'beta')

      do id1 = 1, nefr4Out%blockSize(1)
        do id2 = 1, nefr4Out%blockSize(2)
          do id3 = 1, nefr4Out%blockSize(3)
            do id4 = 1, nefr4Out%blockSize(4)
              out_index = NEFR4Out%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
              do idiab = 1, nefr4B%blockSize(ib)
                do idjab = 1, nefr4B%blockSize(jb)
                  A_index = NEFR4A%blockIndex([ka,id1],[la,id2],[ia,idiab],[ja,idjab])
                  B_index = NEFR4B%blockIndex([kb,id3],[lb,id4],[ib,idiab],[jb,idjab])
                  if (NEFR4A%Zero_Flag(A_index).or.NEFR4B%Zero_Flag(B_index)) cycle
                  if(NEFR4Out%Zero_Flag(out_index)) then
                    NEFR4Out%Main_R4(out_index) =  &
                      mqc_eri_eri_partial_contraction_2eri(&
                      NEFR4A%Main_R4(A_index),NEFR4B%Main_R4(B_index),my_index1,my_index2)
                    NEFR4Out%Zero_Flag(out_index) = .false.
                  else
                    NEFR4Out%Main_R4(out_index) =  NEFR4Out%Main_R4(out_index)+&
                      mqc_eri_eri_partial_contraction_2eri(&
                      NEFR4A%Main_R4(A_index),NEFR4B%Main_R4(B_index),my_index1,my_index2)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_nef_r4_partial_contraction_2nef_r4
!
!
!     PROCEDURE MQC_NEF_R4_ERI_Partial_Contraction_2R4
!
!>    \brief <b> MQC_NEF_R4_ERI_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC 2ERI </b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_ERI_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC ERI
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] ERIB
!>    \verbatim
!>        ERIB is type(mqc_twoeris) 
!>        The second MQC 2ERI to be contracted.
!>    \endverbatim
!>
!>    \param[in] Index1
!>    \verbatim
!>        Index1 is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIs1 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>        Ex: If Index1=[3,2,4,1], the third and second indexes of ERIs1 are contracted.
!>            The fourth and first are multiplied to form the first and second indices
!>            of the result matrix, respectivelly.
!>    \endverbatim
!>
!>    \param[in] Index2
!>    \verbatim
!>        Index2 is Integer(kind=int64),dimension(4),optional
!>        Contain the indexes of ERIs2 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>
!>        Ex: index1=[3,2,4,1] and index2=[2,1,4,3]
!>        sum_{i,j} R4A_{pjiq} ERIB_{jisr} = R4Out_{pqrs}
!>
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_eri_partial_contraction_2nef_r4(NEFR4A,ERIB,index1,index2) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A
      type(mqc_twoeris),intent(in)::ERIB
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_r4)::NEFR4Out
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      type(mqc_scalar)::val
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      Integer(kind=int64)::idiab,idjab,A_index,B_index,out_index,N1234
!
      if(present(index1)) then
        ia=index1(1)
        ja=index1(2)
        ka=index1(3)
        la=index1(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ka', ka,'la',la)
        my_index1 = index1
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_index1 = [1,2,3,4]
      endif
      if(present(index2)) then
        ib=index2(1)
        jb=index2(2)
        kb=index2(3)
        lb=index2(4)
        if (ib.eq.jb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ib', ib,'jb',jb)
        if (ib.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ib', ib,'kb',kb)
        if (ib.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'ib', ib,'lb',lb)
        if (jb.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'jb', jb,'kb',kb)
        if (jb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'jb', jb,'lb',lb)
        if (kb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_ERI_Partial_Contraction_2R4',&
                                        6,'kb', kb,'lb',lb)
        my_index2 = index2
      else
        ib=1
        jb=2
        kb=3
        lb=4
        my_index2 = [1,2,3,4]
      endif
      if (ia.gt.4.or.ia.le.0.or.ib.gt.4.or.ib.le.0.or.&
          ja.gt.4.or.ja.le.0.or.jb.gt.4.or.jb.le.0.or.&
          ka.gt.4.or.ka.le.0.or.kb.gt.4.or.kb.le.0.or.&
          la.gt.4.or.la.le.0.or.lb.gt.4.or.lb.le.0.) Call MQC_Error(&
           'Index out of range in MQC_NEF_R4_ERI_Partial_Contraction_2R4', 6)

      if (nefr4A%blockSize(ia).ne.1.or.nefr4A%blockSize(ja).ne.1) &
          call mqc_error_i('NEFR4 not correctly dimensioned in mqc_nef_r4_eri_contraction_2R4',6,&
           'dim1',nefr4A%blockSize(ia),'dim2',nefr4A%blockSize(ja))
      if(NEFR4A%blockSize(ia,1,'alpha').ne.size(ERIB,'alpha',ib).or. &
         NEFR4A%blockSize(ia,1,'beta').ne.size(ERIB,'beta',ib)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nefNEF_R4_ERI_Partial_Contraction_2R4',6,&
        'NEFR4A%blockSize(ia,1,alpha)',NEFR4A%blockSize(ia,1,'alpha'),&
        'size(ERIB,alpha,ib)',size(ERIB,'alpha',ib),&
        'NEFR4A%blockSize(ia,1,beta)',NEFR4A%blockSize(ia,1,'beta'),&
        'size(ERIB,beta,ib)',size(ERIB,'beta',ib))
      if(NEFR4A%blockSize(ja,1,'alpha').ne.size(ERIB,'alpha',jb).or. &
         NEFR4A%blockSize(ja,1,'beta').ne.size(ERIB,'beta',jb)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_nefNEF_R4_ERI_Partial_Contraction_2R4',6,&
        'NEFR4A%blockSize(ja,1,alpha)',NEFR4A%blockSize(ja,1,'alpha'),&
        'size(ERIB,alpha,jb)',size(ERIB,'alpha',jb),&
        'NEFR4A%blockSize(ja,1,beta)',NEFR4A%blockSize(ja,1,'beta'),&
        'size(ERIB,beta,jb)',size(ERIB,'beta',jb))

      NEFR4Out%NumSpace_1 = NEFR4A%blockSize(ka)
      NEFR4Out%NumSpace_2 = NEFR4A%blockSize(la)
      NEFR4Out%NumSpace_3 = 1 
      NEFR4Out%NumSpace_4 = 1
      N1234 = NEFR4Out%NumSpace_1*NEFR4Out%NumSpace_2

      allocate(NEFR4Out%Main_R4(N1234))
      allocate(NEFR4Out%Zero_Flag(N1234))
      NEFR4Out%Zero_Flag = .true.
      NEFR4Out%Map_Matrix_1 = NEFR4A%blockMap(ka)
      NEFR4Out%Map_Matrix_2 = NEFR4A%blockMap(la)
      if(allocated(NEFR4Out%Map_Matrix_3)) deallocate(NEFR4Out%Map_Matrix_3)
      allocate(NEFR4Out%Map_Matrix_3(1))
      if(allocated(NEFR4Out%Map_Matrix_4)) deallocate(NEFR4Out%Map_Matrix_4)
      allocate(NEFR4Out%Map_Matrix_4(1))
      call NEFR4Out%Map_Matrix_3(1)%init(size(ERIB,'alpha',kb),scalar=1,lengthBeta=size(ERIB,'beta',kb))
      call NEFR4Out%Map_Matrix_4(1)%init(size(ERIB,'alpha',lb),scalar=1,lengthBeta=size(ERIB,'beta',lb))
      do id1 = 1, size(ERIB,'alpha',kb)
        val = id1
        call  NEFR4Out%Map_Matrix_3(1)%put(val,id1,'alpha')
      end do
      do id1 = 1, size(ERIB,'beta',kb)
        val = id1
        call  NEFR4Out%Map_Matrix_3(1)%put(val,id1,'beta')
      end do
      do id1 = 1, size(ERIB,'alpha',lb)
        val = id1
        call  NEFR4Out%Map_Matrix_4(1)%put(val,id1,'alpha')
      end do
      do id1 = 1, size(ERIB,'beta',lb)
        val = id1
        call  NEFR4Out%Map_Matrix_4(1)%put(val,id1,'beta')
      end do
    
      NEFR4Out%Size_vector_alpha_1 = NEFR4A%blockSizeVec(ka,'alpha') 
      NEFR4Out%Size_vector_alpha_2 = NEFR4A%blockSizeVec(la,'alpha') 
      NEFR4Out%Size_vector_beta_1 = NEFR4A%blockSizeVec(ka,'beta')  
      NEFR4Out%Size_vector_beta_2 = NEFR4A%blockSizeVec(la,'beta')  
      if(allocated(NEFR4Out%Size_vector_alpha_3)) deallocate(NEFR4Out%Size_vector_alpha_3)
      allocate(NEFR4Out%Size_vector_alpha_3(1))
      if(allocated(NEFR4Out%Size_vector_alpha_4)) deallocate(NEFR4Out%Size_vector_alpha_4)
      allocate(NEFR4Out%Size_vector_alpha_4(1))
      if(allocated(NEFR4Out%Size_vector_beta_3)) deallocate(NEFR4Out%Size_vector_beta_3)
      allocate(NEFR4Out%Size_vector_beta_3(1))
      if(allocated(NEFR4Out%Size_vector_beta_4)) deallocate(NEFR4Out%Size_vector_beta_4)
      allocate(NEFR4Out%Size_vector_beta_4(1))
      NEFR4Out%Size_vector_alpha_3(1) = size(ERIB,'alpha',kb)
      NEFR4Out%Size_vector_alpha_4(1) = size(ERIB,'alpha',lb)
      NEFR4Out%Size_vector_beta_3(1) = size(ERIB,'beta',kb)
      NEFR4Out%Size_vector_beta_4(1) = size(ERIB,'beta',lb)

      do id1 = 1, nefr4Out%blockSize(1)
        do id2 = 1, nefr4Out%blockSize(2)
          out_index = NEFR4Out%blockIndex([1,id1],[2,id2],[3,1],[4,1])
          A_index = NEFR4A%blockIndex([ka,id1],[la,id2],[ia,1],[ja,1])
          if (.not.NEFR4A%Zero_Flag(A_index)) then
            NEFR4Out%Main_R4(out_index) =  &
              mqc_eri_eri_partial_contraction_2eri(&
              NEFR4A%Main_R4(A_index),ERIB,my_index1,my_index2)
            NEFR4Out%Zero_Flag(out_index) = .false.
          end if
        end do
      end do
!
      end function mqc_nef_r4_eri_partial_contraction_2nef_r4
!
!
!     PROCEDURE MQC_ERI_NEF_R4_Partial_Contraction_2R4
!
!>    \brief <b> MQC_ERI_NEF_R4_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC 2ERI </b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_ERI_NEF_R4_Contraction_2R4 is used to return the contraction of 
!>    MQC NEF R4 with a MQC ERI
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] ERIA
!>    \verbatim
!>        ERIA is type(mqc_twoeris) 
!>        The second MQC 2ERI to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] Index1
!>    \verbatim
!>        Index1 is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIs1 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>        Ex: If Index1=[3,2,4,1], the third and second indexes of ERIs1 are contracted.
!>            The fourth and first are multiplied to form the first and second indices
!>            of the result matrix, respectivelly.
!>    \endverbatim
!>
!>    \param[in] Index2
!>    \verbatim
!>        Index2 is Integer(kind=int64),dimension(4),optional
!>        Contain the indexes of ERIs2 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>
!>        Ex: index1=[3,2,4,1] and index2=[2,1,4,3]
!>        sum_{i,j} R4A_{pjiq} ERIB_{jisr} = R4Out_{pqrs}
!>
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_eri_nef_r4_partial_contraction_2nef_r4(ERIA,NEFR4B,index1,index2) result(NEFR4Out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4B
      type(mqc_twoeris),intent(in)::ERIA
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_r4)::NEFR4Out
      type(mqc_scalar)::val
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      Integer(kind=int64)::idiab,idjab,A_index,B_index,out_index,N1234
!
      if(present(index1)) then
        ia=index1(1)
        ja=index1(2)
        ka=index1(3)
        la=index1(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ka', ka,'la',la)
        my_index1 = index1
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_index1 = [1,2,3,4]
      endif
      if(present(index2)) then
        ib=index2(1)
        jb=index2(2)
        kb=index2(3)
        lb=index2(4)
        if (ib.eq.jb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'jb',jb)
        if (ib.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'kb',kb)
        if (ib.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'ib', ib,'lb',lb)
        if (jb.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'jb', jb,'kb',kb)
        if (jb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'jb', jb,'lb',lb)
        if (kb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_ERI_NEF_R4_Partial_Contraction_2R4',&
                                        6,'kb', kb,'lb',lb)
        my_index2 = index2
      else
        ib=1
        jb=2
        kb=3
        lb=4
        my_index2 = [1,2,3,4]
      endif
      if (ia.gt.4.or.ia.le.0.or.ib.gt.4.or.ib.le.0.or.&
          ja.gt.4.or.ja.le.0.or.jb.gt.4.or.jb.le.0.or.&
          ka.gt.4.or.ka.le.0.or.kb.gt.4.or.kb.le.0.or.&
          la.gt.4.or.la.le.0.or.lb.gt.4.or.lb.le.0.) Call MQC_Error(&
           'Index out of range in MQC_ERI_NEF_R4_Partial_Contraction_2R4', 6)

      if (nefr4B%blockSize(ib).ne.1.or.nefr4B%blockSize(jb).ne.1) &
          call mqc_error_i('NEFR4 not correctly dimensioned in mqc_eri_nef_r4_contraction_2R4',6,&
           'dim1',nefr4B%blockSize(ib),'dim2',nefr4B%blockSize(jb))
      if(NEFR4B%blockSize(ib,1,'alpha').ne.size(ERIA,'alpha',ia).or. &
         NEFR4B%blockSize(ib,1,'beta').ne.size(ERIA,'beta',ia)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_ERI_NEF_R4_Partial_Contraction_2R4',6,&
        'NEFR4B%blockSize(ib,1,alpha)',NEFR4B%blockSize(ib,1,'alpha'),&
        'size(ERIA,alpha,ia)',size(ERIA,'alpha',ia),&
        'NEFR4B%blockSize(ib,1,beta)',NEFR4B%blockSize(ib,1,'beta'),&
        'size(ERIA,beta,ia)',size(ERIA,'beta',ia))
      if(NEFR4B%blockSize(jb,1,'alpha').ne.size(ERIA,'alpha',ja).or. &
         NEFR4B%blockSize(jb,1,'beta').ne.size(ERIA,'beta',ja)) &
        call mqc_error_i('NEF subspace dimensions are &
        &not conformable in mqc_ERI_NEF_R4_Partial_Contraction_2R4',6,&
        'NEFR4B%blockSize(jb,1,alpha)',NEFR4B%blockSize(jb,1,'alpha'),&
        'size(ERIA,alpha,ja)',size(ERIA,'alpha',ja),&
        'NEFR4B%blockSize(jb,1,beta)',NEFR4B%blockSize(jb,1,'beta'),&
        'size(ERIA,beta,ja)',size(ERIA,'beta',ja))

      NEFR4Out%NumSpace_1 = 1
      NEFR4Out%NumSpace_2 = 1
      NEFR4Out%NumSpace_3 = NEFR4B%blockSize(kb) 
      NEFR4Out%NumSpace_4 = NEFR4B%blockSize(lb)
      N1234 = NEFR4Out%NumSpace_3*NEFR4Out%NumSpace_4

      allocate(NEFR4Out%Main_R4(N1234))
      allocate(NEFR4Out%Zero_Flag(N1234))
      NEFR4Out%Zero_Flag = .true.
      NEFR4Out%Map_Matrix_3 = NEFR4B%blockMap(kb)
      NEFR4Out%Map_Matrix_4 = NEFR4B%blockMap(lb)
      if(allocated(NEFR4Out%Map_Matrix_1)) deallocate(NEFR4Out%Map_Matrix_1)
      allocate(NEFR4Out%Map_Matrix_1(1))
      if(allocated(NEFR4Out%Map_Matrix_2)) deallocate(NEFR4Out%Map_Matrix_2)
      allocate(NEFR4Out%Map_Matrix_2(1))
      call NEFR4Out%Map_Matrix_1(1)%init(size(ERIA,'alpha',ka),scalar=1,lengthBeta=size(ERIA,'beta',ka))
      call NEFR4Out%Map_Matrix_2(1)%init(size(ERIA,'alpha',la),scalar=1,lengthBeta=size(ERIA,'beta',la))
      do id1 = 1, size(ERIA,'alpha',ka)
        val = id1
        call  NEFR4Out%Map_Matrix_1(1)%put(val,id1,'alpha')
      end do
      do id1 = 1, size(ERIA,'beta',ka)
        val = id1
        call  NEFR4Out%Map_Matrix_1(1)%put(val,id1,'beta')
      end do
      do id1 = 1, size(ERIA,'alpha',la)
        val = id1
        call  NEFR4Out%Map_Matrix_2(1)%put(val,id1,'alpha')
      end do
      do id1 = 1, size(ERIA,'beta',la)
        val = id1
        call  NEFR4Out%Map_Matrix_2(1)%put(val,id1,'beta')
      end do
    
      NEFR4Out%Size_vector_alpha_3 = NEFR4B%blockSizeVec(kb,'alpha') 
      NEFR4Out%Size_vector_alpha_4 = NEFR4B%blockSizeVec(lb,'alpha') 
      NEFR4Out%Size_vector_beta_3 = NEFR4B%blockSizeVec(kb,'beta')  
      NEFR4Out%Size_vector_beta_4 = NEFR4B%blockSizeVec(lb,'beta')  
      if(allocated(NEFR4Out%Size_vector_alpha_1)) deallocate(NEFR4Out%Size_vector_alpha_1)
      allocate(NEFR4Out%Size_vector_alpha_1(1))
      if(allocated(NEFR4Out%Size_vector_alpha_2)) deallocate(NEFR4Out%Size_vector_alpha_2)
      allocate(NEFR4Out%Size_vector_alpha_2(1))
      if(allocated(NEFR4Out%Size_vector_beta_1)) deallocate(NEFR4Out%Size_vector_beta_1)
      allocate(NEFR4Out%Size_vector_beta_1(1))
      if(allocated(NEFR4Out%Size_vector_beta_2)) deallocate(NEFR4Out%Size_vector_beta_2)
      allocate(NEFR4Out%Size_vector_beta_2(1))
      NEFR4Out%Size_vector_alpha_1(1) = size(ERIA,'alpha',ka)
      NEFR4Out%Size_vector_alpha_2(1) = size(ERIA,'alpha',la)
      NEFR4Out%Size_vector_beta_1(1) = size(ERIA,'beta',ka)
      NEFR4Out%Size_vector_beta_2(1) = size(ERIA,'beta',la)

      do id1 = 1, nefr4Out%blockSize(3)
        do id2 = 1, nefr4Out%blockSize(4)
          out_index = NEFR4Out%blockIndex([1,1],[2,1],[3,id1],[4,id2])
          A_index = NEFR4B%blockIndex([kb,id1],[lb,id2],[ib,1],[jb,1])
          if (.not.NEFR4B%Zero_Flag(A_index)) then
            NEFR4Out%Main_R4(out_index) =  &
              mqc_eri_eri_partial_contraction_2eri(&
              ERIA,NEFR4B%Main_R4(A_index),my_index1,my_index2)
            NEFR4Out%Zero_Flag(out_index) = .false.
          end if
        end do
      end do
!
      end function mqc_eri_nef_r4_partial_contraction_2nef_r4
!
!
!     PROCEDURE MQC_NEF_R4_NEF_R4_Partial_Contraction
!
!>    \brief <b> MQC_NEF_R4_NEF_R4_Contraction is used to return the contraction of 
!>    MQC NEF R4 with a MQC NEF R4 </b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_NEF_R4_Contraction is used to return the contraction of 
!>    MQC NEF R4 with a MQC NEF R4 
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4) 
!>        The first MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] NEFR4B
!>    \verbatim
!>        NEFR4B is type(mqc_nef_r4) 
!>        The second MQC NEF R4 to be contracted.
!>    \endverbatim
!>
!>    \param[in] Index1
!>    \verbatim
!>        Index1 is integer(kind=int64),dimension(4),optional 
!>        Contain the indexes of ERIs1 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>        Ex: If Index1=[3,2,4,1], the third and second indexes of ERIs1 are contracted.
!>            The fourth and first are multiplied to form the first and second indices
!>            of the result matrix, respectivelly.
!>    \endverbatim
!>
!>    \param[in] Index2
!>    \verbatim
!>        Index2 is Integer(kind=int64),dimension(4),optional
!>        Contain the indexes of ERIs2 to be contracted.
!>        The two first integers in the object are the ones to be contracted,
!>        while the two last are the one to be multiplied.
!>
!>        Ex: index1=[3,2,4,1] and index2=[2,1,4,3]
!>        sum_{i,j} R4A_{ljik} R4B_{jilk} = MatOut_{kl}
!>
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_nef_r4_partial_contraction(NEFR4A,NEFR4B,index1,index2) result(MatOut)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A,NEFR4B
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_matrix)::MatOut
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4,mat_index,nef_bidx,nef_bidx2
!
      if(present(index1)) then
        ia=index1(1)
        ja=index1(2)
        ka=index1(3)
        la=index1(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ka', ka,'la',la)
        my_index1 = index1
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_index1 = [1,2,3,4]
      endif
      if(present(index2)) then
        ib=index2(1)
        jb=index2(2)
        kb=index2(3)
        lb=index2(4)
        if (ib.eq.jb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ib', ib,'jb',jb)
        if (ib.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ib', ib,'kb',kb)
        if (ib.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'ib', ib,'lb',lb)
        if (jb.eq.kb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'jb', jb,'kb',kb)
        if (jb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'jb', jb,'lb',lb)
        if (kb.eq.lb) Call MQC_Error_I(&
               'Duplicated index in MQC_NEF_R4_NEF_R4_Partial_Contraction',&
                                        6,'kb', kb,'lb',lb)
        my_index2 = index2
      else
        ib=1
        jb=2
        kb=3
        lb=4
        my_index2 = [1,2,3,4]
      endif
      if (ia.gt.4.or.ia.le.0.or.ib.gt.4.or.ib.le.0.or.&
          ja.gt.4.or.ja.le.0.or.jb.gt.4.or.jb.le.0.or.&
          ka.gt.4.or.ka.le.0.or.kb.gt.4.or.kb.le.0.or.&
          la.gt.4.or.la.le.0.or.lb.gt.4.or.lb.le.0.) Call MQC_Error(&
           'Index out of range in MQC_NEF_R4_NEF_R4_Partial_Contraction', 6)

      if (nefr4A%blockSize(ia).ne.nefr4B%blockSize(ib).or.&
          nefr4A%blockSize(ja).ne.nefr4B%blockSize(jb).or.&
          nefr4A%blockSize(ka).ne.nefr4B%blockSize(kb).or.&
          nefr4A%blockSize(la).ne.nefr4B%blockSize(lb))&
          call mqc_error_i('NEFR4 not correctly dimensioned in mqc_nef_r4_nef_r4_contraction',6,&
           'dim1',nefr4A%blockSize(ia)-nefr4B%blockSize(ib),&
           'dim2',nefr4A%blockSize(ja)-nefr4B%blockSize(jb),&
           'dim3',nefr4A%blockSize(ka)-nefr4B%blockSize(kb),&
           'dim4',nefr4A%blockSize(la)-nefr4B%blockSize(lb))

      Map_Matrix_Col = NEFR4A%blockMap(la)
      Map_Matrix_Row = NEFR4A%blockMap(ka)
      call MatOut%init(Map_Matrix_Row,initial_map_matrix_col=Map_Matrix_Col)
      do id3 = 1, nefr4A%blockSize(ka)
        do id4 = 1, nefr4A%blockSize(la)
          mat_index = MatOut%blockIndex([1,id3],[2,id4])
          do id1 = 1, nefr4A%blockSize(ia)
            do id2 = 1, nefr4A%blockSize(ja)
              nef_bidx = NEFR4A%blockIndex([ka,id3],[la,id4],[ia,id1],[ja,id2])
              nef_bidx2 = NEFR4B%blockIndex([kb,id3],[lb,id4],[ib,id1],[jb,id2])
              if (NEFR4A%Zero_Flag(nef_bidx).or.NEFR4B%Zero_Flag(nef_bidx2)) cycle
              if (MatOut%Zero_Flag(mat_index)) then
                MatOut%Main_Matrix(mat_index) = &
                  partialcontraction(NEFR4A%Main_R4(nef_bidx),&
                  NEFR4B%Main_R4(nef_bidx2),my_index1,my_index2)
                MatOut%Zero_Flag(mat_index) = .false.
              else
                MatOut%Main_Matrix(mat_index) = MatOut%Main_Matrix(mat_index) +&
                  partialcontraction(NEFR4A%Main_R4(nef_bidx),&
                  NEFR4B%Main_R4(nef_bidx2),my_index1,my_index2)
              end if
            end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_nef_r4_partial_contraction
!
!
!     PROCEDURE MQC_NEF_R4_Dimension
!
!>    \brief <b> MQC_NEF_R4_Dimension is a function that returns the length of a spin
!>    block axis of an MQC NEF R4 object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Dimension is a function that returns the length of a spin block axis 
!>    of an MQC NEF R4 object. The required argument label specifies which block is queried. 
!>    The optional arguments axis, subaxis and spin specifies which dimension is returned. 
!>    The following options are available:
!>
!>    1. Axis = '1' selects the first dimension (default).
!>    2. Axis = '2' selects the second dimension.
!>    3. Axis = '3' selects the third dimension.
!>    4. Axis = '4' selects the fourth dimension.
!>
!>    1. .not.present(SubAxis) returns the dimension of sub-spaces in <Axis> (default).
!>    2. SubAxis = 'n' selects the n-th dimension of sub-spaces in <Axis>.
!>
!>    1. Spin = returns the 'alpha' dimension of the n-th dimension of  sub-spaces in <Axis> (default).
!>    2. Spin = returns the 'beta' dimension of the n-th dimension of  sub-spaces in <Axis>.
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The name of the MQC NEF R4 variable.
!>    \endverbatim
!>
!>    \param[in] Axis
!>    \verbatim
!>        Axis is Integer(kind=int64),Optional
!>        The spin block axis to return.
!>        = 1: number of first dimension (default).
!>        = 2: number of second dimension.
!>        = 3: number of third dimension.
!>        = 4: number of fourth dimension.
!>    \endverbatim
!>
!>    \param[in] SubAxis
!>    \verbatim
!>        SubAxis is Integer(kind=int64),Optional
!>       .not.present(SubAxis) returns the dimension of sub-spaces in <Axis> (default).
!>       SubAxis = 'n' selects the n-th dimension of sub-spaces in <Axis>.
!>    \endverbatim
!>
!>    \param[in] 
!>    \verbatim Spin
!>        Spin is Character(Len=*),optional
!>        The spin block to query.
!>        = 'alpha','a': alpha axis dimension.
!>        = 'beta','b':  beta axis dimension.
!>        Returns an error if present while subaxis is not present
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_dimension(nefr4,axis,subaxis,spin) result(dimBlock)
!
      implicit none
      class(mqc_nef_r4),intent(in)::nefr4
      integer(kind=int64)::dimBlock
      Character(Len=*),intent(in),optional::spin
      integer(kind=int64),intent(in),optional::axis,subaxis
      Character(Len=64)::my_spin
      integer(kind=int64)::my_axis,my_sub_axis
      logical::do_sub
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      endIf
      if(my_axis.le.0.or.my_axis.gt.4) &
        call mqc_error_I('Unknown axis specified in MQC_NEF_R4_Dimension',6,'my_axis',my_axis)

      if(present(subaxis)) then
        my_sub_axis = subaxis
        do_sub=.true.
      else
        do_sub=.false.
      endIf
      if(present(spin)) then
        if (.not.do_sub) &
           call mqc_error('Spin specified without subaxis in MQC_NEF_R4_Dimension',6)
        call string_change_case(spin,'l',my_spin)
      end if
     
      if (do_sub) then
        select case(my_spin)
        case('alpha','a')
          if(my_axis.eq.1) then
            if (my_sub_axis.gt.NEFR4%NumSpace_1) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_1',NEFR4%NumSpace_1,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_alpha_1(my_sub_axis)
          else if(my_axis.eq.2) then
            if (my_sub_axis.gt.NEFR4%NumSpace_2) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_2',NEFR4%NumSpace_2,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_alpha_2(my_sub_axis)
          else if(my_axis.eq.3) then
            if (my_sub_axis.gt.NEFR4%NumSpace_3) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_3',NEFR4%NumSpace_3,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_alpha_3(my_sub_axis)
          else if(my_axis.eq.4) then
            if (my_sub_axis.gt.NEFR4%NumSpace_4) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_4',NEFR4%NumSpace_4,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_alpha_4(my_sub_axis)
          end if
        case('beta','b')
          if(my_axis.eq.1) then
            if (my_sub_axis.gt.NEFR4%NumSpace_1) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_1',NEFR4%NumSpace_1,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_beta_1(my_sub_axis)
          else if(my_axis.eq.2) then
            if (my_sub_axis.gt.NEFR4%NumSpace_2) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_2',NEFR4%NumSpace_2,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_beta_2(my_sub_axis)
          else if(my_axis.eq.3) then
            if (my_sub_axis.gt.NEFR4%NumSpace_3) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_3',NEFR4%NumSpace_3,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_beta_3(my_sub_axis)
          else if(my_axis.eq.4) then
            if (my_sub_axis.gt.NEFR4%NumSpace_4) &
               call mqc_error_I('Subaxis out of bounds in MQC_NEF_R4_Dimension',6,&
                'NEFR4%NumSpace_4',NEFR4%NumSpace_4,'my_sub_axis',my_sub_axis)
            dimBlock = NEFR4%Size_vector_beta_4(my_sub_axis)
          end if
        case default
          call mqc_error_A('Spin label not valid in mqc_nef_r4_dimension', 6, &
               'my_spin', my_spin )
        end select
      else
        if(my_axis.eq.1) then
          dimBlock = NEFR4%NumSpace_1
        else if(my_axis.eq.2) then
          dimBlock = NEFR4%NumSpace_2
        else if(my_axis.eq.3) then
          dimBlock = NEFR4%NumSpace_3
        else if(my_axis.eq.4) then
          dimBlock = NEFR4%NumSpace_4
        end if
      end if
!
      end function mqc_nef_r4_dimension
!
!
!     PROCEDURE MQC_NEF_R4_Dimension_vector
!
!>    \brief <b> MQC_NEF_R4_Dimension is a function that returns the size_vector_spin of a
!>    block axis of an MQC NEF R4 object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Dimension is a function that returns the size_vector_spin block axis 
!>    of an MQC NEF R4 object. The required argument label specifies which block is queried. 
!>    The optional arguments axis, subaxis and spin specifies which dimension is returned. 
!>    The following options are available:
!>
!>    1. Axis = '1' selects the first dimension (default).
!>    2. Axis = '2' selects the second dimension.
!>    3. Axis = '3' selects the third dimension.
!>    4. Axis = '4' selects the fourth dimension.
!>
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The name of the MQC NEF R4 variable.
!>    \endverbatim
!>
!>    \param[in] Axis
!>    \verbatim
!>        Axis is Integer(kind=int64),Optional
!>        The spin block axis to return.
!>        = 1: number of first dimension (default).
!>        = 2: number of second dimension.
!>        = 3: number of third dimension.
!>        = 4: number of fourth dimension.
!>    \endverbatim
!>
!>    \param[in] 
!>    \verbatim Spin
!>        Spin is Character(Len=*),optional
!>        The spin block to query.
!>        = 'alpha': alpha axis dimension.
!>        = 'beta':  beta axis dimension.
!>        Returns an error if present while subaxis is not present
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_dimension_vector(nefr4,axis,spin) result(dimVec)
!
      implicit none
      class(mqc_nef_r4),intent(in)::nefr4
      integer,dimension(:),allocatable::dimVec
      Character(Len=*),intent(in),optional::spin
      integer(kind=int64),intent(in),optional::axis
      Character(Len=64)::my_spin
      integer(kind=int64)::my_axis,my_sub_axis
      logical::do_sub=.false.
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      endIf
      if(my_axis.le.0.or.my_axis.gt.4) &
        call mqc_error_I('Unknown axis specified in MQC_NEF_R4_Dimension',6,'my_axis',my_axis)
      if(present(spin)) then
        my_spin = spin
      else
        my_spin = 'alpha'
      endIf

      select case(my_spin)
      case('alpha')
        if(my_axis.eq.1) dimVec = NEFR4%Size_vector_alpha_1
        if(my_axis.eq.2) dimVec = NEFR4%Size_vector_alpha_2
        if(my_axis.eq.3) dimVec = NEFR4%Size_vector_alpha_3
        if(my_axis.eq.4) dimVec = NEFR4%Size_vector_alpha_4
      case('beta')
        if(my_axis.eq.1) dimVec = NEFR4%Size_vector_beta_1
        if(my_axis.eq.2) dimVec = NEFR4%Size_vector_beta_2
        if(my_axis.eq.3) dimVec = NEFR4%Size_vector_beta_3
        if(my_axis.eq.4) dimVec = NEFR4%Size_vector_beta_4
      case default
        call mqc_error_A('Unknown spin specified in MQC_NEF_R4_Dimension',6,'my_spin',my_spin)
      end select
!
      end function mqc_nef_r4_dimension_vector
!
!
!     PROCEDURE MQC_NEF_R4_Map
!
!>    \brief <b> MQC_NEF_R4_Map is a function that returns the Map_Matrix
!>    of an MQC NEF R4 object</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Map is a function that returns the Map_Matrix
!>    of an MQC NEF R4 object
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The name of the MQC NEF R4 variable.
!>    \endverbatim
!>
!>    \param[in] Axis
!>    \verbatim
!>        Axis is Integer(kind=int64),Optional
!>        The spin block axis to return.
!>        = 1: number of first dimension (default).
!>        = 2: number of second dimension.
!>        = 3: number of third dimension.
!>        = 4: number of fourth dimension.
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_map(nefr4,axis) result(MapOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::nefr4
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::MapOut
      integer(kind=int64),intent(in),optional::axis
      integer(kind=int64)::my_axis
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      endIf
      if(my_axis.le.0.or.my_axis.gt.4) &
        call mqc_error_I('Unknown axis specified in MQC_NEF_R4_Dimension',6,'my_axis',my_axis)

      if(my_axis.eq.1) then
        MapOut = NEFR4%Map_Matrix_1
      else if(my_axis.eq.2) then
        MapOut = NEFR4%Map_Matrix_2
      else if(my_axis.eq.3) then
        MapOut = NEFR4%Map_Matrix_3
      else if(my_axis.eq.4) then
        MapOut = NEFR4%Map_Matrix_4
      end if
!
      end function mqc_nef_r4_map
!
!
!     PROCEDURE MQC_NEF_R4_Index
!
!>    \brief <b> MQC_NEF_R4_Index is a function that returns the main index
!>    of an MQC NEF R4 object based on sub-indexes</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Index is a function that returns the main index
!>    of an MQC NEF R4 object based on sub-indexes. This function is used when
!>    the main index is required, but the order of the axis is variable.
!>
!>    Example: given the four sets:
!>      set1 = [1,4]; set2 = [3,1]; set3 = [4,10] and set4 = [2,3]
!>      This function returns the correspondent combined-index for ERI(4,3,1,10)
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The name of the MQC NEF R4 variable.
!>    \endverbatim
!>
!>    \param[in] set1
!>    \verbatim
!>        set1 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!>
!>    \param[in] set2
!>    \verbatim
!>        set2 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!>
!>    \param[in] set3
!>    \verbatim
!>        set3 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!>
!>    \param[in] set4
!>    \verbatim
!>        set4 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_index(nefr4,set1,set2,set3,set4) result(indexOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::nefr4
      integer(kind=int64),dimension(2),intent(in)::set1,set2,set3,set4
      integer(kind=int64)::indexOut
      integer(kind=int64)::id1,id2,id3,id4
      integer(kind=int64)::N1,N12,N123,N1234
!
      if (set1(1).le.0.or.set1(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,&
        'set1(1)',set1(1))
      if (set2(1).le.0.or.set2(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,&
        'set2(1)',set2(1))
      if (set3(1).le.0.or.set3(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,&
        'set3(1)',set3(1))
      if (set4(1).le.0.or.set4(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,&
        'set4(1)',set4(1))
    
      id1=-1
      id2=-1
      id3=-1
      id4=-1

      if (set1(1).eq.1) id1 = set1(2)
      if (set2(1).eq.1) then
        if (id1.gt.0) call mqc_error(&
          'Axis 1 index defined more than once in mqc_nef_r4_index',6)
        id1 = set2(2)
      endIf
      if (set3(1).eq.1) then
        if (id1.gt.0) call mqc_error(&
          'Axis 1 index defined more than once in mqc_nef_r4_index',6)
        id1 = set3(2)
      endIf
      if (set4(1).eq.1) then
        if (id1.gt.0) call mqc_error(&
          'Axis 1 index defined more than once in mqc_nef_r4_index',6)
        id1 = set4(2)
      endIf
      if (set1(1).eq.2) id2 = set1(2)
      if (set2(1).eq.2) then
        if (id2.gt.0) call mqc_error(&
          'Axis 2 index defined more than once in mqc_nef_r4_index',6)
        id2 = set2(2)
      endIf
      if (set3(1).eq.2) then
        if (id2.gt.0) call mqc_error(&
          'Axis 2 index defined more than once in mqc_nef_r4_index',6)
        id2 = set3(2)
      endIf
      if (set4(1).eq.2) then
        if (id2.gt.0) call mqc_error(&
          'Axis 2 index defined more than once in mqc_nef_r4_index',6)
        id2 = set4(2)
      endIf
      if (set1(1).eq.3) id3 = set1(2)
      if (set2(1).eq.3) then
        if (id3.gt.0) call mqc_error(&
          'Axis 3 index defined more than once in mqc_nef_r4_index',6)
        id3 = set2(2)
      endIf
      if (set3(1).eq.3) then
        if (id3.gt.0) call mqc_error(&
          'Axis 3 index defined more than once in mqc_nef_r4_index',6)
        id3 = set3(2)
      endIf
      if (set4(1).eq.3) then
        if (id3.gt.0) call mqc_error(&
          'Axis 3 index defined more than once in mqc_nef_r4_index',6)
        id3 = set4(2)
      endIf
      if (set1(1).eq.4) id4 = set1(2)
      if (set2(1).eq.4) then
        if (id4.gt.0) call mqc_error(&
          'Axis 4 index defined more than once in mqc_nef_r4_index',6)
        id4 = set2(2)
      endIf
      if (set3(1).eq.4) then
        if (id4.gt.0) call mqc_error(&
          'Axis 4 index defined more than once in mqc_nef_r4_index',6)
        id4 = set3(2)
      endIf
      if (set4(1).eq.4) then
        if (id4.gt.0) call mqc_error(&
          'Axis 4 index defined more than once in mqc_nef_r4_index',6)
        id4 = set4(2)
      endIf
      N1 = NEFR4%NumSpace_4
      N12 = N1*NEFR4%NumSpace_3
      N123 = N12*NEFR4%NumSpace_2
      N1234 = N123*NEFR4%NumSpace_1
      indexOut = (id1-1)*N123+(id2-1)*N12+(id3-1)*N1+id4
      !write(*,*) 'Sets = ', set1,set2,set3,set4
      !write(*,*) 'indexIn = ', id1,id2,id3,id4
      !write(*,*) 'indexOut = ', indexOut
!
      end function mqc_nef_r4_index
!
!
!     PROCEDURE MQC_NEF_Matrix_Index
!
!>    \brief <b> MQC_NEF_Matrix_Index is a function that returns the main index
!>    of an MQC NEF Matrix object based on sub-indexes</b>
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Index is a function that returns the main index
!>    of an MQC NEF Matrix object based on sub-indexes. This function is used when
!>    the main index is required, but the order of the axis is variable.
!>
!>    Example: given the four sets:
!>      set1 = [2,4]; set2 = [1,2]
!>      This function returns the correspondent combined-index for ERI(2,4)
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The name of the MQC NEF Matrix variable.
!>    \endverbatim
!>
!>    \param[in] set1
!>    \verbatim
!>        set1 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!>
!>    \param[in] set2
!>    \verbatim
!>        set2 is Integer(kind=int64),dimension(2)
!>        A set containing axis and sub-axis-index.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_index(nefmat,set1,set2) result(indexOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::nefmat
      integer(kind=int64),dimension(2),intent(in)::set1,set2
      integer(kind=int64)::indexOut
      integer(kind=int64)::id1,id2
!
      if (set1(1).le.0.or.set1(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nefmatrix_index',6,&
        'set1(1)',set1(1))
      if (set2(1).le.0.or.set2(1).gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nefmatrix_index',6,&
        'set2(1)',set2(1))
    
      id1=-1
      id2=-1

      if (set1(1).eq.1) id1 = set1(2)
      if (set2(1).eq.1) then
        if (id1.gt.0) call mqc_error(&
          'Axis 1 index defined more than once in mqc_nefmatrix_index',6)
        id1 = set2(2)
      endIf
      if (set1(1).eq.2) id2 = set1(2)
      if (set2(1).eq.2) then
        if (id2.gt.0) call mqc_error(&
          'Axis 2 index defined more than once in mqc_nefmatrix_index',6)
        id2 = set2(2)
      endIf
      indexOut = (id1-1)*NEFMat%NumSpace_col+id2
      !write(*,*) 'Sets = ', set1,set2
      !write(*,*) 'indexIn = ', id1,id2
      !write(*,*) 'indexOut = ', indexOut
!
      end function mqc_nef_matrix_index
!
!
!
!     PROCEDURE MQC_NEF_Matrix_At
!
!>    \brief <b> MQC_NEF_Matrix_At is a function that returns the value of an element 
!>    of a MQC NEF Matrix variable</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_Matrix_At is a function that returns the value of (I,J)th element of a 
!>    (p,q) block of a MQC NEF Matrix variable as an MQC scalar.  
!>    Optional argument SpinBlock specifies the spin block of MQC NEF Matrix.
!>    are available:
!>
!>    1. SpinBlock = 'alpha' uses the alpha spin block.
!>    2. SpinBlock = 'beta' uses the beta spin block.
!>    3. SpinBlock = 'alphaBeta' uses the alpha-beta spin block.
!>    4. SpinBlock = 'betaAlpha' uses the beta-alpha spin block.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] nefM
!>    \verbatim
!>        nefM is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to return the value of the (I,J)th for block (p,q)
!>        element.
!>    \endverbatim
!>    
!>    \param[in] I  
!>    \verbatim
!>        I is Integer(kind=int64)
!>        The first dimension of the element in MQC NEF Matrix
!>          If I>0 row count is from first index
!>    \endverbatim
!>    
!>    \param[in] J  
!>    \verbatim
!>        J is Integer(kind=int64)
!>        The second dimension of the element in MQC NEF Matrix
!>          If J>0 col count is from first index
!>    \endverbatim
!>    
!>    \param[in] p 
!>    \verbatim
!>        p is Integer(kind=int64)
!>        The first excitation dimension of the element in MQC NEF Matrix
!>          If p>0 row count is from first index
!>    \endverbatim
!>    
!>    \param[in] q 
!>    \verbatim
!>        q is Integer(kind=int64)
!>        The second excitation dimension of the element in MQC NEF Matrix
!>          If q>0 row count is from first index
!>    \endverbatim
!>
!>    \param[in] spinBlockIn
!>    \verbatim
!>        spinBlockIn is Character(Len=*),Optional 
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_at(nefM,i,j,p,q,spinBlockIn) result(element)
!
      implicit none
      class(mqc_nef_matrix)::nefM
      class(*),intent(in)::i,j,p,q
      integer(kind=int64)::myI,myJ,myp,myq,nef_bidx
      type(mqc_scalar)::element
      character(len=*),optional,intent(in)::spinBlockIn
      character(len=64)::spinBlock
      real(kind=real64)::zero=0.0d0
!
      element = zero
      select type (i)
      type is (integer)
        myI = i
      type is (mqc_scalar)
        myI = i
      class default
        call mqc_error('index type not defined in MQC_NEF_Matrix_at')
      end select
      select type (j)
      type is (integer)
        myJ = j
      type is (mqc_scalar)
        myJ = j
      class default
        call mqc_error('index type not defined in MQC_NEF_Matrix_at')
      end select
      select type (p)
      type is (integer)
        myp = p
      type is (mqc_scalar)
        myp = p
      class default
        call mqc_error('index type not defined in MQC_NEF_Matrix_at')
      end select
      select type (q)
      type is (integer)
        myq = q
      type is (mqc_scalar)
        myq = q
      class default
        call mqc_error('index type not defined in MQC_NEF_Matrix_at')
      end select
!
      if(myp.gt.nefM%NumSpace_row) &
        call mqc_error('Requested p dimension out of bounds in mqc_NEF_Matrix_at')
      if(myq.gt.nefM%NumSpace_col) &
        call mqc_error('Requested q dimension out of bounds in mqc_NEF_Matrix_at')
      if(myI.gt.(nefM%Size_vector_alpha_row(myp)+nefM%Size_vector_beta_row(myp))) &
        call mqc_error('Requested I dimension out of bounds in mqc_NEF_Matrix_at')
      if(myJ.gt.(nefM%Size_vector_alpha_col(myq)+nefM%Size_vector_beta_col(myq))) &
        call mqc_error('Requested J dimension out of bounds in mqc_NEF_Matrix_at')
!  
      if(present(spinBlockIn)) then
        call string_change_case(spinBlockIn,'L',spinBlock)
      else
        spinBlock = 'alpha'
      end if
      nef_bidx = nefM%blockIndex([1,myp],[2,myq])
      if (.not.nefM%Zero_Flag(nef_bidx)) &
        element = nefM%Main_Matrix(nef_bidx)%at(myI,myJ,spinBlock)
!
      end function mqc_nef_matrix_at
!
!
!
!     PROCEDURE MQC_NEF_R4_At
!
!>    \brief <b> MQC_NEF_R4_At is a function that returns the value of an element 
!>    of a MQC NEF R4 variable</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>    MQC_NEF_R4_At is a function that returns the value of (I,J,K,L)th element of a 
!>    (p,q,r,s) block of a MQC NEF R4 variable as an MQC scalar.  
!>    Optional argument SpinBlock specifies the spin block of MQC NEF R4.
!>    Each index is represented by a set = [dI,I,p,s], where s is 1 for alpha or -1 for beta
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] nefR4
!>    \verbatim
!>        nefR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to return the value 
!>        element.
!>    \endverbatim
!>    
!>    \param[in] Set1  
!>    \verbatim
!>        set1 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>    
!>    \param[in] Set2  
!>    \verbatim
!>        set2 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>
!>    \param[in] Set3  
!>    \verbatim
!>        set3 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>
!>    \param[in] Set4  
!>    \verbatim
!>        set4 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_at(nefR4,set1,set2,set3,set4) result(element)
!
      implicit none
      integer(kind=int64) :: nef_bidx
      class(mqc_nef_r4)::nefR4
      Integer(kind=int64),dimension(4),intent(in)::set1,set2,set3,set4
      integer(kind=int64)::I,J,K,L,myp,myq,myr,mys,sp,sq,sr,ss,id1,id2,id3,id4
      integer(kind=int64)::dI,dJ,dK,dL,r,s,p,q
      type(mqc_scalar)::element
      character(len=1)::sI,sJ,sK,sL
      character(len=4)::spinBlock
      real(kind=real64)::zero=0.0d0
!
      element = zero
!
      dI = set1(1)
      dJ = set2(1)
      dK = set3(1)
      dL = set4(1)
      I = set1(2)
      J = set2(2)
      K = set3(2)
      L = set4(2)
      p = set1(3)
      q = set2(3)
      r = set3(3)
      s = set4(3)
      sp = set1(4)
      sq = set2(4)
      sr = set3(4)
      ss = set4(4)
!
      if (dI.le.0.or.dI.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,'dI',dI)
      if (dJ.le.0.or.dJ.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,'dJ',dJ)
      if (dK.le.0.or.dK.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,'dK',dK)
      if (dL.le.0.or.dL.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_index',6,'dL',dL)
!
      if (dI.eq.dJ) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dI',dI,'dJ',dJ)
      if (dI.eq.dK) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dI',dI,'dK',dK)
      if (dI.eq.dL) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dI',dI,'dL',dL)
      if (dJ.eq.dK) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dJ',dJ,'dK',dK)
      if (dJ.eq.dL) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dJ',dJ,'dL',dL)
      if (dK.eq.dL) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_index',6,'dK',dK,'dL',dL)
!
      if(I.gt.nefR4%blockSize(dI)) &
        call mqc_error('Requested I dimension out of bounds in mqc_NEF_R4_at')
      if(J.gt.nefR4%blockSize(dJ)) &
        call mqc_error('Requested J dimension out of bounds in mqc_NEF_R4_at')
      if(K.gt.nefR4%blockSize(dK)) &
        call mqc_error('Requested K dimension out of bounds in mqc_NEF_R4_at')
      if(L.gt.nefR4%blockSize(dL)) &
        call mqc_error('Requested L dimension out of bounds in mqc_NEF_R4_at')
!
      if (sp.eq.-1) then
        sI = 'b'
      else if (sp.eq.1) then
        sI = 'a'
      else
        call mqc_error_i('Requested I spin not recognized in mqc_NEF_R4_at',6,&
                       'sp',sp)
      end if
      if (sq.eq.-1) then
        sJ = 'b'
      else if (sq.eq.1) then
        sJ = 'a'
      else
        call mqc_error_i('Requested K spin not recognized in mqc_NEF_R4_at',6,&
                       'sq',sq)
      end if
      if (sr.eq.-1) then
        sK = 'b'
      else if (sr.eq.1) then
        sK = 'a'
      else
        call mqc_error_i('Requested K spin not recognized in mqc_NEF_R4_at',6,&
                       'sr',sr)
      end if
      if (ss.eq.-1) then
        sL = 'b'
      else if (ss.eq.1) then
        sL = 'a'
      else
        call mqc_error_i('Requested K spin not recognized in mqc_NEF_R4_at',6,&
                       'ss',ss)
      end if
!
      spinBlock=''
      if (dI.eq.1) then
        id1 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.1) then
        id1 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.1) then
        id1 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else 
        id1 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
      if (dI.eq.2) then
        id2 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.2) then
        id2 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.2) then
        id2 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else 
        id2 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
      if (dI.eq.3) then
        id3 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.3) then
        id3 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.3) then
        id3 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else 
        id3 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
      if (dI.eq.4) then
        id4 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.4) then
        id4 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.4) then
        id4 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else 
        id4 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
!  
      nef_bidx = nefR4%blockIndex(&
                          [dI,I],[dJ,J],[dK,K],[dL,L])
      if (.not.nefR4%Zero_Flag(nef_bidx)) &
        element = nefR4%Main_R4(nef_bidx)%at(&
                           id1,id2,id3,id4,spinBlock)
!
      end function mqc_nef_r4_at
!
!
!     PROCEDURE MQC_NEF_R4_Put
!
!
!>    \brief <b> MQC_NEF_R4_Put sets the value of the (I,J,K,L)th element
!>    of a (p,q,r,s) block of a MQC NEF R4 variable from an MQC scalar.
!>    </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Put is the write counterpart of MQC_NEF_R4_At (same
!>    tagged-axis calling convention: each index is a set = [dI,I,p,s],
!>    where s is 1 for alpha or -1 for beta). The target block
!>    (NEFR4%Main_R4(nef_bidx)) must already be allocated/sized (e.g. via
!>    its own init) before calling put -- put only writes into an
!>    existing block, it does not allocate one. Depends on MQC_TwoERIs
!>    gaining a put(element,i,j,k,l,spinBlockIn) counterpart to its
!>    existing at(i,j,k,l,spinBlockIn), mirroring MQC_SCF_Integral's
!>    put/at pair (mqc_est.F03).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] nefR4
!>    \verbatim
!>        nefR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to write the value into.
!>    \endverbatim
!>
!>    \param[in] Set1
!>    \verbatim
!>        set1 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>
!>    \param[in] Set2
!>    \verbatim
!>        set2 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>
!>    \param[in] Set3
!>    \verbatim
!>        set3 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>
!>    \param[in] Set4
!>    \verbatim
!>        set4 is Integer(kind=int64),dimension(4)
!>        A set containing axis, sub-axis-index, orbital and spin.
!>    \endverbatim
!>
!>    \param[in] element
!>    \verbatim
!>        element is Type(MQC_Scalar)
!>        The value to store.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_r4_put(nefR4,set1,set2,set3,set4,element)
!
      implicit none
      integer(kind=int64) :: nef_bidx
      class(mqc_nef_r4),intent(inOut)::nefR4
      Integer(kind=int64),dimension(4),intent(in)::set1,set2,set3,set4
      type(mqc_scalar),intent(in)::element
      integer(kind=int64)::I,J,K,L,myp,myq,myr,mys,sp,sq,sr,ss,id1,id2,id3,id4
      integer(kind=int64)::dI,dJ,dK,dL,r,s,p,q
      character(len=1)::sI,sJ,sK,sL
      character(len=4)::spinBlock
!
      dI = set1(1)
      dJ = set2(1)
      dK = set3(1)
      dL = set4(1)
      I = set1(2)
      J = set2(2)
      K = set3(2)
      L = set4(2)
      p = set1(3)
      q = set2(3)
      r = set3(3)
      s = set4(3)
      sp = set1(4)
      sq = set2(4)
      sr = set3(4)
      ss = set4(4)
!
      if (dI.le.0.or.dI.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_put',6,'dI',dI)
      if (dJ.le.0.or.dJ.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_put',6,'dJ',dJ)
      if (dK.le.0.or.dK.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_put',6,'dK',dK)
      if (dL.le.0.or.dL.gt.4) call mqc_error_i(&
        'Index out-of-bounds in mqc_nef_r4_put',6,'dL',dL)
!
      if (dI.eq.dJ) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_put',6,'dI',dI,'dJ',dJ)
      if (dI.eq.dK) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_put',6,'dI',dI,'dK',dK)
      if (dI.eq.dL) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_put',6,'dI',dI,'dL',dL)
      if (dJ.eq.dK) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_put',6,'dJ',dJ,'dK',dK)
      if (dJ.eq.dL) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_put',6,'dJ',dJ,'dL',dL)
      if (dK.eq.dL) call mqc_error_i(&
        'Duplicated index in mqc_nef_r4_put',6,'dK',dK,'dL',dL)
!
      if(I.gt.nefR4%blockSize(dI)) &
        call mqc_error('Requested I dimension out of bounds in mqc_NEF_R4_put')
      if(J.gt.nefR4%blockSize(dJ)) &
        call mqc_error('Requested J dimension out of bounds in mqc_NEF_R4_put')
      if(K.gt.nefR4%blockSize(dK)) &
        call mqc_error('Requested K dimension out of bounds in mqc_NEF_R4_put')
      if(L.gt.nefR4%blockSize(dL)) &
        call mqc_error('Requested L dimension out of bounds in mqc_NEF_R4_put')
!
      if (sp.eq.-1) then
        sI = 'b'
      else if (sp.eq.1) then
        sI = 'a'
      else
        call mqc_error_i('Requested I spin not recognized in mqc_NEF_R4_put',6,&
                       'sp',sp)
      end if
      if (sq.eq.-1) then
        sJ = 'b'
      else if (sq.eq.1) then
        sJ = 'a'
      else
        call mqc_error_i('Requested K spin not recognized in mqc_NEF_R4_put',6,&
                       'sq',sq)
      end if
      if (sr.eq.-1) then
        sK = 'b'
      else if (sr.eq.1) then
        sK = 'a'
      else
        call mqc_error_i('Requested K spin not recognized in mqc_NEF_R4_put',6,&
                       'sr',sr)
      end if
      if (ss.eq.-1) then
        sL = 'b'
      else if (ss.eq.1) then
        sL = 'a'
      else
        call mqc_error_i('Requested K spin not recognized in mqc_NEF_R4_put',6,&
                       'ss',ss)
      end if
!
      spinBlock=''
      if (dI.eq.1) then
        id1 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.1) then
        id1 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.1) then
        id1 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else
        id1 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
      if (dI.eq.2) then
        id2 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.2) then
        id2 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.2) then
        id2 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else
        id2 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
      if (dI.eq.3) then
        id3 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.3) then
        id3 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.3) then
        id3 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else
        id3 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
      if (dI.eq.4) then
        id4 = p
        spinBlock=trim(spinBlock)//trim(sI)
      else if (dJ.eq.4) then
        id4 = q
        spinBlock=trim(spinBlock)//trim(sJ)
      else if (dK.eq.4) then
        id4 = r
        spinBlock=trim(spinBlock)//trim(sK)
      else
        id4 = s
        spinBlock=trim(spinBlock)//trim(sL)
      end if
!
      nef_bidx = nefR4%blockIndex(&
                          [dI,I],[dJ,J],[dK,K],[dL,L])
      call nefR4%Main_R4(nef_bidx)%put(element,id1,id2,id3,id4,spinBlock)
      nefR4%Zero_Flag(nef_bidx) = .false.
!
      end subroutine mqc_nef_r4_put
!
!
!     PROCEDURE MQC_NEF_Vector_Dimension_Expand
!
!
!>    \brief <b> MQC_NEF_Vector_Dimension_Expand is used to return an expanded
!>     NEF Matrix object based on a NEF Vector object. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Dimension_Expand is used to return an expanded
!>     NEF Matrix object based on a NEF Vector object.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFV
!>    \verbatim
!>        NEFV is type(mqc_nef_vector) 
!>        NEF Vector Object to be expanded.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        spinType is Character(Len=*),Optional 
!>        spinType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin 
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_dimention_expand(&
         NEFV,map_matrix_row,map_matrix_col,spinType) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_vector),intent(in)::NEFV
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,intent(in)::map_matrix_row,map_matrix_col
      character(len=*),intent(in),optional::spinType
      character(len=64)::mySpinType
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,i,j,id
!
      if(present(spinType)) then
        mySpinType = spinType 
      else
        mySpinType = 'all' 
      end if
!
      call NEFMatOut%init(map_matrix_row,initial_map_matrix_col=map_matrix_col)
      select case(mySpinType)
      case ('singlet')
        do row = 1, NEFMatOut%NumSpace_row
          do col = 1, NEFMatOut%NumSpace_col
            id = NEFMatOut%blockIndex([1,row],[2,col])
            if (NEFV%Zero_Flag(id)) cycle
            call NEFMatOut%Main_Matrix(id)%init(&
                   NEFMatOut%Size_vector_alpha_row(row),&
                   NEFMatOut%Size_vector_beta_row(row),&
                   nAlpha2=NEFMatOut%Size_vector_alpha_col(col),&
                   nBeta2=NEFMatOut%Size_vector_beta_col(col))
            NEFMatOut%Zero_Flag(id) = .false.
            do i = 1, NEFMatOut%Size_vector_alpha_row(row)
              do j = 1, NEFMatOut%Size_vector_alpha_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','a',row,col,'singlet')),'alpha'),&
                    i,j,'alpha')
              end do
            end do
            do i = 1, nefMatOut%Size_vector_beta_row(row)
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','b',row,col,'singlet')),'beta'),&
                    i,j,'beta')
              end do
            end do
          end do
        end do
      case ('triplet')
        do row = 1, NEFMatOut%NumSpace_row
          do col = 1, NEFMatOut%NumSpace_col
            id = nefMatOut%blockIndex([1,row],[2,col])
            if (NEFV%Zero_Flag(id)) cycle
            call NEFMatOut%Main_Matrix(id)%init(&
                   NEFMatOut%Size_vector_alpha_row(row),&
                   NEFMatOut%Size_vector_beta_row(row),&
                   nAlpha2=NEFMatOut%Size_vector_alpha_col(col),&
                   nBeta2=NEFMatOut%Size_vector_beta_col(col))
            NEFMatOut%Zero_Flag(id) = .false.
            do i = 1, nefMatOut%Size_vector_alpha_row(row)
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','b',row,col,'triplet')),'alpha'),&
                    i,j,'alphabeta')
              end do
            end do
            do i = 1, nefMatOut%Size_vector_beta_row(row)
              do j = 1, nefMatOut%Size_vector_alpha_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','a',row,col,'triplet')),'beta'),&
                    i,j,'betaalpha')
              end do
            end do
          end do
        end do
      case ('all','full')
        do row = 1, NEFMatOut%NumSpace_row
          do col = 1, NEFMatOut%NumSpace_col
            id = nefMatOut%blockIndex([1,row],[2,col])
            if (NEFV%Zero_Flag(id)) cycle
            call NEFMatOut%Main_Matrix(id)%init(&
                   NEFMatOut%Size_vector_alpha_row(row),&
                   NEFMatOut%Size_vector_beta_row(row),&
                   nAlpha2=NEFMatOut%Size_vector_alpha_col(col),&
                   nBeta2=NEFMatOut%Size_vector_beta_col(col))
            NEFMatOut%Zero_Flag(id) = .false.
            do i = 1, nefMatOut%Size_vector_alpha_row(row)
              do j = 1, nefMatOut%Size_vector_alpha_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','a',row,col,'full')),'alpha'),&
                    i,j,'alpha')
              end do
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','b',row,col,'full')),'alpha'),&
                    i,j,'alphabeta')
              end do
            end do
            do i = 1, nefMatOut%Size_vector_beta_row(row)
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','b',row,col,'full')),'beta'),&
                    i,j,'beta')
              end do
              do j = 1, nefMatOut%Size_vector_alpha_col(col)
                call NEFMatOut%Main_Matrix(id)%put(&
                  NEFV%Main_Vector(id)%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','a',row,col,'full')),'beta'),&
                    i,j,'betaalpha')
              end do
            end do
          end do
        end do
      case default
        call mqc_error_a('mqc_nef_vector_dimention_excpanded: Unrecognized excitation type',6,&
          'spinType',spinType)
      end select  
!
      end function mqc_nef_vector_dimention_expand
!      
!     PROCEDURE MQC_NEF_Vector_Insert
!
!
!>    \brief <b> MQC_NEF_Vector_Insert is used to insert eigenvalues into a NEF vector </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Insert is used to insert eigenvalues into a NEF vector
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] VecIn
!>    \verbatim
!>        VecIn is type(mqc_nef_vector) 
!>        VecIn is the vector you initially have, and insert the eigenvalues into
!>    \endverbatim
!>
!>    \param[in] idI
!>    \verbatim
!>            idI is the index or location that you would like to place the eigenvalues into
!>    \endverbatim!>    
!>    
!>    \param[in] ExtMap
!>    \verbatim
!>            ExtMap is the map matrix of the extra block of eigenvalues being added, it is optional
!>    \endverbatim
!
!>    \param[in] ExtVec
!>    \verbatim
!>            ExtVec is the block of extra eigenvalues you want to add, it is optional
!>    \endverbatim
!
!
!
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      Subroutine MQC_NEF_Vector_Insert(VecIn,idI,ExtMap,ExtVec)
!
!     Variable Declarations.
      implicit none
      class(mqc_nef_vector),intent(inOut)::VecIn
      type(mqc_scf_eigenvalues)::ExtI
      integer::idI
      integer::i
      integer::zero
      type(mqc_scf_eigenvalues),Optional::ExtMap,ExtVec
      type(mqc_scf_eigenvalues)::EmptyBlock
      type(mqc_scf_eigenvalues),dimension(:),allocatable::tempMain,tempMap
      logical,dimension(:),allocatable::tempZeroFlag
      logical::extIsZero
      zero=0
      allocate(tempMain(VecIn%Numspace+1),tempMap(VecIn%Numspace+1))
      allocate(tempZeroFlag(VecIn%Numspace+1))
      if(present(ExtMap)) then
        EmptyBlock=ExtMap
        if(present(ExtVec)) then
          !! Both Vec and Map are given
          ExtI=ExtVec
          extIsZero=.false.

          if ((ExtI%blockSize("alpha")).ne.(EmptyBlock%blockSize("alpha"))) &
            call mqc_error_i("Dimensionality between Map and Input Block not conformable in mqc_nef_vector_insert",&
            6,"ExtI%blockSize(alpha)",ExtI%blockSize("alpha"),"EmptyBlock%blockSize(alpha)",EmptyBlock%blockSize("alpha"))
          if ((ExtI%blockSize("beta")).ne.(EmptyBlock%blockSize("beta"))) &
            call mqc_error_i("Dimensionality between Map and Input Block not conformable in mqc_nef_vector_insert",&
            6,"ExtI%blockSize(beta)",ExtI%blockSize("beta"),"EmptyBlock%blockSize(beta)",EmptyBlock%blockSize("beta"))
        else
          !! Only Map is given
          if(EmptyBlock%type().eq."space") then
            call ExtI%init(EmptyBlock%blockSize("alpha"),zero)
          else
            call ExtI%init(EmptyBlock%blockSize("alpha"),zero,EmptyBlock%blockSize("beta"))
          end if
          extIsZero=.true.
        end if
      else
        if(present(ExtVec)) then
          !! Only Vec is given
          ExtI=ExtVec
          extIsZero=.false.
          if(ExtI%type().eq."space") then
            call EmptyBlock%init(ExtI%blockSize("alpha"),zero)
          else
            call EmptyBlock%init(ExtI%blockSize("alpha"),zero,ExtI%blockSize("beta"))
          end if
        else
          !! Nothing is given
          call mqc_error("No input given to mqc_nef_vector_insert",6)
        end if
      end if
      do i=1,idI-1
        tempMain(i)=VecIn%Main_vector(i)
        tempMap(i)=VecIn%Map_matrix(i)
        tempZeroFlag(i)=VecIn%Zero_Flag(i)
      end do
      tempMain(idI)=ExtI
      tempMap(idI)=EmptyBlock
      tempZeroFlag(idI)=extIsZero
      do i=idI+1,(VecIn%Numspace)+1
        tempMain(i)=VecIn%Main_vector(i-1)
        tempMap(i)=VecIn%Map_matrix(i-1)
        tempZeroFlag(i)=VecIn%Zero_Flag(i-1)
      end do
!
      call VecIn%init(tempMap)
      VecIn%Main_vector=tempMain
      VecIn%Zero_Flag=tempZeroFlag
!
      end subroutine MQC_NEF_Vector_Insert
!
!     PROCEDURE MQC_NEF_Matrix_Map_Insert
!
!
!>    \brief <b> MQC_NEF_Matrix_Map_Insert is used to insert eigenvalues or NEF vector into a NEF map matrix </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Map_Insert is used to insert eigenvalues or NEF vector into a NEF map matrix
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] MatIn
!>    \verbatim
!>        MatIn is type(mqc_nef_matrix) 
!>        MatIn is the map matrix you initially have, and insert the new blocks into
!>    \endverbatim
!>
!>    \param[in] idIRow
!>    \verbatim
!>            idI is the index or location on the row axis that you would like to place the block(s) into
!>    \endverbatim!>   
!>
!>    \param[in] idICol
!>    \verbatim
!>            idI is the index or location on the column axis that you would like to place the block(s) into
!>    \endverbatim!>  
!>
!>    \param[in] ExtMapC
!>    \verbatim
!>            ExtMapC is a column of a map matrix that you want to be added to MatIn, it is optional
!>    \endverbatim
!>
!>    \param[in] ExtMapR
!>    \verbatim
!>            ExtMapR is a row of a map matrix that you want to be added to MatIn, it is optional
!>    \endverbatim
!>
!>    \param[in] ExtSCFI
!>    \verbatim
!>            ExtSCFI is a block of SCF integrals, it is used to &
!>            determine dimensions of empty blocks, it is optional
!>    \endverbatim
!
!
!
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026
!
      Subroutine MQC_NEF_Matrix_Map_Insert(MatIn,idIRow,idICol,ExtSCFI,ExtMapC,ExtMapR)
!
!     Variable Declarations.
      implicit none
      integer(kind=int64) :: nef_bidx
      class(mqc_nef_matrix),intent(inOut)::MatIn
      type(mqc_scf_integral),allocatable::ExtI
      integer,optional,intent(in)::idIRow,idICol
      integer::myidIRow,myidICol
      integer::i,j,zero
      integer::newNrow,newNcol,extendRow,extendCol
      integer::costRow,costCol
      type(mqc_scf_integral)::EmptyBlockIntegral
      type(mqc_scf_integral),optional,intent(in)::ExtSCFI
      type(mqc_scf_eigenvalues),optional,intent(in)::ExtMapC,ExtMapR
      type(mqc_scf_eigenvalues)::EmptyBlockRow,EmptyBlockCol
      type(mqc_scf_integral),dimension(:),allocatable::tempMain
      type(mqc_scf_eigenvalues),dimension(:),allocatable::tempMapRow,tempMapCol
      zero = 0
!
      if (present(ExtMapR) .or. present(ExtMapC)) then
        if (present(ExtMapR)) then
          EmptyBlockRow = ExtMapR
          extendRow = 1
          newNrow = (mqc_nef_matrix_dimension(MatIn,2))+1
        else
          extendRow = 0
          newNrow = mqc_nef_matrix_dimension(MatIn,2)
        end if
        if (present(ExtMapC)) then
          EmptyBlockCol = ExtMapC
          extendCol = 1
          newNcol = (mqc_nef_matrix_dimension(MatIn,1))+1
        else
          extendCol = 0
          newNcol = mqc_nef_matrix_dimension(MatIn,1)  
        end if
!        
      allocate(tempMain(newNrow*newNcol))
      allocate(tempMapRow(newNrow))
      allocate(tempMapCol(newNcol))
      allocate(ExtI)
!
        if (present(ExtSCFI)) then
          ExtI = ExtSCFI
          if ((ExtI%blockSize("alpha")) .ne. (EmptyBlockRow%blockSize("alpha"))) &
            call mqc_error_i("Dimensionality between Map and Input Block not conformable" // &
            " in mqc_nef_matrix_map_insert",6, &
            "ExtI%blockSize(alpha)",ExtI%blockSize("alpha"), &
            "EmptyBlockRow%blockSize(alpha)",EmptyBlockRow%blockSize("alpha"))
          if ((ExtI%blockSize("beta")) .ne. (EmptyBlockRow%blockSize("beta"))) &
            call mqc_error_i("Dimensionality between Map and Input Block not conformable" // &
            " in mqc_nef_matrix_map_insert",6, &
            "ExtI%blockSize(beta)",ExtI%blockSize("beta"), &
            "EmptyBlockRow%blockSize(beta)",EmptyBlockRow%blockSize("beta"))
        else
          if (extendRow .eq. 1) then
          if (EmptyBlockRow%type() .eq. "space") then
            call ExtI%init(EmptyBlockRow%blockSize("alpha"),zero)
          else
            call ExtI%init(EmptyBlockRow%blockSize("alpha"),zero, &
                EmptyBlockRow%blockSize("beta"))
          end if
          end if
          if (extendCol .eq. 1) then
          if (EmptyBlockCol%type() .eq. "space") then
            call ExtI%init(EmptyBlockCol%blockSize("alpha"),zero)
          else
            call ExtI%init(EmptyBlockCol%blockSize("alpha"),zero, &
                EmptyBlockCol%blockSize("beta"))
                
          end if
         end if 
        end if
      else
        if (present(ExtSCFI)) then
          ExtI = ExtSCFI
          if (ExtI%type() .eq. "space") then
            call EmptyBlockRow%init(ExtI%blockSize("alpha"),zero)
            call EmptyBlockCol%init(ExtI%blockSize("alpha"),zero)
          else
            call EmptyBlockRow%init(ExtI%blockSize("alpha"),zero, &
                ExtI%blockSize("beta"))
            call EmptyBlockCol%init(ExtI%blockSize("alpha"),zero, &
                ExtI%blockSize("beta"))
          end if
        else
          call mqc_error("No input given to mqc_nef_matrix_map_insert",6)
        end if
      end if
!
      if (extendRow .eq. 1) then
        if (idIRow>1) then
          do i = 1,idIRow-1
            tempMapRow(i) = MatIn%Map_Matrix_Row(i)
          end do
        else
            tempMapRow(1) = MatIn%Map_Matrix_Row(idIRow)  
        end if
        tempMapRow(idIRow) = EmptyBlockRow
        do i = idIRow+1,newNrow
          tempMapRow(i) = MatIn%Map_Matrix_Row(i-1)
        end do
        if (extendCol .eq. 0) then
          do j = 1,newNcol
            tempMapCol(j) = MatIn%Map_Matrix_Col(j)
          end do
        end if
      end if
!

      if (extendCol .eq. 1) then
       if (idICol>1) then      
         do j = 1,idICol-1
           tempMapCol(j) = MatIn%Map_Matrix_Col(j)
         end do
       else
         tempMapCol(1)=MatIn%Map_Matrix_Col(idICol)
       end if  
         tempMapCol(idICol) = EmptyBlockCol
         do j = idICol+1,newNcol
           tempMapCol(j) = MatIn%Map_Matrix_Col(j-1)
         end do
        if (extendRow .eq. 0) then 
           do i = 1,newNRow
             tempMapRow(i) = MatIn%Map_Matrix_Row(i)
           end do
        end if   
      end if
 



  do i=1, newNrow
    do j=1,newNcol
      if (i==idIRow .or. j==idICol) then
        call EmptyBlockIntegral%init(size(tempMapRow(i),'alpha'),size(tempMapRow(i),'beta'),&
             0,'',size(tempMapCol(j),'alpha'), size(tempMapCol(j),'beta'))
        tempMain(((i-1)*newNcol)+j) = &
           EmptyBlockIntegral
      else if (i > idIRow .and. idIRow>0) then
        if (j>idICol .and. idICol>0) then
          nef_bidx = MatIn%blockIndex([1,i-1],[2,j-1])
          tempMain(((i-1)*newNcol)+j) = &
                 MatIn%Main_matrix(nef_bidx)
        else
          nef_bidx = MatIn%blockIndex([1,i-1],[2,j])
          tempMain(((i-1)*newNcol)+j) = &
                 MatIn%Main_matrix(nef_bidx)
        end if
      else
        if (j>idICol .and. idICol>0) then
          nef_bidx = MatIn%blockIndex([1,i],[2,j-1])
          tempMain(((i-1)*newNcol)+j) = &
                 MatIn%Main_matrix(nef_bidx)
        else
          nef_bidx = MatIn%blockIndex([1,i],[2,j])
          tempMain(((i-1)*newNcol)+j) = &
                 MatIn%Main_matrix(nef_bidx)
        end if
      end if
    end do 
  end do   
      call MatIn%init(tempMapRow, initial_map_matrix_col=tempMapCol)
      MatIn%main_matrix = tempMain
!
      end subroutine MQC_NEF_Matrix_Map_Insert
!    
!     PROCEDURE MQC_NEF_Matrix_Main_Insert
!
!
!>    \brief <b> MQC_NEF_Matrix_Main_Insert is used to insert SCF integral block(s) into a NEF matrix </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Main_Insert is used to insert SCF integral block(s) into a NEF matrix
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] MatIn
!>    \verbatim
!>        MatIn is type(mqc_nef_vector) 
!>        MatIn is the map matrix you initially have, and insert the new blocks into
!>    \endverbatim
!>
!>    \param[in] idIRow
!>    \verbatim
!>            idI is the index or location on the row axis that you would like to place the block(s) into
!>    \endverbatim!>   
!>
!>    \param[in] idICol
!>    \verbatim
!>            idI is the index or location on the column axis that you would like to place the block(s) into
!>    \endverbatim!>  
!>
!>    \param[in] ExtMapC
!>    \verbatim
!>            ExtMapC is a column of a map matrix that you want to be added to MatIn, it is optional
!>    \endverbatim
!>
!>    \param[in] ExtMapR
!>    \verbatim
!>            ExtMapR is a row of a map matrix that you want to be added to MatIn, it is optional
!>    \endverbatim
!>
!>    \param[in] ExtSCFI
!>    \verbatim
!>            ExtSCFI is a block of SCF integrals, it is used to &
!>            determine dimensions of empty blocks, it is optional
!>    \endverbatim
!
!
!
!
!     Authors:
!     ========
!>    \author N. J. Woolery
!>    \date 2026

!      subroutine MQC_NEF_Matrix_Main_Insert
              
!
!
!     PROCEDURE MQC_NEF_Matrix_Dimension_Collapse
!
!
!>    \brief <b> MQC_NEF_Matrix_Dimension_Collapse is used to return a collapse
!>    two indexes of the input NEF Matrix object into one. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Dimension_Collapse is used to return a collapse two indexes
!>    of the input NEF Matrix object into one.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFM
!>    \verbatim
!>        NEFM is type(mqc_nef_matrix) 
!>        NEF Matrix Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        spinType is Character(Len=*),Optional 
!>        spinType modiffies the excitations considered
!>          Options:
!>            'singlet' => return same spin 
!>            'triplet  => return opposite spin
!>            'all'     => return both types (default)
!>    \endverbatim
!>
!>    \param[in] swapOrder
!>    \verbatim
!>        swapOrder is logical,Optional
!>        swapOrder defines the excitation order
!>          Options:
!>            .false.   => row -> col (default) 
!>            .true.    => col -> row 
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_dimention_collapse(NEFM,spinType,swapOrder) result(NEFVecOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFM
      character(len=*),intent(in),optional::spinType
      logical,intent(in),optional::swapOrder
      character(len=64)::mySpinType
      type(mqc_nef_vector)::NEFVecOut
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,i,j,id,outId
      logical::mySwapOrder
!
      if(present(spinType)) then
        mySpinType = spinType 
      else
        mySpinType = 'all' 
      end if
!
      if(present(swapOrder)) then
        mySwapOrder = swapOrder 
      else
        mySwapOrder = .false.
      end if
!
      call NEFVecOut%init(NEFM%collapseMap(mySpinType,swapOrder=mySwapOrder))
      select case(mySpinType)
      case ('singlet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            if (NEFM%Zero_Flag(id)) cycle
            if (mySwapOrder) then
              outId = (col-1)*NEFM%NumSpace_row + row
            else
              outId = id
            end if
            call NEFVecOut%Main_vector(outId)%init(NEFVecOut%Size_vector_alpha(outId),&
                   lengthBeta=NEFVecOut%Size_vector_beta(outId))
            NEFVecOut%Zero_Flag(outId) = .false.
            do i = 1, nefM%Size_vector_alpha_row(row)
              do j = 1, nefM%Size_vector_alpha_col(col)
                call NEFVecOut%Main_vector(outId)%put(&
                  NEFM%at(i,j,row,col,'alpha'),&
                  int(NEFM%collapseIndex(i,j,'a','a',row,col,'singlet',swapOrder=mySwapOrder)),&
                  'alpha')
              end do
            end do
            do i = 1, nefM%Size_vector_beta_row(row)
              do j = 1, nefM%Size_vector_beta_col(col)
                call NEFVecOut%Main_vector(outId)%put(&
                  NEFM%at(i,j,row,col,'beta'),&
                  int(NEFM%collapseIndex(i,j,'b','b',row,col,'singlet',swapOrder=mySwapOrder)),&
                  'beta')
              end do
            end do
          end do
        end do
      case ('triplet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            if (NEFM%Zero_Flag(id)) cycle
            if (mySwapOrder) then
              outId = (col-1)*NEFM%NumSpace_row + row
            else
              outId = id
            end if
            call NEFVecOut%Main_vector(outId)%init(NEFVecOut%Size_vector_alpha(outId),&
                   lengthBeta=NEFVecOut%Size_vector_beta(outId))
            NEFVecOut%Zero_Flag(outId) = .false.
            do i = 1, nefM%Size_vector_alpha_row(row)
              do j = 1, nefM%Size_vector_beta_col(col)
                call NEFVecOut%Main_vector(outId)%put(&
                  NEFM%at(i,j,row,col,'alphabeta'),&
                  int(NEFM%collapseIndex(i,j,'a','b',row,col,'triplet',swapOrder=mySwapOrder)),&
                  'alpha')
              end do
            end do
            do i = 1, nefM%Size_vector_beta_row(row)
              do j = 1, nefM%Size_vector_alpha_col(col)
                call NEFVecOut%Main_vector(outId)%put(&
                  NEFM%at(i,j,row,col,'betaalpha'),&
                  int(NEFM%collapseIndex(i,j,'b','a',row,col,'triplet',swapOrder=mySwapOrder)),&
                  'beta')
              end do
            end do
          end do
        end do
      case ('all')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            if (NEFM%Zero_Flag(id)) cycle
            if (mySwapOrder) then
              outId = (col-1)*NEFM%NumSpace_row + row
            else
              outId = id
            end if
            call NEFVecOut%Main_vector(outId)%init(NEFVecOut%Size_vector_alpha(outId),&
                   lengthBeta=NEFVecOut%Size_vector_beta(outId))
            NEFVecOut%Zero_Flag(outId) = .false.
            do i = 1, nefM%Size_vector_alpha_row(row)
              do j = 1, nefM%Size_vector_alpha_col(col)
                call NEFVecOut%Main_vector(outId)%put(&
                  NEFM%at(i,j,row,col,'alpha'),&
                  int(NEFM%collapseIndex(i,j,'a','a',row,col,'all',swapOrder=mySwapOrder)),&
                  'alpha')
              end do
              do j = 1, nefM%Size_vector_beta_col(col)
                call NEFVecOut%Main_vector(outId)%put(&
                  NEFM%at(i,j,row,col,'alphabeta'),&
                  int(NEFM%collapseIndex(i,j,'a','b',row,col,'all',swapOrder=mySwapOrder)),&
                  'alpha')
              end do
            end do
            do i = 1, nefM%Size_vector_beta_row(row)
              do j = 1, nefM%Size_vector_beta_col(col)
                call NEFVecOut%Main_vector(outId)%put(&
                  NEFM%at(i,j,row,col,'beta'),&
                  int(NEFM%collapseIndex(i,j,'b','b',row,col,'all',swapOrder=mySwapOrder)),&
                  'beta')
              end do
              do j = 1, nefM%Size_vector_alpha_col(col)
                call NEFVecOut%Main_vector(outId)%put(&
                  NEFM%at(i,j,row,col,'betaalpha'),&
                  int(NEFM%collapseIndex(i,j,'b','a',row,col,'all',swapOrder=mySwapOrder)),&
                  'beta')
              end do
            end do
          end do
        end do
      case default
        call mqc_error_a('mqc_nef_matrix_dimention_collapse: Unrecognized excitation type',6,&
          'spinType',spinType)
      end select  
!
      end function mqc_nef_matrix_dimention_collapse
!
!
!     PROCEDURE MQC_NEF_R4_Dimension_Collapse
!
!
!>    \brief <b> MQC_NEF_R4_Dimension_Collapse is used to return a collapse
!>    two indexes of the input NEF R4 object into one. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Dimension_Collapse is used to return a collapse two indexes
!>    of the input NEF R4 object into one.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFr4 is type(mqc_nef_r4) 
!>        NEF R4 Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        spinType is Character(Len=*),Optional 
!>        spinType modiffies the excitations considered
!>          Options:
!>            'singlet' => return same spin 
!>            'triplet  => return opposite spin
!>            'all'     => return both types (default)
!>    \endverbatim
!>
!>    \param[in] spinTypeKet
!>    \verbatim
!>        spinTypeKet is Character(Len=*),Optional 
!>        spinType modiffies the excitations considered over ket
!>          Options:
!>            'singlet' => return same spin 
!>            'triplet  => return opposite spin
!>            'all'     => return both types (default)
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Order is integer(kind=int64),dimension(4),intent(in),optional
!>        Order contains the indexes to be collapsed following:
!>        Order = [i,j,k,l] => i->j = p and  k->l = q
!>        (default [1,2,3,4])
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_dimention_collapse(NEFr4,spinType,spinTypeKet,Order) result(NEFMatOut)
!
      implicit none
      integer(kind=int64) :: nef_bidx
      class(mqc_nef_r4),intent(in)::NEFr4
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      integer(kind=int64),dimension(4)::my_order
      type(mqc_nef_matrix)::NEFMatOut
      character(len=64)::mySpinType,mySpinTypeKet
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,i,j,k,l,id,idij,idlk,fullid,finalid
      integer(kind=int64)::ia,ja,ka,la
      integer(kind=int64)::dimi,dimj,dimk,diml
!
      if(present(spinType)) then
        mySpinType = spinType 
        if(present(spinTypeKet)) then
          mySpinTypeKet = spinTypeKet 
        else
          mySpinTypeKet = spinType 
        end if
      else
        if(present(spinTypeKet)) then
          mySpinTypeKet = spinTypeKet 
          mySpinType = spinTypeKet
        else
          mySpinTypeKet = 'all'
          mySpinType = 'all' 
        end if
      end if
!
      if(present(Order)) then
        ia=Order(1)
        ja=Order(2)
        ka=Order(3)
        la=Order(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ka', ka,'la',la)
        my_order = Order
      else
        ia=1
        ja=2
        ka=3
        la=4
        my_order = [1,2,3,4]
      endif
!
      call NEFMatOut%init(NEFr4%collapseMap(mySpinType,[ia,ja]),&
         initial_map_matrix_col=NEFr4%collapseMap(mySpinTypeKet,[ka,la]))
      select case(mySpinType)
      case ('singlet')
        do dimi = 1, NEFr4%blockSize(ia)
         do dimj = 1, NEFr4%blockSize(ja)
          idij = NEFr4%blockSize(ja)*(dimi-1) + dimj
          do dimk = 1, NEFr4%blockSize(ka)
           do diml = 1, NEFr4%blockSize(la)
            idlk = NEFr4%blockSize(la)*(dimk-1) + diml
            finalid = nefMatOut%blockIndex([1,idij],[2,idlk])
            nef_bidx = nefr4%blockIndex([ia,dimi],[ja,dimj],[la,diml],[ka,dimk])
            if (NEFr4%Zero_Flag(nef_bidx)) cycle
            call NEFMatOut%Main_Matrix(finalid)%init(&
                   NEFMatOut%Size_vector_alpha_row(idij),&
                   NEFMatOut%Size_vector_beta_row(idij),&
                   nAlpha2=NEFMatOut%Size_vector_alpha_col(idlk),&
                   nBeta2=NEFMatOut%Size_vector_beta_col(idlk))
            NEFMatOut%Zero_Flag(finalid) = .false.
            do i = 1, NEFr4%blockSize(ia,dimi,'alpha')
             do j = 1, NEFr4%blockSize(ja,dimj,'alpha')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,1],[ka,dimk,k,1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'singlet')),&
                 int(NEFr4%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'singlet')),&
                 'alpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,1],[ka,dimk,k,-1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'singlet')),&
                 int(NEFr4%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'singlet')),&
                 'betaalpha')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'beta')
             do j = 1, NEFr4%blockSize(ja,dimj,'beta')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,-1],[ka,dimk,k,1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'singlet')),&
                 int(NEFr4%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'singlet')),&
                 'alphabeta')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,-1],[ka,dimk,k,-1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'singlet')),&
                 int(NEFr4%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'singlet')),&
                 'beta')
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
      case ('triplet')
        do dimi = 1, NEFr4%blockSize(ia)
         do dimj = 1, NEFr4%blockSize(ja)
          idij = nefr4%blockIndex([1,1],[2,1],[3,dimi],[4,dimj])
          do dimk = 1, NEFr4%blockSize(ka)
           do diml = 1, NEFr4%blockSize(la)
            idlk = nefr4%blockIndex([1,1],[2,1],[3,diml],[4,dimk])
            finalid = nefMatOut%blockIndex([1,idij],[2,idlk])
            nef_bidx = nefr4%blockIndex([ia,dimi],[ja,dimj],[la,diml],[ka,dimk])
            if (NEFr4%Zero_Flag(nef_bidx)) cycle
            call NEFMatOut%Main_Matrix(finalid)%init(&
                   NEFMatOut%Size_vector_alpha_row(idij),&
                   NEFMatOut%Size_vector_beta_row(idij),&
                   nAlpha2=NEFMatOut%Size_vector_alpha_col(idlk),&
                   nBeta2=NEFMatOut%Size_vector_beta_col(idlk))
            NEFMatOut%Zero_Flag(finalid) = .false.
            do i = 1, NEFr4%blockSize(ia,dimi,'alpha')
             do j = 1, NEFr4%blockSize(ja,dimj,'beta')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,-1],[ka,dimk,k,1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'triplet')),&
                 int(NEFr4%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'triplet')),&
                 'alpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,-1],[ka,dimk,k,-1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'triplet')),&
                 int(NEFr4%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'triplet')),&
                 'alphabeta')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'beta')
             do j = 1, NEFr4%blockSize(ja,dimj,'alpha')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,1],[ka,dimk,k,1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'triplet')),&
                 int(NEFr4%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'triplet')),&
                 'betaalpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,1],[ka,dimk,k,-1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'triplet')),&
                 int(NEFr4%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'triplet')),&
                 'beta')
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
      case ('all')
        do dimi = 1, NEFr4%blockSize(ia)
         do dimj = 1, NEFr4%blockSize(ja)
          idij = nefr4%blockIndex([1,1],[2,1],[3,dimi],[4,dimj])
          do dimk = 1, NEFr4%blockSize(ka)
           do diml = 1, NEFr4%blockSize(la)
            idlk = nefr4%blockIndex([1,1],[2,1],[3,diml],[4,dimk])
            finalid = nefMatOut%blockIndex([1,idij],[2,idlk])
            nef_bidx = nefr4%blockIndex([ia,dimi],[ja,dimj],[la,diml],[ka,dimk])
            if (NEFr4%Zero_Flag(nef_bidx)) cycle
            call NEFMatOut%Main_Matrix(finalid)%init(&
                   NEFMatOut%Size_vector_alpha_row(idij),&
                   NEFMatOut%Size_vector_beta_row(idij),&
                   nAlpha2=NEFMatOut%Size_vector_alpha_col(idlk),&
                   nBeta2=NEFMatOut%Size_vector_beta_col(idlk))
            NEFMatOut%Zero_Flag(finalid) = .false.
            do i = 1, NEFr4%blockSize(ia,dimi,'alpha')
             do j = 1, NEFr4%blockSize(ja,dimj,'alpha')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,1],[ka,dimk,k,1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'all')),&
                 'alpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,1],[ka,dimk,k,-1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'all')),&
                 'alphabeta')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'beta')
             do j = 1, NEFr4%blockSize(ja,dimj,'beta')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,-1],[ka,dimk,k,1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'all')),&
                 'betaalpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,-1],[ka,dimk,k,-1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'all')),&
                 'beta')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'alpha')
             do j = 1, NEFr4%blockSize(ja,dimj,'beta')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,-1],[ka,dimk,k,1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'all')),&
                 'alpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,1],[ja,dimj,j,-1],[ka,dimk,k,-1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'all')),&
                 'alphabeta')
               end do
              end do
             end do
            end do
            do i = 1, NEFr4%blockSize(ia,dimi,'beta')
             do j = 1, NEFr4%blockSize(ja,dimj,'alpha')
              do k = 1, NEFr4%blockSize(ka,dimk,'alpha')
               do l = 1, NEFr4%blockSize(la,diml,'beta')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,1],[ka,dimk,k,1],[la,diml,l,-1]),&
                 int(NEFr4%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'all')),&
                 'betaalpha')
               end do
              end do
              do k = 1, NEFr4%blockSize(ka,dimk,'beta')
               do l = 1, NEFr4%blockSize(la,diml,'alpha')
                call NEFMatOut%Main_Matrix(finalid)%put(&
                 NEFR4%at(&
                     [ia,dimi,i,-1],[ja,dimj,j,1],[ka,dimk,k,-1],[la,diml,l,1]),&
                 int(NEFr4%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'all')),&
                 int(NEFr4%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'all')),&
                 'beta')
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
      case default
        call mqc_error_a('mqc_nef_r4_dimention_collapse: Unrecognized excitation type',6,&
          'spinType',spinType)
      end select  
!
      end function mqc_nef_r4_dimention_collapse
!
!
!     PROCEDURE MQC_NEF_R4_Block_Diagonal_PInv
!
!>    \brief <b> MQC_NEF_R4_Block_Diagonal_PInv is used to return the block-diagonal
!>    pseudo-inverse of a NEF R4 object, expressed as a collapsed NEF Matrix. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Block_Diagonal_PInv collapses the input NEF R4 object (using the
!>    same collapse convention as MQC_NEF_R4_Dimension_Collapse) and computes the
!>    pseudo-inverse of each of its diagonal blocks -- the blocks whose collapsed
!>    row space-index equals its collapsed column space-index. Off-diagonal blocks
!>    of the result are left zero-flagged, so the returned object is block-diagonal
!>    by construction, regardless of any off-diagonal content in the collapsed
!>    input. Zero diagonal blocks are left zero rather than pseudo-inverted.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFr4 is type(mqc_nef_r4)
!>        NEF R4 Object whose block-diagonal pseudo-inverse is to be returned.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        spinType is Character(Len=*),Optional
!>        spinType modiffies the excitations considered, applied identically to
!>        bra and ket (the collapsed row and column maps always match)
!>          Options:
!>            'singlet' => return same spin
!>            'triplet  => return opposite spin
!>            'all'     => return both types (default)
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Order is integer(kind=int64),dimension(4),intent(in),optional
!>        Order contains the indexes to be collapsed following:
!>        Order = [i,j,k,l] => i->j = p and  k->l = q
!>        (default [1,2,3,4])
!>    \endverbatim
!>
!>    \param[in] thresh
!>    \verbatim
!>        thresh is real(kind=real64),optional
!>        Singular-value threshold passed to each diagonal block's pseudo-inverse
!>        (default 1.0d-8)
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_block_diagonal_pinv(NEFr4,spinType,Order,thresh) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFr4
      character(len=*),intent(in),optional::spinType
      integer(kind=int64),dimension(4),intent(in),optional::Order
      real(kind=real64),intent(in),optional::thresh
      type(mqc_nef_matrix)::NEFMatOut
      type(mqc_nef_matrix)::CollapsedMat
      real(kind=real64)::myThresh
      integer(kind=int64)::N,id
!
      if(present(thresh)) then
        myThresh = thresh
      else
        myThresh = 1.0d-8
      end if
!
      CollapsedMat = NEFr4%collapse(spinType=spinType,Order=Order)
!
      call NEFMatOut%init(CollapsedMat%Map_Matrix_Row,&
           initial_map_matrix_col=CollapsedMat%Map_Matrix_Col)
!
      do N = 1, CollapsedMat%NumSpace_row
        id = CollapsedMat%blockIndex([1,N],[2,N])
        if (.not.CollapsedMat%Zero_Flag(id)) then
          NEFMatOut%Main_Matrix(id) = CollapsedMat%Main_Matrix(id)%pinv(myThresh)
          NEFMatOut%Zero_Flag(id) = .false.
        end if
      end do
!
      end function mqc_nef_r4_block_diagonal_pinv
!
!
!     PROCEDURE MQC_NEF_Matrix_Collapse_Index
!
!
!>    \brief <b> MQC_NEF_Matrix_Collapse_Index is used to return a collapsed
!>    two index excitation. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Collapse_Index is used to return a collapsed two index
!>    two index excitation (i,j->r).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFM
!>    \verbatim
!>        NEFM is type(mqc_nef_matrix) 
!>        NEF Matrix Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] I  
!>    \verbatim
!>        I is Integer(kind=int64)
!>        The first dimension of the element in MQC NEF Matrix
!>          If I>0 row count is from first index
!>    \endverbatim
!>    
!>    \param[in] J  
!>    \verbatim
!>        J is Integer(kind=int64)
!>        The second dimension of the element in MQC NEF Matrix
!>          If J>0 col count is from first index
!>    \endverbatim
!>    
!>    \param[in] i_spin 
!>    \verbatim
!>        i_spin is character(len=1)
!>        Spin of the annihilation operator natorb.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>        Options:
!>        a                 =>    alpha spin
!>        b                 =>    beta spin
!>    \endverbatim
!>
!>    \param[in] j_spin
!>    \verbatim
!>        j_spin is character(len=1)
!>        Spin of the creation operator natorb.
!>        If a_index or a_index is negative, this option does not 
!>        change the result.
!>        Options:
!>        a                 =>    alpha spin
!>        b                 =>    beta spin
!>    \endverbatim
!>
!>    \param[in] i_space
!>    \verbatim
!>        i_space is an integer
!>        Index of the annihilation operator natorb subspace.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>    \endverbatim
!>
!>    \param[in] j_space
!>    \verbatim
!>        j_space is an integer
!>        Index of the creation operator natorb subspace.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>    \endverbatim
!>
!>    \param[in] excType
!>    \verbatim
!>        excType is Character(Len=*),Optional
!>        excType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!>
!>    \param[in] swapOrder
!>    \verbatim
!>        swapOrder is logical,Optional
!>        swapOrder defines the excitation order
!>          Options:
!>            .false.   => row -> col (default), row-major index
!>            .true.    => col -> row, column-major index
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_collapse_index(NEFM,i,j,i_spin,j_spin,i_space,j_space,excType,&
          swapOrder) result(idOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFM
      integer::idOut
      integer(kind=int64)::i,j,i_space,j_space
      character(len=1)::i_spin,j_spin
      character(len=*),intent(in),optional::excType
      logical,intent(in),optional::swapOrder
      character(len=64)::selectedExcType
      logical::mySwapOrder
!
    if(present(swapOrder)) then
      mySwapOrder = swapOrder
    else
      mySwapOrder = .false.
    end if
!
    if(i_spin.ne.'a'.and.i_spin.ne.'b') &
      call mqc_Error_a('mqc_nef_matrix_collapse_index: Unrecognized spin',&
           6,'i_spin',i_spin)
    if(j_spin.ne.'a'.and.j_spin.ne.'b') &
      call mqc_Error_a('mqc_nef_matrix_collapse_index: Unrecognized spin',&
            6,'j_spin',j_spin)
    if(i_space.gt.nefM%NumSpace_row.or.i_space.lt.1) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: Space index out-of-bound',&
            6,'i_space',i_space)
    if(j_space.gt.nefM%NumSpace_col.or.j_space.lt.1)&
      call mqc_Error_i('mqc_nef_matrix_collapse_index: Space index out-of-bound',&
            6,'j_space',j_space)
!
    idOut = 0
    if ((i_spin.eq.'b'.and.nefM%Size_vector_beta_row(i_space).eq.0).or.&
        (i_spin.eq.'a'.and.nefM%Size_vector_alpha_row(i_space).eq.0).or.&
        (j_spin.eq.'b'.and.nefM%Size_vector_beta_col(j_space).eq.0).or.& 
        (j_spin.eq.'a'.and.nefM%Size_vector_alpha_col(j_space).eq.0)) return
!
    if (i_spin.eq.'b'.and.i.gt.nefM%Size_vector_beta_row(i_space).or.i.lt.0) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: i beta index out-of-bound',&
            6,'i',i,'nefM%Size_vector_beta_row(i_space)',nefM%Size_vector_beta_row(i_space))
    if (i_spin.eq.'a'.and.i.gt.nefM%Size_vector_alpha_row(i_space).or.i.lt.0) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: i alpha index out-of-bound',&
            6,'i',i,'nefM%Size_vector_alpha_row(i_space)',nefM%Size_vector_alpha_row(i_space))
!
    if (j_spin.eq.'b'.and.j.gt.nefM%Size_vector_beta_col(j_space).or.j.lt.0) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: j beta index out-of-bound',&
            6,'j',j,'nefM%Size_vector_beta_col(j_space)',nefM%Size_vector_beta_col(j_space))
    if (j_spin.eq.'a'.and.j.gt.nefM%Size_vector_alpha_col(j_space).or.j.lt.0) &
      call mqc_Error_i('mqc_nef_matrix_collapse_index: j alpha index out-of-bound',&
            6,'j',j,'nefM%Size_vector_alpha_col(j_space)',nefM%Size_vector_alpha_col(j_space))
!
    if(present(excType)) then
      selectedExcType = trim(excType) !singlet, triplet, all
    else
      selectedExcType = 'singlet'
    end if
    if (selectedExcType.eq.'singlet'.and.i_spin.ne.j_spin) &
      call mqc_Error('mqc_nef_matrix_collapse_index: wrong spin for triplet',6)
    if (selectedExcType.eq.'triplet'.and.i_spin.eq.j_spin) &
      call mqc_Error('mqc_nef_matrix_collapse_index: wrong spin for triplet',6)
!
    idOut = 0
    select case(selectedExcType)
    case ('singlet')
      if (mySwapOrder) then
        if (i_spin.eq.'b') then
          idOut = (j-1)*nefM%Size_vector_beta_row(i_space) + i
        else
          idOut = (j-1)*nefM%Size_vector_alpha_row(i_space) + i
        end if
      else
        if (j_spin.eq.'b') then
          idOut = (i-1)*nefM%Size_vector_beta_col(j_space) + j
        else
          idOut = (i-1)*nefM%Size_vector_alpha_col(j_space) + j
        end if
      end if
    case ('triplet')
      if (mySwapOrder) then
        if (i_spin.eq.'b') then
          idOut = (j-1)*nefM%Size_vector_alpha_row(i_space) + i
        else
          idOut = (j-1)*nefM%Size_vector_beta_row(i_space) + i
        end if
      else
        if (j_spin.eq.'b') then
          idOut = (i-1)*nefM%Size_vector_alpha_col(j_space) + j
        else
          idOut = (i-1)*nefM%Size_vector_beta_col(j_space) + j
        end if
      end if
    case ('all','full')
      if (mySwapOrder) then
        if (i_spin.eq.'b') then
          idOut = (j-1)*nefM%Size_vector_beta_row(i_space) + i
        else
          idOut = (j-1)*(nefM%Size_vector_beta_row(i_space)+nefM%Size_vector_alpha_row(i_space)) + &
                i + nefM%Size_vector_beta_row(i_space)
        end if
      else
        if (j_spin.eq.'b') then
          idOut = (i-1)*nefM%Size_vector_beta_col(j_space) + j
        else
          idOut = (i-1)*(nefM%Size_vector_beta_col(j_space)+nefM%Size_vector_alpha_col(j_space)) + &
                j + nefM%Size_vector_beta_col(j_space)
        end if
      end if
    case default
      call mqc_error_a('mqc_nef_matrix_collapse_index: Unrecognized',6,'selectedExcType',&
        selectedExcType)
    end select
!
      end function mqc_nef_matrix_collapse_index
!
!
!     PROCEDURE MQC_NEF_R4_Collapse_Index
!
!
!>    \brief <b> MQC_NEF_R4_Collapse_Index is used to return a collapsed
!>    two index excitation. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Collapse_Index is used to return a collapsed two index
!>    two index excitation.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is type(mqc_nef_r4) 
!>        NEF R4 Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] I  
!>    \verbatim
!>        I is Integer(kind=int64)
!>        The dimi of the element in MQC NEF R4
!>          If I>0 row count is from first index
!>    \endverbatim
!>    
!>    \param[in] J  
!>    \verbatim
!>        J is Integer(kind=int64)
!>        The dimj of the element in MQC NEF MR4
!>          If J>0 col count is from first index
!>    \endverbatim
!>    
!>    \param[in] i_spin 
!>    \verbatim
!>        i_spin is character(len=1)
!>        Spin of the annihilation operator natorb.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>        Options:
!>        a                 =>    alpha spin
!>        b                 =>    beta spin
!>    \endverbatim
!>
!>    \param[in] j_spin
!>    \verbatim
!>        j_spin is character(len=1)
!>        Spin of the creation operator natorb.
!>        If a_index or a_index is negative, this option does not 
!>        change the result.
!>        Options:
!>        a                 =>    alpha spin
!>        b                 =>    beta spin
!>    \endverbatim
!>
!>    \param[in] i_space
!>    \verbatim
!>        i_space is an integer
!>        Index of the annihilation operator natorb subspace.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>    \endverbatim
!>
!>    \param[in] j_space
!>    \verbatim
!>        j_space is an integer
!>        Index of the creation operator natorb subspace.
!>        If i_index or a_index is negative, this option does not 
!>        change the result.
!>    \endverbatim
!>
!>    \param[in] dimi
!>    \verbatim
!>        dimi is an integer
!>        Index of the creation operator natorb subspace.
!>    \endverbatim
!>
!>    \param[in] dimj
!>    \verbatim
!>        dimj is an integer
!>        Index of the creation operator natorb subspace.
!>    \endverbatim
!>
!>    \param[in] excType
!>    \verbatim
!>        excType is Character(Len=*),Optional 
!>        excType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin 
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_collapse_index(NEFM,i,j,i_spin,j_spin,&
                                          i_space,j_space,dimi,dimj,excType) result(idOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFM
      integer::idOut
      integer(kind=int64)::i,j,i_space,j_space,dimi,dimj
      character(len=1)::i_spin,j_spin
      character(len=*),intent(in),optional::excType
      character(len=64)::selectedExcType
!
    if(dimi.le.0.or.dimi.gt.4) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: dimention index out-of-bound',&
           6,'dimi',dimi)
    if(dimj.le.0.or.dimj.gt.4) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: dimention index out-of-bound',&
           6,'dimj',dimj)
!
    if(i_spin.ne.'a'.and.i_spin.ne.'b') &
      call mqc_Error_a('mqc_nef_r4_collapse_index: Unrecognized spin',&
           6,'i_spin',i_spin)
    if(j_spin.ne.'a'.and.j_spin.ne.'b') &
      call mqc_Error_a('mqc_nef_r4_collapse_index: Unrecognized spin',&
            6,'j_spin',j_spin)
    if(i_space.gt.nefM%blockSize(dimi).or.i_space.lt.1) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: Space index out-of-bound',&
            6,'i_space',i_space)
    if(j_space.gt.nefM%blockSize(dimj).or.j_space.lt.1)&
      call mqc_Error_i('mqc_nef_r4_collapse_index: Space index out-of-bound',&
            6,'j_space',j_space)
!
    idOut = 0
    if ((i_spin.eq.'b'.and.nefM%blockSize(dimi,i_space,'beta').eq.0).or.&
        (i_spin.eq.'a'.and.nefM%blockSize(dimi,i_space,'alpha').eq.0).or.&
        (j_spin.eq.'b'.and.nefM%blockSize(dimj,j_space,'beta').eq.0).or.& 
        (j_spin.eq.'a'.and.nefM%blockSize(dimj,j_space,'alpha').eq.0)) return
!
    if (i_spin.eq.'b'.and.i.gt.nefM%blockSize(dimi,i_space,'beta').or.i.lt.0) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: i beta index out-of-bound',&
            6,'i',i,'nefM%blockSize(dimi,i_space,beta)',nefM%blockSize(dimi,i_space,'beta'))
    if (i_spin.eq.'a'.and.i.gt.nefM%blockSize(dimi,i_space,'alpha').or.i.lt.0) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: i alpha index out-of-bound',&
            6,'i',i,'nefM%blockSize(dimi,i_space,alpha)',nefM%blockSize(dimi,i_space,'alpha'))
!
    if (j_spin.eq.'b'.and.j.gt.nefM%blockSize(dimj,j_space,'beta').or.j.lt.0) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: j beta index out-of-bound',&
            6,'j',j,'nefM%blockSize(dimj,j_space,beta)',nefM%blockSize(dimj,j_space,'beta'))
    if (j_spin.eq.'a'.and.j.gt.nefM%blockSize(dimj,j_space,'alpha').or.j.lt.0) &
      call mqc_Error_i('mqc_nef_r4_collapse_index: j alpha index out-of-bound',&
            6,'j',j,'nefM%blockSize(dimj,j_space,alpha)',nefM%blockSize(dimj,j_space,'alpha'))
!
    if(present(excType)) then
      selectedExcType = trim(excType) !singlet, triplet, all
    else
      selectedExcType = 'singlet'
    end if
    if (selectedExcType.eq.'singlet'.and.i_spin.ne.j_spin) &
      call mqc_Error('mqc_nef_r4_collapse_index: wrong spin for triplet',6)
    if (selectedExcType.eq.'triplet'.and.i_spin.eq.j_spin) &
      call mqc_Error('mqc_nef_r4_collapse_index: wrong spin for triplet',6)
!
    idOut = 0
    select case(selectedExcType)
    case ('singlet')
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'beta') + j
      else
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'alpha') + j
      end if 
    case ('triplet')
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'alpha') + j
      else
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'beta') + j
      end if 
    case ('all','full')
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'beta') + j
      else
        idOut = (i-1)*(nefM%blockSize(dimj,j_space,'beta') + &
                       nefM%blockSize(dimj,j_space,'alpha')) + &
                   j + nefM%blockSize(dimj,j_space,'beta')
      end if 
    case default
      call mqc_error_a('mqc_nef_r4_collapse_index: Unrecognized',6,'selectedExcType',&
        selectedExcType)
    end select
!
      end function mqc_nef_r4_collapse_index
!
!
!     PROCEDURE MQC_NEF_Matrix_Collapse_Map
!
!
!>    \brief <b> MQC_NEF_Matrix_Collapse_Map is used to return a collapsed
!>    two Map_Vectors. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Collapse_Map is used to return a collapsed two Map_Vectors.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFM
!>    \verbatim
!>        NEFM is type(mqc_nef_matrix) 
!>        NEF Matrix Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] excType
!>    \verbatim
!>        excType is Character(Len=*),Optional 
!>        excType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!>
!>    \param[in] swapOrder
!>    \verbatim
!>        swapOrder is logical,Optional
!>        swapOrder defines the excitation order, and must match whatever
!>        value is later passed to dimention_collapse/collapseIndex for this
!>        same NEFM, so that each block's pre-computed size lands in the same
!>        slot its data is later written to.
!>          Options:
!>            .false.   => row -> col (default), row-major index
!>            .true.    => col -> row, column-major index
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_collapse_map(NEFM,excType,swapOrder) result(MapOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFM
      character(len=*),intent(in),optional::excType
      logical,intent(in),optional::swapOrder
      character(len=64)::selectedExcType
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::MapOut
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map_Col,New_Map_Row
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,id,outId,idOld,r0
      logical::full_row,full_col
      logical::mySwapOrder
!
      if(present(excType)) then
        selectedExcType = trim(excType) !singlet, triplet, all
      else
        selectedExcType = 'singlet'
      end if
!
      if(present(swapOrder)) then
        mySwapOrder = swapOrder
      else
        mySwapOrder = .false.
      end if
!
      if(allocated(MapOut)) deallocate(MapOut)
      allocate(MapOut(NEFM%NumSpace_row*NEFM%NumSpace_col))
!
      select case(selectedExcType)
      case ('singlet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            if (mySwapOrder) then
              outId = (col-1)*NEFM%NumSpace_row + row
            else
              outId = id
            end if
            call MapOut(outId)%init(&
               nefM%Size_vector_alpha_col(col)*nefM%Size_vector_alpha_row(row),0,&
               nefM%Size_vector_beta_col(col)*nefM%Size_vector_beta_row(row))
          end do
        end do
      case ('triplet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            if (mySwapOrder) then
              outId = (col-1)*NEFM%NumSpace_row + row
            else
              outId = id
            end if
            call MapOut(outId)%init(&
               nefM%Size_vector_alpha_col(col)*nefM%Size_vector_beta_row(row),0,&
               nefM%Size_vector_beta_col(col)*nefM%Size_vector_alpha_row(row))
          end do
        end do
      case ('all','full')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            id = nefM%blockIndex([1,row],[2,col])
            if (mySwapOrder) then
              outId = (col-1)*NEFM%NumSpace_row + row
            else
              outId = id
            end if
            call MapOut(outId)%init(&
               nefM%Size_vector_alpha_col(col)*nefM%Size_vector_alpha_row(row)+&
               nefM%Size_vector_alpha_col(col)*nefM%Size_vector_beta_row(row),0,&
               nefM%Size_vector_beta_col(col)*nefM%Size_vector_beta_row(row)+&
               nefM%Size_vector_beta_col(col)*nefM%Size_vector_alpha_row(row))
          end do
        end do
      case default
        call mqc_error_a('mqc_nef_matrix_collapse_map: Unrecognized excitation type',6,&
          'selectedExcType',selectedExcType)
      end select  
!
      end function mqc_nef_matrix_collapse_map
!
!
!     PROCEDURE MQC_NEF_R4_Collapse_Map
!
!
!>    \brief <b> MQC_NEF_R4_Collapse_Map is used to return a collapsed
!>    two Map_Vectors. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Collapse_Map is used to return a collapsed two Map_Vectors.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFr4
!>    \verbatim
!>        NEFr4 is type(mqc_nef_r4) 
!>        NEF R4 Object to be collapsed.
!>    \endverbatim
!>
!>    \param[in] excType
!>    \verbatim
!>        excType is Character(Len=*),Optional 
!>        excType modiffies the excitations considered
!>          Options:
!>            'singlet'    => return same spin 
!>            'triplet     => return opposite spin
!>            'all','full' => return both types (default)
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Order is integer(kind=int64),dimension(2),intent(in),optional
!>        Order contains the indexes to be collapsed following:
!>        Order = [i,j] => i->j = p, where i and j are the dimentions
!>        (default [1,2])
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_collapse_map(NEFr4,excType,Order) result(MapOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFr4
      character(len=*),intent(in),optional::excType
      character(len=64)::selectedExcType
      integer(kind=int64),dimension(2),intent(in),optional::Order
      integer(kind=int64),dimension(2)::my_order
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::MapOut
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::New_Map_Col,New_Map_Row
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::ia,ja,ka,la,row,col,id
      logical::full_row,full_col
!
      if(present(excType)) then
        selectedExcType = trim(excType) !singlet, triplet, all
      else
        selectedExcType = 'singlet'
      end if
!
      if(present(Order)) then
        ia=Order(1)
        ja=Order(2)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in MQC_nef_matrix_dimention_collapse',&
                                        6,'ia', ia,'ja',ja)
        my_order = Order
      else
        ia=1
        ja=2
        my_order = [1,2]
      endif
!
      if(allocated(MapOut)) deallocate(MapOut)
      allocate(MapOut(int(NEFr4%blockSize(ia)*NEFr4%blockSize(ja))))
!
      select case(selectedExcType)
      case ('singlet')
        do row = 1, NEFr4%blockSize(ia)
          do col = 1, NEFr4%blockSize(ja)
            id = nefR4%blockIndex([1,1],[2,1],[3,row],[4,col])
            call MapOut(id)%init(&
               nefR4%blockSize(ia,row,'alpha')*nefR4%blockSize(ja,col,'alpha'),0,&
               nefR4%blockSize(ia,row,'beta')*nefR4%blockSize(ja,col,'beta'))
          end do
        end do
      case ('triplet')
        do row = 1, NEFr4%blockSize(ia)
          do col = 1, NEFr4%blockSize(ja)
            id = nefr4%blockIndex([1,1],[2,1],[3,row],[4,col])
            call MapOut(id)%init(&
               nefR4%blockSize(ia,row,'alpha')*nefR4%blockSize(ja,col,'beta'),0,&
               nefR4%blockSize(ia,row,'beta')*nefR4%blockSize(ja,col,'alpha'))
          end do
        end do
      case ('all','full')
        do row = 1, NEFr4%blockSize(ia)
          do col = 1, NEFr4%blockSize(ja)
            id = nefr4%blockIndex([1,1],[2,1],[3,row],[4,col])
            call MapOut(id)%init(&
               nefR4%blockSize(ia,row,'alpha')*(nefR4%blockSize(ja,col,'beta')+&
                                                nefR4%blockSize(ja,col,'alpha')),0,&
               nefR4%blockSize(ia,row,'beta')*(nefR4%blockSize(ja,col,'alpha')+&
                                               nefR4%blockSize(ja,col,'beta')))
          end do
        end do
      case default
        call mqc_error_a('mqc_nef_r4_collapse_map: Unrecognized excitation type',6,&
          'selectedExcType',selectedExcType)
      end select  
!
      end function mqc_nef_r4_collapse_map
!
!
!     PROCEDURE mqc_nef_vector_output_block_vector
!
!>    \brief <b> mqc_nef_vector_output_block_vector is used to output a MQC Vector of an MQC 
!>    NEF Vector type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    mqc_nef_vector_output_block_vector is used to output a MQC Vector of an MQC 
!>    NEF Vector type variable
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] nefV
!>    \verbatim
!>        nefV is Class(MQC_NEF_Vector)
!>        The MQC_NEF_Vecotr variable to output.
!>    \endverbatim
!>    
!>    \param[in] indexes 
!>    \verbatim
!>        indexes is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the blocks to be returned.
!>        Default: all indexes
!>    \endverbatim
!>    
!>    \param[in] BlockName
!>    \verbatim
!>        BlockName is character(len=*),optional
!>        = 'full':        outputs the full MQC_NEF_Vector.(default) 
!>        = 'alpha','a':   outputs the alpha spin block.
!>        = 'beta','b':    outputs the beta spin block. 
!>    \endverbatim
!>    
!>    \param[in] BlockOrder
!>    \verbatim
!>        BlockOrder is character(len=*),optional
!>        = 'beta': orders the output as b->b, b->a, a->b, a->a (default)
!>        = 'alpha':orders the output as a->b, a->a, b->b, b->a
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_output_block_vector(nefV,indexes,blockName,BlockOrder) result(vectorOut)
!
      implicit none
      class(mqc_nef_vector)::nefV
      integer(kind=int64),dimension(:),optional::indexes
      character(len=*),optional::blockOrder
      character(len=*),optional::blockName
      type(mqc_vector)::vectorOut
      character(len=64)::myBlockName,myBlockOrder
      integer(kind=int64),dimension(:),allocatable::myId
      logical::full
      integer(kind=int64)::fsize,i,j
!
      if(present(blockOrder)) then
        call string_change_case(blockOrder,'l',myBlockOrder)
      else
        myBlockOrder = 'beta'
      endIf
      if(myBlockOrder.ne.'alpha'.and.myBlockOrder.ne.'beta') &
        call mqc_error_A('BlockOrder inrecognized in mqc_nef_vector_output_block_vector',6,&
           'myBlockOrder',myBlockOrder)
      if(present(blockName)) then
        call string_change_case(blockName,'l',myBlockName)
      else
        myBlockName = 'full'
      endIf
      if(myBlockName.ne.'alpha'.and.myBlockName.ne.'beta'.and.myBlockName.ne.'full') &
        call mqc_error_A('BlockName inrecognized in mqc_nef_vector_output_block_vector',6,&
           'myBlockName',myBlockName)
      if(present(indexes)) then
        myid = indexes
        full=.false.
      else
        full=.true.
      end if
!
      fsize = 0
      if (full) then
        do i = 1, nefV%NumSpace
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsize = fsize + nefV%Size_vector_alpha(i)
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsize = fsize + nefV%Size_vector_beta(i)
        end do
      else
        do i = 1, size(myid)
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsize = fsize + nefV%Size_vector_alpha(myid(i))
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsize = fsize + nefV%Size_vector_beta(myid(i))
        end do
      end if
      call vectorOut%init(fsize)
      fsize = 1
      if (full) then
        do i = 1, nefV%NumSpace
          if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then
            if (.not.nefV%Zero_Flag(i).and.(nefV%Size_vector_alpha(i)+nefV%Size_vector_beta(i)).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('full'),fsize)
            fsize = fsize + (nefV%Size_vector_alpha(i)+nefV%Size_vector_beta(i))
          else if (myBlockName.eq.'full') then
            if (.not.nefV%Zero_Flag(i).and.nefV%Size_vector_beta(i).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('beta'),fsize)
            fsize = fsize + nefV%Size_vector_beta(i)
            if (.not.nefV%Zero_Flag(i).and.nefV%Size_vector_alpha(i).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('alpha'),fsize)
            fsize = fsize + nefV%Size_vector_alpha(i)
          else if (myBlockName.eq.'alpha') then
            if (.not.nefV%Zero_Flag(i).and.nefV%Size_vector_alpha(i).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('alpha'),fsize)
            fsize = fsize + nefV%Size_vector_alpha(i)
          else
            if (.not.nefV%Zero_Flag(i).and.nefV%Size_vector_beta(i).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('beta'),fsize)
            fsize = fsize + nefV%Size_vector_beta(i)
          end if
        end do
      else
        do j = 1, size(myid)
          i = myid(j)
          if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then
            if (.not.nefV%Zero_Flag(i).and.(nefV%Size_vector_alpha(i)+nefV%Size_vector_beta(i)).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('full'),fsize)
            fsize = fsize + (nefV%Size_vector_alpha(i)+nefV%Size_vector_beta(i))
          else if (myBlockName.eq.'full') then
            if (.not.nefV%Zero_Flag(i).and.nefV%Size_vector_beta(i).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('beta'),fsize)
            fsize = fsize + nefV%Size_vector_beta(i)
            if (.not.nefV%Zero_Flag(i).and.nefV%Size_vector_alpha(i).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('alpha'),fsize)
            fsize = fsize + nefV%Size_vector_alpha(i)
          else if (myBlockName.eq.'alpha') then
            if (.not.nefV%Zero_Flag(i).and.nefV%Size_vector_alpha(i).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('alpha'),fsize)
            fsize = fsize + nefV%Size_vector_alpha(i)
          else
            if (.not.nefV%Zero_Flag(i).and.nefV%Size_vector_beta(i).gt.0) &
              call vectorOut%vput(nefV%Main_vector(i)%getBlock('beta'),fsize)
            fsize = fsize + nefV%Size_vector_beta(i)
          end if
        end do
      end if
!
      end function mqc_nef_vector_output_block_vector
!
!
!     PROCEDURE mqc_nef_matrix_output_block_matrix
!
!>    \brief <b> mqc_nef_matrix_output_block_matrix is used to output a MQC Matrix of an MQC 
!>    NEF Matrix type variable</b> 
!
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    mqc_nef_matrix_output_block_matrix is used to output a MQC Matrix of an MQC 
!>    NEF Matrix type variable
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] nefMat
!>    \verbatim
!>        nefMat is Class(MQC_NEF_Matrix)
!>        The MQC_NEF_Matrix variable to output.
!>    \endverbatim
!>    
!>    \param[in] indexes 
!>    \verbatim
!>        indexes is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the blocks to be returned.
!>        Default: all indexes
!>    \endverbatim
!>    
!>    \param[in] indexesKet
!>    \verbatim
!>        indexesKet is integer(kind=int64),dimension(:),optional
!>        Contains the indexes of the blocks to be returned.
!>        Default: all indexes
!>    \endverbatim
!>    
!>    \param[in] BlockName
!>    \verbatim
!>        BlockName is character(len=*),optional
!>        = 'full'          :   outputs the full MQC_NEF_Vector.(default) 
!>        = 'alpha','a','aa':   outputs the alpha spin block.
!>        = 'beta', 'b','bb':   outputs the beta spin block. 
!>        = 'alphabeta','ab':   outputs the alphabeta spin block.
!>        = 'betaalpha','ba':   outputs the betaalpha spin block. 
!>    \endverbatim
!>    
!>    \param[in] BlockOrder
!>    \verbatim
!>        BlockOrder is character(len=*),optional
!>        = 'beta': orders the output as b->b, b->a, a->b, a->a (default)
!>        = 'alpha':orders the output as a->b, a->a, b->b, b->a
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_output_block_matrix(nefMat,indexes,indexesKet,blockName,&
                                                  BlockOrder) result(MatrixOut)
!
      implicit none
      class(mqc_nef_matrix)::nefMat
      integer(kind=int64),dimension(:),optional::indexes,indexesKet
      character(len=*),optional::blockOrder
      character(len=*),optional::blockName
      type(mqc_Matrix)::MatrixOut
      type(mqc_matrix)::zeroBlock
      character(len=64)::myBlockName,myBlockOrder
      integer(kind=int64),dimension(:),allocatable::myId,myIdKet
      logical::full,fullKet
      integer(kind=int64)::fsize,fsizeKet,i,j,k,p,size1,size2,size3,size4,id
!
      if(present(blockOrder)) then
        call string_change_case(blockOrder,'l',myBlockOrder)
      else
        myBlockOrder = 'beta'
      endIf
      if(myBlockOrder.ne.'alpha'.and.myBlockOrder.ne.'beta') &
        call mqc_error_A('BlockOrder inrecognized in mqc_nef_matrix_output_block_matrix',6,&
           'myBlockOrder',myBlockOrder)
      if(present(blockName)) then
        call string_change_case(blockName,'l',myBlockName)
      else
        myBlockName = 'full'
      endIf
      if(myBlockName.ne.'alpha'.and.myBlockName.ne.'beta'.and.myBlockName.ne.'full') &
        call mqc_error_A('BlockName inrecognized in mqc_nef_matrix_output_block_matrix',6,&
           'myBlockName',myBlockName)
      if(present(indexes)) then
        myid = indexes
        full = .false.
      else
        full = .true.
      end if
      if(present(indexesKet)) then
        myidKet = indexesKet
        fullKet = .false.
      else
        fullKet = .true.
      end if
!
      fsizeKet = 0
      if (full) then
        do i = 1, nefMat%NumSpace_col
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsizeKet = fsizeKet + nefMat%Size_vector_alpha_col(i)
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsizeKet = fsizeKet + nefMat%Size_vector_beta_col(i)
        end do
      else
        do i = 1, size(myidKet)
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsizeKet = fsizeKet + nefMat%Size_vector_alpha_col(myidKet(i))
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsizeKet = fsizeKet + nefMat%Size_vector_beta_col(myidKet(i))
        end do
      end if
      fsize = 0
      if (fullKet) then
        do i = 1, nefMat%NumSpace_row
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsize = fsize + nefMat%Size_vector_alpha_row(i)
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsize = fsize + nefMat%Size_vector_beta_row(i)
        end do
      else
        do i = 1, size(myid)
          if(myBlockName.eq.'alpha'.or.myBlockName.eq.'full') &
            fsize = fsize + nefMat%Size_vector_alpha_row(myid(i))
          if(myBlockName.eq.'beta'.or.myBlockName.eq.'full') &
            fsize = fsize + nefMat%Size_vector_beta_row(myid(i))
        end do
      end if
      call matrixOut%init(fsize,fsizeKet)
      if (full) then
        fsize = 1
        do i = 1, nefMat%NumSpace_row
          if (fullKet) then
            fsizeKet = 1
            do j = 1, nefMat%NumSpace_col
              id = nefMat%blockIndex([1,i],[2,j])
              if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then
                size1 = nefMat%Size_vector_alpha_row(i)+nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)+nefMat%Size_vector_beta_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                           [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('full'),&
                                           [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'full') then
                size1 = nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_beta_col(j)
                size3 = nefMat%Size_vector_alpha_row(i)
                size4 = nefMat%Size_vector_alpha_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                       [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                       [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                if (size3.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size3,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                                 [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                if (size1.gt.0.and.size4.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size4,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                                 [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  end if
                end if
                if (size3.gt.0.and.size4.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size3,size4,0.0)
                    call MatrixOut%mput(zeroBlock,&
                     [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                     [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2 + size4
                size1 = size1 + size3
              else if (myBlockName.eq.'alpha') then
                size1 = nefMat%Size_vector_alpha_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'beta') then
                size1 = nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_beta_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'alphabeta') then
                size1 = nefMat%Size_vector_alpha_row(i)
                size2 = nefMat%Size_vector_beta_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'betaalpha') then
                size1 = nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('full'),&
                                 [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'full') then
                size1 = nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_beta_col(j)
                size3 = nefMat%Size_vector_alpha_row(i)
                size4 = nefMat%Size_vector_alpha_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                       [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('beta'),&
                                       [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                if (size3.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size3,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                                 [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                if (size1.gt.0.and.size4.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size4,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alphabeta'),&
                                 [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  end if
                end if
                if (size3.gt.0.and.size4.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size3,size4,0.0)
                    call MatrixOut%mput(zeroBlock,&
                     [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                     [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2 + size4
                size1 = size1 + size3
              else if (myBlockName.eq.'alpha') then
                size1 = nefMat%Size_vector_alpha_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (nefMat%Zero_Flag(id)) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('alpha'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'beta') then
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'betaalpha') then
                size1 = nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)
                if (nefMat%Zero_Flag(id)) then
                  call zeroBlock%init(size1,size2,0.0)
                  call MatrixOut%mput(zeroBlock,&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                else
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                end if
                fsizeKet = fsizeKet + size2
                if (nefMat%Zero_Flag(id)) then
                  call zeroBlock%init(size1,size2,0.0)
                  call MatrixOut%mput(zeroBlock,&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                else
                  call MatrixOut%mput(nefMat%Main_matrix(id)%getBlock('betaalpha'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                end if
                fsizeKet = fsizeKet + size2
              else
                call mqc_error_a('block name and order pair unrecognized in&
                                 & mqc_nef_matrix_output_block_matrix',6,&
                     'myBlockName',myBlockName,'myBlockOrder',myBlockOrder)
              end if
            end do
          end if
          fsize = fsize + size1
        end do
      end if
!
      end function mqc_nef_matrix_output_block_matrix
!
!
!     PROCEDURE MQC_NEF_Matrix_Get_Map
!
!
!>    \brief <b> MQC_NEF_Matrix_Get_Map returns one of a MQC_NEF_Matrix's
!>    internal block maps. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Get_Map returns a copy of either Map_Matrix_Row or
!>    Map_Matrix_Col, selected via the optional whichIn argument
!>    (1=row, 2=col). Defaults to 1 (row) when whichIn is absent.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to query.
!>    \endverbatim
!>
!>    \param[in] whichIn
!>    \verbatim
!>        whichIn is Integer,Optional
!>        1 (default) for row, 2 for col.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_get_map(NEFMat,whichIn) result(mapOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      integer,intent(in),optional::whichIn
      integer::which
      type(mqc_scf_eigenvalues),dimension(:),allocatable::mapOut
!
      if(present(whichIn)) then
        which = whichIn
      else
        which = 1
      end if
!
      select case(which)
      case(1)
        mapOut = NEFMat%Map_Matrix_row
      case(2)
        mapOut = NEFMat%Map_Matrix_col
      case default
        call mqc_error_I('mqc_nef_matrix_get_map: whichIn must be 1 (row) or 2 (col)',&
             6,'whichIn',which)
      end select
!
      end function mqc_nef_matrix_get_map
!
!
!     PROCEDURE MQC_NEF_Matrix_Block_Diagonal_PInv
!
!>    \brief <b> MQC_NEF_Matrix_Block_Diagonal_PInv is used to return the
!>    block-diagonal pseudo-inverse of a NEF Matrix object. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Block_Diagonal_PInv computes the pseudo-inverse of
!>    each of the input NEF Matrix's diagonal blocks -- the blocks whose
!>    row space-index equals its column space-index. Off-diagonal blocks
!>    of the result are left zero-flagged, so the returned object is
!>    block-diagonal by construction, regardless of any off-diagonal
!>    content in the input. Zero diagonal blocks are left zero rather
!>    than pseudo-inverted. Same structure as MQC_NEF_R4_Block_Diagonal_
!>    PInv, but operating directly on the Matrix's own row/col block
!>    structure -- no R4-to-Matrix index collapse is needed first.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        NEF Matrix Object whose block-diagonal pseudo-inverse is to be
!>        returned.
!>    \endverbatim
!>
!>    \param[in] thresh
!>    \verbatim
!>        thresh is real(kind=real64),optional
!>        Singular-value threshold passed to each diagonal block's
!>        pseudo-inverse (default 1.0d-8)
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_block_diagonal_pinv(NEFMat,thresh) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      real(kind=real64),intent(in),optional::thresh
      type(mqc_nef_matrix)::NEFMatOut
      real(kind=real64)::myThresh
      integer(kind=int64)::N,id
!
      if(present(thresh)) then
        myThresh = thresh
      else
        myThresh = 1.0d-8
      end if
!
      call NEFMatOut%init(NEFMat%Map_Matrix_Row,&
           initial_map_matrix_col=NEFMat%Map_Matrix_Col)
!
      do N = 1, NEFMat%NumSpace_row
        id = NEFMat%blockIndex([1,N],[2,N])
        if (.not.NEFMat%Zero_Flag(id)) then
          NEFMatOut%Main_Matrix(id) = NEFMat%Main_Matrix(id)%pinv(myThresh)
          NEFMatOut%Zero_Flag(id) = .false.
        end if
      end do
!
      end function mqc_nef_matrix_block_diagonal_pinv
!
!
!     PROCEDURE MQC_NEF_R4_Dimension_Swap
!
!
!>    \brief <b> MQC_NEF_R4_Dimension_Swap is used to return a two-electron integral
!>    after swap two of its dimensions </b> 
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Dimension_Swap is used to return a two-electron integral
!>    after swap two of its dimensions
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is type(mqc_nef_r4) 
!>        The MQC NEF R4 before swapping.
!>    \endverbatim
!>
!>    \param[in] Dim1
!>    \verbatim
!>        Dim1 is integer, between 1 and 4.
!>    \endverbatim
!>
!>    \param[in] Dim2
!>    \verbatim
!>        Dim2 is integer, between 1 and 4.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_dimension_swap(NEFR4,dim1,dim2) result(outR4)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      integer(kind=int64),intent(in)::dim1,dim2
      type(mqc_nef_r4)::outR4
      integer(kind=int64)::mydim1,mydim2,id1,id2,id3,id4,idOld,idNew
!
      if (dim1.le.0.or.dim1.gt.4) call mqc_error_I( &
          'First dimension out of bounds in mqc_nef_r4_dimension_swap',6,&
          'dim1',dim1)
      if (dim2.le.0.or.dim2.gt.4) call mqc_error_I( &
          'Second dimension out of bounds in mqc_nef_r4_dimension_swap',6,&
          'dim2',dim2)
!
      if (dim1.eq.dim2) then 
        outR4 = NEFR4
        return
      else if (dim1.gt.dim2) then
        mydim1 = dim2
        mydim2 = dim1
      else
        mydim1 = dim1
        mydim2 = dim2
      end if
!
      select case (mydim1)
      case (1)
        select case (mydim2)
        case (2)
          outR4%NumSpace_1 = NEFR4%NumSpace_2
          outR4%NumSpace_2 = NEFR4%NumSpace_1
          outR4%NumSpace_3 = NEFR4%NumSpace_3
          outR4%NumSpace_4 = NEFR4%NumSpace_4
          outR4%Map_Matrix_1 = NEFR4%Map_Matrix_2
          outR4%Map_Matrix_2 = NEFR4%Map_Matrix_1
          outR4%Map_Matrix_3 = NEFR4%Map_Matrix_3
          outR4%Map_Matrix_4 = NEFR4%Map_Matrix_4
          outR4%Size_vector_alpha_1 = NEFR4%Size_vector_alpha_2
          outR4%Size_vector_alpha_2 = NEFR4%Size_vector_alpha_1
          outR4%Size_vector_alpha_3 = NEFR4%Size_vector_alpha_3
          outR4%Size_vector_alpha_4 = NEFR4%Size_vector_alpha_4
          outR4%Size_vector_beta_1 = NEFR4%Size_vector_beta_2
          outR4%Size_vector_beta_2 = NEFR4%Size_vector_beta_1
          outR4%Size_vector_beta_3 = NEFR4%Size_vector_beta_3
          outR4%Size_vector_beta_4 = NEFR4%Size_vector_beta_4
          allocate(outR4%Main_R4(size(NEFR4%Main_R4)))
          allocate(outR4%Zero_Flag(size(NEFR4%Main_R4)))
          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  idOld = NEFR4%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
                  idNew = NEFR4%blockIndex([2,id1],[1,id2],[3,id3],[4,id4])
                  outR4%Zero_Flag(idNew) = NEFR4%Zero_Flag(idOld)
                  if (.not.NEFR4%Zero_Flag(idOld)) &
                    outR4%Main_R4(idNew) = NEFR4%Main_R4(idOld)
                end do
              end do
            end do
          end do
        case (3)
          outR4%NumSpace_1 = NEFR4%NumSpace_3
          outR4%NumSpace_2 = NEFR4%NumSpace_2
          outR4%NumSpace_3 = NEFR4%NumSpace_1
          outR4%NumSpace_4 = NEFR4%NumSpace_4
          outR4%Map_Matrix_1 = NEFR4%Map_Matrix_3
          outR4%Map_Matrix_2 = NEFR4%Map_Matrix_2
          outR4%Map_Matrix_3 = NEFR4%Map_Matrix_1
          outR4%Map_Matrix_4 = NEFR4%Map_Matrix_4
          outR4%Size_vector_alpha_1 = NEFR4%Size_vector_alpha_3
          outR4%Size_vector_alpha_2 = NEFR4%Size_vector_alpha_2
          outR4%Size_vector_alpha_3 = NEFR4%Size_vector_alpha_1
          outR4%Size_vector_alpha_4 = NEFR4%Size_vector_alpha_4
          outR4%Size_vector_beta_1 = NEFR4%Size_vector_beta_3
          outR4%Size_vector_beta_2 = NEFR4%Size_vector_beta_2
          outR4%Size_vector_beta_3 = NEFR4%Size_vector_beta_1
          outR4%Size_vector_beta_4 = NEFR4%Size_vector_beta_4
          allocate(outR4%Main_R4(size(NEFR4%Main_R4)))
          allocate(outR4%Zero_Flag(size(NEFR4%Main_R4)))
          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  idOld = NEFR4%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
                  idNew = NEFR4%blockIndex([3,id1],[2,id2],[1,id3],[4,id4])
                  outR4%Zero_Flag(idNew) = NEFR4%Zero_Flag(idOld)
                  if (.not.NEFR4%Zero_Flag(idOld)) &
                    outR4%Main_R4(idNew) = NEFR4%Main_R4(idOld)
                end do
              end do
            end do
          end do
        case (4)
          outR4%NumSpace_1 = NEFR4%NumSpace_4
          outR4%NumSpace_2 = NEFR4%NumSpace_2
          outR4%NumSpace_3 = NEFR4%NumSpace_3
          outR4%NumSpace_4 = NEFR4%NumSpace_1
          outR4%Map_Matrix_1 = NEFR4%Map_Matrix_4
          outR4%Map_Matrix_2 = NEFR4%Map_Matrix_2
          outR4%Map_Matrix_3 = NEFR4%Map_Matrix_3
          outR4%Map_Matrix_4 = NEFR4%Map_Matrix_1
          outR4%Size_vector_alpha_1 = NEFR4%Size_vector_alpha_4
          outR4%Size_vector_alpha_2 = NEFR4%Size_vector_alpha_2
          outR4%Size_vector_alpha_3 = NEFR4%Size_vector_alpha_3
          outR4%Size_vector_alpha_4 = NEFR4%Size_vector_alpha_1
          outR4%Size_vector_beta_1 = NEFR4%Size_vector_beta_4
          outR4%Size_vector_beta_2 = NEFR4%Size_vector_beta_2
          outR4%Size_vector_beta_3 = NEFR4%Size_vector_beta_3
          outR4%Size_vector_beta_4 = NEFR4%Size_vector_beta_1
          allocate(outR4%Main_R4(size(NEFR4%Main_R4)))
          allocate(outR4%Zero_Flag(size(NEFR4%Main_R4)))
          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  idOld = NEFR4%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
                  idNew = NEFR4%blockIndex([4,id1],[2,id2],[3,id3],[1,id4])
                  outR4%Zero_Flag(idNew) = NEFR4%Zero_Flag(idOld)
                  if (.not.NEFR4%Zero_Flag(idOld)) &
                    outR4%Main_R4(idNew) = NEFR4%Main_R4(idOld)
                end do
              end do
            end do
          end do
        end select
      case (2)
        select case (mydim2)
        case (3)
          outR4%NumSpace_1 = NEFR4%NumSpace_1
          outR4%NumSpace_2 = NEFR4%NumSpace_3
          outR4%NumSpace_3 = NEFR4%NumSpace_2
          outR4%NumSpace_4 = NEFR4%NumSpace_4
          outR4%Map_Matrix_1 = NEFR4%Map_Matrix_1
          outR4%Map_Matrix_2 = NEFR4%Map_Matrix_3
          outR4%Map_Matrix_3 = NEFR4%Map_Matrix_2
          outR4%Map_Matrix_4 = NEFR4%Map_Matrix_4
          outR4%Size_vector_alpha_1 = NEFR4%Size_vector_alpha_1
          outR4%Size_vector_alpha_2 = NEFR4%Size_vector_alpha_3
          outR4%Size_vector_alpha_3 = NEFR4%Size_vector_alpha_2
          outR4%Size_vector_alpha_4 = NEFR4%Size_vector_alpha_4
          outR4%Size_vector_beta_1 = NEFR4%Size_vector_beta_1
          outR4%Size_vector_beta_2 = NEFR4%Size_vector_beta_3
          outR4%Size_vector_beta_3 = NEFR4%Size_vector_beta_2
          outR4%Size_vector_beta_4 = NEFR4%Size_vector_beta_4
          allocate(outR4%Main_R4(size(NEFR4%Main_R4)))
          allocate(outR4%Zero_Flag(size(NEFR4%Main_R4)))
          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  idOld = NEFR4%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
                  idNew = NEFR4%blockIndex([1,id1],[3,id2],[2,id3],[4,id4])
                  outR4%Zero_Flag(idNew) = NEFR4%Zero_Flag(idOld)
                  if (.not.NEFR4%Zero_Flag(idOld)) &
                    outR4%Main_R4(idNew) = NEFR4%Main_R4(idOld)
                end do
              end do
            end do
          end do
        case (4)
          outR4%NumSpace_1 = NEFR4%NumSpace_1
          outR4%NumSpace_2 = NEFR4%NumSpace_4
          outR4%NumSpace_3 = NEFR4%NumSpace_3
          outR4%NumSpace_4 = NEFR4%NumSpace_2
          outR4%Map_Matrix_1 = NEFR4%Map_Matrix_1
          outR4%Map_Matrix_2 = NEFR4%Map_Matrix_4
          outR4%Map_Matrix_3 = NEFR4%Map_Matrix_3
          outR4%Map_Matrix_4 = NEFR4%Map_Matrix_2
          outR4%Size_vector_alpha_1 = NEFR4%Size_vector_alpha_1
          outR4%Size_vector_alpha_2 = NEFR4%Size_vector_alpha_4
          outR4%Size_vector_alpha_3 = NEFR4%Size_vector_alpha_3
          outR4%Size_vector_alpha_4 = NEFR4%Size_vector_alpha_2
          outR4%Size_vector_beta_1 = NEFR4%Size_vector_beta_1
          outR4%Size_vector_beta_2 = NEFR4%Size_vector_beta_4
          outR4%Size_vector_beta_3 = NEFR4%Size_vector_beta_3
          outR4%Size_vector_beta_4 = NEFR4%Size_vector_beta_2
          allocate(outR4%Main_R4(size(NEFR4%Main_R4)))
          allocate(outR4%Zero_Flag(size(NEFR4%Main_R4)))
          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  idOld = NEFR4%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
                  idNew = NEFR4%blockIndex([1,id1],[4,id2],[3,id3],[2,id4])
                  outR4%Zero_Flag(idNew) = NEFR4%Zero_Flag(idOld)
                  if (.not.NEFR4%Zero_Flag(idOld)) &
                    outR4%Main_R4(idNew) = NEFR4%Main_R4(idOld)
                end do
              end do
            end do
          end do
        end select
      case (3)
        outR4%NumSpace_1 = NEFR4%NumSpace_1
        outR4%NumSpace_2 = NEFR4%NumSpace_2
        outR4%NumSpace_3 = NEFR4%NumSpace_4
        outR4%NumSpace_4 = NEFR4%NumSpace_3
        outR4%Map_Matrix_1 = NEFR4%Map_Matrix_1
        outR4%Map_Matrix_2 = NEFR4%Map_Matrix_2
        outR4%Map_Matrix_3 = NEFR4%Map_Matrix_4
        outR4%Map_Matrix_4 = NEFR4%Map_Matrix_3
        outR4%Size_vector_alpha_1 = NEFR4%Size_vector_alpha_1
        outR4%Size_vector_alpha_2 = NEFR4%Size_vector_alpha_2
        outR4%Size_vector_alpha_3 = NEFR4%Size_vector_alpha_4
        outR4%Size_vector_alpha_4 = NEFR4%Size_vector_alpha_3
        outR4%Size_vector_beta_1 = NEFR4%Size_vector_beta_1
        outR4%Size_vector_beta_2 = NEFR4%Size_vector_beta_2
        outR4%Size_vector_beta_3 = NEFR4%Size_vector_beta_4
        outR4%Size_vector_beta_4 = NEFR4%Size_vector_beta_3
        allocate(outR4%Main_R4(size(NEFR4%Main_R4)))
        allocate(outR4%Zero_Flag(size(NEFR4%Main_R4)))
        do id1 = 1, NEFR4%NumSpace_1
          do id2 = 1, NEFR4%NumSpace_2
            do id3 = 1, NEFR4%NumSpace_3
              do id4 = 1, NEFR4%NumSpace_4
                idOld = NEFR4%blockIndex([1,id1],[2,id2],[3,id3],[4,id4])
                idNew = NEFR4%blockIndex([1,id1],[2,id2],[4,id3],[3,id4])
                outR4%Zero_Flag(idNew) = NEFR4%Zero_Flag(idOld)
                if (.not.NEFR4%Zero_Flag(idOld)) &
                  outR4%Main_R4(idNew) = NEFR4%Main_R4(idOld)
              end do
            end do
          end do
        end do
      end select
!
      end function mqc_nef_r4_dimension_swap
!
!----------------------------------------------------------------
!                                                               |
!     MQC_NEF_SET PROCEDURES                                   |
!                                                               |
!----------------------------------------------------------------
!
!
!     PROCEDURE MQC_NEF_Set_Initialize
!
!
!>    \brief <b> MQC_NEF_Set_Initialize is a subroutine that initializes an
!>    MQC_NEF_Set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Initialize records the (I,J) grid dimensions of a
!>    MQC_NEF_Set and marks it 'open' -- no payload array (Vector, Matrix,
!>    or R4) is allocated yet. The first call to put (putVector/putMatrix/
!>    putR4) determines which kind the set holds, allocates the matching
!>    payload array to the recorded dimensions, and fixes the set's kind
!>    for the rest of its lifetime.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to initialize.
!>    \endverbatim
!>
!>    \param[in] nI
!>    \verbatim
!>        nI is Integer
!>        Size of the set's first (I) grid dimension.
!>    \endverbatim
!>
!>    \param[in] nJ
!>    \verbatim
!>        nJ is Integer
!>        Size of the set's second (J) grid dimension.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_set_initialize(NEFSet,nI,nJ)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      integer,intent(in)::nI,nJ
!
      NEFSet%NumI = nI
      NEFSet%NumJ = nJ
      NEFSet%Array_Type = 'open'
!
      end subroutine mqc_nef_set_initialize
!
!
!     PROCEDURE MQC_NEF_Set_Is_Open
!
!
!>    \brief <b> MQC_NEF_Set_Is_Open returns whether a MQC_NEF_Set has not
!>    yet had an element written to it. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Is_Open checks whether the set's Array_Type is still
!>    'open' (i.e. init has been called but no put has happened yet).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to query.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_is_open(NEFSet) result(isOpen)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      logical::isOpen
!
      isOpen = (trim(NEFSet%Array_Type).eq.'open')
!
      end function mqc_nef_set_is_open
!
!
!     PROCEDURE MQC_NEF_Set_Is_Zero
!
!
!>    \brief <b> MQC_NEF_Set_Is_Zero returns whether a MQC_NEF_Set's
!>    (I,J) cell has never been written. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Is_Zero checks the set's per-cell Zero_Flag(I,J),
!>    which starts .true. for every cell when the set is first
!>    allocated (on the first put) and is flipped to .false. the
!>    moment that specific (I,J) cell is written via put. A set that
!>    is still 'open' (Zero_Flag not yet allocated) reports every
!>    cell as zero.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to query.
!>    \endverbatim
!>
!>    \param[in] I
!>    \verbatim
!>        I is Integer
!>        First grid index.
!>    \endverbatim
!>
!>    \param[in] J
!>    \verbatim
!>        J is Integer
!>        Second grid index.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_is_zero(NEFSet,I,J) result(isZero)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      integer,intent(in)::I,J
      logical::isZero
!
      isZero = .true.
      if(allocated(NEFSet%Zero_Flag)) isZero = NEFSet%Zero_Flag(I,J)
!
      end function mqc_nef_set_is_zero
!
!
!     PROCEDURE MQC_NEF_Set_Is_Vector
!
!
!>    \brief <b> MQC_NEF_Set_Is_Vector returns whether a MQC_NEF_Set
!>    currently holds MQC_NEF_Vector elements. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Is_Vector checks whether the set's Array_Type is
!>    'vector'.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to query.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_is_vector(NEFSet) result(isVector)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      logical::isVector
!
      isVector = (trim(NEFSet%Array_Type).eq.'vector')
!
      end function mqc_nef_set_is_vector
!
!
!     PROCEDURE MQC_NEF_Set_Is_Matrix
!
!
!>    \brief <b> MQC_NEF_Set_Is_Matrix returns whether a MQC_NEF_Set
!>    currently holds MQC_NEF_Matrix elements. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Is_Matrix checks whether the set's Array_Type is
!>    'matrix'.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to query.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_is_matrix(NEFSet) result(isMatrix)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      logical::isMatrix
!
      isMatrix = (trim(NEFSet%Array_Type).eq.'matrix')
!
      end function mqc_nef_set_is_matrix
!
!
!     PROCEDURE MQC_NEF_Set_Is_R4
!
!
!>    \brief <b> MQC_NEF_Set_Is_R4 returns whether a MQC_NEF_Set currently
!>    holds MQC_NEF_R4 elements. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Is_R4 checks whether the set's Array_Type is 'r4'.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to query.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_is_r4(NEFSet) result(isR4)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      logical::isR4
!
      isR4 = (trim(NEFSet%Array_Type).eq.'r4')
!
      end function mqc_nef_set_is_r4
!
!
!     PROCEDURE MQC_NEF_Set_Array_Type
!
!
!>    \brief <b> MQC_NEF_Set_Array_Type returns the kind tag of a
!>    MQC_NEF_Set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Array_Type returns the set's Array_Type string
!>    ('open', 'vector', 'matrix', or 'r4'), mirroring
!>    mqc_integral_array_type's role for mqc_scf_integral.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to query.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_array_type(NEFSet) result(arrayType)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      character(len=64)::arrayType
!
      arrayType = NEFSet%Array_Type
!
      end function mqc_nef_set_array_type
!
!
!     PROCEDURE MQC_NEF_Set_Allocate
!
!
!>    \brief <b> MQC_NEF_Set_Allocate is a subroutine that allocates a
!>    MQC_NEF_Set's payload array for a given kind and records its
!>    dimension signature. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Allocate allocates the payload array (Set_Vector,
!>    Set_Matrix, or Set_R4) matching kind, sized from the dimensions
!>    recorded by init, and tags the set accordingly. Exactly which of the
!>    eight Size_vector_alpha/beta_1..4 arguments are required depends on
!>    kind: 'vector' requires only axis 1 (size_vector_alpha_1,
!>    size_vector_beta_1) and rejects any axis-2/3/4 argument being
!>    present; 'matrix' requires axes 1 and 2 (row/col) and rejects
!>    axis-3/4; 'r4' requires all eight (axes 1-4). Within every required
!>    axis, the alpha and beta arrays must have the same length. That
!>    length becomes the axis's NumSpace entry; axes not used by kind stay
!>    NumSpace=-1. This is called once, from put, the first time an
!>    element is written to an 'open' set.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to allocate.
!>    \endverbatim
!>
!>    \param[in] kind
!>    \verbatim
!>        kind is Character(len=*)
!>        Which payload to allocate: 'vector', 'matrix', or 'r4'.
!>    \endverbatim
!>
!>    \param[in] size_vector_alpha_1, size_vector_beta_1
!>    \verbatim
!>        Integer,Dimension(:),Optional
!>        Axis-1 per-space sizes. Required for all three kinds (Vector's
!>        only axis; Matrix's row; R4's axis 1). Must have equal length.
!>    \endverbatim
!>
!>    \param[in] size_vector_alpha_2, size_vector_beta_2
!>    \verbatim
!>        Integer,Dimension(:),Optional
!>        Axis-2 per-space sizes. Required for 'matrix' (column) and 'r4'
!>        (axis 2); rejected for 'vector'. Must have equal length.
!>    \endverbatim
!>
!>    \param[in] size_vector_alpha_3, size_vector_beta_3
!>    \verbatim
!>        Integer,Dimension(:),Optional
!>        Axis-3 per-space sizes. Required for 'r4' only; rejected for
!>        'vector'/'matrix'. Must have equal length.
!>    \endverbatim
!>
!>    \param[in] size_vector_alpha_4, size_vector_beta_4
!>    \verbatim
!>        Integer,Dimension(:),Optional
!>        Axis-4 per-space sizes. Required for 'r4' only; rejected for
!>        'vector'/'matrix'. Must have equal length.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_set_allocate(NEFSet,kind,&
           size_vector_alpha_1,size_vector_beta_1,&
           size_vector_alpha_2,size_vector_beta_2,&
           size_vector_alpha_3,size_vector_beta_3,&
           size_vector_alpha_4,size_vector_beta_4)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      character(len=*),intent(in)::kind
      integer,dimension(:),intent(in),optional::size_vector_alpha_1,size_vector_beta_1
      integer,dimension(:),intent(in),optional::size_vector_alpha_2,size_vector_beta_2
      integer,dimension(:),intent(in),optional::size_vector_alpha_3,size_vector_beta_3
      integer,dimension(:),intent(in),optional::size_vector_alpha_4,size_vector_beta_4
!
      select case(trim(kind))
      case('vector','Vector')
        if(.not.(present(size_vector_alpha_1).and.present(size_vector_beta_1))) &
          call mqc_error('mqc_nef_set_allocate: vector requires '//&
               'size_vector_alpha_1 and size_vector_beta_1',6)
        if(present(size_vector_alpha_2).or.present(size_vector_beta_2).or.&
           present(size_vector_alpha_3).or.present(size_vector_beta_3).or.&
           present(size_vector_alpha_4).or.present(size_vector_beta_4)) &
          call mqc_error('mqc_nef_set_allocate: vector only accepts axis-1 size_vectors',6)
        if(size(size_vector_alpha_1).ne.size(size_vector_beta_1)) &
          call mqc_error('mqc_nef_set_allocate: size_vector_alpha_1 and '//&
               'size_vector_beta_1 must have the same dimension',6)
        allocate(NEFSet%Set_Vector(NEFSet%NumI,NEFSet%NumJ))
        NEFSet%Size_vector_alpha_1 = size_vector_alpha_1
        NEFSet%Size_vector_beta_1 = size_vector_beta_1
        NEFSet%NumSpace(1) = size(size_vector_alpha_1)
      case('matrix','Matrix')
        if(.not.(present(size_vector_alpha_1).and.present(size_vector_beta_1).and.&
                 present(size_vector_alpha_2).and.present(size_vector_beta_2))) &
          call mqc_error('mqc_nef_set_allocate: matrix requires '//&
               'axis-1 and axis-2 size_vectors',6)
        if(present(size_vector_alpha_3).or.present(size_vector_beta_3).or.&
           present(size_vector_alpha_4).or.present(size_vector_beta_4)) &
          call mqc_error('mqc_nef_set_allocate: matrix only accepts '//&
               'axis-1 and axis-2 size_vectors',6)
        if(size(size_vector_alpha_1).ne.size(size_vector_beta_1)) &
          call mqc_error('mqc_nef_set_allocate: size_vector_alpha_1 and '//&
               'size_vector_beta_1 must have the same dimension',6)
        if(size(size_vector_alpha_2).ne.size(size_vector_beta_2)) &
          call mqc_error('mqc_nef_set_allocate: size_vector_alpha_2 and '//&
               'size_vector_beta_2 must have the same dimension',6)
        allocate(NEFSet%Set_Matrix(NEFSet%NumI,NEFSet%NumJ))
        NEFSet%Size_vector_alpha_1 = size_vector_alpha_1
        NEFSet%Size_vector_beta_1 = size_vector_beta_1
        NEFSet%NumSpace(1) = size(size_vector_alpha_1)
        NEFSet%Size_vector_alpha_2 = size_vector_alpha_2
        NEFSet%Size_vector_beta_2 = size_vector_beta_2
        NEFSet%NumSpace(2) = size(size_vector_alpha_2)
      case('r4','R4')
        if(.not.(present(size_vector_alpha_1).and.present(size_vector_beta_1).and.&
                 present(size_vector_alpha_2).and.present(size_vector_beta_2).and.&
                 present(size_vector_alpha_3).and.present(size_vector_beta_3).and.&
                 present(size_vector_alpha_4).and.present(size_vector_beta_4))) &
          call mqc_error('mqc_nef_set_allocate: r4 requires all eight size_vector arguments',6)
        if(size(size_vector_alpha_1).ne.size(size_vector_beta_1)) &
          call mqc_error('mqc_nef_set_allocate: size_vector_alpha_1 and '//&
               'size_vector_beta_1 must have the same dimension',6)
        if(size(size_vector_alpha_2).ne.size(size_vector_beta_2)) &
          call mqc_error('mqc_nef_set_allocate: size_vector_alpha_2 and '//&
               'size_vector_beta_2 must have the same dimension',6)
        if(size(size_vector_alpha_3).ne.size(size_vector_beta_3)) &
          call mqc_error('mqc_nef_set_allocate: size_vector_alpha_3 and '//&
               'size_vector_beta_3 must have the same dimension',6)
        if(size(size_vector_alpha_4).ne.size(size_vector_beta_4)) &
          call mqc_error('mqc_nef_set_allocate: size_vector_alpha_4 and '//&
               'size_vector_beta_4 must have the same dimension',6)
        allocate(NEFSet%Set_R4(NEFSet%NumI,NEFSet%NumJ))
        NEFSet%Size_vector_alpha_1 = size_vector_alpha_1
        NEFSet%Size_vector_beta_1 = size_vector_beta_1
        NEFSet%NumSpace(1) = size(size_vector_alpha_1)
        NEFSet%Size_vector_alpha_2 = size_vector_alpha_2
        NEFSet%Size_vector_beta_2 = size_vector_beta_2
        NEFSet%NumSpace(2) = size(size_vector_alpha_2)
        NEFSet%Size_vector_alpha_3 = size_vector_alpha_3
        NEFSet%Size_vector_beta_3 = size_vector_beta_3
        NEFSet%NumSpace(3) = size(size_vector_alpha_3)
        NEFSet%Size_vector_alpha_4 = size_vector_alpha_4
        NEFSet%Size_vector_beta_4 = size_vector_beta_4
        NEFSet%NumSpace(4) = size(size_vector_alpha_4)
      case default
        call mqc_error_A('mqc_nef_set_allocate: unknown kind',6,'kind',trim(kind))
      end select
      allocate(NEFSet%Zero_Flag(NEFSet%NumI,NEFSet%NumJ))
      NEFSet%Zero_Flag = .true.
      NEFSet%Array_Type = trim(kind)
!
      end subroutine mqc_nef_set_allocate
!
!
!     PROCEDURE MQC_NEF_Set_Put_Vector
!
!
!>    \brief <b> MQC_NEF_Set_Put_Vector writes a MQC_NEF_Vector element
!>    into a MQC_NEF_Set at (I,J). </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Put_Vector allocates the set's Vector payload array via
!>    mqc_nef_set_allocate (sized from the dimensions recorded by init,
!>    tagged 'vector') if the set is still 'open'. If the set already
!>    holds a different kind, an error is raised. If the set is already
!>    'vector', the incoming element's Size_vector_alpha/beta are checked
!>    against the set's recorded dimensions and an error is raised on a
!>    mismatch. Otherwise, the element is stored at (I,J). Bound to the
!>    set's Generic put, so callers normally invoke this through
!>    NEFSet%put(element,I,J).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to write into.
!>    \endverbatim
!>
!>    \param[in] I
!>    \verbatim
!>        I is Integer
!>        First grid index.
!>    \endverbatim
!>
!>    \param[in] J
!>    \verbatim
!>        J is Integer
!>        Second grid index.
!>    \endverbatim
!>
!>    \param[in] element
!>    \verbatim
!>        element is Type(MQC_NEF_Vector)
!>        The value to store at (I,J).
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_set_put_vector(NEFSet,element,I,J)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      integer,intent(in)::I,J
      type(mqc_nef_vector),intent(in)::element
!
      if(NEFSet%isOpen()) then
        call mqc_nef_set_allocate(NEFSet,'vector',&
             size_vector_alpha_1=element%Size_vector_alpha,&
             size_vector_beta_1=element%Size_vector_beta)
      else if(.not.NEFSet%isVector()) then
        call mqc_error_A('mqc_nef_set_put_vector: NEF_Set already holds a different kind',&
             6,'Array_Type',NEFSet%Array_Type)
      else if(size(NEFSet%Size_vector_alpha_1).ne.size(element%Size_vector_alpha).or.&
              size(NEFSet%Size_vector_beta_1).ne.size(element%Size_vector_beta)) then
        call mqc_error('mqc_nef_set_put_vector: element dimensions do not '//&
             'match the set''s established dimensions',6)
      else if(any(NEFSet%Size_vector_alpha_1.ne.element%Size_vector_alpha).or.&
              any(NEFSet%Size_vector_beta_1.ne.element%Size_vector_beta)) then
        call mqc_error('mqc_nef_set_put_vector: element dimensions do not '//&
             'match the set''s established dimensions',6)
      end if
      NEFSet%Set_Vector(I,J) = element
      NEFSet%Zero_Flag(I,J) = .false.
!
      end subroutine mqc_nef_set_put_vector
!
!
!     PROCEDURE MQC_NEF_Set_Put_Matrix
!
!
!>    \brief <b> MQC_NEF_Set_Put_Matrix writes a MQC_NEF_Matrix element
!>    into a MQC_NEF_Set at (I,J). </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Put_Matrix allocates the set's Matrix payload array via
!>    mqc_nef_set_allocate (sized from the dimensions recorded by init,
!>    tagged 'matrix') if the set is still 'open'. If the set already
!>    holds a different kind, an error is raised. If the set is already
!>    'matrix', the incoming element's row/col Size_vector arrays are
!>    checked against the set's recorded dimensions and an error is
!>    raised on a mismatch. Otherwise, the element is stored at (I,J).
!>    Bound to the set's Generic put, so callers normally invoke this
!>    through NEFSet%put(element,I,J).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to write into.
!>    \endverbatim
!>
!>    \param[in] I
!>    \verbatim
!>        I is Integer
!>        First grid index.
!>    \endverbatim
!>
!>    \param[in] J
!>    \verbatim
!>        J is Integer
!>        Second grid index.
!>    \endverbatim
!>
!>    \param[in] element
!>    \verbatim
!>        element is Type(MQC_NEF_Matrix)
!>        The value to store at (I,J).
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_set_put_matrix(NEFSet,element,I,J)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      integer,intent(in)::I,J
      type(mqc_nef_matrix),intent(in)::element
!
      if(NEFSet%isOpen()) then
        call mqc_nef_set_allocate(NEFSet,'matrix',&
             size_vector_alpha_1=element%Size_vector_alpha_row,&
             size_vector_beta_1=element%Size_vector_beta_row,&
             size_vector_alpha_2=element%Size_vector_alpha_col,&
             size_vector_beta_2=element%Size_vector_beta_col)
      else if(.not.NEFSet%isMatrix()) then
        call mqc_error_A('mqc_nef_set_put_matrix: NEF_Set already holds a different kind',&
             6,'Array_Type',NEFSet%Array_Type)
      else if(size(NEFSet%Size_vector_alpha_1).ne.size(element%Size_vector_alpha_row).or.&
              size(NEFSet%Size_vector_beta_1).ne.size(element%Size_vector_beta_row).or.&
              size(NEFSet%Size_vector_alpha_2).ne.size(element%Size_vector_alpha_col).or.&
              size(NEFSet%Size_vector_beta_2).ne.size(element%Size_vector_beta_col)) then
        call mqc_error('mqc_nef_set_put_matrix: element dimensions do not '//&
             'match the set''s established dimensions',6)
      else if(any(NEFSet%Size_vector_alpha_1.ne.element%Size_vector_alpha_row).or.&
              any(NEFSet%Size_vector_beta_1.ne.element%Size_vector_beta_row).or.&
              any(NEFSet%Size_vector_alpha_2.ne.element%Size_vector_alpha_col).or.&
              any(NEFSet%Size_vector_beta_2.ne.element%Size_vector_beta_col)) then
        call mqc_error('mqc_nef_set_put_matrix: element dimensions do not '//&
             'match the set''s established dimensions',6)
      end if
      NEFSet%Set_Matrix(I,J) = element
      NEFSet%Zero_Flag(I,J) = .false.
!
      end subroutine mqc_nef_set_put_matrix
!
!
!     PROCEDURE MQC_NEF_Set_Put_R4
!
!
!>    \brief <b> MQC_NEF_Set_Put_R4 writes a MQC_NEF_R4 element into a
!>    MQC_NEF_Set at (I,J). </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Put_R4 allocates the set's R4 payload array via
!>    mqc_nef_set_allocate (sized from the dimensions recorded by init,
!>    tagged 'r4') if the set is still 'open'. If the set already holds a
!>    different kind, an error is raised. If the set is already 'r4', the
!>    incoming element's four axes' Size_vector arrays are checked against
!>    the set's recorded dimensions and an error is raised on a mismatch.
!>    Otherwise, the element is stored at (I,J). Bound to the set's
!>    Generic put, so callers normally invoke this through
!>    NEFSet%put(element,I,J).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to write into.
!>    \endverbatim
!>
!>    \param[in] I
!>    \verbatim
!>        I is Integer
!>        First grid index.
!>    \endverbatim
!>
!>    \param[in] J
!>    \verbatim
!>        J is Integer
!>        Second grid index.
!>    \endverbatim
!>
!>    \param[in] element
!>    \verbatim
!>        element is Type(MQC_NEF_R4)
!>        The value to store at (I,J).
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_set_put_r4(NEFSet,element,I,J)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      integer,intent(in)::I,J
      type(mqc_nef_r4),intent(in)::element
!
      if(NEFSet%isOpen()) then
        call mqc_nef_set_allocate(NEFSet,'r4',&
             size_vector_alpha_1=element%Size_vector_alpha_1,&
             size_vector_beta_1=element%Size_vector_beta_1,&
             size_vector_alpha_2=element%Size_vector_alpha_2,&
             size_vector_beta_2=element%Size_vector_beta_2,&
             size_vector_alpha_3=element%Size_vector_alpha_3,&
             size_vector_beta_3=element%Size_vector_beta_3,&
             size_vector_alpha_4=element%Size_vector_alpha_4,&
             size_vector_beta_4=element%Size_vector_beta_4)
      else if(.not.NEFSet%isR4()) then
        call mqc_error_A('mqc_nef_set_put_r4: NEF_Set already holds a different kind',&
             6,'Array_Type',NEFSet%Array_Type)
      else if(size(NEFSet%Size_vector_alpha_1).ne.size(element%Size_vector_alpha_1).or.&
              size(NEFSet%Size_vector_beta_1).ne.size(element%Size_vector_beta_1).or.&
              size(NEFSet%Size_vector_alpha_2).ne.size(element%Size_vector_alpha_2).or.&
              size(NEFSet%Size_vector_beta_2).ne.size(element%Size_vector_beta_2).or.&
              size(NEFSet%Size_vector_alpha_3).ne.size(element%Size_vector_alpha_3).or.&
              size(NEFSet%Size_vector_beta_3).ne.size(element%Size_vector_beta_3).or.&
              size(NEFSet%Size_vector_alpha_4).ne.size(element%Size_vector_alpha_4).or.&
              size(NEFSet%Size_vector_beta_4).ne.size(element%Size_vector_beta_4)) then
        call mqc_error('mqc_nef_set_put_r4: element dimensions do not '//&
             'match the set''s established dimensions',6)
      else if(any(NEFSet%Size_vector_alpha_1.ne.element%Size_vector_alpha_1).or.&
              any(NEFSet%Size_vector_beta_1.ne.element%Size_vector_beta_1).or.&
              any(NEFSet%Size_vector_alpha_2.ne.element%Size_vector_alpha_2).or.&
              any(NEFSet%Size_vector_beta_2.ne.element%Size_vector_beta_2).or.&
              any(NEFSet%Size_vector_alpha_3.ne.element%Size_vector_alpha_3).or.&
              any(NEFSet%Size_vector_beta_3.ne.element%Size_vector_beta_3).or.&
              any(NEFSet%Size_vector_alpha_4.ne.element%Size_vector_alpha_4).or.&
              any(NEFSet%Size_vector_beta_4.ne.element%Size_vector_beta_4)) then
        call mqc_error('mqc_nef_set_put_r4: element dimensions do not '//&
             'match the set''s established dimensions',6)
      end if
      NEFSet%Set_R4(I,J) = element
      NEFSet%Zero_Flag(I,J) = .false.
!
      end subroutine mqc_nef_set_put_r4
!
!
!     PROCEDURE MQC_NEF_Set_At_Vector
!
!
!>    \brief <b> MQC_NEF_Set_At_Vector returns the MQC_NEF_Vector element
!>    stored at (I,J) in a MQC_NEF_Set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_At_Vector returns a copy of the element stored at
!>    (I,J). An error is raised if the set does not currently hold
!>    vector elements (including if it is still 'open').
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to read from.
!>    \endverbatim
!>
!>    \param[in] I
!>    \verbatim
!>        I is Integer
!>        First grid index.
!>    \endverbatim
!>
!>    \param[in] J
!>    \verbatim
!>        J is Integer
!>        Second grid index.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_at_vector(NEFSet,I,J) result(element)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      integer,intent(in)::I,J
      type(mqc_nef_vector)::element
!
      if(.not.NEFSet%isVector()) call mqc_error_A(&
           'mqc_nef_set_at_vector: NEF_Set does not hold vector elements (or is still open)',&
           6,'Array_Type',NEFSet%Array_Type)
      element = NEFSet%Set_Vector(I,J)
!
      end function mqc_nef_set_at_vector
!
!
!     PROCEDURE MQC_NEF_Set_At_Matrix
!
!
!>    \brief <b> MQC_NEF_Set_At_Matrix returns the MQC_NEF_Matrix element
!>    stored at (I,J) in a MQC_NEF_Set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_At_Matrix returns a copy of the element stored at
!>    (I,J). An error is raised if the set does not currently hold
!>    matrix elements (including if it is still 'open').
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to read from.
!>    \endverbatim
!>
!>    \param[in] I
!>    \verbatim
!>        I is Integer
!>        First grid index.
!>    \endverbatim
!>
!>    \param[in] J
!>    \verbatim
!>        J is Integer
!>        Second grid index.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_at_matrix(NEFSet,I,J) result(element)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      integer,intent(in)::I,J
      type(mqc_nef_matrix)::element
!
      if(.not.NEFSet%isMatrix()) call mqc_error_A(&
           'mqc_nef_set_at_matrix: NEF_Set does not hold matrix elements (or is still open)',&
           6,'Array_Type',NEFSet%Array_Type)
      element = NEFSet%Set_Matrix(I,J)
!
      end function mqc_nef_set_at_matrix
!
!
!     PROCEDURE MQC_NEF_Set_At_R4
!
!
!>    \brief <b> MQC_NEF_Set_At_R4 returns the MQC_NEF_R4 element stored
!>    at (I,J) in a MQC_NEF_Set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_At_R4 returns a copy of the element stored at (I,J).
!>    An error is raised if the set does not currently hold R4 elements
!>    (including if it is still 'open').
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to read from.
!>    \endverbatim
!>
!>    \param[in] I
!>    \verbatim
!>        I is Integer
!>        First grid index.
!>    \endverbatim
!>
!>    \param[in] J
!>    \verbatim
!>        J is Integer
!>        Second grid index.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_at_r4(NEFSet,I,J) result(element)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      integer,intent(in)::I,J
      type(mqc_nef_r4)::element
!
      if(.not.NEFSet%isR4()) call mqc_error_A(&
           'mqc_nef_set_at_r4: NEF_Set does not hold r4 elements (or is still open)',&
           6,'Array_Type',NEFSet%Array_Type)
      element = NEFSet%Set_R4(I,J)
!
      end function mqc_nef_set_at_r4
!
!
!     PROCEDURE MQC_NEF_Set_Get_Map
!
!
!>    \brief <b> MQC_NEF_Set_Get_Map returns one of a matrix-kind
!>    MQC_NEF_Set's cells' internal block maps. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Get_Map requires NEFSet to hold matrix elements, and
!>    returns the (I,J) cell's Map_Matrix_Row or Map_Matrix_Col (via that
!>    cell's own getMap), selected via the optional whichIn argument
!>    (1=row, 2=col, defaulting to 1). Maps are stored per-cell (not at
!>    the grid level), so (I,J) selects which cell's map is returned;
!>    every cell shares the same axis sizes by construction, but only a
!>    specific cell's map is a real MQC_SCF_Eigenvalues object.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The (I,J) grid to query. Must hold matrix elements.
!>    \endverbatim
!>
!>    \param[in] I
!>    \verbatim
!>        I is Integer
!>        First grid index.
!>    \endverbatim
!>
!>    \param[in] J
!>    \verbatim
!>        J is Integer
!>        Second grid index.
!>    \endverbatim
!>
!>    \param[in] whichIn
!>    \verbatim
!>        whichIn is Integer,Optional
!>        1 (default) for row, 2 for col.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_get_map(NEFSet,I,J,whichIn) result(mapOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      integer,intent(in)::I,J
      integer,intent(in),optional::whichIn
      type(mqc_scf_eigenvalues),dimension(:),allocatable::mapOut
!
      if(trim(NEFSet%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nef_set_get_map: NEF_Set must hold matrix elements',&
           6,'Array_Type',NEFSet%Type())
!
      if(present(whichIn)) then
        mapOut = mqc_nef_matrix_get_map(NEFSet%atMatrix(I,J),whichIn)
      else
        mapOut = mqc_nef_matrix_get_map(NEFSet%atMatrix(I,J))
      end if
!
      end function mqc_nef_set_get_map
!
!
!     PROCEDURE MQC_NEFSet_Vector_Multiply
!
!
!>    \brief <b> MQC_NEFSet_Vector_Multiply contracts a MQC_NEF_Set's
!>    second (J) axis against a plain mqc_vector. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFSet_Vector_Multiply computes, for every I,
!>    NEFSetOut(I,1) = Sum_J NEFSet(I,J)*Vec%at(J), where Vec is a plain
!>    mqc_vector (not a MQC_NEF_Vector). The result is a MQC_NEF_Set of
!>    the same kind as NEFSetA, with NumJ collapsed to 1. NEFSetA must
!>    already hold an element (not be 'open'). Arguments are named
!>    NEFSetA/VecB (rather than NEFSet/Vec) to avoid name collisions with
!>    the reversed-operand MatMul(Vec,NEFSet) procedure.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>        The (I,J) grid to contract.
!>    \endverbatim
!>
!>    \param[in] VecB
!>    \verbatim
!>        VecB is Type(MQC_Vector)
!>        The plain vector to contract against NEFSetA's J axis. Its
!>        length must equal NEFSetA's NumJ.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_vector_multiply(NEFSetA,VecB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_vector),intent(in)::VecB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
      type(mqc_nef_vector)::tmpVector
      type(mqc_nef_matrix)::tmpMatrix
      type(mqc_nef_r4)::tmpR4
      logical::first
!
      call NEFSetOut%init(NEFSetA%NumI,1)
!
      select case(trim(NEFSetA%Type()))
      case('vector')
        do I = 1, NEFSetA%NumI
          first = .true.
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J)) cycle
            if(first) then
              tmpVector = NEFSetA%atVector(I,J)*VecB%at(J)
              first = .false.
            else
              tmpVector = tmpVector + NEFSetA%atVector(I,J)*VecB%at(J)
            end if
          end do
          if(.not.first) call NEFSetOut%put(tmpVector,I,1)
        end do
      case('matrix')
        do I = 1, NEFSetA%NumI
          first = .true.
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J)) cycle
            if(first) then
              tmpMatrix = NEFSetA%atMatrix(I,J)*VecB%at(J)
              first = .false.
            else
              tmpMatrix = tmpMatrix + NEFSetA%atMatrix(I,J)*VecB%at(J)
            end if
          end do
          if(.not.first) call NEFSetOut%put(tmpMatrix,I,1)
        end do
      case('r4')
        do I = 1, NEFSetA%NumI
          first = .true.
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J)) cycle
            if(first) then
              tmpR4 = NEFSetA%atR4(I,J)*VecB%at(J)
              first = .false.
            else
              tmpR4 = tmpR4 + NEFSetA%atR4(I,J)*VecB%at(J)
            end if
          end do
          if(.not.first) call NEFSetOut%put(tmpR4,I,1)
        end do
      case default
        call mqc_error_A('mqc_nefset_vector_multiply: NEF_Set is open (holds no elements)',&
             6,'Array_Type',NEFSetA%Type())
      end select
!
      end function mqc_nefset_vector_multiply
!
!
!     PROCEDURE MQC_Vector_NEFSet_Multiply
!
!
!>    \brief <b> MQC_Vector_NEFSet_Multiply contracts a MQC_NEF_Set's
!>    first (I) axis against a plain mqc_vector. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Vector_NEFSet_Multiply computes, for every J,
!>    NEFSetOut(1,J) = Sum_I VecA%at(I)*NEFSetB(I,J), where VecA is a
!>    plain mqc_vector (not a MQC_NEF_Vector). The result is a
!>    MQC_NEF_Set of the same kind as NEFSetB, with NumI collapsed to 1.
!>    NEFSetB must already hold an element (not be 'open'). Arguments are
!>    named VecA/NEFSetB (rather than Vec/NEFSet) to avoid name
!>    collisions with the reversed-operand
!>    MatMul(NEFSetA,VecB) => MQC_NEFSet_Vector_Multiply procedure.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] VecA
!>    \verbatim
!>        VecA is Type(MQC_Vector)
!>        The plain vector to contract against NEFSetB's I axis. Its
!>        length must equal NEFSetB's NumI.
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        The (I,J) grid to contract.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_vector_nefset_multiply(VecA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_vector),intent(in)::VecA
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
      type(mqc_nef_vector)::tmpVector
      type(mqc_nef_matrix)::tmpMatrix
      type(mqc_nef_r4)::tmpR4
      logical::first
!
      call NEFSetOut%init(1,NEFSetB%NumJ)
!
      select case(trim(NEFSetB%Type()))
      case('vector')
        do J = 1, NEFSetB%NumJ
          first = .true.
          do I = 1, NEFSetB%NumI
            if(NEFSetB%isZero(I,J)) cycle
            if(first) then
              tmpVector = VecA%at(I)*NEFSetB%atVector(I,J)
              first = .false.
            else
              tmpVector = tmpVector + VecA%at(I)*NEFSetB%atVector(I,J)
            end if
          end do
          if(.not.first) call NEFSetOut%put(tmpVector,1,J)
        end do
      case('matrix')
        do J = 1, NEFSetB%NumJ
          first = .true.
          do I = 1, NEFSetB%NumI
            if(NEFSetB%isZero(I,J)) cycle
            if(first) then
              tmpMatrix = VecA%at(I)*NEFSetB%atMatrix(I,J)
              first = .false.
            else
              tmpMatrix = tmpMatrix + VecA%at(I)*NEFSetB%atMatrix(I,J)
            end if
          end do
          if(.not.first) call NEFSetOut%put(tmpMatrix,1,J)
        end do
      case('r4')
        do J = 1, NEFSetB%NumJ
          first = .true.
          do I = 1, NEFSetB%NumI
            if(NEFSetB%isZero(I,J)) cycle
            if(first) then
              tmpR4 = VecA%at(I)*NEFSetB%atR4(I,J)
              first = .false.
            else
              tmpR4 = tmpR4 + VecA%at(I)*NEFSetB%atR4(I,J)
            end if
          end do
          if(.not.first) call NEFSetOut%put(tmpR4,1,J)
        end do
      case default
        call mqc_error_A('mqc_vector_nefset_multiply: NEF_Set is open (holds no elements)',&
             6,'Array_Type',NEFSetB%Type())
      end select
!
      end function mqc_vector_nefset_multiply
!
!
!     PROCEDURE MQC_NEFSet_Matrix_Multiply
!
!
!>    \brief <b> MQC_NEFSet_Matrix_Multiply contracts a MQC_NEF_Set's
!>    second (J) axis against a plain mqc_matrix. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFSet_Matrix_Multiply computes, for every I and every column K
!>    of MatB, NEFSetOut(I,K) = Sum_J NEFSetA(I,J)*MatB%at(J,K), where
!>    MatB is a plain mqc_matrix (not a MQC_NEF_Matrix). The result is a
!>    MQC_NEF_Set of the same kind as NEFSetA, with shape (NumI,
!>    size(MatB,2)) -- i.e. NEFSetA is treated as an (I,J) matrix of NEF
!>    elements and MatB as an ordinary (J,K) matrix, contracted over J
!>    exactly like a standard matrix product. MQC_NEFSet_Vector_Multiply
!>    is the special case where MatB has a single column. NEFSetA must
!>    already hold an element (not be 'open'). Arguments are named
!>    NEFSetA/MatB (rather than NEFSet/Mat) to avoid name collisions with
!>    the reversed-operand MatMul(MatA,NEFSet) procedure.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>        The (I,J) grid to contract.
!>    \endverbatim
!>
!>    \param[in] MatB
!>    \verbatim
!>        MatB is Type(MQC_Matrix)
!>        The plain matrix to contract against NEFSetA's J axis. Its
!>        number of rows must equal NEFSetA's NumJ.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_matrix_multiply(NEFSetA,MatB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_matrix),intent(in)::MatB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J,K,NumK
      type(mqc_nef_vector)::tmpVector
      type(mqc_nef_matrix)::tmpMatrix
      type(mqc_nef_r4)::tmpR4
      logical::first
!
      NumK = size(MatB,2)
      call NEFSetOut%init(NEFSetA%NumI,NumK)
!
      select case(trim(NEFSetA%Type()))
      case('vector')
        do I = 1, NEFSetA%NumI
          do K = 1, NumK
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J)) cycle
              if(first) then
                tmpVector = NEFSetA%atVector(I,J)*MatB%at(J,K)
                first = .false.
              else
                tmpVector = tmpVector + NEFSetA%atVector(I,J)*MatB%at(J,K)
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpVector,I,K)
          end do
        end do
      case('matrix')
        do I = 1, NEFSetA%NumI
          do K = 1, NumK
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J)) cycle
              if(first) then
                tmpMatrix = NEFSetA%atMatrix(I,J)*MatB%at(J,K)
                first = .false.
              else
                tmpMatrix = tmpMatrix + NEFSetA%atMatrix(I,J)*MatB%at(J,K)
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpMatrix,I,K)
          end do
        end do
      case('r4')
        do I = 1, NEFSetA%NumI
          do K = 1, NumK
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J)) cycle
              if(first) then
                tmpR4 = NEFSetA%atR4(I,J)*MatB%at(J,K)
                first = .false.
              else
                tmpR4 = tmpR4 + NEFSetA%atR4(I,J)*MatB%at(J,K)
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpR4,I,K)
          end do
        end do
      case default
        call mqc_error_A('mqc_nefset_matrix_multiply: NEF_Set is open (holds no elements)',&
             6,'Array_Type',NEFSetA%Type())
      end select
!
      end function mqc_nefset_matrix_multiply
!
!
!     PROCEDURE MQC_Matrix_NEFSet_Multiply
!
!
!>    \brief <b> MQC_Matrix_NEFSet_Multiply contracts a MQC_NEF_Set's
!>    first (I) axis against a plain mqc_matrix. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Matrix_NEFSet_Multiply computes, for every row K of MatA and
!>    every J, NEFSetOut(K,J) = Sum_I MatA%at(K,I)*NEFSetB(I,J), where
!>    MatA is a plain mqc_matrix (not a MQC_NEF_Matrix). The result is a
!>    MQC_NEF_Set of the same kind as NEFSetB, with shape (size(MatA,1),
!>    NumJ) -- i.e. MatA is an ordinary (K,I) matrix and NEFSetB is
!>    treated as an (I,J) matrix of NEF elements, contracted over I
!>    exactly like a standard matrix product. MQC_Vector_NEFSet_Multiply
!>    is the special case where MatA has a single row. NEFSetB must
!>    already hold an element (not be 'open'). Arguments are named
!>    MatA/NEFSetB (rather than Mat/NEFSet) to avoid name collisions with
!>    the reversed-operand MatMul(NEFSetA,MatB) procedure.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] MatA
!>    \verbatim
!>        MatA is Type(MQC_Matrix)
!>        The plain matrix to contract against NEFSetB's I axis. Its
!>        number of columns must equal NEFSetB's NumI.
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        The (I,J) grid to contract.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_matrix_nefset_multiply(MatA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_matrix),intent(in)::MatA
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J,K,NumK
      type(mqc_nef_vector)::tmpVector
      type(mqc_nef_matrix)::tmpMatrix
      type(mqc_nef_r4)::tmpR4
      logical::first
!
      NumK = size(MatA,1)
      call NEFSetOut%init(NumK,NEFSetB%NumJ)
!
      select case(trim(NEFSetB%Type()))
      case('vector')
        do J = 1, NEFSetB%NumJ
          do K = 1, NumK
            first = .true.
            do I = 1, NEFSetB%NumI
              if(NEFSetB%isZero(I,J)) cycle
              if(first) then
                tmpVector = MatA%at(K,I)*NEFSetB%atVector(I,J)
                first = .false.
              else
                tmpVector = tmpVector + MatA%at(K,I)*NEFSetB%atVector(I,J)
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpVector,K,J)
          end do
        end do
      case('matrix')
        do J = 1, NEFSetB%NumJ
          do K = 1, NumK
            first = .true.
            do I = 1, NEFSetB%NumI
              if(NEFSetB%isZero(I,J)) cycle
              if(first) then
                tmpMatrix = MatA%at(K,I)*NEFSetB%atMatrix(I,J)
                first = .false.
              else
                tmpMatrix = tmpMatrix + MatA%at(K,I)*NEFSetB%atMatrix(I,J)
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpMatrix,K,J)
          end do
        end do
      case('r4')
        do J = 1, NEFSetB%NumJ
          do K = 1, NumK
            first = .true.
            do I = 1, NEFSetB%NumI
              if(NEFSetB%isZero(I,J)) cycle
              if(first) then
                tmpR4 = MatA%at(K,I)*NEFSetB%atR4(I,J)
                first = .false.
              else
                tmpR4 = tmpR4 + MatA%at(K,I)*NEFSetB%atR4(I,J)
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpR4,K,J)
          end do
        end do
      case default
        call mqc_error_A('mqc_matrix_nefset_multiply: NEF_Set is open (holds no elements)',&
             6,'Array_Type',NEFSetB%Type())
      end select
!
      end function mqc_matrix_nefset_multiply
!
!
!     PROCEDURE MQC_NEFSet_NEFMatrix_Multiply
!
!
!>    \brief <b> MQC_NEFSet_NEFMatrix_Multiply matrix-multiplies every
!>    (I,J) entry of a MQC_NEF_Set by a MQC_NEF_Matrix for each (I,J).</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFSet_NEFMatrix_Multiply computes
!>    NEFSetOut(I,J) = MatMul(NEFSetA(I,J),NEFMatB), where NEFMatB
!>    is a single MQC_NEF_Matrix shared across every grid cell.
!>    NEFSetA must hold matrix elements, and every cell's column 
!>    space must be compatible with NEFMatB's row space 
!>    (guaranteed uniform across cells by MQC_NEF_Set's own put-time 
!>    dimension check).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>        The (I,J) grid of MQC_NEF_Matrix elements to contract.
!>    \endverbatim
!>
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is Type(MQC_NEF_Matrix)
!>        The matrix multiplied into every (I,J) cell before summing.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_nefmatrix_multiply(NEFSetA,NEFMatB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(trim(NEFSetA%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefset_nefmatrix_multiply: NEF_Set must hold matrix elements',&
           6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ)
!
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) cycle
          call NEFSetOut%put(MatMul(NEFSetA%atMatrix(I,J),NEFMatB),I,J)
        end do
      end do
!
      end function mqc_nefset_nefmatrix_multiply
!
!
!     PROCEDURE MQC_NEFMatrix_NEFSet_Multiply
!
!
!>    \brief <b> MQC_NEFMatrix_NEFSet_Multiply matrix-multiplies a
!>    (I,J) entry of a MQC_NEF_Set by a MQC_NEF_Matrix for each (I,J).</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_NEFSet_Multiply computes
!>    NEFSetOut(I,J) = MatMul(NEFMatA,NEFSetB(I,J)), where NEFMatA
!>    is a single MQC_NEF_Matrix shared across every grid cell.
!>    NEFSetB must hold matrix elements, and every cell's column 
!>    space must be compatible with NEFMatB's row space 
!>    (guaranteed uniform across cells by MQC_NEF_Set's own put-time 
!>    dimension check).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is Type(MQC_NEF_Matrix)
!>        The matrix multiplied into every (I,J) cell before summing.
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        The (I,J) grid of MQC_NEF_Matrix elements to contract.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_nefset_multiply(NEFMatA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(trim(NEFSetB%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefmatrix_nefset_multiply: NEF_Set must hold matrix elements',&
           6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ)
!
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) cycle
          call NEFSetOut%put(MatMul(NEFMatA,NEFSetB%atMatrix(I,J)),I,J)
        end do
      end do
!
      end function mqc_nefmatrix_nefset_multiply
!
!
!     PROCEDURE MQC_NEFSet_NEFSet_Multiply
!
!
!>    \brief <b> MQC_NEFSet_NEFSet_Multiply matrix-multiplies two
!>    MQC_NEF_Set "matrices of objects" at the grid level, combining
!>    contracted-index pairs via an outer product that raises rank. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFSet_NEFSet_Multiply treats NEFSetA (shape I x J) and
!>    NEFSetB (shape J x K) as ordinary matrices whose "elements" are
!>    NEF objects, and computes the standard matrix product at the grid
!>    level:
!>      NEFOut(I,K) = Sum_J outer(NEFSetA(I,J),NEFSetB(J,K))
!>    where "outer" is the per-cell outer product appropriate to
!>    NEFSetA/NEFSetB's shared kind:
!>      'vector' => MQC_NEF_Vector%outer2Matrix, NEFOut holds Matrix
!>                  elements
!>      'matrix' => MQC_NEF_Matrix%outer2R4, NEFOut holds R4 elements
!>    The result NEFOut has shape (NEFSetA%NumI, NEFSetB%NumJ), exactly
!>    like an ordinary (I,J)x(J,K)=(I,K) matrix product -- only the
!>    shared J index is contracted (summed); I and K remain as the
!>    output's own grid axes. NEFSetA and NEFSetB must hold the same
!>    kind ('vector' or 'matrix'; 'r4' is not supported, as there is no
!>    rank-6 NEF type to hold the per-cell outer product), and
!>    NEFSetA%NumJ must equal NEFSetB%NumI.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>        Shape (I,J).
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        Shape (J,K). NEFSetB%NumI must equal NEFSetA%NumJ.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_nefset_multiply(NEFSetA,NEFSetB) result(NEFOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA,NEFSetB
      type(mqc_nef_set)::NEFOut
      type(mqc_nef_matrix)::tmpMatrix
      type(mqc_nef_r4)::tmpR4
      integer::I,J,K
      logical::first
!
      if(trim(NEFSetA%Type()).ne.trim(NEFSetB%Type())) call mqc_error_A(&
           'mqc_nefset_nefset_multiply: NEF_Sets must hold the same kind',&
           6,'NEFSetA%Type',NEFSetA%Type())
      if(NEFSetA%NumJ.ne.NEFSetB%NumI) call mqc_error_i(&
           'mqc_nefset_nefset_multiply: NEFSetA%NumJ must equal NEFSetB%NumI',&
           6,'NEFSetA%NumJ',int(NEFSetA%NumJ,kind=int64),&
           'NEFSetB%NumI',int(NEFSetB%NumI,kind=int64))
!
      call NEFOut%init(NEFSetA%NumI,NEFSetB%NumJ)
!
      select case(trim(NEFSetA%Type()))
      case('vector')
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                tmpMatrix = mqc_nefvector_nefvector_outer_product(&
                     NEFSetA%atVector(I,J),NEFSetB%atVector(J,K))
                first = .false.
              else
                tmpMatrix = tmpMatrix + mqc_nefvector_nefvector_outer_product(&
                     NEFSetA%atVector(I,J),NEFSetB%atVector(J,K))
              end if
            end do
            if(.not.first) call NEFOut%put(tmpMatrix,I,K)
          end do
        end do
      case('matrix')
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                tmpR4 = mqc_nefmatrix_nefmatrix_outer_product(&
                     NEFSetA%atMatrix(I,J),NEFSetB%atMatrix(J,K))
                first = .false.
              else
                tmpR4 = tmpR4 + mqc_nefmatrix_nefmatrix_outer_product(&
                     NEFSetA%atMatrix(I,J),NEFSetB%atMatrix(J,K))
              end if
            end do
            if(.not.first) call NEFOut%put(tmpR4,I,K)
          end do
        end do
      case default
        call mqc_error_A('mqc_nefset_nefset_multiply: unsupported NEF_Set kind '//&
             '(must be vector or matrix; r4 has no rank-8 outer-product target)',&
             6,'Array_Type',NEFSetA%Type())
      end select
!
      end function mqc_nefset_nefset_multiply
!
!
!     PROCEDURE MQC_NEF_R4_Apply_Order
!
!
!>    \brief <b> MQC_NEF_R4_Apply_Order permutes a MQC_NEF_R4's four axes
!>    into a caller-chosen order. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Apply_Order returns R4in with its four axes permuted so
!>    that physical axis p (p=1..4) of the output holds what was axis
!>    Order(p) of the input. If Order is absent, R4in is returned
!>    unchanged. The permutation is realized as at most 3 pairwise swaps
!>    via R4in%swap(dim1,dim2) (a standard permutation-via-transpositions
!>    decomposition): a "current" bookkeeping array tracks which original
!>    axis occupies each position, and for each position p=1..3 in turn,
!>    if the wrong axis sits there, the position q>p holding the axis
!>    that belongs at p is swapped into place. This is a private helper,
!>    not part of any public generic interface.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] R4in
!>    \verbatim
!>        R4in is type(mqc_nef_r4)
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Order is integer(kind=int64),dimension(4),optional
!>        A permutation of [1,2,3,4]. Absent leaves R4in unchanged.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_apply_order(R4in,Order) result(R4out)
!
      implicit none
      type(mqc_nef_r4),intent(in)::R4in
      integer(kind=int64),dimension(4),optional,intent(in)::Order
      type(mqc_nef_r4)::R4out
      integer(kind=int64),dimension(4)::myOrder,current
      integer(kind=int64)::p,q,tmp
!
      R4out = R4in
      if(.not.present(Order)) return
      myOrder = Order
      current = [1,2,3,4]
      do p = 1, 3
        if(current(p).ne.myOrder(p)) then
          do q = p+1, 4
            if(current(q).eq.myOrder(p)) exit
          end do
          R4out = R4out%swap(p,q)
          tmp = current(p) ; current(p) = current(q) ; current(q) = tmp
        end if
      end do
!
      end function mqc_nef_r4_apply_order
!
!
!     PROCEDURE MQC_NEFSet_NEFSet_Contraction
!
!
!>    \brief <b> MQC_NEFSet_NEFSet_Contraction contracts two MQC_NEF_Set
!>    "matrices of objects" at the grid level, reducing each contracted
!>    cell pair to a scalar. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFSet_NEFSet_Contraction treats NEFSetA (shape I x J) and
!>    NEFSetB (shape J x K) as ordinary matrices whose "elements" are NEF
!>    objects, and computes
!>      MatOut(I,K) = Sum_J Contraction(NEFSetA(I,J),NEFSetB(J,K))
!>    exactly like an ordinary (I,J)x(J,K)=(I,K) matrix product -- only
!>    the shared J index is contracted (summed); I and K remain as
!>    MatOut's own two axes. Since the per-cell Contraction already
!>    reduces each pair to a scalar, MatOut is a plain mqc_matrix (not a
!>    MQC_NEF_Set). NEFSetA and NEFSetB must hold the same kind ('vector',
!>    'matrix', or 'r4'), and NEFSetA%NumJ must equal NEFSetB%NumI.
!>
!>    For matrix-kind cells, transposeB (default false) contracts against
!>    Transpose(NEFSetB(J,K)) instead of NEFSetB(J,K) directly -- the
!>    (I,J,u,v)x(J,K,v,u)->(I,K) variant. For r4-kind cells, label is
!>    passed through to the underlying MQC_NEF_R4,MQC_NEF_R4 Contraction
!>    ('coulomb'/'exchange'/'doublebar'), defaulting here to 'coulomb'
!>    (NOT that procedure's own 'doublebar' default) since the exchange
!>    diagram imposes extra axis-conformability requirements (it swaps
!>    the role of two axes internally) that a caller-chosen Order may
!>    not satisfy; pass label='doublebar' explicitly if exchange is
!>    wanted and the axes are known to support it. Order (if present)
!>    permutes NEFSetB(J,K)'s four axes via MQC_NEF_R4_Apply_Order
!>    before contracting, realizing a contraction against a
!>    caller-chosen axis order.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>        Shape (I,J).
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        Shape (J,K). NEFSetB%NumI must equal NEFSetA%NumJ.
!>    \endverbatim
!>
!>    \param[in] transposeB
!>    \verbatim
!>        transposeB is logical,optional
!>        Matrix-kind only. Default false.
!>    \endverbatim
!>
!>    \param[in] label
!>    \verbatim
!>        label is character(len=*),optional
!>        R4-kind only, passed through to the per-cell Contraction.
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Order is integer(kind=int64),dimension(4),optional
!>        R4-kind only. Permutes each NEFSetB cell's axes before
!>        contracting.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_nefset_contraction(NEFSetA,NEFSetB,transposeB,label,Order) &
           result(MatOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA,NEFSetB
      logical,optional,intent(in)::transposeB
      character(len=*),optional,intent(in)::label
      integer(kind=int64),dimension(4),optional,intent(in)::Order
      type(mqc_matrix)::MatOut
      logical::myTransposeB
      character(len=64)::myLabel
      type(mqc_scalar)::tmpScalar
      integer::I,J,K
      logical::first
!
      if(trim(NEFSetA%Type()).ne.trim(NEFSetB%Type())) call mqc_error_A(&
           'mqc_nefset_nefset_contraction: NEF_Sets must hold the same kind',&
           6,'NEFSetA%Type',NEFSetA%Type())
      if(NEFSetA%NumJ.ne.NEFSetB%NumI) call mqc_error_i(&
           'mqc_nefset_nefset_contraction: NEFSetA%NumJ must equal NEFSetB%NumI',&
           6,'NEFSetA%NumJ',int(NEFSetA%NumJ,kind=int64),&
           'NEFSetB%NumI',int(NEFSetB%NumI,kind=int64))
!
      myTransposeB = .false.
      if(present(transposeB)) myTransposeB = transposeB
      myLabel = 'coulomb'
      if(present(label)) myLabel = label
!
      call MatOut%init(int(NEFSetA%NumI,kind=int64),int(NEFSetB%NumJ,kind=int64))
!
      select case(trim(NEFSetA%Type()))
      case('vector')
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                tmpScalar = Contraction(NEFSetA%atVector(I,J),NEFSetB%atVector(J,K))
                first = .false.
              else
                tmpScalar = tmpScalar + Contraction(NEFSetA%atVector(I,J),NEFSetB%atVector(J,K))
              end if
            end do
            if(.not.first) call MatOut%put(tmpScalar,I,K)
          end do
        end do
      case('matrix')
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                if(myTransposeB) then
                  tmpScalar = Contraction(NEFSetA%atMatrix(I,J),Transpose(NEFSetB%atMatrix(J,K)))
                else
                  tmpScalar = Contraction(NEFSetA%atMatrix(I,J),NEFSetB%atMatrix(J,K))
                end if
                first = .false.
              else
                if(myTransposeB) then
                  tmpScalar = tmpScalar + &
                       Contraction(NEFSetA%atMatrix(I,J),Transpose(NEFSetB%atMatrix(J,K)))
                else
                  tmpScalar = tmpScalar + Contraction(NEFSetA%atMatrix(I,J),NEFSetB%atMatrix(J,K))
                end if
              end if
            end do
            if(.not.first) call MatOut%put(tmpScalar,I,K)
          end do
        end do
      case('r4')
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                tmpScalar = Contraction(NEFSetA%atR4(I,J),&
                     mqc_nef_r4_apply_order(NEFSetB%atR4(J,K),Order),myLabel)
                first = .false.
              else
                tmpScalar = tmpScalar + Contraction(NEFSetA%atR4(I,J),&
                     mqc_nef_r4_apply_order(NEFSetB%atR4(J,K),Order),myLabel)
              end if
            end do
            if(.not.first) call MatOut%put(tmpScalar,I,K)
          end do
        end do
      case default
        call mqc_error_A('mqc_nefset_nefset_contraction: NEF_Set is open '//&
             '(holds no elements)',6,'Array_Type',NEFSetA%Type())
      end select
!
      end function mqc_nefset_nefset_contraction
!
!
!     PROCEDURE MQC_NEFSet_NEFSet_Full_Contraction
!
!
!>    \brief <b> MQC_NEFSet_NEFSet_Full_Contraction fully contracts two
!>    matching MQC_NEF_Set objects down to a single scalar. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFSet_NEFSet_Full_Contraction requires NEFSetA and NEFSetB to
!>    share the exact same (NumI,NumJ) grid shape (mirroring
!>    Contraction's own "exact shape, fully reduce" philosophy one level
!>    up onto the grid), and computes
!>      output = Sum_I Sum_J Contraction(NEFSetA(I,J),NEFSetB(I,J))
!>    Unlike MQC_NEFSet_NEFSet_Contraction (which keeps the outer (I,K)
!>    grid, contracting only the shared J axis), this contracts every
!>    index -- both outer grid axes and every inner axis of each cell --
!>    down to a single mqc_scalar. NEFSetA and NEFSetB must hold the same
!>    kind ('vector', 'matrix', or 'r4').
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        Must share NEFSetA's exact (NumI,NumJ) grid shape.
!>    \endverbatim
!>
!>    \param[in] label
!>    \verbatim
!>        label is character(len=*),optional
!>        R4-kind only, passed through to the per-cell Contraction.
!>        Defaults to 'coulomb' here (not that procedure's own
!>        'doublebar' default) since the exchange diagram imposes
!>        extra axis-conformability requirements.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_nefset_full_contraction(NEFSetA,NEFSetB,label) result(output)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA,NEFSetB
      character(len=*),optional,intent(in)::label
      character(len=64)::myLabel
      type(mqc_scalar)::output
      integer::I,J
!
      if(trim(NEFSetA%Type()).ne.trim(NEFSetB%Type())) call mqc_error_A(&
           'mqc_nefset_nefset_full_contraction: NEF_Sets must hold the same kind',&
           6,'NEFSetA%Type',NEFSetA%Type())
      if(NEFSetA%NumI.ne.NEFSetB%NumI.or.NEFSetA%NumJ.ne.NEFSetB%NumJ) call mqc_error_i(&
           'mqc_nefset_nefset_full_contraction: NEFSetA and NEFSetB must have '//&
           'identical grid shape',6,&
           'NEFSetA%NumI',int(NEFSetA%NumI,kind=int64),&
           'NEFSetB%NumI',int(NEFSetB%NumI,kind=int64))
!
      myLabel = 'coulomb'
      if(present(label)) myLabel = label
      output = 0.0
!
      select case(trim(NEFSetA%Type()))
      case('vector')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(I,J)) cycle
            output = output + Contraction(NEFSetA%atVector(I,J),NEFSetB%atVector(I,J))
          end do
        end do
      case('matrix')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(I,J)) cycle
            output = output + Contraction(NEFSetA%atMatrix(I,J),NEFSetB%atMatrix(I,J))
          end do
        end do
      case('r4')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(I,J)) cycle
            output = output + Contraction(NEFSetA%atR4(I,J),NEFSetB%atR4(I,J),myLabel)
          end do
        end do
      case default
        call mqc_error_A('mqc_nefset_nefset_full_contraction: NEF_Set is open '//&
             '(holds no elements)',6,'Array_Type',NEFSetA%Type())
      end select
!
      end function mqc_nefset_nefset_full_contraction
!
!
!     PROCEDURE MQC_NEFSet_NEFSet_Partial_Contraction
!
!
!>    \brief <b> MQC_NEFSet_NEFSet_Partial_Contraction partially
!>    contracts two MQC_NEF_Set objects at the grid level, reducing each
!>    contracted cell pair's rank by the shared inner axes. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFSet_NEFSet_Partial_Contraction treats NEFSetA (shape I x J)
!>    and NEFSetB (shape J x K) as ordinary matrices whose "elements" are
!>    NEF objects, and computes
!>      NEFSetOut(I,K) = Sum_J <cell-op>(NEFSetA(I,J),NEFSetB(J,K))
!>    exactly like an ordinary (I,J)x(J,K)=(I,K) matrix product -- only
!>    the shared J index is contracted (summed); I and K remain as
!>    NEFSetOut's own two grid axes. <cell-op> and NEFSetOut's resulting
!>    kind depend on NEFSetA/NEFSetB's own kinds:
!>
!>      'r4','r4' => the underlying MQC_NEF_R4,MQC_NEF_R4
!>        PartialContraction, using index1/index2 (role convention:
!>        contracted axes first, surviving axes last -- see that
!>        procedure's own documentation). NEFSetOut holds Matrix.
!>
!>      'matrix','matrix' => a combined-index contraction using ONLY
!>        index1=[a,b,c,d], a permutation of [1,2,3,4] where 1,2 address
!>        MatA's (row,col) and 3,4 address MatB's (row,col). [a,b] is the
!>        contracted pair and MUST come from distinct matrices (one of
!>        1,2 and one of 3,4) -- checked, errors otherwise. [c,d] is the
!>        surviving pair (one index left over from each matrix) IN
!>        OUTPUT ORDER: c becomes NEFSetOut's row-contributing index, d
!>        its col-contributing index (regardless of which matrix each
!>        came from). Implemented as ordinary MatMul(MatA,MatB) with
!>        Transpose applied to whichever operand's contracted axis isn't
!>        already where MatMul itself contracts (its col for the left
!>        operand, its row for the right), plus a final Transpose of the
!>        result if [c,d] asks for the swapped output order -- reusing
!>        the already-correct bare-object MatMul/Transpose rather than a
!>        new primitive. Default index1=[2,3,1,4] (contract MatA's col
!>        against MatB's row, survivors in (MatA-row,MatB-col) order)
!>        reduces exactly to ordinary MatMul(MatA,MatB). NEFSetOut holds
!>        Matrix. index2 is ignored for this kind combination.
!>
!>      'matrix','vector' => index1(1:3)=[a,b,c], a permutation of
!>        [1,2,3] where 1,2 address MatA's (row,col) and 3 addresses
!>        VecB's single axis. [a,b] contracted, must be one of 1,2 and
!>        3 (VecB is always fully consumed, contributing nothing to the
!>        output). c is MatA's surviving axis (the only one, so no
!>        output-order choice is needed). Default index1=[2,3,1]
!>        (contract MatA's col against VecB) reduces exactly to ordinary
!>        MatMul(MatA,VecB). NEFSetOut holds Vector.
!>
!>      'vector','matrix' => index1(1:3)=[a,b,c], 1 addresses VecA's
!>        single axis, 2,3 address MatB's (row,col). [a,b] contracted,
!>        must be 1 and one of 2,3 (VecA always fully consumed). c is
!>        MatB's surviving axis. Default index1=[1,3,2] (contract VecA
!>        against MatB's col) reduces exactly to ordinary
!>        MatMul(MatB,VecA). NEFSetOut holds Vector.
!>
!>      'r4','matrix' => the underlying MQC_NEF_R4,MQC_NEF_Matrix
!>        PartialContraction (same procedure as MQC_NEF_R4_NEF_Matrix_
!>        Contraction, registered under both Contraction and
!>        PartialContraction), using index1/label (index1 plays the same
!>        role as that procedure's own indexA: index1's first two
!>        entries are the R4 axes contracted against NEFSetB's row/col,
!>        the last two survive as the per-cell output Matrix's row/col;
!>        label selects which coulomb/exchange-type diagrams are
!>        summed). Default index1=[1,2,3,4]. NEFSetOut holds Matrix.
!>        index2 is ignored for this kind combination.
!>
!>    NEFSetA%NumJ must equal NEFSetB%NumI in every case.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>        Shape (I,J).
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        Shape (J,K). NEFSetB%NumI must equal NEFSetA%NumJ.
!>    \endverbatim
!>
!>    \param[in] index1
!>    \verbatim
!>        index1 is integer(kind=int64),dimension(4),optional
!>        r4,r4: role array passed through to the per-cell
!>        PartialContraction. matrix,matrix: combined [a,b,c,d] (see
!>        above). matrix,vector / vector,matrix: combined [a,b,c] in
!>        index1(1:3) (index1(4) unused). r4,matrix: plays the role of
!>        MQC_NEF_R4_NEF_Matrix_Contraction's own indexA.
!>    \endverbatim
!>
!>    \param[in] index2
!>    \verbatim
!>        index2 is integer(kind=int64),dimension(4),optional
!>        r4,r4 only, passed through to the per-cell PartialContraction.
!>        Ignored for every other kind combination.
!>    \endverbatim
!>
!>    \param[in] label
!>    \verbatim
!>        label is character(len=*),optional
!>        r4,matrix only. Passed through to the per-cell
!>        MQC_NEF_R4_NEF_Matrix_Contraction ('coulomb'/'exchange'/
!>        'doublebar'), defaulting here to 'coulomb' (not that
!>        procedure's own 'doublebar' default) since the exchange
!>        diagram imposes extra axis-conformability requirements (it
!>        swaps the role of two axes internally) that a caller-chosen
!>        index1 may not satisfy; pass label='doublebar' explicitly if
!>        exchange is wanted and the axes are known to support it.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_nefset_partial_contraction(NEFSetA,NEFSetB,index1,index2,&
           label) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA,NEFSetB
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      character(len=*),optional,intent(in)::label
      character(len=64)::myLabel
      type(mqc_nef_set)::NEFSetOut
      type(mqc_nef_matrix)::tmpMatrix
      type(mqc_nef_vector)::tmpVector
      integer(kind=int64)::a,b,c,d,contractedA,contractedB
      logical::myTransposeA,myTransposeB,swapOut,first
      integer::I,J,K
!
      if(NEFSetA%NumJ.ne.NEFSetB%NumI) call mqc_error_i(&
           'mqc_nefset_nefset_partial_contraction: NEFSetA%NumJ must equal NEFSetB%NumI',&
           6,'NEFSetA%NumJ',int(NEFSetA%NumJ,kind=int64),&
           'NEFSetB%NumI',int(NEFSetB%NumI,kind=int64))
!
      myLabel = 'coulomb'
      if(present(label)) myLabel = label
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetB%NumJ)
!
      select case(trim(NEFSetA%Type())//'_'//trim(NEFSetB%Type()))
      case('r4_r4')
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                tmpMatrix = PartialContraction(NEFSetA%atR4(I,J),NEFSetB%atR4(J,K),index1,index2)
                first = .false.
              else
                tmpMatrix = tmpMatrix + &
                     PartialContraction(NEFSetA%atR4(I,J),NEFSetB%atR4(J,K),index1,index2)
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpMatrix,I,K)
          end do
        end do
      case('matrix_matrix')
        if(present(index1)) then
          a=index1(1) ; b=index1(2) ; c=index1(3) ; d=index1(4)
        else
          a=2 ; b=3 ; c=1 ; d=4
        end if
        if(a.lt.1.or.a.gt.4.or.b.lt.1.or.b.gt.4.or.c.lt.1.or.c.gt.4.or.d.lt.1.or.d.gt.4.or.&
           a.eq.b.or.a.eq.c.or.a.eq.d.or.b.eq.c.or.b.eq.d.or.c.eq.d) call mqc_error_i(&
             'mqc_nefset_nefset_partial_contraction: index1 must be a permutation '//&
             'of [1,2,3,4] for matrix,matrix',6,'a',a,'b',b,'c',c,'d',d)
        if((a.le.2.and.b.le.2).or.(a.gt.2.and.b.gt.2)) call mqc_error_i(&
             'mqc_nefset_nefset_partial_contraction: index1(1:2) (the contracted '//&
             'pair) must come from distinct matrices (one of 1,2 and one of 3,4)',&
             6,'a',a,'b',b)
        if(a.le.2) then
          contractedA = a ; contractedB = b
        else
          contractedA = b ; contractedB = a
        end if
        if(.not.((c.eq.(3-contractedA).and.d.eq.(7-contractedB)).or.&
                 (c.eq.(7-contractedB).and.d.eq.(3-contractedA)))) call mqc_error_i(&
             'mqc_nefset_nefset_partial_contraction: index1(3:4) must be exactly '//&
             'the surviving indexes',6,'c',c,'d',d)
        swapOut = (c.eq.(7-contractedB))
        myTransposeA = (contractedA.eq.1)
        myTransposeB = (contractedB.eq.4)
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                if(myTransposeA.and.myTransposeB) then
                  tmpMatrix = MatMul(Transpose(NEFSetA%atMatrix(I,J)),Transpose(NEFSetB%atMatrix(J,K)))
                else if(myTransposeA) then
                  tmpMatrix = MatMul(Transpose(NEFSetA%atMatrix(I,J)),NEFSetB%atMatrix(J,K))
                else if(myTransposeB) then
                  tmpMatrix = MatMul(NEFSetA%atMatrix(I,J),Transpose(NEFSetB%atMatrix(J,K)))
                else
                  tmpMatrix = MatMul(NEFSetA%atMatrix(I,J),NEFSetB%atMatrix(J,K))
                end if
                first = .false.
              else
                if(myTransposeA.and.myTransposeB) then
                  tmpMatrix = tmpMatrix + &
                       MatMul(Transpose(NEFSetA%atMatrix(I,J)),Transpose(NEFSetB%atMatrix(J,K)))
                else if(myTransposeA) then
                  tmpMatrix = tmpMatrix + &
                       MatMul(Transpose(NEFSetA%atMatrix(I,J)),NEFSetB%atMatrix(J,K))
                else if(myTransposeB) then
                  tmpMatrix = tmpMatrix + &
                       MatMul(NEFSetA%atMatrix(I,J),Transpose(NEFSetB%atMatrix(J,K)))
                else
                  tmpMatrix = tmpMatrix + MatMul(NEFSetA%atMatrix(I,J),NEFSetB%atMatrix(J,K))
                end if
              end if
            end do
            if(.not.first) then
              if(swapOut) tmpMatrix = Transpose(tmpMatrix)
              call NEFSetOut%put(tmpMatrix,I,K)
            end if
          end do
        end do
      case('matrix_vector')
        if(present(index1)) then
          a=index1(1) ; b=index1(2) ; c=index1(3)
        else
          a=2 ; b=3 ; c=1
        end if
        if(a.lt.1.or.a.gt.3.or.b.lt.1.or.b.gt.3.or.c.lt.1.or.c.gt.3.or.&
           a.eq.b.or.a.eq.c.or.b.eq.c) call mqc_error_i(&
             'mqc_nefset_nefset_partial_contraction: index1(1:3) must be a '//&
             'permutation of [1,2,3] for matrix,vector',6,'a',a,'b',b,'c',c)
        if((a.le.2.and.b.le.2).or.(a.eq.3.and.b.eq.3)) call mqc_error_i(&
             'mqc_nefset_nefset_partial_contraction: index1(1:2) (the contracted '//&
             'pair) must come from distinct objects (one of 1,2 and 3)',6,'a',a,'b',b)
        if(a.le.2) then
          contractedA = a
        else
          contractedA = b
        end if
        if(c.ne.(3-contractedA)) call mqc_error_i(&
             'mqc_nefset_nefset_partial_contraction: index1(3) must be the '//&
             'surviving MatA index',6,'c',c)
        myTransposeA = (contractedA.eq.1)
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                if(myTransposeA) then
                  tmpVector = MatMul(Transpose(NEFSetA%atMatrix(I,J)),NEFSetB%atVector(J,K))
                else
                  tmpVector = MatMul(NEFSetA%atMatrix(I,J),NEFSetB%atVector(J,K))
                end if
                first = .false.
              else
                if(myTransposeA) then
                  tmpVector = tmpVector + &
                       MatMul(Transpose(NEFSetA%atMatrix(I,J)),NEFSetB%atVector(J,K))
                else
                  tmpVector = tmpVector + MatMul(NEFSetA%atMatrix(I,J),NEFSetB%atVector(J,K))
                end if
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpVector,I,K)
          end do
        end do
      case('vector_matrix')
        if(present(index1)) then
          a=index1(1) ; b=index1(2) ; c=index1(3)
        else
          a=1 ; b=3 ; c=2
        end if
        if(a.lt.1.or.a.gt.3.or.b.lt.1.or.b.gt.3.or.c.lt.1.or.c.gt.3.or.&
           a.eq.b.or.a.eq.c.or.b.eq.c) call mqc_error_i(&
             'mqc_nefset_nefset_partial_contraction: index1(1:3) must be a '//&
             'permutation of [1,2,3] for vector,matrix',6,'a',a,'b',b,'c',c)
        if(.not.((a.eq.1.and.b.ge.2).or.(b.eq.1.and.a.ge.2))) call mqc_error_i(&
             'mqc_nefset_nefset_partial_contraction: index1(1:2) (the contracted '//&
             'pair) must come from distinct objects (1 and one of 2,3)',6,'a',a,'b',b)
        if(a.eq.1) then
          contractedB = b
        else
          contractedB = a
        end if
        if(c.ne.(5-contractedB)) call mqc_error_i(&
             'mqc_nefset_nefset_partial_contraction: index1(3) must be the '//&
             'surviving MatB index',6,'c',c)
        myTransposeB = (contractedB.eq.2)
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                if(myTransposeB) then
                  tmpVector = MatMul(Transpose(NEFSetB%atMatrix(J,K)),NEFSetA%atVector(I,J))
                else
                  tmpVector = MatMul(NEFSetB%atMatrix(J,K),NEFSetA%atVector(I,J))
                end if
                first = .false.
              else
                if(myTransposeB) then
                  tmpVector = tmpVector + &
                       MatMul(Transpose(NEFSetB%atMatrix(J,K)),NEFSetA%atVector(I,J))
                else
                  tmpVector = tmpVector + MatMul(NEFSetB%atMatrix(J,K),NEFSetA%atVector(I,J))
                end if
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpVector,I,K)
          end do
        end do
      case('r4_matrix')
        do I = 1, NEFSetA%NumI
          do K = 1, NEFSetB%NumJ
            first = .true.
            do J = 1, NEFSetA%NumJ
              if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
              if(first) then
                tmpMatrix = PartialContraction(NEFSetA%atR4(I,J),NEFSetB%atMatrix(J,K),myLabel,index1)
                first = .false.
              else
                tmpMatrix = tmpMatrix + &
                     PartialContraction(NEFSetA%atR4(I,J),NEFSetB%atMatrix(J,K),myLabel,index1)
              end if
            end do
            if(.not.first) call NEFSetOut%put(tmpMatrix,I,K)
          end do
        end do
      case default
        call mqc_error_A('mqc_nefset_nefset_partial_contraction: unsupported NEF_Set '//&
             'kind combination',6,'NEFSetA%Type',NEFSetA%Type())
      end select
!
      end function mqc_nefset_nefset_partial_contraction
!
!
!     PROCEDURE MQC_NEFSet_NEFSet_Partial_Contraction_2NEF_R4
!
!
!>    \brief <b> MQC_NEFSet_NEFSet_Partial_Contraction_2NEF_R4 partially
!>    contracts two r4-kind MQC_NEF_Set objects at the grid level,
!>    keeping each contracted cell pair's R4,R4 result at rank 4. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFSet_NEFSet_Partial_Contraction_2NEF_R4 treats NEFSetA (shape
!>    I x J) and NEFSetB (shape J x K) as ordinary matrices whose
!>    "elements" are MQC_NEF_R4 objects, and computes
!>      NEFSetOut(I,K) = Sum_J PartialContraction2R4(NEFSetA(I,J),&
!>                                                    NEFSetB(J,K),&
!>                                                    index1,index2)
!>    exactly like an ordinary (I,J)x(J,K)=(I,K) matrix product -- only
!>    the shared J index is contracted (summed); I and K remain as
!>    NEFSetOut's own two grid axes. index1/index2 carry the same role
!>    convention as the underlying MQC_NEF_R4,MQC_NEF_R4
!>    PartialContraction2R4 (2 of A's 4 axes and 2 of B's 4 axes
!>    contracted away pairwise, A's other 2 and B's other 2 surviving
!>    independently as the 4 axes of the per-cell R4 result; see that
!>    procedure's own documentation). Since the per-cell
!>    PartialContraction2R4 keeps each R4,R4 pair at rank 4, NEFSetOut
!>    holds r4-kind elements. NEFSetA and NEFSetB must both hold r4
!>    elements, and NEFSetA%NumJ must equal NEFSetB%NumI.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>        Shape (I,J). Must hold r4 elements.
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        Shape (J,K). Must hold r4 elements. NEFSetB%NumI must equal
!>        NEFSetA%NumJ.
!>    \endverbatim
!>
!>    \param[in] index1
!>    \verbatim
!>        index1 is integer(kind=int64),dimension(4),optional
!>        Passed through to the per-cell PartialContraction2R4.
!>    \endverbatim
!>
!>    \param[in] index2
!>    \verbatim
!>        index2 is integer(kind=int64),dimension(4),optional
!>        Passed through to the per-cell PartialContraction2R4.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_nefset_partial_contraction_2nef_r4(NEFSetA,NEFSetB,index1,index2) &
           result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA,NEFSetB
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      type(mqc_nef_set)::NEFSetOut
      type(mqc_nef_r4)::tmpR4
      integer::I,J,K
      logical::first
!
      if(trim(NEFSetA%Type()).ne.'r4'.or.trim(NEFSetB%Type()).ne.'r4') call mqc_error_A(&
           'mqc_nefset_nefset_partial_contraction_2nef_r4: both NEF_Sets must hold '//&
           'r4 elements',6,'NEFSetA%Type',NEFSetA%Type())
      if(NEFSetA%NumJ.ne.NEFSetB%NumI) call mqc_error_i(&
           'mqc_nefset_nefset_partial_contraction_2nef_r4: NEFSetA%NumJ must equal '//&
           'NEFSetB%NumI',6,'NEFSetA%NumJ',int(NEFSetA%NumJ,kind=int64),&
           'NEFSetB%NumI',int(NEFSetB%NumI,kind=int64))
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetB%NumJ)
!
      do I = 1, NEFSetA%NumI
        do K = 1, NEFSetB%NumJ
          first = .true.
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).or.NEFSetB%isZero(J,K)) cycle
            if(first) then
              tmpR4 = PartialContraction2R4(NEFSetA%atR4(I,J),NEFSetB%atR4(J,K),index1,index2)
              first = .false.
            else
              tmpR4 = tmpR4 + &
                   PartialContraction2R4(NEFSetA%atR4(I,J),NEFSetB%atR4(J,K),index1,index2)
            end if
          end do
          if(.not.first) call NEFSetOut%put(tmpR4,I,K)
        end do
      end do
!
      end function mqc_nefset_nefset_partial_contraction_2nef_r4
!
!
!     PROCEDURE MQC_NEF_Set_Expand
!
!
!>    \brief <b> MQC_NEF_Set_Expand is the counterpart of MQC_NEF_Set%
!>    collapse: it applies the underlying element's expand to every
!>    (I,J) entry of a MQC_NEF_Set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Expand returns a new MQC_NEF_Set of the same (I,J)
!>    shape as NEFSet, where every entry is the underlying element's
!>    expand (every cell shares the same target map(s), since
!>    MQC_NEF_Set%put already guarantees identical internal dimensions
!>    on every axis for every cell). If NEFSet holds Vector elements,
!>    expand raises rank 1 to 2 (the result set holds Matrix elements)
!>    and requires map_matrix_1 (row) and map_matrix_2 (col), matching
!>    MQC_NEF_Vector's own expand. If NEFSet holds Matrix elements,
!>    expand raises rank 2 to 4 (the result set holds R4 elements) and
!>    requires map_matrix_1..4, matching MQC_NEF_Matrix's own expand.
!>    NEFSet holding R4 elements (already the highest rank handled here)
!>    or being 'open' is an error.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The set whose elements are each expanded.
!>    \endverbatim
!>
!>    \param[in] map_matrix_1,map_matrix_2
!>    \verbatim
!>        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Optional
!>        Required (together) when NEFSet holds Vector elements (row and
!>        col, respectively); forwarded to the underlying Vector expand.
!>        Also used as the first two of four axes when NEFSet holds
!>        Matrix elements.
!>    \endverbatim
!>
!>    \param[in] map_matrix_3,map_matrix_4
!>    \verbatim
!>        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Optional
!>        Required (together with map_matrix_1,map_matrix_2) when NEFSet
!>        holds Matrix elements; forwarded to the underlying Matrix
!>        expand.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        Character(len=*),Optional
!>        Forwarded to the underlying Vector/Matrix expand.
!>    \endverbatim
!>
!>    \param[in] spinTypeKet
!>    \verbatim
!>        Character(len=*),Optional
!>        Forwarded to the underlying Matrix expand only.
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Integer(kind=int64),Dimension(4),Optional
!>        Forwarded to the underlying Matrix expand only.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_expand(NEFSet,map_matrix_1,map_matrix_2,&
           map_matrix_3,map_matrix_4,spinType,spinTypeKet,Order) result(NEFSetOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,intent(in),optional::&
           map_matrix_1,map_matrix_2,map_matrix_3,map_matrix_4
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      select case(trim(NEFSet%Type()))
      case('vector')
        if(.not.(present(map_matrix_1).and.present(map_matrix_2))) call mqc_error(&
             'mqc_nef_set_expand: vector-kind NEF_Set requires map_matrix_1 '//&
             'and map_matrix_2',6)
        call NEFSetOut%init(NEFSet%NumI,NEFSet%NumJ)
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(NEFSet%Set_Vector(I,J)%expand(&
                 map_matrix_1,map_matrix_2,spinType=spinType),I,J)
          end do
        end do
      case('matrix')
        if(.not.(present(map_matrix_1).and.present(map_matrix_2).and.&
                 present(map_matrix_3).and.present(map_matrix_4))) call mqc_error(&
             'mqc_nef_set_expand: matrix-kind NEF_Set requires map_matrix_1, '//&
             'map_matrix_2, map_matrix_3 and map_matrix_4',6)
        call NEFSetOut%init(NEFSet%NumI,NEFSet%NumJ)
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(NEFSet%Set_Matrix(I,J)%expand(&
                 map_matrix_1,map_matrix_2,map_matrix_3,map_matrix_4,&
                 spinType=spinType,spinTypeKet=spinTypeKet,Order=Order),I,J)
          end do
        end do
      case default
        call mqc_error_A('mqc_nef_set_expand: NEF_Set kind does not support '//&
             'expand (must be vector or matrix)',6,'Array_Type',NEFSet%Type())
      end select
!
      end function mqc_nef_set_expand
!
!
!     PROCEDURE MQC_NEF_Set_Collapse
!
!
!>    \brief <b> MQC_NEF_Set_Collapse applies the underlying element's
!>    collapse to every (I,J) entry of a MQC_NEF_Set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Collapse returns a new MQC_NEF_Set of the same (I,J)
!>    shape as NEFSet, where every entry is NEFSet(I,J)%collapse(...).
!>    If NEFSet holds Matrix elements, collapse reduces rank 2 to 1 (the
!>    result set holds Vector elements) and accepts spinType/swapOrder,
!>    matching MQC_NEF_Matrix's own collapse. If NEFSet holds R4
!>    elements, collapse reduces rank 4 to 2 (the result set holds Matrix
!>    elements) and accepts spinType/spinTypeKet/Order, matching
!>    MQC_NEF_R4's own collapse. NEFSet holding Vector elements (already
!>    rank 1) or being 'open' is an error, since MQC_NEF_Vector has no
!>    collapse of its own.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The set whose elements are each collapsed.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        Character(len=*),Optional
!>        Forwarded to the underlying Matrix/R4 collapse.
!>    \endverbatim
!>
!>    \param[in] swapOrder
!>    \verbatim
!>        Logical,Optional
!>        Forwarded to the underlying Matrix collapse only.
!>    \endverbatim
!>
!>    \param[in] spinTypeKet
!>    \verbatim
!>        Character(len=*),Optional
!>        Forwarded to the underlying R4 collapse only.
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Integer(kind=int64),Dimension(4),Optional
!>        Forwarded to the underlying R4 collapse only.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_collapse(NEFSet,spinType,swapOrder,spinTypeKet,Order) result(NEFSetOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      character(len=*),intent(in),optional::spinType,spinTypeKet
      logical,intent(in),optional::swapOrder
      integer(kind=int64),dimension(4),intent(in),optional::Order
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      select case(trim(NEFSet%Type()))
      case('matrix')
        call NEFSetOut%init(NEFSet%NumI,NEFSet%NumJ)
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            call NEFSetOut%put(NEFSet%Set_Matrix(I,J)%collapse(spinType=spinType,swapOrder=swapOrder),I,J)
          end do
        end do
      case('r4')
        call NEFSetOut%init(NEFSet%NumI,NEFSet%NumJ)
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            call NEFSetOut%put(NEFSet%Set_R4(I,J)%collapse(&
                 spinType=spinType,spinTypeKet=spinTypeKet,Order=Order),I,J)
          end do
        end do
      case default
        call mqc_error_A('mqc_nef_set_collapse: NEF_Set kind does not support '//&
             'collapse (must be matrix or r4)',6,'Array_Type',NEFSet%Type())
      end select
!
      end function mqc_nef_set_collapse
!
!
!     PROCEDURE MQC_NEF_Set_Outer_Collapse
!
!
!>    \brief <b> MQC_NEF_Set_Outer_Collapse folds a MQC_NEF_Set's (I,J)
!>    grid indices into the stored element's own internal axis structure,
!>    producing a single combined element wrapped in a (1,1) set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Outer_Collapse combines the outer (I,J) grid indexing of
!>    a MQC_NEF_Set with the internal space-block structure of whichever
!>    kind it currently holds, producing one combined NEF element of the
!>    same kind, wrapped in a (1,1) MQC_NEF_Set.
!>
!>    If NEFSet holds Vector elements and jointIJ is true (the default
!>    for vector), I and J first jointly fold into one outer index,
!>    outerIdx = (I-1)*NumJ + J, which then combines with the Vector's
!>    own single axis: new block = (outerIdx-1)*NumSpace + origBlock. The
!>    result is a single MQC_NEF_Vector. If jointIJ is false, I instead
!>    merges with the Vector's axis to form a new ROW (new block =
!>    (I-1)*NumSpace + origBlock), and J becomes a plain new COLUMN axis
!>    (NumJ trivial one-element spaces); the result is a single
!>    MQC_NEF_Matrix.
!>
!>    If NEFSet holds Matrix elements and jointIJ is false (the default
!>    for matrix), I folds into the row axis and J folds into the column
!>    axis independently: new row block = (I-1)*NumSpace_row +
!>    origRowBlock, new col block = (J-1)*NumSpace_col + origColBlock. If
!>    jointIJ is true, I and J first jointly fold into one outer index,
!>    G = (I-1)*NumJ + J, which then folds into a single caller-chosen
!>    axis (jointAxis, 1=row/2=col, default 1); the other axis is left
!>    untouched (copied from a representative cell, valid because
!>    MQC_NEF_Set%put already guarantees identical internal dimensions on
!>    every axis for every cell). Either way, the result is a single
!>    MQC_NEF_Matrix.
!>
!>    If NEFSet holds R4 elements and jointIJ is false (the default for
!>    r4), I folds into axis IAxisR4 and J folds into axis JAxisR4
!>    (default axis 1 and axis 4 respectively), independently; the other
!>    two axes are untouched. If jointIJ is true, I and J first jointly
!>    fold into G = (I-1)*NumJ + J, which then folds into a single
!>    caller-chosen axis (jointAxis, 1-4, default 1); the other three
!>    axes are untouched. Either way, the result is a single MQC_NEF_R4.
!>
!>    Converted indexing per case, (grid)*(original block) -> (new
!>    block), u/v/w/x denoting an original per-axis block index:
!>
!>      vector, jointIJ=T:  (I,J)*(u)     -> (IJu)
!>                           IJu = ((I-1)*NumJ+(J-1))*NumSpace + u
!>      vector, jointIJ=F:  (I,J)*(u)     -> (Iu,J)
!>                           Iu  = (I-1)*NumSpace + u
!>      matrix, jointIJ=F:  (I,J)*(u,v)   -> (Iu,Jv)
!>                           Iu  = (I-1)*NumSpace_row + u
!>                           Jv  = (J-1)*NumSpace_col + v
!>      matrix, jointIJ=T,
!>        jointAxis=1:      (I,J)*(u,v)   -> (IJu,v)
!>                           IJu = ((I-1)*NumJ+(J-1))*NumSpace_row + u
!>        jointAxis=2:      (I,J)*(u,v)   -> (u,IJv)
!>                           IJv = ((I-1)*NumJ+(J-1))*NumSpace_col + v
!>      r4, jointIJ=F:      (I,J)*(u,v,w,x) -> fold I into axis IAxisR4,
!>                           J into axis JAxisR4, e.g. IAxisR4=1,JAxisR4=4:
!>                           (I,J)*(u,v,w,x) -> (Iu,v,w,Jx)
!>                           Iu = (I-1)*NumSpace_1 + u
!>                           Jx = (J-1)*NumSpace_4 + x
!>                           (v,w unchanged; generalizes to any axis pair)
!>      r4, jointIJ=T,
!>        jointAxis=k:      (I,J)*(u,v,w,x) -> fold G=(I-1)*NumJ+J into
!>                           axis k only, e.g. jointAxis=1:
!>                           (I,J)*(u,v,w,x) -> (Gu,v,w,x)
!>                           Gu = (G-1)*NumSpace_1 + u
!>                           (v,w,x unchanged; generalizes to any axis)
!>
!>    Because folding a grid index into an axis only changes how many
!>    blocks exist along that axis (never a block's own size), every
!>    surviving block is copied verbatim from its source cell -- no
!>    element-wise re-indexing within a block is needed, unlike collapse.
!>    The one exception is Vector with jointIJ=false, where each source
!>    Vector block (one axis) is transcribed element-by-element into the
!>    corresponding single-column Matrix block (two axes), since the
!>    underlying payload type itself changes from MQC_SCF_Eigenvalues to
!>    MQC_SCF_Integral.
!>
!>    NEFSet must not be 'open'.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The (I,J) grid whose indices are folded into the element's own
!>        axis structure. Must not be 'open'.
!>    \endverbatim
!>
!>    \param[in] jointIJ
!>    \verbatim
!>        jointIJ is Logical,Optional
!>        Whether I and J fold jointly into one combined index (true) or
!>        independently into two separate axes (false). Default: true
!>        when NEFSet holds Vector elements, false when it holds Matrix
!>        or R4 elements.
!>    \endverbatim
!>
!>    \param[in] jointAxis
!>    \verbatim
!>        jointAxis is Integer,Optional
!>        Only used when jointIJ is true and NEFSet holds Matrix (1=row,
!>        2=col) or R4 (1-4) elements. Which single axis the combined
!>        (I,J) index folds into. Default 1.
!>    \endverbatim
!>
!>    \param[in] IAxisR4
!>    \verbatim
!>        IAxisR4 is Integer,Optional
!>        Only used when jointIJ is false and NEFSet%Type()=='r4'. The
!>        physical R4 axis (1-4) that the grid's I index folds into.
!>        Default 1. Must differ from JAxisR4.
!>    \endverbatim
!>
!>    \param[in] JAxisR4
!>    \verbatim
!>        JAxisR4 is Integer,Optional
!>        Only used when jointIJ is false and NEFSet%Type()=='r4'. The
!>        physical R4 axis (1-4) that the grid's J index folds into.
!>        Default 4. Must differ from IAxisR4.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_outer_collapse(NEFSet,jointIJ,jointAxis,IAxisR4,JAxisR4) result(NEFSetOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      logical,intent(in),optional::jointIJ
      integer,intent(in),optional::jointAxis
      integer,intent(in),optional::IAxisR4,JAxisR4
      type(mqc_nef_set)::NEFSetOut
      logical::myJointIJ
      integer::myJointAxis,myIAxisR4,myJAxisR4,axisNum
      integer(kind=int64)::I,J,myNumI,myNumJ,G
      integer(kind=int64)::b,origNumSpace,outerIdx,newIdx,ka
      integer(kind=int64)::r,c,origNumSpace_row,origNumSpace_col
      integer(kind=int64)::newRowId,newColId,origId,newId
      integer(kind=int64)::b1,b2,b3,b4,nb1,nb2,nb3,nb4
      integer(kind=int64)::origNumSpace1,origNumSpace2,origNumSpace3,origNumSpace4
      integer(kind=int64)::axisNumSpace(4),foldGridSize(4),kk
      type(mqc_nef_vector)::cellVector,NEFVecOut
      type(mqc_nef_matrix)::cellMatrix,NEFMatOut
      type(mqc_nef_r4)::cellR4,NEFR4Out
      type(mqc_scf_eigenvalues),dimension(:),allocatable::New_Map_1,New_Map_2,New_Map_3,New_Map_4
      type(mqc_scalar)::tmpScalar
!
      if(present(jointIJ)) then
        myJointIJ = jointIJ
      else
        myJointIJ = (trim(NEFSet%Type()).eq.'vector')
      end if
      if(present(jointAxis)) then
        myJointAxis = jointAxis
      else
        myJointAxis = 1
      end if
      if(present(IAxisR4)) then
        myIAxisR4 = IAxisR4
      else
        myIAxisR4 = 1
      end if
      if(present(JAxisR4)) then
        myJAxisR4 = JAxisR4
      else
        myJAxisR4 = 4
      end if
!
      myNumI = NEFSet%NumI
      myNumJ = NEFSet%NumJ
!
      select case(trim(NEFSet%Type()))
      case('vector')
        origNumSpace = NEFSet%NumSpace(1)
        if(myJointIJ) then
          allocate(New_Map_1(myNumI*myNumJ*origNumSpace))
          do outerIdx = 1, myNumI*myNumJ
            do b = 1, origNumSpace
              call New_Map_1((outerIdx-1)*origNumSpace+b)%init(&
                   NEFSet%Size_vector_alpha_1(b),0,lengthBeta=NEFSet%Size_vector_beta_1(b))
            end do
          end do
          call NEFVecOut%init(New_Map_1)
          do I = 1, myNumI
            do J = 1, myNumJ
              cellVector = NEFSet%atVector(int(I),int(J))
              outerIdx = (I-1)*myNumJ+J
              do b = 1, origNumSpace
                if(cellVector%Zero_Flag(b)) cycle
                newIdx = (outerIdx-1)*origNumSpace+b
                NEFVecOut%Main_Vector(newIdx) = cellVector%Main_Vector(b)
                NEFVecOut%Zero_Flag(newIdx) = .false.
              end do
            end do
          end do
          call NEFSetOut%init(1,1)
          call NEFSetOut%put(NEFVecOut,1,1)
        else
          allocate(New_Map_1(myNumI*origNumSpace))
          do I = 1, myNumI
            do b = 1, origNumSpace
              call New_Map_1((I-1)*origNumSpace+b)%init(&
                   NEFSet%Size_vector_alpha_1(b),0,lengthBeta=NEFSet%Size_vector_beta_1(b))
            end do
          end do
          allocate(New_Map_2(myNumJ))
          do J = 1, myNumJ
            call New_Map_2(J)%init(1,0,lengthBeta=1)
          end do
          call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2)
          do I = 1, myNumI
            do J = 1, myNumJ
              cellVector = NEFSet%atVector(int(I),int(J))
              do b = 1, origNumSpace
                if(cellVector%Zero_Flag(b)) cycle
                newRowId = (I-1)*origNumSpace+b
                newId = NEFMatOut%blockIndex([1,newRowId],[2,J])
                call NEFMatOut%Main_Matrix(newId)%init(&
                     NEFSet%Size_vector_alpha_1(b),NEFSet%Size_vector_beta_1(b),&
                     nAlpha2=1,nBeta2=1)
                NEFMatOut%Zero_Flag(newId) = .false.
                do ka = 1, NEFSet%Size_vector_alpha_1(b)
                  tmpScalar = cellVector%Main_Vector(b)%at(ka,'alpha')
                  call NEFMatOut%Main_Matrix(newId)%put(tmpScalar,ka,1_int64,'aa')
                end do
                do ka = 1, NEFSet%Size_vector_beta_1(b)
                  tmpScalar = cellVector%Main_Vector(b)%at(ka,'beta')
                  call NEFMatOut%Main_Matrix(newId)%put(tmpScalar,ka,1_int64,'bb')
                end do
              end do
            end do
          end do
          call NEFSetOut%init(1,1)
          call NEFSetOut%put(NEFMatOut,1,1)
        end if
!
      case('matrix')
        origNumSpace_row = NEFSet%NumSpace(1)
        origNumSpace_col = NEFSet%NumSpace(2)
        if(.not.myJointIJ) then
          allocate(New_Map_1(myNumI*origNumSpace_row))
          do I = 1, myNumI
            do r = 1, origNumSpace_row
              call New_Map_1((I-1)*origNumSpace_row+r)%init(&
                   NEFSet%Size_vector_alpha_1(r),0,lengthBeta=NEFSet%Size_vector_beta_1(r))
            end do
          end do
          allocate(New_Map_2(myNumJ*origNumSpace_col))
          do J = 1, myNumJ
            do c = 1, origNumSpace_col
              call New_Map_2((J-1)*origNumSpace_col+c)%init(&
                   NEFSet%Size_vector_alpha_2(c),0,lengthBeta=NEFSet%Size_vector_beta_2(c))
            end do
          end do
        else
          if(myJointAxis.ne.1.and.myJointAxis.ne.2) call mqc_error_i(&
               'mqc_nef_set_outer_collapse: jointAxis out-of-bounds for matrix',6,&
               'jointAxis',myJointAxis)
          if(myJointAxis.eq.1) then
            allocate(New_Map_1(myNumI*myNumJ*origNumSpace_row))
            do G = 1, myNumI*myNumJ
              do r = 1, origNumSpace_row
                call New_Map_1((G-1)*origNumSpace_row+r)%init(&
                     NEFSet%Size_vector_alpha_1(r),0,lengthBeta=NEFSet%Size_vector_beta_1(r))
              end do
            end do
            allocate(New_Map_2(origNumSpace_col))
            do c = 1, origNumSpace_col
              call New_Map_2(c)%init(&
                   NEFSet%Size_vector_alpha_2(c),0,lengthBeta=NEFSet%Size_vector_beta_2(c))
            end do
          else
            allocate(New_Map_1(origNumSpace_row))
            do r = 1, origNumSpace_row
              call New_Map_1(r)%init(&
                   NEFSet%Size_vector_alpha_1(r),0,lengthBeta=NEFSet%Size_vector_beta_1(r))
            end do
            allocate(New_Map_2(myNumI*myNumJ*origNumSpace_col))
            do G = 1, myNumI*myNumJ
              do c = 1, origNumSpace_col
                call New_Map_2((G-1)*origNumSpace_col+c)%init(&
                     NEFSet%Size_vector_alpha_2(c),0,lengthBeta=NEFSet%Size_vector_beta_2(c))
              end do
            end do
          end if
        end if
        call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2)
        do I = 1, myNumI
          do J = 1, myNumJ
            cellMatrix = NEFSet%atMatrix(int(I),int(J))
            G = (I-1)*myNumJ+J
            do r = 1, origNumSpace_row
              do c = 1, origNumSpace_col
                origId = cellMatrix%blockIndex([1,r],[2,c])
                if(cellMatrix%Zero_Flag(origId)) cycle
                if(.not.myJointIJ) then
                  newRowId = (I-1)*origNumSpace_row+r
                  newColId = (J-1)*origNumSpace_col+c
                else if(myJointAxis.eq.1) then
                  newRowId = (G-1)*origNumSpace_row+r
                  newColId = c
                else
                  newRowId = r
                  newColId = (G-1)*origNumSpace_col+c
                end if
                newId = NEFMatOut%blockIndex([1,newRowId],[2,newColId])
                NEFMatOut%Main_Matrix(newId) = cellMatrix%Main_Matrix(origId)
                NEFMatOut%Zero_Flag(newId) = .false.
              end do
            end do
          end do
        end do
        call NEFSetOut%init(1,1)
        call NEFSetOut%put(NEFMatOut,1,1)
!
      case('r4')
        origNumSpace1 = NEFSet%NumSpace(1)
        origNumSpace2 = NEFSet%NumSpace(2)
        origNumSpace3 = NEFSet%NumSpace(3)
        origNumSpace4 = NEFSet%NumSpace(4)
        axisNumSpace = [origNumSpace1,origNumSpace2,origNumSpace3,origNumSpace4]
        foldGridSize = 1
        if(.not.myJointIJ) then
          if(myIAxisR4.lt.1.or.myIAxisR4.gt.4) call mqc_error_i(&
               'mqc_nef_set_outer_collapse: IAxisR4 out-of-bounds',6,'IAxisR4',myIAxisR4)
          if(myJAxisR4.lt.1.or.myJAxisR4.gt.4) call mqc_error_i(&
               'mqc_nef_set_outer_collapse: JAxisR4 out-of-bounds',6,'JAxisR4',myJAxisR4)
          if(myIAxisR4.eq.myJAxisR4) call mqc_error_i(&
               'mqc_nef_set_outer_collapse: IAxisR4 and JAxisR4 must differ',6,&
               'IAxisR4',myIAxisR4)
          foldGridSize(myIAxisR4) = myNumI
          foldGridSize(myJAxisR4) = myNumJ
        else
          if(myJointAxis.lt.1.or.myJointAxis.gt.4) call mqc_error_i(&
               'mqc_nef_set_outer_collapse: jointAxis out-of-bounds for r4',6,&
               'jointAxis',myJointAxis)
          foldGridSize(myJointAxis) = myNumI*myNumJ
        end if
!
        NEFR4Out%NumSpace_1 = axisNumSpace(1)*foldGridSize(1)
        NEFR4Out%NumSpace_2 = axisNumSpace(2)*foldGridSize(2)
        NEFR4Out%NumSpace_3 = axisNumSpace(3)*foldGridSize(3)
        NEFR4Out%NumSpace_4 = axisNumSpace(4)*foldGridSize(4)
        allocate(NEFR4Out%Main_R4(NEFR4Out%NumSpace_1*NEFR4Out%NumSpace_2*&
                                   NEFR4Out%NumSpace_3*NEFR4Out%NumSpace_4))
        allocate(NEFR4Out%Zero_Flag(size(NEFR4Out%Main_R4)))
        NEFR4Out%Zero_Flag = .true.
!
        allocate(New_Map_1(axisNumSpace(1)*foldGridSize(1)))
        do kk = 1, foldGridSize(1)
          do b1 = 1, axisNumSpace(1)
            call New_Map_1((kk-1)*axisNumSpace(1)+b1)%init(&
                 NEFSet%Size_vector_alpha_1(b1),0,lengthBeta=NEFSet%Size_vector_beta_1(b1))
          end do
        end do
        NEFR4Out%Map_Matrix_1 = New_Map_1
        NEFR4Out%Size_vector_alpha_1 = &
             [(NEFSet%Size_vector_alpha_1(int(mod(kk-1,axisNumSpace(1)))+1),kk=1,NEFR4Out%NumSpace_1)]
        NEFR4Out%Size_vector_beta_1 = &
             [(NEFSet%Size_vector_beta_1(int(mod(kk-1,axisNumSpace(1)))+1),kk=1,NEFR4Out%NumSpace_1)]
!
        allocate(New_Map_2(axisNumSpace(2)*foldGridSize(2)))
        do kk = 1, foldGridSize(2)
          do b2 = 1, axisNumSpace(2)
            call New_Map_2((kk-1)*axisNumSpace(2)+b2)%init(&
                 NEFSet%Size_vector_alpha_2(b2),0,lengthBeta=NEFSet%Size_vector_beta_2(b2))
          end do
        end do
        NEFR4Out%Map_Matrix_2 = New_Map_2
        NEFR4Out%Size_vector_alpha_2 = &
             [(NEFSet%Size_vector_alpha_2(int(mod(kk-1,axisNumSpace(2)))+1),kk=1,NEFR4Out%NumSpace_2)]
        NEFR4Out%Size_vector_beta_2 = &
             [(NEFSet%Size_vector_beta_2(int(mod(kk-1,axisNumSpace(2)))+1),kk=1,NEFR4Out%NumSpace_2)]
!
        allocate(New_Map_3(axisNumSpace(3)*foldGridSize(3)))
        do kk = 1, foldGridSize(3)
          do b3 = 1, axisNumSpace(3)
            call New_Map_3((kk-1)*axisNumSpace(3)+b3)%init(&
                 NEFSet%Size_vector_alpha_3(b3),0,lengthBeta=NEFSet%Size_vector_beta_3(b3))
          end do
        end do
        NEFR4Out%Map_Matrix_3 = New_Map_3
        NEFR4Out%Size_vector_alpha_3 = &
             [(NEFSet%Size_vector_alpha_3(int(mod(kk-1,axisNumSpace(3)))+1),kk=1,NEFR4Out%NumSpace_3)]
        NEFR4Out%Size_vector_beta_3 = &
             [(NEFSet%Size_vector_beta_3(int(mod(kk-1,axisNumSpace(3)))+1),kk=1,NEFR4Out%NumSpace_3)]
!
        allocate(New_Map_4(axisNumSpace(4)*foldGridSize(4)))
        do kk = 1, foldGridSize(4)
          do b4 = 1, axisNumSpace(4)
            call New_Map_4((kk-1)*axisNumSpace(4)+b4)%init(&
                 NEFSet%Size_vector_alpha_4(b4),0,lengthBeta=NEFSet%Size_vector_beta_4(b4))
          end do
        end do
        NEFR4Out%Map_Matrix_4 = New_Map_4
        NEFR4Out%Size_vector_alpha_4 = &
             [(NEFSet%Size_vector_alpha_4(int(mod(kk-1,axisNumSpace(4)))+1),kk=1,NEFR4Out%NumSpace_4)]
        NEFR4Out%Size_vector_beta_4 = &
             [(NEFSet%Size_vector_beta_4(int(mod(kk-1,axisNumSpace(4)))+1),kk=1,NEFR4Out%NumSpace_4)]
!
        do I = 1, myNumI
          do J = 1, myNumJ
            cellR4 = NEFSet%atR4(int(I),int(J))
            G = (I-1)*myNumJ+J
            do b1 = 1, origNumSpace1
              do b2 = 1, origNumSpace2
                do b3 = 1, origNumSpace3
                  do b4 = 1, origNumSpace4
                    origId = cellR4%blockIndex([1,b1],[2,b2],[3,b3],[4,b4])
                    if(cellR4%Zero_Flag(origId)) cycle
                    nb1=b1; nb2=b2; nb3=b3; nb4=b4
                    if(.not.myJointIJ) then
                      if(myIAxisR4.eq.1) then
                        nb1=(I-1)*origNumSpace1+b1
                      else if(myIAxisR4.eq.2) then
                        nb2=(I-1)*origNumSpace2+b2
                      else if(myIAxisR4.eq.3) then
                        nb3=(I-1)*origNumSpace3+b3
                      else
                        nb4=(I-1)*origNumSpace4+b4
                      end if
                      if(myJAxisR4.eq.1) then
                        nb1=(J-1)*origNumSpace1+b1
                      else if(myJAxisR4.eq.2) then
                        nb2=(J-1)*origNumSpace2+b2
                      else if(myJAxisR4.eq.3) then
                        nb3=(J-1)*origNumSpace3+b3
                      else
                        nb4=(J-1)*origNumSpace4+b4
                      end if
                    else
                      if(myJointAxis.eq.1) then
                        nb1=(G-1)*origNumSpace1+b1
                      else if(myJointAxis.eq.2) then
                        nb2=(G-1)*origNumSpace2+b2
                      else if(myJointAxis.eq.3) then
                        nb3=(G-1)*origNumSpace3+b3
                      else
                        nb4=(G-1)*origNumSpace4+b4
                      end if
                    end if
                    newId = NEFR4Out%blockIndex([1,nb1],[2,nb2],[3,nb3],[4,nb4])
                    NEFR4Out%Main_R4(newId) = cellR4%Main_R4(origId)
                    NEFR4Out%Zero_Flag(newId) = .false.
                  end do
                end do
              end do
            end do
          end do
        end do
        call NEFSetOut%init(1,1)
        call NEFSetOut%put(NEFR4Out,1,1)
!
      case default
        call mqc_error_A('mqc_nef_set_outer_collapse: NEF_Set is open (holds no elements)',&
             6,'Array_Type',NEFSet%Type())
      end select
!
      end function mqc_nef_set_outer_collapse
!
!
!     PROCEDURE MQC_NEF_Set_Outer_Expand
!
!
!>    \brief <b> MQC_NEF_Set_Outer_Expand is the counterpart of
!>    MQC_NEF_Set%outerCollapse: it unfolds a (1,1)-wrapped combined NEF
!>    element back into a full (NumI,NumJ) MQC_NEF_Set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Outer_Expand inverts outerCollapse exactly, covering
!>    every jointIJ/jointAxis/IAxisR4/JAxisR4 mode outerCollapse supports.
!>    Because outerCollapse discards NumI, NumJ, and the original
!>    per-(I,J)-cell axis block sizes when it folds everything into one
!>    combined element, outerExpand needs that information supplied back
!>    explicitly -- either via templateSet (the exact NEF_Set that was
!>    given to outerCollapse, whose NumI/NumJ/per-axis block sizes are
!>    read directly) or via NumI/NumJ plus the explicit size_vector_
!>    alpha/beta_1..4 arrays (mirroring MQC_NEF_Set's own internal
!>    bookkeeping). origKind tells outerExpand what kind the ORIGINAL
!>    (pre-outerCollapse) grid cells held -- this cannot be inferred from
!>    NEFSet%Type() alone, since e.g. a Matrix-kind combined object can
!>    result from either "vector, jointIJ=false" or "matrix, any
!>    jointIJ" in outerCollapse, each requiring a different unfold.
!>
!>    jointIJ/jointAxis/IAxisR4/JAxisR4 must be passed the SAME values
!>    used in the original outerCollapse call (defaults match
!>    outerCollapse's own: jointIJ true for origKind=vector, false
!>    otherwise; jointAxis 1; IAxisR4 1; JAxisR4 4).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The (1,1)-wrapped combined element to unfold (outerCollapse's
!>        result).
!>    \endverbatim
!>
!>    \param[in] origKind
!>    \verbatim
!>        origKind is Character(len=*)
!>        'vector', 'matrix', or 'r4' -- the kind the original grid
!>        cells held before outerCollapse.
!>    \endverbatim
!>
!>    \param[in] NumI,NumJ
!>    \verbatim
!>        Integer,Optional
!>        The original grid shape. Required (together) when templateSet
!>        is absent; ignored (templateSet's own NumI/NumJ used instead)
!>        when templateSet is present.
!>    \endverbatim
!>
!>    \param[in] size_vector_alpha_1,size_vector_beta_1 ... _4
!>    \verbatim
!>        Integer,Dimension(:),Optional
!>        The original per-cell per-axis block sizes. Required (the
!>        subset matching origKind: axis 1 only for vector, 1&2 for
!>        matrix, all four for r4) when templateSet is absent; ignored
!>        when templateSet is present.
!>    \endverbatim
!>
!>    \param[in] templateSet
!>    \verbatim
!>        Type(MQC_NEF_Set),Optional
!>        The exact NEF_Set originally passed to outerCollapse. When
!>        present, NumI/NumJ and all size_vector_alpha/beta arrays are
!>        read from it directly instead of from the explicit arguments
!>        above.
!>    \endverbatim
!>
!>    \param[in] jointIJ
!>    \verbatim
!>        Logical,Optional
!>        Must match the value used in the original outerCollapse call.
!>    \endverbatim
!>
!>    \param[in] jointAxis
!>    \verbatim
!>        Integer,Optional
!>        Must match the value used in the original outerCollapse call.
!>    \endverbatim
!>
!>    \param[in] IAxisR4,JAxisR4
!>    \verbatim
!>        Integer,Optional
!>        Must match the values used in the original outerCollapse call.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_outer_expand(NEFSet,origKind,NumI,NumJ,&
           size_vector_alpha_1,size_vector_beta_1,size_vector_alpha_2,size_vector_beta_2,&
           size_vector_alpha_3,size_vector_beta_3,size_vector_alpha_4,size_vector_beta_4,&
           templateSet,jointIJ,jointAxis,IAxisR4,JAxisR4) result(NEFSetOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      character(len=*),intent(in)::origKind
      integer,intent(in),optional::NumI,NumJ
      integer,dimension(:),intent(in),optional::size_vector_alpha_1,size_vector_beta_1
      integer,dimension(:),intent(in),optional::size_vector_alpha_2,size_vector_beta_2
      integer,dimension(:),intent(in),optional::size_vector_alpha_3,size_vector_beta_3
      integer,dimension(:),intent(in),optional::size_vector_alpha_4,size_vector_beta_4
      type(mqc_nef_set),intent(in),optional::templateSet
      logical,intent(in),optional::jointIJ
      integer,intent(in),optional::jointAxis,IAxisR4,JAxisR4
      type(mqc_nef_set)::NEFSetOut
!
      logical::myJointIJ
      integer::myJointAxis,myIAxisR4,myJAxisR4
      integer(kind=int64)::myNumI,myNumJ,I,J,G,outerIdx,b,origNumSpace
      integer(kind=int64)::r,c,origNumSpace_row,origNumSpace_col
      integer(kind=int64)::b1,b2,b3,b4,nb1,nb2,nb3,nb4
      integer(kind=int64)::origNumSpace1,origNumSpace2,origNumSpace3,origNumSpace4
      integer(kind=int64)::origId,newId,newRowId,newColId,newIdx,ka
      integer,dimension(:),allocatable::mySizeAlpha1,mySizeBeta1,mySizeAlpha2,mySizeBeta2
      integer,dimension(:),allocatable::mySizeAlpha3,mySizeBeta3,mySizeAlpha4,mySizeBeta4
      type(mqc_scf_eigenvalues),dimension(:),allocatable::cellMap1,cellMap2,cellMap3,cellMap4
      type(mqc_nef_vector)::combVector,cellVector
      type(mqc_nef_matrix)::combMatrix,cellMatrix
      type(mqc_nef_r4)::combR4,cellR4
      type(mqc_scalar)::tmpScalar
!
!     --- resolve grid shape ---
      if(present(templateSet)) then
        myNumI = templateSet%NumI
        myNumJ = templateSet%NumJ
      else
        if(.not.(present(NumI).and.present(NumJ))) call mqc_error(&
             'mqc_nef_set_outer_expand: NumI and NumJ are required when '//&
             'templateSet is absent',6)
        myNumI = NumI
        myNumJ = NumJ
      end if
!
!     --- resolve per-original-axis block sizes ---
      select case(trim(origKind))
      case('vector')
        if(present(templateSet)) then
          mySizeAlpha1 = templateSet%Size_vector_alpha_1
          mySizeBeta1 = templateSet%Size_vector_beta_1
        else
          if(.not.(present(size_vector_alpha_1).and.present(size_vector_beta_1))) call mqc_error(&
               'mqc_nef_set_outer_expand: size_vector_alpha_1 and '//&
               'size_vector_beta_1 are required for origKind=vector when '//&
               'templateSet is absent',6)
          mySizeAlpha1 = size_vector_alpha_1
          mySizeBeta1 = size_vector_beta_1
        end if
      case('matrix')
        if(present(templateSet)) then
          mySizeAlpha1 = templateSet%Size_vector_alpha_1
          mySizeBeta1 = templateSet%Size_vector_beta_1
          mySizeAlpha2 = templateSet%Size_vector_alpha_2
          mySizeBeta2 = templateSet%Size_vector_beta_2
        else
          if(.not.(present(size_vector_alpha_1).and.present(size_vector_beta_1).and.&
                   present(size_vector_alpha_2).and.present(size_vector_beta_2))) call mqc_error(&
               'mqc_nef_set_outer_expand: size_vector_alpha_1/2 and '//&
               'size_vector_beta_1/2 are required for origKind=matrix when '//&
               'templateSet is absent',6)
          mySizeAlpha1 = size_vector_alpha_1
          mySizeBeta1 = size_vector_beta_1
          mySizeAlpha2 = size_vector_alpha_2
          mySizeBeta2 = size_vector_beta_2
        end if
      case('r4')
        if(present(templateSet)) then
          mySizeAlpha1 = templateSet%Size_vector_alpha_1
          mySizeBeta1 = templateSet%Size_vector_beta_1
          mySizeAlpha2 = templateSet%Size_vector_alpha_2
          mySizeBeta2 = templateSet%Size_vector_beta_2
          mySizeAlpha3 = templateSet%Size_vector_alpha_3
          mySizeBeta3 = templateSet%Size_vector_beta_3
          mySizeAlpha4 = templateSet%Size_vector_alpha_4
          mySizeBeta4 = templateSet%Size_vector_beta_4
        else
          if(.not.(present(size_vector_alpha_1).and.present(size_vector_beta_1).and.&
                   present(size_vector_alpha_2).and.present(size_vector_beta_2).and.&
                   present(size_vector_alpha_3).and.present(size_vector_beta_3).and.&
                   present(size_vector_alpha_4).and.present(size_vector_beta_4))) call mqc_error(&
               'mqc_nef_set_outer_expand: all eight size_vector_alpha/beta_1..4 '//&
               'arguments are required for origKind=r4 when templateSet is absent',6)
          mySizeAlpha1 = size_vector_alpha_1 ; mySizeBeta1 = size_vector_beta_1
          mySizeAlpha2 = size_vector_alpha_2 ; mySizeBeta2 = size_vector_beta_2
          mySizeAlpha3 = size_vector_alpha_3 ; mySizeBeta3 = size_vector_beta_3
          mySizeAlpha4 = size_vector_alpha_4 ; mySizeBeta4 = size_vector_beta_4
        end if
      case default
        call mqc_error_A('mqc_nef_set_outer_expand: origKind must be vector, '//&
             'matrix, or r4',6,'origKind',origKind)
      end select
!
!     --- resolve mode selectors (same defaults as outerCollapse) ---
      if(present(jointIJ)) then
        myJointIJ = jointIJ
      else
        myJointIJ = (trim(origKind).eq.'vector')
      end if
      if(present(jointAxis)) then
        myJointAxis = jointAxis
      else
        myJointAxis = 1
      end if
      if(present(IAxisR4)) then
        myIAxisR4 = IAxisR4
      else
        myIAxisR4 = 1
      end if
      if(present(JAxisR4)) then
        myJAxisR4 = JAxisR4
      else
        myJAxisR4 = 4
      end if
!
      select case(trim(origKind))
      case('vector')
        origNumSpace = size(mySizeAlpha1)
        allocate(cellMap1(origNumSpace))
        do b = 1, origNumSpace
          call cellMap1(b)%init(mySizeAlpha1(b),0,lengthBeta=mySizeBeta1(b))
        end do
        call NEFSetOut%init(int(myNumI),int(myNumJ))
        if(myJointIJ) then
          combVector = NEFSet%Set_Vector(1,1)
          do I = 1, myNumI
            do J = 1, myNumJ
              outerIdx = (I-1)*myNumJ+J
              call cellVector%init(cellMap1)
              do b = 1, origNumSpace
                newIdx = (outerIdx-1)*origNumSpace+b
                if(combVector%Zero_Flag(newIdx)) cycle
                cellVector%Main_Vector(b) = combVector%Main_Vector(newIdx)
                cellVector%Zero_Flag(b) = .false.
              end do
              call NEFSetOut%put(cellVector,int(I),int(J))
            end do
          end do
        else
          combMatrix = NEFSet%Set_Matrix(1,1)
          do I = 1, myNumI
            do J = 1, myNumJ
              call cellVector%init(cellMap1)
              do b = 1, origNumSpace
                newRowId = (I-1)*origNumSpace+b
                newId = combMatrix%blockIndex([1,newRowId],[2,J])
                if(combMatrix%Zero_Flag(newId)) cycle
                call cellVector%Main_Vector(b)%init(&
                     mySizeAlpha1(b),0,lengthBeta=mySizeBeta1(b))
                do ka = 1, mySizeAlpha1(b)
                  tmpScalar = combMatrix%Main_Matrix(newId)%at(ka,1_int64,'aa')
                  call cellVector%Main_Vector(b)%put(tmpScalar,ka,'alpha')
                end do
                do ka = 1, mySizeBeta1(b)
                  tmpScalar = combMatrix%Main_Matrix(newId)%at(ka,1_int64,'bb')
                  call cellVector%Main_Vector(b)%put(tmpScalar,ka,'beta')
                end do
                cellVector%Zero_Flag(b) = .false.
              end do
              call NEFSetOut%put(cellVector,int(I),int(J))
            end do
          end do
        end if
!
      case('matrix')
        origNumSpace_row = size(mySizeAlpha1)
        origNumSpace_col = size(mySizeAlpha2)
        allocate(cellMap1(origNumSpace_row))
        do r = 1, origNumSpace_row
          call cellMap1(r)%init(mySizeAlpha1(r),0,lengthBeta=mySizeBeta1(r))
        end do
        allocate(cellMap2(origNumSpace_col))
        do c = 1, origNumSpace_col
          call cellMap2(c)%init(mySizeAlpha2(c),0,lengthBeta=mySizeBeta2(c))
        end do
        if(myJointIJ.and.myJointAxis.ne.1.and.myJointAxis.ne.2) call mqc_error_i(&
             'mqc_nef_set_outer_expand: jointAxis out-of-bounds for matrix',6,&
             'jointAxis',myJointAxis)
        combMatrix = NEFSet%Set_Matrix(1,1)
        call NEFSetOut%init(int(myNumI),int(myNumJ))
        do I = 1, myNumI
          do J = 1, myNumJ
            G = (I-1)*myNumJ+J
            call cellMatrix%init(cellMap1,initial_map_matrix_col=cellMap2)
            do r = 1, origNumSpace_row
              do c = 1, origNumSpace_col
                if(.not.myJointIJ) then
                  newRowId = (I-1)*origNumSpace_row+r
                  newColId = (J-1)*origNumSpace_col+c
                else if(myJointAxis.eq.1) then
                  newRowId = (G-1)*origNumSpace_row+r
                  newColId = c
                else
                  newRowId = r
                  newColId = (G-1)*origNumSpace_col+c
                end if
                newId = combMatrix%blockIndex([1,newRowId],[2,newColId])
                if(combMatrix%Zero_Flag(newId)) cycle
                origId = cellMatrix%blockIndex([1,r],[2,c])
                cellMatrix%Main_Matrix(origId) = combMatrix%Main_Matrix(newId)
                cellMatrix%Zero_Flag(origId) = .false.
              end do
            end do
            call NEFSetOut%put(cellMatrix,int(I),int(J))
          end do
        end do
!
      case('r4')
        origNumSpace1 = size(mySizeAlpha1)
        origNumSpace2 = size(mySizeAlpha2)
        origNumSpace3 = size(mySizeAlpha3)
        origNumSpace4 = size(mySizeAlpha4)
        allocate(cellMap1(origNumSpace1))
        do b1 = 1, origNumSpace1
          call cellMap1(b1)%init(mySizeAlpha1(b1),0,lengthBeta=mySizeBeta1(b1))
        end do
        allocate(cellMap2(origNumSpace2))
        do b2 = 1, origNumSpace2
          call cellMap2(b2)%init(mySizeAlpha2(b2),0,lengthBeta=mySizeBeta2(b2))
        end do
        allocate(cellMap3(origNumSpace3))
        do b3 = 1, origNumSpace3
          call cellMap3(b3)%init(mySizeAlpha3(b3),0,lengthBeta=mySizeBeta3(b3))
        end do
        allocate(cellMap4(origNumSpace4))
        do b4 = 1, origNumSpace4
          call cellMap4(b4)%init(mySizeAlpha4(b4),0,lengthBeta=mySizeBeta4(b4))
        end do
        if(.not.myJointIJ) then
          if(myIAxisR4.lt.1.or.myIAxisR4.gt.4.or.myJAxisR4.lt.1.or.myJAxisR4.gt.4.or.&
             myIAxisR4.eq.myJAxisR4) call mqc_error_i(&
               'mqc_nef_set_outer_expand: IAxisR4 and JAxisR4 must be '//&
               'distinct values in [1,4]',6,'IAxisR4',myIAxisR4,'JAxisR4',myJAxisR4)
        else
          if(myJointAxis.lt.1.or.myJointAxis.gt.4) call mqc_error_i(&
               'mqc_nef_set_outer_expand: jointAxis out-of-bounds for r4',6,&
               'jointAxis',myJointAxis)
        end if
        combR4 = NEFSet%Set_R4(1,1)
        call NEFSetOut%init(int(myNumI),int(myNumJ))
        do I = 1, myNumI
          do J = 1, myNumJ
            G = (I-1)*myNumJ+J
            if(allocated(cellR4%Main_R4)) deallocate(cellR4%Main_R4)
            if(allocated(cellR4%Zero_Flag)) deallocate(cellR4%Zero_Flag)
            cellR4%NumSpace_1 = origNumSpace1
            cellR4%NumSpace_2 = origNumSpace2
            cellR4%NumSpace_3 = origNumSpace3
            cellR4%NumSpace_4 = origNumSpace4
            allocate(cellR4%Main_R4(origNumSpace1*origNumSpace2*origNumSpace3*origNumSpace4))
            allocate(cellR4%Zero_Flag(size(cellR4%Main_R4)))
            cellR4%Zero_Flag = .true.
            cellR4%Map_Matrix_1 = cellMap1
            cellR4%Map_Matrix_2 = cellMap2
            cellR4%Map_Matrix_3 = cellMap3
            cellR4%Map_Matrix_4 = cellMap4
            cellR4%Size_vector_alpha_1 = mySizeAlpha1 ; cellR4%Size_vector_beta_1 = mySizeBeta1
            cellR4%Size_vector_alpha_2 = mySizeAlpha2 ; cellR4%Size_vector_beta_2 = mySizeBeta2
            cellR4%Size_vector_alpha_3 = mySizeAlpha3 ; cellR4%Size_vector_beta_3 = mySizeBeta3
            cellR4%Size_vector_alpha_4 = mySizeAlpha4 ; cellR4%Size_vector_beta_4 = mySizeBeta4
            do b1 = 1, origNumSpace1
              do b2 = 1, origNumSpace2
                do b3 = 1, origNumSpace3
                  do b4 = 1, origNumSpace4
                    nb1=b1 ; nb2=b2 ; nb3=b3 ; nb4=b4
                    if(.not.myJointIJ) then
                      if(myIAxisR4.eq.1) then
                        nb1=(I-1)*origNumSpace1+b1
                      else if(myIAxisR4.eq.2) then
                        nb2=(I-1)*origNumSpace2+b2
                      else if(myIAxisR4.eq.3) then
                        nb3=(I-1)*origNumSpace3+b3
                      else
                        nb4=(I-1)*origNumSpace4+b4
                      end if
                      if(myJAxisR4.eq.1) then
                        nb1=(J-1)*origNumSpace1+b1
                      else if(myJAxisR4.eq.2) then
                        nb2=(J-1)*origNumSpace2+b2
                      else if(myJAxisR4.eq.3) then
                        nb3=(J-1)*origNumSpace3+b3
                      else
                        nb4=(J-1)*origNumSpace4+b4
                      end if
                    else
                      if(myJointAxis.eq.1) then
                        nb1=(G-1)*origNumSpace1+b1
                      else if(myJointAxis.eq.2) then
                        nb2=(G-1)*origNumSpace2+b2
                      else if(myJointAxis.eq.3) then
                        nb3=(G-1)*origNumSpace3+b3
                      else
                        nb4=(G-1)*origNumSpace4+b4
                      end if
                    end if
                    newId = combR4%blockIndex([1,nb1],[2,nb2],[3,nb3],[4,nb4])
                    if(combR4%Zero_Flag(newId)) cycle
                    origId = cellR4%blockIndex([1,b1],[2,b2],[3,b3],[4,b4])
                    cellR4%Main_R4(origId) = combR4%Main_R4(newId)
                    cellR4%Zero_Flag(origId) = .false.
                  end do
                end do
              end do
            end do
            call NEFSetOut%put(cellR4,int(I),int(J))
          end do
        end do
!
      end select
!
      end function mqc_nef_set_outer_expand
!
!
!     PROCEDURE MQC_NEF_Set_Full_Collapse
!
!
!>    \brief <b> MQC_NEF_Set_Full_Collapse is equivalent to collapse
!>    followed by outerCollapse, condensed into a single transformation
!>    that never materializes the intermediate per-(I,J) collapsed grid.
!>    </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Full_Collapse produces the same result as calling
!>    collapse (per-cell rank reduction: Matrix->Vector or R4->Matrix)
!>    and then outerCollapse (folding the (I,J) grid into the collapsed
!>    element's own axes) in sequence, but does so in one pass: each
!>    cell's collapse() is computed locally, immediately folded into the
!>    final combined object, and discarded -- no full (NumI,NumJ)-sized
!>    intermediate MQC_NEF_Set of collapsed elements is ever built, so
!>    the index reordering that collapse and outerCollapse each perform
!>    is only ever walked once per cell, not twice.
!>
!>    If NEFSet holds Matrix elements, each cell collapses to a Vector
!>    (via spinType/swapOrder, as MQC_NEF_Matrix%collapse), then that
!>    Vector folds into the grid exactly as MQC_NEF_Set%outerCollapse
!>    does for a Vector-holding set (jointIJ, default true). If NEFSet
!>    holds R4 elements, each cell collapses to a Matrix (via
!>    spinType/spinTypeKet/Order, as MQC_NEF_R4%collapse), then that
!>    Matrix folds into the grid exactly as outerCollapse does for a
!>    Matrix-holding set (jointIJ, default false; jointAxis, default 1).
!>
!>    NEFSet must hold Matrix or R4 elements (not Vector, not 'open') --
!>    the same restriction collapse itself already has, since
!>    MQC_NEF_Vector has no collapse of its own.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The (I,J) grid to collapse per-cell and then outer-fold. Must
!>        hold Matrix or R4 elements.
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        Character(len=*),Optional
!>        Forwarded to the underlying Matrix/R4 collapse.
!>    \endverbatim
!>
!>    \param[in] swapOrder
!>    \verbatim
!>        Logical,Optional
!>        Forwarded to the underlying Matrix collapse only.
!>    \endverbatim
!>
!>    \param[in] spinTypeKet
!>    \verbatim
!>        Character(len=*),Optional
!>        Forwarded to the underlying R4 collapse only.
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Integer(kind=int64),Dimension(4),Optional
!>        Forwarded to the underlying R4 collapse only.
!>    \endverbatim
!>
!>    \param[in] jointIJ
!>    \verbatim
!>        jointIJ is Logical,Optional
!>        Forwarded to the outerCollapse step, applied to the collapsed
!>        element (Vector if NEFSet held Matrix; Matrix if NEFSet held
!>        R4). Default: true for a Matrix-holding NEFSet (collapsed
!>        result is Vector), false for an R4-holding NEFSet (collapsed
!>        result is Matrix) -- matching outerCollapse's own per-kind
!>        defaults, evaluated against the post-collapse kind.
!>    \endverbatim
!>
!>    \param[in] jointAxis
!>    \verbatim
!>        jointAxis is Integer,Optional
!>        Only used when NEFSet holds R4 elements and jointIJ is true
!>        (1=row/2=col of the collapsed Matrix). Default 1.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_full_collapse(NEFSet,spinType,swapOrder,spinTypeKet,Order,&
           jointIJ,jointAxis) result(NEFSetOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      character(len=*),intent(in),optional::spinType,spinTypeKet
      logical,intent(in),optional::swapOrder
      integer(kind=int64),dimension(4),intent(in),optional::Order
      logical,intent(in),optional::jointIJ
      integer,intent(in),optional::jointAxis
      type(mqc_nef_set)::NEFSetOut
      logical::myJointIJ
      integer::myJointAxis
      integer(kind=int64)::I,J,myNumI,myNumJ,G
      integer(kind=int64)::b,origNumSpace,outerIdx,newIdx,ka
      integer(kind=int64)::r,c,origNumSpace_row,origNumSpace_col
      integer(kind=int64)::newRowId,newColId,origId,newId
      type(mqc_nef_vector)::cellVector,NEFVecOut
      type(mqc_nef_matrix)::cellMatrix,NEFMatOut
      type(mqc_scf_eigenvalues),dimension(:),allocatable::New_Map_1,New_Map_2
      type(mqc_scalar)::tmpScalar
!
      myNumI = NEFSet%NumI
      myNumJ = NEFSet%NumJ
!
      select case(trim(NEFSet%Type()))
      case('vector')
        NEFSetOut = NEFSet%outerCollapse(jointIJ=jointIJ)
      case('matrix')
        if(present(jointIJ)) then
          myJointIJ = jointIJ
        else
          myJointIJ = .true.
        end if
        cellVector = NEFSet%Set_Matrix(1,1)%collapse(spinType=spinType,swapOrder=swapOrder)
        origNumSpace = cellVector%NumSpace
        if(myJointIJ) then
          allocate(New_Map_1(myNumI*myNumJ*origNumSpace))
          do outerIdx = 1, myNumI*myNumJ
            do b = 1, origNumSpace
              call New_Map_1((outerIdx-1)*origNumSpace+b)%init(&
                   cellVector%Size_vector_alpha(b),0,lengthBeta=cellVector%Size_vector_beta(b))
            end do
          end do
          call NEFVecOut%init(New_Map_1)
          do I = 1, myNumI
            do J = 1, myNumJ
              cellVector = NEFSet%Set_Matrix(I,J)%collapse(spinType=spinType,swapOrder=swapOrder)
              outerIdx = (I-1)*myNumJ+J
              do b = 1, origNumSpace
                if(cellVector%Zero_Flag(b)) cycle
                newIdx = (outerIdx-1)*origNumSpace+b
                NEFVecOut%Main_Vector(newIdx) = cellVector%Main_Vector(b)
                NEFVecOut%Zero_Flag(newIdx) = .false.
              end do
            end do
          end do
          call NEFSetOut%init(1,1)
          call NEFSetOut%put(NEFVecOut,1,1)
        else
          allocate(New_Map_1(myNumI*origNumSpace))
          do I = 1, myNumI
            do b = 1, origNumSpace
              call New_Map_1((I-1)*origNumSpace+b)%init(&
                   cellVector%Size_vector_alpha(b),0,lengthBeta=cellVector%Size_vector_beta(b))
            end do
          end do
          allocate(New_Map_2(myNumJ))
          do J = 1, myNumJ
            call New_Map_2(J)%init(1,0,lengthBeta=1)
          end do
          call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2)
          do I = 1, myNumI
            do J = 1, myNumJ
              cellVector = NEFSet%Set_Matrix(I,J)%collapse(spinType=spinType,swapOrder=swapOrder)
              do b = 1, origNumSpace
                if(cellVector%Zero_Flag(b)) cycle
                newRowId = (I-1)*origNumSpace+b
                newId = NEFMatOut%blockIndex([1,newRowId],[2,J])
                call NEFMatOut%Main_Matrix(newId)%init(&
                     cellVector%Size_vector_alpha(b),cellVector%Size_vector_beta(b),&
                     nAlpha2=1,nBeta2=1)
                NEFMatOut%Zero_Flag(newId) = .false.
                do ka = 1, cellVector%Size_vector_alpha(b)
                  tmpScalar = cellVector%Main_Vector(b)%at(ka,'alpha')
                  call NEFMatOut%Main_Matrix(newId)%put(tmpScalar,ka,1_int64,'aa')
                end do
                do ka = 1, cellVector%Size_vector_beta(b)
                  tmpScalar = cellVector%Main_Vector(b)%at(ka,'beta')
                  call NEFMatOut%Main_Matrix(newId)%put(tmpScalar,ka,1_int64,'bb')
                end do
              end do
            end do
          end do
          call NEFSetOut%init(1,1)
          call NEFSetOut%put(NEFMatOut,1,1)
        end if
      case('r4')
        if(present(jointIJ)) then
          myJointIJ = jointIJ
        else
          myJointIJ = .false.
        end if
        if(present(jointAxis)) then
          myJointAxis = jointAxis
        else
          myJointAxis = 1
        end if
        cellMatrix = NEFSet%Set_R4(1,1)%collapse(&
             spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
        origNumSpace_row = cellMatrix%NumSpace_row
        origNumSpace_col = cellMatrix%NumSpace_col
        if(.not.myJointIJ) then
          allocate(New_Map_1(myNumI*origNumSpace_row))
          do I = 1, myNumI
            cellMatrix = NEFSet%Set_R4(I,1)%collapse(&
                 spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
            do r = 1, origNumSpace_row
              call New_Map_1((I-1)*origNumSpace_row+r)%init(&
                   cellMatrix%Size_vector_alpha_row(r),0,lengthBeta=cellMatrix%Size_vector_beta_row(r))
            end do
          end do
          allocate(New_Map_2(myNumJ*origNumSpace_col))
          do J = 1, myNumJ
            cellMatrix = NEFSet%Set_R4(1,J)%collapse(&
                 spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
            do c = 1, origNumSpace_col
              call New_Map_2((J-1)*origNumSpace_col+c)%init(&
                   cellMatrix%Size_vector_alpha_col(c),0,lengthBeta=cellMatrix%Size_vector_beta_col(c))
            end do
          end do
        else
          if(myJointAxis.ne.1.and.myJointAxis.ne.2) call mqc_error_i(&
               'mqc_nef_set_full_collapse: jointAxis out-of-bounds for r4-collapsed matrix',6,&
               'jointAxis',myJointAxis)
          if(myJointAxis.eq.1) then
            allocate(New_Map_1(myNumI*myNumJ*origNumSpace_row))
            do I = 1, myNumI
              do J = 1, myNumJ
                cellMatrix = NEFSet%Set_R4(I,J)%collapse(&
                     spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
                G = (I-1)*myNumJ+J
                do r = 1, origNumSpace_row
                  call New_Map_1((G-1)*origNumSpace_row+r)%init(&
                       cellMatrix%Size_vector_alpha_row(r),0,lengthBeta=cellMatrix%Size_vector_beta_row(r))
                end do
              end do
            end do
            allocate(New_Map_2(origNumSpace_col))
            do c = 1, origNumSpace_col
              call New_Map_2(c)%init(&
                   cellMatrix%Size_vector_alpha_col(c),0,lengthBeta=cellMatrix%Size_vector_beta_col(c))
            end do
          else
            allocate(New_Map_1(origNumSpace_row))
            do r = 1, origNumSpace_row
              call New_Map_1(r)%init(&
                   cellMatrix%Size_vector_alpha_row(r),0,lengthBeta=cellMatrix%Size_vector_beta_row(r))
            end do
            allocate(New_Map_2(myNumI*myNumJ*origNumSpace_col))
            do I = 1, myNumI
              do J = 1, myNumJ
                cellMatrix = NEFSet%Set_R4(I,J)%collapse(&
                     spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
                G = (I-1)*myNumJ+J
                do c = 1, origNumSpace_col
                  call New_Map_2((G-1)*origNumSpace_col+c)%init(&
                       cellMatrix%Size_vector_alpha_col(c),0,lengthBeta=cellMatrix%Size_vector_beta_col(c))
                end do
              end do
            end do
          end if
        end if
        call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2)
        do I = 1, myNumI
          do J = 1, myNumJ
            cellMatrix = NEFSet%Set_R4(I,J)%collapse(&
                 spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
            G = (I-1)*myNumJ+J
            do r = 1, origNumSpace_row
              do c = 1, origNumSpace_col
                origId = cellMatrix%blockIndex([1,r],[2,c])
                if(cellMatrix%Zero_Flag(origId)) cycle
                if(.not.myJointIJ) then
                  newRowId = (I-1)*origNumSpace_row+r
                  newColId = (J-1)*origNumSpace_col+c
                else if(myJointAxis.eq.1) then
                  newRowId = (G-1)*origNumSpace_row+r
                  newColId = c
                else
                  newRowId = r
                  newColId = (G-1)*origNumSpace_col+c
                end if
                newId = NEFMatOut%blockIndex([1,newRowId],[2,newColId])
                NEFMatOut%Main_Matrix(newId) = cellMatrix%Main_Matrix(origId)
                NEFMatOut%Zero_Flag(newId) = .false.
              end do
            end do
          end do
        end do
        call NEFSetOut%init(1,1)
        call NEFSetOut%put(NEFMatOut,1,1)
      case default
        call mqc_error_A('mqc_nef_set_full_collapse: NEF_Set kind does not support '//&
             'fullCollapse (must be matrix or r4)',6,'Array_Type',NEFSet%Type())
      end select
!
      end function mqc_nef_set_full_collapse
!
!
!     PROCEDURE MQC_NEF_Set_Full_Expand
!
!
!>    \brief <b> MQC_NEF_Set_Full_Expand is the counterpart of
!>    MQC_NEF_Set%fullCollapse: outerExpand followed by expand. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Full_Expand inverts fullCollapse. Since fullCollapse
!>    itself is just collapse followed by outerCollapse fused into one
!>    pass purely for performance (per its own doc comment), its exact
!>    mathematical inverse is outerExpand followed by expand -- this
!>    function is a thin composition of the two (not a further fused
!>    single-pass implementation).
!>
!>    templateSet%Type() gives the kind the ORIGINAL (pre-fullCollapse)
!>    grid cells held: 'matrix' (whose fullCollapse intermediate was
!>    Vector) or 'r4' (whose fullCollapse intermediate was Matrix) --
!>    the same restriction collapse/fullCollapse themselves have.
!>
!>    The outerExpand stage's shape/size information (NumI, NumJ, and
!>    the intermediate kind's per-cell axis sizes) is only accepted here
!>    via templateSet (the exact NEF_Set originally passed to
!>    fullCollapse) -- unlike outerExpand itself, there is no explicit-
!>    size-array alternative here, since those sizes would describe the
!>    INTERMEDIATE collapsed kind, not the final target kind that
!>    map_matrix_1..4 describe, and reusing the same argument names for
!>    both would be ambiguous. Call outerExpand and expand separately if
!>    explicit (non-template) control over the intermediate stage is
!>    needed.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The (1,1)-wrapped combined element to unfold and expand
!>        (fullCollapse's result).
!>    \endverbatim
!>
!>    \param[in] templateSet
!>    \verbatim
!>        Type(MQC_NEF_Set)
!>        The exact NEF_Set originally passed to fullCollapse. Required
!>        -- its NumI/NumJ and per-cell axis sizes drive the outerExpand
!>        stage, and its Type() ('matrix' or 'r4') tells this function
!>        what kind the original grid cells held before fullCollapse
!>        (no separate origKind argument is needed).
!>    \endverbatim
!>
!>    \param[in] jointIJ
!>    \verbatim
!>        Logical,Optional
!>        Forwarded to the outerExpand stage. Must match the value used
!>        in the original fullCollapse call.
!>    \endverbatim
!>
!>    \param[in] jointAxis
!>    \verbatim
!>        Integer,Optional
!>        Forwarded to the outerExpand stage. Must match the value used
!>        in the original fullCollapse call.
!>    \endverbatim
!>
!>    \param[in] map_matrix_1,map_matrix_2,map_matrix_3,map_matrix_4
!>    \verbatim
!>        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Optional
!>        Forwarded to the final expand stage. map_matrix_1/2 required
!>        when templateSet holds Matrix elements; all four required when
!>        templateSet holds R4 elements.
!>    \endverbatim
!>
!>    \param[in] spinType,spinTypeKet,Order
!>    \verbatim
!>        Forwarded to the final expand stage.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_full_expand(NEFSet,templateSet,jointIJ,jointAxis,&
           map_matrix_1,map_matrix_2,map_matrix_3,map_matrix_4,&
           spinType,spinTypeKet,Order) result(NEFSetOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      type(mqc_nef_set),intent(in)::templateSet
      logical,intent(in),optional::jointIJ
      integer,intent(in),optional::jointAxis
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,intent(in),optional::&
           map_matrix_1,map_matrix_2,map_matrix_3,map_matrix_4
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      type(mqc_nef_set)::NEFSetOut
      type(mqc_nef_set)::intermediateSet
      character(len=64)::intermediateKind
!
      select case(trim(templateSet%Type()))
      case('matrix')
        intermediateKind = 'vector'
      case('r4')
        intermediateKind = 'matrix'
      case default
        call mqc_error_A('mqc_nef_set_full_expand: templateSet must hold '//&
             'matrix or r4 elements',6,'Array_Type',templateSet%Type())
      end select
!
      intermediateSet = mqc_nef_set_outer_expand(NEFSet,trim(intermediateKind),&
           templateSet=templateSet,jointIJ=jointIJ,jointAxis=jointAxis)
      NEFSetOut = mqc_nef_set_expand(intermediateSet,&
           map_matrix_1=map_matrix_1,map_matrix_2=map_matrix_2,&
           map_matrix_3=map_matrix_3,map_matrix_4=map_matrix_4,&
           spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
!
      end function mqc_nef_set_full_expand
!
!
!     PROCEDURE MQC_Scalar_NEFSet_Multiply
!
!
!>    \brief <b> MQC_Scalar_NEFSet_Multiply scales every (I,J) cell of a
!>    MQC_NEF_Set by a scalar. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Scalar_NEFSet_Multiply computes NEFSetOut(I,J) =
!>    scalarIn*NEFSet(I,J) for every (I,J) cell, returning a MQC_NEF_Set
!>    of the same shape and kind as NEFSetB. Zero-flagged cells are left
!>    unset in the output (consistent with every other NEF_Set grid
!>    operation). If scalarA is exactly zero, the output is left fully
!>    'open'/zero-flagged -- no per-cell multiply is performed.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] scalarA
!>    \verbatim
!>        scalarA is Class(*)
!>        Integer, real, complex, or mqc_scalar.
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        The (I,J) grid to scale.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_scalar_nefset_multiply(scalarA,NEFSetB) result(NEFSetOut)
!
      implicit none
      class(*),intent(in)::scalarA
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_scalar)::scalar
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      select type(scalarA)
      type is (integer)
        scalar = scalarA
      type is (real)
        scalar = scalarA
      type is (complex)
        scalar = scalarA
      type is (mqc_scalar)
        scalar = scalarA
      class default
        call mqc_error_I('scalarA type not defined in mqc_scalar_nefset_multiply',6)
      end select
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ)
!
      if(abs(scalar).ne.0.0) then
        select case(trim(NEFSetB%Type()))
        case('vector')
          do I = 1, NEFSetB%NumI
            do J = 1, NEFSetB%NumJ
              if(NEFSetB%isZero(I,J)) cycle
              call NEFSetOut%put(scalar*NEFSetB%atVector(I,J),I,J)
            end do
          end do
        case('matrix')
          do I = 1, NEFSetB%NumI
            do J = 1, NEFSetB%NumJ
              if(NEFSetB%isZero(I,J)) cycle
              call NEFSetOut%put(scalar*NEFSetB%atMatrix(I,J),I,J)
            end do
          end do
        case('r4')
          do I = 1, NEFSetB%NumI
            do J = 1, NEFSetB%NumJ
              if(NEFSetB%isZero(I,J)) cycle
              call NEFSetOut%put(scalar*NEFSetB%atR4(I,J),I,J)
            end do
          end do
        case default
          call mqc_error_A('mqc_scalar_nefset_multiply: NEF_Set is open '//&
               '(holds no elements)',6,'Array_Type',NEFSetB%Type())
        end select
      end if
!
      end function mqc_scalar_nefset_multiply
!
!
!     PROCEDURE MQC_NEFSet_Scalar_Multiply
!
!
!>    \brief <b> MQC_NEFSet_Scalar_Multiply scales every (I,J) cell of a
!>    MQC_NEF_Set by a scalar (reversed-operand form of
!>    MQC_Scalar_NEFSet_Multiply). </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    Identical to MQC_Scalar_NEFSet_Multiply, with NEFSetA*scalarB
!>    operand order instead of scalarA*NEFSetB.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>        The (I,J) grid to scale.
!>    \endverbatim
!>
!>    \param[in] scalarB
!>    \verbatim
!>        scalarB is Class(*)
!>        Integer, real, complex, or mqc_scalar.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_scalar_multiply(NEFSetA,scalarB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      class(*),intent(in)::scalarB
      type(mqc_nef_set)::NEFSetOut
!
      NEFSetOut = mqc_scalar_nefset_multiply(scalarB,NEFSetA)
!
      end function mqc_nefset_scalar_multiply
!
!
!     PROCEDURE MQC_NEFSet_Transpose
!
!
!>    \brief <b> MQC_NEFSet_Transpose transposes a MQC_NEF_Set's outer
!>    (I,J) grid and, for matrix-kind sets, each cell's own inner row/col
!>    axes. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    For matrix-kind sets, computes NEFSetOut(J,I) =
!>    Transpose(NEFSet(I,J)) -- both the outer grid position and each
!>    cell's own row/col axes are transposed. For vector-kind sets, only
!>    the outer grid position is swapped (NEFSetOut(J,I) = NEFSet(I,J));
!>    a vector cell has no row/col structure of its own to transpose.
!>    r4-kind sets are NOT yet supported here: neither MQC_TwoERIs nor
!>    MQC_R4Tensor (the types backing MQC_NEF_R4) has a plain,
!>    non-conjugated Transpose anywhere in mqcPack, only Dagger -- adding
!>    one would require mirroring mqc_twoeris_conjugate_transpose's full
!>    spin-block relabeling logic without conjugation, at both the
!>    MQC_R4Tensor and MQC_TwoERIs levels, which is separate follow-up
!>    work.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Type(MQC_NEF_Set)
!>        The (I,J) grid to transpose.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_transpose(NEFSet) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSet
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      call NEFSetOut%init(NEFSet%NumJ,NEFSet%NumI)
!
      select case(trim(NEFSet%Type()))
      case('vector')
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(NEFSet%atVector(I,J),J,I)
          end do
        end do
      case('matrix')
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(Transpose(NEFSet%atMatrix(I,J)),J,I)
          end do
        end do
      case('r4')
        call mqc_error('mqc_nefset_transpose: r4-kind NEF_Set transpose is '//&
             'not yet implemented (MQC_TwoERIs/MQC_R4Tensor have no plain '//&
             'Transpose, only Dagger)',6)
      case default
        call mqc_error_A('mqc_nefset_transpose: NEF_Set is open (holds no elements)',&
             6,'Array_Type',NEFSet%Type())
      end select
!
      end function mqc_nefset_transpose
!
!
!     PROCEDURE MQC_NEFSet_Conjugate_Transpose
!
!
!>    \brief <b> MQC_NEFSet_Conjugate_Transpose computes the conjugate
!>    transpose of a MQC_NEF_Set: the outer (I,J) grid is transposed and,
!>    for matrix- and r4-kind sets, each cell is replaced by its own
!>    Dagger. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    For matrix-kind sets: NEFSetOut(J,I) = Dagger(NEFSet(I,J)). For
!>    r4-kind sets: NEFSetOut(J,I) = Dagger(NEFSet(I,J)), reusing the
!>    existing MQC_NEFR4_Conjugate_Transpose directly -- unlike plain
!>    Transpose, Dagger is already fully implemented for r4-kind. For
!>    vector-kind sets, only the outer grid position is swapped
!>    (NEFSetOut(J,I) = NEFSet(I,J)) -- the cell's own values are NOT
!>    conjugated, since MQC_NEF_Vector has no Conjg/Dagger primitive
!>    anywhere in mqcPack to call.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Type(MQC_NEF_Set)
!>        The (I,J) grid to conjugate-transpose.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_conjugate_transpose(NEFSet) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSet
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      call NEFSetOut%init(NEFSet%NumJ,NEFSet%NumI)
!
      select case(trim(NEFSet%Type()))
      case('vector')
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(NEFSet%atVector(I,J),J,I)
          end do
        end do
      case('matrix')
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(Dagger(NEFSet%atMatrix(I,J)),J,I)
          end do
        end do
      case('r4')
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(Dagger(NEFSet%atR4(I,J)),J,I)
          end do
        end do
      case default
        call mqc_error_A('mqc_nefset_conjugate_transpose: NEF_Set is open '//&
             '(holds no elements)',6,'Array_Type',NEFSet%Type())
      end select
!
      end function mqc_nefset_conjugate_transpose
!
!
!     PROCEDURE MQC_NEFSet_NEFVector_Sum
!
!>    \brief <b> MQC_NEFSet_NEFVector_Sum adds a plain MQC_NEF_Vector to
!>    every (I,J) cell of a vector-kind MQC_NEF_Set. </b>
!
      function mqc_nefset_nefvector_sum(NEFSetA,NEFVecB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_nef_vector),intent(in)::NEFVecB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'vector') call mqc_error_A(&
           'mqc_nefset_nefvector_sum: NEF_Set must hold vector elements (or be open)',&
           6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ)
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call NEFSetOut%put(NEFVecB,I,J)
          else
            call NEFSetOut%put(NEFSetA%Set_Vector(I,J)+NEFVecB,I,J)
          end if
        end do
      end do
!
      end function mqc_nefset_nefvector_sum
!
!
!     PROCEDURE MQC_NEFVector_NEFSet_Sum
!
!>    \brief <b> MQC_NEFVector_NEFSet_Sum adds a plain MQC_NEF_Vector to
!>    every (I,J) cell of a vector-kind MQC_NEF_Set (reversed-operand
!>    form of MQC_NEFSet_NEFVector_Sum). </b>
!
      function mqc_nefvector_nefset_sum(NEFVecA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecA
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'vector') call mqc_error_A(&
           'mqc_nefvector_nefset_sum: NEF_Set must hold vector elements (or be open)',&
           6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ)
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call NEFSetOut%put(NEFVecA,I,J)
          else
            call NEFSetOut%put(NEFVecA+NEFSetB%Set_Vector(I,J),I,J)
          end if
        end do
      end do
!
      end function mqc_nefvector_nefset_sum
!
!
!     PROCEDURE MQC_NEFSet_NEFVector_Difference
!
!>    \brief <b> MQC_NEFSet_NEFVector_Difference subtracts a plain
!>    MQC_NEF_Vector from every (I,J) cell of a vector-kind
!>    MQC_NEF_Set. </b>
!
      function mqc_nefset_nefvector_difference(NEFSetA,NEFVecB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_nef_vector),intent(in)::NEFVecB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'vector') call mqc_error_A(&
           'mqc_nefset_nefvector_difference: NEF_Set must hold vector elements '//&
           '(or be open)',6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ)
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call NEFSetOut%put(NEFVecB*(-1.0),I,J)
          else
            call NEFSetOut%put(NEFSetA%Set_Vector(I,J)-NEFVecB,I,J)
          end if
        end do
      end do
!
      end function mqc_nefset_nefvector_difference
!
!
!     PROCEDURE MQC_NEFVector_NEFSet_Difference
!
!>    \brief <b> MQC_NEFVector_NEFSet_Difference subtracts every (I,J)
!>    cell of a vector-kind MQC_NEF_Set from a plain MQC_NEF_Vector
!>    (reversed-operand form of MQC_NEFSet_NEFVector_Difference). </b>
!
      function mqc_nefvector_nefset_difference(NEFVecA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_vector),intent(in)::NEFVecA
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'vector') call mqc_error_A(&
           'mqc_nefvector_nefset_difference: NEF_Set must hold vector elements '//&
           '(or be open)',6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ)
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call NEFSetOut%put(NEFVecA,I,J)
          else
            call NEFSetOut%put(NEFVecA-NEFSetB%Set_Vector(I,J),I,J)
          end if
        end do
      end do
!
      end function mqc_nefvector_nefset_difference
!
!
!     PROCEDURE MQC_NEFSet_NEFMatrix_Sum
!
!>    \brief <b> MQC_NEFSet_NEFMatrix_Sum adds a plain MQC_NEF_Matrix to
!>    every (I,J) cell of a matrix-kind MQC_NEF_Set. </b>
!
      function mqc_nefset_nefmatrix_sum(NEFSetA,NEFMatB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefset_nefmatrix_sum: NEF_Set must hold matrix elements (or be open)',&
           6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ)
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call NEFSetOut%put(NEFMatB,I,J)
          else
            call NEFSetOut%put(NEFSetA%Set_Matrix(I,J)+NEFMatB,I,J)
          end if
        end do
      end do
!
      end function mqc_nefset_nefmatrix_sum
!
!
!     PROCEDURE MQC_NEFMatrix_NEFSet_Sum
!
!>    \brief <b> MQC_NEFMatrix_NEFSet_Sum adds a plain MQC_NEF_Matrix to
!>    every (I,J) cell of a matrix-kind MQC_NEF_Set (reversed-operand
!>    form of MQC_NEFSet_NEFMatrix_Sum). </b>
!
      function mqc_nefmatrix_nefset_sum(NEFMatA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefmatrix_nefset_sum: NEF_Set must hold matrix elements (or be open)',&
           6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ)
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call NEFSetOut%put(NEFMatA,I,J)
          else
            call NEFSetOut%put(NEFMatA+NEFSetB%Set_Matrix(I,J),I,J)
          end if
        end do
      end do
!
      end function mqc_nefmatrix_nefset_sum
!
!
!     PROCEDURE MQC_NEFSet_NEFMatrix_Difference
!
!>    \brief <b> MQC_NEFSet_NEFMatrix_Difference subtracts a plain
!>    MQC_NEF_Matrix from every (I,J) cell of a matrix-kind
!>    MQC_NEF_Set. </b>
!
      function mqc_nefset_nefmatrix_difference(NEFSetA,NEFMatB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefset_nefmatrix_difference: NEF_Set must hold matrix elements '//&
           '(or be open)',6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ)
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call NEFSetOut%put(NEFMatB*(-1.0),I,J)
          else
            call NEFSetOut%put(NEFSetA%Set_Matrix(I,J)-NEFMatB,I,J)
          end if
        end do
      end do
!
      end function mqc_nefset_nefmatrix_difference
!
!
!     PROCEDURE MQC_NEFMatrix_NEFSet_Difference
!
!>    \brief <b> MQC_NEFMatrix_NEFSet_Difference subtracts every (I,J)
!>    cell of a matrix-kind MQC_NEF_Set from a plain MQC_NEF_Matrix
!>    (reversed-operand form of MQC_NEFSet_NEFMatrix_Difference). </b>
!
      function mqc_nefmatrix_nefset_difference(NEFMatA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefmatrix_nefset_difference: NEF_Set must hold matrix elements '//&
           '(or be open)',6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ)
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call NEFSetOut%put(NEFMatA,I,J)
          else
            call NEFSetOut%put(NEFMatA-NEFSetB%Set_Matrix(I,J),I,J)
          end if
        end do
      end do
!
      end function mqc_nefmatrix_nefset_difference
!
!
!     PROCEDURE MQC_NEFSet_NEFR4_Sum
!
!>    \brief <b> MQC_NEFSet_NEFR4_Sum adds a plain MQC_NEF_R4 to every
!>    (I,J) cell of an r4-kind MQC_NEF_Set. </b>
!
      function mqc_nefset_nefr4_sum(NEFSetA,NEFR4B) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_nef_r4),intent(in)::NEFR4B
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'r4') call mqc_error_A(&
           'mqc_nefset_nefr4_sum: NEF_Set must hold r4 elements (or be open)',&
           6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ)
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call NEFSetOut%put(NEFR4B,I,J)
          else
            call NEFSetOut%put(NEFSetA%Set_R4(I,J)+NEFR4B,I,J)
          end if
        end do
      end do
!
      end function mqc_nefset_nefr4_sum
!
!
!     PROCEDURE MQC_NEFR4_NEFSet_Sum
!
!>    \brief <b> MQC_NEFR4_NEFSet_Sum adds a plain MQC_NEF_R4 to every
!>    (I,J) cell of an r4-kind MQC_NEF_Set (reversed-operand form of
!>    MQC_NEFSet_NEFR4_Sum). </b>
!
      function mqc_nefr4_nefset_sum(NEFR4A,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'r4') call mqc_error_A(&
           'mqc_nefr4_nefset_sum: NEF_Set must hold r4 elements (or be open)',&
           6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ)
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call NEFSetOut%put(NEFR4A,I,J)
          else
            call NEFSetOut%put(NEFR4A+NEFSetB%Set_R4(I,J),I,J)
          end if
        end do
      end do
!
      end function mqc_nefr4_nefset_sum
!
!
!     PROCEDURE MQC_NEFSet_NEFR4_Difference
!
!>    \brief <b> MQC_NEFSet_NEFR4_Difference subtracts a plain
!>    MQC_NEF_R4 from every (I,J) cell of an r4-kind MQC_NEF_Set. </b>
!
      function mqc_nefset_nefr4_difference(NEFSetA,NEFR4B) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_nef_r4),intent(in)::NEFR4B
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'r4') call mqc_error_A(&
           'mqc_nefset_nefr4_difference: NEF_Set must hold r4 elements (or be open)',&
           6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ)
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call NEFSetOut%put(NEFR4B*(-1.0),I,J)
          else
            call NEFSetOut%put(NEFSetA%Set_R4(I,J)-NEFR4B,I,J)
          end if
        end do
      end do
!
      end function mqc_nefset_nefr4_difference
!
!
!     PROCEDURE MQC_NEFR4_NEFSet_Difference
!
!>    \brief <b> MQC_NEFR4_NEFSet_Difference subtracts every (I,J) cell
!>    of an r4-kind MQC_NEF_Set from a plain MQC_NEF_R4 (reversed-
!>    operand form of MQC_NEFSet_NEFR4_Difference). </b>
!
      function mqc_nefr4_nefset_difference(NEFR4A,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_r4),intent(in)::NEFR4A
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'r4') call mqc_error_A(&
           'mqc_nefr4_nefset_difference: NEF_Set must hold r4 elements (or be open)',&
           6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ)
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call NEFSetOut%put(NEFR4A,I,J)
          else
            call NEFSetOut%put(NEFR4A-NEFSetB%Set_R4(I,J),I,J)
          end if
        end do
      end do
!
      end function mqc_nefr4_nefset_difference
!
!
!     PROCEDURE MQC_NEFSet_NEFSet_Sum
!
!>    \brief <b> MQC_NEFSet_NEFSet_Sum adds two identically-shaped,
!>    same-kind MQC_NEF_Set grids cell-by-cell. </b>
!
      function mqc_nefset_nefset_sum(NEFSetA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA,NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(trim(NEFSetA%Type()).ne.trim(NEFSetB%Type())) call mqc_error_A(&
           'mqc_nefset_nefset_sum: NEF_Sets must hold the same kind',&
           6,'NEFSetA%Type',NEFSetA%Type())
      if(NEFSetA%NumI.ne.NEFSetB%NumI.or.NEFSetA%NumJ.ne.NEFSetB%NumJ) call mqc_error_i(&
           'mqc_nefset_nefset_sum: NEFSetA and NEFSetB must have identical '//&
           'grid shape',6,'NEFSetA%NumI',int(NEFSetA%NumI,kind=int64),&
           'NEFSetB%NumI',int(NEFSetB%NumI,kind=int64))
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ)
!
      select case(trim(NEFSetA%Type()))
      case('vector')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call NEFSetOut%put(NEFSetB%Set_Vector(I,J),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call NEFSetOut%put(NEFSetA%Set_Vector(I,J),I,J)
            else
              call NEFSetOut%put(NEFSetA%Set_Vector(I,J)+NEFSetB%Set_Vector(I,J),I,J)
            end if
          end do
        end do
      case('matrix')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call NEFSetOut%put(NEFSetB%Set_Matrix(I,J),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call NEFSetOut%put(NEFSetA%Set_Matrix(I,J),I,J)
            else
              call NEFSetOut%put(NEFSetA%Set_Matrix(I,J)+NEFSetB%Set_Matrix(I,J),I,J)
            end if
          end do
        end do
      case('r4')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call NEFSetOut%put(NEFSetB%Set_R4(I,J),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call NEFSetOut%put(NEFSetA%Set_R4(I,J),I,J)
            else
              call NEFSetOut%put(NEFSetA%Set_R4(I,J)+NEFSetB%Set_R4(I,J),I,J)
            end if
          end do
        end do
      case default
        call mqc_error_A('mqc_nefset_nefset_sum: NEF_Set is open (holds no elements)',&
             6,'Array_Type',NEFSetA%Type())
      end select
!
      end function mqc_nefset_nefset_sum
!
!
!     PROCEDURE MQC_NEFSet_NEFSet_Difference
!
!>    \brief <b> MQC_NEFSet_NEFSet_Difference subtracts two identically-
!>    shaped, same-kind MQC_NEF_Set grids cell-by-cell. </b>
!
      function mqc_nefset_nefset_difference(NEFSetA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA,NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      if(trim(NEFSetA%Type()).ne.trim(NEFSetB%Type())) call mqc_error_A(&
           'mqc_nefset_nefset_difference: NEF_Sets must hold the same kind',&
           6,'NEFSetA%Type',NEFSetA%Type())
      if(NEFSetA%NumI.ne.NEFSetB%NumI.or.NEFSetA%NumJ.ne.NEFSetB%NumJ) call mqc_error_i(&
           'mqc_nefset_nefset_difference: NEFSetA and NEFSetB must have '//&
           'identical grid shape',6,'NEFSetA%NumI',int(NEFSetA%NumI,kind=int64),&
           'NEFSetB%NumI',int(NEFSetB%NumI,kind=int64))
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ)
!
      select case(trim(NEFSetA%Type()))
      case('vector')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call NEFSetOut%put(NEFSetB%Set_Vector(I,J)*(-1.0),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call NEFSetOut%put(NEFSetA%Set_Vector(I,J),I,J)
            else
              call NEFSetOut%put(NEFSetA%Set_Vector(I,J)-NEFSetB%Set_Vector(I,J),I,J)
            end if
          end do
        end do
      case('matrix')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call NEFSetOut%put(NEFSetB%Set_Matrix(I,J)*(-1.0),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call NEFSetOut%put(NEFSetA%Set_Matrix(I,J),I,J)
            else
              call NEFSetOut%put(NEFSetA%Set_Matrix(I,J)-NEFSetB%Set_Matrix(I,J),I,J)
            end if
          end do
        end do
      case('r4')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call NEFSetOut%put(NEFSetB%Set_R4(I,J)*(-1.0),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call NEFSetOut%put(NEFSetA%Set_R4(I,J),I,J)
            else
              call NEFSetOut%put(NEFSetA%Set_R4(I,J)-NEFSetB%Set_R4(I,J),I,J)
            end if
          end do
        end do
      case default
        call mqc_error_A('mqc_nefset_nefset_difference: NEF_Set is open '//&
             '(holds no elements)',6,'Array_Type',NEFSetA%Type())
      end select
!
      end function mqc_nefset_nefset_difference
!
!
!     PROCEDURE MQC_NEF_Matrix_Dimention_Expand
!
!
!>    \brief <b> MQC_NEF_Matrix_Dimention_Expand is used to return an
!>    expanded NEF R4 object based on a NEF Matrix object -- the reverse
!>    of MQC_NEF_R4%collapse. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Dimention_Expand is used to return an expanded NEF
!>    R4 object based on a NEF Matrix object, given the four original
!>    per-axis maps (map_matrix_1..4, one per physical R4 axis role,
!>    matching whatever Order was used when the Matrix was built by
!>    MQC_NEF_R4%collapse) -- the same relationship
!>    MQC_NEF_Vector_Dimension_Expand has to MQC_NEF_Matrix%collapse, one
!>    rank higher. Mirrors mqc_nef_r4_dimention_collapse's exact
!>    singlet/triplet/all case structure in reverse: for every element
!>    collapse would have read from R4 and written to Matrix, expand
!>    reads the same position from Matrix and writes it to R4, using the
!>    same collapseIndex formulas (now evaluated on the R4 object being
!>    built, since its Map/size metadata is already set from the
!>    supplied maps before the main loop runs).
!>
!>    Depends on MQC_NEF_R4%put (mqc_nef_r4_put), which itself depends on
!>    MQC_TwoERIs gaining a put(element,i,j,k,l,spinBlockIn) counterpart
!>    to its existing at, mirroring MQC_SCF_Integral's put/at pair.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is type(mqc_nef_matrix)
!>        NEF Matrix object to be expanded.
!>    \endverbatim
!>
!>    \param[in] map_matrix_1
!>    \verbatim
!>        map_matrix_1 is type(mqc_scf_eigenvalues),dimension(:),allocatable
!>        Original per-space map for the role that occupies physical
!>        axis Order(1) (default physical axis 1).
!>    \endverbatim
!>
!>    \param[in] map_matrix_2
!>    \verbatim
!>        map_matrix_2 is type(mqc_scf_eigenvalues),dimension(:),allocatable
!>        Original per-space map for the role that occupies physical
!>        axis Order(2) (default physical axis 2).
!>    \endverbatim
!>
!>    \param[in] map_matrix_3
!>    \verbatim
!>        map_matrix_3 is type(mqc_scf_eigenvalues),dimension(:),allocatable
!>        Original per-space map for the role that occupies physical
!>        axis Order(3) (default physical axis 3).
!>    \endverbatim
!>
!>    \param[in] map_matrix_4
!>    \verbatim
!>        map_matrix_4 is type(mqc_scf_eigenvalues),dimension(:),allocatable
!>        Original per-space map for the role that occupies physical
!>        axis Order(4) (default physical axis 4).
!>    \endverbatim
!>
!>    \param[in] spinType
!>    \verbatim
!>        spinType is Character(Len=*),Optional
!>        Forwarded exactly as in MQC_NEF_R4%collapse: 'singlet','triplet',
!>        or 'all'/'full' (default).
!>    \endverbatim
!>
!>    \param[in] spinTypeKet
!>    \verbatim
!>        spinTypeKet is Character(Len=*),Optional
!>        Forwarded exactly as in MQC_NEF_R4%collapse.
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Order is Integer(kind=int64),Dimension(4),Optional
!>        Forwarded exactly as in MQC_NEF_R4%collapse: Order=[i,j,k,l] =>
!>        i->j = row-pair, k->l = col-pair. Default [1,2,3,4].
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_dimention_expand(NEFMat,map_matrix_1,map_matrix_2,&
           map_matrix_3,map_matrix_4,spinType,spinTypeKet,Order) result(NEFR4Out)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,intent(in)::&
           map_matrix_1,map_matrix_2,map_matrix_3,map_matrix_4
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      type(mqc_nef_r4)::NEFR4Out
      character(len=64)::mySpinType,mySpinTypeKet
      integer(kind=int64)::ia,ja,ka,la
      integer(kind=int64)::dimi,dimj,dimk,diml,i,j,k,l,idij,idlk,finalid,nef_bidx
      integer(kind=int64)::nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4
!
      if(present(spinType)) then
        mySpinType = spinType
        if(present(spinTypeKet)) then
          mySpinTypeKet = spinTypeKet
        else
          mySpinTypeKet = spinType
        end if
      else
        if(present(spinTypeKet)) then
          mySpinTypeKet = spinTypeKet
          mySpinType = spinTypeKet
        else
          mySpinTypeKet = 'all'
          mySpinType = 'all'
        end if
      end if
!
      if(present(Order)) then
        ia=Order(1)
        ja=Order(2)
        ka=Order(3)
        la=Order(4)
        if (ia.eq.ja) Call MQC_Error_I(&
               'Duplicated index in mqc_nef_matrix_dimention_expand',&
                                        6,'ia', ia,'ja',ja)
        if (ia.eq.ka) Call MQC_Error_I(&
               'Duplicated index in mqc_nef_matrix_dimention_expand',&
                                        6,'ia', ia,'ka',ka)
        if (ia.eq.la) Call MQC_Error_I(&
               'Duplicated index in mqc_nef_matrix_dimention_expand',&
                                        6,'ia', ia,'la',la)
        if (ja.eq.ka) Call MQC_Error_I(&
               'Duplicated index in mqc_nef_matrix_dimention_expand',&
                                        6,'ja', ja,'ka',ka)
        if (ja.eq.la) Call MQC_Error_I(&
               'Duplicated index in mqc_nef_matrix_dimention_expand',&
                                        6,'ja', ja,'la',la)
        if (ka.eq.la) Call MQC_Error_I(&
               'Duplicated index in mqc_nef_matrix_dimention_expand',&
                                        6,'ka', ka,'la',la)
      else
        ia=1
        ja=2
        ka=3
        la=4
      endif
!
!     Role (1=ia,2=ja,3=ka,4=la) -> physical axis assignment of the four
!     supplied maps, mirroring collapse's own Order convention.
      if(ia.eq.1) then
        NEFR4Out%Map_Matrix_1 = map_matrix_1
      else if(ja.eq.1) then
        NEFR4Out%Map_Matrix_1 = map_matrix_2
      else if(ka.eq.1) then
        NEFR4Out%Map_Matrix_1 = map_matrix_3
      else
        NEFR4Out%Map_Matrix_1 = map_matrix_4
      end if
      if(ia.eq.2) then
        NEFR4Out%Map_Matrix_2 = map_matrix_1
      else if(ja.eq.2) then
        NEFR4Out%Map_Matrix_2 = map_matrix_2
      else if(ka.eq.2) then
        NEFR4Out%Map_Matrix_2 = map_matrix_3
      else
        NEFR4Out%Map_Matrix_2 = map_matrix_4
      end if
      if(ia.eq.3) then
        NEFR4Out%Map_Matrix_3 = map_matrix_1
      else if(ja.eq.3) then
        NEFR4Out%Map_Matrix_3 = map_matrix_2
      else if(ka.eq.3) then
        NEFR4Out%Map_Matrix_3 = map_matrix_3
      else
        NEFR4Out%Map_Matrix_3 = map_matrix_4
      end if
      if(ia.eq.4) then
        NEFR4Out%Map_Matrix_4 = map_matrix_1
      else if(ja.eq.4) then
        NEFR4Out%Map_Matrix_4 = map_matrix_2
      else if(ka.eq.4) then
        NEFR4Out%Map_Matrix_4 = map_matrix_3
      else
        NEFR4Out%Map_Matrix_4 = map_matrix_4
      end if
!
      NEFR4Out%NumSpace_1 = size(NEFR4Out%Map_Matrix_1)
      NEFR4Out%NumSpace_2 = size(NEFR4Out%Map_Matrix_2)
      NEFR4Out%NumSpace_3 = size(NEFR4Out%Map_Matrix_3)
      NEFR4Out%NumSpace_4 = size(NEFR4Out%Map_Matrix_4)
      allocate(NEFR4Out%Size_vector_alpha_1(NEFR4Out%NumSpace_1))
      allocate(NEFR4Out%Size_vector_beta_1(NEFR4Out%NumSpace_1))
      allocate(NEFR4Out%Size_vector_alpha_2(NEFR4Out%NumSpace_2))
      allocate(NEFR4Out%Size_vector_beta_2(NEFR4Out%NumSpace_2))
      allocate(NEFR4Out%Size_vector_alpha_3(NEFR4Out%NumSpace_3))
      allocate(NEFR4Out%Size_vector_beta_3(NEFR4Out%NumSpace_3))
      allocate(NEFR4Out%Size_vector_alpha_4(NEFR4Out%NumSpace_4))
      allocate(NEFR4Out%Size_vector_beta_4(NEFR4Out%NumSpace_4))
      do i = 1, NEFR4Out%NumSpace_1
        NEFR4Out%Size_vector_alpha_1(i) = Size(NEFR4Out%Map_Matrix_1(i),'alpha')
        NEFR4Out%Size_vector_beta_1(i) = Size(NEFR4Out%Map_Matrix_1(i),'beta')
      end do
      do i = 1, NEFR4Out%NumSpace_2
        NEFR4Out%Size_vector_alpha_2(i) = Size(NEFR4Out%Map_Matrix_2(i),'alpha')
        NEFR4Out%Size_vector_beta_2(i) = Size(NEFR4Out%Map_Matrix_2(i),'beta')
      end do
      do i = 1, NEFR4Out%NumSpace_3
        NEFR4Out%Size_vector_alpha_3(i) = Size(NEFR4Out%Map_Matrix_3(i),'alpha')
        NEFR4Out%Size_vector_beta_3(i) = Size(NEFR4Out%Map_Matrix_3(i),'beta')
      end do
      do i = 1, NEFR4Out%NumSpace_4
        NEFR4Out%Size_vector_alpha_4(i) = Size(NEFR4Out%Map_Matrix_4(i),'alpha')
        NEFR4Out%Size_vector_beta_4(i) = Size(NEFR4Out%Map_Matrix_4(i),'beta')
      end do
      allocate(NEFR4Out%Main_R4(NEFR4Out%NumSpace_1*NEFR4Out%NumSpace_2*&
                                 NEFR4Out%NumSpace_3*NEFR4Out%NumSpace_4))
      allocate(NEFR4Out%Zero_Flag(size(NEFR4Out%Main_R4)))
      NEFR4Out%Zero_Flag = .true.
!
      select case(mySpinType)
      case ('singlet')
        do dimi = 1, NEFR4Out%blockSize(ia)
         do dimj = 1, NEFR4Out%blockSize(ja)
          idij = NEFR4Out%blockSize(ja)*(dimi-1) + dimj
          do dimk = 1, NEFR4Out%blockSize(ka)
           do diml = 1, NEFR4Out%blockSize(la)
            idlk = NEFR4Out%blockSize(la)*(dimk-1) + diml
            finalid = NEFMat%blockIndex([1,idij],[2,idlk])
            nef_bidx = NEFR4Out%blockIndex([ia,dimi],[ja,dimj],[la,diml],[ka,dimk])
            if (NEFMat%Zero_Flag(finalid)) cycle
            if(ia.eq.1) then
              nA1=NEFR4Out%blockSize(ia,dimi,'alpha'); nB1=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.1) then
              nA1=NEFR4Out%blockSize(ja,dimj,'alpha'); nB1=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.1) then
              nA1=NEFR4Out%blockSize(ka,dimk,'alpha'); nB1=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA1=NEFR4Out%blockSize(la,diml,'alpha'); nB1=NEFR4Out%blockSize(la,diml,'beta')
            end if
            if(ia.eq.2) then
              nA2=NEFR4Out%blockSize(ia,dimi,'alpha'); nB2=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.2) then
              nA2=NEFR4Out%blockSize(ja,dimj,'alpha'); nB2=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.2) then
              nA2=NEFR4Out%blockSize(ka,dimk,'alpha'); nB2=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA2=NEFR4Out%blockSize(la,diml,'alpha'); nB2=NEFR4Out%blockSize(la,diml,'beta')
            end if
            if(ia.eq.3) then
              nA3=NEFR4Out%blockSize(ia,dimi,'alpha'); nB3=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.3) then
              nA3=NEFR4Out%blockSize(ja,dimj,'alpha'); nB3=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.3) then
              nA3=NEFR4Out%blockSize(ka,dimk,'alpha'); nB3=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA3=NEFR4Out%blockSize(la,diml,'alpha'); nB3=NEFR4Out%blockSize(la,diml,'beta')
            end if
            if(ia.eq.4) then
              nA4=NEFR4Out%blockSize(ia,dimi,'alpha'); nB4=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.4) then
              nA4=NEFR4Out%blockSize(ja,dimj,'alpha'); nB4=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.4) then
              nA4=NEFR4Out%blockSize(ka,dimk,'alpha'); nB4=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA4=NEFR4Out%blockSize(la,diml,'alpha'); nB4=NEFR4Out%blockSize(la,diml,'beta')
            end if
            call NEFR4Out%Main_R4(nef_bidx)%init(nA1,nB1,nAlpha2=nA2,nBeta2=nB2,&
                 nAlpha3=nA3,nBeta3=nB3,nAlpha4=nA4,nBeta4=nB4)
            NEFR4Out%Zero_Flag(nef_bidx) = .false.
            do i = 1, NEFR4Out%blockSize(ia,dimi,'alpha')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'alpha')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'singlet')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'singlet')),&
                       'alpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,-1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'singlet')),&
                       int(NEFR4Out%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'singlet')),&
                       'betaalpha'))
               end do
              end do
             end do
            end do
            do i = 1, NEFR4Out%blockSize(ia,dimi,'beta')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'beta')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'singlet')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'singlet')),&
                       'alphabeta'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,-1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'singlet')),&
                       int(NEFR4Out%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'singlet')),&
                       'beta'))
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
      case ('triplet')
        do dimi = 1, NEFR4Out%blockSize(ia)
         do dimj = 1, NEFR4Out%blockSize(ja)
          idij = NEFR4Out%blockIndex([1,1],[2,1],[3,dimi],[4,dimj])
          do dimk = 1, NEFR4Out%blockSize(ka)
           do diml = 1, NEFR4Out%blockSize(la)
            idlk = NEFR4Out%blockIndex([1,1],[2,1],[3,diml],[4,dimk])
            finalid = NEFMat%blockIndex([1,idij],[2,idlk])
            nef_bidx = NEFR4Out%blockIndex([ia,dimi],[ja,dimj],[la,diml],[ka,dimk])
            if (NEFMat%Zero_Flag(finalid)) cycle
            if(ia.eq.1) then
              nA1=NEFR4Out%blockSize(ia,dimi,'alpha'); nB1=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.1) then
              nA1=NEFR4Out%blockSize(ja,dimj,'alpha'); nB1=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.1) then
              nA1=NEFR4Out%blockSize(ka,dimk,'alpha'); nB1=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA1=NEFR4Out%blockSize(la,diml,'alpha'); nB1=NEFR4Out%blockSize(la,diml,'beta')
            end if
            if(ia.eq.2) then
              nA2=NEFR4Out%blockSize(ia,dimi,'alpha'); nB2=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.2) then
              nA2=NEFR4Out%blockSize(ja,dimj,'alpha'); nB2=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.2) then
              nA2=NEFR4Out%blockSize(ka,dimk,'alpha'); nB2=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA2=NEFR4Out%blockSize(la,diml,'alpha'); nB2=NEFR4Out%blockSize(la,diml,'beta')
            end if
            if(ia.eq.3) then
              nA3=NEFR4Out%blockSize(ia,dimi,'alpha'); nB3=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.3) then
              nA3=NEFR4Out%blockSize(ja,dimj,'alpha'); nB3=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.3) then
              nA3=NEFR4Out%blockSize(ka,dimk,'alpha'); nB3=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA3=NEFR4Out%blockSize(la,diml,'alpha'); nB3=NEFR4Out%blockSize(la,diml,'beta')
            end if
            if(ia.eq.4) then
              nA4=NEFR4Out%blockSize(ia,dimi,'alpha'); nB4=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.4) then
              nA4=NEFR4Out%blockSize(ja,dimj,'alpha'); nB4=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.4) then
              nA4=NEFR4Out%blockSize(ka,dimk,'alpha'); nB4=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA4=NEFR4Out%blockSize(la,diml,'alpha'); nB4=NEFR4Out%blockSize(la,diml,'beta')
            end if
            call NEFR4Out%Main_R4(nef_bidx)%init(nA1,nB1,nAlpha2=nA2,nBeta2=nB2,&
                 nAlpha3=nA3,nBeta3=nB3,nAlpha4=nA4,nBeta4=nB4)
            NEFR4Out%Zero_Flag(nef_bidx) = .false.
            do i = 1, NEFR4Out%blockSize(ia,dimi,'alpha')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'beta')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,-1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'triplet')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'triplet')),&
                       'alpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'triplet')),&
                       int(NEFR4Out%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'triplet')),&
                       'alphabeta'))
               end do
              end do
             end do
            end do
            do i = 1, NEFR4Out%blockSize(ia,dimi,'beta')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'alpha')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,-1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'triplet')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'triplet')),&
                       'betaalpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'triplet')),&
                       int(NEFR4Out%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'triplet')),&
                       'beta'))
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
      case ('all')
        do dimi = 1, NEFR4Out%blockSize(ia)
         do dimj = 1, NEFR4Out%blockSize(ja)
          idij = NEFR4Out%blockIndex([1,1],[2,1],[3,dimi],[4,dimj])
          do dimk = 1, NEFR4Out%blockSize(ka)
           do diml = 1, NEFR4Out%blockSize(la)
            idlk = NEFR4Out%blockIndex([1,1],[2,1],[3,diml],[4,dimk])
            finalid = NEFMat%blockIndex([1,idij],[2,idlk])
            nef_bidx = NEFR4Out%blockIndex([ia,dimi],[ja,dimj],[la,diml],[ka,dimk])
            if (NEFMat%Zero_Flag(finalid)) cycle
            if(ia.eq.1) then
              nA1=NEFR4Out%blockSize(ia,dimi,'alpha'); nB1=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.1) then
              nA1=NEFR4Out%blockSize(ja,dimj,'alpha'); nB1=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.1) then
              nA1=NEFR4Out%blockSize(ka,dimk,'alpha'); nB1=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA1=NEFR4Out%blockSize(la,diml,'alpha'); nB1=NEFR4Out%blockSize(la,diml,'beta')
            end if
            if(ia.eq.2) then
              nA2=NEFR4Out%blockSize(ia,dimi,'alpha'); nB2=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.2) then
              nA2=NEFR4Out%blockSize(ja,dimj,'alpha'); nB2=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.2) then
              nA2=NEFR4Out%blockSize(ka,dimk,'alpha'); nB2=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA2=NEFR4Out%blockSize(la,diml,'alpha'); nB2=NEFR4Out%blockSize(la,diml,'beta')
            end if
            if(ia.eq.3) then
              nA3=NEFR4Out%blockSize(ia,dimi,'alpha'); nB3=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.3) then
              nA3=NEFR4Out%blockSize(ja,dimj,'alpha'); nB3=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.3) then
              nA3=NEFR4Out%blockSize(ka,dimk,'alpha'); nB3=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA3=NEFR4Out%blockSize(la,diml,'alpha'); nB3=NEFR4Out%blockSize(la,diml,'beta')
            end if
            if(ia.eq.4) then
              nA4=NEFR4Out%blockSize(ia,dimi,'alpha'); nB4=NEFR4Out%blockSize(ia,dimi,'beta')
            else if(ja.eq.4) then
              nA4=NEFR4Out%blockSize(ja,dimj,'alpha'); nB4=NEFR4Out%blockSize(ja,dimj,'beta')
            else if(ka.eq.4) then
              nA4=NEFR4Out%blockSize(ka,dimk,'alpha'); nB4=NEFR4Out%blockSize(ka,dimk,'beta')
            else
              nA4=NEFR4Out%blockSize(la,diml,'alpha'); nB4=NEFR4Out%blockSize(la,diml,'beta')
            end if
            call NEFR4Out%Main_R4(nef_bidx)%init(nA1,nB1,nAlpha2=nA2,nBeta2=nB2,&
                 nAlpha3=nA3,nBeta3=nB3,nAlpha4=nA4,nBeta4=nB4)
            NEFR4Out%Zero_Flag(nef_bidx) = .false.
            do i = 1, NEFR4Out%blockSize(ia,dimi,'alpha')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'alpha')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'all')),&
                       'alpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,-1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'all')),&
                       'alphabeta'))
               end do
              end do
             end do
            end do
            do i = 1, NEFR4Out%blockSize(ia,dimi,'beta')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'beta')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'all')),&
                       'betaalpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,-1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'b','b',dimk,diml,ka,la,'all')),&
                       'beta'))
               end do
              end do
             end do
            end do
            do i = 1, NEFR4Out%blockSize(ia,dimi,'alpha')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'beta')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,-1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'all')),&
                       'alpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'all')),&
                       'alphabeta'))
               end do
              end do
             end do
            end do
            do i = 1, NEFR4Out%blockSize(ia,dimi,'beta')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'alpha')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,-1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'all')),&
                       'betaalpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,1_int64],&
                     NEFMat%Main_Matrix(finalid)%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'b','a',dimk,diml,ka,la,'all')),&
                       'beta'))
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do
      case default
        call mqc_error_a('mqc_nef_matrix_dimention_expand: Unrecognized excitation type',6,&
          'spinType',mySpinType)
      end select
!
      end function mqc_nef_matrix_dimention_expand
!
      End Module MQC_NEF  
