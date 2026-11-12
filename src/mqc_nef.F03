
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
      Private::mqc_nef_r4_allocate_slots,mqc_nef_matrix_allocate_blocks
      Private::mqc_nef_symmetry_label,mqc_nef_scalar_symmetry
      Private::mqc_nef_maps_equal,mqc_nef_r4_resolve,mqc_nef_r4_full_storage
      Private::mqc_nef_set_full_grid,mqc_nef_set_rebuild_cell_matrix
      Private::mqc_nef_set_put_built_vector,mqc_nef_set_put_built_matrix,mqc_nef_set_put_built_r4
      Private::mqc_nef_pair_order_label,mqc_nef_matrix_equal,mqc_nef_plain_matrix_label
      Private::mqc_nef_r4_bref,mqc_nef_r4_bref_pairs
      Private::mqc_nef_r4_collapse_part,mqc_nef_r4_pair_spins,mqc_nef_r4_collapse_channel
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
        Character(Len=16),Private::Symmetry='asymmetric'   ! operator label: hermitian/antihermitian/asymmetric
        Logical,Private::sameBraKet=.false.   ! .true. iff the object is <I|O|I>
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
        Procedure, Public::iput => mqc_nef_vector_eigen_put
        Procedure, Public::iat => mqc_nef_vector_eigen_at
        Procedure, Public::iref => mqc_nef_vector_eigen_ref
        Procedure, Public::isBlockZero => mqc_nef_vector_is_block_zero
        Procedure, Public::isSameBraKet => mqc_nef_vector_is_same_braket
        Procedure, Public::getSymm => mqc_nef_vector_get_symmetry
      End Type
!
!        !!!! TODO:
!        !!! Usually a NEF_Matrix is not symmetrical, but their NEF_Set is. 
!        !!! An outerCollapsed Hermitian NEF_Set should be a Hermitian NEF_Matrix
!        !!!
!        !!! <New Space-type variable: orbital, single, double
!        !!!   collapse: orbital x orbital -> single
!        !!!             single  x single  -> double
!        !!!             single  x double  -> triple
!        !!!             single  x triple  -> triple
!        !!!             double  x double  -> quadruple
!        !!!
!        !!!             orbital x single  -> ?
!        !!!!!!
!         The symmetry of the operator (not the object!) is storaged at
!         - Symmetry = asymmetric, hermitian or antihermitian
!         In this context, an operator is only considered hermitian or antihermitian if
!         its complex-conjugated returns an operator with the same 'shape'.
!         Example:
!           dagger(<I|a_p^dagger a_q|J>) = <J|a_q^dagger a_p|I>        :  hermitian
!           dagger(<I|a_p^dagger a_q h |J>) != <J|a_q^dagger a_p h |I> :  asymetric
!         - Symmetries are only taken in consideration when Bra and Ket are the same (i.e. <I|O|I>)
!>    \brief <b> </b>
      Type MQC_NEF_Matrix
        Type(MQC_SCF_Integral),Dimension(:),Allocatable,Private::Main_Matrix
        Logical,Dimension(:),Allocatable,Private::Zero_Flag
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Map_Matrix_Col,Map_Matrix_Row
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_col,Size_vector_beta_col
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_row,Size_vector_beta_row
        Integer,Private::NumSpace_col=0,NumSpace_row=0
        Character(Len=16),Private::Symmetry='asymmetric'   ! operator label: hermitian/antihermitian/asymmetric, set at init only (I1)
        Logical,Private::sameBraKet=.false.   ! .true. iff the object is <I|O|I>, set at init only (I6)
        Logical,Private::Mirrored=.false.   ! sameBraKet .and. Symmetry/='asymmetric' .and. equal row/col maps: one block per dagger pair (I1)
      Contains
        Procedure, Public::print => mqc_print_nef_matrix
        Procedure, Public::init => mqc_nef_matrix_initialize
        Procedure, public::collapseIndex => mqc_nef_matrix_collapse_index
        Procedure, public::collapseMap => mqc_nef_matrix_collapse_map
        Procedure, public::blockIndex => mqc_nef_matrix_index
        Procedure, public::at => mqc_nef_matrix_at
        Procedure, Public::imp => mqc_nef_matrix_import_integral
        Procedure, Public::iput => mqc_nef_matrix_integral_put
        Procedure, Public::iat => mqc_nef_matrix_integral_at
        Procedure, Public::export => mqc_nef_matrix_export_integral
        Procedure, Public::slice => mqc_nef_matrix_dimention_slice
        Procedure, Public::collapse => mqc_nef_matrix_dimention_collapse
        Procedure, Public::expand => mqc_nef_matrix_dimention_expand
        Procedure, Public::dimen => mqc_nef_matrix_dimension
        Procedure, Public::getSymm => mqc_NEF_Matrix_Get_Symmetry
        Procedure, Public::isMapSymm => mqc_NEF_Matrix_Test_Map_Symmetry
!>    \brief <b> Outer Product of two variables into a two NEF R4</b>
        Procedure, Public::outer2R4 => mqc_NEFMatrix_NEFMatrix_outer_product
!>    \brief <b> Outer Product of two variables into a NEF R4</b>
        Procedure, Public::crossouter2R4 => mqc_NEFMatrix_NEFMatrix_crossouter_product
        Procedure, public::getBlockMat => mqc_nef_matrix_output_block_matrix
        Procedure, public::getMap => mqc_nef_matrix_get_map
        Procedure, Public::pinv => mqc_nef_matrix_block_diagonal_pinv
        Procedure, Public::getRotation => mqc_nef_matrix_get_generator_rotation
        Procedure, Public::getBlockRotation => mqc_nef_matrix_get_block_rotation
        Procedure, Public::exp => mqc_nef_matrix_exp
        Procedure, Public::isBlockZero => mqc_nef_matrix_is_block_zero
        Procedure, Public::isSameBraKet => mqc_nef_matrix_is_same_braket
      End Type
!
!>    \brief <b> MQC_NEF_R4_Slot holds one NEF_R4 space block (P,Q,R,S). It is a way to
!>    reduce memory usage of zero blocks: an unallocated R4Slot costs 8 B instead of a
!>    full MQC_TwoERIs (10 240 B). Accessed only through the NEF_R4 iput/iat/iref/
!>    isBlockZero/put/at procedures; iput and put allocate it. </b>
      Type,Private::MQC_NEF_R4_Slot
        Type(MQC_TwoERIs),Allocatable::R4Slot
      End Type
!
!>    \brief <b> </b>
      Type MQC_NEF_R4
        Type(MQC_NEF_R4_Slot),Dimension(:),Allocatable,Private::Main_R4
        Logical,Dimension(:),Allocatable,Private::Zero_Flag      ! kept during the I0 transition
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Map_Matrix_1,Map_Matrix_2,&
                                                                    Map_Matrix_3,Map_Matrix_4
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_1,Size_vector_beta_1
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_2,Size_vector_beta_2
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_3,Size_vector_beta_3
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_4,Size_vector_beta_4
        Integer,Private::NumSpace_1=0,NumSpace_2=0,NumSpace_3=0,NumSpace_4=0
        Character(Len=16),Private::Symmetry='asymmetric'   ! operator label: hermitian/antihermitian/asymmetric, set at init only (I1)
        Logical,Private::sameBraKet=.false.   ! .true. iff the object is <I|O|I>, set at init only (I6)
        Logical,Private::Mirrored=.false.   ! sameBraKet .and. Symmetry/='asymmetric' .and. Map_1=Map_4 .and. Map_2=Map_3: one block per dagger pair (I2)
      Contains
        Procedure, Public::print => mqc_print_nef_r4
        Procedure, Public::init => mqc_nef_r4_initialize
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
        Procedure, Public::transform => mqc_nef_r4_full_transform
        Procedure, Private::mqc_nef_r4_eris_put,mqc_nef_r4_eris_put_pairs
        Generic, Public::iput => mqc_nef_r4_eris_put,mqc_nef_r4_eris_put_pairs
        Procedure, Private::mqc_nef_r4_eris_at,mqc_nef_r4_eris_at_pairs
        Generic, Public::iat => mqc_nef_r4_eris_at,mqc_nef_r4_eris_at_pairs
        Procedure, Private::mqc_nef_r4_eris_ref,mqc_nef_r4_eris_ref_pairs
        Generic, Public::iref => mqc_nef_r4_eris_ref,mqc_nef_r4_eris_ref_pairs
        Procedure, Private::mqc_nef_r4_is_block_zero,mqc_nef_r4_is_block_zero_pairs
        Generic, Public::isBlockZero => mqc_nef_r4_is_block_zero,mqc_nef_r4_is_block_zero_pairs
        Procedure, Public::isSameBraKet => mqc_nef_r4_is_same_braket
        Procedure, Public::getSymm => mqc_nef_r4_get_symmetry
        Procedure, Public::isMirrorSlot => mqc_nef_r4_is_mirror_slot
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
        Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,Private::Map_Matrix_1,Map_Matrix_2,&
                                                                    Map_Matrix_3,Map_Matrix_4
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_1,Size_vector_beta_1
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_2,Size_vector_beta_2
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_3,Size_vector_beta_3
        Integer,Dimension(:),Allocatable,Private::Size_vector_alpha_4,Size_vector_beta_4
        Integer,Dimension(4)::NumSpace = -1
        Logical,Private::sameBraKetSpace=.false.   ! .true. iff the (I,J) grid is bra x ket over one determinant list; only when NumI = NumJ (I6)
        Character(Len=16),Private::Symmetry='asymmetric'   ! operator label of every cell (hermitian/antihermitian/asymmetric), set at init only (P6)
        Logical,Private::GridMirrored=.false.   ! sameBraKetSpace .and. Symmetry/='asymmetric' .and. paired cell maps: only I>=J cells stored (P6)
      Contains
        Procedure, Public::print => mqc_print_nef_set
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
        Procedure, Public::mergeSpace => mqc_nef_set_outer_merge_space
        Procedure, Public::separateSpace => mqc_nef_set_outer_unmerge_space
        Procedure, Public::fullCollapse => mqc_nef_set_full_collapse
        Procedure, Public::fullExpand => mqc_nef_set_full_expand
        Procedure, Public::getMap => mqc_nef_set_get_map
        Procedure, Public::pinvR4 => mqc_nef_set_block_diagonal_pinv
        Procedure, Public::isSameBraKetSpace => mqc_nef_set_is_same_braket_space
        Procedure, Public::getSymm => mqc_nef_set_get_symmetry
        Procedure, Public::isMirrorCell => mqc_nef_set_is_mirror_cell
      End Type
!
!----------------------------------------------------------------
!                                                               |
!     PROCEDURE INTERFACES                                      |
!                                                               |
!----------------------------------------------------------------
!
!
!>    \brief <b> Prints an object</b>
      interface mqc_print
        module procedure mqc_print_nef_vector
        module procedure mqc_print_nef_matrix
        module procedure mqc_print_nef_r4
        module procedure mqc_print_nef_set
      end interface
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
        Module Procedure MQC_NEF_Vector_Conjugate
        Module Procedure MQC_NEFSet_Conjugate
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
!>    \brief <b> Hermitian sandwich B^dagger A B (extends the HermMatMul generic of
!>    mqc_algebra / mqc_est) </b>
      Interface HermMatMul
        Module Procedure MQC_Matrix_NEFSet_Herm_Sandwich
        Module Procedure MQC_NEFSet_Matrix_Herm_Sandwich
        Module Procedure MQC_NEFMatrix_Herm_Sandwich
        Module Procedure MQC_Integral_NEFMatrix_Herm_Sandwich
      End Interface
!
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
!     PROCEDURE MQC_NEF_Vector_Initialize
!
!>    \brief <b> MQC_NEF_Vector_Initialize is a subroutine that initializes a
!>    NEF Vector object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Initialize initializes a MQC_NEF_Vector over the spaces
!>    of initial_map_matrix. Default element values are 0.0 (all blocks zero),
!>    or every element is set to the optional scalar. The operator label and
!>    the sameBraKet flag are set here only; on a single NEF_Vector they have
!>    no in-object effect.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFVec
!>    \verbatim
!>        NEFVec is Class(MQC_NEF_Vector)
!>        The MQC NEF Vector to initialize.
!>    \endverbatim
!>
!>    \param[in] initial_map_matrix
!>    \verbatim
!>        initial_map_matrix is Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable
!>        Maps each element from an original basis into the nef basis.
!>    \endverbatim
!>
!>    \param[in] scalar
!>    \verbatim
!>        scalar is Class(*),Optional
!>        Value of every element. If not present, the value is 0.0.
!>    \endverbatim
!>
!>    \param[in] symmetry
!>    \verbatim
!>        symmetry is Character(len=*), Optional
!>        Operator label given by the builder: 'hermitian', 'antihermitian'
!>        or 'asymmetric' (default). Never changed after init.
!>    \endverbatim
!>
!>    \param[in] sameBraKet
!>    \verbatim
!>        sameBraKet is Logical, Optional
!>        .true. iff the object is <I|O|I>. Default .false.. Never changed
!>        after init.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_vector_initialize(NEFVec,initial_map_matrix,scalar,symmetry,sameBraKet)
!
        implicit none
        class(mqc_nef_vector),intent(inOut)::NEFVec
        type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in)::initial_map_matrix
        Class(*),Optional::scalar
        character(len=*),intent(in),optional::symmetry
        logical,intent(in),optional::sameBraKet
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
        NEFVec%Symmetry = mqc_nef_symmetry_label(symmetry,'mqc_nef_vector_initialize')
        NEFVec%sameBraKet = .false.
        if(present(sameBraKet)) NEFVec%sameBraKet = sameBraKet
!
      end subroutine mqc_nef_vector_initialize
!
!
!     PROCEDURE MQC_NEF_Symmetry_Label
!
!>    \brief <b> MQC_NEF_Symmetry_Label (private) validates the operator symmetry
!>    label given to a NEF object init. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Symmetry_Label returns the operator label of a NEF object:
!>    'hermitian', 'antihermitian' or 'asymmetric' (lower case). The label
!>    is a property of the operator O and is given by the builder at init;
!>    NEF never derives it. An absent label gives 'asymmetric'; any other
!>    value raises an error.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] symmetry
!>    \verbatim
!>        symmetry is Character(len=*), Optional
!>        The requested operator label.
!>    \endverbatim
!>
!>    \param[in] caller
!>    \verbatim
!>        caller is Character(len=*)
!>        Name of the calling init, used in the error message.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_symmetry_label(symmetry,caller) result(label)
!
      implicit none
      character(len=*),intent(in),optional::symmetry
      character(len=*),intent(in)::caller
      character(len=16)::label
!
      label = 'asymmetric'
      if(present(symmetry)) call string_change_case(symmetry,'l',label)
      select case(trim(label))
      case('hermitian','antihermitian','asymmetric')
        continue
      case default
        call mqc_error_a(trim(caller)//': invalid NEF symmetry label (hermitian, '//&
          'antihermitian or asymmetric)',6,'symmetry',label)
      end select
!
      end function mqc_nef_symmetry_label
!
!
!     PROCEDURE MQC_NEF_Scalar_Symmetry
!
!>    \brief <b> MQC_NEF_Scalar_Symmetry (private) returns the operator label of
!>    c*O from the label of O. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Scalar_Symmetry propagates the operator label through a
!>    scalar multiply: a real c keeps it, a purely imaginary c swaps
!>    'hermitian' and 'antihermitian', any other c gives 'asymmetric'.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] label
!>    \verbatim
!>        label is Character(len=*)
!>        The operator label of O.
!>    \endverbatim
!>
!>    \param[in] scalar
!>    \verbatim
!>        scalar is Type(MQC_Scalar)
!>        The scalar c.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_scalar_symmetry(label,scalar) result(outSymm)
!
      implicit none
      character(len=*),intent(in)::label
      type(mqc_scalar),intent(in)::scalar
      character(len=16)::outSymm
      Real(kind=real64)::Thresh
      Complex(kind=real64)::c
!
      If(MQC_very_small.gt.0.0) then
        Thresh = MQC_very_small
      Else
        Thresh = 1e-10
      EndIf
      c = scalar%cval()
      outSymm = 'asymmetric'
      if(abs(aimag(c)).lt.Thresh) then
        outSymm = label
      else if(abs(real(c)).lt.Thresh) then
        select case(trim(label))
        case('hermitian')
          outSymm = 'antihermitian'
        case('antihermitian')
          outSymm = 'hermitian'
        end select
      end if
!
      end function mqc_nef_scalar_symmetry
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
!
      write(6,*) 'Symmetry = ',NEFMat%getSymm()
      if(allocated(NEFMat%Zero_Flag)) then
        if(all(NEFMat%Zero_Flag)) then
          write(iout,1000) '(all zero) dims: row alpha='//trim(axisDimsStr(NEFMat%Size_vector_alpha_row))//&
               ' beta='//trim(axisDimsStr(NEFMat%Size_vector_beta_row))//&
               ' col alpha='//trim(axisDimsStr(NEFMat%Size_vector_alpha_col))//&
               ' beta='//trim(axisDimsStr(NEFMat%Size_vector_beta_col))
          return
        end if
      end if
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
              if (NEFMat%isBlockZero(j,i)) then
                do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                  Line = trim(Line)//'     '//num2char(0.00000)
                  if (sub_i.ne.NEFMat%Size_vector_alpha_col(i)) Line = trim(Line)//' '
                end do
              else
              tmp_int = NEFMat%iat(j,i)
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
              if (NEFMat%isBlockZero(j,i)) then
                do sub_i = 1, NEFMat%Size_vector_beta_col(i)
                  Line = trim(Line)//'     '//num2char(0.00000)
                  if (sub_i.ne.NEFMat%Size_vector_beta_col(i)) Line = trim(Line)//' '
                end do
              else
              tmp_int = NEFMat%iat(j,i)
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
              if (NEFMat%isBlockZero(j,i)) then
                do sub_i = 1, NEFMat%Size_vector_alpha_col(i)
                  Line = trim(Line)//'     '//num2char(0.00000)
                  if (sub_i.ne.NEFMat%Size_vector_beta_col(i)) Line = trim(Line)//' '
                end do
              else
              tmp_int = NEFMat%iat(j,i)
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
              if (NEFMat%isBlockZero(j,i)) then
                do sub_i = 1, NEFMat%Size_vector_beta_col(i)
                  Line = trim(Line)//'     '//num2char(0.00000)
                  if (sub_i.ne.NEFMat%Size_vector_alpha_col(i)) Line = trim(Line)//' '
                end do
              else
              tmp_int = NEFMat%iat(j,i)
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
            if (NEFMat%isBlockZero(j,i)) then
              write(iOut,1000) 'Sub-Space '//trim(num2char(j))//','//trim(num2char(i))//&
                   ' (all zero) dims: row alpha='//trim(num2char(NEFMat%Size_vector_alpha_row(j)))//&
                   ' beta='//trim(num2char(NEFMat%Size_vector_beta_row(j)))//&
                   ' col alpha='//trim(num2char(NEFMat%Size_vector_alpha_col(i)))//&
                   ' beta='//trim(num2char(NEFMat%Size_vector_beta_col(i)))
            else
              call mqc_print(NEFMat%iat(j,i),iOut,&
                    'Sub-Space '//trim(num2char(j))//','//trim(num2char(i)))
            end if
          end do
        end do
      end if
!
      contains
!
        function axisDimsStr(sizes) result(str)
        implicit none
        integer,dimension(:),intent(in)::sizes
        character(len=256)::str
        integer::k
!
        str = '['
        do k = 1, size(sizes)
          if(k.gt.1) str = trim(str)//','
          str = trim(str)//trim(num2char(sizes(k)))
        end do
        str = trim(str)//']'
!
        end function axisDimsStr
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
!>    scalar. The operator label and the sameBraKet flag are set here only. When
!>    sameBraKet is .true., the label is not 'asymmetric' and the row and column
!>    maps are equal, the object symmetry holds: one block per dagger pair is
!>    stored and import/iput enforce the symmetry (a scalar fill must obey it).
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
!>
!>    \param[in] symmetry
!>    \verbatim
!>        symmetry is Character(len=*), Optional
!>        Operator label given by the builder: 'hermitian', 'antihermitian'
!>        or 'asymmetric' (default). Never changed after init.
!>    \endverbatim
!>
!>    \param[in] sameBraKet
!>    \verbatim
!>        sameBraKet is Logical, Optional
!>        .true. iff the object is <I|O|I>. Default .false.. Never changed
!>        after init.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_matrix_initialize(NEFMat,initial_map_matrix_row,scalar,initial_map_matrix_col,&
           symmetry,sameBraKet)
!
        implicit none
        integer(kind=int64) :: nef_bidx
        class(mqc_nef_matrix),intent(inOut)::NEFMat
        type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in)::initial_map_matrix_row
        type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in),optional::initial_map_matrix_col
        type(mqc_scf_eigenvalues),dimension(:),allocatable::active_map_matrix_col
        Class(*),Optional::scalar
        character(len=*),intent(in),optional::symmetry
        logical,intent(in),optional::sameBraKet
        logical::same_vec
        integer::i,j
        type(mqc_scalar)::scalarCheck
        Real(kind=real64)::Thresh
        logical::scalarIsZero,scalarIsReal,scalarIsImag
!
        If(MQC_very_small.gt.0.0) then
          Thresh = MQC_very_small
        Else
          Thresh = 1e-10
        EndIf
!
        scalarIsReal = .false.
        scalarIsImag = .false.
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
          scalarIsZero = (abs(scalarCheck).le.Thresh)
          if (.not.scalarIsZero) then
            scalarIsReal = (abs(aimag(scalarCheck)).le.Thresh)
            scalarIsImag = (abs(real(scalarCheck)).le.Thresh)
          end if
        else
          scalarCheck = 0.0
          scalarIsZero = .true.
        end if
!
        NEFMat%NumSpace_row = size(initial_map_matrix_row)
        NEFMat%Map_Matrix_row = initial_map_matrix_row
        if(present(initial_map_matrix_col)) then
          NEFMat%NumSpace_col = size(initial_map_matrix_col)
          active_map_matrix_col = initial_map_matrix_col
          NEFMat%Map_Matrix_col = initial_map_matrix_col
          same_vec = mqc_NEF_Matrix_Test_Map_Symmetry(NEFMat)
        else
          NEFMat%NumSpace_col = NEFMat%NumSpace_row
          active_map_matrix_col = initial_map_matrix_row
          NEFMat%Map_Matrix_col = initial_map_matrix_row
          same_vec = .true.
        end if
!
!       Operator label and bra-ket flag are given by the builder (never derived).
!       One block per dagger pair is stored only when the object symmetry holds.
        NEFMat%Symmetry = mqc_nef_symmetry_label(symmetry,'mqc_nef_matrix_initialize')
        NEFMat%sameBraKet = .false.
        if(present(sameBraKet)) NEFMat%sameBraKet = sameBraKet
        NEFMat%Mirrored = (NEFMat%sameBraKet.and.NEFMat%Symmetry.ne.'asymmetric'.and.same_vec)
!       Every element equals the scalar c: the filled object obeys the label only
!       for a real c (hermitian) or an imaginary c (antihermitian).
        if(NEFMat%Mirrored.and..not.scalarIsZero) then
          if((NEFMat%Symmetry.eq.'hermitian'.and..not.scalarIsReal).or.&
             (NEFMat%Symmetry.eq.'antihermitian'.and..not.scalarIsImag)) &
            call mqc_error_a('mqc_nef_matrix_initialize: scalar does not obey the symmetry',&
              6,'symmetry',NEFMat%Symmetry)
        end if
!
!       !!!Bookkeeping!!!
        if(allocated(NEFMat%Size_vector_alpha_row)) deallocate(NEFMat%Size_vector_alpha_row)
        if(allocated(NEFMat%Size_vector_alpha_col)) deallocate(NEFMat%Size_vector_alpha_col)
        if(allocated(NEFMat%Size_vector_beta_row)) deallocate(NEFMat%Size_vector_beta_row)
        if(allocated(NEFMat%Size_vector_beta_col)) deallocate(NEFMat%Size_vector_beta_col)
        allocate(NEFMat%Size_vector_alpha_row(NEFMat%NumSpace_row))
        allocate(NEFMat%Size_vector_beta_row(NEFMat%NumSpace_row))
        allocate(NEFMat%Size_vector_alpha_col(NEFMat%NumSpace_col))
        allocate(NEFMat%Size_vector_beta_col(NEFMat%NumSpace_col))
!       !!!MainData!!!
        if(allocated(NEFMat%Main_matrix)) deallocate(NEFMat%Main_matrix)
        if(allocated(NEFMat%Zero_Flag)) deallocate(NEFMat%Zero_Flag)
        allocate(NEFMat%Main_matrix(NEFMat%NumSpace_row*NEFMat%NumSpace_col))
        allocate(NEFMat%Zero_Flag(NEFMat%NumSpace_row*NEFMat%NumSpace_col))
!
        NEFMat%Zero_Flag=.true.
        do i = 1, NEFMat%NumSpace_row
          NEFMat%Size_vector_alpha_row(i) = Size(initial_map_matrix_row(i),'alpha')
          NEFMat%Size_vector_beta_row(i) = Size(initial_map_matrix_row(i),'beta')
          do j = 1, NEFMat%NumSpace_col
            if(i.eq.1) then
              NEFMat%Size_vector_alpha_col(j) = Size(active_map_matrix_col(j),'alpha')
              NEFMat%Size_vector_beta_col(j) = Size(active_map_matrix_col(j),'beta')
            end if
            if (i.lt.j.and.NEFMat%Mirrored) cycle
            if (same_vec.and.i.eq.j) then
              nef_bidx = NEFMat%blockIndex([1,i],[2,i])
              NEFMat%Zero_Flag(nef_bidx) = scalarIsZero
              !if (.not.scalarIsZero) then
                call NEFMat%Main_matrix(nef_bidx)%init(&
                       Size(initial_map_matrix_row(i),'alpha'),&
                       Size(initial_map_matrix_row(i),'beta'),&
                       scalar=scalarCheck)
              !end if
            else
              nef_bidx = NEFMat%blockIndex([1,i],[2,j])
              NEFMat%Zero_Flag(nef_bidx) = scalarIsZero
              !if (.not.scalarIsZero) then
                call NEFMat%Main_matrix(nef_bidx)%init(&
                       Size(initial_map_matrix_row(i),'alpha'),&
                       Size(initial_map_matrix_row(i),'beta'),&
                       scalar=scalarCheck,&
                       nAlpha2=Size(active_map_matrix_col(j),'alpha'),&
                       nBeta2=Size(active_map_matrix_col(j),'beta'))
              !end if
            end if
          end do
        end do
!
      end subroutine mqc_nef_matrix_initialize
!
!
!     PROCEDURE MQC_NEF_Matrix_Test_Map_Symmetry
!
!>    \brief <b> MQC_NEF_Matrix_Test_Map_Symmetry is a function that check if 
!>      the NEF_Matrix row and column maps are the same</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to be queried.
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!
      function mqc_NEF_Matrix_Test_Map_Symmetry(NEFMat) result(same)
!
        implicit none
        class(mqc_nef_matrix),intent(in)::NEFMat
        logical::same
        integer::i
!
        same = (NEFMat%NumSpace_col.eq.NEFMat%NumSpace_row)
        if(same) then
          do i = 1, NEFMat%NumSpace_col 
            if(NEFMat%Map_Matrix_Col(i).ne.NEFMat%Map_Matrix_Row(i)) then
              same = .false.
              exit
            end if
          end do
        end if
!
      end function mqc_NEF_Matrix_Test_Map_Symmetry    
!
!     PROCEDURE MQC_NEF_Matrix_Get_Symmetry
!
!>    \brief <b> MQC_NEF_Matrix_Get_Symmetry is a function that retrun a 
!>     NEF_Matrix symmetry</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>       
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to be queried.
!>    \endverbatim
!>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_NEF_Matrix_Get_Symmetry(NEFMat) result(symm)
!
        implicit none
        class(mqc_nef_matrix),intent(in)::NEFMat
        Character(Len=16)::symm
!
        symm = NEFMat%Symmetry
!
      end function mqc_NEF_Matrix_Get_Symmetry
!
!
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
        type(mqc_vector)::row_alpha,row_beta,col_alpha,col_beta
        integer::row,col,rowSize,colSize
!
!       The object symmetry is fixed at init; with one block per dagger pair the
!       imported integral must obey it (the stored blocks are the lower ones).
        if(NEFMat%Mirrored) then
          if(.not.integral%testSymm(NEFMat%Symmetry,math=.true.)) &
            call mqc_error_a('mqc_nef_matrix_import_integral: integral does not obey '//&
              'the NEF_Matrix symmetry',6,'NEFMat%Symmetry',NEFMat%Symmetry,&
              'integral%getSymm()',integral%getSymm())
        end if
!
        do row = 1, NEFMat%NumSpace_row
          row_alpha = NEFMat%Map_Matrix_Row(row)%getBlock('alpha')
          row_beta = NEFMat%Map_Matrix_Row(row)%getBlock('beta')
          rowSize = size(row_alpha)+size(row_beta)
          do col = 1, NEFMat%NumSpace_col
            col_alpha = NEFMat%Map_Matrix_Col(col)%getBlock('alpha')
            col_beta = NEFMat%Map_Matrix_Col(col)%getBlock('beta')
            colSize = size(col_alpha)+size(col_beta)
            if(row.lt.col.and.NEFMat%Mirrored) cycle
            nef_bidx = NEFMat%blockIndex([1,row],[2,col])
            if(size(row_beta).eq.0.and.size(col_beta).gt.0) then
              tmp_integral = &
                 integral%orbitals('useStrings',&
                  alphaOrbsIn = col_alpha%getArrayI(),&
                   BetaOrbsIn = col_beta%getArrayI(),axis=2)
              NEFMat%Main_matrix(nef_bidx)=&
                  tmp_integral%orbitals('useStrings',&
                alphaOrbsIn = row_alpha%getArrayI(),&
                 BetaOrbsIn = row_beta%getArrayI(),axis=1)
            else
              tmp_integral = &
                 integral%orbitals('useStrings',&
                  alphaOrbsIn = row_alpha%getArrayI(),&
                   BetaOrbsIn = row_beta%getArrayI(),axis=1)
              NEFMat%Main_matrix(nef_bidx)=&
                  tmp_integral%orbitals('useStrings',&
                alphaOrbsIn = col_alpha%getArrayI(),&
                 BetaOrbsIn = col_beta%getArrayI(),axis=2)
            end if
            NEFMat%Zero_Flag(nef_bidx) = (rowSize.eq.0.or.colSize.eq.0)
          end do
        end do
!
      end subroutine mqc_nef_matrix_import_integral
!
!
!     PROCEDURE MQC_NEF_Matrix_Integral_Put
!
!>    \brief <b> MQC_NEF_Matrix_Integral_Put adds an MQC SCF Integral into
!>    one of the Main_Matrix entries of an initialized NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Integral_Put stores an MQC SCF Integral as the
!>    (row,col) space-pair block of an already-initialized NEF Matrix
!>    object. The integral's alpha/beta row and column dimensions must
!>    match the target block's dimensions exactly. When the object
!>    symmetry holds (sameBraKet, symmetry not 'asymmetric' and equal
!>    row/column maps, fixed at init) the put enforces it: a diagonal
!>    block must obey the symmetry (error otherwise) and a block above
!>    the diagonal is stored as its (+/-) dagger in the shared slot.
!>    Otherwise the block is stored as given. A zero off-diagonal
!>    integral is stored and flagged zero.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to update.
!>    \endverbatim
!>
!>    \param[in] integral
!>    \verbatim
!>        integral is type(mqc_scf_integral)
!>        The value to add into the target block.
!>    \endverbatim
!>
!>    \param[in] row
!>    \verbatim
!>        row is integer
!>        The row space index of the target block.
!>    \endverbatim
!>
!>    \param[in] col
!>    \verbatim
!>        col is integer
!>        The column space index of the target block.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_matrix_integral_put(NEFMat,integral,row,col)
!
      implicit none
      class(mqc_nef_matrix),intent(inOut)::NEFMat
      type(mqc_scf_integral),intent(in)::integral
      integer,intent(in)::row,col
      integer(kind=int64)::nef_bidx
      logical::zeroBlock
!
      if(integral%blockSize('alpha',1).ne.NEFMat%Size_vector_alpha_row(row).or.&
         integral%blockSize('beta',1).ne.NEFMat%Size_vector_beta_row(row).or.&
         integral%blockSize('alpha',2).ne.NEFMat%Size_vector_alpha_col(col).or.&
         integral%blockSize('beta',2).ne.NEFMat%Size_vector_beta_col(col)) &
        call mqc_error_i('mqc_nef_matrix_integral_put: integral dimensions do not '//&
          'match the target block dimensions at (row,col)',6,&
          'row',row,'col',col,&
          'integral_alphaRow',integral%blockSize('alpha',1),&
          'block_alphaRow',NEFMat%Size_vector_alpha_row(row),&
          'integral_betaRow',integral%blockSize('beta',1),&
          'block_betaRow',NEFMat%Size_vector_beta_row(row))
!
!     A zero off-diagonal block is flagged zero; a diagonal block is always stored.
      zeroBlock = (row.ne.col.and.integral%getSymm().eq.'fullsymmetric')
      nef_bidx = NEFMat%blockIndex([1,row],[2,col])
      if(.not.NEFMat%Mirrored.or.row.gt.col) then
        NEFMat%Main_Matrix(nef_bidx) = integral
      else if(row.eq.col) then
!       Object symmetry: a diagonal block carries the operator symmetry in-block
!       (mathematical test: math=.true. also checks the diagonal elements).
        if(.not.integral%testSymm(NEFMat%Symmetry,math=.true.)) &
          call mqc_error_a('mqc_nef_matrix_integral_put: diagonal block does not obey '//&
            'the NEF_Matrix symmetry',6,'NEFMat%Symmetry',NEFMat%Symmetry,&
            'integral%getSymm()',integral%getSymm())
        NEFMat%Main_Matrix(nef_bidx) = integral
      else if(NEFMat%Symmetry.eq.'hermitian') then
!       Object symmetry, row < col: the stored (col,row) slot gets the mirror.
        NEFMat%Main_Matrix(nef_bidx) = dagger(integral)
      else
        NEFMat%Main_Matrix(nef_bidx) = (-1)*dagger(integral)
      end if
      NEFMat%Zero_Flag(nef_bidx) = zeroBlock
!
      end subroutine mqc_nef_matrix_integral_put
!
!
!     PROCEDURE MQC_NEF_Matrix_Integral_At
!
!>    \brief <b> MQC_NEF_Matrix_Integral_At returns one of the Main_Matrix
!>    entries of an initialized NEF Matrix object as an MQC SCF Integral</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Integral_At returns the (row,col) space-pair block of
!>    an already-initialized NEF Matrix object as an MQC SCF Integral. This
!>    is the read-side counterpart to MQC_NEF_Matrix_Integral_Put. When the
!>    block is flagged zero, a properly-dimensioned zero-valued integral is
!>    built from Map_Matrix_Row/Map_Matrix_Col instead of returning the
!>    (never initialized) Main_Matrix entry.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to read from.
!>    \endverbatim
!>
!>    \param[in] row
!>    \verbatim
!>        row is integer
!>        The row space index of the target block.
!>    \endverbatim
!>
!>    \param[in] col
!>    \verbatim
!>        col is integer
!>        The column space index of the target block.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_integral_at(NEFMat,row,col) result(integralOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      integer,intent(in)::row,col
      type(mqc_scf_integral)::integralOut
      integer(kind=int64)::nef_bidx
!
      nef_bidx = NEFMat%blockIndex([1,row],[2,col])
      if(.not.allocated(NEFMat%Zero_Flag)) call mqc_error(&
           'mqc_nef_matrix_integral_at: NEF_Matrix is open (never initialized)',6)
      if(NEFMat%Zero_Flag(nef_bidx)) then
!       The returned block is always the (row,col) block, also above the block
!       diagonal (row<col): rows from space row, columns from space col.
        call integralOut%init(&
               Size(NEFMat%Map_Matrix_Row(row),'alpha'),&
               Size(NEFMat%Map_Matrix_Row(row),'beta'),&
               nAlpha2=Size(NEFMat%Map_Matrix_Col(col),'alpha'),&
               nBeta2=Size(NEFMat%Map_Matrix_Col(col),'beta'))
      else if(.not.NEFMat%Mirrored.or.row.ge.col) then
        integralOut = NEFMat%Main_Matrix(nef_bidx)
      else if(NEFMat%Symmetry.eq.'hermitian') then
        integralOut = dagger(NEFMat%Main_Matrix(nef_bidx))
      else
        integralOut = (-1)*dagger(NEFMat%Main_Matrix(nef_bidx))
      end if
!
      end function mqc_nef_matrix_integral_at
!
!
!     PROCEDURE mqc_nef_matrix_is_block_zero
!
!>    \brief <b> mqc_nef_matrix_is_block_zero returns the Zero_Flag
!>    entry of an initialized NEF Matrix object</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    mqc_nef_matrix_is_block_zero returns the Zero_Flag
!>    entry of an initialized NEF Matrix object
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The MQC NEF Matrix to read from.
!>    \endverbatim
!>
!>    \param[in] row
!>    \verbatim
!>        row is integer
!>        The row space index of the target block.
!>    \endverbatim
!>
!>    \param[in] col
!>    \verbatim
!>        col is integer
!>        The column space index of the target block.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_is_block_zero(NEFMat,row,col) result(iszero)

!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      integer,intent(in)::row,col
      logical::isZero
      integer(kind=int64)::nef_bidx
!
      isZero = NEFMat%Zero_Flag(NEFMat%blockIndex([1,row],[2,col]))
!
      end function mqc_nef_matrix_is_block_zero
!
!
!     PROCEDURE MQC_NEF_Matrix_Export_Integral
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
      !!!! TODO: make the inetgral with symmetry based on the NEF Symm!!
!
        implicit none
        class(mqc_nef_matrix),intent(in)::NEFMat
        type(mqc_scf_integral)::integral
        type(mqc_scf_integral)::tmp_integral
        type(mqc_vector)::tmp_alpha,tmp_beta
        integer::ia,ja,ib,jb,row,col,maxColAlpha,maxColBeta,maxRowAlpha,maxRowBeta
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
            if(NEFMat%isBlockZero(row,col)) cycle
            tmp_integral = NEFMat%iat(row,col)
            do ia = 1, NEFMat%Size_vector_alpha_row(row)
              do ja = 1, NEFMat%Size_vector_alpha_col(col)
                   call integral%put(tmp_integral%at(ia,ja,'alpha'),&
                              NEFMat%Map_Matrix_Row(row)%at(ia,'alpha'),&
                              NEFMat%Map_Matrix_Col(col)%at(ja,'alpha'),'alpha')
              end do
              do ja = 1, NEFMat%Size_vector_beta_col(col)
                call integral%put(tmp_integral%at(ia,ja,'betaalpha'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'alpha'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'beta'),'betaalpha')
              end do
            end do
            do ia = 1, NEFMat%Size_vector_beta_row(row)
              do ja = 1, NEFMat%Size_vector_alpha_col(col)
                call integral%put(tmp_integral%at(ia,ja,'alphabeta'),&
                           NEFMat%Map_Matrix_Row(row)%at(ia,'beta'),&
                           NEFMat%Map_Matrix_Col(col)%at(ja,'alpha'),'alphabeta')
              end do
              do ja = 1, NEFMat%Size_vector_beta_col(col)
                call integral%put(tmp_integral%at(ia,ja,'beta'),&
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
      integer(kind=int64)::rowA,colB,colA
      type(mqc_scf_integral)::blockOut
      logical::started
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
          initial_map_matrix_col=NEFMatrixB%Map_Matrix_Col,&
           sameBraKet=(NEFMatrixA%sameBraKet.and.NEFMatrixB%sameBraKet))
      do rowA = 1, NEFMatrixA%NumSpace_row
        do colB = 1, NEFMatrixB%NumSpace_col
          started = .false.
          do colA = 1, NEFMatrixA%NumSpace_col
            if (NEFMatrixA%isBlockZero(rowA,colA).or.NEFMatrixB%isBlockZero(colA,colB)) cycle
            if (.not.started) then
              blockOut = matmul(NEFMatrixA%iat(rowA,colA),NEFMatrixB%iat(colA,colB))
              started = .true.
            else
              blockOut = blockOut + matmul(NEFMatrixA%iat(rowA,colA),NEFMatrixB%iat(colA,colB))
            end if
          end do
          if (started) call NEFMatrixOut%iput(blockOut,rowA,colB)
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
      type(mqc_nef_vector),intent(in),target::NEFVectorB
      type(mqc_nef_vector)::NEFVectorOut
      integer(kind=int64)::rowA,colA
      type(mqc_scf_integral)::blockA
      type(mqc_scf_eigenvalues)::acc
      type(mqc_vector)::padAlpha,padBeta
      integer(kind=int64)::nRowAlpha,nRowBeta
      logical::started
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

      call NEFVectorOut%init(NEFMatrixA%Map_Matrix_Row,&
           sameBraKet=(NEFMatrixA%sameBraKet.and.NEFVectorB%sameBraKet))
      do rowA = 1, NEFMatrixA%NumSpace_row
        started = .false.
        do colA = 1, NEFMatrixA%NumSpace_col
          if (NEFMatrixA%isBlockZero(rowA,colA).or.NEFVectorB%isBlockZero(colA)) cycle
          blockA = NEFMatrixA%iat(rowA,colA)
          if (.not.started) then
            acc = (blockA.dot.NEFVectorB%iref(colA))
            started = .true.
          else
            acc = acc + (blockA.dot.NEFVectorB%iref(colA))
          end if
        end do
        if (started) then
!         a spin part missing from every term is zero: give it the row-space size
          nRowAlpha = NEFMatrixA%Size_vector_alpha_row(rowA)
          nRowBeta = NEFMatrixA%Size_vector_beta_row(rowA)
          if (size(acc,'alpha').ne.nRowAlpha.or.size(acc,'beta').ne.nRowBeta) then
            if (size(acc,'alpha').eq.nRowAlpha) then
              padAlpha = acc%getBlock('alpha')
            else
              call padAlpha%init(nRowAlpha)
            end if
            if (size(acc,'beta').eq.nRowBeta) then
              padBeta = acc%getBlock('beta')
            else
              call padBeta%init(nRowBeta)
            end if
            call mqc_eigenvalues_allocate(acc,'','spin',padAlpha,padBeta)
          end if
          call NEFVectorOut%iput(acc,rowA)
        end if
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
        if (.not.NEFMatA%isBlockZero(row,1)) &
          call NEFMatOut%iput(matmul(NEFMatA%iat(row,1),IntegralB),row,1)
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
        if (.not.NEFMatB%isBlockZero(1,col)) &
          call NEFMatOut%iput(matmul(IntegralA,NEFMatB%iat(1,col)),1,col)
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
      character(len=16)::outSymm
      type(mqc_nef_vector),intent(in),target::NEFVecA,NEFVecB
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
!     Operator label: equal labels are kept, mixed labels give 'asymmetric'.
      outSymm = NEFVecA%Symmetry
      if(NEFVecB%Symmetry.ne.outSymm) outSymm = 'asymmetric'
      call NEFVecOut%init(NEFVecA%Map_Matrix,&
           symmetry=outSymm,sameBraKet=(NEFVecA%sameBraKet.and.NEFVecB%sameBraKet))
      do col = 1, NEFVecB%NumSpace
        if (NEFVecA%isBlockZero(col).and.NEFVecB%isBlockZero(col)) then
          continue
        else if (NEFVecA%isBlockZero(col)) then
          call NEFVecOut%iput(NEFVecB%iref(col),col)
        else if (NEFVecB%isBlockZero(col)) then
          call NEFVecOut%iput(NEFVecA%iref(col),col)
        else
          call NEFVecOut%iput(NEFVecA%iref(col) + NEFVecB%iref(col),col)
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
        character(len=16)::outSymm
        integer::col,row
        type(mqc_scf_eigenvalues),dimension(:),allocatable::Map1,Map2
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
!      
!       Operator label: equal labels are kept, mixed labels give 'asymmetric'.
        outSymm = NEFMatA%Symmetry
        if(NEFMatB%Symmetry.ne.outSymm) outSymm = 'asymmetric'
        Map1 = NEFMatA%getMap(1)
        Map2 = NEFMatA%getMap(2)
        call NEFMatOut%init(Map1,initial_map_matrix_col=Map2,symmetry=outSymm,&
          sameBraKet=(NEFMatA%sameBraKet.and.NEFMatB%sameBraKet))
        do row = 1, NEFMatA%NumSpace_row
          do col = 1, NEFMatA%NumSpace_col
            if(NEFMatOut%Mirrored.and.col.gt.row) cycle
            if( .not.NEFMatA%isBlockZero(row,col).and..not.NEFMatB%isBlockZero(row,col)) then
              call NEFMatOut%iput(NEFMatA%iat(row,col)+NEFMatB%iat(row,col),row,col)
            elseif(NEFMatA%isBlockZero(row,col).and..not.NEFMatB%isBlockZero(row,col)) then
              call NEFMatOut%iput(NEFMatB%iat(row,col),row,col)
            elseif(NEFMatB%isBlockZero(row,col).and..not.NEFMatA%isBlockZero(row,col)) then
              call NEFMatOut%iput(NEFMatA%iat(row,col),row,col)
            end if
          end do
        end do
!
      end function mqc_nefmatrix_sum
!
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
      type(mqc_nef_r4),intent(in),target::NEFR4A,NEFR4B
      type(mqc_nef_r4)::NEFR4Out
      integer::id1,id2,id3,id4
      character(len=16)::outSymm
!
      type(mqc_nef_r4),target::full_NEFR4A
      class(mqc_nef_r4),pointer::ref_NEFR4A
      type(mqc_nef_r4),target::full_NEFR4B
      class(mqc_nef_r4),pointer::ref_NEFR4B
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


!     Operator label: equal labels are kept, mixed labels give 'asymmetric'.
      outSymm = NEFR4A%Symmetry
      if(NEFR4B%Symmetry.ne.outSymm) outSymm = 'asymmetric'
      call NEFR4Out%init(NEFR4A%Map_Matrix_1,NEFR4A%Map_Matrix_2,&
           NEFR4A%Map_Matrix_3,NEFR4A%Map_Matrix_4,&
           symmetry=outSymm,sameBraKet=(NEFR4A%sameBraKet.and.NEFR4B%sameBraKet))
!     P2: inputs share the maps of a mirrored output, so its stored slots are read
!     in place; a mirrored input feeding a full-storage output is read from a
!     full-storage copy (iref holds no mirror block).
      ref_NEFR4A => NEFR4A
      if(NEFR4A%Mirrored.and..not.NEFR4Out%Mirrored) then
        full_NEFR4A = mqc_nef_r4_full_storage(NEFR4A)
        ref_NEFR4A => full_NEFR4A
      end if
      ref_NEFR4B => NEFR4B
      if(NEFR4B%Mirrored.and..not.NEFR4Out%Mirrored) then
        full_NEFR4B = mqc_nef_r4_full_storage(NEFR4B)
        ref_NEFR4B => full_NEFR4B
      end if
      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              if(NEFR4Out%isMirrorSlot(id1,id2,id3,id4)) cycle
              if (NEFR4A%isBlockZero(id1,id2,id3,id4).and.NEFR4B%isBlockZero(id1,id2,id3,id4)) then
                continue
              else if (NEFR4A%isBlockZero(id1,id2,id3,id4)) then
                call NEFR4Out%iput(ref_NEFR4B%iref(id1,id2,id3,id4),id1,id2,id3,id4)
              else if (NEFR4B%isBlockZero(id1,id2,id3,id4)) then
                call NEFR4Out%iput(ref_NEFR4A%iref(id1,id2,id3,id4),id1,id2,id3,id4)
              else
                call NEFR4Out%iput(ref_NEFR4A%iref(id1,id2,id3,id4) + ref_NEFR4B%iref(id1,id2,id3,id4),&
                     id1,id2,id3,id4)
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
      integer(kind=int64)::row,col
!
!     Operator label and bra-ket flag are unchanged.
      call NEFMatOut%init(NEFMat%Map_Matrix_row,&
           initial_map_matrix_col=NEFMat%Map_Matrix_col,&
           symmetry=NEFMat%Symmetry,sameBraKet=NEFMat%sameBraKet)
      do row = 1, NEFMatOut%NumSpace_row
        do col = 1, NEFMatOut%NumSpace_col
          if(NEFMatOut%Mirrored.and.col.gt.row) cycle
          if (.not.NEFMat%isBlockZero(row,col)) &
            call NEFMatOut%iput(conjg(NEFMat%iat(row,col)),row,col)
        end do
      end do
!
      end function mqc_NEFMatrix_conjugate
!
!     PROCEDURE MQC_NEF_Vector_Conjugate
!
!>    \brief <b> MQC_NEF_Vector_Conjugate is used to compute the conjugate
!>    of a MQC NEF Vector type variable.
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVec
!>    \verbatim
!>        NEFVec is type(mqc_nef_vector)
!>        The MQC NEF Vector type variable to conjugate.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_conjugate(NEFVec) result(NEFVecOut)
!
      implicit none
      type(mqc_nef_vector),intent(in),target::NEFVec
      type(mqc_nef_vector)::NEFVecOut
      integer::id
!
      call NEFVecOut%init(NEFVec%Map_Matrix,symmetry=NEFVec%Symmetry,&
           sameBraKet=NEFVec%sameBraKet)
      do id = 1, NEFVec%NumSpace
        if (.not.NEFVec%isBlockZero(id)) then
          call NEFVecOut%iput(conjg(NEFVec%iref(id)),id)
        end if
      end do
!
      end function mqc_nef_vector_conjugate
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
      integer(kind=int64)::row,col
!
!     Operator label and bra-ket flag are unchanged.
      call NEFMatOut%init(NEFMat%Map_Matrix_col,&
           initial_map_matrix_col=NEFMat%Map_Matrix_row,&
           symmetry=NEFMat%Symmetry,sameBraKet=NEFMat%sameBraKet)
      do row = 1, NEFMatOut%NumSpace_row
        do col = 1, NEFMatOut%NumSpace_col
          if(NEFMatOut%Mirrored.and.col.gt.row) cycle
          if (.not.NEFMat%isBlockZero(col,row)) &
            call NEFMatOut%iput(dagger(NEFMat%iat(col,row)),row,col)
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
      type(mqc_nef_r4),intent(in),target::NEFR4
      type(mqc_nef_r4)::NEFR4Out
      integer(kind=int64)::id1,id2,id3,id4
!
      type(mqc_nef_r4),target::full_NEFR4
      class(mqc_nef_r4),pointer::ref_NEFR4
!     P2: with the object symmetry dagger(O) = +/- O (block (S,R,Q,P) is the +/- dagger
!     of block (P,Q,R,S)): a scaled copy of the stored slots, no per-block dagger.
      if(NEFR4%Mirrored) then
        if(NEFR4%Symmetry.eq.'hermitian') then
          NEFR4Out = 1.0*NEFR4
        else
          NEFR4Out = (-1.0)*NEFR4
        end if
      else
  !     Mirror slots of a symmetric input are read from a full-storage copy (iref
  !     holds no mirror block); full-storage inputs are read in place.
        ref_NEFR4 => NEFR4
        if(NEFR4%Mirrored) then
          full_NEFR4 = mqc_nef_r4_full_storage(NEFR4)
          ref_NEFR4 => full_NEFR4
        end if
        call NEFR4Out%init(NEFR4%Map_Matrix_4,NEFR4%Map_Matrix_3,&
             NEFR4%Map_Matrix_2,NEFR4%Map_Matrix_1,&
             symmetry=NEFR4%Symmetry,sameBraKet=NEFR4%sameBraKet)

        do id1 = 1, NEFR4Out%NumSpace_1
          do id2 = 1, NEFR4Out%NumSpace_2
            do id3 = 1, NEFR4Out%NumSpace_3
              do id4 = 1, NEFR4Out%NumSpace_4
                if (.not.NEFR4%isBlockZero(id4,id3,id2,id1)) &
                  call NEFR4Out%iput(dagger(ref_NEFR4%iref(id4,id3,id2,id1)),id1,id2,id3,id4)
              end do
            end do
          end do
        end do
      end if
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
      integer(kind=int64)::row,col
!
!     Operator label and bra-ket flag are unchanged.
      call NEFMatOut%init(NEFMat%Map_Matrix_col,&
           initial_map_matrix_col=NEFMat%Map_Matrix_row,&
           symmetry=NEFMat%Symmetry,sameBraKet=NEFMat%sameBraKet)
      do row = 1, NEFMatOut%NumSpace_row
        do col = 1, NEFMatOut%NumSpace_col
          if(NEFMatOut%Mirrored.and.col.gt.row) cycle
          if (.not.NEFMat%isBlockZero(col,row)) &
            call NEFMatOut%iput(transpose(NEFMat%iat(col,row)),row,col)
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
      type(mqc_nef_vector),intent(in),target::NEFVecB
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
!     Operator label: real c keeps it, imaginary c swaps herm/antiherm.
      call NEFVecOut%init(NEFVecB%Map_Matrix,&
           symmetry=mqc_nef_scalar_symmetry(NEFVecB%Symmetry,scalar),sameBraKet=NEFVecB%sameBraKet)
      if (abs(scalar).ne.0.0) then
        do row = 1, NEFVecB%NumSpace
          if (.not.NEFVecB%isBlockZero(row)) &
            call NEFVecOut%iput(NEFVecB%iref(row)*scalar,row)
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
      type(mqc_nef_vector),intent(in),target::NEFVecA
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
!     Operator label: real c keeps it, imaginary c swaps herm/antiherm.
      call NEFVecOut%init(NEFVecA%Map_Matrix,&
           symmetry=mqc_nef_scalar_symmetry(NEFVecA%Symmetry,scalar),sameBraKet=NEFVecA%sameBraKet)
      if (abs(scalar).ne.0.0) then
        do row = 1, NEFVecA%NumSpace
          if (.not.NEFVecA%isBlockZero(row)) &
            call NEFVecOut%iput(scalar*NEFVecA%iref(row),row)
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
        integer::col,row,id,maxCol
        Real(kind=real64)::Thresh
        type(mqc_scf_eigenvalues),dimension(:),allocatable::Map1,Map2
!
        If(MQC_very_small.gt.0.0) then
          Thresh = MQC_very_small
        Else
          Thresh = 1e-10
        EndIf
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
!       Operator label: real c keeps it, imaginary c swaps herm/antiherm.
        Map1 = NEFMat%getMap(1)
        Map2 = NEFMat%getMap(2)
        call NEFMatOut%init(Map1,initial_map_matrix_col=Map2,&
          symmetry=mqc_nef_scalar_symmetry(NEFMat%Symmetry,scalar),&
          sameBraKet=NEFMat%sameBraKet)
        if(abs(scalar).ge.Thresh) then
          do row = 1, NEFMat%NumSpace_row
            do col = 1, NEFMat%NumSpace_col
              if(NEFMatOut%Mirrored.and.col.gt.row) cycle
              if(NEFMat%isBlockZero(row,col)) cycle
              call NEFMatOut%iput(scalar*NEFMat%iat(row,col),row,col)
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
      NEFMatOut = scalar*NEFMat
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
      type(mqc_nef_r4),intent(in),target::NEFR4
      class(*),intent(in)::scalarIn
      type(mqc_scalar)::scalar
      type(mqc_nef_r4)::NEFR4Out
      integer::id1,id2,id3,id4
!
      type(mqc_nef_r4),target::full_NEFR4
      class(mqc_nef_r4),pointer::ref_NEFR4
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
!     Operator label: real c keeps it, imaginary c swaps herm/antiherm.
      call NEFR4Out%init(NEFR4%Map_Matrix_1,NEFR4%Map_Matrix_2,&
           NEFR4%Map_Matrix_3,NEFR4%Map_Matrix_4,&
           symmetry=mqc_nef_scalar_symmetry(NEFR4%Symmetry,scalar),&
           sameBraKet=NEFR4%sameBraKet)
!     P2: inputs share the maps of a mirrored output, so its stored slots are read
!     in place; a mirrored input feeding a full-storage output is read from a
!     full-storage copy (iref holds no mirror block).
      ref_NEFR4 => NEFR4
      if(NEFR4%Mirrored.and..not.NEFR4Out%Mirrored) then
        full_NEFR4 = mqc_nef_r4_full_storage(NEFR4)
        ref_NEFR4 => full_NEFR4
      end if
      if (abs(scalar).ne.0.0) then
        do id1 = 1, NEFR4%NumSpace_1
          do id2 = 1, NEFR4%NumSpace_2
            do id3 = 1, NEFR4%NumSpace_3
              do id4 = 1, NEFR4%NumSpace_4
                if(NEFR4Out%isMirrorSlot(id1,id2,id3,id4)) cycle
                if (.not.NEFR4%isBlockZero(id1,id2,id3,id4)) &
                  call NEFR4Out%iput(scalar*ref_NEFR4%iref(id1,id2,id3,id4),id1,id2,id3,id4)
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
      type(mqc_nef_r4),intent(in),target::NEFR4
      class(*),intent(in)::scalarB
      type(mqc_scalar)::scalar
      type(mqc_nef_r4)::NEFR4Out
      integer::id1,id2,id3,id4
!
      type(mqc_nef_r4),target::full_NEFR4
      class(mqc_nef_r4),pointer::ref_NEFR4
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
!     Operator label: real c keeps it, imaginary c swaps herm/antiherm.
      call NEFR4Out%init(NEFR4%Map_Matrix_1,NEFR4%Map_Matrix_2,&
           NEFR4%Map_Matrix_3,NEFR4%Map_Matrix_4,&
           symmetry=mqc_nef_scalar_symmetry(NEFR4%Symmetry,scalar),&
           sameBraKet=NEFR4%sameBraKet)
!     P2: inputs share the maps of a mirrored output, so its stored slots are read
!     in place; a mirrored input feeding a full-storage output is read from a
!     full-storage copy (iref holds no mirror block).
      ref_NEFR4 => NEFR4
      if(NEFR4%Mirrored.and..not.NEFR4Out%Mirrored) then
        full_NEFR4 = mqc_nef_r4_full_storage(NEFR4)
        ref_NEFR4 => full_NEFR4
      end if
      if (abs(scalar).ne.0.0) then
        do id1 = 1, NEFR4%NumSpace_1
          do id2 = 1, NEFR4%NumSpace_2
            do id3 = 1, NEFR4%NumSpace_3
              do id4 = 1, NEFR4%NumSpace_4
                if(NEFR4Out%isMirrorSlot(id1,id2,id3,id4)) cycle
                if (.not.NEFR4%isBlockZero(id1,id2,id3,id4)) &
                  call NEFR4Out%iput(ref_NEFR4%iref(id1,id2,id3,id4)*scalar,id1,id2,id3,id4)
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
      character(len=16)::outSymm
      type(mqc_nef_vector),intent(in),target::NEFVecA,NEFVecB
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
!     Operator label: equal labels are kept, mixed labels give 'asymmetric'.
      outSymm = NEFVecA%Symmetry
      if(NEFVecB%Symmetry.ne.outSymm) outSymm = 'asymmetric'
      call NEFVecOut%init(NEFVecA%Map_Matrix,&
           symmetry=outSymm,sameBraKet=(NEFVecA%sameBraKet.and.NEFVecB%sameBraKet))
      do col = 1, NEFVecB%NumSpace
        if (NEFVecA%isBlockZero(col).and.NEFVecB%isBlockZero(col)) then
          continue
        else if (NEFVecA%isBlockZero(col)) then
          call NEFVecOut%iput((-1.0)*NEFVecB%iref(col),col)
        else if (NEFVecB%isBlockZero(col)) then
          call NEFVecOut%iput(NEFVecA%iref(col),col)
        else
          call NEFVecOut%iput(NEFVecA%iref(col) - NEFVecB%iref(col),col)
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
        integer::col,row
        character(len=16)::outSymm
        type(mqc_scf_eigenvalues),dimension(:),allocatable::Map1,Map2
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
!       Operator label: equal labels are kept, mixed labels give 'asymmetric'.
        outSymm = NEFMatA%Symmetry
        if(NEFMatB%Symmetry.ne.outSymm) outSymm = 'asymmetric'
        Map1 = NEFMatA%getMap(1)
        Map2 = NEFMatA%getMap(2)
        call NEFMatOut%init(Map1,initial_map_matrix_col=Map2,symmetry=outSymm,&
          sameBraKet=(NEFMatA%sameBraKet.and.NEFMatB%sameBraKet))
        do row = 1, NEFMatA%NumSpace_row
          do col = 1, NEFMatA%NumSpace_col
            if(NEFMatOut%Mirrored.and.col.gt.row) cycle
            if( .not.NEFMatA%isBlockZero(row,col).and..not.NEFMatB%isBlockZero(row,col)) then
              call NEFMatOut%iput(NEFMatA%iat(row,col)-NEFMatB%iat(row,col),row,col)
            elseif(NEFMatA%isBlockZero(row,col).and..not.NEFMatB%isBlockZero(row,col)) then
              call NEFMatOut%iput((-1)*NEFMatB%iat(row,col),row,col)
            elseif(NEFMatB%isBlockZero(row,col).and..not.NEFMatA%isBlockZero(row,col)) then
              call NEFMatOut%iput(NEFMatA%iat(row,col),row,col)
            end if
          end do
        end do
!
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
      type(mqc_nef_r4),intent(in),target::NEFR4A,NEFR4B
      type(mqc_nef_r4)::NEFR4Out
      integer::id1,id2,id3,id4
      character(len=16)::outSymm
!
      type(mqc_nef_r4),target::full_NEFR4A
      class(mqc_nef_r4),pointer::ref_NEFR4A
      type(mqc_nef_r4),target::full_NEFR4B
      class(mqc_nef_r4),pointer::ref_NEFR4B
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


!     Operator label: equal labels are kept, mixed labels give 'asymmetric'.
      outSymm = NEFR4A%Symmetry
      if(NEFR4B%Symmetry.ne.outSymm) outSymm = 'asymmetric'
      call NEFR4Out%init(NEFR4A%Map_Matrix_1,NEFR4A%Map_Matrix_2,&
           NEFR4A%Map_Matrix_3,NEFR4A%Map_Matrix_4,&
           symmetry=outSymm,sameBraKet=(NEFR4A%sameBraKet.and.NEFR4B%sameBraKet))
!     P2: inputs share the maps of a mirrored output, so its stored slots are read
!     in place; a mirrored input feeding a full-storage output is read from a
!     full-storage copy (iref holds no mirror block).
      ref_NEFR4A => NEFR4A
      if(NEFR4A%Mirrored.and..not.NEFR4Out%Mirrored) then
        full_NEFR4A = mqc_nef_r4_full_storage(NEFR4A)
        ref_NEFR4A => full_NEFR4A
      end if
      ref_NEFR4B => NEFR4B
      if(NEFR4B%Mirrored.and..not.NEFR4Out%Mirrored) then
        full_NEFR4B = mqc_nef_r4_full_storage(NEFR4B)
        ref_NEFR4B => full_NEFR4B
      end if
      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              if(NEFR4Out%isMirrorSlot(id1,id2,id3,id4)) cycle
              if (NEFR4A%isBlockZero(id1,id2,id3,id4).and.NEFR4B%isBlockZero(id1,id2,id3,id4)) then
                continue
              else if (NEFR4A%isBlockZero(id1,id2,id3,id4)) then
                call NEFR4Out%iput((-1.0)*ref_NEFR4B%iref(id1,id2,id3,id4),id1,id2,id3,id4)
              else if (NEFR4B%isBlockZero(id1,id2,id3,id4)) then
                call NEFR4Out%iput(ref_NEFR4A%iref(id1,id2,id3,id4),id1,id2,id3,id4)
              else
                call NEFR4Out%iput(ref_NEFR4A%iref(id1,id2,id3,id4) - ref_NEFR4B%iref(id1,id2,id3,id4),&
                     id1,id2,id3,id4)
              end if
           end do
          end do
        end do
      end do
!
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
      class(mqc_nef_r4),target::NEFR4
      integer(kind=int64),intent(in)::iOut
      character(len=*),intent(in)::header
      character(len=64)::arrayType
      logical,intent(In),optional::blank_at_top,blank_at_bottom
      character(len=*),intent(in),optional::formatStr
      Real(kind=real64),intent(in),optional::short_print_threshold
      Real(kind=real64)::print_threshold
      integer::id1,id2,id3,id4
      type(mqc_twoeris),pointer::blk
!
 type(mqc_nef_r4),target::full_NEFR4
 class(mqc_nef_r4),pointer::ref_NEFR4
!     Mirror slots of a symmetric input are read from a full-storage copy (iref
!     holds no mirror block); full-storage inputs are read in place.
 ref_NEFR4 => NEFR4
 if(NEFR4%Mirrored) then
   full_NEFR4 = mqc_nef_r4_full_storage(NEFR4)
   ref_NEFR4 => full_NEFR4
 end if
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
              if (NEFR4%isBlockZero(id1,id2,id3,id4)) then
                write(iOut,1000) &
   'Sub-Space ('//trim(num2char(id1))//trim(num2char(id2))//'|'//trim(num2char(id3))//trim(num2char(id4))//') '&
               &' (all zero)'
              else
              blk => ref_NEFR4%iref(id1,id2,id3,id4)
              call mqc_print(blk,iOut,&
   'Sub-Space ('//trim(num2char(id1))//trim(num2char(id2))//'|'//trim(num2char(id3))//trim(num2char(id4))//') '&
               &' '//trim(blk%storType())//' '//trim(blk%type()),&
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

      integer::id,id1,id2,id3,id4,N1,N12,N123,N1234
      type(mqc_scf_integral)::blockA,blockB

      call NEFR4Out%init(NEFMatA%Map_Matrix_Row,NEFMatA%Map_Matrix_Col,&
           NEFMatB%Map_Matrix_Row,NEFMatB%Map_Matrix_Col,&
           symmetry='asymmetric',sameBraKet=(NEFMatA%sameBraKet.and.NEFMatB%sameBraKet))
      N1 = NEFR4Out%NumSpace_4
      N12 = N1*NEFR4Out%NumSpace_3
      N123 = N12*NEFR4Out%NumSpace_2
      N1234 = N123*NEFR4Out%NumSpace_1

      do id1 = 1, NEFR4Out%NumSpace_1
        do id2 = 1, NEFR4Out%NumSpace_2
          if (NEFMatA%isBlockZero(id1,id2)) cycle
          blockA = NEFMatA%iat(id1,id2)
          do id3 = 1, NEFR4Out%NumSpace_3
            do id4 = 1, NEFR4Out%NumSpace_4
              if (.not.NEFMatB%isBlockZero(id3,id4)) then
                blockB = NEFMatB%iat(id3,id4)
                call NEFR4Out%iput(blockA%outer2ERI(blockB),id1,id2,id3,id4)
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
      class(mqc_nef_vector),intent(in),target::NEFVecA
      type(mqc_nef_vector),intent(in),target::NEFVecB
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64)::id1,id2
      type(mqc_scf_eigenvalues),pointer::blkA
!
      NEFMatOut%NumSpace_row = NEFVecA%NumSpace
      NEFMatOut%NumSpace_col = NEFVecB%NumSpace
      call mqc_nef_matrix_allocate_blocks(NEFMatOut)
!     Built here without init (no zero blocks allocated): the builder sets the
!     flags at construction; an operator product is 'asymmetric' (default label).
      NEFMatOut%sameBraKet = (NEFVecA%sameBraKet.and.NEFVecB%sameBraKet)
      NEFMatOut%Map_Matrix_Row = NEFVecA%Map_Matrix
      NEFMatOut%Map_Matrix_Col = NEFVecB%Map_Matrix
      NEFMatOut%Size_vector_alpha_row = NEFVecA%Size_vector_alpha
      NEFMatOut%Size_vector_beta_row = NEFVecA%Size_vector_beta
      NEFMatOut%Size_vector_alpha_col = NEFVecB%Size_vector_alpha
      NEFMatOut%Size_vector_beta_col = NEFVecB%Size_vector_beta
!
      do id1 = 1, NEFMatOut%NumSpace_row
        do id2 = 1, NEFMatOut%NumSpace_col
          if (.not.(NEFVecA%isBlockZero(id1).or.NEFVecB%isBlockZero(id2))) then
            blkA => NEFVecA%iref(id1)
            call NEFMatOut%iput(blkA%outer(NEFVecB%iref(id2)),id1,id2)
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

      integer::id,id1,id2,id3,id4,N1,N12,N123,N1234
      type(mqc_scf_integral)::blockA,blockB

      call NEFR4Out%init(NEFMatA%Map_Matrix_Row,NEFMatB%Map_Matrix_Col,&
           NEFMatB%Map_Matrix_Row,NEFMatA%Map_Matrix_Col,&
           symmetry='asymmetric',sameBraKet=(NEFMatA%sameBraKet.and.NEFMatB%sameBraKet))
      N1 = NEFR4Out%NumSpace_4
      N12 = N1*NEFR4Out%NumSpace_3
      N123 = N12*NEFR4Out%NumSpace_2
      N1234 = N123*NEFR4Out%NumSpace_1

      do id1 = 1, NEFR4Out%NumSpace_1
        do id4 = 1, NEFR4Out%NumSpace_4
          if (NEFMatA%isBlockZero(id1,id4)) cycle
          blockA = NEFMatA%iat(id1,id4)
          do id3 = 1, NEFR4Out%NumSpace_3
            do id2 = 1, NEFR4Out%NumSpace_2
              if (.not.NEFMatB%isBlockZero(id3,id2)) then
                blockB = NEFMatB%iat(id3,id2)
                call NEFR4Out%iput(blockA%crossouter2ERI(blockB),id1,id2,id3,id4)
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
      type(mqc_nef_vector),intent(in),target::NEFVecA,NEFVecB
      type(mqc_scalar)::output
      integer::id
      type(mqc_scf_eigenvalues),pointer::blkA,blkB
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
        if(NEFVecA%isBlockZero(id).or.NEFVecB%isBlockZero(id)) cycle
        blkA => NEFVecA%iref(id)
        blkB => NEFVecB%iref(id)
        output = output + (Transpose(blkA%getBlock()).dot.blkB%getBlock())
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
          if (NEFMatA%isBlockZero(row,col).or.NEFMatB%isBlockZero(row,col)) cycle
          output = output + Contraction(NEFMatA%iat(row,col),NEFMatB%iat(row,col))
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
      type(mqc_nef_r4),intent(in),target::NEFR4A,NEFR4B
      character(len=*),optional,intent(in)::label
      character(len=64)::myLabel
      type(mqc_scalar)::output
      type(mqc_scalar)::cpair,sAB
      logical::pairTrick
      integer::id1,id2,id3,id4
      type(mqc_twoeris):: debbug_ERI
      type(mqc_twoeris),pointer::blkB
!
      type(mqc_twoeris),target::scr_NEFR4A
      type(mqc_twoeris),target::scr_NEFR4B
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
!     C2 applies when both inputs are mirrored with the same pairing (equal maps, so the
!     same stored slots); sAB = +/-1 from the two labels.
      pairTrick = (NEFR4A%Mirrored.and.NEFR4B%Mirrored)
      if(pairTrick) pairTrick = mqc_nef_maps_equal(NEFR4A%Map_Matrix_1,NEFR4B%Map_Matrix_1).and.&
                                mqc_nef_maps_equal(NEFR4A%Map_Matrix_2,NEFR4B%Map_Matrix_2)
      sAB = 1.0
      if(NEFR4A%Symmetry.ne.NEFR4B%Symmetry) sAB = -1.0
      output = 0.0
      do id1 = 1, NEFR4A%NumSpace_1
        do id2 = 1, NEFR4A%NumSpace_2
          do id3 = 1, NEFR4A%NumSpace_3
            do id4 = 1, NEFR4A%NumSpace_4
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
                if (pairTrick) then
!                 C2: both mirrored with the same pairing: a dagger pair (Q, Q*) gives
!                 c + sA*sB*conj(c) from the stored blocks; self-dagger slots give c
                  if (.not.(NEFR4A%isMirrorSlot(id1,id2,id3,id4).or.&
                            NEFR4A%isBlockZero(id1,id2,id3,id4).or.NEFR4B%isBlockZero(id1,id2,id3,id4))) then
                    cpair = contraction(NEFR4A%iref(id1,id2,id3,id4),NEFR4B%iref(id1,id2,id3,id4),'coulomb')
                    if (id1.eq.id4.and.id2.eq.id3) then
                      output = output + cpair
                    else
                      output = output + cpair + sAB*conjg(cpair)
                    end if
                  end if
                else if (.not.(NEFR4A%isBlockZero(id1,id2,id3,id4).or.NEFR4B%isBlockZero(id1,id2,id3,id4))) then
                  output = output + contraction(mqc_nef_r4_bref(NEFR4A,id1,id2,id3,id4,scr_NEFR4A),&
                                                mqc_nef_r4_bref(NEFR4B,id1,id2,id3,id4,scr_NEFR4B),'coulomb')
                end if
              end if
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
                 if (.not.(NEFR4A%isBlockZero(id1,id4,id3,id2).or.NEFR4B%isBlockZero(id1,id2,id3,id4))) then
                   blkB => mqc_nef_r4_bref(NEFR4B,id1,id2,id3,id4,scr_NEFR4B)
                   output = output - &
                     contraction(mqc_nef_r4_bref(NEFR4A,id1,id4,id3,id2,scr_NEFR4A),blkB%swap(2,4),'coulomb')
                 end if
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
      type(mqc_nef_r4),intent(in),target::NEFR4A
      type(mqc_nef_matrix),intent(in)::NEFMatB
      character(len=*),optional,intent(in)::label
      integer(kind=int64),dimension(4),intent(in),optional::indexA
      integer(kind=int64),dimension(4)::my_indexA
      character(len=64)::myLabel
      type(mqc_nef_matrix)::output
      integer::id1,id2,id3,id4
      type(mqc_twoeris),pointer::blkA
      type(mqc_scf_integral)::blockB,blockOut
      logical::started
      type(mqc_twoeris):: debbug_ERI
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la
!
      type(mqc_twoeris),target::scr_NEFR4A
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
      call output%init(Map_Matrix_Row,initial_map_matrix_col=Map_Matrix_Col,&
           sameBraKet=(NEFR4A%sameBraKet.and.NEFMatB%sameBraKet))
      do id4 = 1, NEFR4A%blockSize(la)
        do id3 = 1, NEFR4A%blockSize(ka)
          started = .false.
          do id2 = 1, NEFR4A%blockSize(ja)
            do id1 = 1, NEFR4A%blockSize(ia)
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'coulomb') then
                if ((NEFR4A%blockSize(ia,id1,'alpha')+NEFR4A%blockSize(ia,id1,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ja,id2,'alpha')+NEFR4A%blockSize(ja,id2,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ka,id3,'alpha')+NEFR4A%blockSize(ka,id3,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(la,id4,'alpha')+NEFR4A%blockSize(la,id4,'beta')).eq.0) cycle
                if (.not.(NEFR4A%isBlockZero([ia,id1],[ja,id2],[ka,id3],[la,id4]).or.&
                          NEFMatB%isBlockZero(id1,id2))) then
                  blockB = NEFMatB%iat(id1,id2)
                  blkA => mqc_nef_r4_bref_pairs(NEFR4A,[ia,id1],[ja,id2],[ka,id3],[la,id4],scr_NEFR4A)
                  if (ka.gt.la) then
                    if (.not.started) then
                      blockOut = &
                          transpose(contraction(blkA,blockB,[ia,ja]))
                    else
                      blockOut = blockOut &
                          + transpose(contraction(blkA,blockB,[ia,ja]))
                    end if
                  else
                    if (.not.started) then
                      blockOut = &
                          contraction(blkA,blockB,[ia,ja])
                    else
                      blockOut = blockOut &
                          + contraction(blkA,blockB,[ia,ja])
                    end if
                  end if
                  started = .true.
                end if
              end if
              if (myLabel.eq.'doublebar'.or.myLabel.eq.'exchange') then
                if ((NEFR4A%blockSize(ia,id1,'alpha')+NEFR4A%blockSize(ia,id1,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(la,id2,'alpha')+NEFR4A%blockSize(la,id2,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ka,id3,'alpha')+NEFR4A%blockSize(ka,id3,'beta')).eq.0.or.&
                    (NEFR4A%blockSize(ja,id4,'alpha')+NEFR4A%blockSize(ja,id4,'beta')).eq.0) cycle
                if (.not.(NEFR4A%isBlockZero([ia,id1],[la,id2],[ka,id3],[ja,id4]).or.&
                          NEFMatB%isBlockZero(id1,id2))) then
                  blockB = NEFMatB%iat(id1,id2)
                  blkA => mqc_nef_r4_bref_pairs(NEFR4A,[ia,id1],[la,id2],[ka,id3],[ja,id4],scr_NEFR4A)
                  if (.not.started) then
                    blockOut = &
                        (-1.0)*contraction(blkA,blockB,[ia,la])
                  else
                    blockOut = blockOut &
                        - contraction(blkA,blockB,[ia,la])
                  end if
                  started = .true.
                end if
              end if
            end do
          end do
          if (started) call output%iput(blockOut,id3,id4)
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
      integer(kind=int64)::row,col,id,maxCol
      integer(kind=int64),dimension(:),Allocatable::rowSel,colSel
      logical::full_row,full_col,lowerOnly
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
!     Operator label and bra-ket flag are kept (O8: index changes never affect them).
      call NEFMatOut%init(New_Map_Row,&
          initial_map_matrix_col=New_Map_Col,&
          symmetry=NEFMat%Symmetry,sameBraKet=NEFMat%sameBraKet)
      allocate(rowSel(size(New_Map_Row)),colSel(size(New_Map_Col)))
      if(full_row) then
        rowSel = [(id, id = 1, size(New_Map_Row))]
      else
        rowSel = my_row
      end if
      if(full_col) then
        colSel = [(id, id = 1, size(New_Map_Col))]
      else
        colSel = my_col
      end if
      lowerOnly = NEFMatOut%Mirrored
      do row = 1, size(rowSel)
        maxCol = size(colSel)
        if(lowerOnly) maxCol = row
        do col = 1, maxCol
          if (NEFMat%isBlockZero(rowSel(row),colSel(col))) cycle
          call NEFMatOut%iput(NEFMat%iat(rowSel(row),colSel(col)),row,col)
        end do
      end do
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
      type(mqc_nef_r4),intent(in),target::NEFR4A,NEFR4B
      type(mqc_twoeris)::acc
      logical::started
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_r4)::NEFR4Out
      type(mqc_scf_eigenvalues),dimension(:),allocatable::outMap1,outMap2,outMap3,outMap4
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      Integer(kind=int64)::idiab,idjab,A_index,B_index,out_index,N1234
!
      type(mqc_twoeris),target::scr_NEFR4A
      type(mqc_twoeris),target::scr_NEFR4B
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

      outMap1 = NEFR4A%blockMap(ka)
      outMap2 = NEFR4A%blockMap(la)
      outMap3 = NEFR4B%blockMap(kb)
      outMap4 = NEFR4B%blockMap(lb)
      call NEFR4Out%init(outMap1,outMap2,outMap3,outMap4,&
           symmetry='asymmetric',sameBraKet=(NEFR4A%sameBraKet.and.NEFR4B%sameBraKet))
      N1234 = NEFR4Out%NumSpace_1*NEFR4Out%NumSpace_2*NEFR4Out%NumSpace_3*NEFR4Out%NumSpace_4

      do id1 = 1, nefr4Out%blockSize(1)
        do id2 = 1, nefr4Out%blockSize(2)
          do id3 = 1, nefr4Out%blockSize(3)
            do id4 = 1, nefr4Out%blockSize(4)
              started = .false.
              do idiab = 1, nefr4B%blockSize(ib)
                do idjab = 1, nefr4B%blockSize(jb)
                  if (NEFR4A%isBlockZero([ka,id1],[la,id2],[ia,idiab],[ja,idjab]).or.&
                      NEFR4B%isBlockZero([kb,id3],[lb,id4],[ib,idiab],[jb,idjab])) cycle
                  if(.not.started) then
                    acc = mqc_eri_eri_partial_contraction_2eri(&
                      mqc_nef_r4_bref_pairs(NEFR4A,[ka,id1],[la,id2],[ia,idiab],[ja,idjab],scr_NEFR4A),&
                      mqc_nef_r4_bref_pairs(NEFR4B,[kb,id3],[lb,id4],[ib,idiab],[jb,idjab],scr_NEFR4B),my_index1,my_index2)
                    started = .true.
                  else
                    acc = acc + mqc_eri_eri_partial_contraction_2eri(&
                      mqc_nef_r4_bref_pairs(NEFR4A,[ka,id1],[la,id2],[ia,idiab],[ja,idjab],scr_NEFR4A),&
                      mqc_nef_r4_bref_pairs(NEFR4B,[kb,id3],[lb,id4],[ib,idiab],[jb,idjab],scr_NEFR4B),my_index1,my_index2)
                  end if
                end do
              end do
              if(started) call NEFR4Out%iput(acc,id1,id2,id3,id4)
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
      type(mqc_nef_r4),intent(in),target::NEFR4A
      type(mqc_twoeris),intent(in)::ERIB
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_r4)::NEFR4Out
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      type(mqc_scalar)::val
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      Integer(kind=int64)::idiab,idjab,A_index,B_index,out_index,N1234
!
      type(mqc_twoeris),target::scr_NEFR4A
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

      call mqc_nef_r4_allocate_slots(NEFR4Out)
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
          if (.not.NEFR4A%isBlockZero([ka,id1],[la,id2],[ia,1],[ja,1])) &
            call NEFR4Out%iput(mqc_eri_eri_partial_contraction_2eri(&
              mqc_nef_r4_bref_pairs(NEFR4A,[ka,id1],[la,id2],[ia,1],[ja,1],scr_NEFR4A),ERIB,my_index1,my_index2),&
              id1,id2,1,1)
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
      type(mqc_nef_r4),intent(in),target::NEFR4B
      type(mqc_twoeris),intent(in)::ERIA
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_r4)::NEFR4Out
      type(mqc_scalar)::val
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      Integer(kind=int64)::idiab,idjab,A_index,B_index,out_index,N1234
!
      type(mqc_twoeris),target::scr_NEFR4B
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

      call mqc_nef_r4_allocate_slots(NEFR4Out)
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
          if (.not.NEFR4B%isBlockZero([kb,id1],[lb,id2],[ib,1],[jb,1])) &
            call NEFR4Out%iput(mqc_eri_eri_partial_contraction_2eri(&
              ERIA,mqc_nef_r4_bref_pairs(NEFR4B,[kb,id1],[lb,id2],[ib,1],[jb,1],scr_NEFR4B),my_index1,my_index2),&
              1,1,id1,id2)
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
      type(mqc_nef_r4),intent(in),target::NEFR4A,NEFR4B
      integer(kind=int64),dimension(4),intent(in),optional::index1,index2
      integer(kind=int64),dimension(4)::my_index1,my_index2
      type(mqc_nef_matrix)::MatOut
      type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::Map_Matrix_Col,Map_Matrix_Row
      Integer(kind=int64)::ia,ja,ka,la,ib,jb,kb,lb,i,j,k,l,id1,id2,id3,id4
      type(mqc_scf_integral)::blockOut
      logical::started
!
      type(mqc_twoeris),target::scr_NEFR4A
      type(mqc_twoeris),target::scr_NEFR4B
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
      call MatOut%init(Map_Matrix_Row,initial_map_matrix_col=Map_Matrix_Col,&
           sameBraKet=(NEFR4A%sameBraKet.and.NEFR4B%sameBraKet))
      do id3 = 1, nefr4A%blockSize(ka)
        do id4 = 1, nefr4A%blockSize(la)
          started = .false.
          do id1 = 1, nefr4A%blockSize(ia)
            do id2 = 1, nefr4A%blockSize(ja)
              if (NEFR4A%isBlockZero([ka,id3],[la,id4],[ia,id1],[ja,id2]).or.&
                  NEFR4B%isBlockZero([kb,id3],[lb,id4],[ib,id1],[jb,id2])) cycle
              if (.not.started) then
                blockOut = &
                  partialcontraction(mqc_nef_r4_bref_pairs(NEFR4A,[ka,id3],[la,id4],[ia,id1],[ja,id2],scr_NEFR4A),&
                  mqc_nef_r4_bref_pairs(NEFR4B,[kb,id3],[lb,id4],[ib,id1],[jb,id2],scr_NEFR4B),my_index1,my_index2)
                started = .true.
              else
                blockOut = blockOut +&
                  partialcontraction(mqc_nef_r4_bref_pairs(NEFR4A,[ka,id3],[la,id4],[ia,id1],[ja,id2],scr_NEFR4A),&
                  mqc_nef_r4_bref_pairs(NEFR4B,[kb,id3],[lb,id4],[ib,id1],[jb,id2],scr_NEFR4B),my_index1,my_index2)
              end if
            end do
          end do
          if (started) call MatOut%iput(blockOut,id3,id4)
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
!     PROCEDURE MQC_NEF_R4_ERIs_Put
!
!>    \brief <b> MQC_NEF_R4_ERIs_Put (iput) stores the MQC_TwoERIs block of the space
!>    block (P,Q,R,S) of an MQC NEF R4 variable. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    iput is the only block writer of MQC_NEF_R4: it allocates the block's
!>    MQC_NEF_R4_Slot on the first store and marks the block non-zero. The slot
!>    array itself must exist (mqc_nef_r4_allocate_slots). The block indices are
!>    given either as P,Q,R,S (the space indices on axes 1-4) or as four
!>    [axis,index] pairs, as blockIndex takes them.
!>
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_r4_eris_put(NEFR4,ERIs,P,Q,R,S)
!
      implicit none
      class(mqc_nef_r4),intent(inOut)::NEFR4
      type(mqc_twoeris),intent(in)::ERIs
      integer,intent(in)::P,Q,R,S
!
      call mqc_nef_r4_eris_put_pairs(NEFR4,ERIs,[1,P],[2,Q],[3,R],[4,S])
!
      end subroutine mqc_nef_r4_eris_put
!
      subroutine mqc_nef_r4_eris_put_pairs(NEFR4,ERIs,set1,set2,set3,set4)
!
      implicit none
      class(mqc_nef_r4),intent(inOut)::NEFR4
      type(mqc_twoeris),intent(in)::ERIs
      integer(kind=int64),dimension(2),intent(in)::set1,set2,set3,set4
      integer(kind=int64)::nef_bidx
      logical::isMirror
!
      if(.not.allocated(NEFR4%Zero_Flag)) call mqc_error(&
           'mqc_nef_r4_eris_put: NEF_R4 is open (never initialized)',6)
      call mqc_nef_r4_resolve(NEFR4,set1,set2,set3,set4,nef_bidx,isMirror)
      if(.not.isMirror) then
        NEFR4%Main_R4(nef_bidx)%R4Slot = ERIs
      else if(NEFR4%Symmetry.eq.'hermitian') then
!       Object symmetry: the shared slot holds the dagger-paired block.
        NEFR4%Main_R4(nef_bidx)%R4Slot = dagger(ERIs)
      else
        NEFR4%Main_R4(nef_bidx)%R4Slot = (-1.0)*dagger(ERIs)
      end if
      NEFR4%Zero_Flag(nef_bidx) = .false.
!
      end subroutine mqc_nef_r4_eris_put_pairs
!
!
!     PROCEDURE MQC_NEF_R4_ERIs_At
!
!>    \brief <b> MQC_NEF_R4_ERIs_At (iat) returns a copy of the MQC_TwoERIs block of
!>    the space block (P,Q,R,S) of an MQC NEF R4 variable. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    iat returns a copy of the stored block, or a zero block with the block
!>    dimensions (Size_vector_*) when the block is zero. Internal loops use iref
!>    (no copy); iat is meant for callers outside MQC_NEF.
!>
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_eris_at(NEFR4,P,Q,R,S) result(ERIsOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      integer,intent(in)::P,Q,R,S
      type(mqc_twoeris)::ERIsOut
!
      ERIsOut = mqc_nef_r4_eris_at_pairs(NEFR4,[1,P],[2,Q],[3,R],[4,S])
!
      end function mqc_nef_r4_eris_at
!
      function mqc_nef_r4_eris_at_pairs(NEFR4,set1,set2,set3,set4) result(ERIsOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      integer(kind=int64),dimension(2),intent(in)::set1,set2,set3,set4
      type(mqc_twoeris)::ERIsOut
      integer(kind=int64),dimension(4)::blk
      integer(kind=int64)::nef_bidx
      logical::isMirror
!
      if(.not.allocated(NEFR4%Zero_Flag)) call mqc_error(&
           'mqc_nef_r4_eris_at: NEF_R4 is open (never initialized)',6)
      call mqc_nef_r4_resolve(NEFR4,set1,set2,set3,set4,nef_bidx,isMirror)
      if(NEFR4%Zero_Flag(nef_bidx)) then
        blk(set1(1)) = set1(2)
        blk(set2(1)) = set2(2)
        blk(set3(1)) = set3(2)
        blk(set4(1)) = set4(2)
        call ERIsOut%init(NEFR4%Size_vector_alpha_1(blk(1)),NEFR4%Size_vector_beta_1(blk(1)),&
             nAlpha2=NEFR4%Size_vector_alpha_2(blk(2)),nBeta2=NEFR4%Size_vector_beta_2(blk(2)),&
             nAlpha3=NEFR4%Size_vector_alpha_3(blk(3)),nBeta3=NEFR4%Size_vector_beta_3(blk(3)),&
             nAlpha4=NEFR4%Size_vector_alpha_4(blk(4)),nBeta4=NEFR4%Size_vector_beta_4(blk(4)))
      else if(.not.isMirror) then
        ERIsOut = NEFR4%Main_R4(nef_bidx)%R4Slot
      else if(NEFR4%Symmetry.eq.'hermitian') then
        ERIsOut = dagger(NEFR4%Main_R4(nef_bidx)%R4Slot)
      else
        ERIsOut = (-1.0)*dagger(NEFR4%Main_R4(nef_bidx)%R4Slot)
      end if
!
      end function mqc_nef_r4_eris_at_pairs
!
!
!     PROCEDURE MQC_NEF_R4_ERIs_Ref
!
!>    \brief <b> MQC_NEF_R4_ERIs_Ref (iref) returns a read-only reference to the
!>    stored MQC_TwoERIs block of the space block (P,Q,R,S) of an MQC NEF R4
!>    variable, null when the block is zero. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    iref avoids the block copy of iat in hot loops. The NEF_R4 argument must
!>    have the TARGET attribute; the pointer is valid only while the calling
!>    procedure runs and the NEF_R4 is not modified. Never write through it
!>    (use iput). Test isBlockZero first: a zero block gives a null pointer.
!>
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_eris_ref(NEFR4,P,Q,R,S) result(ERIsPtr)
!
      implicit none
      class(mqc_nef_r4),intent(in),target::NEFR4
      integer,intent(in)::P,Q,R,S
      type(mqc_twoeris),pointer::ERIsPtr
!
      ERIsPtr => mqc_nef_r4_eris_ref_pairs(NEFR4,[1,P],[2,Q],[3,R],[4,S])
!
      end function mqc_nef_r4_eris_ref
!
      function mqc_nef_r4_eris_ref_pairs(NEFR4,set1,set2,set3,set4) result(ERIsPtr)
!
      implicit none
      class(mqc_nef_r4),intent(in),target::NEFR4
      integer(kind=int64),dimension(2),intent(in)::set1,set2,set3,set4
      type(mqc_twoeris),pointer::ERIsPtr
      integer(kind=int64)::nef_bidx
      logical::isMirror
!
      if(.not.allocated(NEFR4%Zero_Flag)) call mqc_error(&
           'mqc_nef_r4_eris_ref: NEF_R4 is open (never initialized)',6)
      call mqc_nef_r4_resolve(NEFR4,set1,set2,set3,set4,nef_bidx,isMirror)
!     A mirror slot holds no block of its own (its block is the dagger of the
!     stored one): read it with iat.
      if(isMirror) call mqc_error('mqc_nef_r4_eris_ref: block is a mirror slot of a '//&
        'symmetric NEF_R4, use iat (or test isMirrorSlot)',6)
      ERIsPtr => null()
      if(.not.NEFR4%Zero_Flag(nef_bidx)) ERIsPtr => NEFR4%Main_R4(nef_bidx)%R4Slot
!
      end function mqc_nef_r4_eris_ref_pairs
!
!
!     PROCEDURE MQC_NEF_R4_Is_Block_Zero
!
!>    \brief <b> MQC_NEF_R4_Is_Block_Zero (isBlockZero) reports whether the space
!>    block (P,Q,R,S) of an MQC NEF R4 variable is zero (not stored). </b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_is_block_zero(NEFR4,P,Q,R,S) result(isZero)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      integer,intent(in)::P,Q,R,S
      logical::isZero
!
      isZero = mqc_nef_r4_is_block_zero_pairs(NEFR4,[1,P],[2,Q],[3,R],[4,S])
!
      end function mqc_nef_r4_is_block_zero
!
      function mqc_nef_r4_is_block_zero_pairs(NEFR4,set1,set2,set3,set4) result(isZero)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      integer(kind=int64),dimension(2),intent(in)::set1,set2,set3,set4
      logical::isZero
      integer(kind=int64)::nef_bidx
      logical::isMirror
!
      call mqc_nef_r4_resolve(NEFR4,set1,set2,set3,set4,nef_bidx,isMirror)
      isZero = NEFR4%Zero_Flag(nef_bidx)
!
      end function mqc_nef_r4_is_block_zero_pairs
!
!
!     PROCEDURE MQC_NEF_R4_Is_Mirror_Slot
!
!>    \brief <b> MQC_NEF_R4_Is_Mirror_Slot (isMirrorSlot) reports whether the space
!>    block (P,Q,R,S) of an MQC NEF R4 variable is a mirror slot. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Is_Mirror_Slot returns .true. iff the object symmetry holds
!>    (fixed at init) and (P,Q,R,S) is the member of its dagger pair
!>    (P,Q,R,S)/(S,R,Q,P) that is not stored. Such a block is read with iat
!>    (the +/- TwoERIs dagger of the stored block); iref raises an error.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to query.
!>    \endverbatim
!>
!>    \param[in] P, Q, R, S
!>    \verbatim
!>        P, Q, R, S are Integer
!>        Space indices of the block on axes 1-4.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_is_mirror_slot(NEFR4,P,Q,R,S) result(isMirror)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      integer,intent(in)::P,Q,R,S
      logical::isMirror
      integer(kind=int64)::nef_bidx
!
      call mqc_nef_r4_resolve(NEFR4,[1,P],[2,Q],[3,R],[4,S],nef_bidx,isMirror)
!
      end function mqc_nef_r4_is_mirror_slot
!
!
!     PROCEDURE MQC_NEF_R4_Resolve
!
!>    \brief <b> MQC_NEF_R4_Resolve (private) returns the stored slot of a space
!>    block of an MQC NEF R4 variable. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Resolve converts the four [axis,index] pairs to the block
!>    (P,Q,R,S) and returns the slot that holds it. Without object symmetry
!>    this is the block's own slot. With object symmetry (Mirrored, fixed at
!>    init) the pair (P,Q,R,S)/(S,R,Q,P) shares the slot with the smaller
!>    block index; isMirror is .true. when (P,Q,R,S) is the other member.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>    \endverbatim
!>
!>    \param[in] set1, set2, set3, set4
!>    \verbatim
!>        set1-4 are Integer(kind=int64),dimension(2)
!>        [axis,space index] pairs, one per axis.
!>    \endverbatim
!>
!>    \param[out] nef_bidx
!>    \verbatim
!>        nef_bidx is Integer(kind=int64)
!>        Index of the slot that holds the block.
!>    \endverbatim
!>
!>    \param[out] isMirror
!>    \verbatim
!>        isMirror is Logical
!>        .true. iff the block is the (+/-) dagger of the stored one.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_r4_resolve(NEFR4,set1,set2,set3,set4,nef_bidx,isMirror)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      integer(kind=int64),dimension(2),intent(in)::set1,set2,set3,set4
      integer(kind=int64),intent(out)::nef_bidx
      logical,intent(out)::isMirror
      integer(kind=int64),dimension(4)::blk
      integer(kind=int64)::mirror_bidx
!
      nef_bidx = NEFR4%blockIndex(set1,set2,set3,set4)
      isMirror = .false.
      if(NEFR4%Mirrored) then
        blk(set1(1)) = set1(2)
        blk(set2(1)) = set2(2)
        blk(set3(1)) = set3(2)
        blk(set4(1)) = set4(2)
        mirror_bidx = NEFR4%blockIndex([1_int64,blk(4)],[2_int64,blk(3)],&
                                       [3_int64,blk(2)],[4_int64,blk(1)])
        if(mirror_bidx.lt.nef_bidx) then
          nef_bidx = mirror_bidx
          isMirror = .true.
        end if
      end if
!
      end subroutine mqc_nef_r4_resolve
!
!
!     PROCEDURE MQC_NEF_Maps_Equal
!
!>    \brief <b> MQC_NEF_Maps_Equal (private) tests whether two NEF space-map
!>    arrays are equal. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Maps_Equal returns .true. iff both arrays have the same number
!>    of spaces and every pair of maps is equal (MQC_SCF_Eigenvalues .eq.).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] mapA, mapB
!>    \verbatim
!>        mapA, mapB are Type(MQC_SCF_Eigenvalues),Dimension(:)
!>        The space maps to compare.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_maps_equal(mapA,mapB) result(same)
!
      implicit none
      type(mqc_scf_eigenvalues),dimension(:),intent(in)::mapA,mapB
      logical::same
      integer::i
!
      same = (size(mapA).eq.size(mapB))
      if(same) then
        do i = 1, size(mapA)
          if(mapA(i).ne.mapB(i)) then
            same = .false.
            exit
          end if
        end do
      end if
!
      end function mqc_nef_maps_equal
!
!
!     PROCEDURE MQC_NEF_Pair_Order_Label
!
!>    \brief <b> MQC_NEF_Pair_Order_Label (private) returns the operator label of an
!>    R4 <-> NEF_Matrix collapse/expand. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Pair_Order_Label applies the collapse/expand rule (user,
!>    2026-11-09): for R(p,q,r,s) = <I|a+_p a_q O a+_r a_s|J> the collapsed
!>    element ((i,j),(k,l)) takes axis Order(1..4) = i,j,k,l; it is the
!>    dagger-consistent image only when the row pair is the bra excitation
!>    and the column pair the ket excitation, i.e. Order(3) = 5-Order(1) and
!>    Order(4) = 5-Order(2) (e.g. [2,1,3,4]), with the same spin type on bra
!>    and ket (checked for singlet, triplet and all collapses; the triplet
!>    and all column pair was fixed 2026-11-09). Defaults as in collapse/expand:
!>    Order [1,2,3,4]; spinTypeKet = spinType; both 'all' when absent.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] label
!>    \verbatim
!>        label is Character(len=*)
!>        Operator label of the input.
!>    \endverbatim
!>
!>    \param[in] spinType, spinTypeKet
!>    \verbatim
!>        spinType, spinTypeKet are Character(len=*), Optional
!>        Bra and ket spin types of the collapse/expand.
!>    \endverbatim
!>
!>    \param[in] Order
!>    \verbatim
!>        Order is Integer(kind=int64),Dimension(4), Optional
!>        Axis order of the collapse/expand.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_pair_order_label(label,spinType,spinTypeKet,Order) result(outSymm)
!
      implicit none
      character(len=*),intent(in)::label
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      character(len=16)::outSymm
      character(len=64)::mySpinType,mySpinTypeKet
      integer(kind=int64),dimension(4)::myOrder
!
      mySpinType = 'all'
      mySpinTypeKet = 'all'
      if(present(spinType)) then
        mySpinType = spinType
        mySpinTypeKet = spinType
        if(present(spinTypeKet)) mySpinTypeKet = spinTypeKet
      else if(present(spinTypeKet)) then
        mySpinType = spinTypeKet
        mySpinTypeKet = spinTypeKet
      end if
      myOrder = [1,2,3,4]
      if(present(Order)) myOrder = Order
      outSymm = 'asymmetric'
      if(myOrder(3).eq.5-myOrder(1).and.myOrder(4).eq.5-myOrder(2).and.&
         mySpinType.eq.mySpinTypeKet) outSymm = label
!
      end function mqc_nef_pair_order_label
!
!
!     PROCEDURE MQC_NEF_Matrix_Equal
!
!>    \brief <b> MQC_NEF_Matrix_Equal (private) tests whether two MQC NEF Matrix
!>    variables hold the same content. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Equal returns .true. iff both matrices have equal row and
!>    column maps and the Frobenius norm of their difference is at most
!>    1e-12 relative to the norm of the first. Used to decide whether the
!>    dagger-paired axes of an orbital transform use paired coefficients.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA, NEFMatB
!>    \verbatim
!>        NEFMatA, NEFMatB are Type(MQC_NEF_Matrix)
!>        The matrices to compare.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_equal(NEFMatA,NEFMatB) result(same)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA,NEFMatB
      logical::same
!
      same = mqc_nef_maps_equal(NEFMatA%Map_Matrix_Row,NEFMatB%Map_Matrix_Row).and.&
             mqc_nef_maps_equal(NEFMatA%Map_Matrix_Col,NEFMatB%Map_Matrix_Col)
      if(same) same = (mqc_nef_matrix_frobenius(NEFMatA-NEFMatB).le.&
                       1.0e-12*(1.0+mqc_nef_matrix_frobenius(NEFMatA)))
!
      end function mqc_nef_matrix_equal
!
!
!     PROCEDURE MQC_NEF_Plain_Matrix_Label
!
!>    \brief <b> MQC_NEF_Plain_Matrix_Label (private) returns the operator label of
!>    a plain MQC matrix used as the middle factor of a hermitian sandwich. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Plain_Matrix_Label returns 'hermitian' or 'antihermitian' when the
!>    matrix passes the mathematical test (MQC_Matrix_Test_Symmetric with
!>    math=.true.; real symmetric counts as hermitian), and 'asymmetric'
!>    otherwise.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] Mat
!>    \verbatim
!>        Mat is Type(MQC_Matrix)
!>        The matrix to classify.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_plain_matrix_label(Mat) result(label)
!
      implicit none
      type(mqc_matrix),intent(in)::Mat
      character(len=16)::label
!
      label = 'asymmetric'
      if(MQC_Matrix_Test_Symmetric(Mat,'hermitian',.true.)) then
        label = 'hermitian'
      else if(MQC_Matrix_Test_Symmetric(Mat,'antihermitian',.true.)) then
        label = 'antihermitian'
      end if
!
      end function mqc_nef_plain_matrix_label
!
!
!     PROCEDURE MQC_Matrix_NEFSet_Herm_Sandwich
!
!>    \brief <b> MQC_Matrix_NEFSet_Herm_Sandwich (HermMatMul) computes the grid-level
!>    hermitian sandwich B^dagger A B of a NEF_Set B and a plain matrix A. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Matrix_NEFSet_Herm_Sandwich computes, for a NEF_Set B with Vector or
!>    Matrix cells (grid NumI x NumJ) and a plain NumI x NumI matrix A,
!>      Out(I,K) = Sum_{J,L} A(J,L) outer( op(B(J,I)), B(L,K) )
!>    with op = conj (Vector cells) or dagger (Matrix cells); Vector cells give
!>    NEF_Matrix cells and Matrix cells give NEF_R4 cells (as in the NEF_Set
!>    product). The output grid is NumJ x NumJ over B's column list
!>    (sameBraKetSpace = .true.). Its label is the label of A (mathematical
!>    test): B^dagger A B is hermitian / antihermitian by construction, so the
!>    output is then grid-symmetric and only the I >= K cells are computed.
!>    The double sum is factored, T(J,K) = Sum_L A(J,L) B(L,K) (scalar x cell
!>    sums) and Out(I,K) = Sum_J outer(op(B(J,I)), T(J,K)), so each output
!>    cell takes NumI outer products instead of NumI**2; zero cells and zero
!>    A entries are skipped. The set builds its cells with its own flags.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] MatA
!>    \verbatim
!>        MatA is Type(MQC_Matrix)
!>        The middle factor A (NumI x NumI, over B's row list).
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        The outer factor B (Vector or Matrix cells).
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_matrix_nefset_herm_sandwich(MatA,NEFSetB) result(NEFSetOut)
!
      implicit none
      type(mqc_matrix),intent(in)::MatA
      type(mqc_nef_set),intent(in)::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      type(mqc_nef_vector),dimension(:,:),allocatable::TVec
      type(mqc_nef_matrix),dimension(:,:),allocatable::TMat
      logical,dimension(:,:),allocatable::TZero
      type(mqc_nef_matrix)::cellMat
      type(mqc_nef_r4)::cellR4
      type(mqc_scalar)::aJL
      character(len=16)::outSymm
      integer::nDet,nOut,I,J,K,L
      logical::first,halfGrid
!
      nDet = NEFSetB%NumI
      nOut = NEFSetB%NumJ
      if(MQC_Matrix_Rows(MatA).ne.nDet.or.MQC_Matrix_Columns(MatA).ne.nDet) &
        call mqc_error_i('HermMatMul(matrix,NEF_Set): A must be NumI x NumI of the set',6,&
          'rows(A)',MQC_Matrix_Rows(MatA),'columns(A)',MQC_Matrix_Columns(MatA),&
          'NumI',int(nDet,kind=int64))
      outSymm = mqc_nef_plain_matrix_label(MatA)
      halfGrid = (outSymm.ne.'asymmetric')
      call NEFSetOut%init(nOut,nOut,sameBraKetSpace=.true.,symmetry=outSymm)
      allocate(TZero(nDet,nOut))
      TZero = .true.
      select case(trim(NEFSetB%Type()))
      case('vector')
        allocate(TVec(nDet,nOut))
        do K = 1, nOut
          do L = 1, nDet
            if(NEFSetB%isZero(L,K)) cycle
            do J = 1, nDet
              aJL = MatA%at(J,L)
              if(abs(aJL).eq.0.0) cycle
              if(TZero(J,K)) then
                TVec(J,K) = aJL*NEFSetB%atVector(L,K)
                TZero(J,K) = .false.
              else
                TVec(J,K) = TVec(J,K) + aJL*NEFSetB%atVector(L,K)
              end if
            end do
          end do
        end do
        do I = 1, nOut
          do K = 1, nOut
            if(halfGrid.and.I.lt.K) cycle
            first = .true.
            do J = 1, nDet
              if(NEFSetB%isZero(J,I).or.TZero(J,K)) cycle
              if(first) then
                cellMat = mqc_nefvector_nefvector_outer_product(conjg(NEFSetB%atVector(J,I)),TVec(J,K))
                first = .false.
              else
                cellMat = cellMat + mqc_nefvector_nefvector_outer_product(conjg(NEFSetB%atVector(J,I)),&
                          TVec(J,K))
              end if
            end do
            if(.not.first) call mqc_nef_set_put_built_matrix(NEFSetOut,cellMat,I,K)
          end do
        end do
      case('matrix')
        allocate(TMat(nDet,nOut))
        do K = 1, nOut
          do L = 1, nDet
            if(NEFSetB%isZero(L,K)) cycle
            do J = 1, nDet
              aJL = MatA%at(J,L)
              if(abs(aJL).eq.0.0) cycle
              if(TZero(J,K)) then
                TMat(J,K) = aJL*NEFSetB%atMatrix(L,K)
                TZero(J,K) = .false.
              else
                TMat(J,K) = TMat(J,K) + aJL*NEFSetB%atMatrix(L,K)
              end if
            end do
          end do
        end do
        do I = 1, nOut
          do K = 1, nOut
            if(halfGrid.and.I.lt.K) cycle
            first = .true.
            do J = 1, nDet
              if(NEFSetB%isZero(J,I).or.TZero(J,K)) cycle
              if(first) then
                cellR4 = mqc_nefmatrix_nefmatrix_outer_product(dagger(NEFSetB%atMatrix(J,I)),TMat(J,K))
                first = .false.
              else
                cellR4 = cellR4 + mqc_nefmatrix_nefmatrix_outer_product(dagger(NEFSetB%atMatrix(J,I)),&
                         TMat(J,K))
              end if
            end do
            if(.not.first) call mqc_nef_set_put_built_r4(NEFSetOut,cellR4,I,K)
          end do
        end do
      case default
        call mqc_error_A('HermMatMul(matrix,NEF_Set): the set must hold vector or matrix '//&
             'cells (R4 cells have no rank-8 target)',6,'Array_Type',NEFSetB%Type())
      end select
!
      end function mqc_matrix_nefset_herm_sandwich
!
!
!     PROCEDURE MQC_NEFSet_Matrix_Herm_Sandwich
!
!>    \brief <b> MQC_NEFSet_Matrix_Herm_Sandwich (HermMatMul) computes the grid-level
!>    hermitian sandwich C^dagger S C of a NEF_Set S and a plain coefficient
!>    matrix C. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFSet_Matrix_Herm_Sandwich computes, for a square NEF_Set S (N x N,
!>    any cell kind) and a plain N x M matrix C,
!>      Out(a,b) = Sum_{I,J} conj(C(I,a)) S(I,J) C(J,b)
!>    (cells are linear combinations of S's cells, same kind). The output grid
!>    is M x M (sameBraKetSpace = .true.). If S is grid-symmetric (hermitian or
!>    antihermitian over one determinant list), the output is too and keeps
!>    S's label; only the a >= b cells are computed. Factored as
!>    T(I,b) = Sum_J S(I,J) C(J,b), Out(a,b) = Sum_I conj(C(I,a)) T(I,b); zero
!>    cells and zero coefficients are skipped. The set builds its cells with
!>    its own flags.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSetA
!>    \verbatim
!>        NEFSetA is Type(MQC_NEF_Set)
!>        The middle factor S (square grid).
!>    \endverbatim
!>
!>    \param[in] MatB
!>    \verbatim
!>        MatB is Type(MQC_Matrix)
!>        The outer factor C (N x M).
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_matrix_herm_sandwich(NEFSetA,MatB) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSetA
      type(mqc_matrix),intent(in)::MatB
      type(mqc_nef_set)::NEFSetOut
      type(mqc_nef_vector),dimension(:,:),allocatable::TVec
      type(mqc_nef_matrix),dimension(:,:),allocatable::TMat
      type(mqc_nef_r4),dimension(:,:),allocatable::TR4
      logical,dimension(:,:),allocatable::TZero
      type(mqc_nef_vector)::cellVec
      type(mqc_nef_matrix)::cellMat
      type(mqc_nef_r4)::cellR4
      type(mqc_scalar)::cJb,cIa
      character(len=16)::outSymm
      character(len=64)::kindS
      integer::N,M,I,J,a,b
      logical::first,halfGrid
!
      N = NEFSetA%NumI
      M = MQC_Matrix_Columns(MatB)
      if(NEFSetA%NumJ.ne.N.or.MQC_Matrix_Rows(MatB).ne.N) &
        call mqc_error_i('HermMatMul(NEF_Set,matrix): S must be square and C must have '//&
          'NumI rows',6,'NumI',int(N,kind=int64),'NumJ',int(NEFSetA%NumJ,kind=int64),&
          'rows(C)',MQC_Matrix_Rows(MatB))
      outSymm = 'asymmetric'
      if(NEFSetA%GridMirrored) outSymm = NEFSetA%Symmetry
      halfGrid = (outSymm.ne.'asymmetric')
      call NEFSetOut%init(M,M,sameBraKetSpace=.true.,symmetry=outSymm)
      kindS = NEFSetA%Type()
      allocate(TZero(N,M))
      TZero = .true.
      select case(trim(kindS))
      case('vector')
        allocate(TVec(N,M))
      case('matrix')
        allocate(TMat(N,M))
      case('r4')
        allocate(TR4(N,M))
      case default
        call mqc_error_A('HermMatMul(NEF_Set,matrix): NEF_Set is open (holds no elements)',&
             6,'Array_Type',NEFSetA%Type())
      end select
!     T(I,b) = Sum_J S(I,J) C(J,b)
      do b = 1, M
        do J = 1, N
          cJb = MatB%at(J,b)
          if(abs(cJb).eq.0.0) cycle
          do I = 1, N
            if(NEFSetA%isZero(I,J)) cycle
            select case(trim(kindS))
            case('vector')
              if(TZero(I,b)) then
                TVec(I,b) = cJb*NEFSetA%atVector(I,J)
              else
                TVec(I,b) = TVec(I,b) + cJb*NEFSetA%atVector(I,J)
              end if
            case('matrix')
              if(TZero(I,b)) then
                TMat(I,b) = cJb*NEFSetA%atMatrix(I,J)
              else
                TMat(I,b) = TMat(I,b) + cJb*NEFSetA%atMatrix(I,J)
              end if
            case('r4')
              if(TZero(I,b)) then
                TR4(I,b) = cJb*NEFSetA%atR4(I,J)
              else
                TR4(I,b) = TR4(I,b) + cJb*NEFSetA%atR4(I,J)
              end if
            end select
            TZero(I,b) = .false.
          end do
        end do
      end do
!     Out(a,b) = Sum_I conj(C(I,a)) T(I,b)
      do a = 1, M
        do b = 1, M
          if(halfGrid.and.a.lt.b) cycle
          first = .true.
          do I = 1, N
            if(TZero(I,b)) cycle
            cIa = conjg(MatB%at(I,a))
            if(abs(cIa).eq.0.0) cycle
            select case(trim(kindS))
            case('vector')
              if(first) then
                cellVec = cIa*TVec(I,b)
              else
                cellVec = cellVec + cIa*TVec(I,b)
              end if
            case('matrix')
              if(first) then
                cellMat = cIa*TMat(I,b)
              else
                cellMat = cellMat + cIa*TMat(I,b)
              end if
            case('r4')
              if(first) then
                cellR4 = cIa*TR4(I,b)
              else
                cellR4 = cellR4 + cIa*TR4(I,b)
              end if
            end select
            first = .false.
          end do
          if(first) cycle
          select case(trim(kindS))
          case('vector')
            call mqc_nef_set_put_built_vector(NEFSetOut,cellVec,a,b)
          case('matrix')
            call mqc_nef_set_put_built_matrix(NEFSetOut,cellMat,a,b)
          case('r4')
            call mqc_nef_set_put_built_r4(NEFSetOut,cellR4,a,b)
          end select
        end do
      end do
!
      end function mqc_nefset_matrix_herm_sandwich
!
!
!     PROCEDURE MQC_NEFMatrix_Herm_Sandwich
!
!>    \brief <b> MQC_NEFMatrix_Herm_Sandwich (HermMatMul) computes the hermitian
!>    sandwich B^dagger A B of two MQC NEF Matrix variables. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_Herm_Sandwich computes Out = B^dagger A B (e.g. an operator
!>    transformed to another orbital basis). The output keeps A's operator
!>    label (B^dagger A B of a (anti)hermitian A is (anti)hermitian by
!>    construction) and A's sameBraKet flag; its row and column maps are B's
!>    column maps. T = A B is formed once (NEF matmul); with object symmetry
!>    only the lower blocks Out(r,c) = Sum_p B(p,r)^dagger T(p,c), r >= c, are
!>    computed (the rest is the mirror).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is Type(MQC_NEF_Matrix)
!>        The middle factor A.
!>    \endverbatim
!>
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is Type(MQC_NEF_Matrix)
!>        The outer factor B (rows over A's spaces).
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefmatrix_herm_sandwich(NEFMatA,NEFMatB) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMatA,NEFMatB
      type(mqc_nef_matrix)::NEFMatOut,T
      type(mqc_scf_eigenvalues),dimension(:),allocatable::outMap
      type(mqc_scf_integral)::blk
      integer::r,c,p
      logical::first
!
      T = MatMul(NEFMatA,NEFMatB)
      outMap = NEFMatB%Map_Matrix_Col
      call NEFMatOut%init(outMap,initial_map_matrix_col=outMap,&
           symmetry=NEFMatA%Symmetry,sameBraKet=NEFMatA%sameBraKet)
      do r = 1, NEFMatOut%NumSpace_row
        do c = 1, NEFMatOut%NumSpace_col
          if(NEFMatOut%Mirrored.and.c.gt.r) cycle
          first = .true.
          do p = 1, NEFMatB%NumSpace_row
            if(NEFMatB%isBlockZero(p,r).or.T%isBlockZero(p,c)) cycle
            if(first) then
              blk = matmul(dagger(NEFMatB%iat(p,r)),T%iat(p,c))
              first = .false.
            else
              blk = blk + matmul(dagger(NEFMatB%iat(p,r)),T%iat(p,c))
            end if
          end do
          if(.not.first) call NEFMatOut%iput(blk,r,c)
        end do
      end do
!
      end function mqc_nefmatrix_herm_sandwich
!
!
!     PROCEDURE MQC_Integral_NEFMatrix_Herm_Sandwich
!
!>    \brief <b> MQC_Integral_NEFMatrix_Herm_Sandwich (HermMatMul) computes the
!>    hermitian sandwich B^dagger A B of an MQC SCF Integral A and an MQC NEF
!>    Matrix B. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Integral_NEFMatrix_Herm_Sandwich computes Out = B^dagger A B, e.g. an
!>    AO operator (core Hamiltonian, overlap, dipole, CAP) transformed to the
!>    NEF orbital spaces of B. B must have one row space spanning A's basis
!>    (same conformability as MQC_Integral_NEFMatrix_Multiply). Row and column
!>    maps of Out are B's column maps; sameBraKet is .true. (the same B on both
!>    sides). The operator label comes from A's integral symmetry flag:
!>      'hermitian', 'realdiag'      -> 'hermitian'
!>      'antihermitian', 'imagdiag'  -> 'antihermitian'
!>      anything else                -> 'asymmetric'
!>    With object symmetry (Mirrored) only the blocks r >= c are computed.
!>    Diagonal blocks use the integral HermMatMul (mqc_est, packed output);
!>    off-diagonal blocks are dagger(B(1,r)) T_c with T_c = A B(1,c) formed once
!>    per column space.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] IntegralA
!>    \verbatim
!>        IntegralA is Type(MQC_SCF_Integral)
!>        The middle factor A (square, in the basis of B's row space).
!>    \endverbatim
!>
!>    \param[in] NEFMatB
!>    \verbatim
!>        NEFMatB is Type(MQC_NEF_Matrix)
!>        The outer factor B (one row space).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_integral_nefmatrix_herm_sandwich(IntegralA,NEFMatB) result(NEFMatOut)
!
      implicit none
      type(mqc_scf_integral),intent(in)::IntegralA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_matrix)::NEFMatOut
      type(mqc_scf_eigenvalues),dimension(:),allocatable::outMap
      type(mqc_scf_integral)::T
      character(len=16)::outSymm
      integer::r,c
      logical::haveT
!
      if(NEFMatB%NumSpace_row.ne.1) &
        call mqc_error_i('HermMatMul(integral,NEF_Matrix): B must have one row space',6,&
        'NEFMatB%NumSpace_row',NEFMatB%NumSpace_row)
      if(size(IntegralA,'alpha',1).ne.size(IntegralA,'alpha',2).or.&
         size(IntegralA,'beta',1).ne.size(IntegralA,'beta',2)) &
        call mqc_error_i('HermMatMul(integral,NEF_Matrix): A must be square',6,&
        'size(IntegralA,alpha,1)',size(IntegralA,'alpha',1),&
        'size(IntegralA,alpha,2)',size(IntegralA,'alpha',2),&
        'size(IntegralA,beta,1)',size(IntegralA,'beta',1),&
        'size(IntegralA,beta,2)',size(IntegralA,'beta',2))
      if(NEFMatB%Size_vector_alpha_row(1).ne.size(IntegralA,'alpha',2).or.&
         NEFMatB%Size_vector_beta_row(1).ne.size(IntegralA,'beta',2)) &
        call mqc_error_i('HermMatMul(integral,NEF_Matrix): B rows not conformable with A',6,&
        'NEFMatB%Size_vector_alpha_row',NEFMatB%Size_vector_alpha_row(1),&
        'size(IntegralA,alpha,2)',size(IntegralA,'alpha',2),&
        'NEFMatB%Size_vector_beta_row',NEFMatB%Size_vector_beta_row(1),&
        'size(IntegralA,beta,2)',size(IntegralA,'beta',2))
!
      select case(trim(IntegralA%getSymm()))
      case('hermitian','realdiag')
        outSymm = 'hermitian'
      case('antihermitian','imagdiag')
        outSymm = 'antihermitian'
      case default
        outSymm = 'asymmetric'
      end select
!
      outMap = NEFMatB%Map_Matrix_Col
      call NEFMatOut%init(outMap,initial_map_matrix_col=outMap,&
           symmetry=outSymm,sameBraKet=.true.)
      do c = 1, NEFMatOut%NumSpace_col
        if(NEFMatB%isBlockZero(1,c)) cycle
        haveT = .false.
        do r = 1, NEFMatOut%NumSpace_row
          if(NEFMatOut%Mirrored.and.c.gt.r) cycle
          if(NEFMatB%isBlockZero(1,r)) cycle
          if(r.eq.c) then
            call NEFMatOut%iput(HermMatMul(IntegralA,NEFMatB%iat(1,c)),r,c)
          else
            if(.not.haveT) then
              T = matmul(IntegralA,NEFMatB%iat(1,c))
              haveT = .true.
            end if
            call NEFMatOut%iput(matmul(dagger(NEFMatB%iat(1,r)),T),r,c)
          end if
        end do
      end do
!
      end function mqc_integral_nefmatrix_herm_sandwich
!
!
!     PROCEDURE MQC_NEF_Vector_Is_Same_BraKet
!
!>    \brief <b> MQC_NEF_Vector_Is_Same_BraKet (isSameBraKet) returns the
!>    sameBraKet flag of an MQC NEF Vector variable. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Is_Same_BraKet returns .true. iff the object is marked
!>    as <I|O|I> (given at init). Bra-ket swapping
!>    symmetries of the object (e.g. dagger-paired blocks) may only be used
!>    when this flag is .true..
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVec
!>    \verbatim
!>        NEFVec is Class(MQC_NEF_Vector)
!>        The MQC NEF Vector to query.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_is_same_braket(NEFVec) result(sameBraKet)
!
      implicit none
      class(mqc_nef_vector),intent(in)::NEFVec
      logical::sameBraKet
!
      sameBraKet = NEFVec%sameBraKet
!
      end function mqc_nef_vector_is_same_braket
!
!
!     PROCEDURE MQC_NEF_Vector_Get_Symmetry
!
!>    \brief <b> MQC_NEF_Vector_Get_Symmetry (getSymm) returns the operator
!>    symmetry label of an MQC NEF Vector variable. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Vector_Get_Symmetry returns the operator label given at init:
!>    'hermitian', 'antihermitian' or 'asymmetric'. The label is a property
!>    of the operator O; NEF never derives it.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFVec
!>    \verbatim
!>        NEFVec is Class(MQC_NEF_Vector)
!>        The MQC NEF Vector to query.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_get_symmetry(NEFVec) result(symm)
!
      implicit none
      class(mqc_nef_vector),intent(in)::NEFVec
      character(len=16)::symm
!
      symm = NEFVec%Symmetry
!
      end function mqc_nef_vector_get_symmetry
!
!
!     PROCEDURE MQC_NEF_Matrix_Is_Same_BraKet
!
!>    \brief <b> MQC_NEF_Matrix_Is_Same_BraKet (isSameBraKet) returns the
!>    sameBraKet flag of an MQC NEF Matrix variable. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Matrix_Is_Same_BraKet returns .true. iff the object is marked
!>    as <I|O|I> (given at init). Bra-ket swapping
!>    symmetries of the object (e.g. dagger-paired blocks) may only be used
!>    when this flag is .true..
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
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_is_same_braket(NEFMat) result(sameBraKet)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      logical::sameBraKet
!
      sameBraKet = NEFMat%sameBraKet
!
      end function mqc_nef_matrix_is_same_braket
!
!
!     PROCEDURE MQC_NEF_R4_Is_Same_BraKet
!
!>    \brief <b> MQC_NEF_R4_Is_Same_BraKet (isSameBraKet) returns the
!>    sameBraKet flag of an MQC NEF R4 variable. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Is_Same_BraKet returns .true. iff the object is marked
!>    as <I|O|I> (given at init). Bra-ket swapping
!>    symmetries of the object (e.g. dagger-paired blocks) may only be used
!>    when this flag is .true..
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to query.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_is_same_braket(NEFR4) result(sameBraKet)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      logical::sameBraKet
!
      sameBraKet = NEFR4%sameBraKet
!
      end function mqc_nef_r4_is_same_braket
!
!
!     PROCEDURE MQC_NEF_R4_Get_Symmetry
!
!>    \brief <b> MQC_NEF_R4_Get_Symmetry (getSymm) returns the operator
!>    symmetry label of an MQC NEF R4 variable. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Get_Symmetry returns the operator label given at init:
!>    'hermitian', 'antihermitian' or 'asymmetric'. The label is a property
!>    of the operator O; NEF never derives it.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to query.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_get_symmetry(NEFR4) result(symm)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      character(len=16)::symm
!
      symm = NEFR4%Symmetry
!
      end function mqc_nef_r4_get_symmetry
!
!
!     PROCEDURE MQC_NEF_R4_Allocate_Slots
!
!>    \brief <b> MQC_NEF_R4_Allocate_Slots (private) (re)allocates all block slots of
!>    an MQC NEF R4 variable empty (all blocks zero) from NumSpace_1..4. </b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_r4_allocate_slots(NEFR4)
!
      implicit none
      type(mqc_nef_r4),intent(inOut)::NEFR4
!
      if(allocated(NEFR4%Main_R4)) deallocate(NEFR4%Main_R4)
      if(allocated(NEFR4%Zero_Flag)) deallocate(NEFR4%Zero_Flag)
      allocate(NEFR4%Main_R4(NEFR4%NumSpace_1*NEFR4%NumSpace_2*NEFR4%NumSpace_3*NEFR4%NumSpace_4))
      allocate(NEFR4%Zero_Flag(size(NEFR4%Main_R4)))
      NEFR4%Zero_Flag = .true.
!
      end subroutine mqc_nef_r4_allocate_slots
!
!
!     PROCEDURE MQC_NEF_R4_Initialize
!
!>    \brief <b> MQC_NEF_R4_Initialize (init) initializes an MQC NEF R4 variable
!>    with all space blocks zero. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Initialize records the four axis maps of an MQC NEF R4
!>    (axis k runs over the spaces of map_matrix_k), sets the alpha/beta
!>    block sizes of every space from its map, and allocates the empty
!>    space-block slots (all blocks zero). The operator label and the
!>    sameBraKet flag are given here by the builder and are never changed
!>    afterwards; defaults are 'asymmetric' and .false.. When sameBraKet is
!>    .true., the label is not 'asymmetric' and the maps of axes 1/4 and 2/3
!>    are equal, the object symmetry holds: only one block of each dagger
!>    pair (P,Q,R,S)/(S,R,Q,P) is stored (the one with the smaller block
!>    index) and put/iput/at/iat resolve the other through the TwoERIs dagger.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to initialize.
!>    \endverbatim
!>
!>    \param[in] map_matrix_1
!>    \verbatim
!>        map_matrix_1 is Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable
!>        Space maps of axis 1.
!>    \endverbatim
!>
!>    \param[in] map_matrix_2
!>    \verbatim
!>        map_matrix_2 is Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable
!>        Space maps of axis 2.
!>    \endverbatim
!>
!>    \param[in] map_matrix_3
!>    \verbatim
!>        map_matrix_3 is Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable
!>        Space maps of axis 3.
!>    \endverbatim
!>
!>    \param[in] map_matrix_4
!>    \verbatim
!>        map_matrix_4 is Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable
!>        Space maps of axis 4.
!>    \endverbatim
!>
!>    \param[in] symmetry
!>    \verbatim
!>        symmetry is Character(len=*), Optional
!>        Operator label: 'hermitian', 'antihermitian' or 'asymmetric'
!>        (default).
!>    \endverbatim
!>
!>    \param[in] sameBraKet
!>    \verbatim
!>        sameBraKet is Logical, Optional
!>        .true. iff the object is <I|O|I>. Default .false..
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_r4_initialize(NEFR4,map_matrix_1,map_matrix_2,map_matrix_3,&
           map_matrix_4,symmetry,sameBraKet)
!
      implicit none
      class(mqc_nef_r4),intent(inOut)::NEFR4
      type(mqc_scf_eigenvalues),dimension(:),allocatable,intent(in)::map_matrix_1,&
        map_matrix_2,map_matrix_3,map_matrix_4
      character(len=*),intent(in),optional::symmetry
      logical,intent(in),optional::sameBraKet
      integer::i
!
      NEFR4%NumSpace_1 = size(map_matrix_1)
      NEFR4%NumSpace_2 = size(map_matrix_2)
      NEFR4%NumSpace_3 = size(map_matrix_3)
      NEFR4%NumSpace_4 = size(map_matrix_4)
      NEFR4%Map_Matrix_1 = map_matrix_1
      NEFR4%Map_Matrix_2 = map_matrix_2
      NEFR4%Map_Matrix_3 = map_matrix_3
      NEFR4%Map_Matrix_4 = map_matrix_4
      NEFR4%Size_vector_alpha_1 = [(Size(map_matrix_1(i),'alpha'), i = 1, NEFR4%NumSpace_1)]
      NEFR4%Size_vector_beta_1  = [(Size(map_matrix_1(i),'beta'),  i = 1, NEFR4%NumSpace_1)]
      NEFR4%Size_vector_alpha_2 = [(Size(map_matrix_2(i),'alpha'), i = 1, NEFR4%NumSpace_2)]
      NEFR4%Size_vector_beta_2  = [(Size(map_matrix_2(i),'beta'),  i = 1, NEFR4%NumSpace_2)]
      NEFR4%Size_vector_alpha_3 = [(Size(map_matrix_3(i),'alpha'), i = 1, NEFR4%NumSpace_3)]
      NEFR4%Size_vector_beta_3  = [(Size(map_matrix_3(i),'beta'),  i = 1, NEFR4%NumSpace_3)]
      NEFR4%Size_vector_alpha_4 = [(Size(map_matrix_4(i),'alpha'), i = 1, NEFR4%NumSpace_4)]
      NEFR4%Size_vector_beta_4  = [(Size(map_matrix_4(i),'beta'),  i = 1, NEFR4%NumSpace_4)]
      call mqc_nef_r4_allocate_slots(NEFR4)
      NEFR4%Symmetry = mqc_nef_symmetry_label(symmetry,'mqc_nef_r4_initialize')
      NEFR4%sameBraKet = .false.
      if(present(sameBraKet)) NEFR4%sameBraKet = sameBraKet
!     Object symmetry: the dagger (P,Q,R,S) -> (S,R,Q,P) exchanges axes 1<->4 and 2<->3.
      NEFR4%Mirrored = (NEFR4%sameBraKet.and.NEFR4%Symmetry.ne.'asymmetric')
      if(NEFR4%Mirrored) NEFR4%Mirrored = mqc_nef_maps_equal(map_matrix_1,map_matrix_4).and.&
                                          mqc_nef_maps_equal(map_matrix_2,map_matrix_3)
!
      end subroutine mqc_nef_r4_initialize
!
!
!     PROCEDURE MQC_NEF_R4_Full_Storage
!
!>    \brief <b> MQC_NEF_R4_Full_Storage (private) returns a copy of an MQC NEF R4
!>    variable with every space block stored. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_Full_Storage copies an MQC NEF R4 into the full slot layout:
!>    same maps, operator label and sameBraKet flag, but every block of a
!>    dagger pair is stored (mirror blocks are materialized as the +/- TwoERIs
!>    dagger of the stored block). It is an internal read copy used by
!>    routines that read blocks through iref (transitional rule I5); only the
!>    storage layout differs from the input.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to copy.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_full_storage(NEFR4) result(R4Out)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFR4
      type(mqc_nef_r4)::R4Out
      integer::id1,id2,id3,id4
!
      call R4Out%init(NEFR4%Map_Matrix_1,NEFR4%Map_Matrix_2,NEFR4%Map_Matrix_3,&
           NEFR4%Map_Matrix_4,symmetry=NEFR4%Symmetry,sameBraKet=NEFR4%sameBraKet)
!     storage layout only: every slot is its own block
      R4Out%Mirrored = .false.
      do id1 = 1, NEFR4%NumSpace_1
        do id2 = 1, NEFR4%NumSpace_2
          do id3 = 1, NEFR4%NumSpace_3
            do id4 = 1, NEFR4%NumSpace_4
              if (.not.NEFR4%isBlockZero(id1,id2,id3,id4)) &
                call R4Out%iput(NEFR4%iat(id1,id2,id3,id4),id1,id2,id3,id4)
            end do
          end do
        end do
      end do
!
      end function mqc_nef_r4_full_storage
!
!
!     PROCEDURE MQC_NEF_R4_BRef
!
!>    \brief <b> MQC_NEF_R4_BRef (private) returns a pointer to a space block of an
!>    MQC NEF R4 variable, stored or mirrored. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_R4_BRef returns a pointer to the block (P,Q,R,S) for read-only use
!>    in contraction kernels: for a stored slot the pointer is iref (no copy);
!>    for a mirror slot of a symmetric NEF_R4 the +/- TwoERIs dagger of the
!>    stored block is written into the caller's scratch block and the pointer
!>    points to it (one block, instead of a full-storage copy of the object).
!>    Zero blocks are the caller's responsibility (test isBlockZero first).
!>    The scratch block is overwritten by the next call with the same scratch.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to read.
!>    \endverbatim
!>
!>    \param[in] P, Q, R, S
!>    \verbatim
!>        P, Q, R, S are Integer
!>        Space indices of the block.
!>    \endverbatim
!>
!>    \param[in,out] scratch
!>    \verbatim
!>        scratch is Type(MQC_TwoERIs), target
!>        Holds the mirrored block when (P,Q,R,S) is a mirror slot.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_bref(NEFR4,P,Q,R,S,scratch) result(ERIsPtr)
!
      implicit none
      class(mqc_nef_r4),intent(in),target::NEFR4
      integer,intent(in)::P,Q,R,S
      type(mqc_twoeris),intent(inout),target::scratch
      type(mqc_twoeris),pointer::ERIsPtr
!
      ERIsPtr => mqc_nef_r4_bref_pairs(NEFR4,[1_int64,int(P,int64)],[2_int64,int(Q,int64)],&
                 [3_int64,int(R,int64)],[4_int64,int(S,int64)],scratch)
!
      end function mqc_nef_r4_bref
!
!
!     PROCEDURE MQC_NEF_R4_BRef_Pairs
!
!>    \brief <b> MQC_NEF_R4_BRef_Pairs (private) is MQC_NEF_R4_BRef with the block
!>    given as [axis,index] pairs. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    See MQC_NEF_R4_BRef: iref for a stored slot, the +/- dagger of the stored
!>    block in the caller's scratch for a mirror slot.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFR4 is Class(MQC_NEF_R4)
!>        The MQC NEF R4 to read.
!>    \endverbatim
!>
!>    \param[in] set1, set2, set3, set4
!>    \verbatim
!>        set1-4 are Integer(kind=int64),dimension(2)
!>        [axis,space index] pairs, one per axis.
!>    \endverbatim
!>
!>    \param[in,out] scratch
!>    \verbatim
!>        scratch is Type(MQC_TwoERIs), target
!>        Holds the mirrored block when the block is a mirror slot.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_bref_pairs(NEFR4,set1,set2,set3,set4,scratch) result(ERIsPtr)
!
      implicit none
      class(mqc_nef_r4),intent(in),target::NEFR4
      integer(kind=int64),dimension(2),intent(in)::set1,set2,set3,set4
      type(mqc_twoeris),intent(inout),target::scratch
      type(mqc_twoeris),pointer::ERIsPtr
      integer(kind=int64)::nef_bidx
      logical::isMirror
!
      call mqc_nef_r4_resolve(NEFR4,set1,set2,set3,set4,nef_bidx,isMirror)
      if(isMirror) then
        scratch = mqc_nef_r4_eris_at_pairs(NEFR4,set1,set2,set3,set4)
        ERIsPtr => scratch
      else
        ERIsPtr => mqc_nef_r4_eris_ref_pairs(NEFR4,set1,set2,set3,set4)
      end if
!
      end function mqc_nef_r4_bref_pairs
!
!
!     PROCEDURE MQC_NEF_Vector_Eigen_Put
!
!>    \brief <b> MQC_NEF_Vector_Eigen_Put (iput) stores the MQC_SCF_Eigenvalues block P
!>    of an MQC NEF Vector variable and marks it non-zero. </b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_vector_eigen_put(NEFVec,eigen,P)
!
      implicit none
      class(mqc_nef_vector),intent(inOut)::NEFVec
      type(mqc_scf_eigenvalues),intent(in)::eigen
      integer,intent(in)::P
!
      if(.not.allocated(NEFVec%Zero_Flag)) call mqc_error(&
           'mqc_nef_vector_eigen_put: NEF_Vector is open (never initialized)',6)
      NEFVec%Main_Vector(P) = eigen
      NEFVec%Zero_Flag(P) = .false.
!
      end subroutine mqc_nef_vector_eigen_put
!
!
!     PROCEDURE MQC_NEF_Vector_Eigen_At
!
!>    \brief <b> MQC_NEF_Vector_Eigen_At (iat) returns a copy of block P of an MQC NEF
!>    Vector variable; a zero block is returned as zeros of the block size. Internal
!>    loops use iref (no copy). </b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_eigen_at(NEFVec,P) result(eigenOut)
!
      implicit none
      class(mqc_nef_vector),intent(in)::NEFVec
      integer,intent(in)::P
      type(mqc_scf_eigenvalues)::eigenOut
      integer::zero
!
      if(.not.allocated(NEFVec%Zero_Flag)) call mqc_error(&
           'mqc_nef_vector_eigen_at: NEF_Vector is open (never initialized)',6)
      zero = 0
      if(NEFVec%Zero_Flag(P)) then
        if(NEFVec%Map_Matrix(P)%type().eq.'space') then
          call eigenOut%init(NEFVec%Size_vector_alpha(P),zero)
        else
          call eigenOut%init(NEFVec%Size_vector_alpha(P),zero,NEFVec%Size_vector_beta(P))
        end if
      else
        eigenOut = NEFVec%Main_Vector(P)
      end if
!
      end function mqc_nef_vector_eigen_at
!
!
!     PROCEDURE MQC_NEF_Vector_Eigen_Ref
!
!>    \brief <b> MQC_NEF_Vector_Eigen_Ref (iref) returns a read-only reference to the
!>    stored block P of an MQC NEF Vector variable, null when the block is zero. The
!>    NEF_Vector argument must have the TARGET attribute; the pointer is valid only
!>    while the calling procedure runs and the NEF_Vector is not modified. Never
!>    write through it (use iput). </b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_eigen_ref(NEFVec,P) result(eigenPtr)
!
      implicit none
      class(mqc_nef_vector),intent(in),target::NEFVec
      integer,intent(in)::P
      type(mqc_scf_eigenvalues),pointer::eigenPtr
!
      if(.not.allocated(NEFVec%Zero_Flag)) call mqc_error(&
           'mqc_nef_vector_eigen_ref: NEF_Vector is open (never initialized)',6)
      eigenPtr => null()
      if(.not.NEFVec%Zero_Flag(P)) eigenPtr => NEFVec%Main_Vector(P)
!
      end function mqc_nef_vector_eigen_ref
!
!
!     PROCEDURE MQC_NEF_Vector_Is_Block_Zero
!
!>    \brief <b> MQC_NEF_Vector_Is_Block_Zero (isBlockZero) reports whether block P of
!>    an MQC NEF Vector variable is zero (not stored). </b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_is_block_zero(NEFVec,P) result(isZero)
!
      implicit none
      class(mqc_nef_vector),intent(in)::NEFVec
      integer,intent(in)::P
      logical::isZero
!
      isZero = NEFVec%Zero_Flag(P)
!
      end function mqc_nef_vector_is_block_zero
!
!
!     PROCEDURE MQC_NEF_Matrix_Allocate_Blocks
!
!>    \brief <b> MQC_NEF_Matrix_Allocate_Blocks (private) (re)allocates all blocks of an
!>    MQC NEF Matrix variable empty (all zero) from NumSpace_row x NumSpace_col. </b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_matrix_allocate_blocks(NEFMat)
!
      implicit none
      type(mqc_nef_matrix),intent(inOut)::NEFMat
!
      if(allocated(NEFMat%Main_Matrix)) deallocate(NEFMat%Main_Matrix)
      if(allocated(NEFMat%Zero_Flag)) deallocate(NEFMat%Zero_Flag)
      allocate(NEFMat%Main_Matrix(NEFMat%NumSpace_row*NEFMat%NumSpace_col))
      allocate(NEFMat%Zero_Flag(size(NEFMat%Main_Matrix)))
      NEFMat%Zero_Flag = .true.
!
      end subroutine mqc_nef_matrix_allocate_blocks
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
!>    The indexes are shared by symmetric and asymmetric matrices (square or not)
!>
!>    Ex: no-sym
!>     5x5=               3x5=              5x3=
!>     | 1 16 17 19 22|   | 1  7  8 10 13|  | 1 13 14|
!>     | 2  3 18 20 23|   | 2  3  9 11 14|  | 2  3 15|
!>     | 4  5  6 21 24|   | 4  5  6 12 15|  | 4  5  6|
!>     | 7  8  9 10 25|                     | 7  8  9|
!>     |11 12 13 14 15|                     |10 11 12|
!> 
!>     If square and symm, the indexes of the upper-diagonal are returned based on the lower
!>    Ex: sym
!>     5x5=              
!>     | 1  2  4  7 11|  
!>     | 2  3  5  8 12|  
!>     | 4  5  6  9 13|  
!>     | 7  8  9 10 14|  
!>     |11 12 13 14 15|  
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
      if(.not.NEFMat%Mirrored) then
        if (id1.ge.id2) then
          if(NEFMat%NumSpace_row.eq.NEFMat%NumSpace_col.or.&
             id1.lt.NEFMat%NumSpace_col) then
            indexOut = id1*(id1-1)/2+id2
          else if(id1.ge.NEFMat%NumSpace_col) then
            indexOut =  NEFMat%NumSpace_col*(NEFMat%NumSpace_col-1)/2+&
                        (id1-NEFMat%NumSpace_col)*NEFMat%NumSpace_col+id2
          end if
        else
          if(NEFMat%NumSpace_row.eq.NEFMat%NumSpace_col.or.&
             id2.le.NEFMat%NumSpace_row) then
            indexOut = Min(NEFMat%NumSpace_row,NEFMat%NumSpace_col)*&
                      (Min(NEFMat%NumSpace_row,NEFMat%NumSpace_col)+1)/2+&
                       Max(0,NEFMat%NumSpace_row-NEFMat%NumSpace_col)*NEFMat%NumSpace_col
            indexOut = indexOut + (id2-2)*(id2-1)/2+id1
          else if(id2.gt.NEFMat%NumSpace_row) then
            indexOut =  Min(NEFMat%NumSpace_row,NEFMat%NumSpace_col)**2+&
                        (id2-NEFMat%NumSpace_row-1)*NEFMat%NumSpace_row+id1
          end if
        end if
      else
!       One block per dagger pair: (id1,id2) and (id2,id1) share a slot.
        indexOut = Max(id1,id2)*(Max(id1,id2)-1)/2+Min(id1,id2)
      end if
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
      if(.not.allocated(nefM%Zero_Flag)) call mqc_error(&
           'mqc_nef_matrix_at: NEF_Matrix is open (never initialized)',6)
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
      if (.not.nefM%Zero_Flag(nef_bidx)) then
        if (.not.nefM%Mirrored.or.myp.ge.myq) then
          element = nefM%Main_Matrix(nef_bidx)%at(myI,myJ,spinBlock)
        else
!         Object symmetry: the element is read from the stored (q,p) block.
          select case(spinBlock)
          case('alphabeta','ab')
            element = conjg(nefM%Main_Matrix(nef_bidx)%at(myJ,myI,'betaAlpha'))
          case('betaalpha','ba')
            element = conjg(nefM%Main_Matrix(nef_bidx)%at(myJ,myI,'alphaBeta'))
          case('beta','alpha','aa','bb','a','b')
            element = conjg(nefM%Main_Matrix(nef_bidx)%at(myJ,myI,spinBlock))
          case default
            call mqc_error_A('Spin-block unrecognized in mqc_NEF_Matrix_at',6,&
                      'spinBlock',spinBlock)
          end select
          if(nefM%Symmetry.eq.'antihermitian') element = (-1)*element
        end if
      end if
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
      logical::isMirror
      real(kind=real64)::zero=0.0d0
!
      element = zero
      if(.not.allocated(nefR4%Zero_Flag)) call mqc_error(&
           'mqc_nef_r4_at: NEF_R4 is open (never initialized)',6)
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
      call mqc_nef_r4_resolve(nefR4,[dI,I],[dJ,J],[dK,K],[dL,L],nef_bidx,isMirror)
      if (.not.nefR4%Zero_Flag(nef_bidx)) then
        if (.not.isMirror) then
          element = nefR4%Main_R4(nef_bidx)%R4Slot%at(&
                             id1,id2,id3,id4,spinBlock)
        else
!         Object symmetry: element of the dagger of the stored (S,R,Q,P) block.
          element = conjg(nefR4%Main_R4(nef_bidx)%R4Slot%at(id4,id3,id2,id1,&
                    spinBlock(4:4)//spinBlock(3:3)//spinBlock(2:2)//spinBlock(1:1)))
          if (nefR4%Symmetry.eq.'antihermitian') element = (-1)*element
        end if
      end if
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
!>    where s is 1 for alpha or -1 for beta). If the target block is
!>    absent (its MQC_NEF_R4_Slot is unallocated), put allocates it and
!>    initializes it with the block dimensions (Size_vector_*), then
!>    writes the element and marks the block non-zero. Depends on MQC_TwoERIs
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
      integer(kind=int64),dimension(4)::blk
      type(mqc_twoeris)::newBlock
      type(mqc_scalar)::mirElement
      logical::isMirror
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
      call mqc_nef_r4_resolve(nefR4,[dI,I],[dJ,J],[dK,K],[dL,L],nef_bidx,isMirror)
!     An absent (zero) block gets its slot here. Zero_Flag is the slot state:
!     'allocated' is a generic name in mqcPack (MQC_TwoERIs_isAllocated), so it
!     cannot test the allocation of the R4Slot component itself.
      if(nefR4%Zero_Flag(nef_bidx)) then
        blk(dI) = I
        blk(dJ) = J
        blk(dK) = K
        blk(dL) = L
!       the stored block of a mirror slot is (S,R,Q,P)
        if(isMirror) blk = blk(4:1:-1)
        call newBlock%init(&
             nefR4%Size_vector_alpha_1(blk(1)),nefR4%Size_vector_beta_1(blk(1)),&
             nAlpha2=nefR4%Size_vector_alpha_2(blk(2)),nBeta2=nefR4%Size_vector_beta_2(blk(2)),&
             nAlpha3=nefR4%Size_vector_alpha_3(blk(3)),nBeta3=nefR4%Size_vector_beta_3(blk(3)),&
             nAlpha4=nefR4%Size_vector_alpha_4(blk(4)),nBeta4=nefR4%Size_vector_beta_4(blk(4)))
        nefR4%Main_R4(nef_bidx)%R4Slot = newBlock
      end if
      if(.not.isMirror) then
        call nefR4%Main_R4(nef_bidx)%R4Slot%put(element,id1,id2,id3,id4,spinBlock)
      else
!       Object symmetry: the stored (S,R,Q,P) block gets the dagger-paired element.
        mirElement = conjg(element)
        if(nefR4%Symmetry.eq.'antihermitian') mirElement = (-1)*mirElement
        call nefR4%Main_R4(nef_bidx)%R4Slot%put(mirElement,id4,id3,id2,id1,&
             spinBlock(4:4)//spinBlock(3:3)//spinBlock(2:2)//spinBlock(1:1))
      end if
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
!>
!>    \param[in] swapOrder
!>    \verbatim
!>        swapOrder is Logical,Optional (default .false.)
!>        Must match the swapOrder passed to the Matrix collapse that
!>        originally produced this Vector. Collapse's swapOrder affects two
!>        things that expand must both undo: (1) which flat Main_Vector
!>        slot a given (row,col) block was written to -- outId=
!>        (col-1)*NumSpace_row+row (column-major) instead of row-major
!>        blockIndex; and (2) collapseIndex's inner (i,j)->flat-position
!>        formula within that slot (row-major vs column-major). Mismatching
!>        it silently reads the wrong flat positions (typically
!>        zero-flagged ones) instead of raising an error.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_vector_dimention_expand(&
         NEFV,map_matrix_row,map_matrix_col,spinType,swapOrder) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_vector),intent(in),target::NEFV
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,intent(in)::map_matrix_row,map_matrix_col
      character(len=*),intent(in),optional::spinType
      logical,intent(in),optional::swapOrder
      character(len=64)::mySpinType
      logical::mySwapOrder
      type(mqc_nef_matrix)::NEFMatOut
      integer(kind=int64),dimension(:),Allocatable::my_row,my_col
      integer(kind=int64)::row,col,i,j,vecId
      type(mqc_scf_eigenvalues),pointer::blkV
      type(mqc_scf_integral)::blkM,blkMEmpty
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
      call NEFMatOut%init(map_matrix_row,initial_map_matrix_col=map_matrix_col,&
           sameBraKet=NEFV%sameBraKet)
      select case(mySpinType)
      case ('singlet')
        do row = 1, NEFMatOut%NumSpace_row
          do col = 1, NEFMatOut%NumSpace_col
            if (mySwapOrder) then
              vecId = (col-1)*NEFMatOut%NumSpace_row + row
            else
              vecId = (row-1)*NEFMatOut%NumSpace_col + col
            end if
            if (NEFV%isBlockZero(vecId)) cycle
            blkV => NEFV%iref(vecId)
            blkM = blkMEmpty
            call blkM%init(&
                   NEFMatOut%Size_vector_alpha_row(row),&
                   NEFMatOut%Size_vector_beta_row(row),&
                   nAlpha2=NEFMatOut%Size_vector_alpha_col(col),&
                   nBeta2=NEFMatOut%Size_vector_beta_col(col))
            do i = 1, NEFMatOut%Size_vector_alpha_row(row)
              do j = 1, NEFMatOut%Size_vector_alpha_col(col)
                call blkM%put(&
                  blkV%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','a',row,col,'singlet',swapOrder=mySwapOrder)),'alpha'),&
                    i,j,'alpha')
              end do
            end do
            do i = 1, nefMatOut%Size_vector_beta_row(row)
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call blkM%put(&
                  blkV%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','b',row,col,'singlet',swapOrder=mySwapOrder)),'beta'),&
                    i,j,'beta')
              end do
            end do
!           An all-zero block ('fullsymmetric') stays a zero block.
            if (blkM%getSymm().ne.'fullsymmetric') call NEFMatOut%iput(blkM,row,col)
          end do
        end do
      case ('triplet')
        do row = 1, NEFMatOut%NumSpace_row
          do col = 1, NEFMatOut%NumSpace_col
            if (mySwapOrder) then
              vecId = (col-1)*NEFMatOut%NumSpace_row + row
            else
              vecId = (row-1)*NEFMatOut%NumSpace_col + col
            end if
            if (NEFV%isBlockZero(vecId)) cycle
            blkV => NEFV%iref(vecId)
            blkM = blkMEmpty
            call blkM%init(&
                   NEFMatOut%Size_vector_alpha_row(row),&
                   NEFMatOut%Size_vector_beta_row(row),&
                   nAlpha2=NEFMatOut%Size_vector_alpha_col(col),&
                   nBeta2=NEFMatOut%Size_vector_beta_col(col))
            do i = 1, nefMatOut%Size_vector_alpha_row(row)
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call blkM%put(&
                  blkV%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','b',row,col,'triplet',swapOrder=mySwapOrder)),'alpha'),&
                    i,j,'alphabeta')
              end do
            end do
            do i = 1, nefMatOut%Size_vector_beta_row(row)
              do j = 1, nefMatOut%Size_vector_alpha_col(col)
                call blkM%put(&
                  blkV%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','a',row,col,'triplet',swapOrder=mySwapOrder)),'beta'),&
                    i,j,'betaalpha')
              end do
            end do
!           An all-zero block ('fullsymmetric') stays a zero block.
            if (blkM%getSymm().ne.'fullsymmetric') call NEFMatOut%iput(blkM,row,col)
          end do
        end do
      case ('all','full')
        do row = 1, NEFMatOut%NumSpace_row
          do col = 1, NEFMatOut%NumSpace_col
            if (mySwapOrder) then
              vecId = (col-1)*NEFMatOut%NumSpace_row + row
            else
              vecId = (row-1)*NEFMatOut%NumSpace_col + col
            end if
            if (NEFV%isBlockZero(vecId)) cycle
            blkV => NEFV%iref(vecId)
            blkM = blkMEmpty
            call blkM%init(&
                   NEFMatOut%Size_vector_alpha_row(row),&
                   NEFMatOut%Size_vector_beta_row(row),&
                   nAlpha2=NEFMatOut%Size_vector_alpha_col(col),&
                   nBeta2=NEFMatOut%Size_vector_beta_col(col))
            do i = 1, nefMatOut%Size_vector_alpha_row(row)
              do j = 1, nefMatOut%Size_vector_alpha_col(col)
                call blkM%put(&
                  blkV%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','a',row,col,'full',swapOrder=mySwapOrder)),'alpha'),&
                    i,j,'alpha')
              end do
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call blkM%put(&
                  blkV%at(&
                    int(NEFMatOut%collapseIndex(i,j,'a','b',row,col,'full',swapOrder=mySwapOrder)),'alpha'),&
                    i,j,'alphabeta')
              end do
            end do
            do i = 1, nefMatOut%Size_vector_beta_row(row)
              do j = 1, nefMatOut%Size_vector_beta_col(col)
                call blkM%put(&
                  blkV%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','b',row,col,'full',swapOrder=mySwapOrder)),'beta'),&
                    i,j,'beta')
              end do
              do j = 1, nefMatOut%Size_vector_alpha_col(col)
                call blkM%put(&
                  blkV%at(&
                    int(NEFMatOut%collapseIndex(i,j,'b','a',row,col,'full',swapOrder=mySwapOrder)),'beta'),&
                    i,j,'betaalpha')
              end do
            end do
!           An all-zero block ('fullsymmetric') stays a zero block.
            if (blkM%getSymm().ne.'fullsymmetric') call NEFMatOut%iput(blkM,row,col)
          end do
        end do
      case default
        call mqc_error_a('mqc_nef_vector_dimention_expanded: Unrecognized excitation type',6,&
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
      class(mqc_nef_matrix),intent(inOut)::MatIn
      type(mqc_nef_matrix)::oldMat
      integer::oi,oj
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
 



      oldMat = MatIn
      call MatIn%init(tempMapRow, initial_map_matrix_col=tempMapCol)
      do i = 1, newNrow
        do j = 1, newNcol
          if (i==idIRow .or. j==idICol) cycle
          oi = i
          if (i > idIRow .and. idIRow>0) oi = i-1
          oj = j
          if (j > idICol .and. idICol>0) oj = j-1
          if (oldMat%isBlockZero(oi,oj)) cycle
          call MatIn%iput(oldMat%iat(oi,oj),i,j)
        end do
      end do
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
      integer(kind=int64)::row,col,i,j,outId
      type(mqc_scf_eigenvalues)::blkE,blkEEmpty
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
!     The vector index (row,col) is reversed by the dagger: the label is 'asymmetric'.
      call NEFVecOut%init(NEFM%collapseMap(mySpinType,swapOrder=mySwapOrder),&
           sameBraKet=NEFM%sameBraKet)
      select case(mySpinType)
      case ('singlet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            if (NEFM%isBlockZero(row,col)) cycle
            if (mySwapOrder) then
              outId = (col-1)*NEFM%NumSpace_row + row
            else
              outId = (row-1)*NEFM%NumSpace_col + col
            end if
            blkE = blkEEmpty
            call blkE%init(NEFVecOut%Size_vector_alpha(outId),&
                   lengthBeta=NEFVecOut%Size_vector_beta(outId))
            do i = 1, nefM%Size_vector_alpha_row(row)
              do j = 1, nefM%Size_vector_alpha_col(col)
                write(*,*) i,j
                call blkE%put(&
                  NEFM%at(i,j,row,col,'alpha'),&
                  int(NEFM%collapseIndex(i,j,'a','a',row,col,'singlet',swapOrder=mySwapOrder)),&
                  'alpha')
              end do
            end do
            do i = 1, nefM%Size_vector_beta_row(row)
              do j = 1, nefM%Size_vector_beta_col(col)
                call blkE%put(&
                  NEFM%at(i,j,row,col,'beta'),&
                  int(NEFM%collapseIndex(i,j,'b','b',row,col,'singlet',swapOrder=mySwapOrder)),&
                  'beta')
              end do
            end do
            call NEFVecOut%iput(blkE,outId)
          end do
        end do
      case ('triplet')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            if (NEFM%isBlockZero(row,col)) cycle
            if (mySwapOrder) then
              outId = (col-1)*NEFM%NumSpace_row + row
            else
              outId = (row-1)*NEFM%NumSpace_col + col
            end if
            blkE = blkEEmpty
            call blkE%init(NEFVecOut%Size_vector_alpha(outId),&
                   lengthBeta=NEFVecOut%Size_vector_beta(outId))
            do i = 1, nefM%Size_vector_alpha_row(row)
              do j = 1, nefM%Size_vector_beta_col(col)
                call blkE%put(&
                  NEFM%at(i,j,row,col,'alphabeta'),&
                  int(NEFM%collapseIndex(i,j,'a','b',row,col,'triplet',swapOrder=mySwapOrder)),&
                  'alpha')
              end do
            end do
            do i = 1, nefM%Size_vector_beta_row(row)
              do j = 1, nefM%Size_vector_alpha_col(col)
                call blkE%put(&
                  NEFM%at(i,j,row,col,'betaalpha'),&
                  int(NEFM%collapseIndex(i,j,'b','a',row,col,'triplet',swapOrder=mySwapOrder)),&
                  'beta')
              end do
            end do
            call NEFVecOut%iput(blkE,outId)
          end do
        end do
      case ('all')
        do row = 1, NEFM%NumSpace_row
          do col = 1, NEFM%NumSpace_col
            if (NEFM%isBlockZero(row,col)) cycle
            if (mySwapOrder) then
              outId = (col-1)*NEFM%NumSpace_row + row
            else
              outId = (row-1)*NEFM%NumSpace_col + col
            end if
            blkE = blkEEmpty
            call blkE%init(NEFVecOut%Size_vector_alpha(outId),&
                   lengthBeta=NEFVecOut%Size_vector_beta(outId))
            do i = 1, nefM%Size_vector_alpha_row(row)
              do j = 1, nefM%Size_vector_alpha_col(col)
                call blkE%put(&
                  NEFM%at(i,j,row,col,'alpha'),&
                  int(NEFM%collapseIndex(i,j,'a','a',row,col,'all',swapOrder=mySwapOrder)),&
                  'alpha')
              end do
              do j = 1, nefM%Size_vector_beta_col(col)
                call blkE%put(&
                  NEFM%at(i,j,row,col,'alphabeta'),&
                  int(NEFM%collapseIndex(i,j,'a','b',row,col,'all',swapOrder=mySwapOrder)),&
                  'alpha')
              end do
            end do
            do i = 1, nefM%Size_vector_beta_row(row)
              do j = 1, nefM%Size_vector_beta_col(col)
                call blkE%put(&
                  NEFM%at(i,j,row,col,'beta'),&
                  int(NEFM%collapseIndex(i,j,'b','b',row,col,'all',swapOrder=mySwapOrder)),&
                  'beta')
              end do
              do j = 1, nefM%Size_vector_alpha_col(col)
                call blkE%put(&
                  NEFM%at(i,j,row,col,'betaalpha'),&
                  int(NEFM%collapseIndex(i,j,'b','a',row,col,'all',swapOrder=mySwapOrder)),&
                  'beta')
              end do
            end do
            call NEFVecOut%iput(blkE,outId)
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
      class(mqc_nef_r4),intent(in),target::NEFr4
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      integer(kind=int64),dimension(4)::my_order
      type(mqc_nef_matrix)::NEFMatOut
      character(len=64)::mySpinType,mySpinTypeKet
      integer(kind=int64)::idij,idlk
      type(mqc_scf_integral)::blkM
      integer(kind=int64)::ia,ja,ka,la
      integer(kind=int64)::dimi,dimj,dimk,diml
      type(mqc_twoeris),pointer::ERIsP
      type(mqc_twoeris),target::ERIsMirror
      type(mqc_matrix)::MatAlpha,MatBeta,MatAlphaBeta,MatBetaAlpha
      integer(kind=int64),dimension(4)::blk
      integer(kind=int64),dimension(2)::nIs,nJs,nKs,nLs
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
!     The collapsed NEF_Matrix (element ((i,j),(k,l)) = R with axis Order(1..4) = i,j,k,l)
!     is the dagger-consistent image of R(p,q,r,s) = <I|a+_p a_q O a+_r a_s|J> only when the
!     row pair is the excitation q->p of the bra and the column pair r->s of the ket, i.e.
!     the dagger (axis x -> 5-x) maps them onto each other: ka=5-ia, la=5-ja (e.g. Order
!     [2,1,3,4]; user, 2026-11-09), with the same spin type on bra and ket. Checked for the
!     singlet, triplet and all collapses (triplet/all column pair fixed 2026-11-09).
      call NEFMatOut%init(NEFr4%collapseMap(mySpinType,[ia,ja]),&
         initial_map_matrix_col=NEFr4%collapseMap(mySpinTypeKet,[ka,la]),&
           symmetry=merge(NEFr4%Symmetry,'asymmetric      ',ka.eq.5-ia.and.la.eq.5-ja.and.&
           mySpinType.eq.mySpinTypeKet),sameBraKet=NEFr4%sameBraKet)
!     Unknown types stop here (collapseMap above already checked the ket type).
      select case(mySpinType)
      case ('singlet','triplet','all')
        continue
      case default
        call mqc_error_a('mqc_nef_r4_dimention_collapse: Unrecognized excitation type',6,&
          'mySpinType',mySpinType)
      end select
!
!     Block-wise: each stored space block is read once per spin channel as a dense
!     array and copied into the collapsed blocks by array index arithmetic. Pair
!     (i,j) -> row, (k,l) -> column; the row/column spin part is the spin of i/k
!     (collapseMap), so a channel goes to alpha (a,a), betaalpha (row a, column b),
!     alphabeta (row b, column a) or beta (b,b). The second index of a pair has
!     the same spin (singlet), the opposite spin (triplet) or both, beta first
!     (all), as in collapseIndex. spinType sets the rows, spinTypeKet the columns.
      do dimi = 1, NEFr4%blockSize(ia)
       do dimj = 1, NEFr4%blockSize(ja)
        idij = NEFr4%blockSize(ja)*(dimi-1) + dimj
        do dimk = 1, NEFr4%blockSize(ka)
         do diml = 1, NEFr4%blockSize(la)
          idlk = NEFr4%blockSize(la)*(dimk-1) + diml
          if (NEFr4%isBlockZero([ia,dimi],[ja,dimj],[la,diml],[ka,dimk])) cycle
!         whole stored block: pointer for own slots, dagger copy for mirror slots
          blk(ia) = dimi
          blk(ja) = dimj
          blk(ka) = dimk
          blk(la) = diml
          if (NEFr4%isMirrorSlot(int(blk(1)),int(blk(2)),int(blk(3)),int(blk(4)))) then
            ERIsMirror = NEFr4%iat([ia,dimi],[ja,dimj],[ka,dimk],[la,diml])
            ERIsP => ERIsMirror
          else
            ERIsP => NEFr4%iref([ia,dimi],[ja,dimj],[ka,dimk],[la,diml])
          end if
          nIs = [NEFr4%blockSize(ia,dimi,'alpha'),NEFr4%blockSize(ia,dimi,'beta')]
          nJs = [NEFr4%blockSize(ja,dimj,'alpha'),NEFr4%blockSize(ja,dimj,'beta')]
          nKs = [NEFr4%blockSize(ka,dimk,'alpha'),NEFr4%blockSize(ka,dimk,'beta')]
          nLs = [NEFr4%blockSize(la,diml,'alpha'),NEFr4%blockSize(la,diml,'beta')]
          call mqc_nef_r4_collapse_part(ERIsP,my_order,mySpinType,mySpinTypeKet,'a','a',&
            nIs(1),nJs,nKs(1),nLs,MatAlpha)
          call mqc_nef_r4_collapse_part(ERIsP,my_order,mySpinType,mySpinTypeKet,'a','b',&
            nIs(1),nJs,nKs(2),nLs,MatBetaAlpha)
          call mqc_nef_r4_collapse_part(ERIsP,my_order,mySpinType,mySpinTypeKet,'b','a',&
            nIs(2),nJs,nKs(1),nLs,MatAlphaBeta)
          call mqc_nef_r4_collapse_part(ERIsP,my_order,mySpinType,mySpinTypeKet,'b','b',&
            nIs(2),nJs,nKs(2),nLs,MatBeta)
          ERIsP => null()
!         An all-zero block stays a zero block (no matrix allocated).
          if (MQC_Matrix_isAllocated(MatAlpha).or.MQC_Matrix_isAllocated(MatBeta).or.&
              MQC_Matrix_isAllocated(MatAlphaBeta).or.MQC_Matrix_isAllocated(MatBetaAlpha)) then
            call mqc_integral_allocate(blkM,'',&
              merge('general','spin   ',MQC_Matrix_isAllocated(MatAlphaBeta).or.&
                                        MQC_Matrix_isAllocated(MatBetaAlpha)),&
              alpha=MatAlpha,beta=MatBeta,alphaBeta=MatAlphaBeta,betaAlpha=MatBetaAlpha,&
              nDimA1=int(NEFMatOut%Size_vector_alpha_row(idij),kind=int64),&
              nDimB1=int(NEFMatOut%Size_vector_beta_row(idij),kind=int64),&
              nDimA2=int(NEFMatOut%Size_vector_alpha_col(idlk),kind=int64),&
              nDimB2=int(NEFMatOut%Size_vector_beta_col(idlk),kind=int64))
            call NEFMatOut%iput(blkM,idij,idlk)
          end if
         end do
        end do
       end do
      end do
!
      end function mqc_nef_r4_dimention_collapse
!
!
!     PROCEDURE MQC_NEF_R4_Collapse_Part
!
!>    \brief <b> MQC_NEF_R4_Collapse_Part (private) returns one spin part of a
!>    collapsed NEF_R4 space block as a dense MQC matrix. </b>
!>
!>    Rows are the pairs (i,j) with i of spin RowSpin, columns the pairs (k,l) with
!>    k of spin ColSpin; i,j,k,l lie on the storage axes Ax = [ia,ja,ka,la]. The
!>    spins of j (l) follow RowType (ColType): singlet = same, triplet = opposite,
!>    all = beta then alpha. Element ((i-1)*nJ+jj,(k-1)*nL+ll) = T(i,j,k,l), jj (ll)
!>    the position of j (l) in its pair, as in collapseIndex. nI, nK are the sizes
!>    of the i, k spins; nJ, nL = [alpha,beta] sizes of the j, l spaces. As in the
!>    element-wise put (lazy allocation), entries before the first |T| >
!>    MQC_very_small in row-major order are zeroed, and MatOut stays unallocated
!>    when no entry exceeds MQC_very_small. MatOut is not classified (StorFull /
!>    StorNoSym), as an element-wise put leaves it.
      subroutine mqc_nef_r4_collapse_part(ERIs,Ax,RowType,ColType,RowSpin,ColSpin,&
                                          nI,nJ,nK,nL,MatOut)
!
      implicit none
      type(mqc_twoeris),intent(in)::ERIs
      integer(kind=int64),dimension(4),intent(in)::Ax
      character(len=*),intent(in)::RowType,ColType
      character(len=1),intent(in)::RowSpin,ColSpin
      integer(kind=int64),intent(in)::nI,nK
      integer(kind=int64),dimension(2),intent(in)::nJ,nL
      type(mqc_matrix),intent(out)::MatOut
      character(len=1),dimension(2)::jSpins,lSpins
      integer(kind=int64),dimension(2)::jOff,lOff,jLen,lLen
      integer::njs,nls,js,ls
      integer(kind=int64)::jExt,lExt,NRow,NCol,r,c
      real(kind=real64),dimension(:,:,:,:),allocatable::DR,QR
      complex(kind=real64),dimension(:,:,:,:),allocatable::DC,QC
      real(kind=real64),dimension(:,:),allocatable::MR,AbsM
      complex(kind=real64),dimension(:,:),allocatable::MC
      logical::IsCmplx,QIsCmplx,HasData
!
      call mqc_nef_r4_pair_spins(RowType,RowSpin,nJ,jSpins,jOff,jLen,njs,jExt)
      call mqc_nef_r4_pair_spins(ColType,ColSpin,nL,lSpins,lOff,lLen,nls,lExt)
      if (nI*jExt*nK*lExt.gt.0) then
        allocate(QR(jExt,nI,lExt,nK))
        QR = 0.0
        QIsCmplx = .false.
        do js = 1, njs
          do ls = 1, nls
            call mqc_nef_r4_collapse_channel(ERIs,[RowSpin,jSpins(js),ColSpin,lSpins(ls)],&
              Ax,DR,DC,IsCmplx)
            if (IsCmplx) then
              if (any(shape(DC,kind=int64).ne.[jLen(js),nI,lLen(ls),nK])) &
                call mqc_error_a('mqc_nef_r4_collapse_part: spin block size mismatch',6,&
                  'channel (i,j,k,l)',RowSpin//jSpins(js)//ColSpin//lSpins(ls))
            else
              if (any(shape(DR,kind=int64).ne.[jLen(js),nI,lLen(ls),nK])) &
                call mqc_error_a('mqc_nef_r4_collapse_part: spin block size mismatch',6,&
                  'channel (i,j,k,l)',RowSpin//jSpins(js)//ColSpin//lSpins(ls))
            end if
            if (IsCmplx.and..not.QIsCmplx) then
              QC = cmplx(QR,kind=real64)
              deallocate(QR)
              QIsCmplx = .true.
            end if
            if (QIsCmplx) then
              if (IsCmplx) then
                QC(jOff(js)+1:jOff(js)+jLen(js),:,lOff(ls)+1:lOff(ls)+lLen(ls),:) = DC
              else
                QC(jOff(js)+1:jOff(js)+jLen(js),:,lOff(ls)+1:lOff(ls)+lLen(ls),:) = &
                  cmplx(DR,kind=real64)
              end if
            else
              QR(jOff(js)+1:jOff(js)+jLen(js),:,lOff(ls)+1:lOff(ls)+lLen(ls),:) = DR
            end if
          end do
        end do
!       (jj,i,ll,k) flattened column-major: row (i-1)*jExt+jj, column (k-1)*lExt+ll
        NRow = jExt*nI
        NCol = lExt*nK
        if (QIsCmplx) then
          MC = reshape(QC,[NRow,NCol])
          AbsM = abs(MC)
        else
          MR = reshape(QR,[NRow,NCol])
          AbsM = abs(MR)
        end if
        HasData = .false.
        Rows: do r = 1, NRow
          do c = 1, NCol
            if (AbsM(r,c).gt.MQC_very_small) then
              HasData = .true.
              exit Rows
            end if
          end do
        end do Rows
        if (HasData) then
          if (QIsCmplx) then
            MC(1:r-1,:) = (0.0,0.0)
            MC(r,1:c-1) = (0.0,0.0)
            call MQC_Array_To_Matrix(MatOut,MC,Classify=.false.)
          else
            MR(1:r-1,:) = 0.0
            MR(r,1:c-1) = 0.0
            call MQC_Array_To_Matrix(MatOut,MR,Classify=.false.)
          end if
        end if
      end if
!
      end subroutine mqc_nef_r4_collapse_part
!
!
!     PROCEDURE MQC_NEF_R4_Pair_Spins
!
!>    \brief <b> MQC_NEF_R4_Pair_Spins (private) returns the spins of the second index
!>    of a collapsed pair whose first index has spin FirstSpin. </b>
!>
!>    singlet: the same spin; triplet: the opposite spin; all: beta then alpha (the
!>    order of collapseIndex). Off, Len are the offset and size of each spin in the
!>    pair position, Ext their total; nSecond = [alpha,beta] sizes of the second space.
      subroutine mqc_nef_r4_pair_spins(ExcType,FirstSpin,nSecond,Spins,Off,Len,nSpins,Ext)
!
      implicit none
      character(len=*),intent(in)::ExcType
      character(len=1),intent(in)::FirstSpin
      integer(kind=int64),dimension(2),intent(in)::nSecond
      character(len=1),dimension(2),intent(out)::Spins
      integer(kind=int64),dimension(2),intent(out)::Off,Len
      integer,intent(out)::nSpins
      integer(kind=int64),intent(out)::Ext
      integer::m
!
      Spins = ' '
      select case(ExcType)
      case('singlet')
        nSpins = 1
        Spins(1) = FirstSpin
      case('triplet')
        nSpins = 1
        Spins(1) = merge('b','a',FirstSpin.eq.'a')
      case('all','full')
        nSpins = 2
        Spins = ['b','a']
      case default
        call mqc_error_a('mqc_nef_r4_pair_spins: Unrecognized excitation type',6,&
          'ExcType',ExcType)
      end select
      Off = 0
      Len = 0
      Ext = 0
      do m = 1, nSpins
        Off(m) = Ext
        Len(m) = merge(nSecond(1),nSecond(2),Spins(m).eq.'a')
        Ext = Ext + Len(m)
      end do
!
      end subroutine mqc_nef_r4_pair_spins
!
!
!     PROCEDURE MQC_NEF_R4_Collapse_Channel
!
!>    \brief <b> MQC_NEF_R4_Collapse_Channel (private) returns one spin block of a
!>    stored NEF_R4 space block as a dense array in collapse order (j,i,l,k). </b>
!>
!>    Spins(1:4) are the spins of i,j,k,l ('a'/'b'), on the storage axes Ax =
!>    [ia,ja,ka,la]; the spin block is read with getBlock (unfolded with FullView).
!>    The result is in DC (IsCmplx) or DR.
      subroutine mqc_nef_r4_collapse_channel(ERIs,Spins,Ax,DR,DC,IsCmplx)
!
      implicit none
      type(mqc_twoeris),intent(in)::ERIs
      character(len=1),dimension(4),intent(in)::Spins
      integer(kind=int64),dimension(4),intent(in)::Ax
      real(kind=real64),dimension(:,:,:,:),allocatable,intent(out)::DR
      complex(kind=real64),dimension(:,:,:,:),allocatable,intent(out)::DC
      logical,intent(out)::IsCmplx
      type(mqc_r4tensor)::T
      real(kind=real64),dimension(:,:,:,:),allocatable::SR
      complex(kind=real64),dimension(:,:,:,:),allocatable::SC
      character(len=4)::BlockName
      integer(kind=int64),dimension(4)::NewAx,NewShape,ROrder,NStor
      integer::m
!
      do m = 1, 4
        BlockName(Ax(m):Ax(m)) = Spins(m)
      end do
      T = ERIs%getBlock(BlockName)
      IsCmplx = mqc_have_complex(T)
      if (IsCmplx) then
        call T%getArray(SC,FullView=.true.)
        NStor = shape(SC,kind=int64)
      else
        call T%getArray(SR,FullView=.true.)
        NStor = shape(SR,kind=int64)
      end if
!     target axis order (j,i,l,k); ROrder maps each storage axis to its target position
      NewAx = [Ax(2),Ax(1),Ax(4),Ax(3)]
      do m = 1, 4
        NewShape(m) = NStor(NewAx(m))
        ROrder(NewAx(m)) = m
      end do
      if (IsCmplx) then
        DC = reshape(SC,NewShape,order=ROrder)
      else
        DR = reshape(SR,NewShape,order=ROrder)
      end if
!
      end subroutine mqc_nef_r4_collapse_channel
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
!>
!>    \param spaceDim
!>    \verbatim
!>        spaceDim is Integer,Optional
!>        The number of linear-independent vectors.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_block_diagonal_pinv(NEFr4,spinType,Order,thresh,spaceDim) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_r4),intent(in)::NEFr4
      character(len=*),intent(in),optional::spinType
      integer(kind=int64),dimension(4),intent(in),optional::Order
      real(kind=real64),intent(in),optional::thresh
      Integer(kind=int64),optional::spaceDim
      type(mqc_nef_matrix)::NEFMatOut
      type(mqc_nef_matrix)::CollapsedMat
      real(kind=real64)::myThresh
      integer(kind=int64)::N,id
      Integer(kind=int64)::mySpaceDim
!
      if(present(thresh)) then
        myThresh = thresh
      else
        myThresh = 1.0d-8
      end if
!
      CollapsedMat = NEFr4%collapse(spinType=spinType,Order=Order)
!
      NEFMatOut = CollapsedMat%pinv(myThresh,mySpaceDim)
      if(present(spaceDim)) spaceDim = mySpaceDim
!
      end function mqc_nef_r4_block_diagonal_pinv
!
!
!     PROCEDURE MQC_NEF_Set_Block_Diagonal_PInv
!
!>    \brief <b> MQC_NEF_Set_Block_Diagonal_PInv is used to return the block-diagonal
!>    pseudo-inverse of a NEF Set-R4 object, expressed as a collapsed NEF Matrix. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4
!>    \verbatim
!>        NEFSet is type(mqc_nef_set)
!>        NEF Set containing R4 Object whose block-diagonal pseudo-inverse is to be returned.
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
!>
!>    \param spaceDim
!>    \verbatim
!>        spaceDim is Integer,Optional
!>        The number of linear-independent vectors.
!>    \endverbatim
!>
!>    \param[in] wholeMatrix
!>    \verbatim
!>        wholeMatrix is Logical,Optional (default .false.)
!>        .false.: pseudo-invert each diagonal space block of the merged matrix
!>                 (off-diagonal space blocks are ignored).
!>        .true. : pseudo-invert the whole merged matrix, off-diagonal space blocks
!>                 included (mqc_nef_matrix_whole_pinv); exact when the space
!>                 partition does not block-diagonalize the matrix.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_block_diagonal_pinv(NEFSet,spinType,Order,thresh,spaceDim,wholeMatrix) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      character(len=*),intent(in),optional::spinType
      integer(kind=int64),dimension(4),intent(in),optional::Order
      real(kind=real64),intent(in),optional::thresh
      Integer(kind=int64),optional::spaceDim
      logical,intent(in),optional::wholeMatrix
      type(mqc_nef_matrix)::TmpMat
      type(mqc_nef_set)::CollapsedMat,NEFMatOut
      real(kind=real64)::myThresh
      Integer(kind=int64)::mySpaceDim
      logical::myWhole
!
      if(present(thresh)) then
        myThresh = thresh
      else
        myThresh = 1.0d-8
      end if
      myWhole = .false.
      if(present(wholeMatrix)) myWhole = wholeMatrix
!
      if(.not.NEFset%isR4()) &
        call mqc_error_A('mqc_nef_set_block_diagonal_pinv: NEFSet is not R4-type ',6,&
                       'NEFSet%type()',NEFSet%type())
!
      CollapsedMat = NEFSet%collapse(spinType=spinType,Order=Order)
!
      if(CollapsedMat%NumSpace(1).ne.CollapsedMat%NumSpace(2)) &
        call mqc_error_i('mqc_nef_set_block_diagonal_pinv: collapsed set is not square',6,&
          'NumSpace(1)',CollapsedMat%NumSpace(1),'NumSpace(2)',CollapsedMat%NumSpace(2))
!
      CollapsedMat = CollapsedMat%mergeSpace()
!     the pseudo-inverse keeps the label and bra-ket flag of the merged cell
      call NEFMatOut%init(1,1,sameBraKetSpace=CollapsedMat%isSameBraKetSpace(),&
           symmetry=CollapsedMat%getSymm())
      call TmpMat%init(CollapsedMat%Map_Matrix_1,&
                       initial_map_matrix_col=CollapsedMat%Map_Matrix_2,&
                       symmetry=CollapsedMat%getSymm(),&
                       sameBraKet=CollapsedMat%isSameBraKetSpace())
!
      mySpaceDim=0
      if (.not.CollapsedMat%Zero_Flag(1,1)) then
        if (myWhole) then
          TmpMat = mqc_nef_matrix_whole_pinv(CollapsedMat%Set_Matrix(1,1),myThresh,mySpaceDim)
        else
          TmpMat = CollapsedMat%Set_Matrix(1,1)%pinv(myThresh,mySpaceDim)
        end if
      end if
      if (mySpaceDim.ne.0) then
        call NEFMatOut%putMatrix(TmpMat,1,1)
      else
        call mqc_nef_set_allocate(NEFMatOut,'matrix',Map_Matrix_1=CollapsedMat%Map_Matrix_1,&
                                                     Map_Matrix_2=CollapsedMat%Map_Matrix_2)
      end if
      if (present(spaceDim)) spaceDim = mySpaceDim
! 
      end function mqc_nef_set_block_diagonal_pinv
!
!
!     PROCEDURE MQC_NEF_Matrix_Whole_PInv
!
!>    \brief <b> MQC_NEF_Matrix_Whole_PInv (private) returns the pseudo-inverse of
!>    the whole NEF Matrix, all space blocks coupled. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    The space blocks are assembled into one full matrix (per space: alpha then
!>    beta rows / columns, the layout of getBlock('full')), pseudo-inverted with
!>    the singular-value threshold thresh, and split back into 'general' space
!>    blocks (alphaBeta = beta rows x alpha columns, as mqc_integral_split_full).
!>    The output rows follow the input columns and vice versa; it keeps the
!>    operator label and the bra-ket flag of the input (a Mirrored output gets
!>    only its lower blocks). Empty spaces are skipped; a space holding only
!>    alpha or only beta functions is not supported (error).
!>
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_whole_pinv(NEFMat,thresh,spaceDim) result(NEFMatOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::NEFMat
      real(kind=real64),intent(in)::thresh
      integer(kind=int64),intent(out)::spaceDim
      type(mqc_nef_matrix)::NEFMatOut
      type(mqc_scf_integral)::blk,outInt
      type(mqc_matrix)::blkFull,wholeMat,wholePinv,sub,pa,pb,pab,pba
      integer(kind=int64),dimension(:),allocatable::naR,nbR,offR,naC,nbC,offC
      integer(kind=int64)::r,c,i,j,nT,mT,nRow,nCol
!
      allocate(naR(NEFMat%NumSpace_row),nbR(NEFMat%NumSpace_row),offR(NEFMat%NumSpace_row))
      allocate(naC(NEFMat%NumSpace_col),nbC(NEFMat%NumSpace_col),offC(NEFMat%NumSpace_col))
      nT = 0
      do r = 1, NEFMat%NumSpace_row
        naR(r) = NEFMat%Size_vector_alpha_row(r)
        nbR(r) = NEFMat%Size_vector_beta_row(r)
        offR(r) = nT
        nT = nT + naR(r) + nbR(r)
      end do
      mT = 0
      do c = 1, NEFMat%NumSpace_col
        naC(c) = NEFMat%Size_vector_alpha_col(c)
        nbC(c) = NEFMat%Size_vector_beta_col(c)
        offC(c) = mT
        mT = mT + naC(c) + nbC(c)
      end do
!     assemble
      call wholeMat%init(nT,mT,(0.0,0.0))
      do r = 1, NEFMat%NumSpace_row
        do c = 1, NEFMat%NumSpace_col
          if (naR(r)+nbR(r).eq.0.or.naC(c)+nbC(c).eq.0) cycle
          if (NEFMat%isBlockZero(r,c)) cycle
          blk = NEFMat%iat(r,c)
          blkFull = blk%getBlock('full')
          do i = 1, naR(r)+nbR(r)
            do j = 1, naC(c)+nbC(c)
              call wholeMat%put(blkFull%at(i,j),offR(r)+i,offC(c)+j)
            end do
          end do
        end do
      end do
      wholePinv = wholeMat%pinv(thresh,spaceDim=spaceDim)
!     split back: output space (r,c) = rows of input column space r, columns of input row space c
      call NEFMatOut%init(NEFMat%Map_Matrix_Col,initial_map_matrix_col=NEFMat%Map_Matrix_Row,&
           symmetry=NEFMat%Symmetry,sameBraKet=NEFMat%sameBraKet)
      do r = 1, NEFMat%NumSpace_col
        do c = 1, NEFMat%NumSpace_row
          if (NEFMatOut%Mirrored.and.c.gt.r) cycle
          nRow = naC(r) + nbC(r)
          nCol = naR(c) + nbR(c)
          if (nRow.eq.0.or.nCol.eq.0) cycle
          sub = wholePinv%mat([offC(r)+1,offC(r)+nRow],[offR(c)+1,offR(c)+nCol])
          if (sub%norm().eq.0.0) cycle
          if (naC(r).eq.0.or.nbC(r).eq.0.or.naR(c).eq.0.or.nbR(c).eq.0) &
            call mqc_error_i('mqc_nef_matrix_whole_pinv: a space with only alpha or only beta '//&
              'functions is not supported',6,'r',r,'c',c)
          pa = sub%mat([1,naC(r)],[1,naR(c)])
          pb = sub%mat([naC(r)+1,nRow],[naR(c)+1,nCol])
          pab = sub%mat([naC(r)+1,nRow],[1,naR(c)])
          pba = sub%mat([1,naC(r)],[naR(c)+1,nCol])
          call mqc_integral_allocate(outInt,'','general',pa,pb,pab,pba)
          call NEFMatOut%iput(outInt,r,c)
        end do
      end do
!
      end function mqc_nef_matrix_whole_pinv
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
      integer(kind=int64)::jg
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
!       swapOrder: j outer (betas first, then alphas), i inner within the vector spin
!       block of i (was: i_spin switch with a row stride, which overlapped entries)
        if (j_spin.eq.'b') then
          jg = j
        else
          jg = j + nefM%Size_vector_beta_col(j_space)
        end if
        if (i_spin.eq.'b') then
          idOut = (jg-1)*nefM%Size_vector_beta_row(i_space) + i
        else
          idOut = (jg-1)*nefM%Size_vector_alpha_row(i_space) + i
        end if
      else
!       no swap: each i owns nb+na entries (betas first); the beta stride is nb+na too
!       (was nb, which overlapped the alpha entries for i >= 2)
        if (j_spin.eq.'b') then
          idOut = (i-1)*(nefM%Size_vector_beta_col(j_space)+nefM%Size_vector_alpha_col(j_space)) + j
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
    case ('singlet','triplet')
!     row-major over (i,j): the stride is the size of j's own spin (triplet: the
!     alpha part (i a,j b) has n_j^beta columns per i, as in collapseMap; it used
!     the opposite spin, wrong when n_j^alpha /= n_j^beta)
      if (j_spin.eq.'b') then
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'beta') + j
      else
        idOut = (i-1)*nefM%blockSize(dimj,j_space,'alpha') + j
      end if 
    case ('all','full')
!     each i owns nb+na entries (betas first): the beta stride is nb+na too (was nb,
!     which overlapped the alpha entries for i >= 2)
      if (j_spin.eq.'b') then
        idOut = (i-1)*(nefM%blockSize(dimj,j_space,'beta') + &
                       nefM%blockSize(dimj,j_space,'alpha')) + j
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
              outId = (row-1)*NEFM%NumSpace_col + col
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
              outId = (row-1)*NEFM%NumSpace_col + col
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
              outId = (row-1)*NEFM%NumSpace_col + col
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
            id = NEFr4%blockSize(ja)*(row-1) + col
            call MapOut(id)%init(&
               nefR4%blockSize(ia,row,'alpha')*nefR4%blockSize(ja,col,'alpha'),0,&
               nefR4%blockSize(ia,row,'beta')*nefR4%blockSize(ja,col,'beta'))
          end do
        end do
      case ('triplet')
        do row = 1, NEFr4%blockSize(ia)
          do col = 1, NEFr4%blockSize(ja)
            id = NEFr4%blockSize(ja)*(row-1) + col
            call MapOut(id)%init(&
               nefR4%blockSize(ia,row,'alpha')*nefR4%blockSize(ja,col,'beta'),0,&
               nefR4%blockSize(ia,row,'beta')*nefR4%blockSize(ja,col,'alpha'))
          end do
        end do
      case ('all','full')
        do row = 1, NEFr4%blockSize(ia)
          do col = 1, NEFr4%blockSize(ja)
            id = NEFr4%blockSize(ja)*(row-1) + col
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
      if(.not.allocated(nefV%Zero_Flag)) call mqc_error(&
           'mqc_nef_vector_output_block_vector: NEF_Vector is open (never initialized)',6)
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
      type(mqc_scf_integral)::tmpBlock
      character(len=64)::myBlockName,myBlockOrder
      integer(kind=int64),dimension(:),allocatable::myId,myIdKet
      logical::full,fullKet,blockZero
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
              blockZero = nefMat%isBlockZero(i,j)
              if (.not.blockZero) tmpBlock = nefMat%iat(i,j)
              if(myBlockOrder.eq.'alpha'.and.myBlockName.eq.'full') then
                size1 = nefMat%Size_vector_alpha_row(i)+nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)+nefMat%Size_vector_beta_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                           [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('full'),&
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
                  if (blockZero) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                       [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('beta'),&
                                       [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                if (size3.gt.0.and.size2.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size3,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('betaalpha'),&
                                 [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                if (size1.gt.0.and.size4.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size1,size4,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('alphabeta'),&
                                 [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  end if
                end if
                if (size3.gt.0.and.size4.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size3,size4,0.0)
                    call MatrixOut%mput(zeroBlock,&
                     [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('alpha'),&
                     [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2 + size4
                size1 = size1 + size3
              else if (myBlockName.eq.'alpha') then
                size1 = nefMat%Size_vector_alpha_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('alpha'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'beta') then
                size1 = nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_beta_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('beta'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'alphabeta') then
                size1 = nefMat%Size_vector_alpha_row(i)
                size2 = nefMat%Size_vector_beta_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('alphabeta'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'betaalpha') then
                size1 = nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('betaalpha'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                    call MatrixOut%mput(tmpBlock%getBlock('full'),&
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
                  if (blockZero) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                       [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('beta'),&
                                       [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                if (size3.gt.0.and.size2.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size3,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('betaalpha'),&
                                 [fsize+size1,fsize+size1+size3-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                if (size1.gt.0.and.size4.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size1,size4,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                 [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('alphabeta'),&
                                 [fsize,fsize+size1-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  end if
                end if
                if (size3.gt.0.and.size4.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size3,size4,0.0)
                    call MatrixOut%mput(zeroBlock,&
                     [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('alpha'),&
                     [fsize+size1,fsize+size1+size3-1],[fsizeKet+size2,fsizeKet+size2+size4-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2 + size4
                size1 = size1 + size3
              else if (myBlockName.eq.'alpha') then
                size1 = nefMat%Size_vector_alpha_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)
                if (size1.gt.0.and.size2.gt.0) then
                  if (blockZero) then
                    call zeroBlock%init(size1,size2,0.0)
                    call MatrixOut%mput(zeroBlock,&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  else
                    call MatrixOut%mput(tmpBlock%getBlock('alpha'),&
                                     [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                  end if
                end if
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'beta') then
                fsizeKet = fsizeKet + size2
              else if (myBlockName.eq.'betaalpha') then
                size1 = nefMat%Size_vector_beta_row(i)
                size2 = nefMat%Size_vector_alpha_col(j)
                if (blockZero) then
                  call zeroBlock%init(size1,size2,0.0)
                  call MatrixOut%mput(zeroBlock,&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                else
                  call MatrixOut%mput(tmpBlock%getBlock('betaalpha'),&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                end if
                fsizeKet = fsizeKet + size2
                if (blockZero) then
                  call zeroBlock%init(size1,size2,0.0)
                  call MatrixOut%mput(zeroBlock,&
                                   [fsize,fsize+size1-1],[fsizeKet,fsizeKet+size2-1])
                else
                  call MatrixOut%mput(tmpBlock%getBlock('betaalpha'),&
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
!>    Map_Matrix_Col, selected via the optional axis argument
!>    (1=row, 2=col). Defaults to 1 (row) when axis is absent.
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
!>    \param[in] axis
!>    \verbatim
!>        axis is Integer,Optional
!>        1 (default) for row, 2 for col.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_get_map(NEFMat,axis) result(mapOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      integer,intent(in),optional::axis
      integer::my_axis
      type(mqc_scf_eigenvalues),dimension(:),allocatable::mapOut
!
      if(present(axis)) then
        my_axis = axis
      else
        my_axis = 1
      end if
!
      select case(my_axis)
      case(1)
        mapOut = NEFMat%Map_Matrix_row
      case(2)
        mapOut = NEFMat%Map_Matrix_col
      case default
        call mqc_error_I('mqc_nef_matrix_get_map: whichIn must be 1 (row) or 2 (col)',&
             6,'my_axis',my_axis)
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
!>
!>    \param spaceDim
!>    \verbatim
!>        spaceDim is Integer,Optional
!>        The number of linear-independent vectors.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_block_diagonal_pinv(NEFMat,thresh,spaceDim) result(NEFMatOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      real(kind=real64),intent(in),optional::thresh
      Integer(kind=int64),optional::spaceDim
      type(mqc_nef_matrix)::NEFMatOut
      real(kind=real64)::myThresh
      integer(kind=int64)::N
      type(mqc_scf_integral)::tmpBlock,pinvBlock
      Integer(kind=int64)::mySpaceDim,tmpspacedim
!
      if(present(thresh)) then
        myThresh = thresh
      else
        myThresh = 1.0d-8
      end if
!
!     The pseudo-inverse of a (anti)hermitian operator keeps its label.
      call NEFMatOut%init(NEFMat%Map_Matrix_Row,&
           initial_map_matrix_col=NEFMat%Map_Matrix_Col,&
           symmetry=NEFMat%Symmetry,sameBraKet=NEFMat%sameBraKet)
!
      mySpaceDim = 0
      do N = 1, NEFMat%NumSpace_row
        if (.not.NEFMat%isBlockZero(N,N)) then
          tmpBlock = NEFMat%iat(N,N)
          pinvBlock = tmpBlock%pinv(myThresh,spaceDim=tmpSpaceDim)
          if(tmpSpaceDim.ne.0) then
            call NEFMatOut%iput(pinvBlock,N,N)
            mySpaceDim = mySpaceDim + tmpSpaceDim
          end if
        end if
      end do
      if (present(spaceDim)) spaceDim = mySpaceDim
!
      end function mqc_nef_matrix_block_diagonal_pinv
!
!
!     PROCEDURE MQC_NEF_Matrix_Get_Rotation
!
!>    \brief <b> MQC_NEF_Matrix_Get_Rotation builds a NEF rotation matrix
!>    as a product of elementary orbital-pair rotations, one per nonzero
!>    entry of an input NEF_Matrix whose own values and sign pattern
!>    specify both the rotation angle and the generator type at each
!>    orbital pair -- mirroring SNOWMAN_NEF's rotateSolutions/rotateMO
!>    (OperatorsMod.f03). </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    Unlike a diagonalization routine, this does not compute anything
!>    from an input matrix's own structure to drive some entry to zero --
!>    it takes externally-specified rotation angles and builds the
!>    corresponding finite rotation as a product of elementary two-orbital
!>    rotations, exactly as SNOWMAN_NEF's rotateSolutions/rotateMO do when
!>    applying a quasi-Newton orbital-optimization step. SNOWMAN_NEF uses
!>    8 named elementary generators (Ka0,KaZ,Ks0,KsZ,KaX,KsX,KaY,KsY,
!>    GenMod.f03), each an identity matrix save for 4 diagonal cos(theta)
!>    entries and 4 off-diagonal entries equal to sin(theta) times a fixed
!>    +-1/+-i sign pattern.
!>
!>    ThetaMat encodes, for every orbital pair (p in space si, q in space
!>    sj), BOTH the angle and which of the 8 generators to use, via the
!>    *infinitesimal-generator* convention value = theta*K (K in
!>    {+-1,+-i}, |K|=1), read across up to 8 slots per pair:
!>    alpha(p,q), alpha(q,p), beta(p,q), beta(q,p) for the 4 same-spin
!>    generators (Ka0,KaZ,Ks0,KsZ -- off blocks untouched), and
!>    alphaBeta(p,q), alphaBeta(q,p), betaAlpha(p,q), betaAlpha(q,p) for
!>    the 4 spin-flip generators (KaX,KsX,KaY,KsY -- same/beta blocks
!>    untouched). The 8 slots of a pair must form an anti-Hermitian
!>    generator, X(q,p) = -conjg(X(p,q)) (for a real theta this is exactly
!>    one of the 8 templates or a combination of them); anything else is an
!>    error. Each pair is rotated exactly on its disjoint 2x2 planes: for
!>    X(a,b) = z, R(a,a) = R(b,b) = cos|z|, R(a,b) = z/|z| sin|z| and
!>    R(b,a) = -conjg(z)/|z| sin|z|, so complex entries give a unitary
!>    rotation. Same-spin planes (alpha_p,alpha_q),(beta_p,beta_q) form one
!>    factor and spin-flip planes (beta_p,alpha_q),(alpha_p,beta_q) a
!>    second one. A pair with all 8 slots at zero is skipped.
!>
!>    Elementary rotations are applied in ascending (space,orbital)
!>    order -- an arbitrary but fixed choice, since these generators do
!>    not all commute, this is NOT equivalent to rotateSolutions unless
!>    the caller's own encoding order happens to match. The second-order
!>    Zassenhaus/BCH correction rotateSolutions applies when its zaussFlag
!>    is set is NOT implemented here.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] ThetaMat
!>    \verbatim
!>        ThetaMat is type(mqc_nef_matrix)
!>        Encodes the rotation to build, per the convention above. Its
!>        own Map_Matrix_Row/Col define RotOut's structure.
!>    \endverbatim
!>
!>    \param[in] thresh
!>    \verbatim
!>        thresh is real(kind=real64),optional
!>        Numerical tolerance for "is this slot nonzero" / "do these
!>        slots agree on one theta". Default 1.0d-10.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_get_generator_rotation(ThetaMat,thresh) result(RotOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::ThetaMat
      real(kind=real64),intent(in),optional::thresh
      type(mqc_nef_matrix)::RotOut
      real(kind=real64)::myThresh
      integer::row,si,sj
!
      if(present(thresh)) then
        myThresh = thresh
      else
        myThresh = 1.0d-10
      end if
!
!     RotOut starts as the identity NEF_Matrix matching ThetaMat's own
!     structure.
!
      call RotOut%init(ThetaMat%Map_Matrix_Row,&
           initial_map_matrix_col=ThetaMat%Map_Matrix_Col)
      do row = 1, RotOut%NumSpace_row
        call RotOut%iput(mqc_integral_identity_func(&
             RotOut%Size_vector_alpha_row(row),RotOut%Size_vector_beta_row(row),&
             nAlpha2=RotOut%Size_vector_alpha_col(row),&
             nBeta2=RotOut%Size_vector_beta_col(row)),row,row)
      end do
!
      do si = 1, ThetaMat%NumSpace_row
        do sj = si, ThetaMat%NumSpace_row
          call applyPair(si,sj)
        end do
      end do
!
      contains
!
        subroutine buildPlanes(si,p,sj,q,zFirst,zSecond,spinFlip)
!       Exact exponential of the pair generator on two disjoint 2x2 planes:
!       same-spin (alpha_p,alpha_q)<-zFirst, (beta_p,beta_q)<-zSecond;
!       spin-flip (beta_p,alpha_q)<-zFirst, (alpha_p,beta_q)<-zSecond.
        implicit none
        integer,intent(in)::si,sj,p,q
        complex(kind=real64),intent(in)::zFirst,zSecond
        logical,intent(in)::spinFlip
        type(mqc_nef_matrix)::Rk
        integer::rowSeed
!
!       Rk starts as the identity, stored in full layout ('unknown' label) so
!       the (si,sj) and (sj,si) blocks written by putPlane are independent.
!
        call Rk%init(ThetaMat%Map_Matrix_Row,&
             initial_map_matrix_col=ThetaMat%Map_Matrix_Col)
        do rowSeed = 1, Rk%NumSpace_row
          call Rk%iput(mqc_integral_identity_func(&
               Rk%Size_vector_alpha_row(rowSeed),Rk%Size_vector_beta_row(rowSeed),&
               nAlpha2=Rk%Size_vector_alpha_col(rowSeed),&
               nBeta2=Rk%Size_vector_beta_col(rowSeed)),rowSeed,rowSeed)
        end do
!
        if(spinFlip) then
          call putPlane(Rk,si,sj,p,q,zFirst,'beta','alpha','alphaBeta','betaAlpha')
          call putPlane(Rk,si,sj,p,q,zSecond,'alpha','beta','betaAlpha','alphaBeta')
        else
          call putPlane(Rk,si,sj,p,q,zFirst,'alpha','alpha','alpha','alpha')
          call putPlane(Rk,si,sj,p,q,zSecond,'beta','beta','beta','beta')
        end if
!
        RotOut = matmul(RotOut,Rk)
        end subroutine buildPlanes
!
        subroutine putPlane(Rk,si,sj,p,q,z,spinP,spinQ,blockFwd,blockBwd)
!       One 2x2 plane (a,b) of exp(X) with X(a,b) = z, X(b,a) = -conjg(z).
!       Each touched block is read with iat, updated and written back with iput,
!       so si==sj (all four writes on the same diagonal block) stays consistent.
        implicit none
        type(mqc_nef_matrix),intent(inout)::Rk
        integer,intent(in)::si,sj,p,q
        complex(kind=real64),intent(in)::z
        character(len=*),intent(in)::spinP,spinQ,blockFwd,blockBwd
        type(mqc_scalar)::sElement
        real(kind=real64)::r
        complex(kind=real64)::cVal
        type(mqc_scf_integral)::tmpBlock
!
        r = abs(z)
        if(r.gt.myThresh) then
          cVal = cos(r); sElement = cVal
          tmpBlock = Rk%iat(si,si)
          call tmpBlock%put(sElement,p,p,spinP)
          call Rk%iput(tmpBlock,si,si)
          tmpBlock = Rk%iat(sj,sj)
          call tmpBlock%put(sElement,q,q,spinQ)
          call Rk%iput(tmpBlock,sj,sj)
          cVal = z/r*sin(r); sElement = cVal
          tmpBlock = Rk%iat(si,sj)
          call tmpBlock%put(sElement,p,q,blockFwd)
          call Rk%iput(tmpBlock,si,sj)
          cVal = -conjg(z)/r*sin(r); sElement = cVal
          tmpBlock = Rk%iat(sj,si)
          call tmpBlock%put(sElement,q,p,blockBwd)
          call Rk%iput(tmpBlock,sj,si)
        end if
        end subroutine putPlane
!
        subroutine applyPair(si,sj)
        implicit none
        integer,intent(in)::si,sj
        integer::p,q,qStart,sizeI,sizeJ
        type(mqc_scalar)::sElement
        complex(kind=real64),dimension(8)::V
!
        sizeI = ThetaMat%Size_vector_alpha_row(si)
        sizeJ = ThetaMat%Size_vector_alpha_row(sj)
        do p = 1, sizeI
          if(si.eq.sj) then
            qStart = p+1
          else
            qStart = 1
          end if
          do q = qStart, sizeJ
            sElement = ThetaMat%at(p,q,si,sj,'alpha');     V(1) = sElement
            sElement = ThetaMat%at(q,p,sj,si,'alpha');     V(2) = sElement
            sElement = ThetaMat%at(p,q,si,sj,'beta');      V(3) = sElement
            sElement = ThetaMat%at(q,p,sj,si,'beta');      V(4) = sElement
            sElement = ThetaMat%at(p,q,si,sj,'alphaBeta'); V(5) = sElement
            sElement = ThetaMat%at(q,p,sj,si,'alphaBeta'); V(6) = sElement
            sElement = ThetaMat%at(p,q,si,sj,'betaAlpha'); V(7) = sElement
            sElement = ThetaMat%at(q,p,sj,si,'betaAlpha'); V(8) = sElement
            if(maxval(abs(V)).le.myThresh) cycle
            if(abs(V(2)+conjg(V(1))).gt.myThresh.or.abs(V(4)+conjg(V(3))).gt.myThresh.or.&
               abs(V(8)+conjg(V(5))).gt.myThresh.or.abs(V(6)+conjg(V(7))).gt.myThresh) &
              call mqc_error_i('mqc_nef_matrix_get_generator_rotation: generator is not '//&
                 'anti-Hermitian at this orbital pair',6,&
                 'space si',si,'space sj',sj,'orbital p',p,'orbital q',q)
            if(abs(V(1)).gt.myThresh.or.abs(V(3)).gt.myThresh) &
              call buildPlanes(si,p,sj,q,V(1),V(3),.false.)
            if(abs(V(5)).gt.myThresh.or.abs(V(7)).gt.myThresh) &
              call buildPlanes(si,p,sj,q,V(5),V(7),.true.)
          end do
        end do
        end subroutine applyPair
!
      end function mqc_nef_matrix_get_generator_rotation
!
!
!     PROCEDURE MQC_NEF_Generator_Template
!
!>    \brief <b> MQC_NEF_Generator_Template returns the 8-slot +-1/+-i sign
!>    pattern for one of SNOWMAN_NEF's 8 orbital-rotation generators
!>    (Ka0,KaZ,Ks0,KsZ,KaX,KsX,KaY,KsY, GenMod.f03), shared by
!>    MQC_NEF_Matrix_Get_Rotation and MQC_NEF_Matrix_Get_Block_Rotation.
!>    </b>
!
!     Arguments:
!     ==========
!>    \param[in] genIdx
!>    \verbatim
!>        genIdx is integer
!>        1=Ka0, 2=KaZ, 3=Ks0, 4=KsZ, 5=KaX, 6=KsX, 7=KaY, 8=KsY.
!>    \endverbatim
!>
!>    \param[out] K
!>    \verbatim
!>        K is complex(kind=real64),dimension(8)
!>        Slot order: alpha(p,q), alpha(q,p), beta(p,q), beta(q,p),
!>        alphaBeta(p,q), alphaBeta(q,p), betaAlpha(p,q), betaAlpha(q,p).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_generator_template(genIdx,K)
!
      implicit none
      integer,intent(in)::genIdx
      complex(kind=real64),dimension(8),intent(out)::K
      complex(kind=real64),parameter::zero=(0.0d0,0.0d0),one=(1.0d0,0.0d0),&
           ii=(0.0d0,1.0d0)
      select case(genIdx)
      case(1)
        K = [one,-one,one,-one,zero,zero,zero,zero]
      case(2)
        K = [one,-one,-one,one,zero,zero,zero,zero]
      case(3)
        K = [ii,ii,ii,ii,zero,zero,zero,zero]
      case(4)
        K = [ii,ii,-ii,-ii,zero,zero,zero,zero]
      case(5)
        K = [zero,zero,zero,zero,one,-one,one,-one]
      case(6)
        K = [zero,zero,zero,zero,ii,ii,ii,ii]
      case(7)
        K = [zero,zero,zero,zero,ii,-ii,-ii,ii]
      case(8)
        K = [zero,zero,zero,zero,one,one,-one,-one]
      case default
        call mqc_error_i('mqc_nef_generator_template: genIdx must be 1-8 '//&
             '(Ka0/KaZ/Ks0/KsZ/KaX/KsX/KaY/KsY) -- generator 9 (Zero) has no '//&
             'single template, it decomposes into Ka0+KaZ',6,'genIdx',genIdx)
      end select
!
      end subroutine mqc_nef_generator_template
!
!
!     PROCEDURE MQC_NEF_Generator_Name_To_Index
!
!>    \brief <b> MQC_NEF_Generator_Name_To_Index maps one of the 8
!>    generator-type labels ('Ka0','KaZ','Ks0','KsZ','KaX','KsX','KaY',
!>    'KsY', case-insensitive) to its index 1-8, erroring on anything
!>    else. </b>
!
!     Arguments:
!     ==========
!>    \param[in] nameIn
!>    \verbatim
!>        nameIn is character(len=*)
!>    \endverbatim
!>
!>    \param[out] genIdx
!>    \verbatim
!>        genIdx is integer
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_generator_name_to_index(nameIn,genIdx)
!
      implicit none
      character(len=*),intent(in)::nameIn
      integer,intent(out)::genIdx
      character(len=4)::myName
!
      call string_change_case(nameIn,'L',myName)
      select case(trim(myName))
      case('ka0')
        genIdx = 1
      case('kaz')
        genIdx = 2
      case('ks0')
        genIdx = 3
      case('ksz')
        genIdx = 4
      case('kax')
        genIdx = 5
      case('ksx')
        genIdx = 6
      case('kay')
        genIdx = 7
      case('ksy')
        genIdx = 8
      case('zero')
        genIdx = 9
      case default
        call mqc_error_a('mqc_nef_generator_name_to_index: unrecognized generator '//&
             'type (expected one of Ka0/KaZ/Ks0/KsZ/KaX/KsX/KaY/KsY/Zero)',6,&
             'nameIn',trim(nameIn))
      end select
!
      end subroutine mqc_nef_generator_name_to_index
!
!
!     PROCEDURE MQC_NEF_Generator_Write_Slots
!
!>    \brief <b> MQC_NEF_Generator_Write_Slots writes mag*K(slot) into
!>    TargetMat's two cross blocks (idFwd for (p,q), idBwd for (q,p)) at
!>    whichever of the 8 generator slots are nonzero -- the off-diagonal
!>    placement logic shared by MQC_NEF_Matrix_Get_Rotation (mag=sin(theta),
!>    building an elementary rotation) and MQC_NEF_Matrix_Get_Block_Rotation
!>    (mag=theta, building a ThetaMat to feed into getRotation). </b>
!
!     Arguments:
!     ==========
!>    \param[in,out] TargetMat
!>    \verbatim
!>        TargetMat is type(mqc_nef_matrix)
!>        Main_Matrix(idFwd)/(idBwd) must already be allocated/non-zero-
!>        flagged with the correct (p,q)/(q,p) dimensions.
!>    \endverbatim
!>
!>    \param[in] idFwd, idBwd
!>    \verbatim
!>        idFwd, idBwd are integer
!>        Flat block indices for the (si,sj) and (sj,si) blocks (equal
!>        when si==sj).
!>    \endverbatim
!>
!>    \param[in] p, q
!>    \verbatim
!>        p, q are integer
!>        Local orbital indices within spaces si, sj respectively.
!>    \endverbatim
!>
!>    \param[in] mag
!>    \verbatim
!>        mag is complex(kind=real64)
!>    \endverbatim
!>
!>    \param[in] K
!>    \verbatim
!>        K is complex(kind=real64),dimension(8)
!>        One generator's sign template, from mqc_nef_generator_template.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_generator_write_slots(TargetMat,rowF,colF,p,q,mag,K)
!
      implicit none
      type(mqc_nef_matrix),intent(inout)::TargetMat
      integer,intent(in)::rowF,colF,p,q
      complex(kind=real64),intent(in)::mag
      complex(kind=real64),dimension(8),intent(in)::K
      type(mqc_scalar)::sElement
      complex(kind=real64)::cVal
      type(mqc_scf_integral)::tmpBlock
!
!     Forward block (rowF,colF): slots 1,3,5,7. Backward block (colF,rowF): slots
!     2,4,6,8. Each block is read with iat and written back with iput.
!
      tmpBlock = TargetMat%iat(rowF,colF)
      if(abs(K(1)).gt.0.5) then
        cVal = mag*K(1); sElement = cVal
        call tmpBlock%put(sElement,p,q,'alpha')
      end if
      if(abs(K(3)).gt.0.5) then
        cVal = mag*K(3); sElement = cVal
        call tmpBlock%put(sElement,p,q,'beta')
      end if
      if(abs(K(5)).gt.0.5) then
        cVal = mag*K(5); sElement = cVal
        call tmpBlock%put(sElement,p,q,'alphaBeta')
      end if
      if(abs(K(7)).gt.0.5) then
        cVal = mag*K(7); sElement = cVal
        call tmpBlock%put(sElement,p,q,'betaAlpha')
      end if
      call TargetMat%iput(tmpBlock,rowF,colF)
      tmpBlock = TargetMat%iat(colF,rowF)
      if(abs(K(2)).gt.0.5) then
        cVal = mag*K(2); sElement = cVal
        call tmpBlock%put(sElement,q,p,'alpha')
      end if
      if(abs(K(4)).gt.0.5) then
        cVal = mag*K(4); sElement = cVal
        call tmpBlock%put(sElement,q,p,'beta')
      end if
      if(abs(K(6)).gt.0.5) then
        cVal = mag*K(6); sElement = cVal
        call tmpBlock%put(sElement,q,p,'alphaBeta')
      end if
      if(abs(K(8)).gt.0.5) then
        cVal = mag*K(8); sElement = cVal
        call tmpBlock%put(sElement,q,p,'betaAlpha')
      end if
      call TargetMat%iput(tmpBlock,colF,rowF)
!
      end subroutine mqc_nef_generator_write_slots
!
!
!     PROCEDURE MQC_NEF_Generator_Write_AntiHermitian
!
!>    \brief <b> Writes one orbital pair of an anti-Hermitian generator into a
!>    ThetaMat: forward slots get mag*K, each backward slot gets minus the complex
!>    conjugate of its forward partner. Identical to MQC_NEF_Generator_Write_Slots
!>    for a real mag; for a complex mag the generator stays anti-Hermitian, so the
!>    rotation built from it is unitary. </b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_generator_write_antihermitian(TargetMat,rowF,colF,p,q,mag,K)
!
      implicit none
      type(mqc_nef_matrix),intent(inout)::TargetMat
      integer,intent(in)::rowF,colF,p,q
      complex(kind=real64),intent(in)::mag
      complex(kind=real64),dimension(8),intent(in)::K
      complex(kind=real64),dimension(8)::V
      type(mqc_scalar)::sElement
      type(mqc_scf_integral)::tmpBlock
!
      V = mag*K
      V(2) = -conjg(V(1))
      V(4) = -conjg(V(3))
      V(8) = -conjg(V(5))
      V(6) = -conjg(V(7))
!
!     Forward block (rowF,colF) first, then backward block (colF,rowF). When
!     rowF=colF both are the same block; the second iat then already holds the
!     forward writes, and the (p,q)/(q,p) elements never overlap (q>p).
!
      tmpBlock = TargetMat%iat(rowF,colF)
      if(abs(K(1)).gt.0.5) then
        sElement = V(1); call tmpBlock%put(sElement,p,q,'alpha')
      end if
      if(abs(K(3)).gt.0.5) then
        sElement = V(3); call tmpBlock%put(sElement,p,q,'beta')
      end if
      if(abs(K(5)).gt.0.5) then
        sElement = V(5); call tmpBlock%put(sElement,p,q,'alphaBeta')
      end if
      if(abs(K(7)).gt.0.5) then
        sElement = V(7); call tmpBlock%put(sElement,p,q,'betaAlpha')
      end if
      call TargetMat%iput(tmpBlock,rowF,colF)
      tmpBlock = TargetMat%iat(colF,rowF)
      if(abs(K(1)).gt.0.5) then
        sElement = V(2); call tmpBlock%put(sElement,q,p,'alpha')
      end if
      if(abs(K(3)).gt.0.5) then
        sElement = V(4); call tmpBlock%put(sElement,q,p,'beta')
      end if
      if(abs(K(5)).gt.0.5) then
        sElement = V(8); call tmpBlock%put(sElement,q,p,'betaAlpha')
      end if
      if(abs(K(7)).gt.0.5) then
        sElement = V(6); call tmpBlock%put(sElement,q,p,'alphaBeta')
      end if
      call TargetMat%iput(tmpBlock,colF,rowF)
!
      end subroutine mqc_nef_generator_write_antihermitian
!
!
!     PROCEDURE MQC_NEF_Matrix_Get_Block_Rotation
!
!>    \brief <b> MQC_NEF_Matrix_Get_Block_Rotation builds the rotation
!>    corresponding to one off-diagonal block of orbital-rotation angles
!>    (InputNEF's block at InputPair), of a caller-specified generator
!>    type, placed at a given position (TemplatePair) within a larger
!>    output structure (TemplateMat). </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    A single off-diagonal block of values cannot, by itself, say which
!>    of the 8 generators (Ka0,KaZ,Ks0,KsZ,KaX,KsX,KaY,KsY) it represents
!>    -- GeneratorType supplies that missing information. Depending on
!>    GeneratorType, the block is read from InputNEF%Main_Matrix(InputPair)
!>    via %getBlock('alpha') (the 4 same-spin generators) or
!>    %getBlock('alphaBeta') (the 4 spin-flip generators). Each entry
!>    (p,q) of that block is a (generally complex) rotation angle theta;
!>    a fresh ThetaMat, shaped like TemplateMat and zero everywhere except
!>    the two blocks at TemplatePair and its reverse, is populated at
!>    every (p,q) via mqc_nef_generator_write_antihermitian(...,mag=theta,...),
!>    then handed to %getRotation() to build the actual rotation. When
!>    TemplatePair's two space-indices are equal (intra-space), only the
!>    upper triangle (p<q) of the input block is used, matching
!>    getRotation's own intra-space convention.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] InputNEF
!>    \verbatim
!>        InputNEF is type(mqc_nef_matrix)
!>        Holds the (generally complex) rotation-angle values.
!>    \endverbatim
!>
!>    \param[in] TemplateMat
!>    \verbatim
!>        TemplateMat is type(mqc_nef_matrix)
!>        Defines RotOut's full structure (Map_Matrix_Row/Col).
!>    \endverbatim
!>
!>    \param[in] InputPair
!>    \verbatim
!>        InputPair is integer,dimension(2)
!>        The space-index pair of InputNEF holding the block to use.
!>    \endverbatim
!>
!>    \param[in] TemplatePair
!>    \verbatim
!>        TemplatePair is integer,dimension(2)
!>        Where that block's space-indices map to within TemplateMat.
!>    \endverbatim
!>
!>    \param[in] GeneratorType
!>    \verbatim
!>        GeneratorType is character(len=*)
!>        One of 'Ka0','KaZ','Ks0','KsZ','KaX','KsX','KaY','KsY'
!>        (case-insensitive).
!>    \endverbatim
!>
!>    \param[in] thresh
!>    \verbatim
!>        thresh is real(kind=real64),optional
!>        Forwarded to %getRotation(). Default 1.0d-10.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_get_block_rotation(InputNEF,TemplateMap,InputPair,&
           TemplatePair,GeneratorType,thresh) result(RotOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::InputNEF
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,intent(in)::TemplateMap
      integer,dimension(2),intent(in)::InputPair,TemplatePair
      character(len=*),intent(in)::GeneratorType
      real(kind=real64),intent(in),optional::thresh
      type(mqc_nef_matrix)::RotOut
      type(mqc_nef_matrix)::ThetaMat
      type(mqc_matrix)::baseBlock,baseBlockAlpha,baseBlockBeta
      complex(kind=real64),dimension(8)::K,KAlphaOnly,KBetaOnly
      type(mqc_scalar)::sElement
      complex(kind=real64)::theta,thetaAlpha,thetaBeta
      real(kind=real64)::myThresh
      integer::genIdx,p,q,qStart
      type(mqc_scf_integral)::inBlock
      integer::nRowIn,nColIn,nRowTemplate,nColTemplate
!
      if(present(thresh)) then
        myThresh = thresh
      else
        myThresh = 1.0d-10
      end if
!
      call mqc_nef_generator_name_to_index(GeneratorType,genIdx)
      if(genIdx.ne.9) call mqc_nef_generator_template(genIdx,K)
!
      call ThetaMat%init(templateMap)
!
      if(genIdx.le.4.or.genIdx.eq.9) then
        nRowIn = InputNEF%Size_vector_alpha_row(InputPair(1))
        nColIn = InputNEF%Size_vector_alpha_col(InputPair(2))
      else
        nRowIn = InputNEF%Size_vector_beta_row(InputPair(1))
        nColIn = InputNEF%Size_vector_alpha_col(InputPair(2))
      end if
!
      nRowTemplate = ThetaMat%Size_vector_alpha_row(TemplatePair(1))
      nColTemplate = ThetaMat%Size_vector_alpha_col(TemplatePair(2))
      if(nRowIn.ne.nRowTemplate.or.nColIn.ne.nColTemplate) call mqc_error_i(&
        'mqc_nef_matrix_get_block_rotation: input block dimensions do not '//&
        'match the template dimensions at templatePair',6,&
        'nRowIn',nRowIn,'nRowTemplate',nRowTemplate,&
        'nColIn',nColIn,'nColTemplate',nColTemplate)
!
      if(.not.InputNEF%isBlockZero(InputPair(1),InputPair(2))) then
        inBlock = InputNEF%iat(InputPair(1),InputPair(2))
        if(genIdx.eq.9) then
          KAlphaOnly = [(1.0d0,0.0d0),(-1.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),&
               (0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)]
          KBetaOnly = [(0.0d0,0.0d0),(0.0d0,0.0d0),(1.0d0,0.0d0),(-1.0d0,0.0d0),&
               (0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0),(0.0d0,0.0d0)]
          baseBlockAlpha = inBlock%getBlock('alpha')
          baseBlockBeta = inBlock%getBlock('beta')
          do p = 1, nRowIn
            if(TemplatePair(1).eq.TemplatePair(2)) then
              qStart = p+1
            else
              qStart = 1
            end if
            do q = qStart, nColIn
              sElement = baseBlockAlpha%at(p,q); thetaAlpha = sElement
              sElement = baseBlockBeta%at(p,q); thetaBeta = sElement
              if(abs(thetaAlpha).gt.myThresh) &
                call mqc_nef_generator_write_antihermitian(ThetaMat,TemplatePair(1),TemplatePair(2),p,q,thetaAlpha,KAlphaOnly)
              if(abs(thetaBeta).gt.myThresh) &
                call mqc_nef_generator_write_antihermitian(ThetaMat,TemplatePair(1),TemplatePair(2),p,q,thetaBeta,KBetaOnly)
            end do
          end do
        else
          if(genIdx.le.4) then
            baseBlock = inBlock%getBlock('alpha')
          else
            baseBlock = inBlock%getBlock('alphaBeta')
          end if
          do p = 1, nRowIn
            if(TemplatePair(1).eq.TemplatePair(2)) then
              qStart = p+1
            else
              qStart = 1
            end if
            do q = qStart, nColIn
              sElement = baseBlock%at(p,q)
              theta = sElement
              if(abs(theta).le.myThresh) cycle
              call mqc_nef_generator_write_antihermitian(ThetaMat,TemplatePair(1),TemplatePair(2),p,q,theta,K)
            end do
          end do
        end if
      end if
!
      RotOut = ThetaMat%getRotation(thresh=myThresh)
!
      end function mqc_nef_matrix_get_block_rotation
!
!
!     PROCEDURE MQC_NEF_Matrix_Exp
!
!>    \brief <b> MQC_NEF_Matrix_Exp returns the matrix exponential exp(A) of a
!>    square NEF Matrix</b>
!>    \par Purpose:
!     =============
!>    \verbatim
!>    Scaling and squaring with a Taylor series. A is scaled by 2**(-s) so its
!>    Frobenius norm is at most 0.5, the series is summed until the Frobenius
!>    norm of the last term is below thresh, and the sum is squared s times.
!>    Works on the NEF block structure directly (no export); the row and column
!>    maps of A must describe the same spaces. For an anti-Hermitian A
!>    (A = i*t*H, H Hermitian) the result is unitary to roundoff.
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMat
!>    \verbatim
!>        NEFMat is Class(MQC_NEF_Matrix)
!>        The square NEF Matrix A.
!>    \endverbatim
!>
!>    \param[in] thresh
!>    \verbatim
!>        thresh is real(kind=real64),optional
!>        Convergence threshold on the Frobenius norm of the last Taylor term.
!>        Default 1.0e-16.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_exp(NEFMat,thresh) result(ExpOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      real(kind=real64),intent(in),optional::thresh
      type(mqc_nef_matrix)::ExpOut
      type(mqc_nef_matrix)::scaledMat,term
      type(mqc_scalar)::factor
      real(kind=real64)::myThresh,matNorm
      integer::nSquare,k,row
      integer,parameter::maxTerms=40
!
      if(present(thresh)) then
        myThresh = thresh
      else
        myThresh = 1.0e-16
      end if
!
      if(NEFMat%NumSpace_row.ne.NEFMat%NumSpace_col) call mqc_error_i(&
        'mqc_nef_matrix_exp: NEF Matrix is not square',6,&
        'NumSpace_row',NEFMat%NumSpace_row,'NumSpace_col',NEFMat%NumSpace_col)
      do row = 1, NEFMat%NumSpace_row
        if(NEFMat%Size_vector_alpha_row(row).ne.NEFMat%Size_vector_alpha_col(row).or.&
           NEFMat%Size_vector_beta_row(row).ne.NEFMat%Size_vector_beta_col(row)) &
          call mqc_error_i('mqc_nef_matrix_exp: NEF Matrix is not square',6,'space',row)
      end do
!
      matNorm = mqc_nef_matrix_frobenius(NEFMat)
      nSquare = 0
      if(matNorm.gt.0.5) nSquare = ceiling(log(matNorm/0.5)/log(2.0))
      factor = 1.0/(2.0**nSquare)
      scaledMat = NEFMat*factor
!
!     exp(X) is an operator product: label 'asymmetric' (also for exp(0) = 1)
      ExpOut = mqc_nef_matrix_identity(NEFMat,symmetry='asymmetric')
      term = ExpOut
      do k = 1, maxTerms
        factor = 1.0/real(k,kind=real64)
        term = matmul(term,scaledMat)*factor
        ExpOut = ExpOut + term
        if(mqc_nef_matrix_frobenius(term).le.myThresh) exit
      end do
      do k = 1, nSquare
        ExpOut = matmul(ExpOut,ExpOut)
      end do
!
      end function mqc_nef_matrix_exp
!
!
!     PROCEDURE MQC_NEF_Matrix_Identity
!
!>    \brief <b> Identity NEF Matrix with the row/column structure of a square
!>    template NEF Matrix</b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_identity(Template,symmetry) result(IdentOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::Template
      character(len=*),intent(in),optional::symmetry
      character(len=16)::mySymm
      type(mqc_nef_matrix)::IdentOut
      integer::row
!
!     The identity operator is hermitian (default label); the bra-ket flag follows
!     the template. A caller building another operator from it (exp) passes its label.
      mySymm = 'hermitian'
      if(present(symmetry)) mySymm = symmetry
      call IdentOut%init(Template%Map_Matrix_Row,&
           initial_map_matrix_col=Template%Map_Matrix_Col,&
           symmetry=mySymm,sameBraKet=Template%sameBraKet)
      do row = 1, IdentOut%NumSpace_row
        call IdentOut%iput(mqc_integral_identity_func(&
             IdentOut%Size_vector_alpha_row(row),IdentOut%Size_vector_beta_row(row),&
             nAlpha2=IdentOut%Size_vector_alpha_col(row),&
             nBeta2=IdentOut%Size_vector_beta_col(row)),row,row)
      end do
!
      end function mqc_nef_matrix_identity
!
!
!     PROCEDURE MQC_NEF_Matrix_Frobenius
!
!>    \brief <b> Frobenius norm of a NEF Matrix (non-zero blocks only)</b>
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_matrix_frobenius(NEFMat) result(normOut)
!
      implicit none
      class(mqc_nef_matrix),intent(in)::NEFMat
      real(kind=real64)::normOut
      type(mqc_scalar)::blockNorm
      integer::row,col
      type(mqc_scf_integral)::tmpBlock
      normOut = 0.0
      do row = 1, NEFMat%NumSpace_row
        do col = 1, NEFMat%NumSpace_col
          if(NEFMat%isBlockZero(row,col)) cycle
          tmpBlock = NEFMat%iat(row,col)
          blockNorm = tmpBlock%norm('F')
          normOut = normOut + abs(blockNorm%cval())**2
        end do
      end do
      normOut = sqrt(normOut)
!
      end function mqc_nef_matrix_frobenius
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
      class(mqc_nef_r4),intent(in),target::NEFR4
      integer(kind=int64),intent(in)::dim1,dim2
      type(mqc_nef_r4)::outR4
      integer(kind=int64)::mydim1,mydim2,id1,id2,id3,id4
      type(mqc_twoeris),pointer::blk
!
      type(mqc_nef_r4),target::full_NEFR4
      class(mqc_nef_r4),pointer::ref_NEFR4
!     Mirror slots of a symmetric input are read from a full-storage copy (iref
!     holds no mirror block); full-storage inputs are read in place.
      ref_NEFR4 => NEFR4
      if(NEFR4%Mirrored) then
        full_NEFR4 = mqc_nef_r4_full_storage(NEFR4)
        ref_NEFR4 => full_NEFR4
      end if
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
          call outR4%init(NEFR4%Map_Matrix_2,NEFR4%Map_Matrix_1,&
               NEFR4%Map_Matrix_3,NEFR4%Map_Matrix_4,&
               symmetry='asymmetric',sameBraKet=NEFR4%sameBraKet)

          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  if (.not.NEFR4%isBlockZero(id1,id2,id3,id4)) then
                    blk => ref_NEFR4%iref(id1,id2,id3,id4)
                    call outR4%iput(blk%swap(mydim1,mydim2),id2,id1,id3,id4)
                  end if
                end do
              end do
            end do
          end do
        case (3)
          call outR4%init(NEFR4%Map_Matrix_3,NEFR4%Map_Matrix_2,&
               NEFR4%Map_Matrix_1,NEFR4%Map_Matrix_4,&
               symmetry='asymmetric',sameBraKet=NEFR4%sameBraKet)

          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  if (.not.NEFR4%isBlockZero(id1,id2,id3,id4)) then
                    blk => ref_NEFR4%iref(id1,id2,id3,id4)
                    call outR4%iput(blk%swap(mydim1,mydim2),id3,id2,id1,id4)
                  end if
                end do
              end do
            end do
          end do
        case (4)
          call outR4%init(NEFR4%Map_Matrix_4,NEFR4%Map_Matrix_2,&
               NEFR4%Map_Matrix_3,NEFR4%Map_Matrix_1,&
               symmetry='asymmetric',sameBraKet=NEFR4%sameBraKet)

          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  if (.not.NEFR4%isBlockZero(id1,id2,id3,id4)) then
                    blk => ref_NEFR4%iref(id1,id2,id3,id4)
                    call outR4%iput(blk%swap(mydim1,mydim2),id4,id2,id3,id1)
                  end if
                end do
              end do
            end do
          end do
        end select
      case (2)
        select case (mydim2)
        case (3)
          call outR4%init(NEFR4%Map_Matrix_1,NEFR4%Map_Matrix_3,&
               NEFR4%Map_Matrix_2,NEFR4%Map_Matrix_4,&
               symmetry='asymmetric',sameBraKet=NEFR4%sameBraKet)

          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  if (.not.NEFR4%isBlockZero(id1,id2,id3,id4)) then
                    blk => ref_NEFR4%iref(id1,id2,id3,id4)
                    call outR4%iput(blk%swap(mydim1,mydim2),id1,id3,id2,id4)
                  end if
                end do
              end do
            end do
          end do
        case (4)
          call outR4%init(NEFR4%Map_Matrix_1,NEFR4%Map_Matrix_4,&
               NEFR4%Map_Matrix_3,NEFR4%Map_Matrix_2,&
               symmetry='asymmetric',sameBraKet=NEFR4%sameBraKet)

          do id1 = 1, NEFR4%NumSpace_1
            do id2 = 1, NEFR4%NumSpace_2
              do id3 = 1, NEFR4%NumSpace_3
                do id4 = 1, NEFR4%NumSpace_4
                  if (.not.NEFR4%isBlockZero(id1,id2,id3,id4)) then
                    blk => ref_NEFR4%iref(id1,id2,id3,id4)
                    call outR4%iput(blk%swap(mydim1,mydim2),id1,id4,id3,id2)
                  end if
                end do
              end do
            end do
          end do
        end select
      case (3)
        call outR4%init(NEFR4%Map_Matrix_1,NEFR4%Map_Matrix_2,&
             NEFR4%Map_Matrix_4,NEFR4%Map_Matrix_3,&
             symmetry='asymmetric',sameBraKet=NEFR4%sameBraKet)

        do id1 = 1, NEFR4%NumSpace_1
          do id2 = 1, NEFR4%NumSpace_2
            do id3 = 1, NEFR4%NumSpace_3
              do id4 = 1, NEFR4%NumSpace_4
                if (.not.NEFR4%isBlockZero(id1,id2,id3,id4)) then
                  blk => ref_NEFR4%iref(id1,id2,id3,id4)
                  call outR4%iput(blk%swap(mydim1,mydim2),id1,id2,id4,id3)
                end if
              end do
            end do
          end do
        end do
      end select
!
      end function mqc_nef_r4_dimension_swap
!
!
!     PROCEDURE MQC_NEF_R4_Full_Transform
!
!>    \brief <b> MQC_NEF_R4_Full_Transform transforms a MQC NEF R4
!>    two-electron integral set to a new basis using four independent
!>    MQC NEF Matrix coefficients, one per axis, each of which may freely
!>    mix NEF spaces (e.g. a getRotation()-style coefficient). </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    Generalizes TwoERI_Trans_Full (mqc_est.F03) -- which applies four
!>    independent, dense MO coefficient matrices to a flat mqc_twoeris --
!>    to the NEF level, in two ways at once: (1) the four coefficients are
!>    MQC_NEF_Matrix objects that may have nonzero blocks between ANY pair
!>    of NEF spaces, not just the diagonal, so an output space can receive
!>    contributions mixed in from multiple input spaces; and (2) every
!>    nonzero block of NEFR4A and of the four coefficients may be of ANY
!>    type (regular/space/spin/general for NEFR4A; space/spin/general for
!>    the coefficients).
!>
!>    This is achieved by working entirely in the merged ('full',
!>    alpha+beta-combined) representation throughout the transform itself
!>    (via the existing, already-type-agnostic %getBlock('full') on both
!>    the integral and the coefficient blocks), and only splitting the
!>    result back into named alpha/beta/alphaBeta/betaAlpha/... blocks at
!>    the very end, once the output type (space/spin/general) has been
!>    resolved. Implemented as two sequential stages (axes 3,4 then axes
!>    1,2), mirroring TwoERI_Trans_Full's own split, with each stage
!>    summing over every NEF-space combination that is not structurally
!>    zero-flagged.
!>
!>    Output is 'general'-type (16 blocks per space-quadruple) if any
!>    nonzero block of NEFR4A or of the four coefficients is 'general'-
!>    type; else 'spin'-type (4 blocks) if any is 'spin'-type; else
!>    'space'-type (1 block). storage is always 'full'.
!>
!>    Every output space-quadruple is built, even when every contribution
!>    to it is structurally zero -- this does not yet exploit the NEF-level
!>    sparsity a more optimized version could. Known, deferred performance
!>    limitation, not a correctness one.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFR4A
!>    \verbatim
!>        NEFR4A is type(mqc_nef_r4)
!>        The two-electron integrals to transform. Any block type.
!>    \endverbatim
!>
!>    \param[in] NEFCoeff_1
!>    \verbatim
!>        NEFCoeff_1 is type(mqc_nef_matrix)
!>        Coefficients driving axis 1 (P). Required.
!>    \endverbatim
!>
!>    \param[in] NEFCoeff_2
!>    \verbatim
!>        NEFCoeff_2 is type(mqc_nef_matrix),optional
!>        Coefficients driving axis 2 (Q). Defaults to NEFCoeff_1.
!>    \endverbatim
!>
!>    \param[in] NEFCoeff_3
!>    \verbatim
!>        NEFCoeff_3 is type(mqc_nef_matrix),optional
!>        Coefficients driving axis 3 (R). Defaults to NEFCoeff_1.
!>    \endverbatim
!>
!>    \param[in] NEFCoeff_4
!>    \verbatim
!>        NEFCoeff_4 is type(mqc_nef_matrix),optional
!>        Coefficients driving axis 4 (S). Defaults to NEFCoeff_1.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_r4_full_transform(NEFR4A,NEFCoeff_1,NEFCoeff_2,&
           NEFCoeff_3,NEFCoeff_4) result(NEFR4Out)
!
      implicit none
      class(mqc_nef_r4),intent(in),target::NEFR4A
      type(mqc_nef_matrix),intent(in)::NEFCoeff_1
      type(mqc_nef_matrix),optional,intent(in)::NEFCoeff_2,NEFCoeff_3,NEFCoeff_4
      logical::pairedCoeff
      type(mqc_nef_r4)::NEFR4Out
      type(mqc_nef_r4),target::full_NEFR4A
      class(mqc_nef_r4),pointer::ref_NEFR4A
      type(mqc_scf_eigenvalues),dimension(:),allocatable::outMap1,outMap2,outMap3,outMap4
      type(mqc_nef_matrix)::myCoeff_2,myCoeff_3,myCoeff_4
      type(mqc_nef_r4)::stageFull,finalFull
      character(len=15)::storageType,interactionType
      logical::hasSpin,hasGeneral,anySpin,anyGeneral
      integer(kind=int64)::P,Q,R,S
      integer(kind=int64)::nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4
      type(mqc_twoeris)::blkOut,blkEmpty
      type(mqc_r4tensor)::fullTensor
      type(mqc_r4tensor)::tAlpha,tBeta,tAlphaBeta,tBetaAlpha
      type(mqc_r4tensor)::tABAB,tABBA,tBAAB,tBABA,tAAAB,tAABA,tABAA,tBAAA
      type(mqc_r4tensor)::tBBBA,tBBAB,tBABB,tABBB
!
      if(present(NEFCoeff_2)) then
        myCoeff_2 = NEFCoeff_2
      else
        myCoeff_2 = NEFCoeff_1
      end if
      if(present(NEFCoeff_3)) then
        myCoeff_3 = NEFCoeff_3
      else
        myCoeff_3 = NEFCoeff_1
      end if
      if(present(NEFCoeff_4)) then
        myCoeff_4 = NEFCoeff_4
      else
        myCoeff_4 = NEFCoeff_1
      end if
!
      anySpin = .false.
      anyGeneral = .false.
!     Mirror slots of a symmetric input are read from a full-storage copy.
      ref_NEFR4A => NEFR4A
      if(NEFR4A%Mirrored) then
        full_NEFR4A = mqc_nef_r4_full_storage(NEFR4A)
        ref_NEFR4A => full_NEFR4A
      end if
      call detectR4InteractionType(ref_NEFR4A,interactionType)
      call detectR4SpinGeneral(ref_NEFR4A,hasSpin,hasGeneral)
      anySpin = anySpin.or.hasSpin; anyGeneral = anyGeneral.or.hasGeneral
      call detectCoeffSpinGeneral(NEFCoeff_1,hasSpin,hasGeneral)
      anySpin = anySpin.or.hasSpin; anyGeneral = anyGeneral.or.hasGeneral
      call detectCoeffSpinGeneral(myCoeff_2,hasSpin,hasGeneral)
      anySpin = anySpin.or.hasSpin; anyGeneral = anyGeneral.or.hasGeneral
      call detectCoeffSpinGeneral(myCoeff_3,hasSpin,hasGeneral)
      anySpin = anySpin.or.hasSpin; anyGeneral = anyGeneral.or.hasGeneral
      call detectCoeffSpinGeneral(myCoeff_4,hasSpin,hasGeneral)
      anySpin = anySpin.or.hasSpin; anyGeneral = anyGeneral.or.hasGeneral
!
      if(anyGeneral) then
        storageType = 'general'
      else if(anySpin) then
        storageType = 'spin'
      else
        storageType = 'space'
      end if
!
      stageFull = buildStage(ref_NEFR4A,3,4,myCoeff_3,myCoeff_4)
      finalFull = buildStage(stageFull,1,2,NEFCoeff_1,myCoeff_2)
!
!     Final assembly: wrap each output space-quadruple's merged 'full'
!     result into the appropriately-typed mqc_twoeris block.
!
      outMap1 = NEFCoeff_1%getMap(2)
      outMap2 = myCoeff_2%getMap(2)
      outMap3 = myCoeff_3%getMap(2)
      outMap4 = myCoeff_4%getMap(2)
!     The orbital transform keeps the bra-ket flag (O8). It keeps the operator label
!     when the dagger-paired axes use the same coefficients, C1 = C4 and C2 = C3 (T1, user),
!     compared by content, real or complex: the two passes apply conjg(C1), C2, conjg(C3), C4,
!     which keeps R(p,q,r,s) = +/- conj R(s,r,q,p) (B.5, brute-force check).
      pairedCoeff = .true.
      if(present(NEFCoeff_2).or.present(NEFCoeff_3).or.present(NEFCoeff_4)) &
        pairedCoeff = mqc_nef_matrix_equal(NEFCoeff_1,myCoeff_4).and.&
                      mqc_nef_matrix_equal(myCoeff_2,myCoeff_3)
      call NEFR4Out%init(outMap1,outMap2,outMap3,outMap4,&
           symmetry=merge(NEFR4A%Symmetry,'asymmetric      ',pairedCoeff),&
           sameBraKet=NEFR4A%sameBraKet)

!
      do P = 1, NEFR4Out%NumSpace_1
        do Q = 1, NEFR4Out%NumSpace_2
          do R = 1, NEFR4Out%NumSpace_3
            do S = 1, NEFR4Out%NumSpace_4
!             P2: a mirrored output stores one block per dagger pair
              if(NEFR4Out%isMirrorSlot(int(P),int(Q),int(R),int(S))) cycle
              blkOut = blkEmpty
              nA1 = NEFCoeff_1%Size_vector_alpha_col(P); nB1 = NEFCoeff_1%Size_vector_beta_col(P)
              nA2 = myCoeff_2%Size_vector_alpha_col(Q);  nB2 = myCoeff_2%Size_vector_beta_col(Q)
              nA3 = myCoeff_3%Size_vector_alpha_col(R);  nB3 = myCoeff_3%Size_vector_beta_col(R)
              nA4 = myCoeff_4%Size_vector_alpha_col(S);  nB4 = myCoeff_4%Size_vector_beta_col(S)
              call fetchFull(finalFull,P,Q,R,S,nA1+nB1,nA2+nB2,nA3+nB3,nA4+nB4,fullTensor)
              tAlpha = sliceCorner(fullTensor,'a','a','a','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
              select case(storageType)
              case('space')
                call mqc_twoeris_allocate(blkOut,'full-'//trim(interactionType),'space',tAlpha)
              case('spin')
                tBeta = sliceCorner(fullTensor,'b','b','b','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                select case(trim(interactionType))
                case('coulomb')
                  tAlphaBeta = sliceCorner(fullTensor,'a','a','b','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tBetaAlpha = sliceCorner(fullTensor,'b','b','a','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  call mqc_twoeris_allocate(blkOut,'full-coulomb','spin',&
                       tAlpha,tBeta,tAlphaBeta,tBetaAlpha)
                case('exchange')
                  tAlphaBeta = sliceCorner(fullTensor,'a','b','b','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tBetaAlpha = sliceCorner(fullTensor,'b','a','a','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  call mqc_twoeris_allocate(blkOut,'full-exchange','spin',&
                       tAlpha,tBeta,tAlphaBeta,tBetaAlpha)
                case('inverted')
                  tAlphaBeta = sliceCorner(fullTensor,'a','b','a','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tBetaAlpha = sliceCorner(fullTensor,'b','a','b','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  call mqc_twoeris_allocate(blkOut,'full-inverted','spin',&
                       tAlpha,tBeta,tAlphaBeta,tBetaAlpha)
                case('doublebar')
                  tAlphaBeta = sliceCorner(fullTensor,'a','a','b','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tBetaAlpha = sliceCorner(fullTensor,'b','b','a','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tABBA = sliceCorner(fullTensor,'a','b','b','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tBAAB = sliceCorner(fullTensor,'b','a','a','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  call mqc_twoeris_allocate(blkOut,'full-doublebar','spin',&
                       tAlpha,tBeta,tAlphaBeta,tBetaAlpha,abba=tABBA,baab=tBAAB)
                case('ms')
                  tAlphaBeta = sliceCorner(fullTensor,'a','a','b','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tBetaAlpha = sliceCorner(fullTensor,'b','b','a','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tABAB = sliceCorner(fullTensor,'a','b','a','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tABBA = sliceCorner(fullTensor,'a','b','b','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tBAAB = sliceCorner(fullTensor,'b','a','a','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  tBABA = sliceCorner(fullTensor,'b','a','b','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                  call mqc_twoeris_allocate(blkOut,'full-ms','spin',&
                       tAlpha,tBeta,tAlphaBeta,tBetaAlpha,tABAB,tABBA,tBAAB,tBABA)
                case default
                  call mqc_error_a('mqc_nef_r4_full_transform: unsupported interaction type',6,&
                       'interactionType',trim(interactionType))
                end select
              case('general')
                tBeta = sliceCorner(fullTensor,'b','b','b','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tAlphaBeta = sliceCorner(fullTensor,'a','a','b','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tBetaAlpha = sliceCorner(fullTensor,'b','b','a','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tABAB = sliceCorner(fullTensor,'a','b','a','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tABBA = sliceCorner(fullTensor,'a','b','b','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tBAAB = sliceCorner(fullTensor,'b','a','a','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tBABA = sliceCorner(fullTensor,'b','a','b','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tAAAB = sliceCorner(fullTensor,'a','a','a','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tAABA = sliceCorner(fullTensor,'a','a','b','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tABAA = sliceCorner(fullTensor,'a','b','a','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tBAAA = sliceCorner(fullTensor,'b','a','a','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tBBBA = sliceCorner(fullTensor,'b','b','b','a',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tBBAB = sliceCorner(fullTensor,'b','b','a','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tBABB = sliceCorner(fullTensor,'b','a','b','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                tABBB = sliceCorner(fullTensor,'a','b','b','b',nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4)
                call mqc_twoeris_allocate(blkOut,'full-'//trim(interactionType),'general',&
                     tAlpha,tBeta,tAlphaBeta,tBetaAlpha,tABAB,tABBA,tBAAB,tBABA,&
                     tAAAB,tAABA,tABAA,tBAAA,tBBBA,tBBAB,tBABB,tABBB)
              end select
              call NEFR4Out%iput(blkOut,P,Q,R,S)
            end do
          end do
        end do
      end do
!
      Contains
!
        subroutine detectR4InteractionType(R4In,interactionTypeOut)
        implicit none
        class(mqc_nef_r4),intent(in),target::R4In
        character(len=*),intent(out)::interactionTypeOut
        character(len=64)::blockInteractionType
        integer(kind=int64)::i1,i2,i3,i4
        type(mqc_twoeris),pointer::block
        interactionTypeOut = ''
        do i1=1,R4In%blockSize(1)
          do i2=1,R4In%blockSize(2)
            do i3=1,R4In%blockSize(3)
              do i4=1,R4In%blockSize(4)
                if(R4In%isBlockZero(i1,i2,i3,i4)) cycle
                block => R4In%iref(i1,i2,i3,i4)
                blockInteractionType = trim(block%intType())
                if(len_trim(interactionTypeOut).eq.0) then
                  interactionTypeOut = trim(blockInteractionType)
                else if(trim(interactionTypeOut).ne.trim(blockInteractionType)) then
                  call mqc_error('mqc_nef_r4_full_transform: source blocks have mixed interaction types',6)
                end if
              end do
            end do
          end do
        end do
        if(len_trim(interactionTypeOut).eq.0) interactionTypeOut = 'coulomb'
        end subroutine detectR4InteractionType
!
        subroutine detectR4SpinGeneral(R4In,hasSpinOut,hasGeneralOut)
        implicit none
        type(mqc_nef_r4),intent(in),target::R4In
        logical,intent(out)::hasSpinOut,hasGeneralOut
        integer(kind=int64)::i1,i2,i3,i4
        type(mqc_twoeris),pointer::blkIn
        hasSpinOut = .false.
        hasGeneralOut = .false.
        do i1=1,R4In%blockSize(1)
          do i2=1,R4In%blockSize(2)
            do i3=1,R4In%blockSize(3)
              do i4=1,R4In%blockSize(4)
                if(R4In%isBlockZero(i1,i2,i3,i4)) cycle
                blkIn => R4In%iref(i1,i2,i3,i4)
                if(blkIn%type().eq.'general') hasGeneralOut = .true.
                if(blkIn%type().eq.'spin') hasSpinOut = .true.
              end do
            end do
          end do
        end do
        end subroutine detectR4SpinGeneral
!
        subroutine detectCoeffSpinGeneral(CoeffIn,hasSpinOut,hasGeneralOut)
        implicit none
        type(mqc_nef_matrix),intent(in)::CoeffIn
        logical,intent(out)::hasSpinOut,hasGeneralOut
        integer(kind=int64)::i1,i2
        type(mqc_scf_integral)::blkCoeff
        hasSpinOut = .false.
        hasGeneralOut = .false.
        do i1=1,CoeffIn%NumSpace_row
          do i2=1,CoeffIn%NumSpace_col
            if(CoeffIn%isBlockZero(i1,i2)) cycle
            blkCoeff = CoeffIn%iat(i1,i2)
            if(blkCoeff%type().eq.'general') hasGeneralOut = .true.
            if(blkCoeff%type().eq.'spin') hasSpinOut = .true.
          end do
        end do
        end subroutine detectCoeffSpinGeneral
!
        function buildStage(srcR4,axisRow,axisCol,CoeffRow,CoeffCol) result(outR4)
        implicit none
        type(mqc_nef_r4),intent(in),target::srcR4
        integer(kind=int64),intent(in)::axisRow,axisCol
        type(mqc_nef_matrix),intent(in)::CoeffRow,CoeffCol
        type(mqc_nef_r4)::outR4
        integer(kind=int64)::axisO1,axisO2,nRowOut,nColOut,nO1,nO2
        integer(kind=int64)::pOut,qOut,idO1,idO2,iRow,iCol
        type(mqc_twoeris),pointer::blkSrc
        type(mqc_twoeris)::blkStage,blkStageEmpty
        type(mqc_scf_integral)::coeffBlockRow,coeffBlockCol
        type(mqc_r4tensor)::rawSrc,contribution,acc
        type(mqc_matrix)::denseRow,denseCol
        logical::first
!
        if(axisRow.eq.3.and.axisCol.eq.4) then
          axisO1=1; axisO2=2
        else if(axisRow.eq.1.and.axisCol.eq.2) then
          axisO1=3; axisO2=4
        else
          call mqc_error_a('mqc_nef_r4_full_transform: unsupported axis pair in buildStage',6,&
               'axisRow',' ')
        end if
!
        nRowOut = CoeffRow%NumSpace_col
        nColOut = CoeffCol%NumSpace_col
        nO1 = srcR4%blockSize(axisO1)
        nO2 = srcR4%blockSize(axisO2)
!
        call setAxisCount(outR4,axisRow,nRowOut)
        call setAxisCount(outR4,axisCol,nColOut)
        call setAxisCount(outR4,axisO1,nO1)
        call setAxisCount(outR4,axisO2,nO2)
        call mqc_nef_r4_allocate_slots(outR4)
        call setAxisMap(outR4,axisRow,CoeffRow%getMap(2),&
             CoeffRow%Size_vector_alpha_col,CoeffRow%Size_vector_beta_col)
        call setAxisMap(outR4,axisCol,CoeffCol%getMap(2),&
             CoeffCol%Size_vector_alpha_col,CoeffCol%Size_vector_beta_col)
        call setAxisMap(outR4,axisO1,srcR4%blockMap(axisO1),&
             srcR4%blockSizeVec(axisO1,'alpha'),srcR4%blockSizeVec(axisO1,'beta'))
        call setAxisMap(outR4,axisO2,srcR4%blockMap(axisO2),&
             srcR4%blockSizeVec(axisO2,'alpha'),srcR4%blockSizeVec(axisO2,'beta'))
!
        do idO1 = 1, nO1
          do idO2 = 1, nO2
            do pOut = 1, nRowOut
              do qOut = 1, nColOut
                first = .true.
                do iRow = 1, srcR4%blockSize(axisRow)
                  do iCol = 1, srcR4%blockSize(axisCol)
                    if(srcR4%isBlockZero([axisRow,iRow],[axisCol,iCol],&
                         [axisO1,idO1],[axisO2,idO2])) cycle
                    if(CoeffRow%isBlockZero(iRow,pOut)) cycle
                    if(CoeffCol%isBlockZero(iCol,qOut)) cycle
                    blkSrc => srcR4%iref([axisRow,iRow],[axisCol,iCol],&
                         [axisO1,idO1],[axisO2,idO2])
                    rawSrc = blkSrc%getBlock('full')
                    coeffBlockRow = CoeffRow%iat(iRow,pOut)
                    coeffBlockCol = CoeffCol%iat(iCol,qOut)
                    denseRow = coeffBlockRow%getBlock('full')
                    denseCol = coeffBlockCol%getBlock('full')
                    contribution = rawAxisTransform(rawSrc,axisRow,axisCol,denseRow,denseCol)
                    if(first) then
                      acc = contribution
                      first = .false.
                    else
                      acc = acc + contribution
                    end if
                  end do
                end do
                if(.not.first) then
                  blkStage = blkStageEmpty
                  call mqc_twoeris_allocate(blkStage,'full-coulomb','space',acc)
                  call outR4%iput(blkStage,[axisRow,pOut],[axisCol,qOut],&
                       [axisO1,idO1],[axisO2,idO2])
                end if
              end do
            end do
          end do
        end do
        end function buildStage
!
        subroutine setAxisCount(r4Obj,axisNum,nVal)
        implicit none
        type(mqc_nef_r4),intent(inout)::r4Obj
        integer(kind=int64),intent(in)::axisNum,nVal
        if(axisNum.eq.1) then
          r4Obj%NumSpace_1 = nVal
        else if(axisNum.eq.2) then
          r4Obj%NumSpace_2 = nVal
        else if(axisNum.eq.3) then
          r4Obj%NumSpace_3 = nVal
        else
          r4Obj%NumSpace_4 = nVal
        end if
        end subroutine setAxisCount
!
        subroutine setAxisMap(r4Obj,axisNum,mapIn,sizeAlphaIn,sizeBetaIn)
        implicit none
        type(mqc_nef_r4),intent(inout)::r4Obj
        integer(kind=int64),intent(in)::axisNum
        type(mqc_scf_eigenvalues),dimension(:),intent(in)::mapIn
        integer,dimension(:),intent(in)::sizeAlphaIn,sizeBetaIn
        if(axisNum.eq.1) then
          r4Obj%Map_Matrix_1 = mapIn
          r4Obj%Size_vector_alpha_1 = sizeAlphaIn
          r4Obj%Size_vector_beta_1 = sizeBetaIn
        else if(axisNum.eq.2) then
          r4Obj%Map_Matrix_2 = mapIn
          r4Obj%Size_vector_alpha_2 = sizeAlphaIn
          r4Obj%Size_vector_beta_2 = sizeBetaIn
        else if(axisNum.eq.3) then
          r4Obj%Map_Matrix_3 = mapIn
          r4Obj%Size_vector_alpha_3 = sizeAlphaIn
          r4Obj%Size_vector_beta_3 = sizeBetaIn
        else
          r4Obj%Map_Matrix_4 = mapIn
          r4Obj%Size_vector_alpha_4 = sizeAlphaIn
          r4Obj%Size_vector_beta_4 = sizeBetaIn
        end if
        end subroutine setAxisMap
!
        function rawAxisTransform(rawIn,axisRow,axisCol,coeffRowDense,coeffColDense) result(rawOut)
!       Transforms the axes (axisRow,axisCol) = (3,4), otherwise (1,2), of rawIn:
!       out(..p..q..) = sum conjg(coeffRow(x,p)) coeffCol(y,q) in(..x..y..), with
!       R4_Trans_Full (identity on the two other axes, which keeps their symmetry).
        implicit none
        type(mqc_r4tensor),intent(in)::rawIn
        integer(kind=int64),intent(in)::axisRow,axisCol
        type(mqc_matrix),intent(in)::coeffRowDense,coeffColDense
        type(mqc_r4tensor)::rawOut
        type(mqc_matrix)::Id1,Id2
!
        if(axisRow.eq.3.and.axisCol.eq.4) then
          Id1 = identity(size(rawIn,1),size(rawIn,1))
          Id2 = identity(size(rawIn,2),size(rawIn,2))
          rawOut = R4_Trans_Full(rawIn,coeffRowDense,coeffColDense,Id1,Id2,Map=[3,4,1,2], &
            Conj=[.false.,.false.,.true.,.false.])
        else
          Id1 = identity(size(rawIn,3),size(rawIn,3))
          Id2 = identity(size(rawIn,4),size(rawIn,4))
          rawOut = R4_Trans_Full(rawIn,coeffRowDense,coeffColDense,Id1,Id2,Map=[1,2,3,4], &
            Conj=[.true.,.false.,.false.,.false.])
        end if
        end function rawAxisTransform
!
        subroutine fetchFull(R4In,pIn,qIn,rIn,sIn,d1,d2,d3,d4,tOut)
        implicit none
        type(mqc_nef_r4),intent(in),target::R4In
        integer(kind=int64),intent(in)::pIn,qIn,rIn,sIn,d1,d2,d3,d4
        type(mqc_r4tensor),intent(out)::tOut
        type(mqc_twoeris),pointer::blkIn
        if(R4In%isBlockZero(pIn,qIn,rIn,sIn)) then
          call tOut%init(d1,d2,d3,d4)
        else
          blkIn => R4In%iref(pIn,qIn,rIn,sIn)
          tOut = blkIn%getBlock('full')
        end if
        end subroutine fetchFull
!
        function sliceCorner(fullR4,spin1,spin2,spin3,spin4,nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4) &
             result(blockOut)
        implicit none
        type(mqc_r4tensor),intent(in)::fullR4
        character(len=1),intent(in)::spin1,spin2,spin3,spin4
        integer(kind=int64),intent(in)::nA1,nB1,nA2,nB2,nA3,nB3,nA4,nB4
        type(mqc_r4tensor)::blockOut
        integer(kind=int64)::lo1,hi1,lo2,hi2,lo3,hi3,lo4,hi4
!
        if(spin1.eq.'a') then; lo1=1; hi1=nA1; else; lo1=nA1+1; hi1=nA1+nB1; end if
        if(spin2.eq.'a') then; lo2=1; hi2=nA2; else; lo2=nA2+1; hi2=nA2+nB2; end if
        if(spin3.eq.'a') then; lo3=1; hi3=nA3; else; lo3=nA3+1; hi3=nA3+nB3; end if
        if(spin4.eq.'a') then; lo4=1; hi4=nA4; else; lo4=nA4+1; hi4=nA4+nB4; end if
        if(hi1.lt.lo1.or.hi2.lt.lo2.or.hi3.lt.lo3.or.hi4.lt.lo4) then
          call blockOut%init(max(hi1-lo1+1,0),max(hi2-lo2+1,0),max(hi3-lo3+1,0),max(hi4-lo4+1,0))
        else
          blockOut = fullR4%tat([lo1,hi1],[lo2,hi2],[lo3,hi3],[lo4,hi4])
        end if
        end function sliceCorner
!
      end function mqc_nef_r4_full_transform
!
      function mqc_twoeris_nef_r4_transform(ERIs,NEFCoeff_1,NEFCoeff_2,&
           NEFCoeff_3,NEFCoeff_4,symmetry,sameBraKet) result(NEFR4Out)
!
! Rotates all 4 legs of a flat AO-basis mqc_twoERIs into up to 4 independent
! NEF_Matrix bases. Built directly on TwoERI_Trans_Full (mqc_est.F03), which
! already handles 'space'/'spin'/'general' coefficient types against a flat
! 'regular' ERI correctly, with no space->spin duplication needed --
! independent of mqc_nef_r4_full_transform.
!
      implicit none
      type(mqc_twoERIs),intent(in)::ERIs
      type(mqc_nef_matrix),intent(in)::NEFCoeff_1
      type(mqc_nef_matrix),optional,intent(in)::NEFCoeff_2,NEFCoeff_3,NEFCoeff_4
      character(len=*),intent(in),optional::symmetry
      logical,intent(in),optional::sameBraKet
      type(mqc_scf_eigenvalues),dimension(:),allocatable::outMap1,outMap2,outMap3,outMap4
      character(len=16)::outSymm
      logical::pairedCoeff,outSBK
      type(mqc_nef_r4)::NEFR4Out
      type(mqc_nef_matrix)::myCoeff_2,myCoeff_3,myCoeff_4
      type(mqc_scf_integral)::block_1,block_2,block_3,block_4
      type(mqc_twoERIs)::localERIs,MO_ERIs
      integer(kind=int64)::p,q,r,s
!
      if(present(NEFCoeff_2)) then
        myCoeff_2 = NEFCoeff_2
      else
        myCoeff_2 = NEFCoeff_1
      end if
      if(present(NEFCoeff_3)) then
        myCoeff_3 = NEFCoeff_3
      else
        myCoeff_3 = NEFCoeff_1
      end if
      if(present(NEFCoeff_4)) then
        myCoeff_4 = NEFCoeff_4
      else
        myCoeff_4 = NEFCoeff_1
      end if
!
!     The TwoERIs carry no operator label: the caller states it (symmetry, sameBraKet;
!     defaults 'asymmetric' / .false.). As in full_transform, the label is kept for
!     coefficients (real or complex) with C1 = C4 and C2 = C3 (T1, B.5: TwoERI_Trans_Full
!     applies the chemist conjugation conjg(C1), C2, conjg(C3), C4); with it the output is
!     mirrored and only one block per dagger pair is transformed.
      outSymm = mqc_nef_symmetry_label(symmetry,'mqc_twoeris_nef_r4_transform')
      outSBK = .false.
      if(present(sameBraKet)) outSBK = sameBraKet
      if(outSymm.ne.'asymmetric') then
        pairedCoeff = .true.
        if(present(NEFCoeff_2).or.present(NEFCoeff_3).or.present(NEFCoeff_4)) &
          pairedCoeff = mqc_nef_matrix_equal(NEFCoeff_1,myCoeff_4).and.&
                        mqc_nef_matrix_equal(myCoeff_2,myCoeff_3)
        if(.not.pairedCoeff) outSymm = 'asymmetric'
      end if
      outMap1 = NEFCoeff_1%getMap(2)
      outMap2 = myCoeff_2%getMap(2)
      outMap3 = myCoeff_3%getMap(2)
      outMap4 = myCoeff_4%getMap(2)
      call NEFR4Out%init(outMap1,outMap2,outMap3,outMap4,symmetry=outSymm,sameBraKet=outSBK)
!
      do p = 1, NEFR4Out%NumSpace_1
        if(NEFCoeff_1%isBlockZero(1,p)) cycle
        block_1 = NEFCoeff_1%iat(1,p)
        do q = 1, NEFR4Out%NumSpace_2
          if(myCoeff_2%isBlockZero(1,q)) cycle
          block_2 = myCoeff_2%iat(1,q)
          do r = 1, NEFR4Out%NumSpace_3
            if(myCoeff_3%isBlockZero(1,r)) cycle
            block_3 = myCoeff_3%iat(1,r)
            do s = 1, NEFR4Out%NumSpace_4
              if(NEFR4Out%isMirrorSlot(int(p),int(q),int(r),int(s))) cycle
              if(myCoeff_4%isBlockZero(1,s)) cycle
              block_4 = myCoeff_4%iat(1,s)
              localERIs = ERIs
              call TwoERI_Trans_Full(6,0,block_1,block_2,block_3,block_4,localERIs,MO_ERIs)
              call NEFR4Out%iput(MO_ERIs,p,q,r,s)
            end do
          end do
        end do
      end do
!
      end function mqc_twoeris_nef_r4_transform
!
!----------------------------------------------------------------
!                                                               |
!     MQC_NEF_SET PROCEDURES                                   |
!                                                               |
!----------------------------------------------------------------
!
!
!     PROCEDURE MQC_Print_NEF_Set
!
!
!>    \brief <b> MQC_Print_NEF_Set is a subroutine used to print a
!>    MQC_NEF_Set object. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_Print_NEF_Set prints every (I,J) cell of a MQC_NEF_Set,
!>    labeling each with its grid position and delegating to the
!>    cell's own kind-specific print (Vector/Matrix/R4). A zero-flagged
!>    cell is reported as such without calling its print (it was never
!>    put, so there is nothing meaningful to print), but its per-axis
!>    alpha/beta block dimensions (NEFSet%Size_vector_alpha/beta_1..4,
!>    uniform across every cell of the grid) are printed alongside it.
!>
!>    I/J select which (I,J) pairs are printed:
!>      - neither present: every (I,J) pair (the default).
!>      - only I present: every J for that one I.
!>      - only J present: every I for that one J.
!>      - both present: just that single (I,J) cell.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to print.
!>    \endverbatim
!>
!>    \param[in] iOut
!>    \verbatim
!>        iOut is Integer(kind=int64)
!>        The Fortran file number to print to.
!>    \endverbatim
!>
!>    \param[in] header
!>    \verbatim
!>        Header is Character(Len=*)
!>        The title to print along with NEFSet.
!>    \endverbatim
!>
!>    \param[in] I
!>    \verbatim
!>        I is Integer,Optional
!>        Restrict printing to this grid I value. See Purpose.
!>    \endverbatim
!>
!>    \param[in] J
!>    \verbatim
!>        J is Integer,Optional
!>        Restrict printing to this grid J value. See Purpose.
!>    \endverbatim
!>
!>    \param[in] blank_at_top
!>    \verbatim
!>        blank_at_top is Logical,Optional
!>        = .True.:  print blank line above output
!>    \endverbatim
!>
!>    \param[in] blank_at_bottom
!>    \verbatim
!>        blank_at_bottom is Logical,Optional
!>        Unused (kept for signature symmetry with the other NEF print
!>        routines).
!>    \endverbatim
!>
!>    \param[in] formatStr
!>    \verbatim
!>        formatStr is Character(Len=*),Optional
!>        Forwarded to the underlying Vector/Matrix/R4 print.
!>    \endverbatim
!>
!>    \param[in] blockPrint
!>    \verbatim
!>        blockPrint is Logical,Optional
!>        Forwarded to the underlying Matrix print only.
!>    \endverbatim
!>
!>    \param[in] short_print_threshold
!>    \verbatim
!>        short_print_threshold is Real(kind=real64),Optional
!>        Forwarded to the underlying R4 print only.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_print_nef_set(NEFSet,iOut,header,I,J,&
          blank_at_top,blank_at_bottom,formatStr,blockPrint,short_print_threshold)
!
      implicit none
      class(mqc_nef_set),target::NEFSet
      integer(kind=int64),intent(in)::iOut
      character(len=*),intent(in)::header
      integer,intent(in),optional::I,J
      logical,intent(in),optional::blank_at_top,blank_at_bottom
      character(len=*),intent(in),optional::formatStr
      logical,intent(in),optional::blockPrint
      real(kind=real64),intent(in),optional::short_print_threshold
      integer::myI,myJ,loI,hiI,loJ,hiJ
      character(len=64)::cellHeader
      character(len=256)::dimStr
!
      type(mqc_nef_set),target::full_NEFSet
      class(mqc_nef_set),pointer::ref_NEFSet
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSet => NEFSet
      if(NEFSet%isMirrorCell(1,2)) then
        full_NEFSet = mqc_nef_set_full_grid(NEFSet)
        ref_NEFSet => full_NEFSet
      end if
 1000 Format(1x,A)
 1020 Format( " " )
!
      if(present(blank_at_top)) then
        if(blank_at_top) write(iout,1020)
      endif
      write(iout,1000) trim(header)
!
      if(present(I)) then
        loI = I ; hiI = I
      else
        loI = 1 ; hiI = NEFSet%NumI
      end if
      if(present(J)) then
        loJ = J ; hiJ = J
      else
        loJ = 1 ; hiJ = NEFSet%NumJ
      end if
!
      select case(trim(NEFSet%Type()))
      case('vector')
        dimStr = 'dims: alpha='//trim(axisDimsStr(NEFSet%Size_vector_alpha_1))//&
                 ' beta='//trim(axisDimsStr(NEFSet%Size_vector_beta_1))
       do myI = loI, hiI
          do myJ = loJ, hiJ
            cellHeader = '(I,J) = ('//trim(num2char(myI))//','//trim(num2char(myJ))//')'
            if(NEFSet%isZero(myI,myJ)) then
              write(iout,1000) trim(cellHeader)//' (all zero) '//trim(dimStr)
            else
              call ref_NEFSet%Set_Vector(myI,myJ)%print(iOut,trim(cellHeader),formatStr=formatStr)
            end if
          end do
        end do
      case('matrix')
        dimStr = 'dims: row alpha='//trim(axisDimsStr(NEFSet%Size_vector_alpha_1))//&
                 ' beta='//trim(axisDimsStr(NEFSet%Size_vector_beta_1))//&
                 ' col alpha='//trim(axisDimsStr(NEFSet%Size_vector_alpha_2))//&
                 ' beta='//trim(axisDimsStr(NEFSet%Size_vector_beta_2))
        do myI = loI, hiI
          do myJ = loJ, hiJ
            cellHeader = '(I,J) = ('//trim(num2char(myI))//','//trim(num2char(myJ))//')'
            if(NEFSet%isZero(myI,myJ)) then
              write(iout,1000) trim(cellHeader)//' (all zero) '//trim(dimStr)
            else
              call ref_NEFSet%Set_Matrix(myI,myJ)%print(iOut,trim(cellHeader),&
                   formatStr=formatStr,blockPrint=blockPrint)
            end if
          end do
        end do
      case('r4')
        dimStr = 'dims: axis1 alpha='//trim(axisDimsStr(NEFSet%Size_vector_alpha_1))//&
                 ' beta='//trim(axisDimsStr(NEFSet%Size_vector_beta_1))//&
                 ' axis2 alpha='//trim(axisDimsStr(NEFSet%Size_vector_alpha_2))//&
                 ' beta='//trim(axisDimsStr(NEFSet%Size_vector_beta_2))//&
                 ' axis3 alpha='//trim(axisDimsStr(NEFSet%Size_vector_alpha_3))//&
                 ' beta='//trim(axisDimsStr(NEFSet%Size_vector_beta_3))//&
                 ' axis4 alpha='//trim(axisDimsStr(NEFSet%Size_vector_alpha_4))//&
                 ' beta='//trim(axisDimsStr(NEFSet%Size_vector_beta_4))
        do myI = loI, hiI
          do myJ = loJ, hiJ
            cellHeader = '(I,J) = ('//trim(num2char(myI))//','//trim(num2char(myJ))//')'
            if(NEFSet%isZero(myI,myJ)) then
              write(iout,1000) trim(cellHeader)//' (all zero) '//trim(dimStr)
            else
              call ref_NEFSet%Set_R4(myI,myJ)%print(iOut,trim(cellHeader),&
                   formatStr=formatStr,short_print_threshold=short_print_threshold)
            end if
          end do
        end do
      case default
        write(iout,1000) '(open -- holds no elements)'
      end select
!
      contains
!
        function axisDimsStr(sizes) result(str)
        implicit none
        integer,dimension(:),intent(in)::sizes
        character(len=256)::str
        integer::k
!
        str = '['
        do k = 1, size(sizes)
          if(k.gt.1) str = trim(str)//','
          str = trim(str)//trim(num2char(sizes(k)))
        end do
        str = trim(str)//']'
!
        end function axisDimsStr
!
      end subroutine mqc_print_nef_set
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
!>    The optional sameBraKetSpace marks the grid as bra x ket over the
!>    same determinant list (default .false.); it is allowed only when
!>    nI = nJ, otherwise an error is raised. It cannot be changed later.
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
!>
!>    \param[in] sameBraKetSpace
!>    \verbatim
!>        sameBraKetSpace is Logical, Optional
!>        .true. if bra and ket of the grid run over the same determinant
!>        list (requires nI = nJ). Default .false..
!>    \endverbatim
!>
!>    \param[in] symmetry
!>    \verbatim
!>        symmetry is Character(len=*), Optional
!>        Operator label of the set and of every cell: 'hermitian',
!>        'antihermitian' or 'asymmetric' (default). Never changed after
!>        init. With sameBraKetSpace it gives the grid symmetry
!>        Set(J,I) = +/- op(Set(I,J)) (op = conj for Vector cells, dagger
!>        for Matrix/R4 cells); only the I >= J cells are then stored.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_set_initialize(NEFSet,nI,nJ,sameBraKetSpace,symmetry)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      integer,intent(in)::nI,nJ
      logical,intent(in),optional::sameBraKetSpace
      character(len=*),intent(in),optional::symmetry
!
      NEFSet%NumI = nI
      NEFSet%NumJ = nJ
      NEFSet%Array_Type = 'open'
      NEFSet%sameBraKetSpace = .false.
      if(present(sameBraKetSpace)) then
        if(sameBraKetSpace.and.nI.ne.nJ) call mqc_error_i('mqc_nef_set_initialize: '//&
          'sameBraKetSpace requires NumI = NumJ',6,'nI',nI,'nJ',nJ)
        NEFSet%sameBraKetSpace = sameBraKetSpace
      end if
      NEFSet%Symmetry = mqc_nef_symmetry_label(symmetry,'mqc_nef_set_initialize')
!     the grid symmetry also needs the cell maps: decided when they are established
      NEFSet%GridMirrored = .false.
!
      end subroutine mqc_nef_set_initialize
!
!
!     PROCEDURE MQC_NEF_Set_Is_Same_BraKet_Space
!
!>    \brief <b> MQC_NEF_Set_Is_Same_BraKet_Space (isSameBraKetSpace) returns
!>    whether the (I,J) grid of an MQC_NEF_Set runs over the same determinant
!>    list for bra and ket. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Is_Same_BraKet_Space returns the set's sameBraKetSpace
!>    flag. The flag is set only by MQC_NEF_Set_Initialize, default .false.,
!>    and can be .true. only when NumI = NumJ. When .true., every (I,I) cell
!>    must carry sameBraKet = .true. (checked by put; zero cells returned by
!>    at are built with it); all other cells carry .false..
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
      function mqc_nef_set_is_same_braket_space(NEFSet) result(sameBraKetSpace)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      logical::sameBraKetSpace
!
      sameBraKetSpace = NEFSet%sameBraKetSpace
!
      end function mqc_nef_set_is_same_braket_space
!
!
!     PROCEDURE MQC_NEF_Set_Get_Symmetry
!
!>    \brief <b> MQC_NEF_Set_Get_Symmetry (getSymm) returns the operator label of
!>    an MQC_NEF_Set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Get_Symmetry returns the operator label given at init
!>    ('hermitian', 'antihermitian' or 'asymmetric'); every cell carries
!>    the same label (checked by put).
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
      function mqc_nef_set_get_symmetry(NEFSet) result(symm)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      character(len=16)::symm
!
      symm = NEFSet%Symmetry
!
      end function mqc_nef_set_get_symmetry
!
!
!     PROCEDURE MQC_NEF_Set_Is_Mirror_Cell
!
!>    \brief <b> MQC_NEF_Set_Is_Mirror_Cell (isMirrorCell) reports whether the
!>    (I,J) cell of an MQC_NEF_Set is a mirror cell. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Is_Mirror_Cell returns .true. iff the grid symmetry holds
!>    (sameBraKetSpace, label not 'asymmetric', paired cell maps) and I < J:
!>    the cell is not stored and at returns +/- op of the (J,I) cell.
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
!>    \param[in] I, J
!>    \verbatim
!>        I, J are Integer
!>        Grid indices.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_is_mirror_cell(NEFSet,I,J) result(isMirror)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      integer,intent(in)::I,J
      logical::isMirror
!
      isMirror = (NEFSet%GridMirrored.and.I.lt.J)
!
      end function mqc_nef_set_is_mirror_cell
!
!
!     PROCEDURE MQC_NEF_Set_Full_Grid
!
!>    \brief <b> MQC_NEF_Set_Full_Grid (private) returns a copy of an MQC_NEF_Set
!>    with every cell stored. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Full_Grid copies an MQC_NEF_Set into the full grid layout:
!>    same kind, maps, label and sameBraKetSpace, but the mirror cells of a
!>    grid-symmetric set are materialized (+/- op of the stored cell). It is
!>    an internal read copy for routines that read the cells directly
!>    (transitional rule, as MQC_NEF_R4_Full_Storage); only the storage
!>    layout differs from the input.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The MQC NEF Set to copy.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_full_grid(NEFSet) result(SetOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      type(mqc_nef_set)::SetOut
      integer::I,J
!
      SetOut = NEFSet
      if(NEFSet%GridMirrored) then
        SetOut%GridMirrored = .false.
        do I = 1, NEFSet%NumI
          do J = I+1, NEFSet%NumJ
            SetOut%Zero_Flag(I,J) = NEFSet%isZero(I,J)
            if(SetOut%Zero_Flag(I,J)) cycle
            select case(trim(NEFSet%Array_Type))
            case('vector')
              SetOut%Set_Vector(I,J) = NEFSet%atVector(I,J)
            case('matrix')
              SetOut%Set_Matrix(I,J) = NEFSet%atMatrix(I,J)
            case('r4')
              SetOut%Set_R4(I,J) = NEFSet%atR4(I,J)
            end select
          end do
        end do
      end if
!
      end function mqc_nef_set_full_grid
!
!
!     PROCEDURE MQC_NEF_Set_Rebuild_Cell_Matrix
!
!>    \brief <b> MQC_NEF_Set_Rebuild_Cell_Matrix (private) copies a computed
!>    NEF_Matrix cell into a new object built with the cell flags of a set. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Rebuild_Cell_Matrix is used by NEF_Set products, where the
!>    NEF_Set is the builder of its cells: the (I,J) cell of a product is
!>    <I|O1 O2|J>, so its sameBraKet is sameBraKetSpace .and. I = J, which
!>    the NEF-level sum of the contributing products cannot know. A new
!>    NEF_Matrix is initialized with the set's label and that flag and the
!>    blocks of the computed cell are copied into it.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] cell
!>    \verbatim
!>        cell is Type(MQC_NEF_Matrix)
!>        The computed cell.
!>    \endverbatim
!>
!>    \param[in] symmetry
!>    \verbatim
!>        symmetry is Character(len=*)
!>        Operator label of the new cell.
!>    \endverbatim
!>
!>    \param[in] sameBraKet
!>    \verbatim
!>        sameBraKet is Logical
!>        Bra-ket flag of the new cell.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_rebuild_cell_matrix(cell,symmetry,sameBraKet) result(cellOut)
!
      implicit none
      type(mqc_nef_matrix),intent(in)::cell
      character(len=*),intent(in)::symmetry
      logical,intent(in)::sameBraKet
      type(mqc_nef_matrix)::cellOut
      integer::row,col
!
      call cellOut%init(cell%Map_Matrix_Row,initial_map_matrix_col=cell%Map_Matrix_Col,&
           symmetry=symmetry,sameBraKet=sameBraKet)
      do row = 1, cell%NumSpace_row
        do col = 1, cell%NumSpace_col
          if(cellOut%Mirrored.and.col.gt.row) cycle
          if(cell%isBlockZero(row,col)) cycle
          call cellOut%iput(cell%iat(row,col),row,col)
        end do
      end do
!
      end function mqc_nef_set_rebuild_cell_matrix
!
!
!     PROCEDURE MQC_NEF_Set_Put_Built_Vector
!
!>    \brief <b> MQC_NEF_Set_Put_Built_Vector (private) stores a computed cell in
!>    an MQC_NEF_Set, the set acting as the builder of the cell. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Put_Built_Vector is used by NEF_Set operations whose cells are
!>    computed with NEF arithmetic (products, sums with zero cells, set +/-
!>    object). The (I,J) cell of the set must carry the set label and
!>    sameBraKet = sameBraKetSpace .and. I = J. If the computed cell already
!>    does, it is stored as is (no copy); otherwise a new cell is built with
!>    those flags and the computed blocks (the NEF_Set is the builder of its
!>    cells; flags of an existing object are never changed). Used only where
!>    this drops a claim (label to 'asymmetric', sameBraKet to .false.).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The set to write into.
!>    \endverbatim
!>
!>    \param[in] cell
!>    \verbatim
!>        cell is Type(mqc_nef_vector)
!>        The computed cell.
!>    \endverbatim
!>
!>    \param[in] I, J
!>    \verbatim
!>        I, J are Integer
!>        Grid indices.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_set_put_built_vector(NEFSet,cell,I,J)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      type(mqc_nef_vector),intent(in)::cell
      integer,intent(in)::I,J
      type(mqc_nef_vector)::cellOut
      logical::myFlag
      integer::id1,id2,id3,id4
!
      myFlag = (NEFSet%sameBraKetSpace.and.(I.eq.J))
      if(cell%getSymm().eq.NEFSet%Symmetry.and.(cell%isSameBraKet().eqv.myFlag)) then
        call NEFSet%put(cell,I,J)
      else
        call cellOut%init(cell%Map_Matrix,symmetry=NEFSet%Symmetry,sameBraKet=myFlag)
        do id1 = 1, cell%NumSpace
          if(.not.cell%isBlockZero(id1)) call cellOut%iput(cell%iat(id1),id1)
        end do
        call NEFSet%put(cellOut,I,J)
      end if
!
      end subroutine mqc_nef_set_put_built_vector
!
!
!     PROCEDURE MQC_NEF_Set_Put_Built_Matrix
!
!>    \brief <b> MQC_NEF_Set_Put_Built_Matrix (private) stores a computed cell in
!>    an MQC_NEF_Set, the set acting as the builder of the cell. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Put_Built_Matrix is used by NEF_Set operations whose cells are
!>    computed with NEF arithmetic (products, sums with zero cells, set +/-
!>    object). The (I,J) cell of the set must carry the set label and
!>    sameBraKet = sameBraKetSpace .and. I = J. If the computed cell already
!>    does, it is stored as is (no copy); otherwise a new cell is built with
!>    those flags and the computed blocks (the NEF_Set is the builder of its
!>    cells; flags of an existing object are never changed). Used only where
!>    this drops a claim (label to 'asymmetric', sameBraKet to .false.).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The set to write into.
!>    \endverbatim
!>
!>    \param[in] cell
!>    \verbatim
!>        cell is Type(mqc_nef_matrix)
!>        The computed cell.
!>    \endverbatim
!>
!>    \param[in] I, J
!>    \verbatim
!>        I, J are Integer
!>        Grid indices.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_set_put_built_matrix(NEFSet,cell,I,J)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      type(mqc_nef_matrix),intent(in)::cell
      integer,intent(in)::I,J
      type(mqc_nef_matrix)::cellOut
      logical::myFlag
      integer::id1,id2,id3,id4
!
      myFlag = (NEFSet%sameBraKetSpace.and.(I.eq.J))
      if(cell%getSymm().eq.NEFSet%Symmetry.and.(cell%isSameBraKet().eqv.myFlag)) then
        call NEFSet%put(cell,I,J)
      else
        cellOut = mqc_nef_set_rebuild_cell_matrix(cell,NEFSet%Symmetry,myFlag)

        call NEFSet%put(cellOut,I,J)
      end if
!
      end subroutine mqc_nef_set_put_built_matrix
!
!
!     PROCEDURE MQC_NEF_Set_Put_Built_R4
!
!>    \brief <b> MQC_NEF_Set_Put_Built_R4 (private) stores a computed cell in
!>    an MQC_NEF_Set, the set acting as the builder of the cell. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Put_Built_R4 is used by NEF_Set operations whose cells are
!>    computed with NEF arithmetic (products, sums with zero cells, set +/-
!>    object). The (I,J) cell of the set must carry the set label and
!>    sameBraKet = sameBraKetSpace .and. I = J. If the computed cell already
!>    does, it is stored as is (no copy); otherwise a new cell is built with
!>    those flags and the computed blocks (the NEF_Set is the builder of its
!>    cells; flags of an existing object are never changed). Used only where
!>    this drops a claim (label to 'asymmetric', sameBraKet to .false.).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in,out] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The set to write into.
!>    \endverbatim
!>
!>    \param[in] cell
!>    \verbatim
!>        cell is Type(mqc_nef_r4)
!>        The computed cell.
!>    \endverbatim
!>
!>    \param[in] I, J
!>    \verbatim
!>        I, J are Integer
!>        Grid indices.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      subroutine mqc_nef_set_put_built_r4(NEFSet,cell,I,J)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      type(mqc_nef_r4),intent(in)::cell
      integer,intent(in)::I,J
      type(mqc_nef_r4)::cellOut
      logical::myFlag
      integer::id1,id2,id3,id4
!
      myFlag = (NEFSet%sameBraKetSpace.and.(I.eq.J))
      if(cell%getSymm().eq.NEFSet%Symmetry.and.(cell%isSameBraKet().eqv.myFlag)) then
        call NEFSet%put(cell,I,J)
      else
        call cellOut%init(cell%Map_Matrix_1,cell%Map_Matrix_2,cell%Map_Matrix_3,&
             cell%Map_Matrix_4,symmetry=NEFSet%Symmetry,sameBraKet=myFlag)
        do id1 = 1, cell%NumSpace_1
          do id2 = 1, cell%NumSpace_2
            do id3 = 1, cell%NumSpace_3
              do id4 = 1, cell%NumSpace_4
                if(cell%isBlockZero(id1,id2,id3,id4)) cycle
                if(cellOut%isMirrorSlot(id1,id2,id3,id4)) cycle
                call cellOut%iput(cell%iat(id1,id2,id3,id4),id1,id2,id3,id4)
              end do
            end do
          end do
        end do
        call NEFSet%put(cellOut,I,J)
      end if
!
      end subroutine mqc_nef_set_put_built_r4
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
      if(allocated(NEFSet%Zero_Flag)) then
        if(NEFSet%GridMirrored.and.I.lt.J) then
          isZero = NEFSet%Zero_Flag(J,I)
        else
          isZero = NEFSet%Zero_Flag(I,J)
        end if
      end if
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
      subroutine mqc_nef_set_allocate(NEFSet,NEFType,&
           Map_Matrix_1,Map_Matrix_2,Map_Matrix_3,Map_Matrix_4)
!
      implicit none
      class(mqc_nef_set),intent(inOut)::NEFSet
      character(len=*),intent(in)::NEFType
!
      Type(MQC_SCF_Eigenvalues),Dimension(:),&
        Allocatable,intent(in),optional::Map_Matrix_1,Map_Matrix_2,&
                                         Map_Matrix_3,Map_Matrix_4
      
!
      character(len=64)::myType
      integer::i
!
      if(allocated(NEFSet%Set_Vector)) deallocate(NEFSet%Set_Vector)
      if(allocated(NEFSet%Set_Matrix)) deallocate(NEFSet%Set_Matrix)
      if(allocated(NEFSet%Set_R4)) deallocate(NEFSet%Set_R4)
      if(allocated(NEFSet%Size_vector_alpha_1)) deallocate(NEFSet%Size_vector_alpha_1)
      if(allocated(NEFSet%Size_vector_alpha_2)) deallocate(NEFSet%Size_vector_alpha_2)
      if(allocated(NEFSet%Size_vector_alpha_3)) deallocate(NEFSet%Size_vector_alpha_3)
      if(allocated(NEFSet%Size_vector_alpha_4)) deallocate(NEFSet%Size_vector_alpha_4)
      if(allocated(NEFSet%Size_vector_beta_1)) deallocate(NEFSet%Size_vector_beta_1)
      if(allocated(NEFSet%Size_vector_beta_2)) deallocate(NEFSet%Size_vector_beta_2)
      if(allocated(NEFSet%Size_vector_beta_3)) deallocate(NEFSet%Size_vector_beta_3)
      if(allocated(NEFSet%Size_vector_beta_4)) deallocate(NEFSet%Size_vector_beta_4)
!
      call string_change_case(NEFType,'l',myType)
!
      select case(trim(myType))
      case('vector')
        if(.not.(present(Map_Matrix_1))) &
          call mqc_error('mqc_nef_set_allocate: vector requires '//&
               'Map_Matrix_1',6)
        if(present(Map_Matrix_2).or.present(Map_Matrix_3).or.present(Map_Matrix_4)) &
          call mqc_error('mqc_nef_set_allocate: vector only accepts axis-1 Map',6)
        allocate(NEFSet%Set_Vector(NEFSet%NumI,NEFSet%NumJ))
        NEFSet%NumSpace(1) = size(Map_Matrix_1)
        NEFSet%Map_Matrix_1 = Map_Matrix_1
        allocate(NEFSet%Size_vector_alpha_1(NEFSet%NumSpace(1)),&
                 NEFSet%Size_vector_beta_1(NEFSet%NumSpace(1)))
        do i = 1, NEFSet%NumSpace(1)
          NEFSet%Size_vector_alpha_1(i) = Size(Map_Matrix_1(i),'alpha')
          NEFSet%Size_vector_beta_1(i)  = Size(Map_Matrix_1(i),'beta')
        end do
      case('matrix')
        if(.not.(present(Map_Matrix_1).and.present(Map_Matrix_2))) &
          call mqc_error('mqc_nef_set_allocate: matrix requires '//&
               'axis-1 and axis-2 Map',6)
        if(present(Map_Matrix_3).or.present(Map_Matrix_4)) &
          call mqc_error('mqc_nef_set_allocate: matrix only accepts '//&
               'axis-1 and axis-2 Map',6)
        allocate(NEFSet%Set_Matrix(NEFSet%NumI,NEFSet%NumJ))
        NEFSet%NumSpace(1) = size(Map_Matrix_1)
        NEFSet%NumSpace(2) = size(Map_Matrix_2)
        NEFSet%Map_Matrix_1 = Map_Matrix_1
        NEFSet%Map_Matrix_2 = Map_Matrix_2
        allocate(NEFSet%Size_vector_alpha_1(NEFSet%NumSpace(1)),&
                 NEFSet%Size_vector_beta_1(NEFSet%NumSpace(1)))
        allocate(NEFSet%Size_vector_alpha_2(NEFSet%NumSpace(2)),&
                 NEFSet%Size_vector_beta_2(NEFSet%NumSpace(2)))
        do i = 1, NEFSet%NumSpace(1)
          NEFSet%Size_vector_alpha_1(i) = Size(Map_Matrix_1(i),'alpha')
          NEFSet%Size_vector_beta_1(i)  = Size(Map_Matrix_1(i),'beta')
        end do
        do i = 1, NEFSet%NumSpace(2)
          NEFSet%Size_vector_alpha_2(i) = Size(Map_Matrix_2(i),'alpha')
          NEFSet%Size_vector_beta_2(i)  = Size(Map_Matrix_2(i),'beta')
        end do
      case('r4')
        if(.not.(present(Map_Matrix_1).and.present(Map_Matrix_2).and.&
                 present(Map_Matrix_3).and.present(Map_Matrix_4))) &
          call mqc_error('mqc_nef_set_allocate: r4 requires all four Map arguments',6)
        allocate(NEFSet%Set_R4(NEFSet%NumI,NEFSet%NumJ))
        NEFSet%NumSpace(1) = size(Map_Matrix_1)
        NEFSet%NumSpace(2) = size(Map_Matrix_2)
        NEFSet%NumSpace(3) = size(Map_Matrix_3)
        NEFSet%NumSpace(4) = size(Map_Matrix_4)
        NEFSet%Map_Matrix_1 = Map_Matrix_1
        NEFSet%Map_Matrix_2 = Map_Matrix_2
        NEFSet%Map_Matrix_3 = Map_Matrix_3
        NEFSet%Map_Matrix_4 = Map_Matrix_4
        allocate(NEFSet%Size_vector_alpha_1(NEFSet%NumSpace(1)),&
                 NEFSet%Size_vector_beta_1(NEFSet%NumSpace(1)))
        allocate(NEFSet%Size_vector_alpha_2(NEFSet%NumSpace(2)),&
                 NEFSet%Size_vector_beta_2(NEFSet%NumSpace(2)))
        allocate(NEFSet%Size_vector_alpha_3(NEFSet%NumSpace(3)),&
                 NEFSet%Size_vector_beta_3(NEFSet%NumSpace(3)))
        allocate(NEFSet%Size_vector_alpha_4(NEFSet%NumSpace(4)),&
                 NEFSet%Size_vector_beta_4(NEFSet%NumSpace(4)))
        do i = 1, NEFSet%NumSpace(1)
          NEFSet%Size_vector_alpha_1(i) = Size(Map_Matrix_1(i),'alpha')
          NEFSet%Size_vector_beta_1(i)  = Size(Map_Matrix_1(i),'beta')
        end do
        do i = 1, NEFSet%NumSpace(2)
          NEFSet%Size_vector_alpha_2(i) = Size(Map_Matrix_2(i),'alpha')
          NEFSet%Size_vector_beta_2(i)  = Size(Map_Matrix_2(i),'beta')
        end do
        do i = 1, NEFSet%NumSpace(3)
          NEFSet%Size_vector_alpha_3(i) = Size(Map_Matrix_3(i),'alpha')
          NEFSet%Size_vector_beta_3(i)  = Size(Map_Matrix_3(i),'beta')
        end do
        do i = 1, NEFSet%NumSpace(4)
          NEFSet%Size_vector_alpha_4(i) = Size(Map_Matrix_4(i),'alpha')
          NEFSet%Size_vector_beta_4(i)  = Size(Map_Matrix_4(i),'beta')
        end do
      case default
        call mqc_error_A('mqc_nef_set_allocate: unknown kind',6,'NEFType',trim(NEFType))
      end select
      if(allocated(NEFSet%Zero_Flag)) deallocate(NEFSet%Zero_Flag)
      allocate(NEFSet%Zero_Flag(NEFSet%NumI,NEFSet%NumJ))
      NEFSet%Zero_Flag = .true.
      NEFSet%Array_Type = trim(NEFType)
!     Grid symmetry: Set(J,I) = +/- op(Set(I,J)) needs the op to keep the cell maps
!     (Matrix: rows = cols; R4: axes 1/4 and 2/3).
      NEFSet%GridMirrored = (NEFSet%sameBraKetSpace.and.NEFSet%Symmetry.ne.'asymmetric')
      if(NEFSet%GridMirrored) then
        select case(trim(myType))
        case('matrix')
          NEFSet%GridMirrored = mqc_nef_maps_equal(NEFSet%Map_Matrix_1,NEFSet%Map_Matrix_2)
        case('r4')
          NEFSet%GridMirrored = mqc_nef_maps_equal(NEFSet%Map_Matrix_1,NEFSet%Map_Matrix_4).and.&
                                mqc_nef_maps_equal(NEFSet%Map_Matrix_2,NEFSet%Map_Matrix_3)
        end select
      end if
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
!>    The element's sameBraKet flag must be .true. iff the set's
!>    sameBraKetSpace is .true. and I = J; otherwise an error is raised
!>    (the flag is given by the builder at the element's init).
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
        call mqc_nef_set_allocate(NEFSet,'vector',Map_Matrix_1=element%Map_Matrix)
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
      if(.not.(mqc_nef_maps_equal(element%Map_Matrix,NEFSet%Map_Matrix_1))) &
        call mqc_error_i('mqc_nef_set_put_vector: element maps do not match the set''s '//&
          'established maps',6,'I',I,'J',J)
      if(element%isSameBraKet().neqv.(NEFSet%sameBraKetSpace.and.(I.eq.J))) &
        call mqc_error_i('mqc_nef_set_put_vector: element sameBraKet does not match the '//&
          'set (sameBraKetSpace .and. I = J)',6,'I',I,'J',J)
      if(element%getSymm().ne.NEFSet%Symmetry) &
        call mqc_error_a('mqc_nef_set_put_vector: element symmetry does not match the set',&
          6,'element%getSymm()',element%getSymm(),'NEFSet%Symmetry',NEFSet%Symmetry)
      if(NEFSet%GridMirrored.and.I.lt.J) then
!       Grid symmetry: the stored (J,I) cell holds the +/- op of the (I,J) cell.
        if(NEFSet%Symmetry.eq.'hermitian') then
          NEFSet%Set_Vector(J,I) = conjg(element)
        else
          NEFSet%Set_Vector(J,I) = (-1.0)*conjg(element)
        end if
        NEFSet%Zero_Flag(J,I) = .false.
      else
        NEFSet%Set_Vector(I,J) = element
        NEFSet%Zero_Flag(I,J) = .false.
      end if
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
!>    The element's sameBraKet flag must be .true. iff the set's
!>    sameBraKetSpace is .true. and I = J; otherwise an error is raised
!>    (the flag is given by the builder at the element's init).
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
             Map_Matrix_1=element%Map_Matrix_Row,&
             Map_Matrix_2=element%Map_Matrix_Col)
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
      if(.not.(mqc_nef_maps_equal(element%Map_Matrix_Row,NEFSet%Map_Matrix_1).and.&
          mqc_nef_maps_equal(element%Map_Matrix_Col,NEFSet%Map_Matrix_2))) &
        call mqc_error_i('mqc_nef_set_put_matrix: element maps do not match the set''s '//&
          'established maps',6,'I',I,'J',J)
      if(element%isSameBraKet().neqv.(NEFSet%sameBraKetSpace.and.(I.eq.J))) &
        call mqc_error_i('mqc_nef_set_put_matrix: element sameBraKet does not match the '//&
          'set (sameBraKetSpace .and. I = J)',6,'I',I,'J',J)
      if(element%getSymm().ne.NEFSet%Symmetry) &
        call mqc_error_a('mqc_nef_set_put_matrix: element symmetry does not match the set',&
          6,'element%getSymm()',element%getSymm(),'NEFSet%Symmetry',NEFSet%Symmetry)
      if(NEFSet%GridMirrored.and.I.lt.J) then
!       Grid symmetry: the stored (J,I) cell holds the +/- op of the (I,J) cell.
        if(NEFSet%Symmetry.eq.'hermitian') then
          NEFSet%Set_Matrix(J,I) = dagger(element)
        else
          NEFSet%Set_Matrix(J,I) = (-1.0)*dagger(element)
        end if
        NEFSet%Zero_Flag(J,I) = .false.
      else
        NEFSet%Set_Matrix(I,J) = element
        NEFSet%Zero_Flag(I,J) = .false.
      end if
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
!>    The element's sameBraKet flag must be .true. iff the set's
!>    sameBraKetSpace is .true. and I = J; otherwise an error is raised
!>    (the flag is given by the builder at the element's init).
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
             Map_Matrix_1=element%Map_Matrix_1,&
             Map_Matrix_2=element%Map_Matrix_2,&
             Map_Matrix_3=element%Map_Matrix_3,&
             Map_Matrix_4=element%Map_Matrix_4)
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
      if(.not.(mqc_nef_maps_equal(element%Map_Matrix_1,NEFSet%Map_Matrix_1).and.&
          mqc_nef_maps_equal(element%Map_Matrix_2,NEFSet%Map_Matrix_2).and.&
          mqc_nef_maps_equal(element%Map_Matrix_3,NEFSet%Map_Matrix_3).and.&
          mqc_nef_maps_equal(element%Map_Matrix_4,NEFSet%Map_Matrix_4))) &
        call mqc_error_i('mqc_nef_set_put_r4: element maps do not match the set''s '//&
          'established maps',6,'I',I,'J',J)
      if(element%isSameBraKet().neqv.(NEFSet%sameBraKetSpace.and.(I.eq.J))) &
        call mqc_error_i('mqc_nef_set_put_r4: element sameBraKet does not match the '//&
          'set (sameBraKetSpace .and. I = J)',6,'I',I,'J',J)
      if(element%getSymm().ne.NEFSet%Symmetry) &
        call mqc_error_a('mqc_nef_set_put_r4: element symmetry does not match the set',&
          6,'element%getSymm()',element%getSymm(),'NEFSet%Symmetry',NEFSet%Symmetry)
      if(NEFSet%GridMirrored.and.I.lt.J) then
!       Grid symmetry: the stored (J,I) cell holds the +/- op of the (I,J) cell.
        if(NEFSet%Symmetry.eq.'hermitian') then
          NEFSet%Set_R4(J,I) = dagger(element)
        else
          NEFSet%Set_R4(J,I) = (-1.0)*dagger(element)
        end if
        NEFSet%Zero_Flag(J,I) = .false.
      else
        NEFSet%Set_R4(I,J) = element
        NEFSet%Zero_Flag(I,J) = .false.
      end if
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
!>    A zero cell is returned with sameBraKet = .true. iff the set's
!>    sameBraKetSpace is .true. and I = J.
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
      if(NEFSet%isZero(I,J)) then
        call element%init(NEFSet%Map_Matrix_1,symmetry=NEFSet%Symmetry,&
             sameBraKet=(NEFSet%sameBraKetSpace.and.(I.eq.J)))
      else if(NEFSet%GridMirrored.and.I.lt.J) then
!       Grid symmetry: +/- conj of the stored (J,I) cell.
        element = conjg(NEFSet%Set_Vector(J,I))
        if(NEFSet%Symmetry.eq.'antihermitian') element = (-1.0)*element
      else
        element = NEFSet%Set_Vector(I,J)
      end if
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
!>    A zero cell is returned with sameBraKet = .true. iff the set's
!>    sameBraKetSpace is .true. and I = J.
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
      if(NEFSet%isZero(I,J)) then
        call element%init(NEFSet%Map_Matrix_1,initial_map_matrix_col=NEFSet%Map_Matrix_2,&
             symmetry=NEFSet%Symmetry,sameBraKet=(NEFSet%sameBraKetSpace.and.(I.eq.J)))
      else if(NEFSet%GridMirrored.and.I.lt.J) then
!       Grid symmetry: +/- dagger of the stored (J,I) cell.
        element = dagger(NEFSet%Set_Matrix(J,I))
        if(NEFSet%Symmetry.eq.'antihermitian') element = (-1.0)*element
      else
        element = NEFSet%Set_Matrix(I,J)
      end if
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
!>    A zero cell is returned with sameBraKet = .true. iff the set's
!>    sameBraKetSpace is .true. and I = J.
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
      if(NEFSet%isZero(I,J)) then
        call element%init(NEFSet%Map_Matrix_1,NEFSet%Map_Matrix_2,&
             NEFSet%Map_Matrix_3,NEFSet%Map_Matrix_4,symmetry=NEFSet%Symmetry,&
             sameBraKet=(NEFSet%sameBraKetSpace.and.(I.eq.J)))
      else if(NEFSet%GridMirrored.and.I.lt.J) then
!       Grid symmetry: +/- dagger of the stored (J,I) cell.
        element = dagger(NEFSet%Set_R4(J,I))
        if(NEFSet%Symmetry.eq.'antihermitian') element = (-1.0)*element
      else
        element = NEFSet%Set_R4(I,J)
      end if
!
      end function mqc_nef_set_at_r4
!
!
!     PROCEDURE MQC_NEF_Set_Get_Map
!
!
!>    \brief <b> MQC_NEF_Set_Get_Map returns one of a 
!>    MQC_NEF_Set's cells' internal block maps. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Class(MQC_NEF_Set)
!>        The NEF_Set to  
!>    \endverbatim
!>
!>    \param[in] acis
!>    \verbatim
!>        acis is Integer,Optional
!>        1 (default), 2, 3 or 4
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_get_map(NEFSet,axis) result(mapOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      integer,intent(in),optional::axis
      integer::myaxis
      type(mqc_scf_eigenvalues),dimension(:),allocatable::mapOut
!
      if(present(axis)) then
        myaxis = axis
      else
        myaxis = 1
      end if
!
      if(myaxis.le.0.or.myaxis.gt.4) call mqc_error_I(&
           'mqc_nef_set_get_map: Axis out-of range',&
           6,'myaxis',myaxis)
      if(NEFSet%isOpen()) call mqc_error(& 
           'mqc_nef_set_get_map: NEFSet type is open, no axis defined',&
           6)
      if(NEFSet%isVector().and.myaxis.ne.1) call mqc_error_I(& 
           'mqc_nef_set_get_map: NEFSet type is vector, axis out-of-bound',&
           6,'myaxis',myaxis)
      if(NEFSet%isMatrix().and.myaxis.gt.2) call mqc_error_I(& 
           'mqc_nef_set_get_map: NEFSet type is matrix, axis out-of-bound',&
           6,'myaxis',myaxis)
!
      if(myaxis.eq.1) mapOut = NEFSet%Map_Matrix_1
      if(myaxis.eq.2) mapOut = NEFSet%Map_Matrix_2
      if(myaxis.eq.3) mapOut = NEFSet%Map_Matrix_3
      if(myaxis.eq.4) mapOut = NEFSet%Map_Matrix_4
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
          if(.not.first) call mqc_nef_set_put_built_vector(NEFSetOut,tmpVector,I,1)
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
          if(.not.first) call mqc_nef_set_put_built_matrix(NEFSetOut,tmpMatrix,I,1)
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
          if(.not.first) call mqc_nef_set_put_built_r4(NEFSetOut,tmpR4,I,1)
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
          if(.not.first) call mqc_nef_set_put_built_vector(NEFSetOut,tmpVector,1,J)
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
          if(.not.first) call mqc_nef_set_put_built_matrix(NEFSetOut,tmpMatrix,1,J)
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
          if(.not.first) call mqc_nef_set_put_built_r4(NEFSetOut,tmpR4,1,J)
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
            if(.not.first) call mqc_nef_set_put_built_vector(NEFSetOut,tmpVector,I,K)
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
            if(.not.first) call mqc_nef_set_put_built_matrix(NEFSetOut,tmpMatrix,I,K)
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
            if(.not.first) call mqc_nef_set_put_built_r4(NEFSetOut,tmpR4,I,K)
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
            if(.not.first) call mqc_nef_set_put_built_vector(NEFSetOut,tmpVector,K,J)
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
            if(.not.first) call mqc_nef_set_put_built_matrix(NEFSetOut,tmpMatrix,K,J)
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
            if(.not.first) call mqc_nef_set_put_built_r4(NEFSetOut,tmpR4,K,J)
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
          call mqc_nef_set_put_built_matrix(NEFSetOut,MatMul(NEFSetA%atMatrix(I,J),NEFMatB),I,J)
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
!>    MQC_NEF_Matrix into every (I,J) cell of a matrix- or vector-kind
!>    MQC_NEF_Set, broadcasting the same NEFMatA across the whole
!>    grid. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEFMatrix_NEFSet_Multiply computes, for every (I,J),
!>    NEFSetOut(I,J) = MatMul(NEFMatA,NEFSetB(I,J)), where NEFMatA is a
!>    single MQC_NEF_Matrix shared across every grid cell. The result is
!>    a MQC_NEF_Set of the same (I,J) shape as NEFSetB, holding Matrix
!>    elements if NEFSetB holds Matrix elements, or Vector elements if
!>    NEFSetB holds Vector elements (via MQC_NEFMatrix_NEFVector_
!>    Multiply per cell). Zero-flagged cells of NEFSetB are left unset
!>    in the output. Every cell's row/column space must be compatible
!>    with NEFMatA's column space (guaranteed uniform across cells by
!>    MQC_NEF_Set's own put-time dimension check).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFMatA
!>    \verbatim
!>        NEFMatA is Type(MQC_NEF_Matrix)
!>        The matrix multiplied into every (I,J) cell.
!>    \endverbatim
!>
!>    \param[in] NEFSetB
!>    \verbatim
!>        NEFSetB is Type(MQC_NEF_Set)
!>        The (I,J) grid of Matrix or Vector elements to transform.
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
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ)
!
      select case(trim(NEFSetB%Type()))
      case('matrix')
        do I = 1, NEFSetB%NumI
          do J = 1, NEFSetB%NumJ
            if(NEFSetB%isZero(I,J)) cycle
            call mqc_nef_set_put_built_matrix(NEFSetOut,MatMul(NEFMatA,NEFSetB%atMatrix(I,J)),I,J)
          end do
        end do
      case('vector')
        do I = 1, NEFSetB%NumI
          do J = 1, NEFSetB%NumJ
            if(NEFSetB%isZero(I,J)) cycle
            call mqc_nef_set_put_built_vector(NEFSetOut,MatMul(NEFMatA,NEFSetB%atVector(I,J)),I,J)
          end do
        end do
      case default
        call mqc_error_A('mqc_nefmatrix_nefset_multiply: NEF_Set must hold '//&
             'matrix or vector elements',6,'Array_Type',NEFSetB%Type())
      end select
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
           6,'NEFSetA%Type',NEFSetA%Type(),'NEFSetB%Type',NEFSetB%Type())
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
            if(.not.first) call mqc_nef_set_put_built_matrix(NEFOut,tmpMatrix,I,K)
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
            if(.not.first) call mqc_nef_set_put_built_r4(NEFOut,tmpR4,I,K)
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
           6,'NEFSetA%Type',NEFSetA%Type(),'NEFSetB%Type',NEFSetB%Type())
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
      type(mqc_scf_eigenvalues),dimension(:),allocatable::outMap1,outMap2
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
            if(.not.first) call mqc_nef_set_put_built_matrix(NEFSetOut,tmpMatrix,I,K)
          end do
        end do
        if(NEFSetOut%isOpen()) then
          if(present(index1)) then
            c=index1(3) ; d=index1(4)
          else
            c=3 ; d=4
          end if
          select case(c)
          case(1); outMap1 = NEFSetA%Map_Matrix_1
          case(2); outMap1 = NEFSetA%Map_Matrix_2
          case(3); outMap1 = NEFSetA%Map_Matrix_3
          case(4); outMap1 = NEFSetA%Map_Matrix_4
          end select
          select case(d)
          case(1); outMap2 = NEFSetA%Map_Matrix_1
          case(2); outMap2 = NEFSetA%Map_Matrix_2
          case(3); outMap2 = NEFSetA%Map_Matrix_3
          case(4); outMap2 = NEFSetA%Map_Matrix_4
          end select
          call mqc_nef_set_allocate(NEFSetOut,'matrix',Map_Matrix_1=outMap1,Map_Matrix_2=outMap2)
        end if
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
              call mqc_nef_set_put_built_matrix(NEFSetOut,tmpMatrix,I,K)
            end if
          end do
        end do
        if(NEFSetOut%isOpen()) then
          if(myTransposeA) then
            outMap1 = NEFSetA%Map_Matrix_2
          else
            outMap1 = NEFSetA%Map_Matrix_1
          end if
          if(myTransposeB) then
            outMap2 = NEFSetB%Map_Matrix_1
          else
            outMap2 = NEFSetB%Map_Matrix_2
          end if
          if(swapOut) then
            call mqc_nef_set_allocate(NEFSetOut,'matrix',Map_Matrix_1=outMap2,Map_Matrix_2=outMap1)
          else
            call mqc_nef_set_allocate(NEFSetOut,'matrix',Map_Matrix_1=outMap1,Map_Matrix_2=outMap2)
          end if
        end if
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
            if(.not.first) call mqc_nef_set_put_built_vector(NEFSetOut,tmpVector,I,K)
          end do
        end do
        if(NEFSetOut%isOpen()) then
          if(myTransposeA) then
            call mqc_nef_set_allocate(NEFSetOut,'vector',Map_Matrix_1=NEFSetA%Map_Matrix_2)
          else
            call mqc_nef_set_allocate(NEFSetOut,'vector',Map_Matrix_1=NEFSetA%Map_Matrix_1)
          end if
        end if
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
            if(.not.first) call mqc_nef_set_put_built_vector(NEFSetOut,tmpVector,I,K)
          end do
        end do
        if(NEFSetOut%isOpen()) then
          if(myTransposeB) then
            call mqc_nef_set_allocate(NEFSetOut,'vector',Map_Matrix_1=NEFSetB%Map_Matrix_2)
          else
            call mqc_nef_set_allocate(NEFSetOut,'vector',Map_Matrix_1=NEFSetB%Map_Matrix_1)
          end if
        end if
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
            if(.not.first) call mqc_nef_set_put_built_matrix(NEFSetOut,tmpMatrix,I,K)
          end do
        end do
        if(NEFSetOut%isOpen()) then
          if(present(index1)) then
            c=index1(3) ; d=index1(4)
          else
            c=3 ; d=4
          end if
          select case(c)
          case(1); outMap1 = NEFSetA%Map_Matrix_1
          case(2); outMap1 = NEFSetA%Map_Matrix_2
          case(3); outMap1 = NEFSetA%Map_Matrix_3
          case(4); outMap1 = NEFSetA%Map_Matrix_4
          end select
          select case(d)
          case(1); outMap2 = NEFSetA%Map_Matrix_1
          case(2); outMap2 = NEFSetA%Map_Matrix_2
          case(3); outMap2 = NEFSetA%Map_Matrix_3
          case(4); outMap2 = NEFSetA%Map_Matrix_4
          end select
          call mqc_nef_set_allocate(NEFSetOut,'matrix',Map_Matrix_1=outMap1,Map_Matrix_2=outMap2)
        end if
      case default
        call mqc_error_A('mqc_nefset_nefset_partial_contraction: unsupported NEF_Set '//&
             'kind combination',6,'NEFSetA%Type',NEFSetA%Type(),'NEFSetB%Type',NEFSetB%Type())
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
          if(.not.first) call mqc_nef_set_put_built_r4(NEFSetOut,tmpR4,I,K)
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
!>
!>    \param[in] swapOrder
!>    \verbatim
!>        Logical,Optional
!>        Forwarded to the underlying Vector expand only. Must match the
!>        swapOrder used when the Vector-kind cell was originally produced
!>        by a Matrix collapse, otherwise expand reads back the wrong flat
!>        positions (landing on zero-flagged entries instead of the real
!>        data) -- collapse's swapOrder is not auto-detected or stored.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_expand(NEFSet,map_matrix_1,map_matrix_2,&
           map_matrix_3,map_matrix_4,spinType,spinTypeKet,Order,swapOrder) result(NEFSetOut)
!
      implicit none
      class(mqc_nef_set),intent(in),target::NEFSet
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,intent(in),optional::&
           map_matrix_1,map_matrix_2,map_matrix_3,map_matrix_4
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      logical,intent(in),optional::swapOrder
      type(mqc_nef_set)::NEFSetOut
      logical::mySwap
      integer(kind=int64),dimension(4)::myOrder
      integer::I,J
      integer(kind=int64)::ia,ja,ka,la,s1,s2,s3,s4,expectedRow,expectedCol
      character(len=64)::mySpinType,mySpinTypeKet
      type(mqc_scf_eigenvalues),dimension(:),allocatable::outMap1,outMap2,outMap3,outMap4
!
      type(mqc_nef_set),target::full_NEFSet
      class(mqc_nef_set),pointer::ref_NEFSet
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSet => NEFSet
      if(NEFSet%isMirrorCell(1,2)) then
        full_NEFSet = mqc_nef_set_full_grid(NEFSet)
        ref_NEFSet => full_NEFSet
      end if
      select case(trim(NEFSet%Type()))
      case('vector')
        if(.not.(present(map_matrix_1).and.present(map_matrix_2))) call mqc_error(&
             'mqc_nef_set_expand: vector-type NEF_Set requires map_matrix_1 '//&
             'and map_matrix_2',6)
        if(NEFSet%NumSpace(1).ne.size(map_matrix_1)*size(map_matrix_2)) call mqc_error_i(&
             'mqc_nef_set_expand: vector cells'' NumSpace does not match '//&
             'size(map_matrix_1)*size(map_matrix_2) -- supplied maps do not '//&
             'match this NEF_Set''s actual per-cell content',6,&
             'NEFSet%NumSpace(1)',NEFSet%NumSpace(1),&
             'size(map_matrix_1)*size(map_matrix_2)',&
             size(map_matrix_1)*size(map_matrix_2))
!
        if(present(Order))  call mqc_error(&
             'mqc_nef_set_expand: vector-type should not recieve a Order input',6)
        if(present(spinTypeKet))  call mqc_error(&
             'mqc_nef_set_expand: vector-type should not recieve a spinTypeKet input',6)
!
        if(present(swapOrder)) then
          mySwap = swapOrder
        else
          mySwap = .false.
        end if
!
        if(present(spinType)) then
          mySpinType = spinType
        else
          mySpinType = 'all'
        end if
!          
  !     (matrix <-> vector: the dagger reverses the vector index, label 'asymmetric')
      call NEFSetOut%init(NEFSet%NumI,NEFSet%NumJ,sameBraKetSpace=NEFSet%isSameBraKetSpace())
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(ref_NEFSet%Set_Vector(I,J)%expand(&
                 map_matrix_1,map_matrix_2,spinType=mySpinType,swapOrder=mySwap),I,J)
          end do
        end do
        if(NEFSetOut%isOpen()) &
          call mqc_nef_set_allocate(NEFSetOut,'matrix',&
               Map_Matrix_1=map_matrix_1,Map_Matrix_2=map_matrix_2)
      case('matrix')
        if(.not.(present(map_matrix_1).and.present(map_matrix_2).and.&
                 present(map_matrix_3).and.present(map_matrix_4))) call mqc_error(&
             'mqc_nef_set_expand: matrix-type NEF_Set requires map_matrix_1, '//&
             'map_matrix_2, map_matrix_3 and map_matrix_4',6)
!
        if(present(swapOrder))  call mqc_error(&
             'mqc_nef_set_expand: matrix-type should not recieve a swapOrder input',6)
!
        if(present(Order)) then
          ia=Order(1) ; ja=Order(2) ; ka=Order(3) ; la=Order(4)
          myOrder = Order
        else
          ia=1 ; ja=2 ; ka=3 ; la=4
          myOrder = [1,2,3,4]
        end if
        if(ia.eq.1) then
          s1 = size(map_matrix_1)
        else if(ja.eq.1) then
          s1 = size(map_matrix_2)
        else if(ka.eq.1) then
          s1 = size(map_matrix_3)
        else
          s1 = size(map_matrix_4)
        end if
        if(ia.eq.2) then
          s2 = size(map_matrix_1)
        else if(ja.eq.2) then
          s2 = size(map_matrix_2)
        else if(ka.eq.2) then
          s2 = size(map_matrix_3)
        else
          s2 = size(map_matrix_4)
        end if
        if(ia.eq.3) then
          s3 = size(map_matrix_1)
        else if(ja.eq.3) then
          s3 = size(map_matrix_2)
        else if(ka.eq.3) then
          s3 = size(map_matrix_3)
        else
          s3 = size(map_matrix_4)
        end if
        if(ia.eq.4) then
          s4 = size(map_matrix_1)
        else if(ja.eq.4) then
          s4 = size(map_matrix_2)
        else if(ka.eq.4) then
          s4 = size(map_matrix_3)
        else
          s4 = size(map_matrix_4)
        end if
!
        if(present(spinTypeKet)) then
          mySpinTypeKet = spinTypeKet
        else
          mySpinTypeKet = 'all'
        end if
!          
        if(present(spinType)) then
          mySpinType = spinType
        else
          mySpinType = 'all'
        end if
!          
        expectedRow = s1*s2
        expectedCol = s3*s4
        if(NEFSet%NumSpace(1).ne.expectedRow.or.NEFSet%NumSpace(2).ne.expectedCol) &
          call mqc_error_i(&
             'mqc_nef_set_expand: matrix cells'' row/col NumSpace does not '//&
             'match the expected sizes implied by map_matrix_1..4 and Order '//&
             '-- supplied maps do not match this NEF_Set''s actual per-cell '//&
             'content',6,&
             'NEFSet%NumSpace(1)',NEFSet%NumSpace(1),'expectedRow',expectedRow,&
             'NEFSet%NumSpace(2)',NEFSet%NumSpace(2),'expectedCol',expectedCol)
        call NEFSetOut%init(NEFSet%NumI,NEFSet%NumJ,sameBraKetSpace=NEFSet%isSameBraKetSpace(),&
           symmetry=mqc_nef_pair_order_label(NEFSet%getSymm(),spinType,spinTypeKet,Order))
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(ref_NEFSet%Set_Matrix(I,J)%expand(&
                 map_matrix_1,map_matrix_2,map_matrix_3,map_matrix_4,&
                 spinType=mySpinType,spinTypeKet=mySpinTypeKet,Order=myOrder),I,J)
          end do
        end do
        if(NEFSetOut%isOpen()) then
          if(ia.eq.1) then
            outMap1 = map_matrix_1
          else if(ja.eq.1) then
            outMap1 = map_matrix_2
          else if(ka.eq.1) then
            outMap1 = map_matrix_3
          else
            outMap1 = map_matrix_4
          end if
          if(ia.eq.2) then
            outMap2 = map_matrix_1
          else if(ja.eq.2) then
            outMap2 = map_matrix_2
          else if(ka.eq.2) then
            outMap2 = map_matrix_3
          else
            outMap2 = map_matrix_4
          end if
          if(ia.eq.3) then
            outMap3 = map_matrix_1
          else if(ja.eq.3) then
            outMap3 = map_matrix_2
          else if(ka.eq.3) then
            outMap3 = map_matrix_3
          else
            outMap3 = map_matrix_4
          end if
          if(ia.eq.4) then
            outMap4 = map_matrix_1
          else if(ja.eq.4) then
            outMap4 = map_matrix_2
          else if(ka.eq.4) then
            outMap4 = map_matrix_3
          else
            outMap4 = map_matrix_4
          end if
          call mqc_nef_set_allocate(NEFSetOut,'r4',Map_Matrix_1=outMap1,&
               Map_Matrix_2=outMap2,Map_Matrix_3=outMap3,Map_Matrix_4=outMap4)
        end if
      case default
        call mqc_error_A('mqc_nef_set_expand: NEF_Set type does not support '//&
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
      class(mqc_nef_set),intent(in),target::NEFSet
      character(len=*),intent(in),optional::spinType,spinTypeKet
      logical,intent(in),optional::swapOrder
      integer(kind=int64),dimension(4),intent(in),optional::Order
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
      character(len=64)::mySpinType,mySpinTypeKet
      logical::mySwapOrder
      integer::row,col,id,outId,dimi,dimj,idij,ia,ja,ka,la,nRow,nCol
      Type(MQC_SCF_Eigenvalues),dimension(:),allocatable::outMap1,rowMap,colMap
!
      type(mqc_nef_set),target::full_NEFSet
      class(mqc_nef_set),pointer::ref_NEFSet
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSet => NEFSet
      if(NEFSet%isMirrorCell(1,2)) then
        full_NEFSet = mqc_nef_set_full_grid(NEFSet)
        ref_NEFSet => full_NEFSet
      end if
      select case(trim(NEFSet%Type()))
      case('matrix')
        if(present(spinTypeKet)) call mqc_error(&
             'mqc_nef_set_collapse: matrix-type should not recieve a spinTypeKet input',6)
        if(present(Order)) call mqc_error(&
             'mqc_nef_set_collapse: matrix-type should not recieve a Order input',6)
  !     (matrix <-> vector: the dagger reverses the vector index, label 'asymmetric')
      call NEFSetOut%init(NEFSet%NumI,NEFSet%NumJ,sameBraKetSpace=NEFSet%isSameBraKetSpace())
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(ref_NEFSet%Set_Matrix(I,J)%collapse(spinType=spinType,swapOrder=swapOrder),I,J)
          end do
        end do
        if(NEFSetOut%isOpen()) then
          if(present(spinType)) then
            mySpinType = spinType
          else
            mySpinType = 'all'
          end if
          if(present(swapOrder)) then
            mySwapOrder = swapOrder
          else
            mySwapOrder = .false.
          end if
          allocate(outMap1(NEFSet%NumSpace(1)*NEFSet%NumSpace(2)))
          select case(trim(mySpinType))
          case('singlet')
            do row = 1, NEFSet%NumSpace(1)
              do col = 1, NEFSet%NumSpace(2)
                id = (row-1)*NEFSet%NumSpace(2)+col
                if(mySwapOrder) then
                  outId = (col-1)*NEFSet%NumSpace(1)+row
                else
                  outId = id
                end if
                call outMap1(outId)%init(&
                     NEFSet%Size_vector_alpha_2(col)*NEFSet%Size_vector_alpha_1(row),0,&
                     NEFSet%Size_vector_beta_2(col)*NEFSet%Size_vector_beta_1(row))
              end do
            end do
          case('triplet')
            do row = 1, NEFSet%NumSpace(1)
              do col = 1, NEFSet%NumSpace(2)
                id = (row-1)*NEFSet%NumSpace(2)+col
                if(mySwapOrder) then
                  outId = (col-1)*NEFSet%NumSpace(1)+row
                else
                  outId = id
                end if
                call outMap1(outId)%init(&
                     NEFSet%Size_vector_alpha_2(col)*NEFSet%Size_vector_beta_1(row),0,&
                     NEFSet%Size_vector_beta_2(col)*NEFSet%Size_vector_alpha_1(row))
              end do
            end do
          case('all','full')
            do row = 1, NEFSet%NumSpace(1)
              do col = 1, NEFSet%NumSpace(2)
                id = (row-1)*NEFSet%NumSpace(2)+col
                if(mySwapOrder) then
                  outId = (col-1)*NEFSet%NumSpace(1)+row
                else
                  outId = id
                end if
                call outMap1(outId)%init(&
                     NEFSet%Size_vector_alpha_2(col)*NEFSet%Size_vector_alpha_1(row)+&
                     NEFSet%Size_vector_alpha_2(col)*NEFSet%Size_vector_beta_1(row),0,&
                     NEFSet%Size_vector_beta_2(col)*NEFSet%Size_vector_beta_1(row)+&
                     NEFSet%Size_vector_beta_2(col)*NEFSet%Size_vector_alpha_1(row))
              end do
            end do
          case default
            call mqc_error_a('mqc_nef_set_collapse: Unrecognized excitation type',6,&
                 'spinType',mySpinType)
          end select
          call mqc_nef_set_allocate(NEFSetOut,'vector',Map_Matrix_1=outMap1)
          deallocate(outMap1)
        end if
      case('r4')
        if(present(swapOrder)) call mqc_error(&
             'mqc_nef_set_collapse: r4-type should not recieve a swapOrder input',6)
        call NEFSetOut%init(NEFSet%NumI,NEFSet%NumJ,sameBraKetSpace=NEFSet%isSameBraKetSpace(),&
           symmetry=mqc_nef_pair_order_label(NEFSet%getSymm(),spinType,spinTypeKet,Order))
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(ref_NEFSet%Set_R4(I,J)%collapse(&
                 spinType=spinType,spinTypeKet=spinTypeKet,Order=Order),I,J)
          end do
        end do
        if(NEFSetOut%isOpen()) then
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
          if(present(Order)) then
            ia = int(Order(1)) ; ja = int(Order(2)) ; ka = int(Order(3)) ; la = int(Order(4))
          else
            ia = 1 ; ja = 2 ; ka = 3 ; la = 4
          end if
          nRow = axisNumSpace(ia)*axisNumSpace(ja)
          nCol = axisNumSpace(ka)*axisNumSpace(la)
          allocate(rowMap(nRow))
          allocate(colMap(nCol))
          select case(trim(mySpinType))
          case('singlet')
            do dimi = 1, axisNumSpace(ia)
              do dimj = 1, axisNumSpace(ja)
                idij = axisNumSpace(ja)*(dimi-1)+dimj
                call rowMap(idij)%init(&
                     axisAlphaSize(ia,dimi)*axisAlphaSize(ja,dimj),0,&
                     axisBetaSize(ia,dimi)*axisBetaSize(ja,dimj))
              end do
            end do
          case('triplet')
            do dimi = 1, axisNumSpace(ia)
              do dimj = 1, axisNumSpace(ja)
                idij = axisNumSpace(ja)*(dimi-1)+dimj
                call rowMap(idij)%init(&
                     axisAlphaSize(ia,dimi)*axisBetaSize(ja,dimj),0,&
                     axisBetaSize(ia,dimi)*axisAlphaSize(ja,dimj))
              end do
            end do
          case('all','full')
            do dimi = 1, axisNumSpace(ia)
              do dimj = 1, axisNumSpace(ja)
                idij = axisNumSpace(ja)*(dimi-1)+dimj
                call rowMap(idij)%init(&
                     axisAlphaSize(ia,dimi)*(axisBetaSize(ja,dimj)+axisAlphaSize(ja,dimj)),0,&
                     axisBetaSize(ia,dimi)*(axisAlphaSize(ja,dimj)+axisBetaSize(ja,dimj)))
              end do
            end do
          case default
            call mqc_error_a('mqc_nef_set_collapse: Unrecognized excitation type',6,&
                 'spinType',mySpinType)
          end select
          select case(trim(mySpinTypeKet))
          case('singlet')
            do dimi = 1, axisNumSpace(ka)
              do dimj = 1, axisNumSpace(la)
                idij = axisNumSpace(la)*(dimi-1)+dimj
                call colMap(idij)%init(&
                     axisAlphaSize(ka,dimi)*axisAlphaSize(la,dimj),0,&
                     axisBetaSize(ka,dimi)*axisBetaSize(la,dimj))
              end do
            end do
          case('triplet')
            do dimi = 1, axisNumSpace(ka)
              do dimj = 1, axisNumSpace(la)
                idij = axisNumSpace(la)*(dimi-1)+dimj
                call colMap(idij)%init(&
                     axisAlphaSize(ka,dimi)*axisBetaSize(la,dimj),0,&
                     axisBetaSize(ka,dimi)*axisAlphaSize(la,dimj))
              end do
            end do
          case('all','full')
            do dimi = 1, axisNumSpace(ka)
              do dimj = 1, axisNumSpace(la)
                idij = axisNumSpace(la)*(dimi-1)+dimj
                call colMap(idij)%init(&
                     axisAlphaSize(ka,dimi)*(axisBetaSize(la,dimj)+axisAlphaSize(la,dimj)),0,&
                     axisBetaSize(ka,dimi)*(axisAlphaSize(la,dimj)+axisBetaSize(la,dimj)))
              end do
            end do
          case default
            call mqc_error_a('mqc_nef_set_collapse: Unrecognized excitation type',6,&
                 'spinTypeKet',mySpinTypeKet)
          end select
          call mqc_nef_set_allocate(NEFSetOut,'matrix',&
               Map_Matrix_1=rowMap,Map_Matrix_2=colMap)
          deallocate(rowMap,colMap)
        end if
      case default
        call mqc_error_A('mqc_nef_set_collapse: NEF_Set kind does not support '//&
             'collapse (must be matrix or r4)',6,'Array_Type',NEFSet%Type())
      end select
!
      contains
!
        function axisAlphaSize(axis,idx) result(sz)
        implicit none
        integer,intent(in)::axis,idx
        integer::sz
        select case(axis)
        case(1); sz = NEFSet%Size_vector_alpha_1(idx)
        case(2); sz = NEFSet%Size_vector_alpha_2(idx)
        case(3); sz = NEFSet%Size_vector_alpha_3(idx)
        case(4); sz = NEFSet%Size_vector_alpha_4(idx)
        end select
        end function axisAlphaSize
!
        function axisBetaSize(axis,idx) result(sz)
        implicit none
        integer,intent(in)::axis,idx
        integer::sz
        select case(axis)
        case(1); sz = NEFSet%Size_vector_beta_1(idx)
        case(2); sz = NEFSet%Size_vector_beta_2(idx)
        case(3); sz = NEFSet%Size_vector_beta_3(idx)
        case(4); sz = NEFSet%Size_vector_beta_4(idx)
        end select
        end function axisBetaSize
!
        function axisNumSpace(axis) result(n)
        implicit none
        integer,intent(in)::axis
        integer::n
        n = NEFSet%NumSpace(axis)
        end function axisNumSpace
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
      character(len=16)::outSymm
      logical::outSBK
      type(mqc_scf_integral)::blockTmp
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
      type(mqc_nef_vector),target::cellVector
      type(mqc_nef_vector)::NEFVecOut
      type(mqc_scf_eigenvalues),pointer::blkCV
      type(mqc_nef_matrix)::cellMatrix,NEFMatOut
      type(mqc_nef_r4),target::cellR4
      type(mqc_nef_r4)::NEFR4Out
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
          outSymm = 'asymmetric'
          outSBK = .false.
          do I = 1, myNumI
            do J = 1, myNumJ
              if(NEFSet%isZero(I,J)) cycle
              cellVector = NEFSet%atVector(int(I),int(J))
              outerIdx = (I-1)*myNumJ+J
              do b = 1, origNumSpace
                if(cellVector%isBlockZero(b)) cycle
                blkCV => cellVector%iref(b)
                newIdx = (outerIdx-1)*origNumSpace+b
                call NEFVecOut%iput(blkCV,newIdx)
              end do
            end do
          end do
          call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
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
!         rows (I,u), columns J: bra x ket kept, no dagger pairing ('asymmetric')
          outSymm = 'asymmetric'
          outSBK = NEFSet%isSameBraKetSpace()
          call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2,sameBraKet=outSBK)
          do I = 1, myNumI
            do J = 1, myNumJ
              if(NEFSet%isZero(I,J)) cycle
              cellVector = NEFSet%atVector(int(I),int(J))
              do b = 1, origNumSpace
                if(cellVector%isBlockZero(b)) cycle
                blkCV => cellVector%iref(b)
                newRowId = (I-1)*origNumSpace+b
                call blockTmp%init(&
                     NEFSet%Size_vector_alpha_1(b),NEFSet%Size_vector_beta_1(b),&
                     nAlpha2=1,nBeta2=1)
                do ka = 1, NEFSet%Size_vector_alpha_1(b)
                  tmpScalar = blkCV%at(ka,'alpha')
                  call blockTmp%put(tmpScalar,ka,1_int64,'aa')
                end do
                do ka = 1, NEFSet%Size_vector_beta_1(b)
                  tmpScalar = blkCV%at(ka,'beta')
                  call blockTmp%put(tmpScalar,ka,1_int64,'bb')
                end do
                call NEFMatOut%iput(blockTmp,newRowId,J)
              end do
            end do
          end do
          call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
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
!       Grid -> object symmetry (rule G): (I,J)*(u,v) -> (Iu,Jv) keeps the set label and
!       gives a <Psi|O|Psi>-like matrix; a joint (I,J) index does not.
        outSymm = 'asymmetric'
        if(.not.myJointIJ) outSymm = NEFSet%getSymm()
        outSBK = (NEFSet%isSameBraKetSpace().and..not.myJointIJ)
        call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2,symmetry=outSymm,&
             sameBraKet=outSBK)
        do I = 1, myNumI
          do J = 1, myNumJ
            if(NEFSet%isZero(I,J)) cycle
            cellMatrix = NEFSet%atMatrix(int(I),int(J))
            G = (I-1)*myNumJ+J
            do r = 1, origNumSpace_row
              do c = 1, origNumSpace_col
                if(cellMatrix%isBlockZero(r,c)) cycle
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
                call NEFMatOut%iput(cellMatrix%iat(r,c),newRowId,newColId)
              end do
            end do
          end do
        end do
        call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
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
!
        allocate(New_Map_1(axisNumSpace(1)*foldGridSize(1)))
        do kk = 1, foldGridSize(1)
          do b1 = 1, axisNumSpace(1)
            call New_Map_1((kk-1)*axisNumSpace(1)+b1)%init(&
                 NEFSet%Size_vector_alpha_1(b1),0,lengthBeta=NEFSet%Size_vector_beta_1(b1))
          end do
        end do
!
        allocate(New_Map_2(axisNumSpace(2)*foldGridSize(2)))
        do kk = 1, foldGridSize(2)
          do b2 = 1, axisNumSpace(2)
            call New_Map_2((kk-1)*axisNumSpace(2)+b2)%init(&
                 NEFSet%Size_vector_alpha_2(b2),0,lengthBeta=NEFSet%Size_vector_beta_2(b2))
          end do
        end do
!
        allocate(New_Map_3(axisNumSpace(3)*foldGridSize(3)))
        do kk = 1, foldGridSize(3)
          do b3 = 1, axisNumSpace(3)
            call New_Map_3((kk-1)*axisNumSpace(3)+b3)%init(&
                 NEFSet%Size_vector_alpha_3(b3),0,lengthBeta=NEFSet%Size_vector_beta_3(b3))
          end do
        end do
!
        allocate(New_Map_4(axisNumSpace(4)*foldGridSize(4)))
        do kk = 1, foldGridSize(4)
          do b4 = 1, axisNumSpace(4)
            call New_Map_4((kk-1)*axisNumSpace(4)+b4)%init(&
                 NEFSet%Size_vector_alpha_4(b4),0,lengthBeta=NEFSet%Size_vector_beta_4(b4))
          end do
        end do
!
!       Grid -> object symmetry (rule G): I and J folded on the dagger-paired axes
!       (IAxisR4+JAxisR4 = 5) keep the set label; the merged R4 is <Psi|O|Psi>-like
!       (sameBraKet = sameBraKetSpace) unless I,J are folded jointly into one axis.
        outSymm = 'asymmetric'
        if(.not.myJointIJ.and.myIAxisR4+myJAxisR4.eq.5) outSymm = NEFSet%getSymm()
        outSBK = (NEFSet%isSameBraKetSpace().and..not.myJointIJ)
        call NEFR4Out%init(New_Map_1,New_Map_2,New_Map_3,New_Map_4,symmetry=outSymm,&
             sameBraKet=outSBK)
        do I = 1, myNumI
          do J = 1, myNumJ
            if(NEFSet%isZero(I,J)) cycle
            cellR4 = NEFSet%atR4(int(I),int(J))
!           mirror slots of a symmetric cell are read from a full-storage copy
            if(cellR4%Mirrored) cellR4 = mqc_nef_r4_full_storage(cellR4)
            G = (I-1)*myNumJ+J
            do b1 = 1, origNumSpace1
              do b2 = 1, origNumSpace2
                do b3 = 1, origNumSpace3
                  do b4 = 1, origNumSpace4
                    if(cellR4%isBlockZero(b1,b2,b3,b4)) cycle
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
                    call NEFR4Out%iput(cellR4%iref(b1,b2,b3,b4),nb1,nb2,nb3,nb4)
                  end do
                end do
              end do
            end do
          end do
        end do
        call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
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
!
!     PROCEDURE MQC_NEF_Set_Outer_Merge_Space
!
!
!>    \brief <b> MQC_NEF_Set_Outer_Merge_Space
!>     </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
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
      function mqc_nef_set_outer_merge_space(NEFSet,jointIJ,jointAxis,IAxisR4,JAxisR4) result(NEFSetOut)
!
      implicit none
      character(len=16)::outSymm
      logical::outSBK
      type(mqc_scf_integral)::blockTmp
      logical::blockStarted
      class(mqc_nef_set),intent(in)::NEFSet
      logical,intent(in),optional::jointIJ
      integer,intent(in),optional::jointAxis
      integer,intent(in),optional::IAxisR4,JAxisR4
      type(mqc_nef_set)::NEFSetOut
      logical::myJointIJ,started
      integer::myJointAxis,myIAxisR4,myJAxisR4,axisNum
      integer(kind=int64)::I,J,myNumI,myNumJ,G
      integer(kind=int64)::b,origNumSpace,outerIdx,outerIdxb,newIdx,ka
      integer(kind=int64)::r,c,origNumSpace_row,origNumSpace_col
      integer(kind=int64)::newRowId,newColId,origId,newId,newRowIdb,newColIdb
      integer(kind=int64)::b1,b2,b3,b4,nb1,nb2,nb3,nb4
      integer(kind=int64)::origNumSpace1,origNumSpace2,origNumSpace3,origNumSpace4
      integer(kind=int64)::axisNumSpace(4),foldGridSize(4),kk
      type(mqc_nef_vector),target::cellVector
      type(mqc_nef_vector)::NEFVecOut
      type(mqc_scf_eigenvalues),pointer::blkCV
      type(mqc_scf_eigenvalues)::blkE,blkEEmpty
      type(mqc_scf_integral)::blkM,blkMEmpty,cellBlk
      logical::vecStarted,matStarted
      type(mqc_nef_matrix)::cellMatrix,NEFMatOut
      type(mqc_nef_r4),target::cellR4
      type(mqc_nef_r4)::NEFR4Out
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
          allocate(New_Map_1(origNumSpace))
          do b = 1, origNumSpace
            call New_Map_1(b)%init(&
                 NEFSet%Size_vector_alpha_1(b)*myNumI*myNumJ,0,&
                 lengthBeta=NEFSet%Size_vector_beta_1(b)*myNumI*myNumJ)
          end do
          started=.false.
          call NEFVecOut%init(New_Map_1)
          outSymm = 'asymmetric'
          outSBK = .false.
          do b = 1, origNumSpace
            vecStarted = .false.
            do I = 1, myNumI
              do J = 1, myNumJ
                if(NEFSet%isZero(I,J)) cycle
                cellVector = NEFSet%atVector(I,J)
                if (cellVector%isBlockZero(b)) cycle
                blkCV => cellVector%iref(b)
                started=.true.
                if (.not.vecStarted) then
                   blkE = blkEEmpty
                   call blkE%init(NEFVecOut%Size_vector_alpha(b),&
                                      0.0,lengthBeta=NEFVecOut%Size_vector_beta(b))
                end if
                vecStarted = .true.
                !!! THIS BLOCK should be a vput in mqc_est !!!
                outerIdx = ((I-1)*myNumJ+J-1)*NEFSet%Size_vector_alpha_1(b)
                do b1 = 1, NEFSet%Size_vector_alpha_1(b)
                  call blkE%put(&
                       blkCV%at(b1,'alpha'),b1+outerIdx,'alpha')
                end do
                outerIdxb = ((I-1)*myNumJ+J-1)*NEFSet%Size_vector_beta_1(b)
                do b1 = 1, NEFSet%Size_vector_beta_1(b)
                  call blkE%put(&
                       blkCV%at(b1,'beta'),b1+outerIdxb,'beta')
                end do
                !!! !!!
              end do
            end do
            if (vecStarted) call NEFVecOut%iput(blkE,b)
          end do
          call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
          if (started) then 
            call NEFSetOut%put(NEFVecOut,1,1)
          else
            call mqc_nef_set_allocate(NEFSetOut,'vector',Map_Matrix_1=New_Map_1)
          end if
        else
          allocate(New_Map_1(origNumSpace))
          allocate(New_Map_2(1))
          do b = 1, origNumSpace
            call New_Map_1(b)%init(&
                            NEFSet%Size_vector_alpha_1(b)*myNumI,0,&
                 lengthBeta=NEFSet%Size_vector_beta_1(b)*myNumI)
          end do
          call New_Map_2(1)%init(myNumJ,0,lengthBeta=0)
!         rows (I,u), columns J: bra x ket kept, no dagger pairing ('asymmetric')
          outSymm = 'asymmetric'
          outSBK = NEFSet%isSameBraKetSpace()
          call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2,sameBraKet=outSBK)
          started=.false.
          do b = 1, origNumSpace
            blockStarted = .false.
            do I = 1, myNumI
              outerIdx = (I-1)*NEFSet%Size_vector_alpha_1(b)
              outerIdxb = (I-1)*NEFSet%Size_vector_beta_1(b)
              do J = 1, myNumJ
                if(NEFSet%isZero(I,J)) cycle
                cellVector = NEFSet%atVector(I,J)
                if(cellVector%isBlockZero(b)) cycle
                blkCV => cellVector%iref(b)
                started=.true.
                if (.not.blockStarted)&
                   call blockTmp%init(NEFMatOut%Size_vector_alpha_row(b),&
                                      NEFMatOut%Size_vector_beta_row(b),&
                                      0.0,'',myNumJ,0)
                blockStarted = .true.
                !!! THIS BLOCK should be a mput in mqc_est !!!
                do b1 = 1, NEFSet%Size_vector_alpha_1(b)
                  call blockTmp%put(&
                       blkCV%at(b1,'alpha'),b1+outerIdx,J,'alpha')
                end do
                do b1 = 1, NEFSet%Size_vector_beta_1(b)
                  call blockTmp%put(&
                       blkCV%at(b1,'beta'),b1+outerIdxb,J,'alphaBeta')
                end do
                !!! !!!
              end do
            end do
            if (blockStarted) call NEFMatOut%iput(blockTmp,b,1)
          end do
          call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
          if (started) then 
            call NEFSetOut%put(NEFMatOut,1,1)
          else
            call mqc_nef_set_allocate(NEFSetOut,'matrix',Map_Matrix_1=New_Map_1,Map_Matrix_2=New_Map_2)
          end if
        end if
!
      case('matrix')
        origNumSpace_row = NEFSet%NumSpace(1)
        origNumSpace_col = NEFSet%NumSpace(2)
        allocate(New_Map_1(origNumSpace_row))
        allocate(New_Map_2(origNumSpace_col))
        if(.not.myJointIJ) then
          do r = 1, origNumSpace_row
            call New_Map_1(r)%init(&
                            NEFSet%Size_vector_alpha_1(r)*myNumI,0,&
                 lengthBeta=NEFSet%Size_vector_beta_1(r)*myNumI)
          end do
          do c = 1, origNumSpace_col
            call New_Map_2(c)%init(&
                            NEFSet%Size_vector_alpha_2(c)*myNumJ,0,&
                 lengthBeta=NEFSet%Size_vector_beta_2(c)*myNumJ)
          end do
        else
          if(myJointAxis.ne.1.and.myJointAxis.ne.2) call mqc_error_i(&
               'mqc_nef_set_outer_collapse: jointAxis out-of-bounds for matrix',6,&
               'jointAxis',myJointAxis)
          if(myJointAxis.eq.1) then
            do r = 1, origNumSpace_row
              call New_Map_1(r)%init(&
                              NEFSet%Size_vector_alpha_1(r)*myNumI*myNumJ,0,&
                   lengthBeta=NEFSet%Size_vector_beta_1(r)*myNumI*myNumJ)
            end do
            do c = 1, origNumSpace_col
              call New_Map_2(c)%init(&
                   NEFSet%Size_vector_alpha_2(c),0,lengthBeta=NEFSet%Size_vector_beta_2(c))
            end do
          else
            do r = 1, origNumSpace_row
              call New_Map_1(r)%init(&
                   NEFSet%Size_vector_alpha_1(r),0,lengthBeta=NEFSet%Size_vector_beta_1(r))
            end do
            do c = 1, origNumSpace_col
              call New_Map_2(c)%init(&
                              NEFSet%Size_vector_alpha_2(c)*myNumI*myNumJ,0,&
                   lengthBeta=NEFSet%Size_vector_beta_2(c)*myNumI*myNumJ)
            end do
          end if
        end if
!       Grid -> object symmetry (rule G): I folded into the rows and J into the columns
!       keeps the set label and gives a <Psi|O|Psi>-like matrix; a joint index does not.
        outSymm = 'asymmetric'
        if(.not.myJointIJ) outSymm = NEFSet%getSymm()
        outSBK = (NEFSet%isSameBraKetSpace().and..not.myJointIJ)
        call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2,symmetry=outSymm,&
             sameBraKet=outSBK)
        started=.false.
        do r = 1, origNumSpace_row
          do c = 1, origNumSpace_col
            matStarted = .false.
            do I = 1, myNumI
              do J = 1, myNumJ
                if(NEFSet%isZero(I,J)) cycle
                cellMatrix = NEFSet%atMatrix(I,J)
                if(cellMatrix%isBlockZero(r,c)) cycle
                cellBlk = cellMatrix%iat(r,c)
                started=.true.
                if (.not.matStarted) then
                   blkM = blkMEmpty
                   call blkM%init(NEFMatOut%Size_vector_alpha_row(r),&
                                                           NEFMatOut%Size_vector_beta_row(r),&
                                                    0.0,'',NEFMatOut%Size_vector_alpha_col(c),&
                                                           NEFMatOut%Size_vector_beta_col(c))
                end if
                matStarted = .true.
                G = (I-1)*myNumJ+J
                if(.not.myJointIJ) then
                  newRowId = (I-1)*NEFSet%Size_vector_alpha_1(r)
                  newRowIdb = (I-1)*NEFSet%Size_vector_beta_1(r)
                  newColId = (J-1)*NEFSet%Size_vector_alpha_2(c)
                  newColIdb = (J-1)*NEFSet%Size_vector_beta_2(c)
                else if(myJointAxis.eq.1) then
                  newRowId = (G-1)*NEFSet%Size_vector_alpha_1(r)
                  newRowIdb = (G-1)*NEFSet%Size_vector_beta_1(r)
                  newColId = 0
                  newColIdb = 0
                else
                  newRowId = 0
                  newRowIdb = 0
                  newColId = (G-1)*NEFSet%Size_vector_alpha_2(c)
                  newColIdb = (G-1)*NEFSet%Size_vector_beta_2(c)
                end if
                do b1 = 1, NEFSet%Size_vector_alpha_1(r)
                  do b2 = 1, NEFSet%Size_vector_alpha_2(c)
                    call blkM%put(&
                        cellBlk%at(b1,b2,'alpha'),&
                           newRowId+b1,newColId+b2,'alpha')
                  end do
                end do
                do b1 = 1, NEFSet%Size_vector_beta_1(r)
                  do b2 = 1, NEFSet%Size_vector_beta_2(c)
                    call blkM%put(&
                        cellBlk%at(b1,b2,'beta'),&
                           newRowIdb+b1,newColIdb+b2,'beta')
                  end do
                end do
                if(cellBlk%Type().eq.'general') then
                  do b1 = 1, NEFSet%Size_vector_alpha_1(r)
                    do b2 = 1, NEFSet%Size_vector_beta_2(c)
                      call blkM%put(&
                          cellBlk%at(b1,b2,'betaalpha'),&
                             newRowId+b1,newColIdb+b2,'betaalpha')
                    end do
                  end do
                  do b1 = 1, NEFSet%Size_vector_beta_1(r)
                    do b2 = 1, NEFSet%Size_vector_alpha_2(c)
                      call blkM%put(&
                          cellBlk%at(b1,b2,'alphabeta'),&
                             newRowIdb+b1,newColId+b2,'alphabeta')
                    end do
                  end do
                end if
              end do
            end do
!           An all-zero block ('fullsymmetric') stays a zero block.
            if (matStarted) then
              if (blkM%getSymm().ne.'fullsymmetric') call NEFMatOut%iput(blkM,r,c)
            end if
          end do
        end do
        call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
        if (started) then 
          call NEFSetOut%put(NEFMatOut,1,1)
        else
          call mqc_nef_set_allocate(NEFSetOut,'matrix',Map_Matrix_1=New_Map_1,Map_Matrix_2=New_Map_2)
        end if
!
      case('r4')
        call mqc_error('mqc_nef_set_outer_merge_space: r4 case not yet implemented',6)
!
      case default
        call mqc_error_A('mqc_nef_set_outer_collapse: NEF_Set is open (holds no elements)',&
             6,'Array_Type',NEFSet%Type())
      end select
!
      end function mqc_nef_set_outer_merge_space
!
!
!     PROCEDURE MQC_NEF_Set_Outer_Unmerge_Space
!
!>    \brief <b> MQC_NEF_Set_Outer_Unmerge_Space undoes MQC_NEF_Set_Outer_Merge_Space</b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    MQC_NEF_Set_Outer_Unmerge_Space takes the single merged NEF_Vector/
!>    NEF_Matrix produced by MQC_NEF_Set_Outer_Merge_Space and splits it
!>    back into the original (NumI x NumJ) grid of per-cell objects, using
!>    a template NEF_Set (the original, pre-merge NEF_Set) to recover the
!>    grid size and the per-block alpha/beta lengths needed to slice the
!>    merged object apart. jointIJ/jointAxis must match whatever was passed
!>    to the forward merge call.
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] mergedSet
!>    \verbatim
!>        mergedSet is Class(MQC_NEF_Set)
!>        The 1x1 NEF_Set produced by MQC_NEF_Set_Outer_Merge_Space.
!>    \endverbatim
!>
!>    \param[in] template
!>    \verbatim
!>        template is Class(MQC_NEF_Set)
!>        The original, pre-merge NEF_Set (NumI, NumJ, NumSpace and the
!>        per-block Size_vector_* arrays are read from this object).
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_outer_unmerge_space(mergedSet,template,jointIJ,jointAxis) result(NEFSetOut)
!
      implicit none
      character(len=16)::outSymm
      logical::outSBK
      type(mqc_scf_integral)::mergedBlock,cellBlock
      class(mqc_nef_set),intent(in)::mergedSet,template
      logical,intent(in),optional::jointIJ
      integer,intent(in),optional::jointAxis
      type(mqc_nef_set)::NEFSetOut
      logical::myJointIJ
      integer::myJointAxis
      integer(kind=int64)::I,J,myNumI,myNumJ,G
      integer(kind=int64)::b,origNumSpace,outerIdxAlpha,outerIdxBeta,b1,b2
      integer(kind=int64)::r,c,origNumSpace_row,origNumSpace_col
      integer(kind=int64)::newRowId,newRowIdb,newColId,newColIdb,origId
      type(mqc_nef_vector)::cellVector
      type(mqc_nef_vector),target::mergedVec
      type(mqc_scf_eigenvalues),pointer::blkMV
      type(mqc_scf_eigenvalues)::blkE,blkEEmpty
      type(mqc_nef_matrix)::cellMatrix,mergedMat
      type(mqc_scalar)::tmpScalar
!
      if(present(jointIJ)) then
        myJointIJ = jointIJ
      else
        myJointIJ = (trim(template%Type()).eq.'vector')
      end if
      if(present(jointAxis)) then
        myJointAxis = jointAxis
      else
        myJointAxis = 1
      end if
!
      myNumI = template%NumI
      myNumJ = template%NumJ
!
      select case(trim(template%Type()))
      case('vector')
        origNumSpace = template%NumSpace(1)
        if(myJointIJ) then
!         joint (I,J) index: no bra x ket structure
          outSymm = 'asymmetric'
          outSBK = .false.
          call NEFSetOut%init(int(myNumI),int(myNumJ))
          mergedVec = mergedSet%atVector(1,1)
          do I = 1, myNumI
            do J = 1, myNumJ
              G = (I-1)*myNumJ+J
              call cellVector%init(template%Map_Matrix_1)
              do b = 1, origNumSpace
                if(mergedVec%isBlockZero(b)) cycle
                blkMV => mergedVec%iref(b)
                outerIdxAlpha = (G-1)*template%Size_vector_alpha_1(b)
                outerIdxBeta  = (G-1)*template%Size_vector_beta_1(b)
                blkE = blkEEmpty
                call blkE%init(template%Size_vector_alpha_1(b),&
                        0.0,lengthBeta=template%Size_vector_beta_1(b))
                do b1 = 1, template%Size_vector_alpha_1(b)
                  tmpScalar = blkMV%at(b1+outerIdxAlpha,'alpha')
                  call blkE%put(tmpScalar,b1,'alpha')
                end do
                do b1 = 1, template%Size_vector_beta_1(b)
                  tmpScalar = blkMV%at(b1+outerIdxBeta,'beta')
                  call blkE%put(tmpScalar,b1,'beta')
                end do
                call cellVector%iput(blkE,b)
              end do
              call NEFSetOut%put(cellVector,int(I),int(J))
            end do
          end do
        else
          mergedMat = mergedSet%atMatrix(1,1)
!         rows (I,u), columns J: the bra x ket structure follows the merged matrix
          outSymm = 'asymmetric'
          outSBK = (mergedMat%isSameBraKet().and.myNumI.eq.myNumJ)
          call NEFSetOut%init(int(myNumI),int(myNumJ),sameBraKetSpace=outSBK)
          do I = 1, myNumI
            do J = 1, myNumJ
              call cellVector%init(template%Map_Matrix_1,sameBraKet=(outSBK.and.I.eq.J))
              do b = 1, origNumSpace
                if(mergedMat%isBlockZero(b,1)) cycle
                mergedBlock = mergedMat%iat(b,1)
                outerIdxAlpha = (I-1)*template%Size_vector_alpha_1(b)
                outerIdxBeta = (I-1)*template%Size_vector_beta_1(b)
                blkE = blkEEmpty
                call blkE%init(template%Size_vector_alpha_1(b),&
                        0.0,lengthBeta=template%Size_vector_beta_1(b))
                do b1 = 1, template%Size_vector_alpha_1(b)
                  tmpScalar = mergedBlock%at(b1+outerIdxAlpha,J,'alpha')
                  call blkE%put(tmpScalar,b1,'alpha')
                end do
                do b1 = 1, template%Size_vector_beta_1(b)
                  tmpScalar = mergedBlock%at(b1+outerIdxBeta,J,'alphaBeta')
                  call blkE%put(tmpScalar,b1,'beta')
                end do
                call cellVector%iput(blkE,b)
              end do
              call NEFSetOut%put(cellVector,int(I),int(J))
            end do
          end do
        end if
!
      case('matrix')
        origNumSpace_row = template%NumSpace(1)
        origNumSpace_col = template%NumSpace(2)
        mergedMat = mergedSet%atMatrix(1,1)
!       Object -> grid symmetry (inverse of rule G): a non-joint merge keeps the label of
!       the merged matrix; its sameBraKet gives the bra x ket grid (square grid only).
        outSymm = 'asymmetric'
        if(.not.myJointIJ) outSymm = mergedMat%getSymm()
        outSBK = (mergedMat%isSameBraKet().and..not.myJointIJ.and.myNumI.eq.myNumJ)
        call NEFSetOut%init(int(myNumI),int(myNumJ),sameBraKetSpace=outSBK,symmetry=outSymm)
        do I = 1, myNumI
          do J = 1, myNumJ
            G = (I-1)*myNumJ+J
            call cellMatrix%init(template%Map_Matrix_1,initial_map_matrix_col=template%Map_Matrix_2,&
                 symmetry=outSymm,sameBraKet=(outSBK.and.I.eq.J))
            do r = 1, origNumSpace_row
              do c = 1, origNumSpace_col
                if(mergedMat%isBlockZero(r,c)) cycle
                mergedBlock = mergedMat%iat(r,c)
                if(.not.myJointIJ) then
                  newRowId  = (I-1)*template%Size_vector_alpha_1(r)
                  newRowIdb = (I-1)*template%Size_vector_beta_1(r)
                  newColId  = (J-1)*template%Size_vector_alpha_2(c)
                  newColIdb = (J-1)*template%Size_vector_beta_2(c)
                else if(myJointAxis.eq.1) then
                  newRowId  = (G-1)*template%Size_vector_alpha_1(r)
                  newRowIdb = (G-1)*template%Size_vector_beta_1(r)
                  newColId  = 0
                  newColIdb = 0
                else
                  newRowId  = 0
                  newRowIdb = 0
                  newColId  = (G-1)*template%Size_vector_alpha_2(c)
                  newColIdb = (G-1)*template%Size_vector_beta_2(c)
                end if
                call cellBlock%init(template%Size_vector_alpha_1(r),&
                     template%Size_vector_beta_1(r),0.0,'',template%Size_vector_alpha_2(c),&
                     template%Size_vector_beta_2(c))
                do b1 = 1, template%Size_vector_alpha_1(r)
                  do b2 = 1, template%Size_vector_alpha_2(c)
                    tmpScalar = mergedBlock%at(newRowId+b1,newColId+b2,'alpha')
                    call cellBlock%put(tmpScalar,b1,b2,'alpha')
                  end do
                end do
                do b1 = 1, template%Size_vector_beta_1(r)
                  do b2 = 1, template%Size_vector_beta_2(c)
                    tmpScalar = mergedBlock%at(newRowIdb+b1,newColIdb+b2,'beta')
                    call cellBlock%put(tmpScalar,b1,b2,'beta')
                  end do
                end do
                if(cellBlock%Type().eq.'general') then
                  do b1 = 1, template%Size_vector_alpha_1(r)
                    do b2 = 1, template%Size_vector_beta_2(c)
                      tmpScalar = mergedBlock%at(newRowId+b1,newColIdb+b2,'betaalpha')
                      call cellBlock%put(tmpScalar,b1,b2,'betaalpha')
                    end do
                  end do
                  do b1 = 1, template%Size_vector_beta_1(r)
                    do b2 = 1, template%Size_vector_alpha_2(c)
                      tmpScalar = mergedBlock%at(newRowIdb+b1,newColId+b2,'alphabeta')
                      call cellBlock%put(tmpScalar,b1,b2,'alphabeta')
                    end do
                  end do
                end if
                call cellMatrix%iput(cellBlock,r,c)
              end do
            end do
            call NEFSetOut%put(cellMatrix,int(I),int(J))
          end do
        end do
!
      case('r4')
        call mqc_error('mqc_nef_set_outer_unmerge_space: r4 case not yet implemented',6)
      case default
        call mqc_error_A('mqc_nef_set_outer_unmerge_space: template NEF_Set is open (holds no elements)',&
             6,'Array_Type',template%Type())
      end select
!
      end function mqc_nef_set_outer_unmerge_space
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
!>    MQC_NEF_Set_Outer_Expand inverts outerCollapse, covering every
!>    jointIJ/jointAxis/IAxisR4/JAxisR4 mode outerCollapse supports.
!>    outerCollapse discards NumI and NumJ when it folds everything into
!>    one combined (1,1) element, so those must be supplied back via
!>    templateSet (any MQC_NEF_Set with the right NumI/NumJ -- its own
!>    type, per-cell content, and whether any particular cell happens to
!>    be zero-flagged are all irrelevant and never inspected). Every
!>    other quantity needed to unfold NEFSet -- which type its cells
!>    were originally (vector/matrix/r4) and their per-axis block sizes
!>    -- is derived directly from NEFSet's own Type()/NumSpace/
!>    Size_vector_* fields instead, since those are always self-
!>    consistent with what NEFSet actually holds. The one mode
!>    outerCollapse supports that this does NOT invert is "vector,
!>    jointIJ=.false." (it produces a matrix-type combined object
!>    indistinguishable, by NEFSet%Type() alone, from a genuine
!>    matrix-origin fold) -- call this with jointIJ=.true. explicitly
!>    unset only for a genuinely vector-type NEFSet.
!>
!>    jointIJ/jointAxis/IAxisR4/JAxisR4 must be passed the SAME values
!>    used in the original outerCollapse call (defaults: jointIJ true
!>    when NEFSet is vector-type, false otherwise; jointAxis 1; IAxisR4
!>    1; JAxisR4 4).
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
!>    \param[in] templateSet
!>    \verbatim
!>        Type(MQC_NEF_Set)
!>        Any MQC_NEF_Set with the NumI/NumJ of the original, pre-
!>        outerCollapse grid -- that is all this is used for; its own
!>        type and content are never inspected.
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
      function mqc_nef_set_outer_expand(NEFSet,&
           templateSet,jointIJ,jointAxis,IAxisR4,JAxisR4) result(NEFSetOut)
!
      implicit none
      character(len=16)::outSymm
      logical::outSBK
      class(mqc_nef_set),intent(in)::NEFSet
      type(mqc_nef_set),intent(in)::templateSet
      logical,intent(in),optional::jointIJ
      integer,intent(in),optional::jointAxis,IAxisR4,JAxisR4
!
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
      type(mqc_nef_vector),target::combVector
      type(mqc_nef_vector)::cellVector
      type(mqc_nef_matrix)::combMatrix,cellMatrix
      type(mqc_nef_r4),target::combR4
      type(mqc_nef_r4)::cellR4
      type(mqc_scalar)::tmpScalar
      integer(kind=int64),dimension(4)::foldFactor
!
!     --- resolve grid shape (the ONLY thing templateSet is used for --
!     everything about per-cell structure below comes directly from
!     NEFSet's own NumSpace/Size_vector_* fields instead, since those
!     are always self-consistent with NEFSet's actual content
!     regardless of templateSet's own type/shape or whether any
!     particular cell of templateSet happens to be zero-flagged) ---
      myNumI = templateSet%NumI
      myNumJ = templateSet%NumJ
!
!     --- resolve mode selectors (same defaults as outerCollapse) ---
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
      select case(trim(NEFSet%Type()))
      case('vector')
        if(.not.myJointIJ) call mqc_error(&
             'mqc_nef_set_outer_expand: a vector-type NEFSet can only '//&
             'result from jointIJ=.true. folding (jointIJ=.false. for a '//&
             'vector origin instead produces a matrix-type combined '//&
             'object -- call this with that matrix-type NEFSet, not this '//&
             'vector-type one)',6)
        combVector = NEFSet%Set_Vector(1,1)
        if(mod(combVector%NumSpace,myNumI*myNumJ).ne.0) call mqc_error_i(&
             'mqc_nef_set_outer_expand: NEFSet''s combined Vector length is '//&
             'not evenly divisible by templateSet''s NumI*NumJ -- '//&
             'templateSet''s grid shape does not match this NEFSet',6,&
             'combVector%NumSpace',combVector%NumSpace,&
             'myNumI*myNumJ',myNumI*myNumJ)
        origNumSpace = combVector%NumSpace/(myNumI*myNumJ)
        mySizeAlpha1 = NEFSet%Size_vector_alpha_1(1:origNumSpace)
        mySizeBeta1 = NEFSet%Size_vector_beta_1(1:origNumSpace)
        allocate(cellMap1(origNumSpace))
        do b = 1, origNumSpace
          call cellMap1(b)%init(mySizeAlpha1(b),0,lengthBeta=mySizeBeta1(b))
        end do
        call NEFSetOut%init(int(myNumI),int(myNumJ))
        do I = 1, myNumI
          do J = 1, myNumJ
            outerIdx = (I-1)*myNumJ+J
            call cellVector%init(cellMap1)
            do b = 1, origNumSpace
              newIdx = (outerIdx-1)*origNumSpace+b
              if(combVector%isBlockZero(newIdx)) cycle
              call cellVector%iput(combVector%iref(newIdx),b)
            end do
            call NEFSetOut%put(cellVector,int(I),int(J))
          end do
        end do
!
      case('matrix')
        if(myJointIJ.and.myJointAxis.ne.1.and.myJointAxis.ne.2) call mqc_error_i(&
             'mqc_nef_set_outer_expand: jointAxis out-of-bounds for matrix',6,&
             'jointAxis',myJointAxis)
        combMatrix = NEFSet%Set_Matrix(1,1)
        if(.not.myJointIJ) then
          if(mod(combMatrix%NumSpace_row,myNumI).ne.0.or.&
             mod(combMatrix%NumSpace_col,myNumJ).ne.0) call mqc_error_i(&
               'mqc_nef_set_outer_expand: NEFSet''s combined Matrix row/col '//&
               'sizes are not evenly divisible by templateSet''s NumI/NumJ '//&
               '-- templateSet''s grid shape does not match this NEFSet',6,&
               'combMatrix%NumSpace_row',combMatrix%NumSpace_row,'myNumI',myNumI,&
               'combMatrix%NumSpace_col',combMatrix%NumSpace_col,'myNumJ',myNumJ)
          origNumSpace_row = combMatrix%NumSpace_row/myNumI
          origNumSpace_col = combMatrix%NumSpace_col/myNumJ
        else if(myJointAxis.eq.1) then
          if(mod(combMatrix%NumSpace_row,myNumI*myNumJ).ne.0) call mqc_error_i(&
               'mqc_nef_set_outer_expand: NEFSet''s combined Matrix row size '//&
               'is not evenly divisible by templateSet''s NumI*NumJ -- '//&
               'templateSet''s grid shape does not match this NEFSet',6,&
               'combMatrix%NumSpace_row',combMatrix%NumSpace_row,&
               'myNumI*myNumJ',myNumI*myNumJ)
          origNumSpace_row = combMatrix%NumSpace_row/(myNumI*myNumJ)
          origNumSpace_col = combMatrix%NumSpace_col
        else
          if(mod(combMatrix%NumSpace_col,myNumI*myNumJ).ne.0) call mqc_error_i(&
               'mqc_nef_set_outer_expand: NEFSet''s combined Matrix col size '//&
               'is not evenly divisible by templateSet''s NumI*NumJ -- '//&
               'templateSet''s grid shape does not match this NEFSet',6,&
               'combMatrix%NumSpace_col',combMatrix%NumSpace_col,&
               'myNumI*myNumJ',myNumI*myNumJ)
          origNumSpace_row = combMatrix%NumSpace_row
          origNumSpace_col = combMatrix%NumSpace_col/(myNumI*myNumJ)
        end if
        mySizeAlpha1 = NEFSet%Size_vector_alpha_1(1:origNumSpace_row)
        mySizeBeta1 = NEFSet%Size_vector_beta_1(1:origNumSpace_row)
        mySizeAlpha2 = NEFSet%Size_vector_alpha_2(1:origNumSpace_col)
        mySizeBeta2 = NEFSet%Size_vector_beta_2(1:origNumSpace_col)
        allocate(cellMap1(origNumSpace_row))
        do r = 1, origNumSpace_row
          call cellMap1(r)%init(mySizeAlpha1(r),0,lengthBeta=mySizeBeta1(r))
        end do
        allocate(cellMap2(origNumSpace_col))
        do c = 1, origNumSpace_col
          call cellMap2(c)%init(mySizeAlpha2(c),0,lengthBeta=mySizeBeta2(c))
        end do
!       Object -> grid symmetry (inverse of rule G): a non-joint fold keeps the label of
!       the combined matrix; its sameBraKet gives the bra x ket grid (square grid only).
        outSymm = 'asymmetric'
        if(.not.myJointIJ) outSymm = combMatrix%getSymm()
        outSBK = (combMatrix%isSameBraKet().and..not.myJointIJ.and.myNumI.eq.myNumJ)
        call NEFSetOut%init(int(myNumI),int(myNumJ),sameBraKetSpace=outSBK,symmetry=outSymm)
        do I = 1, myNumI
          do J = 1, myNumJ
            G = (I-1)*myNumJ+J
            call cellMatrix%init(cellMap1,initial_map_matrix_col=cellMap2,symmetry=outSymm,&
                 sameBraKet=(outSBK.and.I.eq.J))
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
                if(combMatrix%isBlockZero(newRowId,newColId)) cycle
                call cellMatrix%iput(combMatrix%iat(newRowId,newColId),r,c)
              end do
            end do
            call NEFSetOut%put(cellMatrix,int(I),int(J))
          end do
        end do
!
      case('r4')
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
!       mirror slots of a symmetric cell are read from a full-storage copy
        if(combR4%Mirrored) combR4 = mqc_nef_r4_full_storage(combR4)
        foldFactor = 1
        if(.not.myJointIJ) then
          foldFactor(myIAxisR4) = foldFactor(myIAxisR4)*myNumI
          foldFactor(myJAxisR4) = foldFactor(myJAxisR4)*myNumJ
        else
          foldFactor(myJointAxis) = myNumI*myNumJ
        end if
        if(mod(combR4%NumSpace_1,foldFactor(1)).ne.0.or.&
           mod(combR4%NumSpace_2,foldFactor(2)).ne.0) call mqc_error_i(&
             'mqc_nef_set_outer_expand: NEFSet''s combined R4 axis 1/2 sizes '//&
             'are not evenly divisible by the expected grid fold factor -- '//&
             'templateSet''s grid shape does not match this NEFSet',6,&
             'combR4%NumSpace_1',combR4%NumSpace_1,'foldFactor_1',foldFactor(1),&
             'combR4%NumSpace_2',combR4%NumSpace_2,'foldFactor_2',foldFactor(2))
        if(mod(combR4%NumSpace_3,foldFactor(3)).ne.0.or.&
           mod(combR4%NumSpace_4,foldFactor(4)).ne.0) call mqc_error_i(&
             'mqc_nef_set_outer_expand: NEFSet''s combined R4 axis 3/4 sizes '//&
             'are not evenly divisible by the expected grid fold factor -- '//&
             'templateSet''s grid shape does not match this NEFSet',6,&
             'combR4%NumSpace_3',combR4%NumSpace_3,'foldFactor_3',foldFactor(3),&
             'combR4%NumSpace_4',combR4%NumSpace_4,'foldFactor_4',foldFactor(4))
        origNumSpace1 = combR4%NumSpace_1/foldFactor(1)
        origNumSpace2 = combR4%NumSpace_2/foldFactor(2)
        origNumSpace3 = combR4%NumSpace_3/foldFactor(3)
        origNumSpace4 = combR4%NumSpace_4/foldFactor(4)
        mySizeAlpha1 = NEFSet%Size_vector_alpha_1(1:origNumSpace1)
        mySizeBeta1 = NEFSet%Size_vector_beta_1(1:origNumSpace1)
        mySizeAlpha2 = NEFSet%Size_vector_alpha_2(1:origNumSpace2)
        mySizeBeta2 = NEFSet%Size_vector_beta_2(1:origNumSpace2)
        mySizeAlpha3 = NEFSet%Size_vector_alpha_3(1:origNumSpace3)
        mySizeBeta3 = NEFSet%Size_vector_beta_3(1:origNumSpace3)
        mySizeAlpha4 = NEFSet%Size_vector_alpha_4(1:origNumSpace4)
        mySizeBeta4 = NEFSet%Size_vector_beta_4(1:origNumSpace4)
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
!       Object -> grid symmetry (inverse of rule G): I and J unfolded from the
!       dagger-paired axes (IAxisR4+JAxisR4 = 5) keep the label of the combined R4.
        outSymm = 'asymmetric'
        if(.not.myJointIJ.and.myIAxisR4+myJAxisR4.eq.5) outSymm = combR4%getSymm()
        outSBK = (combR4%isSameBraKet().and..not.myJointIJ.and.myNumI.eq.myNumJ)
        call NEFSetOut%init(int(myNumI),int(myNumJ),sameBraKetSpace=outSBK,symmetry=outSymm)
        do I = 1, myNumI
          do J = 1, myNumJ
            G = (I-1)*myNumJ+J
            call cellR4%init(cellMap1,cellMap2,cellMap3,cellMap4,symmetry=outSymm,&
                 sameBraKet=(outSBK.and.I.eq.J))
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
                    if(combR4%isBlockZero(nb1,nb2,nb3,nb4)) cycle
                    call cellR4%iput(combR4%iref(nb1,nb2,nb3,nb4),b1,b2,b3,b4)
                  end do
                end do
              end do
            end do
            call NEFSetOut%put(cellR4,int(I),int(J))
          end do
        end do
!
      case default
        call mqc_error_A('mqc_nef_set_outer_expand: NEFSet must hold '//&
             'vector, matrix, or r4 elements (it holds whatever '//&
             'outerCollapse actually produced)',6,'NEFSet%Type()',NEFSet%Type())
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
      character(len=16)::outSymm
      logical::outSBK
      type(mqc_scf_integral)::blockTmp
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
      type(mqc_nef_vector),target::cellVector
      type(mqc_nef_vector)::NEFVecOut
      type(mqc_scf_eigenvalues),pointer::blkCV
      type(mqc_nef_matrix)::cellMatrix,NEFMatOut,tmpMatIn
      type(mqc_nef_r4)::tmpR4In
      type(mqc_scf_eigenvalues),dimension(:),allocatable::New_Map_1,New_Map_2
      type(mqc_scalar)::tmpScalar
!
      myNumI = NEFSet%NumI
      myNumJ = NEFSet%NumJ
!
      select case(trim(NEFSet%Type()))
      case('vector')
        if(present(spinType)) call mqc_error(&
             'mqc_nef_set_full_collapse: vector-type should not recieve a spinType input',6)
        if(present(swapOrder)) call mqc_error(&
             'mqc_nef_set_full_collapse: vector-type should not recieve a swapOrder input',6)
        if(present(spinTypeKet)) call mqc_error(&
             'mqc_nef_set_full_collapse: vector-type should not recieve a spinTypeKet input',6)
        if(present(Order)) call mqc_error(&
             'mqc_nef_set_full_collapse: vector-type should not recieve a Order input',6)
        if(present(jointAxis)) call mqc_error(&
             'mqc_nef_set_full_collapse: vector-type should not recieve a jointAxis input',6)
        NEFSetOut = NEFSet%outerCollapse(jointIJ=jointIJ)
      case('matrix')
        if(present(spinTypeKet)) call mqc_error(&
             'mqc_nef_set_full_collapse: matrix-type should not recieve a spinTypeKet input',6)
        if(present(Order)) call mqc_error(&
             'mqc_nef_set_full_collapse: matrix-type should not recieve a Order input',6)
        if(present(jointIJ)) then
          myJointIJ = jointIJ
        else
          myJointIJ = .true.
        end if
        tmpMatIn = NEFSet%atMatrix(1,1)
        cellVector = tmpMatIn%collapse(spinType=spinType,swapOrder=swapOrder)
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
          outSymm = 'asymmetric'
          outSBK = .false.
          do I = 1, myNumI
            do J = 1, myNumJ
              tmpMatIn = NEFSet%atMatrix(I,J)
              cellVector = tmpMatIn%collapse(spinType=spinType,swapOrder=swapOrder)
              outerIdx = (I-1)*myNumJ+J
              do b = 1, origNumSpace
                if(cellVector%isBlockZero(b)) cycle
                blkCV => cellVector%iref(b)
                newIdx = (outerIdx-1)*origNumSpace+b
                call NEFVecOut%iput(blkCV,newIdx)
              end do
            end do
          end do
          call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
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
!         matrix cells collapse to vectors (label 'asymmetric'); rows (I,u), columns J
          outSymm = 'asymmetric'
          outSBK = NEFSet%isSameBraKetSpace()
          call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2,sameBraKet=outSBK)
          do I = 1, myNumI
            do J = 1, myNumJ
              tmpMatIn = NEFSet%atMatrix(I,J)
              cellVector = tmpMatIn%collapse(spinType=spinType,swapOrder=swapOrder)
              do b = 1, origNumSpace
                if(cellVector%isBlockZero(b)) cycle
                blkCV => cellVector%iref(b)
                newRowId = (I-1)*origNumSpace+b
                call blockTmp%init(&
                     cellVector%Size_vector_alpha(b),cellVector%Size_vector_beta(b),&
                     nAlpha2=1,nBeta2=1)
                do ka = 1, cellVector%Size_vector_alpha(b)
                  tmpScalar = blkCV%at(ka,'alpha')
                  call blockTmp%put(tmpScalar,ka,1_int64,'aa')
                end do
                do ka = 1, cellVector%Size_vector_beta(b)
                  tmpScalar = blkCV%at(ka,'beta')
                  call blockTmp%put(tmpScalar,ka,1_int64,'bb')
                end do
                call NEFMatOut%iput(blockTmp,newRowId,J)
              end do
            end do
          end do
          call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
          call NEFSetOut%put(NEFMatOut,1,1)
        end if
      case('r4')
        if(present(swapOrder)) call mqc_error(&
             'mqc_nef_set_full_collapse: r4-type should not recieve a swapOrder input',6)
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
        tmpR4In = NEFSet%atR4(1,1)
        cellMatrix = tmpR4In%collapse(&
             spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
        origNumSpace_row = cellMatrix%NumSpace_row
        origNumSpace_col = cellMatrix%NumSpace_col
        if(.not.myJointIJ) then
          allocate(New_Map_1(myNumI*origNumSpace_row))
          do I = 1, myNumI
            tmpR4In = NEFSet%atR4(I,1)
            cellMatrix = tmpR4In%collapse(&
                 spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
            do r = 1, origNumSpace_row
              call New_Map_1((I-1)*origNumSpace_row+r)%init(&
                   cellMatrix%Size_vector_alpha_row(r),0,lengthBeta=cellMatrix%Size_vector_beta_row(r))
            end do
          end do
          allocate(New_Map_2(myNumJ*origNumSpace_col))
          do J = 1, myNumJ
            tmpR4In = NEFSet%atR4(1,J)
            cellMatrix = tmpR4In%collapse(&
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
                tmpR4In = NEFSet%atR4(I,J)
                cellMatrix = tmpR4In%collapse(&
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
                tmpR4In = NEFSet%atR4(I,J)
                cellMatrix = tmpR4In%collapse(&
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
!       R4 cells collapse to matrices (label by the Order rule); I,J folded into rows and
!       columns keep it and give a <Psi|O|Psi>-like matrix (rule G); a joint index does not.
        outSymm = 'asymmetric'
        if(.not.myJointIJ) outSymm = mqc_nef_pair_order_label(NEFSet%getSymm(),spinType,&
                                     spinTypeKet,Order)
        outSBK = (NEFSet%isSameBraKetSpace().and..not.myJointIJ)
        call NEFMatOut%init(New_Map_1,initial_map_matrix_col=New_Map_2,symmetry=outSymm,&
             sameBraKet=outSBK)
        do I = 1, myNumI
          do J = 1, myNumJ
            tmpR4In = NEFSet%atR4(I,J)
            cellMatrix = tmpR4In%collapse(&
                 spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
            G = (I-1)*myNumJ+J
            do r = 1, origNumSpace_row
              do c = 1, origNumSpace_col
                if(cellMatrix%isBlockZero(r,c)) cycle
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
                call NEFMatOut%iput(cellMatrix%iat(r,c),newRowId,newColId)
              end do
            end do
          end do
        end do
        call NEFSetOut%init(1,1,sameBraKetSpace=outSBK,symmetry=outSymm)
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
!>    templateSet%Type() gives the type the ORIGINAL (pre-fullCollapse)
!>    grid cells held: 'matrix' (whose fullCollapse intermediate was
!>    Vector) or 'r4' (whose fullCollapse intermediate was Matrix) --
!>    the same restriction collapse/fullCollapse themselves have.
!>
!>    outerExpand derives all its per-cell sizing directly from NEFSet's
!>    own NumSpace/Size_vector_* fields (always self-consistent with
!>    what NEFSet actually holds, regardless of any other object's type
!>    or content), so templateSet is only needed here for its NumI/NumJ
!>    (the outer grid shape) at the outerExpand stage, and for its cell
!>    (1,1)'s own real axis maps (read directly, forwarded as the final
!>    expand stage's target maps) -- templateSet is never collapsed or
!>    otherwise required to match NEFSet's type. Call outerExpand and
!>    expand separately if explicit (non-template) control over either
!>    stage is needed.
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
!>        -- its Type() ('matrix' or 'r4') tells this function what type
!>        the original grid cells held before fullCollapse (no separate
!>        origType argument is needed); its own cell (1,1)'s real
!>        Map_Matrix_Row/Col (matrix-type) or Map_Matrix_1..4 (r4-type)
!>        are read directly and forwarded as the final expand stage's
!>        target maps (no separate map_matrix_1..4 argument is needed
!>        either); and its NumI/NumJ drive the outerExpand stage (see
!>        Purpose).
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
!>    \param[in] spinType
!>    \verbatim
!>        Character(len=*),Optional
!>        Forwarded to the final expand stage. Must match the value used
!>        in the original fullCollapse call.
!>    \endverbatim
!>
!>    \param[in] spinTypeKet,Order
!>    \verbatim
!>        Forwarded to the final expand stage. Must match the values
!>        used in the original fullCollapse call.
!>    \endverbatim
!
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nef_set_full_expand(NEFSet,templateSet,jointIJ,jointAxis,&
           spinType,spinTypeKet,Order,swapOrder) result(NEFSetOut)
!
      implicit none
      class(mqc_nef_set),intent(in)::NEFSet
      type(mqc_nef_set),intent(in)::templateSet
      logical,intent(in),optional::jointIJ
      integer,intent(in),optional::jointAxis
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      logical,intent(in),optional::swapOrder
      type(mqc_nef_set)::NEFSetOut
      type(mqc_nef_set)::intermediateSet
      type(mqc_nef_matrix)::tmpMat
      type(mqc_nef_r4)::tmpR4
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable::tMap1,tMap2,tMap3,tMap4
!
!
      select case(trim(templateSet%Type()))
      case('matrix')
        if(present(spinTypeKet)) call mqc_error(&
             'mqc_nef_set_full_expand: matrix templateSet should not recieve a spinTypeKet input',6)
        if(present(Order)) call mqc_error(&
             'mqc_nef_set_full_expand: matrix templateSet should not recieve a Order input',6)
        tmpMat = templateSet%atMatrix(1,1)
        tMap1 = mqc_nef_matrix_get_map(tmpMat,1)
        tMap2 = mqc_nef_matrix_get_map(tmpMat,2)
      case('r4')
        if(present(swapOrder)) call mqc_error(&
             'mqc_nef_set_full_expand: r4 templateSet should not recieve a swapOrder input',6)
        tmpR4 = templateSet%atR4(1,1)
        tMap1 = mqc_nef_r4_map(tmpR4,1)
        tMap2 = mqc_nef_r4_map(tmpR4,2)
        tMap3 = mqc_nef_r4_map(tmpR4,3)
        tMap4 = mqc_nef_r4_map(tmpR4,4)
      case default
        call mqc_error_A('mqc_nef_set_full_expand: templateSet must hold '//&
             'matrix or r4 elements',6,'Array_Type',templateSet%Type())
      end select
!
      intermediateSet = mqc_nef_set_outer_expand(NEFSet,&
                       templateSet,jointIJ=jointIJ,jointAxis=jointAxis)
      if(trim(templateSet%Type()).eq.'matrix') then
        NEFSetOut = mqc_nef_set_expand(intermediateSet,&
             map_matrix_1=tMap1,map_matrix_2=tMap2,&
             spinType=spinType,spinTypeKet=spinTypeKet,swapOrder=swapOrder)
      else
        NEFSetOut = mqc_nef_set_expand(intermediateSet,&
             map_matrix_1=tMap1,map_matrix_2=tMap2,&
             map_matrix_3=tMap3,map_matrix_4=tMap4,&
             spinType=spinType,spinTypeKet=spinTypeKet,Order=Order)
      end if
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
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ,sameBraKetSpace=NEFSetB%isSameBraKetSpace(),&
           symmetry=mqc_nef_scalar_symmetry(NEFSetB%getSymm(),scalar))
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
      else
        select case(trim(NEFSetB%Type()))
        case('vector')
          call mqc_nef_set_allocate(NEFSetOut,'vector',&
               Map_Matrix_1=NEFSetB%Map_Matrix_1)
        case('matrix')
          call mqc_nef_set_allocate(NEFSetOut,'matrix',&
               Map_Matrix_1=NEFSetB%Map_Matrix_1,&
               Map_Matrix_2=NEFSetB%Map_Matrix_2)
        case('r4')
          call mqc_nef_set_allocate(NEFSetOut,'r4',&
               Map_Matrix_1=NEFSetB%Map_Matrix_1,&
               Map_Matrix_2=NEFSetB%Map_Matrix_2,&
               Map_Matrix_3=NEFSetB%Map_Matrix_3,&
               Map_Matrix_4=NEFSetB%Map_Matrix_4)
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
      call NEFSetOut%init(NEFSet%NumJ,NEFSet%NumI,sameBraKetSpace=NEFSet%isSameBraKetSpace(),&
           symmetry=NEFSet%getSymm())
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
!     PROCEDURE MQC_NEFSet_Conjugate
!
!
!>    \brief <b> MQC_NEFSet_Conjugate computes the conjugate of a
!>    MQC_NEF_Set: the outer (I,J) grid is left unchanged and each cell
!>    is replaced by its own Conjg. </b>
!>
!>    \par Purpose:
!     =============
!>
!>    \verbatim
!>
!>    For vector-kind sets: NEFSetOut(I,J) = Conjg(NEFSet(I,J)), using
!>    MQC_NEF_Vector_Conjugate. For matrix-kind sets:
!>    NEFSetOut(I,J) = Conjg(NEFSet(I,J)). R4-kind sets are not yet
!>    supported (no plain Conjg primitive exists for MQC_NEF_R4).
!>
!>    \endverbatim
!
!     Arguments:
!     ==========
!>    \param[in] NEFSet
!>    \verbatim
!>        NEFSet is Type(MQC_NEF_Set)
!>        The (I,J) grid to conjugate.
!>    \endverbatim
!
!     Authors:
!     ========
!>    \author M. M. F. Moraes
!>    \date 2026
!
      function mqc_nefset_conjugate(NEFSet) result(NEFSetOut)
!
      implicit none
      type(mqc_nef_set),intent(in)::NEFSet
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      call NEFSetOut%init(NEFSet%NumI,NEFSet%NumJ,sameBraKetSpace=NEFSet%isSameBraKetSpace(),&
           symmetry=NEFSet%getSymm())
!
      select case(trim(NEFSet%Type()))
      case('vector')
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(conjg(NEFSet%atVector(I,J)),I,J)
          end do
        end do
      case('matrix')
        do I = 1, NEFSet%NumI
          do J = 1, NEFSet%NumJ
            if(NEFSet%isZero(I,J)) cycle
            call NEFSetOut%put(conjg(NEFSet%atMatrix(I,J)),I,J)
          end do
        end do
      case default
        call mqc_error_A('mqc_nefset_conjugate: NEF_Set kind not supported '//&
             '(only vector/matrix implemented)',6,'Array_Type',NEFSet%Type())
      end select
!
      end function mqc_nefset_conjugate
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
      call NEFSetOut%init(NEFSet%NumJ,NEFSet%NumI,sameBraKetSpace=NEFSet%isSameBraKetSpace(),&
           symmetry=NEFSet%getSymm())
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
      type(mqc_nef_set),intent(in),target::NEFSetA
      type(mqc_nef_vector),intent(in)::NEFVecB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetA
      class(mqc_nef_set),pointer::ref_NEFSetA
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetA => NEFSetA
      if(NEFSetA%isMirrorCell(1,2)) then
        full_NEFSetA = mqc_nef_set_full_grid(NEFSetA)
        ref_NEFSetA => full_NEFSetA
      end if
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'vector') call mqc_error_A(&
           'mqc_nefset_nefvector_sum: NEF_Set must hold vector elements (or be open)',&
           6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ,&
           sameBraKetSpace=(NEFSetA%isSameBraKetSpace().and.NEFVecB%isSameBraKet()),&
           symmetry=merge(NEFSetA%getSymm(),'asymmetric      ',NEFSetA%getSymm().eq.NEFVecB%getSymm()))
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call mqc_nef_set_put_built_vector(NEFSetOut,NEFSetA%atVector(I,J)+NEFVecB,I,J)
          else
            call mqc_nef_set_put_built_vector(NEFSetOut,ref_NEFSetA%Set_Vector(I,J)+NEFVecB,I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetB
      class(mqc_nef_set),pointer::ref_NEFSetB
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetB => NEFSetB
      if(NEFSetB%isMirrorCell(1,2)) then
        full_NEFSetB = mqc_nef_set_full_grid(NEFSetB)
        ref_NEFSetB => full_NEFSetB
      end if
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'vector') call mqc_error_A(&
           'mqc_nefvector_nefset_sum: NEF_Set must hold vector elements (or be open)',&
           6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ,&
           sameBraKetSpace=(NEFSetB%isSameBraKetSpace().and.NEFVecA%isSameBraKet()),&
           symmetry=merge(NEFSetB%getSymm(),'asymmetric      ',NEFSetB%getSymm().eq.NEFVecA%getSymm()))
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call mqc_nef_set_put_built_vector(NEFSetOut,NEFVecA+NEFSetB%atVector(I,J),I,J)
          else
            call mqc_nef_set_put_built_vector(NEFSetOut,NEFVecA+ref_NEFSetB%Set_Vector(I,J),I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetA
      type(mqc_nef_vector),intent(in)::NEFVecB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetA
      class(mqc_nef_set),pointer::ref_NEFSetA
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetA => NEFSetA
      if(NEFSetA%isMirrorCell(1,2)) then
        full_NEFSetA = mqc_nef_set_full_grid(NEFSetA)
        ref_NEFSetA => full_NEFSetA
      end if
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'vector') call mqc_error_A(&
           'mqc_nefset_nefvector_difference: NEF_Set must hold vector elements '//&
           '(or be open)',6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ,&
           sameBraKetSpace=(NEFSetA%isSameBraKetSpace().and.NEFVecB%isSameBraKet()),&
           symmetry=merge(NEFSetA%getSymm(),'asymmetric      ',NEFSetA%getSymm().eq.NEFVecB%getSymm()))
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call mqc_nef_set_put_built_vector(NEFSetOut,NEFSetA%atVector(I,J)-NEFVecB,I,J)
          else
            call mqc_nef_set_put_built_vector(NEFSetOut,ref_NEFSetA%Set_Vector(I,J)-NEFVecB,I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetB
      class(mqc_nef_set),pointer::ref_NEFSetB
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetB => NEFSetB
      if(NEFSetB%isMirrorCell(1,2)) then
        full_NEFSetB = mqc_nef_set_full_grid(NEFSetB)
        ref_NEFSetB => full_NEFSetB
      end if
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'vector') call mqc_error_A(&
           'mqc_nefvector_nefset_difference: NEF_Set must hold vector elements '//&
           '(or be open)',6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ,&
           sameBraKetSpace=(NEFSetB%isSameBraKetSpace().and.NEFVecA%isSameBraKet()),&
           symmetry=merge(NEFSetB%getSymm(),'asymmetric      ',NEFSetB%getSymm().eq.NEFVecA%getSymm()))
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call mqc_nef_set_put_built_vector(NEFSetOut,NEFVecA-NEFSetB%atVector(I,J),I,J)
          else
            call mqc_nef_set_put_built_vector(NEFSetOut,NEFVecA-ref_NEFSetB%Set_Vector(I,J),I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetA
      class(mqc_nef_set),pointer::ref_NEFSetA
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetA => NEFSetA
      if(NEFSetA%isMirrorCell(1,2)) then
        full_NEFSetA = mqc_nef_set_full_grid(NEFSetA)
        ref_NEFSetA => full_NEFSetA
      end if
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefset_nefmatrix_sum: NEF_Set must hold matrix elements (or be open)',&
           6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ,&
           sameBraKetSpace=(NEFSetA%isSameBraKetSpace().and.NEFMatB%isSameBraKet()),&
           symmetry=merge(NEFSetA%getSymm(),'asymmetric      ',NEFSetA%getSymm().eq.NEFMatB%getSymm()))
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call mqc_nef_set_put_built_matrix(NEFSetOut,NEFSetA%atMatrix(I,J)+NEFMatB,I,J)
          else
            call mqc_nef_set_put_built_matrix(NEFSetOut,ref_NEFSetA%Set_Matrix(I,J)+NEFMatB,I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetB
      class(mqc_nef_set),pointer::ref_NEFSetB
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetB => NEFSetB
      if(NEFSetB%isMirrorCell(1,2)) then
        full_NEFSetB = mqc_nef_set_full_grid(NEFSetB)
        ref_NEFSetB => full_NEFSetB
      end if
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefmatrix_nefset_sum: NEF_Set must hold matrix elements (or be open)',&
           6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ,&
           sameBraKetSpace=(NEFSetB%isSameBraKetSpace().and.NEFMatA%isSameBraKet()),&
           symmetry=merge(NEFSetB%getSymm(),'asymmetric      ',NEFSetB%getSymm().eq.NEFMatA%getSymm()))
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call mqc_nef_set_put_built_matrix(NEFSetOut,NEFMatA+NEFSetB%atMatrix(I,J),I,J)
          else
            call mqc_nef_set_put_built_matrix(NEFSetOut,NEFMatA+ref_NEFSetB%Set_Matrix(I,J),I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetA
      type(mqc_nef_matrix),intent(in)::NEFMatB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetA
      class(mqc_nef_set),pointer::ref_NEFSetA
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetA => NEFSetA
      if(NEFSetA%isMirrorCell(1,2)) then
        full_NEFSetA = mqc_nef_set_full_grid(NEFSetA)
        ref_NEFSetA => full_NEFSetA
      end if
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefset_nefmatrix_difference: NEF_Set must hold matrix elements '//&
           '(or be open)',6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ,&
           sameBraKetSpace=(NEFSetA%isSameBraKetSpace().and.NEFMatB%isSameBraKet()),&
           symmetry=merge(NEFSetA%getSymm(),'asymmetric      ',NEFSetA%getSymm().eq.NEFMatB%getSymm()))
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call mqc_nef_set_put_built_matrix(NEFSetOut,NEFSetA%atMatrix(I,J)-NEFMatB,I,J)
          else
            call mqc_nef_set_put_built_matrix(NEFSetOut,ref_NEFSetA%Set_Matrix(I,J)-NEFMatB,I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetB
      class(mqc_nef_set),pointer::ref_NEFSetB
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetB => NEFSetB
      if(NEFSetB%isMirrorCell(1,2)) then
        full_NEFSetB = mqc_nef_set_full_grid(NEFSetB)
        ref_NEFSetB => full_NEFSetB
      end if
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'matrix') call mqc_error_A(&
           'mqc_nefmatrix_nefset_difference: NEF_Set must hold matrix elements '//&
           '(or be open)',6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ,&
           sameBraKetSpace=(NEFSetB%isSameBraKetSpace().and.NEFMatA%isSameBraKet()),&
           symmetry=merge(NEFSetB%getSymm(),'asymmetric      ',NEFSetB%getSymm().eq.NEFMatA%getSymm()))
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call mqc_nef_set_put_built_matrix(NEFSetOut,NEFMatA-NEFSetB%atMatrix(I,J),I,J)
          else
            call mqc_nef_set_put_built_matrix(NEFSetOut,NEFMatA-ref_NEFSetB%Set_Matrix(I,J),I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetA
      type(mqc_nef_r4),intent(in)::NEFR4B
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetA
      class(mqc_nef_set),pointer::ref_NEFSetA
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetA => NEFSetA
      if(NEFSetA%isMirrorCell(1,2)) then
        full_NEFSetA = mqc_nef_set_full_grid(NEFSetA)
        ref_NEFSetA => full_NEFSetA
      end if
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'r4') call mqc_error_A(&
           'mqc_nefset_nefr4_sum: NEF_Set must hold r4 elements (or be open)',&
           6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ,&
           sameBraKetSpace=(NEFSetA%isSameBraKetSpace().and.NEFR4B%isSameBraKet()),&
           symmetry=merge(NEFSetA%getSymm(),'asymmetric      ',NEFSetA%getSymm().eq.NEFR4B%getSymm()))
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call mqc_nef_set_put_built_r4(NEFSetOut,NEFSetA%atR4(I,J)+NEFR4B,I,J)
          else
            call mqc_nef_set_put_built_r4(NEFSetOut,ref_NEFSetA%Set_R4(I,J)+NEFR4B,I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetB
      class(mqc_nef_set),pointer::ref_NEFSetB
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetB => NEFSetB
      if(NEFSetB%isMirrorCell(1,2)) then
        full_NEFSetB = mqc_nef_set_full_grid(NEFSetB)
        ref_NEFSetB => full_NEFSetB
      end if
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'r4') call mqc_error_A(&
           'mqc_nefr4_nefset_sum: NEF_Set must hold r4 elements (or be open)',&
           6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ,&
           sameBraKetSpace=(NEFSetB%isSameBraKetSpace().and.NEFR4A%isSameBraKet()),&
           symmetry=merge(NEFSetB%getSymm(),'asymmetric      ',NEFSetB%getSymm().eq.NEFR4A%getSymm()))
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call mqc_nef_set_put_built_r4(NEFSetOut,NEFR4A+NEFSetB%atR4(I,J),I,J)
          else
            call mqc_nef_set_put_built_r4(NEFSetOut,NEFR4A+ref_NEFSetB%Set_R4(I,J),I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetA
      type(mqc_nef_r4),intent(in)::NEFR4B
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetA
      class(mqc_nef_set),pointer::ref_NEFSetA
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetA => NEFSetA
      if(NEFSetA%isMirrorCell(1,2)) then
        full_NEFSetA = mqc_nef_set_full_grid(NEFSetA)
        ref_NEFSetA => full_NEFSetA
      end if
      if(.not.NEFSetA%isOpen().and.trim(NEFSetA%Type()).ne.'r4') call mqc_error_A(&
           'mqc_nefset_nefr4_difference: NEF_Set must hold r4 elements (or be open)',&
           6,'Array_Type',NEFSetA%Type())
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ,&
           sameBraKetSpace=(NEFSetA%isSameBraKetSpace().and.NEFR4B%isSameBraKet()),&
           symmetry=merge(NEFSetA%getSymm(),'asymmetric      ',NEFSetA%getSymm().eq.NEFR4B%getSymm()))
      do I = 1, NEFSetA%NumI
        do J = 1, NEFSetA%NumJ
          if(NEFSetA%isZero(I,J)) then
            call mqc_nef_set_put_built_r4(NEFSetOut,NEFSetA%atR4(I,J)-NEFR4B,I,J)
          else
            call mqc_nef_set_put_built_r4(NEFSetOut,ref_NEFSetA%Set_R4(I,J)-NEFR4B,I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetB
      class(mqc_nef_set),pointer::ref_NEFSetB
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetB => NEFSetB
      if(NEFSetB%isMirrorCell(1,2)) then
        full_NEFSetB = mqc_nef_set_full_grid(NEFSetB)
        ref_NEFSetB => full_NEFSetB
      end if
      if(.not.NEFSetB%isOpen().and.trim(NEFSetB%Type()).ne.'r4') call mqc_error_A(&
           'mqc_nefr4_nefset_difference: NEF_Set must hold r4 elements (or be open)',&
           6,'Array_Type',NEFSetB%Type())
!
      call NEFSetOut%init(NEFSetB%NumI,NEFSetB%NumJ,&
           sameBraKetSpace=(NEFSetB%isSameBraKetSpace().and.NEFR4A%isSameBraKet()),&
           symmetry=merge(NEFSetB%getSymm(),'asymmetric      ',NEFSetB%getSymm().eq.NEFR4A%getSymm()))
      do I = 1, NEFSetB%NumI
        do J = 1, NEFSetB%NumJ
          if(NEFSetB%isZero(I,J)) then
            call mqc_nef_set_put_built_r4(NEFSetOut,NEFR4A-NEFSetB%atR4(I,J),I,J)
          else
            call mqc_nef_set_put_built_r4(NEFSetOut,NEFR4A-ref_NEFSetB%Set_R4(I,J),I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetA,NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetA
      class(mqc_nef_set),pointer::ref_NEFSetA
      type(mqc_nef_set),target::full_NEFSetB
      class(mqc_nef_set),pointer::ref_NEFSetB
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetA => NEFSetA
      if(NEFSetA%isMirrorCell(1,2)) then
        full_NEFSetA = mqc_nef_set_full_grid(NEFSetA)
        ref_NEFSetA => full_NEFSetA
      end if
      ref_NEFSetB => NEFSetB
      if(NEFSetB%isMirrorCell(1,2)) then
        full_NEFSetB = mqc_nef_set_full_grid(NEFSetB)
        ref_NEFSetB => full_NEFSetB
      end if
      if(trim(NEFSetA%Type()).ne.trim(NEFSetB%Type())) call mqc_error_A(&
           'mqc_nefset_nefset_sum: NEF_Sets must hold the same kind',&
           6,'NEFSetA%Type',NEFSetA%Type(),'NEFSetB%Type',NEFSetB%Type())
      if(NEFSetA%NumI.ne.NEFSetB%NumI.or.NEFSetA%NumJ.ne.NEFSetB%NumJ) call mqc_error_i(&
           'mqc_nefset_nefset_sum: NEFSetA and NEFSetB must have identical '//&
           'grid shape',6,'NEFSetA%NumI',int(NEFSetA%NumI,kind=int64),&
           'NEFSetB%NumI',int(NEFSetB%NumI,kind=int64))
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ,&
           sameBraKetSpace=(NEFSetA%isSameBraKetSpace().and.NEFSetB%isSameBraKetSpace()),&
           symmetry=merge(NEFSetA%getSymm(),'asymmetric      ',NEFSetA%getSymm().eq.NEFSetB%getSymm()))
!
      select case(trim(NEFSetA%Type()))
      case('vector')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call mqc_nef_set_put_built_vector(NEFSetOut,ref_NEFSetB%Set_Vector(I,J),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call mqc_nef_set_put_built_vector(NEFSetOut,ref_NEFSetA%Set_Vector(I,J),I,J)
            else
              call mqc_nef_set_put_built_vector(NEFSetOut,ref_NEFSetA%Set_Vector(I,J)+ref_NEFSetB%Set_Vector(I,J),I,J)
            end if
          end do
        end do
      case('matrix')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call mqc_nef_set_put_built_matrix(NEFSetOut,ref_NEFSetB%Set_Matrix(I,J),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call mqc_nef_set_put_built_matrix(NEFSetOut,ref_NEFSetA%Set_Matrix(I,J),I,J)
            else
              call mqc_nef_set_put_built_matrix(NEFSetOut,ref_NEFSetA%Set_Matrix(I,J)+ref_NEFSetB%Set_Matrix(I,J),I,J)
            end if
          end do
        end do
      case('r4')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call mqc_nef_set_put_built_r4(NEFSetOut,ref_NEFSetB%Set_R4(I,J),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call mqc_nef_set_put_built_r4(NEFSetOut,ref_NEFSetA%Set_R4(I,J),I,J)
            else
              call mqc_nef_set_put_built_r4(NEFSetOut,ref_NEFSetA%Set_R4(I,J)+ref_NEFSetB%Set_R4(I,J),I,J)
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
      type(mqc_nef_set),intent(in),target::NEFSetA,NEFSetB
      type(mqc_nef_set)::NEFSetOut
      integer::I,J
!
      type(mqc_nef_set),target::full_NEFSetA
      class(mqc_nef_set),pointer::ref_NEFSetA
      type(mqc_nef_set),target::full_NEFSetB
      class(mqc_nef_set),pointer::ref_NEFSetB
!     Mirror cells of a grid-symmetric input are read from a full-grid copy;
!     other inputs are read in place.
      ref_NEFSetA => NEFSetA
      if(NEFSetA%isMirrorCell(1,2)) then
        full_NEFSetA = mqc_nef_set_full_grid(NEFSetA)
        ref_NEFSetA => full_NEFSetA
      end if
      ref_NEFSetB => NEFSetB
      if(NEFSetB%isMirrorCell(1,2)) then
        full_NEFSetB = mqc_nef_set_full_grid(NEFSetB)
        ref_NEFSetB => full_NEFSetB
      end if
      if(trim(NEFSetA%Type()).ne.trim(NEFSetB%Type())) call mqc_error_A(&
           'mqc_nefset_nefset_difference: NEF_Sets must hold the same kind',&
           6,'NEFSetA%Type',NEFSetA%Type())
      if(NEFSetA%NumI.ne.NEFSetB%NumI.or.NEFSetA%NumJ.ne.NEFSetB%NumJ) call mqc_error_i(&
           'mqc_nefset_nefset_difference: NEFSetA and NEFSetB must have '//&
           'identical grid shape',6,'NEFSetA%NumI',int(NEFSetA%NumI,kind=int64),&
           'NEFSetB%NumI',int(NEFSetB%NumI,kind=int64))
!
      call NEFSetOut%init(NEFSetA%NumI,NEFSetA%NumJ,&
           sameBraKetSpace=(NEFSetA%isSameBraKetSpace().and.NEFSetB%isSameBraKetSpace()),&
           symmetry=merge(NEFSetA%getSymm(),'asymmetric      ',NEFSetA%getSymm().eq.NEFSetB%getSymm()))
!
      select case(trim(NEFSetA%Type()))
      case('vector')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call mqc_nef_set_put_built_vector(NEFSetOut,ref_NEFSetB%Set_Vector(I,J)*(-1.0),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call mqc_nef_set_put_built_vector(NEFSetOut,ref_NEFSetA%Set_Vector(I,J),I,J)
            else
              call mqc_nef_set_put_built_vector(NEFSetOut,ref_NEFSetA%Set_Vector(I,J)-ref_NEFSetB%Set_Vector(I,J),I,J)
            end if
          end do
        end do
      case('matrix')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call mqc_nef_set_put_built_matrix(NEFSetOut,ref_NEFSetB%Set_Matrix(I,J)*(-1.0),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call mqc_nef_set_put_built_matrix(NEFSetOut,ref_NEFSetA%Set_Matrix(I,J),I,J)
            else
              call mqc_nef_set_put_built_matrix(NEFSetOut,ref_NEFSetA%Set_Matrix(I,J)-ref_NEFSetB%Set_Matrix(I,J),I,J)
            end if
          end do
        end do
      case('r4')
        do I = 1, NEFSetA%NumI
          do J = 1, NEFSetA%NumJ
            if(NEFSetA%isZero(I,J).and.NEFSetB%isZero(I,J)) then
              cycle
            else if(NEFSetA%isZero(I,J)) then
              call mqc_nef_set_put_built_r4(NEFSetOut,ref_NEFSetB%Set_R4(I,J)*(-1.0),I,J)
            else if(NEFSetB%isZero(I,J)) then
              call mqc_nef_set_put_built_r4(NEFSetOut,ref_NEFSetA%Set_R4(I,J),I,J)
            else
              call mqc_nef_set_put_built_r4(NEFSetOut,ref_NEFSetA%Set_R4(I,J)-ref_NEFSetB%Set_R4(I,J),I,J)
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
      type(mqc_scf_integral)::matBlock
      class(mqc_nef_matrix),intent(in)::NEFMat
      Type(MQC_SCF_Eigenvalues),Dimension(:),Allocatable,intent(in)::&
           map_matrix_1,map_matrix_2,map_matrix_3,map_matrix_4
      character(len=*),intent(in),optional::spinType,spinTypeKet
      integer(kind=int64),dimension(4),intent(in),optional::Order
      type(mqc_nef_r4)::NEFR4Out
      type(mqc_scf_eigenvalues),dimension(:),allocatable::axisMap1,axisMap2,axisMap3,axisMap4
      character(len=64)::mySpinType,mySpinTypeKet
      integer(kind=int64)::ia,ja,ka,la
      integer(kind=int64)::dimi,dimj,dimk,diml,i,j,k,l,idij,idlk
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
        axisMap1 = map_matrix_1
      else if(ja.eq.1) then
        axisMap1 = map_matrix_2
      else if(ka.eq.1) then
        axisMap1 = map_matrix_3
      else
        axisMap1 = map_matrix_4
      end if
      if(ia.eq.2) then
        axisMap2 = map_matrix_1
      else if(ja.eq.2) then
        axisMap2 = map_matrix_2
      else if(ka.eq.2) then
        axisMap2 = map_matrix_3
      else
        axisMap2 = map_matrix_4
      end if
      if(ia.eq.3) then
        axisMap3 = map_matrix_1
      else if(ja.eq.3) then
        axisMap3 = map_matrix_2
      else if(ka.eq.3) then
        axisMap3 = map_matrix_3
      else
        axisMap3 = map_matrix_4
      end if
      if(ia.eq.4) then
        axisMap4 = map_matrix_1
      else if(ja.eq.4) then
        axisMap4 = map_matrix_2
      else if(ka.eq.4) then
        axisMap4 = map_matrix_3
      else
        axisMap4 = map_matrix_4
      end if
!
!     The expanded R4 (inverse of collapse) is the dagger-consistent image of R(p,q,r,s) = <I|a+_p a_q O a+_r a_s|J> only when the
!     row pair is the excitation q->p of the bra and the column pair r->s of the ket, i.e.
!     the dagger (axis x -> 5-x) maps them onto each other: ka=5-ia, la=5-ja (e.g. Order
!     [2,1,3,4]; user, 2026-11-09), with the same spin type on bra and ket. Checked for the
!     singlet, triplet and all expands (triplet/all column pair fixed 2026-11-09).
      call NEFR4Out%init(axisMap1,axisMap2,axisMap3,axisMap4,&
           symmetry=merge(NEFMat%Symmetry,'asymmetric      ',ka.eq.5-ia.and.la.eq.5-ja.and.&
           mySpinType.eq.mySpinTypeKet),sameBraKet=NEFMat%sameBraKet)
!
      select case(mySpinType)
      case ('singlet')
        do dimi = 1, NEFR4Out%blockSize(ia)
         do dimj = 1, NEFR4Out%blockSize(ja)
          idij = NEFR4Out%blockSize(ja)*(dimi-1) + dimj
          do dimk = 1, NEFR4Out%blockSize(ka)
           do diml = 1, NEFR4Out%blockSize(la)
            idlk = NEFR4Out%blockSize(la)*(dimk-1) + diml
            if (NEFMat%isBlockZero(idij,idlk)) cycle
            matBlock = NEFMat%iat(idij,idlk)
            do i = 1, NEFR4Out%blockSize(ia,dimi,'alpha')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'alpha')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,1_int64],&
                     matBlock%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'singlet')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'singlet')),&
                       'alpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,-1_int64],&
                     matBlock%at(&
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
                     matBlock%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'singlet')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'singlet')),&
                       'alphabeta'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,-1_int64],&
                     matBlock%at(&
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
          idij = NEFR4Out%blockSize(ja)*(dimi-1) + dimj
          do dimk = 1, NEFR4Out%blockSize(ka)
           do diml = 1, NEFR4Out%blockSize(la)
            idlk = NEFR4Out%blockSize(la)*(dimk-1) + diml
            if (NEFMat%isBlockZero(idij,idlk)) cycle
            matBlock = NEFMat%iat(idij,idlk)
            do i = 1, NEFR4Out%blockSize(ia,dimi,'alpha')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'beta')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,-1_int64],&
                     matBlock%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'triplet')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'triplet')),&
                       'alpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,1_int64],&
                     matBlock%at(&
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
                     matBlock%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'triplet')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'triplet')),&
                       'betaalpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,1_int64],&
                     matBlock%at(&
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
          idij = NEFR4Out%blockSize(ja)*(dimi-1) + dimj
          do dimk = 1, NEFR4Out%blockSize(ka)
           do diml = 1, NEFR4Out%blockSize(la)
            idlk = NEFR4Out%blockSize(la)*(dimk-1) + diml
            if (NEFMat%isBlockZero(idij,idlk)) cycle
            matBlock = NEFMat%iat(idij,idlk)
            do i = 1, NEFR4Out%blockSize(ia,dimi,'alpha')
             do j = 1, NEFR4Out%blockSize(ja,dimj,'alpha')
              do k = 1, NEFR4Out%blockSize(ka,dimk,'alpha')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,1_int64],[la,diml,l,1_int64],&
                     matBlock%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','a',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'all')),&
                       'alpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,-1_int64],&
                     matBlock%at(&
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
                     matBlock%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','b',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','a',dimk,diml,ka,la,'all')),&
                       'betaalpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'beta')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,-1_int64],&
                     matBlock%at(&
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
                     matBlock%at(&
                       int(NEFR4Out%collapseIndex(i,j,'a','b',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'all')),&
                       'alpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,1_int64],[ja,dimj,j,-1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,1_int64],&
                     matBlock%at(&
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
                     matBlock%at(&
                       int(NEFR4Out%collapseIndex(i,j,'b','a',dimi,dimj,ia,ja,'all')),&
                       int(NEFR4Out%collapseIndex(k,l,'a','b',dimk,diml,ka,la,'all')),&
                       'betaalpha'))
               end do
              end do
              do k = 1, NEFR4Out%blockSize(ka,dimk,'beta')
               do l = 1, NEFR4Out%blockSize(la,diml,'alpha')
                call NEFR4Out%put([ia,dimi,i,-1_int64],[ja,dimj,j,1_int64],&
                     [ka,dimk,k,-1_int64],[la,diml,l,1_int64],&
                     matBlock%at(&
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
