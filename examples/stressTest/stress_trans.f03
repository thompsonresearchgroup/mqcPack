program stress_trans
! Stress test of R4_Trans_Full: time to transform an n^4 tensor of a given class.
! Usage: stress_trans <n> <class>
! One MO matrix C (Conj [T,F,T,F]) unless noted:
!   Full   StorFull, real            FullC  StorFull, complex
!   S2     S_2_13_24, real           S2C    S_2_13_24, complex
!   S2P    S_2_12_34 (2143*), complex
!   S4     S_4_1234, complex         Klein  S_4_1234 (Klein), real
!   KleinC S_4_1234_plain, complex; Conj [F,F,F,F] (keeps the plain 2143)
!   R8     R_8_1234, real            R4     R_4_12_34, real
!   R2     R_2_12_3_4, real
!   R8C    R_8_1234, real tensor with a complex C (-> S_4_1234)
!   R8LR   R_8_1234, real, two real matrices, Map [1,2,1,2] (-> S_2_13_24)
! Real classes use a real C, complex classes a complex C.
use mqc_algebra
implicit none
type(mqc_r4tensor)::ao,out
type(mqc_matrix)::c,c2
real(kind=real64),dimension(:,:,:,:),allocatable::x,y
complex(kind=real64),dimension(:,:,:,:),allocatable::xz
real(kind=real64),dimension(:,:),allocatable::a,b
complex(kind=real64),dimension(:,:),allocatable::z
character(len=16)::arg,cls
character(len=64)::lab
integer(kind=int64)::n,cnt0,cnt1,rate,i
integer::nseed
integer,dimension(:),allocatable::seed
logical,dimension(4)::cj
logical::ok,twoMat
call get_command_argument(1,arg); read(arg,*) n
call get_command_argument(2,cls)
call random_seed(size=nseed); allocate(seed(nseed)); seed = [(17+3*i,i=1,nseed)]; call random_seed(put=seed)
allocate(a(n,n),b(n,n)); call random_number(a); call random_number(b)
ok = .true.
twoMat = .false.
cj = [.true.,.false.,.true.,.false.]
lab = ''
select case(trim(cls))
! real tensors
case('Full','S2','R8','Klein','R4','R2','R8C','R8LR')
  allocate(x(n,n,n,n)); call random_number(x)
  select case(trim(cls))
  case('R8','R8C','R8LR')
    x = x+reshape(x,shape(x),order=[2,1,3,4])
    x = x+reshape(x,shape(x),order=[1,2,4,3])
    x = x+reshape(x,shape(x),order=[3,4,1,2])
    lab = 'R_8_1234'
  case('R4')
    x = x+reshape(x,shape(x),order=[2,1,3,4])
    x = x+reshape(x,shape(x),order=[1,2,4,3])
    lab = 'R_4_12_34'
  case('R2')
    x = x+reshape(x,shape(x),order=[2,1,3,4])
    lab = 'R_2_12_3_4'
  case('Klein')
    x = x+reshape(x,shape(x),order=[2,1,4,3])
    x = x+reshape(x,shape(x),order=[3,4,1,2])
    lab = 'Klein'
  case('S2')
    x = x+reshape(x,shape(x),order=[3,4,1,2])
    lab = 'S_2_13_24'
  case('Full')
  case default
    stop 'stress_trans: unknown real class'
  end select
  ao = x; deallocate(x)
  if(len_trim(lab).gt.0) call MQC_R4Tensor_Full2Symm(ao,trim(lab),ok)
  select case(trim(cls))
  case('R8C')
    allocate(z(n,n)); z = cmplx(a,b,kind=real64); c = z
  case('R8LR')
    c = a; c2 = b; twoMat = .true.
  case default
    c = a
  end select
! complex tensors
case('FullC','S4','S2C','S2P','KleinC')
  allocate(x(n,n,n,n),y(n,n,n,n)); call random_number(x); call random_number(y)
  allocate(xz(n,n,n,n)); xz = cmplx(x,y,kind=real64); deallocate(x,y)
  select case(trim(cls))
  case('S4')
    xz = xz+reshape(xz,shape(xz),order=[3,4,1,2])
    xz = xz+conjg(reshape(xz,shape(xz),order=[2,1,4,3]))
    lab = 'S_4_1234'
  case('S2C')
    xz = xz+reshape(xz,shape(xz),order=[3,4,1,2])
    lab = 'S_2_13_24'
  case('S2P')
    xz = xz+conjg(reshape(xz,shape(xz),order=[2,1,4,3]))
    lab = 'S_2_12_34'
  case('KleinC')
    xz = xz+reshape(xz,shape(xz),order=[2,1,4,3])
    xz = xz+reshape(xz,shape(xz),order=[3,4,1,2])
    lab = 'S_4_1234_plain'
    cj = .false.
  case('FullC')
  case default
    stop 'stress_trans: unknown complex class'
  end select
  ao = xz; deallocate(xz)
  if(len_trim(lab).gt.0) call MQC_R4Tensor_Full2Symm(ao,trim(lab),ok)
  allocate(z(n,n)); z = cmplx(a,b,kind=real64); c = z
case default
  stop 'class must be Full, FullC, S2, S2C, S2P, S4, Klein, KleinC, R8, R4, R2, R8C or R8LR'
end select
if(.not.ok) stop 'Full2Symm failed'
! peak memory of the transform only (input build excluded)
call reset_peak()
call system_clock(cnt0,rate)
if(twoMat) then
  out = R4_Trans_Full(ao,c,c2,Map=[1,2,1,2])
else
  out = R4_Trans_Full(ao,c,Conj=cj)
endif
call system_clock(cnt1)
write(*,'(a6,i5,f12.3,a,1x,a,1x,a,a,i8)') trim(cls),n,real(cnt1-cnt0)/real(rate),' s', &
  trim(ao%getSymm()),'->'//trim(out%getSymm()),'  peakMB=',peak_mb()
contains
  subroutine reset_peak()
! reset VmHWM to the current resident size (Linux: write 5 to /proc/self/clear_refs),
! so peak_mb() counts only what follows; no effect where the file cannot be written
  integer::u,ios
  open(newunit=u,file='/proc/self/clear_refs',status='old',action='write',iostat=ios)
  if(ios.eq.0) then
    write(u,'(a)',iostat=ios) '5'
    close(u)
  endif
  end subroutine reset_peak
  integer(kind=int64) function peak_mb()
  character(len=256)::line
  integer::u,ios
  peak_mb = -1
  open(newunit=u,file='/proc/self/status',status='old',action='read',iostat=ios)
  if(ios.ne.0) return
  do
    read(u,'(a)',iostat=ios) line
    if(ios.ne.0) exit
    if(line(1:6).eq.'VmHWM:') then
      read(line(7:),*) peak_mb
      peak_mb = peak_mb/1024
      exit
    endif
  enddo
  close(u)
  end function peak_mb
end program stress_trans
