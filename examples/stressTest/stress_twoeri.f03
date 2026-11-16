program stress_twoeri
! Stress test of TwoERI_Trans (-> TwoERI_Trans_Full) on a real R_8_1234 'regular' AO ERI
! (or StorFull for the *_full cases). Usage: stress_twoeri <n> <case>
! Optional third argument: a file name; all computed blocks are then written there (values).
! case: spaceR       space, real MOs                         Map [1,1,1,1]
!       spaceC       space, complex MOs                      Map [1,1,1,1]
!       spaceC_LR    space, complex bra and ket MOs          Map [1,2,1,2] (Right_MOs)
!       spinR        spin (alpha, beta), real MOs            Map [1,1,1,1]
!       spinC        spin (alpha, beta), complex MOs         Map [1,1,1,1]
!       genC         general (4 blocks), complex MOs         Map [1,1,1,1]
! The *_full cases (spaceR_full, spaceC_full, spinR_full, spinC_full, genC_full) use the
! same MOs on an AO ERI stored StorFull (no AO symmetry): references for the timing.
use mqc_algebra
use mqc_est
implicit none
type(mqc_r4tensor)::ao,blk
type(mqc_twoeris)::eri,moeri
type(mqc_scf_integral)::mo,mo2
type(mqc_matrix),dimension(4)::m
real(kind=real64),dimension(:,:,:,:),allocatable::x
real(kind=real64),dimension(:,:),allocatable::a,b
complex(kind=real64),dimension(:,:),allocatable::z
character(len=16)::arg,cs
character(len=256)::dumpf
character(len=4),dimension(16)::bl
type(mqc_scalar)::sv
integer(kind=int64)::ib,j2,k2,l2,i2
integer::ud
character(len=4)::stor
character(len=64)::ityp,mtyp,lab
integer(kind=int64)::n,cnt0,cnt1,rate,i,iout,iprint
integer::nseed
integer,dimension(:),allocatable::seed
logical::ok,cplx
call get_command_argument(1,arg); read(arg,*) n
call get_command_argument(2,cs)
iout = 6; iprint = 0
call random_seed(size=nseed); allocate(seed(nseed)); seed = [(17+3*i,i=1,nseed)]; call random_seed(put=seed)
! AO ERI: real, 8-fold symmetric values
allocate(x(n,n,n,n)); call random_number(x)
x = x+reshape(x,shape(x),order=[2,1,3,4])
x = x+reshape(x,shape(x),order=[1,2,4,3])
x = x+reshape(x,shape(x),order=[3,4,1,2])
ao = x; deallocate(x)
stor = 'full'
if(index(trim(cs),'_full').eq.0) then
  call MQC_R4Tensor_Full2Symm(ao,'R_8_1234',ok)
  if(.not.ok) stop 'Full2Symm failed'
  stor = 'symm'
endif
ityp = 'regular'
call mqc_twoeris_allocate(eri,stor,ityp,ao)
call MQC_Deallocate_R4Tensor(ao)
! MO coefficient blocks (random, n x n)
cplx = .not.(trim(cs).eq.'spaceR'.or.trim(cs).eq.'spaceR_full'.or.trim(cs).eq.'spinR' &
  .or.trim(cs).eq.'spinR_full')
allocate(a(n,n),b(n,n))
do i = 1, 4
  call random_number(a)
  if(cplx) then
    call random_number(b); allocate(z(n,n)); z = cmplx(a,b,kind=real64); m(i) = z; deallocate(z)
  else
    m(i) = a
  endif
enddo
lab = 'mo'
select case(trim(cs))
case('spaceR_full','spaceC_full','spaceR','spaceC')
  mtyp = 'space'; call mqc_integral_allocate(mo,lab,mtyp,m(1))
case('spaceC_LR')
  mtyp = 'space'; call mqc_integral_allocate(mo,lab,mtyp,m(1))
  call mqc_integral_allocate(mo2,lab,mtyp,m(2))
case('spinR','spinC','spinR_full','spinC_full')
  mtyp = 'spin'; call mqc_integral_allocate(mo,lab,mtyp,m(1),m(2))
case('genC','genC_full')
  mtyp = 'general'; call mqc_integral_allocate(mo,lab,mtyp,m(1),m(2),m(3),m(4))
case default
  stop 'case must be spaceR, spaceC, spaceC_LR, spinR, spinC, genC or a *_full reference &
    &(spaceR_full, spaceC_full, spinR_full, spinC_full, genC_full)'
end select
call system_clock(cnt0,rate)
if(trim(cs).eq.'spaceC_LR') then
  call TwoERI_Trans(iout,iprint,mo,eri,moeri,Right_MOs=mo2)
else
  call TwoERI_Trans(iout,iprint,mo,eri,moeri)
endif
call system_clock(cnt1)
blk = moeri%getBlock('aaaa')
write(*,'(a12,i5,f12.3,a,1x,a,1x,a,1x,a,a,i8)') trim(cs),n,real(cnt1-cnt0)/real(rate),' s', &
  trim(moeri%type()),trim(moeri%storType()),'aaaa='//trim(blk%getSymm()),'  peakMB=',peak_mb()
if(command_argument_count().ge.3) then
  call get_command_argument(3,dumpf)
  bl = [character(len=4)::'aaaa','bbbb','aabb','bbaa','abab','abba','baab','baba', &
    'aaab','aaba','abaa','baaa','bbba','bbab','babb','abbb']
  open(newunit=ud,file=trim(dumpf),status='replace')
  do ib = 1, 16
    if(trim(moeri%type()).eq.'space'.and.ib.gt.1) exit
    if(trim(moeri%type()).eq.'spin'.and.ib.gt.4) exit
    blk = moeri%getBlock(bl(ib))
    if(.not.MQC_R4Tensor_isAllocated(blk)) cycle
    write(ud,'(a,1x,a)') bl(ib),trim(blk%getSymm())
    do l2=1,n; do k2=1,n; do j2=1,n; do i2=1,n
      sv = blk%at(i2,j2,k2,l2)
      write(ud,'(2es24.15)') real(sv%cval()),aimag(sv%cval())
    enddo; enddo; enddo; enddo
  enddo
  close(ud)
endif
contains
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
end program stress_twoeri
