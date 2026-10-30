! Standalone algorithm harness, independent of the Python port and legacy I/O.
program reference
  implicit none
  integer*2 :: ix(5762),mhv(360,2),idm(3,2),lx(1440,2),idum,yjd
  real*4 :: xd,dx,sv(2),snd(2,2),av(2,2),se(2,2),mode(2,2),median(2,2)
  integer*2 :: mhi(2),mlo(2),rng(2),mends(2),ida(10),mqd(10),th(10),vh(10)
  character*2 :: id
  character :: debug,snh
  integer*2 :: ikm,ixy,igs,isp,icsp,igk,mp(8,2),icc(8,2),i3r(8,2),idt,nbh,neh,ntr(2)
  real*4 :: del(8,2),sesrc(8,2),sep(8,2),sepa(8,2)
  real*4 :: xx(124),yy(124),yhat(124,2),z(124,2),b1(124,2),delta(124,2)
  integer*2 :: ns(2),e,nppm,nm,ncp,npts,jda,i,j,quiet(3),idx(9),ki
  real*4 :: k9,fmin9,value,rkv
  common /able/ ix
  common /baker/ mhv,idm,xd
  common /charly/ lx,idum,yjd
  common /dog/ mhi,mlo,id,dx,sv,snd,rng,mends,ida,mqd,av,se,mode,median,th,vh
  common /fox/ debug,ikm,ixy,igs,isp,icsp,igk
  common /golf/ del,sesrc,mp,icc,i3r
  common /hotel/ idt,sep,sepa,snh,nbh,neh,ntr
  common /india/ xx,yy,yhat,z,b1,delta,ns
  integer*2 :: rnd
  read(*,*) nppm,ncp,xd,k9
  read(*,*) quiet
  do e=1,2
    read(*,*) mhi(e),mlo(e),mends(e),rng(e)
    read(*,*) sep(:,e)
    read(*,*) sepa(:,e)
  end do
  nm=1440/nppm
  mhv=32767; idm=0; ix=0; lx=0;xx=0;yy=0;yhat=0;z=0;b1=0;delta=0
  debug='N';snh='N';ida=-99;del=32767;sesrc=0
  open(7,status='scratch')
  do i=1,3
    if(quiet(i)==1) ida(i)=i
  end do
  do jda=1,3
    do e=1,2
      do i=1,1440
        read(*,*) ix((e-1)*1440+i)
      end do
    end do
    call kmean(jda,nppm,nm)
    if(jda<3) call movdat(nm)
  end do
  do e=1,2
    call xyset(e,nm,nppm,npts)
    call gs(e,nm,nppm,npts,ncp)
    call spln2(e,nm)
    call gkr(e)
  end do
  fmin9=k9/2500.
  call curvk(fmin9,idx)
  do e=1,2
    do i=1,8
      value=rnd(del(i,e))
      if(del(i,e)>=32767 .or. value>6000) then
        write(*,'(ES24.16,1X,F8.1)') del(i,e),99.9
      else
        value=rkv(value,fmin9,idx,ki)
        write(*,'(ES24.16,1X,F8.1)') del(i,e),value
      end if
    end do
  end do
end program
