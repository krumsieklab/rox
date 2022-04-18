! This code is adapted from concreg package
! Please refer to https://cran.r-project.org/web/packages/concreg/
SUBROUTINE concregfit(cards, parms, IOARRAY)
!DEC$ ATTRIBUTES DLLEXPORT :: concregfit

IMPLICIT DOUBLE PRECISION (A-H,O-Z)
double precision, dimension (17) :: parms
integer N, IP,JCODE,Irobust,ISEP,ITER,IMAXIT,IMAXHS,Ioffset
double precision, dimension (int(parms(1))) :: T, WT, gt,  Doffset
double precision, dimension (int(parms(1)),int(parms(2))) :: X
!double precision, dimension (int(parms(1)),int(parms(2))) :: X, XMSF, bresx
double precision, dimension (int(parms(2))) :: B, B0, FD, zw1, xx, yy, dinit
double precision, dimension (int(parms(2)),int(parms(2))) :: SD, VM,  fish, rvm
!integer, dimension (int(parms(1))) :: ibresc, IC, ICMSF, patid
integer, dimension (int(parms(1))) :: IC, strata, help
integer, dimension (int(parms(2))) :: IFLAG
!double precision, dimension (int(parms(15)), int(parms(2))) :: dfbeta
double precision, dimension (int(parms(1)),(int(parms(2))+6+int(parms(16)))) :: cards
double precision, dimension (int((3+2*(parms(2)))),int((parms(2)))) :: IOARRAY
!double precision, dimension (14) ::  DER,EREST
logical, dimension (int(parms(2)),int(parms(2))) :: mask
!double precision, dimension (int(parms(1)), int(parms(2))) :: score_weights
integer, dimension (int(parms(14))+1) :: stratamap
integer maxstrata, inpar



INTRINSIC DABS, DSQRT

ifail=0
irobust=0
N=int(parms(1))
IP=int(parms(2))
Irobust=int(parms(3))
imaxit=int(Parms(4))
imaxhs=int(parms(5))
step=parms(6)
xconv=parms(7)
ioffset=int(parms(16))
maxstrata=int(parms(14))
parms(10)=-11
inpar=int(parms(17))
isep=0
ilastlike=0

dinit=ioarray(2,:)
iflag=int(ioarray(1,:))
t=cards(:,ip+1+ioffset)
ic=int(cards(:,ip+2+ioffset))    ! 0...censored, 1... event, 2... competing event
wt=cards(:,ip+3+ioffset)    ! Gewicht zu den Zeitpunkten
Gt=cards(:,ip+4+ioffset)    ! Follow-up-KM f?r competing event- gewichte
strata=int(cards(:,ip+6+ioffset))

! define stratamap: first indices of each stratum plus N+1


stratamap(1)=1
if (maxstrata .gt. 1) then
 do istrat=2,maxstrata
  where (strata<istrat)
   help=1
  elsewhere
   help=0
  end where
  stratamap(istrat)=sum(help)+1
 end do
end if
stratamap(maxstrata+1)=N+1



mask=.FALSE.
do j=1,Ip
 mask(j,j)=.TRUE.
end do

x=cards(:,(ioffset+1):(ioffset+ip))
Doffset = 0.
if (ioffset .eq. 1) then
 Doffset = cards(:,1)
end if
XL=0.
xl0=xl-.000002

b0=0.

where(iflag .eq. 0)
 b0=dinit
elsewhere(iflag .eq. 2)
 b0=dinit
 iflag=1
endwhere

isflag=sum(iflag)
!write(6,*) "ISFLAG", isflag
b(:)=b0(:)

ITER=0
iconv=0
JCODE=0

do while((iconv .eq. 0) .and. (iter .lt. imaxit))
! write(6,*) iter, b
 iter=iter+1
 b0(:)=b(:)
 XL0=XL
 parms(10)=-10
! write(6,*) "Vor 1. LIKE", b
 if (iter .eq. 1) then
  CALL LIKE(N,IP,X,T,IC,XL,FD,SD,VM,B,wt, gt,doffset,irobust,rvm, maxstrata, stratamap, inpar)
  ! LIKE(N,IP,X,T1,t2,IC,XL,FD,vm,B,JCODE,ngv,score_weights,ntde,ft,ftmap,ilastlike,doffset,ainv)
 end if
! write(6,*) "Nach 1. LIKE"

 parms(10)=-9
 parms(8)= real (JCODE)
 IF (JCODE .GE. 1) RETURN
 parms(11)=xl
 parms(10)=iter
 parms(9)=isep
 If (ISFLAG.ne.0) then
  IF(IFAIL.ne.0) then
!   "Save" Variance matrix if INVERT failed
 !  WRITE(6,*) 'Inversion failed', ITER,IFIRTH
   ICONV=0
   parms(8)=3
   return
  else
   DO I=1,(IP)
    IF (IFLAG(I).EQ.1) then
     TT=dot_product(vm(I,:),fd(:)*iflag(:))
     IF (DABS(TT).GT.STEP) TT=TT/DABS(TT)*STEP
     B(I)=B(I)+TT
    end if
   end do
!   half step if new log-likelihood is less than old one
   ICONV=0
   IHS=0

   CALL LIKE(N,IP,X,T,IC,XL,FD,SD,VM,B,wt, gt,doffset,irobust,rvm,  maxstrata, stratamap, inpar)
   do while(((XL .le. XL0) .AND. (ITER.ne.1)) .AND. (ihs .le. imaxhs))
    IHS=IHS+1
    where (iflag .eq. 1)
     b=(b+b0)/2
    end where
    CALL LIKE(N,IP,X,T,IC,XL,FD,SD,VM,B,wt, gt,doffset,irobust,rvm, maxstrata, stratamap, inpar)
   end do
  end if
 end if
 ICONV=1
 if (isflag .gt. 0) then
  XX=dabs(B-B0)
  IF(any(XX.GT.xconv)) ICONV=0
 end if
end do



! robuste Varianz berechnen
irobust=1

call like(N,IP,X,T,IC,XL,FD,SD,VM,B,wt, gt,doffset,irobust,RVM,  maxstrata, stratamap, inpar)

! alles in parms und ioarray eintragen

parms(10)= iter
ioarray(3,:) = B
ioarray(4:3+IP,:) = VM
ioarray((4+IP):(3+2*IP),:) = RVM




fish = RVM
call invrt(fish,ip)


!do j=1,(ip)
! ioarray(3,j)=b(j)
! do k=1,(ip)
!  ioarray(3+j,k)=vm(j,k)
! end do
!end do
!return

!stderr=dsqrt(pack(VM,mask))
!parms(10)=-8



yy=pack(vm,mask)
yy=dabs(yy)
if (any(yy .gt. 10000)) isep=1



zw=0.
zw1=matmul(fish, b)
zw=dot_product(b,zw1)
parms(9)=zw

parms(8)=jcode
parms(11)=xl
parms(10)=iter

!cards(1:numbpatients, 1:(ip)) = dfbeta

!close(unit=6)

RETURN

end

SUBROUTINE INVRT(A,IA)


!
!...original matrix=a inverse matrix =a (on exit)
!...note that a is changed on exit
!
 INTEGER IA,n
! double precision eps
 double precision, dimension (IA,ia) :: A, B, WK
 INTRINSIC DABS

 wk=a

 IFAIL=0
 b=a
 N=ia

 CALL vert(b, IA, N, WK)
 a=b

 RETURN
END

SUBROUTINE INVERT(A,IA,N,B,IFAIL)
!DEC$ ATTRIBUTES DLLEXPORT :: invert


!
!...original matrix=a inverse matrix =b
!...note that a is changed on exit
!...eps is a small quantity used to see if matrix singular
!...ifail on exit ifail=0 ok, ifail=1 matrix nearly singular
!
 INTEGER IA,N,ifail
 !double precision eps
 double precision, dimension (IA,N) :: A, B, WK
 INTRINSIC DABS

 wk=a

 IFAIL=0
 b=a

 CALL vert(b, IA, N, WK)


 RETURN
END
SUBROUTINE LIKE(N,IP,X,T,IC,XL,FD,SD,VM,B,wt, gt,offset,irobust,rvm, maxstrata, stratamap, inpar)
!DEC$ ATTRIBUTES DLLEXPORT :: like

 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 double precision, dimension (IP,IP) :: SD, WK,  vm, rvm
 double precision, dimension (N,IP,IP) :: XXEBX
 double precision, dimension (IP) :: FD, B, paircont
 double precision, dimension (N) :: EBX, BX, T,  offset, wt, gt
 integer, dimension (N) :: IC
 !integer, dimension (IP) :: iflag
 double precision, dimension (N,IP) :: X, XEBX, UCONT
 double precision, dimension (IP) :: score_cont, pseudoxi, pseudoxj
 integer maxstrata
 integer, dimension (maxstrata+1) :: stratamap

 intrinsic dexp, dsqrt

 dlowest=0.000000001
! write(6,*) "in LIKE"
 XL=0.


 ipges=ip

 ! bx=matmul(x,b)
 ! ebx=dexp(bx)


 xl=0.
 fd(:)=0.
 sd(:,:)=0.
 ucont(:,:)=0.


! Likelihood (XL) is only correct if all or none of the variables is weighted


! do i=1,n
!  do j=1,ip
!   xges(i,j)=x(i,j)
!  end do
! end do
if (inpar .eq. 0) then
 bx=matmul(x,b)+offset
 ebx=dexp(bx)
 do i=1,N
  xebx(i,:)=x(i,:)*ebx(i)
  do k=1,ip
   xxebx(i,k,:)=xebx(i,:)*x(i,k)    ! ???
  end do
 end do
end if

do istratum = 1, maxstrata
 if (stratamap(istratum+1)-stratamap(istratum) .GE. 2) then
 do i=stratamap(istratum),(stratamap(istratum+1)-2)
  if (ic(i) .EQ. 1) then
   do j=i+1,(stratamap(istratum+1)-1)
    zeitp=t(i)-0.00001
    if ((ic(j) .NE. 1) .OR. (t(j) .NE. t(i))) then
     if ((ic(j) .EQ. 2) .and. (t(j) .lt. t(i))) then
      wcr = Gt(i) / Gt(j)
     else
      wcr = 1.
     end if
     if (inpar .eq. 1) then
      ! define pseudovalues
      do jj=1,ip
       pseudoxi(jj)=0.5
       pseudoxj(jj)=0.5
       if ((x(i,jj) - x(j,jj)) .GT. 0.0001) then
        pseudoxi(jj)=1.
        pseudoxj(jj)=0.
       end if
       if ((x(i,jj) - x(j, jj)) .LT. -0.0001) then
        pseudoxi(jj)=0.
        pseudoxj(jj)=1.
       end if
      end do
      bx(i)=dot_product(pseudoxi,b)+offset(i)
      bx(j)=dot_product(pseudoxj,b)+offset(j)
      ebx(i)=dexp(bx(i))
      ebx(j)=dexp(bx(j))
      xebx(i,:)=pseudoxi*ebx(i)
      xebx(j,:)=pseudoxj*ebx(j)
      do k=1,ip
       xxebx(i,k,:)=xebx(i,:)*pseudoxi(k)
       xxebx(j,k,:)=xebx(j,:)*pseudoxj(k)
      end do
     else
      pseudoxi=x(i,:)
     end if


     summe = ebx(i)+ ebx(j)*wcr
     score_cont = (xebx(i,:) + xebx(j,:)*wcr)/summe
     xl = xl + wt(i) *(bx(i) - dlog(summe))
     paircont = wt(i) *(pseudoxi - score_cont)
     fd = fd + paircont
     if (irobust .eq. 1) then
      ucont(i,:) = ucont(i,:) + paircont
      ucont(j,:) = ucont(j,:) + paircont
     end if
     do k=1,ip
      wk(k,:) = score_cont * score_cont(k)
     end do
     sd = sd - wt(i) * ( (xxebx(i,:,:)+xxebx(j,:,:)*wcr)/summe - wk )

    end if
   end do
  end if
 end do
 end if
end do

 wk=-sd
 EPS=.000000000001D0
 ifail=0
 if (IPGES .EQ. 1) then
  vm=1/WK
 else
  CALL INVERT(WK,ipges,Ipges,vm,IFAIL)
 end if
 if (irobust .eq. 1) then
  wk = matmul(transpose(ucont),ucont)
  wk = matmul(vm, wk)
  rvm = matmul(wk, vm)
 end if
 RETURN
END



SUBROUTINE vert(v,lv,n,w)

INTEGER, INTENT(IN OUT)                  :: lv
INTEGER, INTENT(IN)                      :: n
double precision, INTENT(IN OUT)                     :: v(lv,N)
double precision, INTENT(OUT)                     :: w(N)
double precision :: s,t
INTEGER :: i,j,k,l,m, p

!Anm GH bei den Dimensionen die Indizes ver?ndert, waren v(lv,1) und w(1) vorher

intrinsic dabs

k = 0
IF ( n == 1 ) GO TO 110
l = 0
m = 1
10    IF ( l == n ) GO TO 90
k = l
l = m
m = m + 1
!     ---------------------------------------
!     |*** FIND PIVOT AND START ROW SWAP ***|
!     ---------------------------------------
p = l
IF ( m > n ) GO TO 30
s = DABS(v(l,l))
DO  i = m,n
  t = DABS(v(i,l))
  IF ( t <= s ) CYCLE
  p = i
  s = t
END DO
w(l) = p
30    s = v(p,l)
v(p,l) = v(l,l)
IF ( s == 0. ) GO TO 120
!     -----------------------------
!     |*** COMPUTE MULTIPLIERS ***|
!     -----------------------------
v(l,l) = -1.
s = 1./s
DO  i = 1,n
  v(i,l) = -s*v(i,l)
END DO
j = l
50    j = j + 1
IF ( j > n ) j = 1
IF ( j == l ) GO TO 10
t = v(p,j)
v(p,j) = v(l,j)
v(l,j) = t
IF ( t == 0. ) GO TO 50
!     ------------------------------
!     |*** ELIMINATE BY COLUMNS ***|
!     ------------------------------
IF ( k == 0 ) GO TO 70
DO  i = 1,k
  v(i,j) = v(i,j) + t*v(i,l)
END DO
70    v(l,j) = s*t
IF ( m > n ) GO TO 50
DO  i = m,n
  v(i,j) = v(i,j) + t*v(i,l)
END DO
GO TO 50
!     -----------------------
!     |*** PIVOT COLUMNS ***|
!     -----------------------
90    l = int(w(k))
DO  i = 1,n
  t = v(i,l)
  v(i,l) = v(i,k)
  v(i,k) = t
END DO
k = k - 1
IF ( k > 0 ) GO TO 90
RETURN
110   IF ( v(1,1) == 0. ) GO TO 120
v(1,1) = 1./v(1,1)
RETURN
120   continue
!WRITE(6,*) 'ERROR: MATRIX HAS NO INVERSE'
END SUBROUTINE vert
