!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Stefan Grimme
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

      subroutine solve(nspins,sigma,jab,nkcp,Ispink,
     .                 dumspin,mfreq,maxline,nline,line)
      implicit none

      real*8 :: jab(nspins*(nspins+1)/2),sigma(nspins),mfreq
      real*8 :: line(2,maxline),Ispink(nspins)
      integer :: nkcp(nspins),dumspin(nspins)
      integer :: nspins,maxline,nline

      integer :: anz_sub, anzp
      integer :: k, isub          
      integer :: i, j, nmt, itmpline
      integer :: anz_mat
      integer :: Gfun,gew
      integer,allocatable :: pascal(:,:)                         
      integer,allocatable :: anz_mt(:),anz_mt_cp(:),anz_k(:)
      real*8 :: norm,dum
      real*8,allocatable :: maxspincp(:),minspincp(:)
      real*8,allocatable :: maxspink(:),minspink(:)
      real*8,allocatable :: isp(:)
      real*8,allocatable :: tmpline(:,:) 
      real*8, allocatable  :: cpspin(:,:)      
      real*8, allocatable  :: spin(:,:)       
      real*4, allocatable  :: A(:,:),U(:,:),E(:)

      allocate(pascal(99,2))
      allocate(anz_mt(nspins),anz_mt_cp(nspins),anz_k(nspins))
      allocate(maxspincp(nspins),minspincp(nspins))
      allocate(maxspink(nspins),minspink(nspins))
      allocate(isp(nspins),tmpline(2,maxline))

      if(nspins.eq.1)then
         nline = nline +1
         line(1,nline) = sigma(1)*mfreq
         line(2,nline) = 1.0d0 * nkcp(1)
         return
      endif

c     scale for speedup
      jab = jab *0.5

      itmpline=0 ! line counter

cProcedure Calc_min_max_cp(anzcp : byte);
      do k=1,nspins
c copy
         anz_k(k)=nkcp  (k)
         isp  (k)=Ispink(k)
         maxspincp(k) = anz_k(k)*isp(k)  
         minspincp(k)=0
         if(mod(anz_k(k),2).ne.0) minspincp(k)=0.5
      enddo
cProcedure Calc_anz_mt(anzk : byte);
      do k=1,nspins
         anz_mt(k)=idint(2.*maxspincp(k))+1
      enddo
cProcedure Calc_anz_mt_cp(anzk : byte);
      do k=1,nspins
         anz_mt_cp(k)=idint(maxspincp(k)-minspincp(k)+1)
      enddo
ccalc_anz_pr(anzahl_kerne, anz_mt_cp);
      anz_sub=1
      do k=1,nspins
         anz_sub=anz_sub*anz_mt_cp(k)
      enddo

      allocate(cpspin(nspins+1,anz_sub))
      call produktfunktionen(nspins, anz_sub, cpspin, 
     .                       maxspincp,
     .                       minspincp, anz_mt_cp)


c     write(*,*) 'number of super blocks ', anz_sub
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do isub=1, anz_sub
cProcedure Calc_Kerne;
         do k=1, nspins
            maxspink(k) =  cpspin(k, isub);
            minspink(k) = -cpspin(k, isub);
         enddo
         do k=1,nspins
            anz_mt(k)=idint(2*maxspink(k)+1)
         enddo
         anzp=1
         do k=1,nspins
            anzp=anzp*anz_mt(k)
         enddo

         allocate(spin(nspins+1,anzp))
         call produktfunktionen(nspins, anzp, spin, maxspink,
     .                          minspink, anz_mt)
c        do j=1,anzp    
c           write(*,'(16F4.1)')(spin(i,j),i=1,nspins+1)
c        enddo
         call pascal_dreieck(nspins,anzp,spin,pascal,anz_mat)
         i=0
         do nmt =1, anz_mat
            if(pascal(nmt,2)-pascal(nmt,1)+1.gt.i)i=
     .         pascal(nmt,2)-pascal(nmt,1)+1
         enddo
         write(*,'(8x,I6,'' product functions'',I6,
     .                   '' Mt blocks, largest is '',I6)')anzp,anz_mat,i

c allocate 
         allocate(A(anzp,anzp),U(anzp,anzp),E(anzp))
c built Mt sub blocks
         A = 0
!$OMP PARALLEL PRIVATE ( nmt ) SHARED ( A )
!$OMP DO 
         do nmt =1, anz_mat
            write(*,'(i3,''('',i4,'')'',$)') 
     .      nmt,pascal(nmt,2)-pascal(nmt,1)+1
            call hmat(pascal(nmt,1),pascal(nmt,2),
     .                nspins, sigma, jab, spin, cpspin, 
     .                anzp, anz_sub, A, mfreq, isub)
         enddo
!$OMP END DO
!$OMP END PARALLEL
         write(*,*)
c diag
         U = 0
         do nmt =1, anz_mat
c           write(*,'(i3,$)') nmt
            call diag(pascal(nmt,1),pascal(nmt,2),anzp,A,U,E)
         enddo
c done
c intensity
         gew=1        
         do i=1,nspins
            gew=gew*Gfun(anz_k(i),cpspin(i,isub),maxspincp(i))
         enddo
         call intens(isub,gew,nspins,spin,cpspin,anz_sub,
     .               pascal,anz_mat,anzp,E,A,U,
     .               tmpline,itmpline,maxline,
     .               dumspin)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         deallocate(A,U,E,spin)
      enddo

c normalize
      norm=0
      do i=1,itmpline
         norm=norm+tmpline(2,i)
      enddo
      dum=0
      do k=1,nspins
         if(dumspin(k).eq.0) dum=dum+nkcp(k)
      enddo
      norm=dum/norm

      do i=1,itmpline
         nline=nline+1
         line(1,nline)=tmpline(1,i)
         line(2,nline)=tmpline(2,i)*norm
      enddo
      
      deallocate(tmpline,isp,minspink,maxspink)
      deallocate(minspincp,maxspincp,anz_k,anz_mt_cp)
      deallocate(anz_mt,pascal)

      end

C------------------------------------------------------------------------------     
C sparse matrix multiply version, SG July 17
C------------------------------------------------------------------------------     
      subroutine intens(isub,gew,nspins,spin,cpspin,anz_sub,
     .                  pascal,anz_mat,anzp,E,A,U,line,nline,
     .                  maxline,dummy)  
      implicit none
      integer isub,gew,anz_mat,nspins,anzp,nline,maxline,targnuc
      integer anz_sub,pascal(99,2)
      integer dummy(nspins)      
      real*8  line(2,maxline)
      real*8  spin  (nspins+1,anzp)                         
      real*8  cpspin(nspins+1,anz_sub)                      
      real*4  A (anzp,anzp)                            
      real*4  U (anzp,anzp)                            
      real*4  E (anzp)                                 

      integer ende2,anfang2,mm,j,i,k,l,ende,anfang,ls
      integer add,add2,n,kern                               
      real*8  moment, xsum, ifak

      character*6 s
      real*4, allocatable  :: val(:)
      integer,allocatable  :: rind(:),cind(:)
      integer nnz
      real*4 dum
      
      open(unit=33,file='tmpanmr_tmat',form='unformatted')

      nnz=0
      do mm = 1,anz_mat-1 
         anfang = pascal(mm,  1) 
         ende =   pascal(mm,  2)
         anfang2 =pascal(mm+1,1)
         ende2  = pascal(mm+1,2)
         do l = anfang, ende 
            do ls = anfang2, ende2 
CCCCC
            add = 0 
            n = 1 
            add2=0
  10        If(spin(n, L) .ne. spin(n, Ls)) add2=add2+1         
            If(spin(n, L) -  spin(n, Ls).eq.1) then 
               add=add+1      
               kern=n        
            endif           
            n=n+1    
            if(add.gt.1.or.add2.gt.1.or.n.gt.nspins) goto 20
            goto 10
  20        continue
c           If (add.eq.1.and.add2.eq.1) Then 
            If (add.eq.1.and.add2.eq.1.and.dummy(kern).eq.0) Then
               dum     =sqrt(cpspin(kern,isub) + spin(kern,l )) *
     .                  sqrt(cpspin(kern,isub) - spin(kern,ls))
               nnz=nnz+1
               write(33) dum,ls,l
            endif

CCCCC
            enddo
         enddo
      enddo

c runtime critical step
c performs  moment_(i->j)=[( U⁺ T U )_ij]²  
c where U are the Eigenvectors and T is the
c dipole transition matrix on file
c U is blocked according to the Pascal triangle, T has entries in between the blocks
c U is less sparse than T
c for larges cases the sparse version is faster by a factor of > 2

      allocate(val(nnz),cind(nnz),rind(nnz))
      rewind 33
      do i=1,nnz
         read(33) val(i),cind(i),rind(i)
      enddo
      close(33,status='delete')

      dum=100.-100.*dble(nnz)/dble(anzp)**2
      write(*,'(''first maxtrix multiply,  sparsity in % '',
     .f9.3,'' ...'')')dum
      s='GxxFxx'
c     val=T is sparse, U is dense
c     A=(TU)
      call mkl_scoomm('N', anzp, anzp, anzp, 1.0, s,           
     .     val, rind, cind, nnz, U, anzp, 0.0, A, anzp)
      deallocate(val,cind,rind)

c     write U to file and allocate
      open(unit=33,file='tmpanmr_tmat',form='unformatted')
      nnz=0
      do i=1,anzp
      do j=1,anzp
         if(abs(u(j,i)).gt.1.d-5) then
            nnz=nnz+1
            write(33) u(j,i),j,i
         endif
      enddo
      enddo

      allocate(val(nnz),cind(nnz),rind(nnz))
      rewind 33
      do i=1,nnz
         read(33) val(i),cind(i),rind(i)
      enddo
      close(33,status='delete')

      dum=100.-100.*dble(nnz)/dble(anzp)**2
      write(*,'(''second maxtrix multiply, sparsity in % '',
     .f9.3,'' ...'')')dum
C     res=U⁺(TU)
c        =U⁺A   
      call mkl_scoomm('N', anzp, anzp, anzp, 1.0, s,           
     .     val, rind, cind, nnz, A, anzp, 0.0, U, anzp)
      deallocate(val,cind,rind)

c     xsum=0
      do j=1,anzp
         do i=1,anzp
            moment=gew*U(i,j)**2
            if(moment.gt.1.d-5)then
               nline=nline+1
               if(nline.gt.maxline) stop 'too many transitions'
               line(1,nline)=E(i)-E(j)  
               line(2,nline)=moment
c              xsum=xsum+line(2,nline)
            endif
         enddo
      enddo
c     write(*,*) 'Intsum = ',xsum,nspins
 
      end


C------------------------------------------------------------------------------     
      subroutine diag(ia,ie,ndim,A,U,E)                              
      implicit none
      integer ia,ie,ndim
      integer n,i,j,k,lwork,info
      real*4 A(ndim,ndim)
      real*4 U(ndim,ndim)
      real*4 E(ndim)
      real*4, allocatable:: eig(:),c(:,:),work(:)

      n=ie-ia+1
c     write(*,*) 'diagonalizing H sub-block of dimension ',n
      if(n.eq.1)then
         E(ia)=A(ia,ie)
         U(ia,ie)=1.0d0
      else
         lwork  = 1 + 6*n + 2*n**2
         allocate(eig(n),c(n,n),work(lwork))
         do i=ia,ie
            do j=ia,i
               c(i-ia+1,j-ia+1) = a(i,j)
               c(j-ia+1,i-ia+1) = a(i,j)
c              H(k)=A(i,j)
            enddo
         enddo
         call ssyev ('V','U',n,c,n,eig,work,lwork,info)
         U(ia:ie,ia:ie)=c(1:n,1:n)
         E(ia:ie)=eig(1:n)
         deallocate(eig,c,work)
      endif

      end

C------------------------------------------------------------------------------     
      subroutine hmat(ia,ie,nspins,sigma,jab,spin,cpspin,
     .                anzp,anz_sub,A,mfreq,isub)
      implicit none
      integer nspins, anzp, ia, ie, isub, anz_sub
      real*8  spin  (nspins+1,anzp)                         
      real*8  cpspin(nspins+1,anz_sub)                      
      real*8 jab(nspins*(nspins+1)/2),sigma(nspins),mfreq
      real*4  A(anzp,anzp)                  

      integer i, j, k, l, iadr, lin 
      integer is, n, add               
      real*8 hilf, hilf2, a1,a2,b,b2

c diagonal
      do k=ia,ie
      hilf = 0
      do l = 1, nspins 
         hilf = hilf+sigma(l)*spin(l, k)*mfreq
      enddo 
      hilf2 = 0
      do i = 1, nspins-1 
         do j = i+1, nspins 
            iadr = lin(i,j)
            hilf2 = hilf2 + jab(iadr) * spin(i, k)*spin(j, k)
         enddo     
      enddo
      A(k,k) = hilf+2.0*hilf2
      enddo

c off-diagonal
      do k=ia,ie
         do L=k+1,ie
CCCCCCCCCCCCCCCCCCCCCCC
          do i=1, nspins                  
             do is=i+1, nspins              
                if(((spin(i,k)-spin(i,L))*
     .              (spin(is,k)-spin(is,L))).eq.-1)then  
                  n=1
                  add=0
ccc
  10              if(spin(n, k).ne.spin(n,L)) add=add+1         
                  n=n+1
                  if(add.le.2.and.n.le.nspins) goto 10
ccc
                  if(add.eq.2) then
                     a1=cpspin(i,isub)*(cpspin(i,isub)+1)
                     b =cpspin(is,isub)*(cpspin(is,isub)+1)
                     a2=spin(i, k) * spin(i, L)
                     b2=spin(is,k) * spin(is, L)
                     iadr = lin(i,is)
                     A(L, k) = jab(iadr) * sqrt(a1-a2) * sqrt(b-b2);
c                    A(k, L) = A(L, k);
                  endif
                endif
             enddo
          enddo
CCCCCCCCCCCCCCCCCCCCCCC
         enddo
      enddo

      end

C-----------------------------------------------------------------------------     
      subroutine produktfunktionen(nspins, anzp, spin, maxspin,
     .                             minspin, anzmt)
      implicit none
      integer nspins, anzp
      real*8  spin(nspins+1,anzp)                                
      real*8  maxspin(nspins),minspin(nspins)
      integer anzmt(nspins)

      integer l, k, i, j ,anzk2, anzk, ii, m, pr, schritte
      real*8 hilf, mi, sc1, pp

      anzk=nspins

      anzk2=anzk+1
      spin = 0
      do i = 1, anzk 
          pr = 1
          do  l = 1, i 
              pr = pr*anzmt(l)
          enddo
          schritte = anzp/pr
          mi = maxspin(i)
          do k = 1, anzp 
             spin(i, k) = mi
             spin(anzk2, k) = mi+spin(anzk2, k) 
             If (mod(k,schritte) .eq. 0) Then 
                mi = mi-1
                If (mi .lt.  minspin(i)) mi = maxspin(i)  
             Endif
          enddo
      Enddo
c sort
      do 140   ii = 2, anzp
         i = ii - 1
         k = i
         pp= spin(anzk2,i)
         do 120   j = ii, anzp
            if (spin(anzk2,j) .lt. pp) go to 120
            k = j
            pp= spin(anzk2,j)
  120    continue
         if (k .eq. i) go to 140
         spin(anzk2,k) = spin(anzk2,i)
         spin(anzk2,i) = pp
         do m=1,anzk
            sc1=spin(m,i)
            spin(m,i)=spin(m,k)
            spin(m,k)=sc1
         enddo
  140 continue

      end

C------------------------------------------------------------------------------     
      subroutine pascal_dreieck(anzk, anzp, spin, 
     .                          pascal, anz_mat)
      implicit none
      integer j, i, k, anzk, anzp, anz_mat
      integer pascal(99,2)
      real*8 spin(anzk+1,anzp)
    
      k = 1      
      j = 1
      pascal(1, 1) = 1
      pascal(1, 2) = 1
      do i = 1, anzp-1 
         If(abs(spin(anzk+1, i)-spin(anzk+1, i+1)).gt.1.d-8) Then 
            pascal(k, 2) = j
            k=k+1  
            j=j+1  
            pascal(k, 1) = j
          else
            j=j+1
          Endif
      enddo
      if(k.gt.99) stop 'internal error'
      pascal(k, 2) = j
      anz_mat = k

      end

C------------------------------------------------------------------------------     
      integer Function fak(x) 
      implicit none
      integer  x,p,z                              
      real*8 integer
      z  = 1 
      If (x .gt. 0) Then 
         do p = 1, x 
            z = z*p
         enddo
      endif
      fak = z
      End

      integer Function pdrei(zeile, element) 
      implicit none
      integer zeile,element,fak
      pdrei=fak(zeile)/(fak(element)*fak(zeile-element))
      End

      integer Function Gfun(ak, iakt, imax)
      implicit none
      integer n1, n2, a, ak, pdrei    
      real*8 iakt, imax
      if(abs(imax-0.5).lt.1.d-8) then 
         a=1 
      Else
         n1 = idint(imax-iakt)
         If(n1.eq.0) Then 
            n2 = 0
         Else 
            n2 = pdrei(ak, n1-1)
         endif
         a = pdrei(ak, n1)-n2
      Endif
      gfun=a
      End



