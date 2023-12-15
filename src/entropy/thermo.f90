! This file is part of xtb, modified for crest
!
! Copyright (C) 2017-2020 Stefan Grimme
! Copyright (C) 2019-2020 Sebastian Ehlert
! Copyright (C)      2020 Philipp Pracht
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!========================================================================!
! Taken from the xtb code and modified for crest
!========================================================================!
module crest_thermo
   use iso_fortran_env, only : wp => real64
   implicit none

   real(wp),private,parameter :: pi = 3.1415926535897932384626433832795029_wp
   real(wp),private,parameter :: twopi = 2.0_wp * pi
   real(wp),private,parameter :: kB = 3.166808578545117e-06_wp

   real(wp),private,parameter :: autoaa = 0.52917726_wp
   real(wp),private,parameter :: aatoau = 1.0_wp/autoaa
   real(wp),private,parameter :: amutokg = 1.660539040e-27_wp
   real(wp),private,parameter :: autokj = 2625.49964038_wp
   real(wp),private,parameter :: autokcal = 627.50947428_wp
   real(wp),private,parameter :: kcaltoau = 1.0_wp/autokcal
   real(wp),private,parameter :: autorcm = 219474.63067_wp
   real(wp),private,parameter :: rcmtoau = 1.0_wp/autorcm
   real(wp),private,parameter :: metokg  = 9.10938356e-31_wp
   real(wp),private,parameter :: kgtome  = 1.0_wp/metokg

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

subroutine thermodyn(iunit,A_rcm,B_rcm,C_rcm,avmom_si,linear,atom,sym,molmass, &
      &              vibs,nvibs,T,sthr_rcm,et,ht,g,ts,zp,pr)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer, intent(in)  :: iunit       !< output_unit
   integer, intent(in)  :: nvibs       !< number of vibrational frequencies
   real(wp),intent(in)  :: A_rcm       !< rotational constants in cm-1
   real(wp),intent(in)  :: B_rcm       !< rotational constants in cm-1
   real(wp),intent(in)  :: C_rcm       !< rotational constants in cm-1
   real(wp),intent(in)  :: avmom_si    !< average moment of inertia in whatever
   real(wp),intent(in)  :: sym         !< symmetry number
   real(wp),intent(in)  :: molmass     !< molecular mass in in amu
   real(wp),intent(in)  :: T           !< temperature in K
   real(wp),intent(in)  :: sthr_rcm    !< rotor cutoff in cm-1
   real(wp),intent(out) :: et          !< enthalpy in Eh
   real(wp),intent(out) :: ht          !< enthalpy in Eh
   real(wp),intent(out) :: g           !< free energy in Eh
   real(wp),intent(out) :: ts          !< entropy in Eh
   real(wp),intent(in)  :: zp          !< zero point vibrational energy in Eh
   real(wp),intent(in)  :: vibs(nvibs) !< vibrational frequencies in cm-1
   logical, intent(in)  :: linear      !< is linear
   logical, intent(in)  :: atom        !< only one atom
   logical, intent(in)  :: pr          !< clutter the screen with printout
   real(wp),parameter :: R = 1.98726D0    ! GAS CONSTANT IN CALORIES/MOLE
   real(wp),parameter :: H = 6.626176D-27 ! PLANCK'S CONSTANT IN ERG-SECONDS
   real(wp),parameter :: AK = 1.3807D-16  ! BOLTZMANN CONSTANT IN ERG/DEGREE
   real(wp),parameter :: conv3 = amutokg*1000 ! 1.6606d-24
   real(wp),parameter :: magic4 = 2.2868d0 ! unknown
   real(wp),parameter :: magic5 = 2.3135d0 ! unknown
   real(wp),parameter :: caltoj = autokj/autokcal

   integer  :: i
   real(wp) :: s_tr,s_rot,s_vib,s_int,s_tot
   real(wp) :: h_tr,h_rot,h_vib,h_int,h_tot
   real(wp) :: q_tr,q_rot,q_vib,q_int
   real(wp) :: cptr,cprot,cpvib,cpint,cptot
   real(wp) :: beta,sthr,avmom,A,B,C
   real(wp) :: ewj,omega,mu
   real(wp) :: wofrot
   real(wp) :: cp_ho,cp_rr
   real(wp) :: sv_ho,sv_rr
   !*******************************************************************

   ! convert EVERYTHING to atomic units NOW and avoid horror and dispair later
   beta=1.0_wp/kB/T ! beta in 1/Eh
   !c1=h*ac/ak/T    ! beta in cm
   sthr  = sthr_rcm * rcmtoau ! sthr in Eh
   avmom = avmom_si*kgtome*aatoau**2*1.0e+20_wp ! in me·α² (au)
   A = A_rcm * rcmtoau
   B = B_rcm * rcmtoau
   C = C_rcm * rcmtoau

   !   ***   INITIALISE SOME VARIABLES   ***
   q_vib = 1.0_wp
   h_vib = 0.0_wp
   cpvib = 0.0_wp
   s_vib = 0.0_wp
   q_rot = 1.0_wp
   h_rot = 0.0_wp
   cprot = 0.0_wp
   s_rot = 0.0_wp
   ! free rotor heat capacity is constant 0.5*R
   ! we use the below to work with Eh² units 
   ! conversion to cal/mol/K is done later
   cp_rr = 0.5_wp/(beta**2) 

   ! construct the frequency dependent parts of partition function
   do i=1,nvibs
      omega=vibs(i)
      ! omega in Eh, beta in 1/Eh
      ewj=exp(-omega*beta)
      q_vib=q_vib/(1.0_wp-ewj)
      ! h_vib in Eh
      h_vib=h_vib+omega*ewj/(1.0_wp-ewj)
      ! cp_ho in Eh²
      cp_ho=omega**2 * ewj/(1.0_wp-ewj)/(1.0_wp-ewj)
      ! replace low-lying vibs for S by rotor approx.
      mu = 0.5_wp / (omega + 1.0e-14_wp)
      ! this reduced moment limits the rotational moment of inertia for
      ! this vibration to that of the total molecule rotation/3
      ! avmom and xmom are in me·α² (au)
      mu = mu*avmom/(mu+avmom)
      if(omega.gt.0)then
         ! sv is S/R which is dimensionless
         ! harm. osc. entropy
         sv_ho = vibs(i)*beta*ewj/(1.0_wp-ewj) - log(1.0_wp-ewj)
         ! free rotor entropy
         ! Cramer, page 328 for one degree of freedom or
         ! http://cccbdb.nist.gov/thermo.asp, eq. 35, sigma=1
         !sv_rr = (0.5_wp+log(sqrt(8.0_wp*pi**3*xmom*sik*t)/sih))
         sv_rr = 0.5_wp + log(sqrt(pi/beta*2.0_wp*mu))
      else
         sv_ho = 0.0_wp
         sv_rr = 0.0_wp
      endif
      ! fermi weigthing
      ! wofrot=1./(1.+exp( (omega-sthr)/20.0 ) )
      ! Head-Gordon weighting
      wofrot=1.0_wp-chg_switching(omega,sthr)
      ! heat capacity (cp_rr is a constant), all in Eh²
      cpvib = cpvib + ((1.0_wp-wofrot)*cp_ho + wofrot*cp_rr)
      ! entropy s_vib is converted to cal/mol/K... by multiplying with R
      s_vib=s_vib+R*((1.0_wp-wofrot)*sv_ho + wofrot*sv_rr)
   enddo
   !   ***   FINISH CALCULATION OF VIBRATIONAL PARTS   ***
   ! now unit conversion again...
   ! h_vib in Eh, beta is in 1/Eh, T is in K, R is in cal/mol/K,
   h_vib=h_vib*R*beta*T
   ! same here
   ! cpvib is in Eh², beta in 1/Eh, R in cal/mol/K
   cpvib=cpvib*R*beta**2
   !   ***   NOW CALCULATE THE ROTATIONAL PARTS  (FIRST LINEAR MOLECULES)
   if (.not.atom) then
      if(linear) then
         ! A is in Eh, beta is in 1/Eh, q_rot is dimensionless
         q_rot=1/(beta*A*sym)
         ! h_rot is in cal/mol
         h_rot=R*T
         ! cprot  is in cal/mol/K
         cprot=R
         ! s_rot is in cal/mol/K
         s_rot=R*((log(1.0_wp/(beta*A*sym)))+1.0_wp)
      else
         ! see above
         q_rot=sqrt(pi/(A*B*C*beta**3))/sym
         h_rot=3.0_wp*R*T/2.0_wp
         cprot=3.0_wp*R/2.0_wp
         s_rot=0.5_wp*R*(-2.0_wp*log(sym)+log(pi/(A*B*C*beta**3))+3.0_wp)
      endif
   endif
   !   ***   CALCULATE INTERNAL CONTRIBUTIONS   ***
   q_int=q_vib*q_rot
   h_int=h_vib+h_rot
   cpint=cpvib+cprot
   s_int=s_vib+s_rot
   !   ***   CONSTRUCT TRANSLATION CONTRIBUTIONS   ***
   q_tr=(sqrt(2.0_wp*pi*molmass*t*ak*conv3)/h)**3
   ! this is 3/2rt+pv=5/2rt
   h_tr=5.0_wp*R*T/2.0_wp
   cptr=5.0_wp*R/2.0_wp
   s_tr=magic4*(5.0_wp*log10(t)+3.0_wp*log10(molmass))-magic5
   !   ***   CONSTRUCT TOTALS   ***
   cptot=cptr+cpint
   s_tot=s_tr+s_int
   h_tot=h_tr+h_int

   if(pr)then
      write(iunit,'(a)')
      write(iunit,'("   temp. (K)  partition function ", &
         & "  enthalpy   heat capacity  entropy")')
      write(iunit,'(  "                                   ", &
         & "cal/mol     cal/K/mol   cal/K/mol   J/K/mol")')
      write(iunit,'(  f7.2,"  VIB ",G10.3,10X,3F11.3)') &
         & T,q_vib,  h_vib,  cpvib,  s_vib
      write(iunit,'(7X,"  ROT ",G10.3,10X,3F11.3)') &
         & q_rot,  h_rot,  cprot,  s_rot
      write(iunit,'(7X,"  INT ",G10.3,10X,3F11.3)') &
         & q_int,h_int,cpint,s_int
      write(iunit,'(7X,"  TR  ",G10.3,10X,3F11.3)') &
         & q_tr, h_tr, cptr, s_tr
      write(iunit,'(7X,"  TOT ",21X,F11.4,3F11.4)') &
         & h_tot,cptot,s_tot,s_tot*caltoj
   endif

   ht=h_tot/1000.0_wp*kcaltoau
   et=ht+zp
   ts=s_tot*t/1000.0_wp*kcaltoau

   g=et-ts

   return
end subroutine thermodyn

pure elemental function lnqrot(temp,f,avmom) result(lnq_r)
   implicit none
   real(wp), parameter :: rcmtoj = rcmtoau*autokj*1000.0_wp
   real(wp), parameter :: avogadro = 6.0221413e23_wp ! 1/mol
   real(wp), parameter :: planck = 6.62606957e-34_wp ! J*s
   real(wp), parameter :: hbar = planck/(2.0_wp*pi)
   real(wp), parameter :: kb = 1.3806488e-23_wp ! J/K
   real(wp), intent(in) :: temp   !< temperature in K
   real(wp), intent(in) :: f      !< vibrational frequency in cm⁻¹
   real(wp), intent(in) :: avmom  !< average moment of inertia in kg·m²
   real(wp):: e
   real(wp):: mu
   real(wp):: lnq_r
   real(wp):: t_rot

   ! moment of inertia corresponding to the rotor with frequency f(ifreq)
   ! convert frequency first from cm⁻¹ to J, we add also a little offset to avoid Infinities
   e = ((f*rcmtoJ)+1.0e-14_wp)/avogadro
   mu = hbar**2/(2.0_wp*e)
   ! the vibrational moment of inertia is now in SI
   ! reduce the moment relative to the total rotational moment
   mu = avmom*mu/(avmom+mu)
   ! now we need a rotational temperature of mu,
   ! since we are SI already no unit conversion needed
   t_rot = hbar**2/(2.0_wp*mu*kb)
   ! ln(q) of the free rotor partition function, we assume σ=1
   lnq_r = log(sqrt(pi*temp/t_rot))

end function lnqrot

pure elemental function lnqvib(temp,f) result(lnq_v)
   implicit none
   real(wp), parameter  :: planck = 6.62606957e-34_wp ! J*s
   real(wp), parameter  :: kb = 1.3806488e-23_wp ! J/K
   real(wp), parameter  :: speed_of_light = 299792458.0_wp ! m/s
   real(wp), parameter  :: factor = planck*100.0_wp*speed_of_light/kb
   real(wp), intent(in) :: temp   !< temperature in K
   real(wp), intent(in) :: f      !< vibrational frequency in cm⁻¹
   real(wp):: lnq_v
   real(wp):: t_vib
   ! get the vibrational temperature (which is in K, BTW)
   t_vib = factor*f
   ! modified oscillator, for sthr = 0 -> HO.
   ! ln(q) of the harmonic oscillator partition function
   lnq_v = - 0.5_wp*t_vib/temp - log(1.0_wp - exp(-t_vib/temp))

end function lnqvib

pure elemental function lnqvibmod(temp,f,sthr,avmom) result(lnq)
   implicit none
   real(wp), intent(in) :: temp   !< temperature in K
   real(wp), intent(in) :: f      !< vibrational frequency in cm⁻¹
   real(wp), intent(in) :: sthr   !< rotor cutoff in cm⁻¹
   real(wp), intent(in) :: avmom  !< average moment of inertia in kg·m²
   real(wp):: fswitch
   real(wp):: lnq_r,lnq_v,lnq
   ! ln(q) of the harmonic oscillator partition function
   lnq_v = lnqvib(temp,f)
   ! ln(q) of the free rotor partition function, we assume σ=1
   lnq_r =  lnqrot(temp,f,avmom)

   ! Chai--Head-Gordon weighting
   fswitch = 1.0_wp - chg_switching(sthr,f)

   ! now final modified vibrational partiation function
   lnq = (1.0_wp-fswitch) * lnq_v + fswitch * lnq_r

end function lnqvibmod

pure elemental function chg_switching(omega,sthr) result(f)
   real(wp),intent(in) :: omega
   real(wp),intent(in) :: sthr
   real(wp) :: f
   if(sthr.ge.0.0_wp) then
      f = 1.0_wp/(1.0_wp+(sthr/omega)**4)
   else
      f = 1.0_wp
   endif
end function chg_switching

pure elemental function chg_inverted(f,sthr) result(omega)
   real(wp),intent(in) :: f
   real(wp),intent(in) :: sthr
   real(wp) :: omega
   omega = sthr/(1.0_wp/f - 1.0_wp)**0.25_wp

end function chg_inverted


subroutine print_thermo_sthr_ts(iunit,nvib,vibs,avmom_si,sthr_rcm,temp)
   implicit none

   integer, intent(in) :: iunit      !< output unit, usually STDOUT
   integer, intent(in) :: nvib       !< number of frequencies
   real(wp),intent(in) :: vibs(nvib) !< frequencies in Eh
   real(wp),intent(in) :: avmom_si   !< average moment
   real(wp),intent(in) :: sthr_rcm   !< rotor cutoff
   real(wp),intent(in) :: temp       !< temperature

   integer  :: i
   real(wp) :: maxfreq,omega,s_r,s_v,fswitch
   real(wp) :: beta,ewj,mu,RT,sthr,avmom
   beta = 1.0_wp/kB/temp ! beta in 1/Eh
   sthr = sthr_rcm * rcmtoau ! sthr in Eh
   RT = kb*temp*autokcal ! RT in kcal/mol for printout
   avmom = avmom_si*kgtome*aatoau**2*1.0e+20_wp ! in me·α²

   write(iunit,'(a)')
   maxfreq = max(300.0_wp,chg_inverted(0.99_wp,sthr_rcm))
   write(iunit,'(1x,a,f6.2,a,/)') 'Frequencies up to ',maxfreq,' cm⁻¹ treated with '// &
   & 'the modified RRHO approximation:'
   write(iunit,'(a8,a14,1x,a27,a27,a12)') &
      "mode","ω/cm⁻¹","T·S(HO)/kcal·mol⁻¹","T·S(FR)/kcal·mol⁻¹","T·S(vib)"
   write(iunit,'(3x,72("-"))')
   do i = 1, nvib
      ! frequency is Eh
      omega=vibs(i)
      ! omega in Eh, beta in 1/Eh
      ewj=exp(-omega*beta)
      ! moment of intertia corresponding to the rotor with frequency omega
      ! mu is in me·α² (au)
      mu = 0.5_wp / (omega+1.0e-14_wp)
      ! this reduced moment limits the rotational moment of inertia for
      ! this vibration to that of the total molecule rotation/3
      ! avmom and mu are in au
      mu=mu*avmom/(mu+avmom)
      !              free rotor entropy
      ! Cramer, page 328 for one degree of freedom or
      ! http://cccbdb.nist.gov/thermo.asp, eq. 35, sigma=1
      !              harm. osc. entropy
      if(omega.gt.0)then
         ! this is S/R which is dimensionless
         s_v = omega*beta*ewj/(1.0_wp-ewj) - log(1.0_wp-ewj)
         s_r = 0.5_wp + log(sqrt(pi/beta*2.0_wp*mu))
      else
         s_v = 0.0_wp
         s_r = 0.0_wp
      endif
      ! Head-Gordon weighting
      fswitch=1.0_wp-chg_switching(omega,sthr)
      if (omega > maxfreq*rcmtoau) exit
      write(iunit,'(i8,f10.2,2(f12.5,1x,"(",f6.2,"%)"),f12.5)') &
         i,omega*autorcm,-RT*s_v,(1.0_wp-fswitch)*100, &
         -RT*s_r,fswitch*100,-RT*((1.0_wp-fswitch) * s_v + fswitch * s_r)
   enddo
   write(iunit,'(3x,72("-"))')

end subroutine print_thermo_sthr_ts

!========================================================================================!
!========================================================================================!
end module crest_thermo
