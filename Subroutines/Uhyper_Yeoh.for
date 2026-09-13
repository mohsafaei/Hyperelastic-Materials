c---------------------------------------------------------------c
c                                                               c
c            UHYPER subroutine: Gent material model             c
c                                                               c
c---------------------------------------------------------------c
      subroutine uhyper(bi1,bi2,aj,u,ui1,ui2,ui3,temp,noel,
     1 cmname,incmpflag,numstatev,statev,numfieldv,fieldv,
     2 fieldvinc,numprops,props)

      implicit none
      ! variables that must be defined within the subroutine
      real*8 u(2),ui1(3),ui2(6),ui3(6),statev(*)
      ! variables that are passed to the subroutine for information
      character*80 cmname
      real*8 bi1,bi2,aj,temp,fieldv(*),fieldvinc(*),props(*)
      integer noel,incmpflag,numstatev,numfieldv,numprops
      
      ! variables that are used locally inside the subroutine
      real*8 mu,kapa,imax
      real*8 zero,half,one,two,three
      parameter(zero=0.d0,half=5.d-1,one=1.d0,two=2.d0,three=3.d0,six=6.d0)
      
      ! get the material properties
      if(numprops.eq.3) then
         c1 = props(1)
         c2 = props(2)
         c3 = props(3)
      else
         print*, '***Error: the number of properties of the Gent
     + model in the uhyper subroutine must be 3.***'
         call xit
      endif 
      
      ! strain energy density function:
      ! deviatoric part of the strain energy density
      u(2) = c1*(bi1-three)+ c2*(bi1 - three)**two + c3*(bi1-three)**three
      ! total strain energy (deviatoric and volumetric)
      !u(1) = u(2) + half*kapa*(aj-one)**two
	  u(1) = u(2)
      
      ! first derivatives of the strain energy with respect to
      !  the stretch invariants
      ! du/di1
      ui1(1) = c1+c2*(two*bi1 - six) + three*c3*(bi1-three)**two
      ! du/di2
      ui1(2) = zero
      ! du/di3
      ui1(3) = zero
      
      ! second derivatives of the strain energy
      ! d2u/d2i1
      ui2(1) = 2*c2 + three*c3*(2*bi1-six)
      ! d2u/d2i2
      ui2(2) = zero
      ! d2u/d2i3
      ui2(3) = zero
      ! d2u/di1di2
      ui2(4) = zero
      ! d2u/di1di3
      ui2(5) = zero
      ! d2u/di2di3
      ui2(6) = zero
      
      ! third derivatives of the strain energy
      ui3 = zero
      
      return
      end
	  
	  
	  # Equilibrium (domain integral)
weq_dom = - test(d(rcur,x))*sig_rr
          + test(rcur) * ( 2*lam_r/rcur ) * ( sig_rr - sig_tt )

# Incompressibility constraint (domain integral)
wJ_dom  = test(pLM) * ( lam_r*lam_t^2 - 1 )