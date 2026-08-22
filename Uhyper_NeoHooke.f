c---------------------------------------------------------------c
c                                                               c
c         UHYPER subroutine: Neo-Hookean material model         c
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
      real*8 mu
      real*8 zero,half,one,two,three
      parameter(zero=0.d0,half=5.d-1,one=1.d0,two=2.d0,three=3.d0)
      ! get the material properties
      if(numprops.eq.1) then
         mu   = props(1)
      else
         print*, '***Error: the number of properties of the Mooney-Rivlin
     + model in the uhyper subroutine must be 1.***'
         call xit
      endif 
	  
	  
	  
      
      
	  ! strain energy density function:
      ! deviatoric part of the strain energy density
	  u(2) = mu* (bi1 - three) / two
      
	  ! total strain energy (deviatoric and volumetric)
      u(1) = u(2)

	  ! first derivatives of the strain energy with respect to
      !  the stretch invariants
      ! du/di1
      ui1(1) = mu/2
      ! du/di2
      ui1(2) = zero
      ! du/di3
      ui1(3) = zero
      
      ! second derivatives of the strain energy
      ! d2u/d2i1
      ui2(1) = zero
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
      
	  
	  
	  statev(1) = bi1
      statev(2) = bi2
	  
	  return
      end