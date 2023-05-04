!!!..............................................................................

module parameters

!!!*****************************************************************************
!!!*****************************************************************************

  implicit none

  real(8),parameter:: Na=0.602214,pi=acos(-1.0d0),pKw=14d0
  real(8),parameter:: cbjrrm=16710.075,eps_sol=78.5

  real(8),allocatable:: R_mg(:),F_mg(:),Q_mg(:),kappa_mg(:),gamma_ads(:), degree(:)


  type cg_unit
     integer:: npKs
     character(6):: nm,vdw
     real(8),allocatable:: fi(:),f(:),pK(:),z(:)
     ! ideal & system's deg of charge, pK
     real(8):: v,Mw
     real(8):: chi,chi_i,chi_f,T_pt,deltaT
  end type cg_unit
 

  type molecule
     integer:: ntyp ! # unit types
     type(cg_unit),allocatable:: cgu(:)
     real(8):: c,rhoq,rhoqbulk,rho,rhobulk
     real(8):: nchains,ns1,nseg ! # chains, # segs in 1 chain, # segs
     integer,allocatable:: nu(:)
     real(8),allocatable:: xbulk(:),x(:),cn(:) ! vol fraction & composition # in %
     character(10):: molname ! molecule name
  end type molecule

  type(molecule):: poly,ads


  integer:: nads,Np,max_mc,eq_mc

  real(8):: vsol,vH,vOH,vpls,vmin
  real(8):: zH,zOH,zpls,zmin
  real(8):: pH,csalt,temp,lb_sol
  real(8):: xsolbulk,xHbulk,xOHbulk,xplsbulk,xminbulk,kappa_sol
  
  real(8):: R0,Vpol,Mw_pol,DRmg,Rmg_max,Rmg_min,a_min,a_max

  real(8):: conc,dens,units





  
contains
  
  
  
  
  
    
!!!..............................................................................
  
  subroutine read_inputs(montecarlo)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none

    integer:: ityp,i,npKs
    
    real(8):: rpls,rmin,z0,xtot
    
    character(4):: chx4
    character(6):: chx6
    character(8):: chx8,chx8v(10)
    character(20):: title,stread
    character(30):: inputnm

    logical,intent(out):: montecarlo
    logical:: loop
    
    
    



    call getarg(1,inputnm)    

    if (trim(adjustl(inputnm))=="") then
       inputnm="inputs.in"       
    endif
    
    

    open(unit=1,file=trim(adjustl(inputnm)))





    write(*,'(/39("*")/a/39("*"))')"****** Microgel Solution: NVT MC ******"
    write(*,'(/4xa,1xa/)')"Reading from file:",trim(adjustl(inputnm))
    write(*,'(a/)')"--- Input variables --- "
    
    

    
!!!!!!!!!! Environment !!!!!!!!!!
    

    rewind(1)
    
    title="Environment"
    
    loop=.true.
    do while(loop)
       read(1,*,end=37)stread
       if (stread==title) loop=.false.
    enddo
    
37  if (loop) then
       write(*,38)
38     format(/1x"Error in input file - No Environment found!"/)
       stop
    endif
       
    read(1,*)
    


    
    read(1,*)pH
    read(1,*)csalt

    read(1,*)temp
    
    temp=temp+273



    write(chx8,"(f5.2)")pH
    write(*,'(4x,a,1xa)')"pH:",trim(adjustl(chx8))

    write(chx8,"(1pe8.2e2)")csalt
    write(*,'(4x,a,1xa,1xa)')"[salt]:",chx8,"M"

    write(chx8,"(f5.1)")temp
    write(chx4,"(f4.1)")temp-273
    write(*,'(4x,a,1xa,1xa,1xa,1xa)')"Temp:",trim(adjustl(chx8)),"K /",&
         trim(adjustl(chx4)),"C"

    write(*,*)
    




!!!!!!!!!! Solvent !!!!!!!!!!


    read(1,*)
    read(1,*)vsol
    read(1,*)vH
    read(1,*)vOH
    read(1,*)
    read(1,*)rpls
    read(1,*)rmin

    vH=vH/vsol
    vOH=vOH/vsol

    vpls=4d0/3d0*pi*rpls**3/vsol
    vmin=4d0/3d0*pi*rmin**3/vsol

    read(1,*)
    read(1,*)zH
    read(1,*)zOH
    read(1,*)zpls
    read(1,*)zmin


    write(*,'(a)')"--- Solvent ---"
    write(*,'(4x,a5,1xa2,1xf6.3,1xa)')"H2O: ","v=",vsol,"nm^3"
    write(*,'(4x,a5,1xa2,1xf6.3,1xa,1xf4.1)')"H3O+:","v=",vsol*vH,"nm^3 - z=",zH
    write(*,'(4x,a5,1xa2,1xf6.3,1xa,1xf4.1)')"OH-: ","v=",vsol*vOH,"nm^3 - z=",zOH
    write(*,'(4x,a5,1xa2,1xf6.3,1xa,1xf4.1)')"+:   ","v=",vsol*vpls,"nm^3 - z=",zpls
    write(*,'(4x,a5,1xa2,1xf6.3,1xa,1xf4.1)')"-:   ","v=",vsol*vmin,"nm^3 - z=",zmin
    write(*,*)
  




    lb_sol=cbjrrm/temp/eps_sol     ! Bjerrum length in nm

    




!!!!!!!!!! microgel !!!!!!!!!!



    rewind(1)
    
    title="Microgel"
    
    loop=.true.
    do while(loop)
       read(1,*,end=39)stread
       if (stread==title) loop=.false.
    enddo
    
39  if (loop) then
       write(*,40)
40     format(/1x"Error in input file - No Microgel found!"/)
       stop
    endif
       

    read(1,*)


    read(1,*)Rmg_min,Rmg_max,DRmg

    if (Rmg_min>Rmg_max) DRmg=-DRmg


    


!!!!!!!!!! network !!!!!!!!!!
 
    read(1,*)
    read(1,*)poly%nchains

    read(1,*)poly%ns1

    poly%nseg=poly%nchains*poly%ns1

    read(1,*)poly%ntyp

    allocate(poly%cgu(poly%ntyp),poly%x(poly%ntyp),poly%cn(poly%ntyp))


    write(*,'(a)')"--- Microgel ---"
    write(*,'(4xa)')"Polymer:"


    write(chx8,"(1pe8.3e1)")poly%nseg
    write(*,'(4x,a,a)')chx8," segments"

    write(chx8,"(1pe8.3e1)")poly%nchains
    write(*,'(4x,a,a)')chx8," chains"

    write(chx8,"(f8.1)")poly%ns1
    write(*,'(4x,a,a)')trim(adjustl(chx8))," segments/chain"


    Vpol=0d0
    Mw_pol=0d0

    do ityp=1,poly%ntyp
       
       read(1,*)
       read(1,*)poly%cgu(ityp)%nm,poly%cn(ityp)
       read(1,*)poly%cgu(ityp)%Mw,poly%cgu(ityp)%v,z0,poly%cgu(ityp)%npKs
       
       npKs=poly%cgu(ityp)%npKs
              
       allocate(poly%cgu(ityp)%z(0:npKs))

       poly%cgu(ityp)%z(0)=z0

       if (npKs>0) then
          
          allocate(poly%cgu(ityp)%fi(npKs),poly%cgu(ityp)%f(npKs),&
               poly%cgu(ityp)%pK(npKs))

          read(1,*)(poly%cgu(ityp)%pK(i),poly%cgu(ityp)%z(i),i=1,npKs)

          poly%cgu(ityp)%fi=0d0
          poly%cgu(ityp)%f=0d0
       endif
       

       read(1,*)poly%cgu(ityp)%vdw

       poly%cgu(ityp)%vdw=to_lowercase(trim(adjustl(poly%cgu(ityp)%vdw)))


       if (poly%cgu(ityp)%vdw=="chi") then
          
          backspace(1)
          read(1,*)chx6,poly%cgu(ityp)%chi

       elseif (poly%cgu(ityp)%vdw=="logi") then

          backspace(1)
          read(1,*)chx6,poly%cgu(ityp)%chi_i,poly%cgu(ityp)%chi_f,&
               poly%cgu(ityp)%T_pt,poly%cgu(ityp)%deltaT

          poly%cgu(ityp)%T_pt=poly%cgu(ityp)%T_pt+273

       endif

                    

       poly%x(ityp)=0d0



       chx8=trim(poly%cgu(ityp)%nm)//":"
       write(chx8v(1),'(f6.1)')100*poly%cn(ityp)
       write(chx8v(2),'(f8.1)')poly%cgu(ityp)%Mw
       write(chx8v(3),'(f8.1)')poly%cgu(ityp)%z(0)

       write(*,'(4x,a7,a,a,1xf6.3,a,1xa,1xa,a,1xa)')chx8,&
            trim(adjustl(chx8v(1))),&
            "% - v=",poly%cgu(ityp)%v," nm^3 - Mw=",trim(adjustl(chx8v(2))),&
            "Da"," - z=",trim(adjustl(chx8v(3)))

       if (npKs>0) then
          do i=1,npKs
             write(*,'(11x,a,1xf4.2,a,1xf4.1)')"pK=",poly%cgu(ityp)%pK(i),&
                  " - z=",poly%cgu(ityp)%z(i)
          enddo
       endif

       Vpol=Vpol+poly%cn(ityp)*poly%cgu(ityp)%v

       Mw_pol=Mw_pol+poly%cn(ityp)*poly%cgu(ityp)%Mw


    enddo

    
    xtot=0d0
    do ityp=1,poly%ntyp
       xtot=xtot+poly%cn(ityp)
    enddo

    if (abs(xtot-1d0)>1d-5) then
       write(*,'(/4xa/)')"ERROR!: network monomers cn's don't add up to 100%..."
       stop
    endif






    Vpol=Vpol*dble(poly%nseg) ! in nm^3

    Mw_pol=Mw_pol*dble(poly%nseg) ! in Da


    write(*,*)



    

    poly%cgu%v=poly%cgu%v/vsol




  
    R0=(3d0*Vpol/4d0/pi)**(1d0/3d0)

    units=4d0/3d0*pi*R0**3


    if (DRmg<0) then
       if (Rmg_max<R0) Rmg_max=R0-DRmg
    else
       if (Rmg_min<R0) Rmg_min=R0+DRmg
    endif




    write(*,'(4xa)')"Elastic model:"
    write(chx8,"(f8.1)")R0
    write(*,'(4xa,1xa,1xa)')"R0:",trim(adjustl(chx8)),"nm"
    write(chx8,"(f8.1)")R0*2d0
    write(*,'(4xa,1xa,1xa)')"(D0:",trim(adjustl(chx8)),"nm)"
    

    write(chx8v(1),"(f8.1)")Rmg_min
    write(chx8v(2),"(f8.1)")Rmg_max
    write(chx8v(3),"(f8.1)")DRmg

    write(*,'(4xa,1x,3(1xa))')"R (nm)=",(trim(adjustl(chx8v(i)(1:6))),i=1,3)



    write(*,'(/4xa/)')"internal energy units (iu): 4/3*pi*R0^3 KBT"












!!!!!!!!!! Adsorbates !!!!!!!!!!



    rewind(1)
    
    title="Adsorbates"
    
    loop=.true.
    do while(loop)
       read(1,*,end=41)stread
       if (stread==title) loop=.false.
    enddo
    
41  if (loop) then
       write(*,42)
42     format(/1x"Error in input file - No Adsorbates found!"/)
       stop
    endif
       


 
    read(1,*)
    read(1,*)nads



    if (nads>0d0) then


       write(*,'(a)')"--- Adsorbates ---"
       
       write(chx8,"(i2)")nads
       write(*,'(4x,a,a)')trim(adjustl(chx8))," adsorbate/s"

       read(1,*)
       read(1,*)ads%molname
       read(1,*)ads%c

       read(1,*)ads%ntyp

       allocate(ads%cgu(ads%ntyp),ads%nu(ads%ntyp))
       allocate(ads%xbulk(0:ads%ntyp),ads%x(0:ads%ntyp))

       ads%x=0d0
       ads%xbulk=0d0
       
       write(*,"(/4x'***'1x,a)")trim(adjustl(ads%molname))
       write(*,'(4x,a,1x1pe9.3e2,1xa)')"concentration:",ads%c,"M"


       do ityp=1,ads%ntyp
          
          read(1,*)ads%cgu(ityp)%nm,ads%nu(ityp),ads%cgu(ityp)%v,z0,&
               ads%cgu(ityp)%npKs
       
          npKs=ads%cgu(ityp)%npKs
              
          allocate(ads%cgu(ityp)%z(0:npKs))

          ads%cgu(ityp)%z(0)=z0

          if (npKs>0) then
             
             allocate(ads%cgu(ityp)%fi(npKs),ads%cgu(ityp)%f(npKs),&
                  ads%cgu(ityp)%pK(npKs))
             
             read(1,*)(ads%cgu(ityp)%pK(i),ads%cgu(ityp)%z(i),i=1,npKs)
             
             ads%cgu(ityp)%fi=0d0
             ads%cgu(ityp)%f=0d0
          endif


          chx8=trim(ads%cgu(ityp)%nm)//":"
          write(chx8v(1),'(i6)')ads%nu(ityp)
          write(chx8v(2),'(f8.1)')ads%cgu(ityp)%z(0)

          write(*,'(4x,a7,a,a,1xf6.3,a,a,1xa)')chx8,trim(adjustl(chx8v(1))),&
               " units - v=",ads%cgu(ityp)%v," nm^3",&
               " - z=",trim(adjustl(chx8v(2)))

          if (npKs>0) then
             do i=1,npKs
                write(*,'(11x,a,1xf4.2,a,1xf4.1)')"pK=",ads%cgu(ityp)%pK(i),&
                     " - z=",ads%cgu(ityp)%z(i)
             enddo
          endif


       enddo

       ads%cgu%v=ads%cgu%v/vsol

       write(*,*)


    else

 
       write(*,'(a/)')"--- NO adsorbates ---"
       
 
    endif


    




!!!!!!!!!! Monte Carlo !!!!!!!!!!


    montecarlo=.true.


    rewind(1)
    
    title="Monte_Carlo"
    
    loop=.true.
    do while(loop)
       read(1,*,end=43)stread
       if (stread==title) loop=.false.
    enddo
    
43  if (loop) then
       write(*,44)
44     format(/1x"Error in input file - No Monte_Carlo found!"/)
       stop
    endif
       
 
    read(1,*)
    read(1,*)Np
    read(1,*)dens

    read(1,*)max_mc
    read(1,*)eq_mc

    read(1,*)
    read(1,*)montecarlo

    

    if (montecarlo) then


       write(*,'(a)')"--- NVT MC ---"
    
       write(chx8,"(i4)")Np
       write(*,'(4x,a,a)')trim(adjustl(chx8))," microgels"


       write(chx8,"(f8.2)")dens
       write(*,'(4x,a,1xa,1xa)')"density:",trim(adjustl(chx8)),"mg/ml"


       write(chx8v(1),"(i8)")max_mc
       write(chx8v(2),"(i8)")eq_mc


       write(*,'(4xa,1xa)')"Max. MC steps:",trim(adjustl(chx8v(1)))
       write(*,'(4xa,1xa,1xa)')"Equilibration:",trim(adjustl(chx8v(2))),"steps"

      
    else
    
       write(*,'(a)')"--- no Monte Carlo simulation ---"
       
    endif


    write(*,*)









  end subroutine read_inputs










  





!!!..............................................................................
  
  subroutine bulk
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none
    
    integer:: ityp,ipK
    real(8):: Qx
    
    
    
    xHbulk=10**(-pH)*Na*vH  ! bulk H+ volume fraction
    
    xOHbulk=10**(pH-pKw)*Na*vOH ! bulk OH- volume fraction
    
    xplsbulk=csalt*Na*vpls
    xminbulk=csalt*Na*vmin
    

    if (nads>0) then
       
       ads%rhobulk=ads%c*Na ! this has units: nm^3
       
       ads%rhoqbulk=0d0

       ads%xbulk(0)=0d0

       do ityp=1,ads%ntyp          
          
          ads%xbulk(ityp)=ads%rhobulk*ads%cgu(ityp)%v*ads%nu(ityp)
         
          ads%xbulk(0)=ads%xbulk(0)+ads%xbulk(ityp)

          
          do ipK=1,ads%cgu(ityp)%npKs
             ads%cgu(ityp)%fi(ipK)=1d0/(1d0+10d0**(sign(1d0,ads%cgu(ityp)%z(ipK))*&
                  (pH-ads%cgu(ityp)%pK(ipK))))
             
             ads%rhoqbulk=ads%rhoqbulk+ads%xbulk(ityp)*&
                  ads%cgu(ityp)%z(ipK)/ads%cgu(ityp)%v*ads%cgu(ityp)%fi(ipK)
          enddo
          
          ads%rhoqbulk=ads%rhoqbulk+&
               ads%xbulk(ityp)*ads%cgu(ityp)%z(0)/ads%cgu(ityp)%v
          
       enddo
       
       
       Qx=xHbulk*zH/vH+xOHbulk*zOH/vOH+&
            xplsbulk*zpls/vpls+xminbulk*zmin/vmin+ads%rhoqbulk
       
       
    else
       
       
       Qx=xHbulk*zH/vH+xOHbulk*zOH/vOH+&
            xplsbulk*zpls/vpls+xminbulk*zmin/vmin
       
       
    endif
    





    if (Qx>0d0) then     ! pH < 7 - Acidic
       xminbulk=xminbulk-Qx/zmin*vmin ! add HCl
       
    elseif (Qx<0d0) then ! pH > 7 - Basic
       xplsbulk=xplsbulk-Qx/zpls*vpls ! add NaOH
    endif


    xHbulk=xHbulk*vsol
    xOHbulk=xOHbulk*vsol
    xplsbulk=xplsbulk*vsol
    xminbulk=xminbulk*vsol



    if (nads>0) then
       
       ads%xbulk=ads%xbulk*vsol
       
       xsolbulk=1d0-xHbulk-xOHbulk-xplsbulk-xminbulk-ads%xbulk(0)
       
    else
       
       xsolbulk=1d0-xHbulk-xOHbulk-xplsbulk-xminbulk

    endif



       
    do ityp=1,poly%ntyp          
       do ipK=1,poly%cgu(ityp)%npKs
          poly%cgu(ityp)%fi(ipK)=1d0/(1d0+10d0**(sign(1d0,poly%cgu(ityp)%z(ipK))*&
               (pH-poly%cgu(ityp)%pK(ipK))))
       enddo
    enddo
       
       
       



    if (nads>0) then
       kappa_sol=xHbulk/vH+xOHbulk/vOH+xplsbulk/vpls+xminbulk/vmin+abs(ads%rhoqbulk) 
    else
       kappa_sol=xHbulk/vH+xOHbulk/vOH+xplsbulk/vpls+xminbulk/vmin
    endif

    kappa_sol=4d0*pi*lb_sol*kappa_sol/vsol
    kappa_sol=sqrt(kappa_sol)
    
    






  end subroutine bulk




  
  
  function to_lowercase(str_in) result(str_out)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
    ! Original author: Clive Page
    
    implicit none
    
    character(len=*), intent(in):: str_in
    character(len=len(str_in)):: str_out
    integer:: i,j
    
    do i=1,len(str_in)
       j=iachar(str_in(i:i))
       if ( j>=iachar("A") .and. j<=iachar("Z") ) then
          str_out(i:i)=achar(iachar(str_in(i:i))+32)
       else
          str_out(i:i)=str_in(i:i)
       endif
    enddo
    
  end function to_lowercase
  
    
    






end module parameters
