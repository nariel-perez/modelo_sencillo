!!!*****************************************************************************
!!!*****************************************************************************

program microgels_MC

!!!*****************************************************************************
!!!*****************************************************************************

!!! NVT Monte Carlo simulation of a microgel solution
!!! microgel's size distribution can evolve according to a simple model
  
  
  use parameters
  use simple_microgel
  use nvt_mc  
  
  implicit none

  logical:: montecarlo
  
  





  call read_inputs(montecarlo)


  call bulk


  call isolated_microgel  



  

  if (montecarlo) then
     call monte_carlo
 
  else
     call print_stuff
  endif
   








contains




!!!.............................................................................
  subroutine print_stuff
!!!*****************************************************************************
!!!*****************************************************************************
    implicit none

    integer:: iamin
    integer,parameter:: outfi=85
    
    real(8):: a0,q,qnet,kappa,zeta
    character(12):: chx12,chx12b
    
    

    open(file="results.dat",unit=outfi)

    iamin=minloc(F_mg,1)
    a0=R_mg(iamin)

     
    q=Q_mg(iamin)

    !kappa=kappa_mg(iamin)
   
    call q_net(iamin,qnet)

    call zeta_pot(iamin,a0,qnet,zeta)

        
    write(chx12,"(f8.1)")a0
    write(chx12b,"(f8.1)")2d0*a0
    write(*,'(4xa,1xa,1xa)')"MG radius:",trim(adjustl(chx12)),"nm"
    write(*,'(4xa,1xa,1xa/)')"diameter:",trim(adjustl(chx12b)),"nm"

    
    write(outfi,*)"Radius=",a0
    write(outfi,*)"zeta_pot=",zeta

    
    write(chx12,"(f10.3)")zeta
    write(*,'(4xa,1xa,1xa/)')"zeta pot:",trim(adjustl(chx12)),"mV"


    if (nads>0) then

       write(chx12,'(1pe12.3)')gamma_ads(iamin)
       write(*,'(4xa,1xa,1xa/)')"adsorption:",trim(adjustl(chx12)),"molec."

       write(outfi,*)"gamma_ads=",gamma_ads(iamin)

    endif



    close(outfi)
    



  end subroutine print_stuff











end program microgels_MC
