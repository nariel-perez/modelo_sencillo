!!!..............................................................................

module simple_microgel

!!!*****************************************************************************
!!!*****************************************************************************


  use parameters

  implicit none

  integer,parameter,private:: neq=2
  integer,private:: itmax,nsig,info
  real(8),private:: fnorm,tol,par(1),x(neq),f(neq)


  real(8),private:: psi,xpi,xsol,xH,xOH,xpls,xmin,xpol,kappa






  
contains
  
  
  
  

  
  
  
!!!*****************************************************************************
  subroutine isolated_microgel
!!!*****************************************************************************
    implicit none

    integer:: iR,nRs,n_spl,ityp,ipK
    integer,parameter:: dn_spl=5 ! spline points between data points

    real(8),allocatable:: tmp(:,:),auxR(:)

    real(8):: Rmg
    real(8):: Df,Fel,Fmt,Ftot

    character(8):: chx8

    logical:: first

    integer:: iostat_adsor

    Real(8), allocatable::degree(:,:)


    

    open(file="degree.dat",unit=7)
    
    
    open(file="DF-Rmg.dat",unit=2)
    write(2,*)"### 1. Rmg ### 2. bFtot/(4/3*pi*R0^3) ### 3." ,&
         "bFmt/(4/3*pi*R0^3) ### 4. bFel/(4/3*pi*R0^3)"



    open(file="isolated-mg.dat",unit=9)
   
    write(9,'(/43("*")/aa/43("*")//6xa/)')"***** Isolated Microgel: ",&
         "Simple Model *****","...calculations:"
       
    
    

    first=.true.
    
    Rmg=Rmg_min


    nRs=int(abs(Rmg_min-Rmg_max)/abs(DRmg))+1

    if (nads>0) open(file="adsor.dat",unit=3, iostat = iostat_adsor)
    if (nads>0) then
       allocate(tmp(nRs,5))
    else
       allocate(tmp(nRs,4))
    endif
    
    allocate(degree(nRs,1))
    iR=0



    do


       if (DRmg>0d0.and.Rmg>Rmg_max) exit
       if (DRmg<0d0.and.Rmg<Rmg_max) exit


       call microgel_size(xpol,Rmg)

       
       write(chx8,"(f8.1)")Rmg  

       write(9,'(43("*")//,4xa,1xa,1xa,1x1pe7.2e1/)')"R=",&
            trim(adjustl(chx8(1:6))),"nm - xp=",xpol   
      


       if (first) then          
          x(1)=xsolbulk-xpol
          x(2)=0d0

          first=.false.
       else
          x(1)=xpi
          x(2)=psi
       endif


       itmax=200
       nsig=6
       par(1)=1d0

       tol=1d-5


       
       call hybrd1(fcn,neq,x,f,tol,info)
     

       fnorm=sqrt(sum(f*f))
       write(9,'(4x,"fnorm=",1PE9.2E2)')fnorm

  
       xpi=x(1)
       psi=x(2)
       call calcular(xpi,psi)

       write(9,*)

       do ityp=1,poly%ntyp
          do ipK=1,poly%cgu(ityp)%npKs
             write(9,'(4xa,i2,":"1x,"f_q=",1PE10.3E2,1x"- f_id=",1PE10.3E2)')&
                  trim(adjustl(poly%cgu(ityp)%nm)),ipK,&
                  poly%cgu(ityp)%f(ipK),poly%cgu(ityp)%fi(ipK)
             write(7,*) Rmg, poly%cgu(ityp)%f(ipK)
          enddo
       enddo




       Df=0d0
       Fmt=0d0
       Fel=0d0
       

       call free_energy_dft(xpi,Df)

       call add_elastic(Rmg,Fel)
       
       call sum_FE(Rmg,Df,Fmt,Fel,Ftot)



       if (info==1) then

          iR=iR+1


          do ityp=1,poly%ntyp
             do ipK=1,poly%cgu(ityp)%npKs
                degree(iR,1) = poly%cgu(ityp)%f(ipK)
                write(*,*)iR, degree(iR,1) 
             enddo
          enddo
          

          write(9,'(/4x,"F_tot/(4*pi*R0^3)=",1PE11.4E2,1x"kT")')Ftot
       
          
          write(2,*)Rmg,Ftot,Fmt,Fel
          


          tmp(iR,1)=Rmg
          tmp(iR,2)=Ftot
          tmp(iR,3)=poly%rhoq*4d0*pi/3d0*Rmg**3/vsol
          tmp(iR,4)=kappa
          
          if (nads>0) tmp(iR,5)=(ads%rho-ads%rhobulk)*4d0*pi/3d0*Rmg**3

          if (nads>0) write(3,*)Rmg, tmp(iR,5)

  
       endif

       write(9,*)

       Rmg=Rmg+DRmg


    enddo

    nRs=iR






!!!....... splines .......!!!


    a_min=minval(tmp(:,1))

    a_max=maxval(tmp(:,1))
    
    n_spl=int((a_max-a_min)*dble(dn_spl))+1



    allocate(auxR(n_spl),R_mg(n_spl),F_mg(n_spl),Q_mg(n_spl),kappa_mg(n_spl))
    
    

    call akima(nRs,tmp(:,1),tmp(:,2),n_spl,R_mg,F_mg) ! Free energy spline
    
    call akima(nRs,tmp(:,1),tmp(:,3),n_spl,auxR,Q_mg) ! polymer charge spline

    call akima(nRs,tmp(:,1),tmp(:,4),n_spl,auxR,kappa_mg) ! kappa spline

    

    if (nads>0) then
   
       allocate(gamma_ads(n_spl))

       call akima(nRs,tmp(:,1),tmp(:,5),n_spl,auxR,gamma_ads) ! kappa spline

    endif



    deallocate(tmp,auxR)
    
    



    if (iostat_adsor == 0) close(3)


    close(2)
    close(9)
    !close(3)
    close(7)




    

    










  end subroutine isolated_microgel
  
  
  
  












!!!..............................................................................
  
  subroutine microgel_size(xpol,Rmg)
    
!!!*****************************************************************************
!!!******************************************************************************

    implicit none
    real(8),intent(in):: Rmg
    real(8),intent(out):: xpol
        
    integer:: ityp

    real(8):: x_tmp


    x_tmp=(poly%nseg*vsol/Rmg**3)*3d0/4d0/pi ! auxiliary

    xpol=0d0

    do ityp=1,poly%ntyp       
       poly%x(ityp)=poly%cn(ityp)*poly%cgu(ityp)%v*x_tmp ! vf of each segment
       xpol=xpol+poly%x(ityp)

    enddo
    
        


  end subroutine microgel_size





  













!!!..............................................................................
  
  subroutine calcular(xpi,psi)
    
!!!*****************************************************************************
!!!******************************************************************************

    implicit none
    real(8),intent(in):: xpi,psi
 
    integer:: ityp,ipK
    real(8):: ypi,epsi,fq,nu,chi
    
    
    ypi=xpi/xsolbulk

    epsi=dexp(-psi)
    
    xH=xHbulk*ypi**vH*epsi**zH
    xOH=xOHbulk*ypi**vOH*epsi**zOH
    xmin=xminbulk*ypi**vmin*epsi**zmin
    xpls=xplsbulk*ypi**vpls*epsi**zpls

    xsol=xpi

    do ityp=1,poly%ntyp

       if (poly%cgu(ityp)%vdw=="nipam") then

          call nipam(chi,poly%x(ityp))


       elseif (poly%cgu(ityp)%vdw=="logi") then

          call logistic_chi(poly%cgu(ityp)%chi_i,poly%cgu(ityp)%chi_f,&
               poly%cgu(ityp)%T_pt,poly%cgu(ityp)%deltaT,chi)


       elseif (poly%cgu(ityp)%vdw=="chi") then

          chi=poly%cgu(ityp)%chi*(dble(298)/Temp)
          
       else
          chi=0d0

       endif
       
       xsol=xsol*exp(-chi*poly%x(ityp))

    enddo





    if (nads>0) then

       ads%x(0)=0d0

       ads%rho=ads%rhobulk
       
       do ityp=1,ads%ntyp          
          
          nu=ads%nu(ityp)
          
          ads%rho=ads%rho*ypi**(nu*ads%cgu(ityp)%v)*epsi**(nu*ads%cgu(ityp)%z(0))
          
          do ipK=1,ads%cgu(ityp)%npKs

             fq=ads%cgu(ityp)%fi(ipK)/(1d0-ads%cgu(ityp)%fi(ipK))*&
                  epsi**ads%cgu(ityp)%z(ipK)
             
             ads%cgu(ityp)%f(ipK)=fq/(1d0+fq)
          
             ads%rho=ads%rho*(ads%cgu(ityp)%fi(ipK)/ads%cgu(ityp)%f(ipK))**nu
             
             ads%rho=ads%rho*epsi**(nu*ads%cgu(ityp)%z(ipK))

          enddo
          
          ads%x(ityp)=nu*ads%cgu(ityp)%v*vsol

       enddo
       
       ads%x=ads%rho*ads%x
       
       ads%x(0)=sum(ads%x(1:ads%ntyp))
       
       ads%rhoq=0d0

       do ityp=1,ads%ntyp          

          ads%rhoq=ads%rhoq+ads%x(ityp)*ads%cgu(ityp)%z(0)/ads%cgu(ityp)%v

          do ipK=1,ads%cgu(ityp)%npKs    

             ads%rhoq=ads%rhoq+ads%x(ityp)*ads%cgu(ityp)%f(ipK)*&
                  ads%cgu(ityp)%z(ipK)/ads%cgu(ityp)%v

          enddo
       enddo
       

    endif


    




    poly%rhoq=0d0

    xpol=0d0

    do ityp=1,poly%ntyp          
       
       poly%rhoq=poly%rhoq+poly%x(ityp)*poly%cgu(ityp)%z(0)/poly%cgu(ityp)%v
       
       do ipK=1,poly%cgu(ityp)%npKs

          fq=poly%cgu(ityp)%fi(ipK)/(1d0-poly%cgu(ityp)%fi(ipK))*&
               epsi**poly%cgu(ityp)%z(ipK)
          
          poly%cgu(ityp)%f(ipK)=fq/(1d0+fq)
          
          poly%rhoq=poly%rhoq+poly%x(ityp)*poly%cgu(ityp)%f(ipK)*&
               poly%cgu(ityp)%z(ipK)/poly%cgu(ityp)%v
          
       enddo

       xpol=xpol+poly%x(ityp)
       
    enddo



    


    
    kappa=0d0

    if (nads>0) then
       kappa=xH/vH+xOH/vOH+xpls/vpls+xmin/vmin+abs(ads%rhoq) 
    else
       kappa=xH/vH+xOH/vOH+xpls/vpls+xmin/vmin
    endif

    !kappa=kappa+abs(poly%rhoq) !!!!!!!!!!!!!

    kappa=4d0*pi*lb_sol*kappa/vsol
    kappa=sqrt(kappa)
    

    






  end subroutine calcular







!!!..............................................................................
  
  subroutine fcn(neq,x,f,iflag)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none
    
    integer,intent(in):: neq,iflag
    real(8),intent(out):: f(neq)
    real(8),intent(inout):: x(neq)

    real(8):: xpi

  
    xpi=x(1)
    psi=x(2)


    call calcular(xpi,psi)

 

    if (nads>0) then
 
       f(1)=xH+xOH+xmin+xpls+xsol+xpol+ads%x(0)-1d0

       f(2)=poly%rhoq+xH*zH/vH+xOH*zOH/vOH+xpls*zpls/vpls+xmin*zmin/vmin+ads%rhoq 
    
    else

       f(1)=xH+xOH+xmin+xpls+xsol+xpol-1d0

       f(2)=poly%rhoq+xH*zH/vH+xOH*zOH/vOH+xpls*zpls/vpls+xmin*zmin/vmin

    endif
 



    
  end subroutine fcn















!!!..............................................................................
  
  subroutine free_energy_dft(xpi,Df)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none

    real(8),intent(in):: xpi
    real(8),intent(out):: Df

    integer:: ityp,ipK
    real(8):: bfbulk,bf




         
    bfbulk=-xHbulk/vH-xOHbulk/vOH-xplsbulk/vpls-xminbulk/vmin-xsolbulk
    bfbulk=bfbulk+dlog(xsolbulk) ! in 1/vsol units

    if (nads>0) then
       bfbulk=bfbulk-ads%rhobulk*vsol ! ads%rhobulk has units!
    endif




    bf=-xH/vH-xOH/vOH-xpls/vpls-xmin/vmin-xsol
    bf=bf+(1d0-xpol)*dlog(xpi) 


    do ityp=1,poly%ntyp          
       do ipK=1,poly%cgu(ityp)%npKs

          bf=bf+poly%x(ityp)/poly%cgu(ityp)%v*dlog(1d0-poly%cgu(ityp)%f(ipK))
          ! in 1/vsol units
       
       enddo
    enddo


    if (nads>0) then
       bf=bf-ads%rho*vsol ! ads%rho has units!
    endif





    Df=bf-bfbulk
    Df=Df/vsol

    
    

  end subroutine free_energy_dft
  
  
  
  




!!!..............................................................................
  
  subroutine add_elastic(Rmg,Fel)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none

    real(8),intent(in):: Rmg
    real(8),intent(out):: Fel




    Fel=(Rmg/R0)**2-dlog(Rmg/R0)-1

    Fel=9d0/8d0/pi*(poly%nchains/R0**3)*Fel ! in units of 4/3*pi*R0^3




  end subroutine add_elastic
  







!!!..............................................................................
  
  subroutine nipam(chi,xnipam)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none

    real(8),intent(in):: xnipam
    real(8),intent(out):: chi

    real(8):: g00,g02,g10,g12,g20,g22,g0,g1,g2

    
    g00=-12.947
    g02=0.044959
    g10=17.920
    g12=-0.056944
    g20=14.814
    g22=-0.051419

    g0=g00+g02*temp
    g1=g10+g12*temp
    g2=g20+g22*temp
    

    chi=g0+g1*xnipam+g2*xnipam**2


  end subroutine nipam







!!!..............................................................................
  
  subroutine logistic_chi(chi_i,chi_f,T_pt,deltaT,chi)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none

    real(8),intent(in):: chi_i,chi_f,T_pt,deltaT
    real(8),intent(out):: chi

    
    
    chi=chi_i+(chi_f-chi_i)/(1d0+dexp((t_pt-temp)/deltat))
    

  end subroutine logistic_chi
  




!!!..............................................................................
  
  subroutine sum_FE(Rmg,Df,Fmt,Fel,Ftot)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none

    real(8),intent(in):: Df,Rmg,Fel
    real(8),intent(out):: Ftot,Fmt



    Fmt=Df*(Rmg/R0)**3 ! in units of 4/3*pi*R0^3

    

    Ftot=Fmt+Fel




  end subroutine sum_FE










!!!..............................................................................
  
  subroutine akima(ndat,xdat_in,ydat_in,n_spl,x_spl,y_spl)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    use fvn_interpol
    
    implicit none

    integer,intent(in):: ndat,n_spl

    real(8),intent(in),dimension(ndat):: xdat_in,ydat_in
    real(8),intent(out),dimension(n_spl):: x_spl,y_spl
    
    real(8),dimension(ndat):: xdat,ydat

    integer:: i
    
    real(8):: breakpts(ndat),coeff_fvn(4,ndat)
    real(8):: x,xmin,xmax,a

    
    logical:: mk(ndat)

    mk=.true.


    
    do i=1,ndat ! reorder data
       
       xdat(i)=minval(xdat_in,mk)

       ydat(i)=ydat_in(minloc(xdat_in,1,mk))

       mk(minloc(xdat_in,mk))=.false.
   
    enddo

    

    ! spline.......................................
    
    call fvn_d_akima(ndat,xdat,ydat,breakpts,coeff_fvn)


    xmin=minval(xdat)
    xmax=maxval(xdat)

    a=(xmax-xmin)/dble(n_spl-1)


    do i=1,n_spl
       
       x=a*dble(i-1)+xmin
       
       x_spl(i)=x
       y_spl(i)=fvn_d_spline_eval(x,ndat-1,breakpts,coeff_fvn)
      
    enddo
    
  
    





    

  end subroutine akima
  






















end module simple_microgel
