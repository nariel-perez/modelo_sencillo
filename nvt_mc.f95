!!!*****************************************************************************
!!!*****************************************************************************

module nvt_mc
  
!!!*****************************************************************************
!!!*****************************************************************************
  
  use parameters
  
  implicit none
  
  
contains
  

  
  
  
!!!..............................................................................
  
  subroutine monte_carlo
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none
    
    integer,allocatable:: seed(:)    
    
    integer:: imc,ip,i,j,ig,ntrials,nplot,nha,iq,nq
    integer:: att,acc,counter,nav,nwidom
    integer:: ianew,iatrial,iamin

    integer,allocatable:: ia(:)    
    
    real(8),allocatable:: g(:),x(:),y(:),z(:),a(:),ha(:),q(:),Sq(:,:)
    
    real(8):: rn,rn4(4)
    real(8):: xnew,ynew,znew,maxd
    real(8):: u,uold,unew,uij,vij,deltau
    real(8):: ratio,wtest,uav,u2av
    real(8):: xtrial,ytrial,ztrial,atrial,utrial,qtrial
    real(8):: v,vb,nid,r,rij,dr,rcut,dq
    real(8):: p,bmu,vol,box
    real(8):: a0,sig,anew,da,ava,a2av,qnew,qaux,maxda
    real(8):: df,bf,bfold,bfnew,duf,bfav,bf2av,bf0
    real(8):: zeta,zetaav,zeta2av,kappa,gamma_av

    
    character(8):: chx8
    character(12):: chx12,chx12b

    integer:: ratio2
    
    


    call random_seed(size=i)
    allocate(seed(i))
    call random_seed(get=seed)

    seed=999
    call random_seed(put=seed)



    allocate(x(Np),y(Np),z(Np),a(Np),q(Np),ia(Np))



!!!....................


    iamin=minloc(F_mg,1)
    a0=R_mg(iamin)
    
    bf0=F_mg(iamin)

    F_mg=F_mg-bf0


!!!....................



    call set_MC_param(box,vol,sig,rcut,dr,nplot,da,nha,maxd,nq,dq,a0,maxda)

    
    kappa_sol=kappa_sol*sig
    
    lb_sol=lb_sol/sig
    


!!!....................

  
    !a=a0/sig
    a=a_min
    
    call get_ia(a(1),sig,i)

    ia=i
   
    call DF_Da(ia(1),bf)

    bf=bf*dble(Np)

    call q_net(ia(1),qaux)

    q=qaux




    !bf=0d0    
    !do ip=1,Np
       !call random_number(rn)

       !a(ip)=a_min+rn*(a_max-a_min)
       
       !call DF_Da(a(ip),sig,df,q(ip))
       
       !bf=bf+df
    !enddo
  


!!!....................



    call initial_coordinates(x,y,z,Np,box)



    call energy(x,y,z,a,q,Np,u,box,rcut)

    
    


    write(*,'(a)')"--- Initial Values:"

    write(chx12,"(1pe11.3e2)")real(u/dble(Np))
    write(*,'(4xa,9xa10)')"Energy pp:",trim(adjustl(chx12))

    write(chx12,"(1pe11.3e2)")real(bf/dble(Np))
    write(*,'(4xa,1xa10)')"MG Free Energy pp:",trim(adjustl(chx12))

    write(chx12,"(1pe11.3e2)")real((bf+u)/dble(Np))
    write(*,'(4xa,13xa10/)')"Total:",trim(adjustl(chx12))



    write(chx12,"(f8.1)")R0
    write(chx8,"(f8.1)")a0

    write(*,'(4xa,1xa,1xa,1xa,1xa)')"MG radius, R0:",trim(adjustl(chx12)),&
         "nm, a0:",trim(adjustl(chx8)),"nm"

    write(chx12,"(f8.1)")2d0*a0
    write(*,'(4xa,1xa,1xa)')"diameter:",trim(adjustl(chx12)),"nm"



    call q_net(iamin,qaux)

    !kappa=kappa_mg(iamin)*sig

    call zeta_pot(iamin,a0/sig,qaux,zeta)

    write(chx12,"(f10.3)")zeta
    write(*,'(4xa,1xa,1xa)')"zeta pot:",trim(adjustl(chx12)),"mV"


    write(chx8,"(f8.1)")sig
    write(*,'(4xa,1xa,1xa)')"sig: 2*R0:",trim(adjustl(chx8)),"nm"


    write(chx12,"(1pe12.3e2)")pi*dble(Np)/6d0/vol
    write(*,'(4xa,1xa)')"vol. fraction:",trim(adjustl(chx12))
    

    write(chx8,"(f8.4)")dens
    write(*,'(4xa,1xa)')"dens*:",trim(adjustl(chx8))
    write(chx8,"(f8.4)")dens*(a0/R0)**3
    write(*,'(4xa,1xa)')"dens* (a0/R0)^3:",trim(adjustl(chx8))

    write(chx8,"(f8.4)")box
    write(*,'(4xa,1xa)')"box*:",trim(adjustl(chx8))
   
    

    write(*,'(/4xa/)')"...starting MC simulation:"



    
    open(unit=100,file="upot.dat")
    open(unit=120,file="traj.xyz")
    
    

    acc=0
    att=0
    counter=0
    nwidom=0


    nav=0
    uav=0d0
    u2av=0d0

    bfav=0d0
    bf2av=0d0

    ava=0d0
    a2av=0d0

    zetaav=0d0
    zeta2av=0d0

    gamma_av = 0d0
    wtest=0d0
    v=0d0
    

    allocate(g(nplot))
    g=0d0


    allocate(Sq(nq,2))
    Sq=0d0

    do iq=1,nq
       Sq(iq,1)=dble(iq-1)*dq
    enddo


    allocate(ha(nha))
    ha=0d0
    
    ratio2 = 0
    do imc=1,max_mc ! loop over the # of MC steps
       
       do ip=1,Np ! MC step --> Np attempts to move randomly chosen particles
          
          call random_number(rn)
          
          i=int(rn*Np)+1
          
          att=att+1
          
          call random_number(rn4)
          
          xnew=x(i)+(rn4(1)-0.5)*maxd ! try new positions
          ynew=y(i)+(rn4(2)-0.5)*maxd
          znew=z(i)+(rn4(3)-0.5)*maxd


          !anew=a_min+rn4(4)*(a_max-a_min)

          anew=a(i)+(rn4(4)-0.5)*maxda
          
          if (anew>a_max) anew=a_max
          if (anew<a_min) anew=a_min
          
          
          
          xnew=xnew-box*anint(xnew/box)
          ynew=ynew-box*anint(ynew/box)
          znew=znew-box*anint(znew/box)


          call DF_Da(ia(i),bfold)

          call get_ia(anew,sig,ianew)
          call DF_Da(ianew,bfnew)
          call q_net(ianew,qnew)

          df=bfnew-bfold


          uold=0d0 ! save the energy of particle i in old position
          unew=0d0

          do j=1,Np
             if (i/=j) then
                
                call pair_pot(x(i),y(i),z(i),a(i),q(i),&
                     x(j),y(j),z(j),a(j),q(j),rij,uij,vij,box,rcut)
                uold=uold+uij
                
                call pair_pot(xnew,ynew,znew,anew,qnew,&
                     x(j),y(j),z(j),a(j),q(j),rij,uij,vij,box,rcut)
                unew=unew+uij

             endif
          enddo

          deltau=unew-uold
  
          duf=df+deltau          

           
          call random_number(rn)

          if (duf<20.) then
             if (duf<=0) then
                u=u+deltau
                bf=bf+df
                acc=acc+1
                x(i)=xnew
                y(i)=ynew
                z(i)=znew
                a(i)=anew
                q(i)=qnew
                ia(i)=ianew
                ratio2 = ratio2 +1
             elseif (dexp(-duf)>rn) then
                u=u+deltau
                bf=bf+df
                acc=acc+1
                x(i)=xnew
                y(i)=ynew
                z(i)=znew
                a(i)=anew
                q(i)=qnew
                ia(i)=ianew
                ratio2 = ratio2 +1
             endif
          endif


       enddo



       write(100,*)imc,(bf+u)/dble(Np),bf/dble(Np),u/dble(Np)



       if (mod(imc,1000)==0) then          
 
          write(120,*)Np
          write(120,*)
          do ip=1,Np
             call pbc(x(ip),y(ip),z(ip),box)
             
             write(120,*)"LJ",x(ip),y(ip),z(ip)
          enddo
          
          
          write(chx12,"(1pe12.3e2)")real(u/dble(Np))
          write(chx12b,"(1pe12.3e2)")real(bf/dble(Np))

          write(*,'(4xa,i9,2(1xa,1xa10))')"step:",imc,"/ u:",&
               trim(adjustl(chx12)),"/ f_MG:",trim(adjustl(chx12b))
          
       endif




       !............................................................

       
       
       if (mod(imc,100)==0) then ! correct maxd

          ratio=dble(acc)/dble(att)
          if (ratio>0.45) then
             maxd=maxd*1.05
             maxda=maxda*1.01
          endif
          if (ratio<0.45) then
             maxd=maxd*0.95
             maxda=maxda*0.99
          endif
          if (maxd>box) maxd=box
          if (maxda>a_max-a_min) maxda=a_max-a_min
       endif
       

       
       !............................................................


       
       if (imc>=eq_mc) then


          if (mod(imc,50)==0) then     ! take values every 50 MC steps
             uav=uav+u
             u2av=u2av+u*u
 
             bfav=bfav+bf
             bf2av=bf2av+bf*bf             
 
             ava=ava+sum(a)
             a2av=a2av+sum(a*a)
             
             nav=nav+1
             
             do i=1,Np
               
                !kappa=kappa_mg(ia(i))*sig
                
                call zeta_pot(ia(i),a(i),q(i),zeta)
                if (nads>0) then
                   gamma_av= gamma_av + gamma_ads(ia(i))
                end if
                
                zetaav=zetaav+zeta
                zeta2av=zeta2av+zeta*zeta
                
             enddo
             
             
          endif
         
         
          
          ! Virial, Pressure & g(r) ...................................
         
          if (mod(imc,50)==0) then
             counter=counter+1
             
             do i=1,Np-1
                do j=i+1,Np
                   
                   call pair_pot(x(i),y(i),z(i),a(i),q(i),x(j),y(j),z(j),&
                        a(j),q(j),rij,uij,vij,box,rcut)
                   
                   v=v+vij

                   if (rij<box/2.) then              
                      ig=int(rij/dr)           
                      g(ig)=g(ig)+2 ! g(r)

                      Sq(2:nq,2)=Sq(2:nq,2)+sin(Sq(2:nq,1)*rij)/Sq(2:nq,1)/rij
                      Sq(1,2)=Sq(1,2)+1d0
                      
                   endif

                enddo
             enddo


             do i=1,Np
                j=int((a(i)-a_min)/da)+1
                ha(j)=ha(j)+1d0
             enddo
             
          endif



          !............................................................





          ! Chemical Potential  .......................................

          ntrials=3
          do i=1,ntrials ! I have no idea if this is right
             
             call random_number(rn4)

             xtrial=(rn4(1)-0.5)*box ! trail positions
             ytrial=(rn4(2)-0.5)*box
             ztrial=(rn4(3)-0.5)*box

             atrial=a_min+rn4(4)*(a_max-a_min)

             call get_ia(atrial,sig,iatrial)

             call Df_Da(iatrial,df)
             call q_net(iatrial,qtrial)

             
             utrial=df ! bf_min=0
             
             do j=1,Np
                
                call pair_pot(xtrial,ytrial,ztrial,atrial,qtrial,&
                     x(j),y(j),z(j),a(j),q(j),rij,uij,vij,box,rcut)  
                
                utrial=utrial+uij

             enddo

             wtest=wtest+dexp(-utrial) 
             nwidom=nwidom+1

          enddo
          

          !............................................................




       endif

    enddo



    close(100)
    close(120)





!!! save final configuration
    
    open(unit=121,file="final_config.xyz")    
    write(121,*)Np
    write(121,*)
    do ip=1,Np
       call pbc(x(ip),y(ip),z(ip),box)
             
       write(121,*)"LJ",x(ip),y(ip),z(ip),a(ip)
    enddo
    close(121)

    write(*,*)'aceptados:',ratio2

    
!!!................................................
!!!................................................



    

    open(unit=111,file="size.dat")
    write(111,'(a,6xa,15xa)')"#","MG radius","probability"
    
    do i=1,nha
       write(111,*)((i-1)*da+a_min)*sig,ha(i)/dble(Np)/dble(counter)/da/sig
    enddo
    
    close(111)




    open(unit=110,file="rdf.dat")

    do i=1,nplot
       r=0
       r=dr*(i+0.5)
       vb=((i+1)**3.-i**3.)*dr**3.
       nid=(4./3.)*pi*vb*dens
       g(i)=g(i)/(dble(counter*Np)*nid)

       !write(110,*)r*sig,g(i)
       write(110,*)r,g(i)
    enddo

    close(110)




   open(unit=112,file="Sq.dat")

    do i=1,nq
       write(112,*)Sq(i,1),1d0+2d0*Sq(i,2)/dble(counter*Np)
    enddo

    close(112)









!!!................................................
!!!................................................
!!!................................................

!!!!!!!!!! check from here !!!!!!!!!!
!!! temp?
!!! cv=.... temp?
!!! bmu=.... temp?

!!! remember that everything (u,bf) is in 4/3*pi*R0^3 kBT units

!!!................................................
!!!................................................
!!!................................................



    wtest=wtest/dble(nwidom)
    !bmu=-temp*dlog(wtest)
    bmu=-dlog(wtest) ! excess 


    v=v/dble(counter)
    p=dens*temp+v/(3d0*vol)


    uav=uav/dble(nav)
    u2av=u2av/dble(nav)


    !Cv=(u2av-uav*uav)/temp**2/dble(Np) ! residual, in kB units


    uav=uav/dble(Np)
    u2av=u2av/dble(Np*Np)

    bfav=bfav/dble(nav*Np)
    bf2av=bf2av/dble(nav*Np*Np)


    ava=ava/dble(nav*Np)*sig
    a2av=a2av/dble(nav*Np)*sig*sig ! notice Np not Np**2 is correct here!
    

    zetaav=zetaav/dble(nav*Np)
    zeta2av=zeta2av/dble(nav*Np)
    gamma_av=gamma_av/dble(nav*Np)



!!!................................................
!!!................................................
!!!................................................



    write(chx12,"(1pe12.3e2)")maxd
    write(chx12b,"(1pe12.3e2)")box
    write(*,'(/4xa,1xa,1x"("a,")")')"max. displacement* (box):",&
         trim(adjustl(chx12)),trim(adjustl(chx12b))

    write(chx12,"(1pe12.3e1)")maxda*sig
    write(chx12b,"(1pe12.3e1)")(a_max-a_min)*sig
    write(*,'(4xa,1xa," nm (",a," nm)")')"max. swelling:",&
         trim(adjustl(chx12)),trim(adjustl(chx12b))

    open(unit=2121,file="averages.dat")
    
    write(*,'(/a/)')"--- Average Values:"
    write(2121,'(/a/)')"--- Average Values:"
   
    write(chx12,"(1pe12.5e2)")uav
    write(*,'(4xa,1xa,1xa)')"u:",trim(adjustl(chx12)),"iu"

    write(chx12,"(1pe12.5e2)")sqrt(u2av-uav*uav)
    write(*,'(4xa,1xa,1xa/)')"Du:",trim(adjustl(chx12)),"iu"



    write(chx12,"(1pe12.5e2)")bfav
    write(*,'(4xa,1xa,1xa)')"f_MG:",trim(adjustl(chx12)),"iu"
    write(2121,'(4xa,1xa,1xa)')"f_MG:",trim(adjustl(chx12)),"iu"

    !write(chx12,"(1pe12.5e2)")bf2av
    !write(*,'(4xa,1xa)')"<F_MG^2>:",trim(adjustl(chx12))

    write(chx12,"(1pe12.5e2)")sqrt(bf2av-bfav*bfav)
    write(*,'(4xa,1xa,1xa/)')"Df_MG:",trim(adjustl(chx12)),"iu"


 
    write(chx12,"(1pe12.5e2)")(bfav+uav)
    write(*,'(4xa,1xa,1xa/)')"u+f_MG:",trim(adjustl(chx12)),"iu"




    write(chx12,"(f8.1)")ava
    write(*,'(4xa,1xa,1xa)')"<R>:",trim(adjustl(chx12)),"nm"
    write(2121,'(4xa,1xa,1xa)')"<R>:",trim(adjustl(chx12)),"nm"
   
    write(chx12,"(f8.1)")sqrt(a2av-ava*ava)
    write(chx12b,"(f8.1)")sqrt(a2av-ava*ava)/ava*100
    write(*,'(4xa,2xa,1xa,1xa,1xa)')"DR:",trim(adjustl(chx12)),&
         "nm -",trim(adjustl(chx12b)),"%"
    
    write(chx12,"(f8.1)")2d0*ava
    write(*,'(4xa,1xa,1xa/)')"diameter:",trim(adjustl(chx12)),"nm"
    write(2121,'(4xa,1xa,1xa/)')"diameter:",trim(adjustl(chx12)),"nm"







    write(chx12,"(f12.3)")zetaav
    write(*,'(4xa,1xa,1xa)')"<zeta>:",trim(adjustl(chx12)),"mV"
    write(2121,'(4xa,1xa,1xa)')"<zeta>:",trim(adjustl(chx12)),"mV"

    if (nads>0) then
       write(chx12,"(f12.3)")gamma_av
       write(*,'(4xa,1xa,1xa)')"<gamma>:",trim(adjustl(chx12)),"molec"
       write(2121,'(4xa,1xa,1xa)')"<gamma>:",trim(adjustl(chx12)),"molec"
    end if

    
    !write(chx12,"(1pe12.2e2)")sqrt(zeta2av-zetaav*zetaav)
    write(chx12,"(f12.3)")sqrt(zeta2av-zetaav*zetaav)
    write(chx12b,"(f8.1)")sqrt(zeta2av-zetaav*zetaav)/abs(zetaav)*100d0
    write(*,'(4xa,1xa,1xa,1xa,1xa/)')"Dzeta:",trim(adjustl(chx12)),&
         "mV -",trim(adjustl(chx12b)),"%"


    

    close(2121)


    





  end subroutine monte_carlo
  
  
  
  
  
  
  
  














!!!..............................................................................
  
  subroutine set_MC_param(box,vol,sig,rcut,dr,nplot,da,nha,maxd,nq,dq,a0,maxda)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none
    
    real(8),intent(in):: a0

    integer,intent(out):: nplot,nha,nq
    
    real(8),intent(out):: box,vol,rcut,dr,maxd,sig,da,dq,maxda



    ! distances in sig units    
    ! energies in 4/3*pi*R0^3 kBT units
    
   
    sig=2d0*R0


    dens=Na*dens/Mw_pol*sig**3

    
    vol=dble(Np)/dens

    box=vol**(1d0/3d0)
    
    maxd=0.1 ! in sig units

    
    a_max=a_max/sig
    a_min=a_min/sig


    maxda=a0/sig/10d0
    


    nha=50
    da=(a_max-a_min)/dble(nha-1)

    

    nplot=500    
    dr=box/(nplot*2.)


    
    nq=500
    dq=2d0*pi/box/10d0
    


    rcut=0d0 !!! obsolete in this code!
    
    


    
    
    
    
  end subroutine set_MC_param

  
  
  








!!!..............................................................................
  
  subroutine initial_coordinates(x,y,z,Np,box)
    
!!!*****************************************************************************
!!!******************************************************************************
    
    implicit none
    
    integer,intent(in):: Np
    
    real(8),intent(in):: box
    real(8),intent(out):: x(Np),y(Np),z(Np)
    
    integer:: i,j,k,l,naux
    real(8):: aux
    
    
    
    naux=int(Np**(1d0/3d0))+1
    aux=box/naux
     
    l=0
    do i=1,naux
       do j=1,naux
          do k=1,naux
             l=l+1

             if (l>Np) cycle

             x(l)=aux/2+aux*(i-1)
             y(l)=aux/2+aux*(j-1)
             z(l)=aux/2+aux*(k-1)

          enddo
       enddo
    enddo


  end subroutine initial_coordinates
  











!!!..............................................................................
   
  subroutine energy(x,y,z,a,q,Np,u,box,rcut)
     
!!!*****************************************************************************
!!!******************************************************************************
!!! This subroutine calculates the total potential energy

    implicit none
    
    integer,intent(in):: Np
    real(8),intent(in):: x(Np),y(Np),z(Np),a(Np),q(Np),box,rcut
     
    real(8),intent(out):: u

    integer:: i,j
    real(8):: rij,uij,vij

    
    
    u=0d0
    
    do i=1,Np-1
       do j=i+1,Np
           
          call pair_pot(x(i),y(i),z(i),a(i),q(i),x(j),y(j),z(j),a(j),q(j),&
               rij,uij,vij,box,rcut)
           
          u=u+uij

       enddo
    enddo
     
     
     
  end subroutine energy












!!!..............................................................................
  
  subroutine pair_pot(x1,y1,z1,a1,q1,x2,y2,z2,a2,q2,r,u12,v12,box,rcut)

!!!*****************************************************************************
!!!******************************************************************************
!!! This subroutine defines the pair potential

    implicit none

    real(8),intent(in):: x1,y1,z1,x2,y2,z2,box,rcut,a1,a2,q1,q2
    real(8),intent(out):: r,u12,v12




 


    call Herts_Yukawa(x1,y1,z1,a1,q1,x2,y2,z2,a2,q2,r,u12,v12,box)

    
    


    
    

    
  end subroutine pair_pot
  
  
  
  
  


   










!!!..............................................................................
  
  subroutine Herts_Yukawa(x1,y1,z1,a1,q1,x2,y2,z2,a2,q2,r,u12,v12,box)
    
!!!*****************************************************************************
!!!******************************************************************************
!!! This subroutine...
    
    implicit none
    
    real(8),intent(in):: x1,y1,z1,x2,y2,z2,a1,a2,q1,q2,box
    real(8),intent(out):: r,u12,v12
    real(8):: dx,dy,dz,r2,b12,nu
    real(8):: u_hz,v_hz,u_ykw,v_ykw
     

    u12=0d0
    v12=0d0


    dx=x2-x1
    dy=y2-y1
    dz=z2-z1

    dx=dx-box*anint(dx/box)
    dy=dy-box*anint(dy/box)
    dz=dz-box*anint(dz/box)

    r2=dx*dx+dy*dy+dz*dz

    r=sqrt(r2)
  

    


!!! Hertz potential ..........


    nu=0.5


    u_hz=0d0
    v_hz=0d0
    

    if (r<a1 +a2) then

       b12=6d0/5d0/pi*poly%nchains/units/(1d0-nu**2)*(a1+a2)**2*sqrt(a1*a2)/&
            (a1**3+a2**3)
       
       u_hz=(1d0-r/(a1+a2))**(5d0/2d0)*b12

       v_hz=5d0/2d0*r/(a1+a2)*(1d0-r/(a1+a2))**(3d0/2d0)*b12

    endif
    





!!! Yukawa potential ..........
    
    u_ykw=0d0
    v_ykw=0d0

    !if (r>a1+a2) then

    u_ykw=(q1*q2)/units*dexp(kappa_sol*(a1+a2-r))/(1d0+kappa_sol*a1)/&
         (1d0+kappa_sol*a2)/r

    v_ykw=u_ykw*(1d0+kappa_sol*r)/r
    
    !endif


!!! .......... ..........


    u12=u_ykw+u_hz

    v12=v_ykw+v_hz


     







  end subroutine Herts_Yukawa








!!!.............................................................................
   
  subroutine get_ia(a,sig,ia)
     
!!!*****************************************************************************
!!!*****************************************************************************
!!! This finds the free energy of a given particle radius
    
    implicit none
     
    real(8),intent(in):: a,sig
    
    integer,intent(out):: ia

    real(8):: R


    R=a*sig
    
    ia=minloc(abs(R_mg-R),1)
    
    
  end subroutine get_ia
  
  





!!!..............................................................................
   
  subroutine Df_Da(ia,df)
     
!!!*****************************************************************************
!!!******************************************************************************
!!! This finds the free energy of a given particle radius
    
    implicit none
     
    integer,intent(in):: ia
    
    real(8),intent(out):: df


    df=F_mg(ia)


  end subroutine Df_Da












!!!..............................................................................
  
  subroutine q_net(ia,qnet)
    
!!!*****************************************************************************
!!!******************************************************************************
!!! This subroutine ... linear-response theory
    
    implicit none

    integer,intent(in):: ia
    real(8),intent(out):: qnet

    real(8):: kappa_a,q
    



    kappa_a=kappa_mg(ia)*R_mg(ia) ! this is kappa*a
    !kappa_a=kappa_sol*R_mg(ia)


    q=Q_mg(ia)



    qnet=(1+kappa_a)*dexp(-kappa_a)*3d0*q/kappa_a**2*&
         (cosh(kappa_a)-sinh(kappa_a)/kappa_a)

    !qnet= 0.1*q
    
    
    
  end subroutine q_net
  





!!!.............................................................................
  
  subroutine zeta_pot(ia,a,qnet,zeta)
    
!!!*****************************************************************************
!!!*****************************************************************************
    implicit none

    integer,intent(in):: ia
    real(8),intent(in):: a,qnet
    real(8),intent(out):: zeta

    real(8),parameter:: c_mV=0.08617 ! kB/e in mV/K

    real(8):: z0,ka



    ka=a*kappa_sol

    !ka=kappa_mg(ia)*R_mg(ia) 


    z0=c_mV*lb_sol*temp/a*qnet




    zeta=z0/(1d0+ka)

    !zeta=z0*dexp(-ka)

    
    
    

  end subroutine zeta_pot






















!!!..............................................................................
  
  subroutine pbc(x,y,z,box)
    
!!!*****************************************************************************
!!!******************************************************************************
!!! ...
    
    implicit none

    real(8),intent(in):: box
    real(8),intent(inout):: x,y,z


    if (x<0d0) x=x+box
    if (y<0d0) y=y+box
    if (z<0d0) z=z+box

    if (x>=box) x=x-box
    if (y>=box) y=y-box
    if (z>=box) z=z-box





  end subroutine pbc











 end module nvt_mc
