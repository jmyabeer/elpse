      program initialize_grid
!
      use module_globalp, ONLY: rk
      use module_all_variable_set, ONLY: contape
      use module_globali, ONLY: geometry, number_of_dimensions, &
                                GEO_PLANAR, &
                                GEO_CYLINDRICAL_FULL, &
                                GEO_CYLINDRICAL_HALF, &
                                GEO_CYLINDRICAL_CONE, &
                                GEO_SPHERICAL, &
                                ie, is, je, js, &
                                number_of_materials, mdx, mdy, mdz
      use module_globalr, ONLY: mc_zero, mc_ohalf, mc_one, mc_two, mc_pi
      use module_eos,     ONLY: rho, tele, tion, mzoneid, &
                                meos_mat, mavz_mat, mopac_mat
      use module_eos_mix, ONLY: gammae_mat, gammai_mat, vof, rho_mat
      use module_eos_tab, ONLY: ideal_gas_mode, thomas_fermi_mode, sesame,  &
                                wisconsin, astro, lanl_fit, fully_ionized,  &
                                non_ionized
      use module_material, ONLY: material_id
      use module_dynamic,  ONLY: xl, yl, vxl, vyl
      use module_initial,  ONLY: grid_type, grid_version, initial_grid_parameters
      use module_initial_grid
      use module_plot

      implicit none
!---------------------------------------------------------------------------
      real*8, parameter :: um=1.0d-4
      integer*4 :: i,j,n,nl,iflag,nm,nnm,mat_req0(10),polar_cells,n2_fine, Immesh0
      real*8  :: dtheta,tmp,x1c,x2c,RingR,RingX,RingY,eps,tan1,tan2,factor,d_0,d_1,d_2
      real*8, allocatable :: x1(:), x2(:), x2x(:)
      real(8) :: R,R1,R2, Rtarget0,delx1,dx10,dx1n,dx1f,dRf,Rm,fabl, dR2,dR1,alpha2,Nf, &
                 Rabl,dRstar,dR3,dx1star, a,z0,z, coef, Vmesh, dRmove, Rmmesh0, Rtarget00
      real(8),dimension(1:30) :: adr,aR,aeps,aN
      integer :: nints,int1,int0
      real*8  :: lambda, h, t1, t2 , lambda_min, lambda_fin, lambda_max
      integer :: irt, lmode0, jcells_per_lambda_min
      real*8  :: ampl_corrugation, lambda_corrugation, width_outer_surface_smear_over, rt, t, teleShell, teleHalo
      integer :: nmShell, nmHalo

!---------------------------------------------------------------------------
! User customizable section begins        
!---------------------------------------------------------------------------
      grid_type               = 'user'
      grid_version            = 'with_grid_parameters'
      number_of_dimensions    = 2        

      geometry = GEO_CYLINDRICAL_CONE !GEO_CYLINDRICAL_HALF

!===========================================================================
	  Rtarget0 = 30.0_rk/(mc_two*mc_pi)*2048*10*3 !in microns
      Rtarget00 = Rtarget0

      jcells_per_lambda_min = 1 !use only enough resolution to resolve 30um

      lambda_min = 30.0_rk*um*1.4/2.0/10.
      lambda_fin = 30.0_rk*um*71./2.0/10.*3.4 !80*1.8 !extend the area to 270um
      lambda_max = 30.0_rk*um*100*4.8 !width of the cone at the Rtarget0
      
      ampl_corrugation = 0.1*um*0  !NO CORRUGATION
      lambda_corrugation = 30*um
      width_outer_surface_smear_over = 0.1*um*1.0e-16      !NO SMOOTHING
!===========================================================================

      roffset                 = 0.0     

!
!     ... Be sure to specify rtarget0 in laserins ...
!
!     ... layers are listed from the inside out ...
!

!     ... CH solid ...
      mat_lay (2) = 110
      xlay    (2) = 700.0_rk ! 362.4_rk
      rho_lay (2) = 1.044_rk  
      tev_lay (2) = 0.0252_rk
      eos_lay (2) = sesame          !thomas_fermi_mode, ideal_gas_mode
      opac_lay(2) = astro           !astro, fully_ionized
      avz_lay (2) = astro           !fully_ionized  !astro           

!     ... DD vapor ...
      mat_lay   (1) = 101
      xlay      (1) = Rtarget0 - xlay(2)
      rho_lay   (1) = 0.625e-4_rk 
      tev_lay   (1) = 0.0252_rk
      eos_lay   (1) = sesame        !thomas_fermi_mode, ideal_gas_mode
      opac_lay  (1) = astro         !astro, fully_ionized
      avz_lay   (1) = astro         !fully_ionized  !astro             

!     ... DD vapor ...
      xlay    (3) = 20000.0_rk
      mat_lay (3) = 101
      rho_lay (3) = 2.09e-04_rk*0.01d0
      tev_lay (3) = 0.025_rk * 100d0
      eos_lay (3) = sesame          !thomas_fermi_mode, ideal_gas_mode
      opac_lay(3) = astro           !astro, fully_ionized
      avz_lay (3) = astro           !fully_ionized  !astro           
                                                                                

!
!---------------------------------------------------------------------------
! You can change everything below this line! 
! Only be sure that this will work properly.
!---------------------------------------------------------------------------
!
      xlay(:)=xlay(:)*um            

!
!     count number of layers:
!
      nl = 1
      do while(mat_lay(nl) /= 0)
         nl = nl + 1
      enddo
      nlay = nl - 1
      if(nlay > num_layers) then
        write(contape,*) "Increase num_layers in module_initial_grid.f90"
        stop
      elseif(nlay == 0) then
        write(contape,*) "No layers specified"
        stop
      endif


!
!     ... set azimthal (x2) grid (uniform) ...
!
!===========================================================================
      Rtarget0 = sum(xlay(1:nlay-1),1)      !initial target radius
      
      lmode0 = 4*nint(0.5_rk*mc_pi*Rtarget0/lambda_max) !to simulate half-wavelength
print*,"Rtarget0=",Rtarget0,"cm"	  
print*,"lmode0=",lmode0," corresponds to lambda=",lambda_max/um,"um"      
      cone_angle_min = 0.0_rk
      cone_angle_max = cone_angle_min + mc_two*180.0/real(lmode0,rk) 
!===========================================================================



      js = 3
      ycells = nint(lambda_max/lambda_min)*jcells_per_lambda_min
      je = js + ycells - 1
      allocate(x2(js:je+1))

      if( geometry /= GEO_CYLINDRICAL_CONE )then
         print*,"ERROR{}: only GEO_CYLINDRICAL_CONE is allowed in this initialize_grid.f90 file!"
         stop
      endif

      x2(js) = cone_angle_min*mc_pi/180.0_rk
      dtheta = mc_pi/180.0_rk*(cone_angle_max - cone_angle_min)/ycells
      eps = (lambda_max - lambda_fin)/(lambda_max - lambda_fin+lambda_min/jcells_per_lambda_min-max(500.0_rk*um,lambda_min/jcells_per_lambda_min)) !feather from dtheta to max of 20um or lambda_min/jcells_per_lambda_min
      j = js+1
      do
        x2(j) = x2(j-1)+dtheta
        if( (x2(j)-x2(js))*Rtarget0 >= lambda_max ) exit
        if( (x2(j)-x2(js))*Rtarget0 >= lambda_fin ) dtheta = dtheta*eps
        j = j + 1
      end do
      je = j -1      
      x2(je+1) = cone_angle_max*mc_pi/180.0_rk
!+ashv the following does not allow creating a very thin j=je zone
if( x2(je+1) - x2(je) < 0.3*(x2(je) - x2(je-1)) ) then
    x2(je) = x2(je+1)
    je = je - 1
endif
!-ashv
      ycells = je - js + 1


!
!     ... add points in x2 grid (Igor's original method)...
!
      if(.false.) then
        n2_fine=8 ! 7
        allocate(x2x(js:je+1))
        x2x=x2
        d_0=x2(je+1)-x2(je)
        d_1=d_0/3.0_rk
        d_2=x2(je+1)-x2(je+1-n2_fine)
        factor=(d_2-d_0)/(d_2-d_1)
        n=1
        tmp=d_0
        d_2=x2(js)
        do while (d_2 < x2(je+1))
          n=n+1
          if(d_2 > x2x(je+1-2*n2_fine)) tmp=max(d_1,tmp*factor)
          d_2=d_2+tmp
        end do
        deallocate(x2)
        allocate(x2(js:js+n-1))
        x2(js)=x2x(js)
        do j=js+1,js+n-1
          if(x2(j-1) > x2x(je+1-2*n2_fine)) d_0=max(d_1,d_0*factor)
          x2(j)=x2(j-1)+d_0
        end do
        tmp=x2x(je+1)
        deallocate(x2x)
        je=n-1
        x2(js:je+1)=x2(js:je+1)*tmp/x2(je+1)
      end if

      !Pats non-uniform theta mesh
      !x2(js) = 0.0d0
      !eps = 0.01067
      !dtheta = mc_pi/180._rk
      !do j=js,je-1
      !  x2(j+1)=x2(j)+dtheta
      !  dtheta=dtheta*(mc_one-eps)
      !end do
      !x2(je+1)=mc_ohalf*mc_pi
!
!     ... set radial (x1) grid ...
!

!
!     ... specify moving mesh parameters ...
!
      Rtarget0 = sum(xlay(1:nlay-1),1)       !initial target radius
      Rmmesh0  = Rtarget0 - 810.0_rk*um      !140.0_rk*um

      Rtarget0 = Rtarget0 - Rmmesh0          !need to redefine Rtarger0 relative to offset Rmmesh0
      dx1n     = 15.0_rk*um                  !size of the last (non-moving) cell of the moving grid
      dx1f     = 1.0_rk*um                   !cell size of the fine-mesh
      Rm       = Rtarget0 + 1000.0_rk*um      !outer radius of the moving mesh minus Rmmesh0
      Vmesh    = 1.0_rk                      !sliding grid velocity relative to ablation-surface velocity 0 <= Vmesh <= 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dRf      = 740.0_rk*um    !range occupied by the fine-mesh
      fabl     = 1.0_rk-10.0_rk*um/dRf            !0<fabl<1, fabl*dRf is the supposed position of the ablation surface within the fine mesh region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      dRmove   = 15.0*um


      is = 3
      n=1218 ! 1536 ! 1910 ! 1470 !1314 !3435 !2656
      allocate(x1(is:is+n))
      x1(:)  = mc_zero
!
!     ... specify inner portion of the static mesh going from 0 to Rmmesh0
!
      Immesh0 = is+1
      x1(Immesh0)=Rmmesh0
      Rmmesh0 = x1(Immesh0)
!
!     ... set up the moving mesh (that follows the ablation surface) parameters ...
!
      Nf = dble(ceiling( (Rm-dRf)/(dx1n-dx1f)*log(dx1n/dx1f) + dRf/dx1f ))
      Rm = (Nf-dRf/dx1f)*(dx1n-dx1f)/log(dx1n/dx1f)+dRf
      alpha2 = (dx1n-dx1f)/(Rm-(Rtarget0+dRf*(1.0_rk-fabl)))

!
!     ... save the moving mesh parameters into initial_grid_parameters varible ...
!
      allocate(initial_grid_parameters(10))
      initial_grid_parameters(1) = 1     !version of the moving mesh algorithm to use for eulerian runs
      initial_grid_parameters(2) = Rtarget0
      initial_grid_parameters(3) = dx1n
      initial_grid_parameters(4) = dx1f
      initial_grid_parameters(5) = dRf
      initial_grid_parameters(6) = Rm
      initial_grid_parameters(7) = fabl  
      initial_grid_parameters(8) = Vmesh
      initial_grid_parameters(9) = dRmove
      initial_grid_parameters(10)= Immesh0

      Rabl = Rtarget0 !150.0e-4 !Rtarget0

      dR1 = Rabl - dRf*fabl
      if(dR1<=0) then
         nints=3
         aR(1:nints) =(/0.0_rk,dRf,Rm/)
         adr(1:nints)=(/dx1f,dx1f,dx1n/)
      else
         dRstar = Rm - dR1 - dRf
         dR3 = (alpha2*dRstar-dx1n+dx1f)/(alpha2-(dx1n-dx1f)/(Rm-dRf))
         dR2 = dRstar - dR3
         dx1star = dx1f + alpha2*dR2

         a = Nf - dRf/dx1f - dR2/(dx1star-dx1f)*log(dx1star/dx1f)
         if( dx1n - dx1star > 0 ) a = a - dR3/(dx1n-dx1star)*log(dx1n/dx1star)
         a = a/dR1*dx1f
         !now we need to solve the equation a*(z-1)-Log(z)=0
         if(a>=1.0_rk) then !only one solution z=1
            z=1.0_rk
         else              !two solutions z=1 and z>1 (we want z>1 root)
            z0=1.0_rk/a
            do
               z=log(z0)/a+1
               if( abs(z-z0) < 1e-13_rk*abs(z0) )exit
               z0=z
            enddo
         endif
         dx10=dx1f*z
         nints=5
         aR(1:nints) =(/0.0_rk,dR1,dR1+dRf,dR1+dRf+dR2,Rm/)
         adr(1:nints)=(/dx10,dx1f,dx1f,dx1star,dx1n/)
      endif

      !to get rid of round off noise
      where(aR>Rm)aR=Rm

      Rmmesh0=x1(Immesh0)
      call piecewise_log_mesh(aR,adr,nints,x1(Immesh0),n)
      i=Immesh0+n-1
      x1(Immesh0:i)=x1(Immesh0:i)+Rmmesh0

!
!     ... Finish with points for upper half of ring and Halo region 
!

      R1=x1(i)
      R2=sum(xlay(1:nlay),1) !cm
      dr1=dx1n 
!dr2 is the size of the outermost cell
      dr2=200.0*um

      eps=(R2-R1)/(R2-R1+dr1-dr2)
      delx1=eps*dr1
      do
         if(x1(i)+delx1>R2)exit
         x1(i+1)=x1(i)+delx1
         i=i+1
         delx1=delx1*eps
      enddo
      ie=i-1
!
      call check_geometry
!
!     count number of materials:
!
      number_of_materials = 0
      do nl = 1, nlay
         iflag=0
         do nm = 1, number_of_materials
            if(mat_lay(nl) == mat_req0(nm)) iflag=1
         enddo
         if(iflag == 0) then
           number_of_materials = number_of_materials + 1
           mat_req0(number_of_materials) = mat_lay(nl)
         endif
      enddo
!
!     allocate arrays:
!
      call set_valid_namelists
      call set_variable_attributes
      call initial_alloc_arrays
      call initial_alloc_eos_arrays
      call initial_alloc_eos_mix_arrays
!
!     set material properties:
!
      nnm = 0
      do nl = 1, nlay
         iflag=0
         avz_lay(nl) = opac_lay(nl)
         do nm = 1, number_of_materials
            if(mat_lay(nl) == material_id(nm)) iflag=1
         enddo
         if(iflag == 0) then
           nnm = nnm + 1
           material_id(nnm) = mat_lay(nl)
           meos_mat(nnm)    = eos_lay(nl)
           mopac_mat(nnm)   = opac_lay(nl)
           mavz_mat(nnm)    = opac_lay(nl)
           if (gammae(nnm) /= mc_zero) gammae_mat(nnm) = gammae(nnm)
           if (gammai(nnm) /= mc_zero) gammai_mat(nnm) = gammai(nnm)
         endif
      enddo

!
!     fill the grid
!

      xl     (1:mdx,1:mdy,1:mdz) = mc_zero
      yl     (1:mdx,1:mdy,1:mdz) = mc_zero
      vxl    (1:mdx,1:mdy,1:mdz) = mc_zero
      vyl    (1:mdx,1:mdy,1:mdz) = mc_zero
      rho    (1:mdx,1:mdy,1:mdz) = mc_zero
      tele   (1:mdx,1:mdy,1:mdz) = mc_zero
      tion   (1:mdx,1:mdy,1:mdz) = mc_zero
      mzoneid(1:mdx,1:mdy,1:mdz) = mc_zero
      rho_mat(:,1:mdx,1:mdy,1:mdz) = mc_zero
      vof    (:,1:mdx,1:mdy,1:mdz) = mc_zero

!
!     ... Spherical layers ...
!

!
!     ... Assign rho, tele, tion, and mzoneid
!
      do j=js,je
      do i=is,ie
         x1c=0.5_rk*(x1(i)+x1(i+1))
         if(j.eq.js) print *,i,x1c,x1(i)
         do nl=1,nlay
           if(x1c<sum(xlay(1:nl),1)) exit
         enddo
!         if (((((x2(j)-x2(js))*Rtarget00*um)**2+(x1c-2.933543911069815e1-0.035)**2) >= 0.035**2) .AND. nl==2) then
         if (((((x2(j)-x2(js))*Rtarget00*um)**2+(x1(i)-2.933543911069815e1+0.035)**2) >= 0.035**2) .AND. nl==2) then
         rho      (i,j,1) = rho_lay(3)
         tele     (i,j,1) = tev_lay(3)
         tion     (i,j,1) = tev_lay(3)
         mzoneid  (i,j,1) = mat_lay(3)
         else
         rho      (i,j,1) = rho_lay(nl)
         tele     (i,j,1) = tev_lay(nl)
         tion     (i,j,1) = tev_lay(nl)
         mzoneid  (i,j,1) = mat_lay(nl)
         endif
             
      end do
      end do

!==============================================================
!      do irt=is,ie
!         x1c=0.5_rk*(x1(irt)+x1(irt+1))
!         if( x1c > Rtarget0+x1(Immesh0) )exit
!      enddo   
!      irt = max(is,irt-1)
!      if( irt-1 < is ) then
!         print*,"ERROR{initialize_grid}: irt-1 < is. Terminating..."; stop
!      endif
!
!      h = x1(irt) - x1(irt-1)
!
!      
!      nints = nint(3.0_rk*mc_two*A0/h - mc_one)
!      nints = max(1,nints)
!
!      if( irt-nints < is ) then
!         print*,"ERROR{initialize_grid}: irt-nints < is. Terminating..."; stop
!      endif
!      if(any( x1(irt-nints+1:irt)-x1(irt-nints:irt-1) /=h ))then
!         print*,"ERROR{initialize_grid}: any( x1(irt-nints+1:irt)-x1(irt-nints:irt-1) /= h ). Terminating..."; stop
!      endif
!      if(any( rho(irt-nints:irt,js:je,1) /= rho(irt,js,1) ))then
!         print*,"ERROR{initialize_grid}: any( rho(irt-nints:irt,js:je,1) /= rho(irt,js,1) ). Terminating..."; stop
!      endif
!
!      do j=js,je
!         t1 = mc_two*mc_pi*real(j  -js,rk)/real(je+1-js,rk)
!         t2 = mc_two*mc_pi*real(j+1-js,rk)/real(je+1-js,rk)
!         a = -A0*(sin(t2)-sin(t1))/(t2-t1)
!         do n = 1,nints
!            i = (irt-nints) + n
!            rho(i,j,1) = rho(i,j,1)*(mc_one+a/h*real(n,rk)/real(nints*(nints+1)/2,rk))
!         enddo
!      end do
!
!      do j=js,je
!         t1 = mc_two*mc_pi*real(j  -js,rk)/real(je+1-js,rk)
!         t2 = mc_two*mc_pi*real(j+1-js,rk)/real(je+1-js,rk)
!         a = -A0*(sin(t2)-sin(t1))/(t2-t1)
!         if( abs( sum( rho(irt-nints+1:irt,j,1)/rho(irt-nints,j,1)-mc_one, 1 ) - (a/h)) > 1.0e-10_rk*(a/h))then
!            print*,"ERROR{initialize_grid}: failed the check of the imposed perturbation"
!         endif
!      end do
!===============================================================

!
!     ... Assign material densities and volume fractions
!
      do j=js,je
      do i=is,ie
         rho_mat(:,i,j,1) = mc_zero
         vof    (:,i,j,1) = mc_zero
         where( material_id(1:number_of_materials) == mzoneid(i,j,1) )
            rho_mat(1:number_of_materials,i,j,1) = rho(i,j,1)
            vof    (1:number_of_materials,i,j,1) = mc_one
         endwhere
      end do
      end do

!+++++++++++++++++++++++ Impose the corrugation ++++++++++++++++++++++++++

GOTO 100
      do irt=is,ie
         x1c=0.5_rk*(x1(irt)+x1(irt+1))
         if( x1c > Rtarget0+x1(Immesh0) )exit
      enddo   
      if( irt < is+1 ) then
         print*,"ERROR{initialize_grid}: irt < is+1. Terminating..."; stop
      endif

! real*8  :: ampl_corrugation, lambda_corrugation, width_outer_surface_smear_over, rt, a, t

      nmShell = maxloc(vof(:,irt-1,js,1),1)
      nmHalo  = maxloc(vof(:,irt,  js,1),1)
	  teleShell = tele(irt-1,js,1)
	  teleHalo  = tele(irt  ,js,1)
      
      do j=js,je
         t = 0.5_rk*( x2(j) + x2(j+1) ) - x2(js)
         rt = x1(irt) - ampl_corrugation*cos(mc_two*mc_pi*t*x1(irt)/lambda_corrugation)
         i = irt
         do i=is,ie
            if( x1(i+1) > rt - ampl_corrugation - width_outer_surface_smear_over )exit
         enddo
         do 
            if( x1(i) > rt + ampl_corrugation + width_outer_surface_smear_over .OR. i == ie )exit
            a = 1.0_rk/width_outer_surface_smear_over
            vof(nmHalo, i,j,1) = ( p(a*x1(i+1)-a*rt) - p(a*x1(i)-a*rt) )/( a*x1(i+1) - a*x1(i) )
            vof(nmShell,i,j,1) = 1.0_rk - vof(nmHalo,i,j,1)
            rho_mat(nmHalo, i,j,1) = rho_mat(nmHalo, irt,  j,1)
            rho_mat(nmShell,i,j,1) = rho_mat(nmShell,irt-1,j,1)
            rho(i,j,1) = sum(rho_mat(1:number_of_materials,i,j,1)*vof(:number_of_materials,i,j,1),1)
            !tele and tion modulations
			tele(i,j,1) = ( teleShell*rho_mat(nmShell, i,j,1)*vof(nmShell, i,j,1) + teleHalo*rho_mat(nmHalo, i,j,1)*vof(nmHalo, i,j,1) )/rho(i,j,1)
			tion(i,j,1) = tele(i,j,1)
            !uniform tele and tion smearing
            a = 1.0_rk/(width_outer_surface_smear_over+ampl_corrugation)
            tele(i,j,1) = teleShell + (( p(a*x1(i+1)-a*x1(irt)) - p(a*x1(i)-a*x1(irt)) )/( a*x1(i+1) - a*x1(i) ))*(teleHalo-teleShell)
			tion(i,j,1) = tele(i,j,1)

            i = i + 1
         enddo
      end do
100 CONTINUE
      
!-------------------------------------------------------------------------
      

!
!     ... Assign node coordinates and velocities
!
      do i=is,ie
        if(number_of_dimensions == 1) then
          xl (i,1,1)=x1(i)
          vxl(i,1,1)=mc_zero
        else
          do j=js,je+1
            xl (i,j,1)=x1(i)*cos(x2(j))
            yl (i,j,1)=x1(i)*sin(x2(j))
            vxl(i,j,1)=mc_zero
            vyl(i,j,1)=mc_zero
          end do
        end if
      end do
      if(number_of_dimensions == 1) then
        xl (ie+1,1,1)=x1(ie+1)
        vxl(ie+1,1,1)=mc_zero
      else
        do j=js,je+1
          xl(ie+1,j,1)=x1(ie+1)*cos(x2(j))
          yl(ie+1,j,1)=x1(ie+1)*sin(x2(j))
          vxl(ie+1,j,1)=mc_zero
          vyl(ie+1,j,1)=mc_zero
        end do
      end if
!
      call check_grid
!
      call write_initial_grid
!
      
      contains
      
      function UnitStep(x)
        real(rk):: unitstep, x
        UnitStep = 0.5_rk*(sign(1.0_rk,x)+1.0_rk)
      end function UnitStep
      
      function p(x)
        real(rk) :: p,x
        p = ((-1 + x)**3*UnitStep(-1 + x) - 2*x**3*UnitStep(x) + (1 + x)**3*UnitStep(1 + x))/6.0_rk
      end function p
      
      end program initialize_grid
