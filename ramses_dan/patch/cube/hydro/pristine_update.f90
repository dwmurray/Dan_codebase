!################################################################
! This module defines the vars that are used to compute an
! update to the pristine fraction
! Rick Sarmento - Jan 2014
!################################################################
module prist_commons
  use amr_commons
  integer::outputcnt = 0       ! Counter used to print output
  integer::ncache,igrid        ! # of grids, index into the grids being processed
  integer::debugPrint
  real(dp)::Pgas               ! The pristrine fraction scalar
  real(dp)::Pbefore
  real(dp)::Pafter
  real(dp)::opf                ! The original pristine fraction scalar - before update
  real(dp)::totEnergy          ! Gas total energy
  real(dp)::kinEnergy          ! Kinetic energy of the gas
  real(dp)::thermEnergy        ! Gas thermal energy
  real(dp)::n                  ! Scale factor in the pristine update eqn
  real(dp)::tau                ! Mixing time scale 
  real(dp)::turbVel            ! Turbulent velocity
  real(dp)::mach               ! Turb mach #
  real(dp)::m                  ! fit value for mach
  real(dp)::temp               ! Gas temperature
  real(dp)::temp2              ! Gas temperature directly...
  real(dp)::tproper            ! Proper time
  real(dp)::Cs                 ! Sound speed
  real(dp)::scale_nH           ! Converts rho in user units into nH in H/cc
  real(dp)::scale_T2           ! Converts (P/rho) in user unit into (T/mu) in Kelvin
  real(dp)::scale_l            ! Converts distance from user units into cm
  real(dp)::scale_d            ! Converts mass density from user units into g/cc
  real(dp)::scale_t            ! Converts time from user units into seconds
  real(dp)::scale_v            ! Converts velocity in user units into cm/s
  real(kind=8)::dx,dx_loc,scale,vol_loc  ! Distance and vol
  real(kind=8)::nH             ! Number density of H
  integer,dimension(1:nvector)::ind_cell,ind_leaf
end module prist_commons

!################################################################
!################################################################
! This subroutine is called every course step
! Rick Sarmento - Jan 2014
!################################################################
!################################################################
subroutine pristine_update
  use prist_commons
  use hydro_commons
  use amr_commons
  use pm_commons ! For star particle stuff...
  integer::ilevel
  !-------------------------------------------------------------------
  ! Update the pristine gas fraction for the cells
  !-------------------------------------------------------------------
  integer::i,ngrid
  integer,dimension(1:nvector),save::ind_grid

  ! Rick Sarmento
  ! See amr_commons for definition of 'active'
  ! active%ngrid field is number of grids - 
  ! each grid can have 8 octs (cells)
  ! active%igrid is array of ptrs to grids
  ! ------------------------------------------------------------------
  ! The outter loop, we look over all levels... 
  ! For the igrid loop, steps over grids in blocks of nvector, 
  ! Defaults to 500; can be set in the Makefile
  ! For each grid, get the indices of the active grids at that level
  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     if(verbose)write(*,111)ilevel
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector ! step size is nvector
        ngrid=MIN(nvector,ncache-igrid+1) ! process nvector @ a time, unless < nvector left
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1) ! array of points to ngrid grids
        end do
        call pristine_up(ind_grid,ngrid,ilevel) ! Process ngrid (=nvector) grids
     end do
  end do
111 format('   Entering pristine_update for level',i2)

end subroutine pristine_update

!################################################################
!################################################################
!################################################################
!################################################################
subroutine pristine_up(ind_grid,ngrid,ilevel)
  use hydro_commons
  use prist_commons
  use pm_commons ! For star particles -- v_t debugging
  implicit none
  integer::ilevel,ngrid ! level of the grids and number of grids
  ! ind_grid is an array of pointers to the grids we need to process
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer,dimension(1:nvector,0:twondim)::ind_father
  integer,dimension(1:nvector,0:twondim)::nbor_level
  real(dp),dimension(1:twotondim,1:3)::xc      ! Cell centers relative to grid center
  integer::i,j,k,ind,indx, iskip,idim,nleaf
  integer::nx_loc,ix,iy,iz,ivar
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:ndim,1:ndim)::Aij,Sij   ! The vel derivatives, & strain tensor
  ! vel's first dim includes 7 items (0-6) for the central cell and the 6 surrounding
  real(dp),dimension(0:twondim,1:ndim)::vel    ! The vel of the gas in cell (0,*) and nbors 
  real(dp)::ev, rightLevel, leftLevel
  real(dp)::tnorm            ! normalize time step: n * dtnew/tau
  real(dp)::dtnorm           ! norm divided into sub-steps
  real(dp)::Zcrit            ! Z_critical as a mass frac (not solar)
  real(dp)::Zave                                 
  real(dp)::x                ! temp variable used with tau
  real(dp)::xcoord,ycoord,zcoord
  integer::nsteps            ! number of sub-cycle steps
  real(dp)::tempDecay        ! DEBUGGING....
  real(dp)::threshold        ! Used to determine if we decay Pgas or increase it
  integer::ind_part

  ! Mesh spacing in that level
  ! Remember, the boxlen = 1.0, so @ level 1
  ! we have 2 boxes (along each axis) each 0.5 on side.
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale  = boxlen/dble(nx_loc) ! scale = 1.0 for my cosmo runs
  dx_loc = dx * scale
  vol_loc=dx_loc**ndim
  ! if ( nx_loc > 0) then
  !    print *,"icoarse_max",icoarse_max
  !    print *,"icoarse_min",icoarse_min
  !    print *,"dx @ ilevel",dx, ilevel
  !    print *,"nx_loc @ icoarse_max, icoarse_min", nx_loc, icoarse_max,icoarse_min
  !    print *,"scale: ",scale
  !    print *,"dx_loc: ", dx_loc
  !    print *,"vol_loc: ", vol_loc
  ! end if

  ! Cells center position relative to grid center position
  do ind=1,twotondim  ! 1 -> 2^3 = 8
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Compute Zcrit as mass fraction
  ! Move to above -- don't need to recompute each time ***
  Zcrit = z_crit * 0.02d0 ! z_crit is in solar units, z_crit in amr_parameters

  ! Loop over cells in each grid
  ! There are 8 cells to a grid (for 3D)
  do ind=1,twotondim
     ! ncoarse = nx.ny.nz,  ! n's are # coarse cells in each dir
     iskip=ncoarse+(ind-1)*ngridmax  
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)  
     end do
     ! At this point we have a list of indices for all cells for our grids.
     ! grids (in array ind_grid.

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i) ! Finding leaf cells - no children
        end if
     end do

     ! Find the neighboring cells/grids for each leaf
     ! ind_leaf is an array of indices of cells we're interested in...  
     ! ind_father is the retuned list such that
     ! ind_father(i,0)=ind_leaf(i)
     ! ind_father(i,j) with j -> {1..6} are the indicies of the nbors of 
     ! ind_leaf(i) [or ind_father of the nbor if a cell doesn't exist] ... 
     ! nbor_level(i,j) is the level of that nboring cell [or grid]
     ! Normally, a nboring cell is returned (at same level), if not,
     ! the next coarser level nbor is returned and ilevel for that 
     ! nbor is set to ilevel - 1 
     call getnborfatherwithlevel(ind_leaf,ind_father,nbor_level,nleaf,ilevel)
!!$     call printNeighbors(ind_father,nbor_level,nleaf,ilevel)

     ! ******************** MAIN LOOP ********************
     ! Update the pristine fraction of gas in the cell: uold(*,iprist)
     ! Mixing formula:
     !   P = P - n/tau * dt * P(1-P^(1/n))
     do i=1,nleaf
        nH=max(uold(ind_leaf(i),1),smallr) ! rho
        if (nH == smallr) print *,"pristine_update: *** nH < smallr *** ",nH
        
        ! ***********************************************************
        ! Need to verify that using delta time new (dtnew) is correct
        ! ***********************************************************
        ! If the time step is 0, no point in doing this computation
        ! However, this should only happen at the start when Pgas
        ! is 1.0 ... hence it would be skipped anyway.
        if (dtnew(ilevel) == 0.0d0) cycle

        ! If the cell has non-zero PGF OR we want to compute/output v_t ... 
        if ((uold(ind_leaf(i),iprist) > 0.0d0) .OR. turbulent_velocity ) then
           
           ! Remove the density factor from the pristine fraction scalar
           ! Note that if the cell's density has been updated (due to
           ! SN or mass-flow from neighboring cells) Pgas might be > 1.0
           ! However, it shouldn't be > 1.0 since we handle a nH
           ! decrease in a cell by correcting iprist in feedback.f90.
           Pgas = uold(ind_leaf(i),iprist)/nH
           ! Check and correct numerical problems that may
           ! cause Pgas to be > 1.0 or < 0.0 
           if (Pgas > 1.0d0) then
              if (Pgas > 1.0001d0) then
                 print *, "Entering pristine_update with Pgas > 1.0, fixing: ", Pgas
              end if
              Pgas = 1.0d0
           else if (Pgas < 0.0d0) then ! Still need to check here since turb_vel may be .true.
              print *, "Entering pristine_update with NEGATIVE Pgas *** *** fixing: ", Pgas
              Pgas = 0.0d0
           endif

           ! ***********************************************************
           ! At this point the pristine fraction has been 'diluted' if
           ! the cell's density has increased since last time step.
           ! ***********************************************************
           
           ! Get the total fluid energy density (still have rho in here)
           totEnergy = uold(ind_leaf(i),ndim+2) ! An energy density (* rho)
           ! Subtract the bulk motion/kinetic energy to get internal energy 
           ! Just the random motions of the gas particles...
           kinEnergy = 0.0d0;
           do idim=1,ndim
              ! kinEnergy = 1/2 * rho * v^2
              ! Note that we're leaving a factor of rho in kinEnergy by squaring
              ! before dividing out nH. 
              kinEnergy = kinEnergy + 0.5d0 * (uold(ind_leaf(i),idim+1)**2)/nH
           end do
           thermEnergy = totEnergy - kinEnergy

           ! Compute temperature 
           ! P = nkT ( or PV = NkT -> P = rho kT w/ rho = n, number density)
           ! T = P/rho 1/k, with (gamma-1) * thermEnergy = P/rho
           ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
           temp = (gamma-1.0d0) * thermEnergy ! Pressure 
           temp = temp/nH * scale_T2 ! Now in terms of K/mu
           temp2 = uold(ind_leaf(i),ndim+2)/uold(ind_leaf(i),1)/uold(ind_leaf(i),1)

           ! divide by nH to get internal energy per unit vol -- thermal energy density
           thermEnergy = thermEnergy/nH

           ! Compute the sound speed
           ! Sqrt( gamma * P/rho ) with (gamma-1) * thermEnergy = P/rho
           Cs = sqrt( gamma*(gamma-1.0d0) * thermEnergy ) * scale_v;

           ! ***********************************************************
           ! ***********************************************************
           ! Compute the jacobian of the velocity across our cell
           ! This matrix describes how each vel component changes in
           ! each of the directions: dvx/dx, dvy/dx, dvz/dx, etc.
           ! ***********************************************************
           ! ***********************************************************
           ! First step, create a 2d array of the velocities for the
           ! nboring cells (and for this cell - cell "0").
           ! ind_father(i,0) is the index of ind_leaf(i) - this cell 
           ! ind_father(i,1->6) are the indices of the neigboring cells
           ! The second index uold(*,2->4) are the velocity components: x,y,z
           ! Remember that the uold velocity info has to be div'd by nH

           ! DEBUG
!!$           if (Pgas < 0.9999d0) then
!!$              print *, "Vels around cell ",ind_father(i,0)
!!$              write (*,2457) "x nbors ", ind_father(i,1), ind_father(i,2)
!!$              write (*,2457) "y nbors ", ind_father(i,3), ind_father(i,4)
!!$              write (*,2457) "z nbors ", ind_father(i,5), ind_father(i,6)
!!$           end if
!!$2457       format(A8,I7,",",I7)
           ! END DEBUG

           ! Build a 2d array of velocities for nboring cells.
           ! First index is cell (0->6), Second is vel component
           ! 1 => x, 2 => y, 3 => z for that cell
           ! Cell ordering map:
           ! vel(1,k) -> +x, vel(2,k) -> -x
           ! vel(3,k) -> +y, vel(4,k) -> -y
           ! vel(5,k) -> +z, vel(6,k) -> -z
           ! uold velocity, x, y, z are in 2nd indices 2, 3, 4
           ! which I copy here into indices 1,2,3 (xyz components
           ! of the vel in that cell)
           ! So, vel(1,1) = v_x for cell in +x direction
           !     vel(1,2) = v_y for cell in +x direction
           !     vel(2,1) = v_x for cell in -x direction 
           !     vel(2,2) = v_y for cell in -x direction 
           !     ...
           !     vel(5,3) = v_z for cell in +z direction
           ! vel(0,*) are velocities for "this" cell, the 
           ! one that is surrounded by the cells in x,y,z
           do j=0,twondim 
              do idim=1,ndim
                 vel(j,idim) = uold(ind_father(i,j),idim+1)/nH
              end do
              ! DEBUG Print *********
!!$              if (Pgas < 0.9999d0) then
!!$                 write (*,2456) j,vel(j,1),vel(j,2),vel(j,3) 
!!$2456             format("Vel(",I1,")=","(",ES12.5e2,",",ES12.5e2,",",ES12.5e2,")")
!!$              end if
              ! END DEBUG ***********
           end do

           ! Form the partials matrix - Aij
           ! Loop only to twondim by 2's since we use the left & right
           ! nbors each loop -- to compute the slope (derivative)
           indx = 1
           do j=1, twondim, 2
              do k=1, ndim
                 ! vel has had rho divided out (above)
                 ! dx_loc and vel are still in user units however
                 ! Interpolate if the nbor cells aren't the same
                 ! size ...
                 rightLevel = nbor_level(i,j);
                 leftLevel  = nbor_level(i,j+1);
                 if (rightLevel == leftLevel) then
                    ! Aij's are the difference of cellToRight - cellToLeft
                    ! for each coordiante direction (x,y,z)
                    Aij(indx,k) = (vel(j,k) - vel(j+1,k)) / (2.0d0 * dx_loc) 
                 else if (rightLevel > leftLevel) then
                    ! Interpolate by weighting the cell at the next level (to right)
                    ! up [vel(j,k)] by 2/3 and the current cell by 1/3
                    Aij(indx,k) = ((vel(j,k) * 2.0d0/3.0d0 + vel(0,k) * 1.0d0/3.0d0) - &
                         vel(j+1,k)) / (2.0d0 * dx_loc)
!!                    print *, "pristine_update - level Right > Left - Aij = ", Aij(indx,k) 
                 else ! left > right 
                    Aij(indx,k) = (vel(j,k) - &
                         (vel(j+1,k) * 2.0d0/3.0d0 + vel(0,k) * 1.0d0/3.0d0)) / (2.0d0 * dx_loc)
!!                    print *, "pristine_update - level Left > Right - Aij = ", Aij(indx,k) 
                    
                 end if
              end do
              indx = indx + 1
           end do

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           ! DEBUG
!!$           if ( Pgas < 0.9999d0) then
!!$              print *,"Aij for cell", ind_leaf(i)
!!$              do j=1, ndim
!!$                 write(*,1225) Aij(j,1),Aij(j,2),Aij(j,3)
!!$              end do
!!$           end if
!!$1225       format(ES12.5e2,",",ES12.5e2,",",ES12.5e2)

! -----------------------------------------------------------------------

           ! Compute 1/2 (Aij + Aij^T) - vel shear tensor
           ! So Sij = 0.5 (Aij + Aij^T)
           ! Also compute the eddy viscosity: ev = sqrt(2 Sij**2)
           ev = 0.0d0
           do j=1, ndim
              do k=1, ndim
                 Sij(j,k) = 0.5d0*(Aij(j,k) + Aij(k,j))
                 ! Save time -- also compute ev =  Sij * Sij
                 ! this is the square of the terms, summed
                 ev = ev + Sij(j,k)**2.0d0
              end do
           end do

           ! We need ev = sqrt(2 Sij Sij). We have Sij Sij from above... 
           ev = sqrt(2.0d0 * ev)

           turbVel = dx_loc * ev 
           
           !! ********* DEBUG *********
           !! turbVel = 1.0d7 / scale_v ! 100 km/s to cm/s
           !! ********* DEBUG *********

           ! Remember, when we save the v_t, it must be in code units
           if (turbulent_velocity) uold(ind_leaf(i),iturbvel) = turbVel * nH 

           turbVel = turbVel * scale_v ! Work in physical units: cm/s
           mach    = turbVel/Cs 

!!$ -----------------------------------------------------------------------
           ! Compute a threshold for P. 
           ! What we want to do is look at the metallicity of the polluted
           ! region: Z = <Z>/(1-P). If Z > Z_crit, we mix otherwise continued
           ! mixing really increases P since we've diluted Z below Z_crit
           ! Another way to write it:
           ! If P > 1 - <Z>/Z_crit we mix.
           ! If P < 1 - <Z>/Z_crit further mixing actually increases P.
           Zave      = uold(ind_leaf(i),imetal)/nH 
           threshold = 1.0d0 - Zave/Zcrit

           ! If we're pristine, nothing in the cell so
           ! exit this loop. 
           if (Pgas >= 1.0d0) then
              uold(ind_leaf(i),iprist) = 1.0d0 * nH
              cycle
           else if (Pgas < 1.0d-10) then
              uold(ind_leaf(i),iprist) = 0.0d0
           end if
           ! DEBUG **** 
           ! if ((Pgas < threshold - 1d-6) .and. (Pgas < 1.0d0)) then
           !    print *,"P < 1-<Z>/Z_crit -- Increase P"
           !    print *,"P=",Pgas
           !    print *,"1 - <Z>/Z_crit=",threshold
           ! end if
                         
           ! **************************************************
           ! **************************************************
           ! Update the pristine gas fraction           
           ! **************************************************
           ! **************************************************
!!           opf = Pgas ! Original Pristine Fraction - before update

           ! If <Z>/f_pol < Z_c futher mixing doesn't increase the polluted fraction
           ! <Z>/f_pol < Z_c => <Z>/Z_c < f_pol => <Z>/Z_c < 1-P => P < 1-<Z>/Z_c
           ! to evolve it back toward pristine. 
           if ((Pgas < 1.0d0) .and. (Pgas < threshold)) then
              ! P needs to increase
              ! Increase as an exponential of the eddy turnover time
              Pgas = 1.0d0 - (1.0d0 - Pgas) * exp(-dtnew(ilevel) * scale_t/ (dx_loc * scale_l / turbVel))
              uold(ind_leaf(i),iprist) = Pgas * nH
              cycle 
           end if

           n = 1.0d0 + 11.0d0 * exp(-mach/3.5d0)

           ! If we made it here compute an updated P
           x = -log10(10.0d0**(-7.0d0) * Zave/Zcrit) / log10(Zave/Zcrit)

           ! We need the log term to make tau smaller as Zave gets bigger
           ! i.e. - a shorter dynamical mixing time for bigger <Z>
           ! 7 Nov 2014 FIX
           if (Pgas >= 0.9d0)then
              tau = (0.225d0 - (0.055d0 * exp(-mach**(3.0d0/2.0d0) / 4.0d0))) * sqrt(x/5.0d0 + 1.0d0)
!!$              print *,"Raw tau (large Pgas): ",tau
           else ! 0 < Pgas < 0.9
              tau = (0.335d0 - (0.095d0 * exp(-mach**(2.0d0) / 4.0d0))) * sqrt(x/3.0d0 + 1.0d0)
!!$              print *,"Raw tau (smaller Pgas): ",tau
           end if
           tau = tau * dx_loc/turbVel * scale_l ! tau in sec. scale dx_loc, turbVel already scaled
!!$           print *,"Scaled tau: ", tau

           ! ------------------------------------
           ! Update the pristine gas fraction
           ! ------------------------------------
           tnorm = n * dtnew(ilevel) * scale_t / tau

           Pbefore = Pgas
           if (tnorm < 0.25d0) then
              Pgas = Pgas * (1.0d0 - tnorm * ( 1.0d0 - Pgas**(1.0d0/n)))
           else
              ! The time step is too large... 
              ! Evolve the pristine fraction over several smaller steps
              nsteps = int(tnorm/0.25) + 1
              dtnorm = tnorm/real(nsteps)

              tempDecay = 0.0d0
              do j=1,nsteps 
                 Pgas = Pgas * (1.0d0 - dtnorm * ( 1.0d0 - Pgas**(1.0d0/n)))
                 tempDecay = tempDecay + Pgas * tnorm * ( 1.0d0 - Pgas**(1.0d0/n))
              end do
!!$              print *, "pristine_update: decay amount: ", tempDecay
           end if

           ! DEBUGGING TEST
           if (Pgas > Pbefore) then
              print *, "pristine_update: ERROR gas got MORE pristine! *** *** *** ", Pbefore, Pgas
           end if

           if (Pgas < 0.0d0) then
              print *,"pristine_update: ERROR after update - Pgas < 0.0, fixing: ",Pgas
              Pgas = 0.0d0
           else if (Pgas > 1.01) then 
              print *,"pristine_update: ERROR gas got MORE pristine! *** ERROR *** ",Pgas
              Pgas = 1.0d0
           endif
           uold(ind_leaf(i),iprist) = Pgas * nH

        else ! Don't let the pristine fraction go below 0.0
           uold(ind_leaf(i),iprist) = 0.0d0;
        endif
        ! DEBUG
!!$        write (*,'("cell ", I7," PGF before ", F8.5," PFG after ", F8.5 )') ind_leaf(i), Pbefore, Pgas
        ! END DEBUG
     end do 
  end do
  ! End loop over cells

end subroutine pristine_up

!################################################################
!################################################################
subroutine printSij(indx,Sij,i,ind_father)
  use hydro_commons
  use prist_commons
  integer::indx,i ! indx is the index number of "this" cell, i is leaf index
  real(dp),dimension(1:ndim,1:ndim)::Sij
  integer,dimension(1:nvector,0:twondim)::ind_father
  
  write (*, '("Sij matrix for cell ",I10)') indx
!  write (*, '("Neighbors indices x ",I7," ",I7)') ind_father(i,1),ind_father(i,2)
!  write (*, '("Neighbors indices y ",I7," ",I7)') ind_father(i,3),ind_father(i,4)
!  write (*, '("Neighbors indices z ",I7," ",I7)') ind_father(i,5),ind_father(i,6)
  write (*,2001) Sij(1,1),Sij(1,2),Sij(1,3)
  write (*,2001) Sij(2,1),Sij(2,2),Sij(2,3)
  write (*,2001) Sij(3,1),Sij(3,2),Sij(3,3)
2001 format ( (ES12.5E2,' ',ES12.5E2,' ',ES12.5E2))
end subroutine printSij
!################################################################
!################################################################
subroutine printNeighbors(igridn,nbor_level,ngrid,ilevel)
  use hydro_commons
  use prist_commons
  integer::ngrid,ilevel
  integer,dimension(1:nvector,0:twondim)::igridn,nbor_level
  integer::i,j,nsmall
  !---------------------------------------------------------
  ! Print out the neigboring cell indices for the cells
  ! of interest.
  !---------------------------------------------------------
  nsmall = min(ngrid,5) ! Don't print out too many cells from igridn

  do i=1,nsmall
     do j=1,twondim
        if (nbor_level(i,j) /= ilevel ) then
           print *,"ilevel = ", ilevel
           print *,"igridn(",i,") = ",igridn(i,0)
           print *,"should also be eq to ind_leaf(",i,") ->",ind_leaf(i)
           print *,"igridn's nbor ", igridn(i,j)
           print *,"igridn's level ",nbor_level(i,j) 
        end if
     end do
  end do
  
end subroutine printNeighbors


!################################################################
!################################################################
subroutine pristine_output(i,ilevel)
  use hydro_commons
  use prist_commons
  integer::i
  integer::ilevel
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  ! Conversion factor from user units to cgs units
!!$  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  print *, "----- Pristine Fraction Update ----- "
  print *, "cell = ", i
  print *, "Pgas start=", opf
  print *, "Pgas decayed=", Pgas
  print *, "Pgas * rho=", uold(ind_leaf(i),iprist)
  ! Print the time in years... 
  print *,"scalet= ", scale_t
  print *,"scalel= ", scale_l
  print *, "scaled=", scale_d
  print *, "scalev=", scale_v
  print *, "scalenH= ", scale_nH
  print *, "rho (_)= ", nH
  print *, "rho (g/cc)= ", nH * scale_d
  print *, "nH (H/cc)= ", nH * scale_nH
  print *, "dtold (_)= ", dtold(ilevel)
  print *, "dtnew (_)= ", dtnew(ilevel)
  print *, "time (_)= ", t 
  print *, "dtold (sec)= ", dtold(ilevel)*scale_t
  print *, "dtnew (sec)= ", dtnew(ilevel)*scale_t
  print *, "dtold (Myr)= ", dtold(ilevel)*scale_t/(365.0 * 24.0 * 3600.0)
  print *, "dtnew (Myr)= ", dtnew(ilevel)*scale_t/(365.0 * 24.0 * 3600.0)
!!$  print *, "aexp_frw= ", aexp_frw
  print *, "aexp_old= ", aexp_old
  print *, "time (yr)= ", t * scale_t/(365.0 * 24.0 * 3600.0)
  call getProperTime(t,tproper)
  print *, "proper time (Myr)", - tproper * scale_t/(365.0 * 24.0 * 3600.0)/(10.0**6)
  print *, "gamma ", gamma
  print *, "gas' Cs (cm/s) = ", Cs
  print *, "turb vel (cm/s)= ", turbVel
  print *, "turb M (#) = ", mach
  print *, "tau (s)= ", tau
  write(*,1215) "n= ",n
  print *, "Total fluid: E * rho (_)= ", uold(ind_leaf(i),ndim+2)
  print *, "E (_)= ", totEnergy/nH
  print *, "E (T/mu)= ", totEnergy/nH * scale_T2
  print *, "thermEnergy (_)= ", thermEnergy
  print *, "thermEnergy (T/mu)= ", thermEnergy * scale_T2
  print *, "thermEnergy * rho (_)= ", thermEnergy * nH
  print *, "kinEnergy (_)= ", kinEnergy/nH
  print *, "kinEnergy (T/mu)= ", kinEnergy/nH * scale_T2
  print *, "kinEnergy * rho (_)= ", kinEnergy
  ! Compute T2=T/mu in Kelvin
  print *,"T/mu (K)= ", temp
  print *,"T (assuming mu = 2)(K)= ", temp * 2
  print *, "vel1 (cm/s)= ", uold(ind_leaf(i),2)/nH * scale_v
  print *, "vel2 (cm/s)= ", uold(ind_leaf(i),3)/nH * scale_v
  print *, "vel3 (cm/s)= ", uold(ind_leaf(i),4)/nH * scale_v
  print *, "dx ()= ", dx 
  print *, "dx (pc)= ", dx * scale_l / (3.086E18)
  print *
1215 format(A5, ES12.5e2)
end subroutine pristine_output

