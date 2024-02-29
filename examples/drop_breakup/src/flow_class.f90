!> Various definitions and tools for running an NGA2 simulation
module flow_class
   use string,            only: str_medium
   use precision,         only: WP
   use config_class,      only: config
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use ccl_class,         only: ccl
   use lpt_class,         only: lpt
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private

   public :: flow

   !> Flow object
   type :: flow

      !> Config
      type(config) :: cfg
   
      !> Two-phase incompressible flow solver, VF solver with CCL, and corresponding time tracker and sgs model
      type(tpns),        public :: fs
      type(vfs),         public :: vf
      type(timetracker), public :: time
      type(sgsmodel),    public :: sgs
      type(ccl),         public :: cc
      type(lpt),         public :: lp
      type(hypre_str),   public :: ps
      type(ddadi),       public :: vs
      
      !> Provide two datafiles and an event tracker for saving restarts
      type(event)    :: save_evt
      type(datafile) :: df
      character(len=str_medium) :: irl_file
      character(len=str_medium) :: lpt_file
      logical :: restarted
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh
      type(partmesh) :: pmesh
      type(ensight)  :: ens_out
      type(event)    :: ens_evt
      
      !> Simulation monitor file
      type(monitor) :: mfile,cflfile,sprayfile
      
      !> Private work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
      real(WP), dimension(:,:,:,:), allocatable :: SR
      
      !> Problem definition
      real(WP), dimension(3) :: center,radii
      
      !> Transfer model parameters
      real(WP) :: filmthickness_over_dx  =5.0e-1_WP
      real(WP) :: min_filmthickness      =1.0e-3_WP
      real(WP) :: diam_over_filmthickness=1.0e+1_WP
      real(WP) :: max_eccentricity       =5.0e-1_WP
      real(WP) :: d_threshold            =1.0e-1_WP
      
      !> SGS surface tension model
      real(WP), dimension(:,:,:), allocatable :: sgsSTx,sgsSTy,sgsSTz


   contains
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
   end type flow


contains
   
   
   !> Function that localizes the left (x-) of the domain
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator
   
   
   !> Function that localizes the right (x+) of the domain
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   
   
   !> Initialization of problem solver
   subroutine init(this)
      use param, only: param_read
      implicit none
      class(flow), intent(inout) :: this
      
      
      ! Handle restart/saves here
      restart_and_save: block
         character(len=str_medium) :: dir_restart
         ! CAREFUL - WE NEED TO CREATE THE TIMETRACKER BEFORE THE EVENT !
         this%time=timetracker(this%cfg%amRoot,name='drop break-up')
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call param_read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call param_read(tag='Restart from',val=dir_restart,short='r',default='')
         this%restarted=.false.; if (len_trim(dir_restart).gt.0) this%restarted=.true.
         if (this%restarted) then
            ! If we are, read the name of the directory
            call param_read('Restart from',dir_restart,'r')
            ! Read the datafile and the name of the IRL file to read later
            this%df=datafile(pg=this%cfg,fdata=trim(adjustl(dir_restart))//'/'//'data')
            this%irl_file=trim(adjustl(dir_restart))//'/'//'data.irl'
            this%lpt_file=trim(adjustl(dir_restart))//'/'//'data.lpt'
         else
            ! If we are not restarting, we will still need datafiles for saving restart files
            this%df=datafile(pg=this%cfg,filename=trim(this%cfg%name),nval=2,nvar=10)
            this%df%valname(1)='t'
            this%df%valname(2)='dt'
            this%df%varname(1)='U'
            this%df%varname(2)='V'
            this%df%varname(3)='W'
            this%df%varname(4)='P'
            this%df%varname(5)='Pjx'
            this%df%varname(6)='Pjy'
            this%df%varname(7)='Pjz'
            this%df%varname(8)='LM'
            this%df%varname(9)='MM'
            this%df%varname(10)='VF'
         end if
      end block restart_and_save
      
      ! Create the flow mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         ! Read in grid definition
         call param_read('Domain Lx',Lx); call param_read('Domain nx',nx); allocate(x(nx+1))
         call param_read('Domain Ly',Ly); call param_read('Domain ny',ny); allocate(y(ny+1))
         call param_read('Domain Lz',Lz); call param_read('Domain nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='Flow')
         ! Read in partition
         call param_read('Whole Domain Partition',partition,short='p')
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config

      ! Allocate work arrays for cfg
      allocate_work_arrays: block
         allocate(this%resU  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%SR  (6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%sgsSTx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%sgsSTx=0.0_WP
         allocate(this%sgsSTy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%sgsSTy=0.0_WP
         allocate(this%sgsSTz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%sgsSTz=0.0_WP
      end block allocate_work_arrays
      
      
      ! Initialize time tracker
      initialize_timetracker: block
         !time=timetracker(cfg%amRoot,name='drop break-up')   !< This is moved up for restarts!
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         call param_read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
         ! Handle restart
         if (this%restarted) then
            call this%df%pullval(name='t' ,val=this%time%t )
            call this%df%pullval(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
         end if
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         ! use vfs_class, only: VFlo,VFhi,lvira,r2p,art,swartz
         use vfs_class, only: VFlo,VFhi,lvira,r2p,swartz
         integer :: i,j,k,si,sj,sk,n
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         ! integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver with LVIRA
         this%vf=vfs(cfg=this%cfg,reconstruction_method=lvira,name='VOF')
         ! Create a VOF solver with R2P
         !vf=vfs(cfg=cfg,reconstruction_method=r2p,name='VOF')
         !vf%VFflot =1.0e-4_WP !< Enables flotsam removal
         !vf%VFsheet=1.0e-2_WP !< Enables sheet removal
         ! Create a VOF solver with ART
         !vf=vfs(cfg=cfg,reconstruction_method=art,name='VOF')
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[this%vf%cfg%x(i+si),this%vf%cfg%y(j+sj),this%vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  ! call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_drop,0.0_WP,amr_ref_lvl)
                  this%vf%VF(i,j,k)=vol/this%vf%cfg%vol(i,j,k)
                  if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
                     this%vf%Lbary(:,i,j,k)=v_cent
                     this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
                  else
                     this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Handle restart - using IRL data - tested but not working
         !if (restarted) then
         !   ! Get the IRL interface
         !   call vf%read_interface(filename=trim(irl_file))
         !   ! Reset moments to guarantee compatibility with interface reconstruction
         !   call vf%reset_volume_moments()
         !end if
         ! Handle restart - using VF data
         if (this%restarted) call this%df%pullvar(name='VF',var=this%vf%VF)
         ! Update the band
         call this%vf%update_band()
         ! Perform interface reconstruction from VOF field
         call this%vf%build_interface()
         ! Set interface planes at the boundaries
         call this%vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call this%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call this%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call this%vf%subcell_vol()
         ! Calculate curvature
         call this%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call this%vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver with bconds
      create_solver: block
         use tpns_class, only: dirichlet,clipped_neumann,neumann
         ! use ils_class,  only: pcg_amg,gmres,gmres_amg
         use hypre_str_class
         use mathtools,       only: Pi
         ! Create a two-phase flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',this%fs%visc_l)
         call param_read('Gas dynamic viscosity'   ,this%fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',this%fs%rho_l)
         call param_read('Gas density'   ,this%fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',this%fs%sigma)
         ! Define boundary conditions: gas inflow on the left and outflow on the right
         call this%fs%add_bcond(name='inflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=xp_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=10
         ! call param_read('Pressure iteration',fs%psolv%maxit)
         ! call param_read('Pressure tolerance',fs%psolv%rcvg)
         call param_read('Pressure iteration',this%ps%maxit)
         call param_read('Pressure tolerance',this%ps%rcvg)
         ! Configure implicit velocity solver
         ! call param_read('Implicit iteration',fs%implicit%maxit)
         ! call param_read('Implicit tolerance',fs%implicit%rcvg)
         this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
         ! Calculate cell-centered velocities and divergence
		   call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
		   call this%fs%get_div()
         ! call fs%setup(pressure_solver=gmres_amg,implicit_solver=gmres_amg)
      end block create_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools,  only: pi
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         real(WP) :: Uin
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Handle restart
         if (this%restarted) then
            call this%df%pullvar(name='U'  ,var=this%fs%U  )
            call this%df%pullvar(name='V'  ,var=this%fs%V  )
            call this%df%pullvar(name='W'  ,var=this%fs%W  )
            call this%df%pullvar(name='P'  ,var=this%fs%P  )
            call this%df%pullvar(name='Pjx',var=this%fs%Pjx)
            call this%df%pullvar(name='Pjy',var=this%fs%Pjy)
            call this%df%pullvar(name='Pjz',var=this%fs%Pjz)
         end if
         ! Apply Dirichlet at inflow
         call param_read('Gas velocity',Uin)
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=Uin
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Compute MFR through all boundary conditions
         call this%fs%get_mfr()
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
      end block initialize_velocity
      
      
      ! Create a connected-component labeling object
      create_and_initialize_ccl: block
         use vfs_class, only: VFlo
         ! Create the CCL object
         this%cc=ccl(cfg=this%cfg,name='CCL')
         this%cc%max_interface_planes=2
         this%cc%VFlo=VFlo
         this%cc%dot_threshold=-0.5_WP
         this%cc%thickness_cutoff=this%filmthickness_over_dx
         ! Perform CCL step
         call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
         call this%cc%film_classify(Lbary=this%vf%Lbary,Gbary=this%vf%Gbary)
         call this%cc%deallocate_lists()
      end block create_and_initialize_ccl
      
      
      ! Create a Lagrangian spray tracker
      create_lpt: block
         ! Create the solver
         this%lp=lpt(cfg=this%cfg,name='spray')
         ! Get particle density from the flow solver
         this%lp%rho=this%fs%rho_l
         ! Handle restarts
         if (this%restarted) call this%lp%read(filename=trim(this%lpt_file))
      end block create_lpt
      
      
      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
         ! Handle restart
         if (this%restarted) then
            call this%df%pullvar(name='LM',var=this%sgs%LM)
            call this%df%pullvar(name='MM',var=this%sgs%MM)
         end if
      end block create_sgs
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         this%smesh=surfmesh(nvar=1,name='plic')
         this%smesh%varname(1)='nplane'
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh(this%smesh)
         ! Also populate nplane variable
         this%smesh%var(1,:)=1.0_WP
         np=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         ! Include an extra variable for droplet diameter
         this%pmesh=partmesh(nvar=1,nvec=0,name='lpt')
         this%pmesh%varname(1)='diameter'
         ! Transfer particles to pmesh
         call this%lp%update_partmesh(this%pmesh)
         ! Also populate diameter variable
         do i=1,this%lp%np_
            this%pmesh%var(1,i)=this%lp%p(i)%d
         end do
      end block create_pmesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(this%cfg,'drop_breakup')
         ! Create event for Ensight output
         this%ens_evt=event(this%time,'Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curv',this%vf%curv)
         call this%ens_out%add_vector('sgsST',this%sgsSTx,this%sgsSTy,this%sgsSTz)
         call this%ens_out%add_scalar('filmThickness',this%cc%film_thickness)
         call this%ens_out%add_surface('vofplic',this%smesh)
         call this%ens_out%add_particle('spray',this%pmesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         call this%lp%get_max()
         ! Create simulation monitor
         this%mfile=monitor(amroot=this%fs%cfg%amRoot,name='simulation')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFmax,'VOF maximum')
         call this%mfile%add_column(this%vf%VFmin,'VOF minimum')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(amroot=this%fs%cfg%amRoot,name='cfl')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%add_column(this%fs%CFLst,'ST CFL')
         call this%cflfile%write()
         ! Create a spray monitor
         this%sprayfile=monitor(amroot=this%lp%cfg%amRoot,name='spray')
         call this%sprayfile%add_column(this%time%n,'Timestep number')
         call this%sprayfile%add_column(this%time%t,'Time')
         call this%sprayfile%add_column(this%time%dt,'Timestep size')
         call this%sprayfile%add_column(this%lp%np,'Droplet number')
         call this%sprayfile%add_column(this%lp%Umin, 'Umin')
         call this%sprayfile%add_column(this%lp%Umax, 'Umax')
         call this%sprayfile%add_column(this%lp%Umean,'Umean')
         call this%sprayfile%add_column(this%lp%Vmin, 'Vmin')
         call this%sprayfile%add_column(this%lp%Vmax, 'Vmax')
         call this%sprayfile%add_column(this%lp%Vmean,'Vmean')
         call this%sprayfile%add_column(this%lp%Wmin, 'Wmin')
         call this%sprayfile%add_column(this%lp%Wmax, 'Wmax')
         call this%sprayfile%add_column(this%lp%Wmean,'Wmean')
         call this%sprayfile%add_column(this%lp%dmin, 'dmin')
         call this%sprayfile%add_column(this%lp%dmax, 'dmax')
         call this%sprayfile%add_column(this%lp%dmean,'dmean')
         call this%sprayfile%write()
      end block create_monitor
      
      
   end subroutine init
   
   
   !> Perform an NGA2 simulation
   subroutine step(this)
      implicit none
      class(flow), intent(inout) :: this
      
      ! Perform time integration - the second solver is the main driver here
      do while (.not.this%time%done())
         
         ! Increment time
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%time%adjust_dt()
         call this%time%increment()
         
         ! Advance our spray
         this%resU=this%fs%rho_g; this%resV=this%fs%visc_g
         call this%lp%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU,visc=this%resV)
         
         ! Remember old VOF
         this%vf%VFold=this%vf%VF
         
         ! Remember old velocity
         this%fs%Uold=this%fs%U
         this%fs%Vold=this%fs%V
         this%fs%Wold=this%fs%W
         
         ! Prepare old staggered density (at n)
         call this%fs%get_olddensity(vf=this%vf)
         
         ! VOF solver step
         call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call this%fs%get_viscosity(vf=this%vf)
         
         ! Turbulence modeling - only work with gas properties here
         sgsmodel: block
            integer :: i,j,k
            ! call fs%get_strainrate(Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
            call this%fs%get_strainrate(SR=this%SR)
            this%resU=this%fs%rho_g
            call this%sgs%get_visc(type=2,dt=this%time%dtold,rho=this%resU,Ui=this%Ui,Vi=this%Vi,Wi=this%Wi,SR=this%SR)
            where (this%sgs%visc.lt.-this%fs%visc_g)
               this%sgs%visc=-this%fs%visc_g
            end where
            do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_
               do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_
                  do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_
                     this%fs%visc(i,j,k)   =this%fs%visc(i,j,k)   +this%sgs%visc(i,j,k)
                     this%fs%visc_xy(i,j,k)=this%fs%visc_xy(i,j,k)+sum(this%fs%itp_xy(:,:,i,j,k)*this%sgs%visc(i-1:i,j-1:j,k))
                     this%fs%visc_yz(i,j,k)=this%fs%visc_yz(i,j,k)+sum(this%fs%itp_yz(:,:,i,j,k)*this%sgs%visc(i,j-1:j,k-1:k))
                     this%fs%visc_zx(i,j,k)=this%fs%visc_zx(i,j,k)+sum(this%fs%itp_xz(:,:,i,j,k)*this%sgs%visc(i-1:i,j,k-1:k))
                  end do
               end do
            end do
         end block sgsmodel
         
         ! Perform sub-iterations
         do while (this%time%it.le.this%time%itmax)
            
            ! Build mid-time velocity
            this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
            this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
            this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call this%fs%prepare_advection_upwind(dt=this%time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
            
            ! ! Add sgs ST model
            ! STmodel_add: block
            !    integer :: i,j,k
            !    do k=vf%cfg%kmin_,vf%cfg%kmax_
            !       do j=vf%cfg%jmin_,vf%cfg%jmax_
            !          do i=vf%cfg%imin_,vf%cfg%imax_
            !             if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+sum(fs%itpi_x(:,i,j,k)*sgsSTx(i-1:i,j,k))
            !             if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+sum(fs%itpi_y(:,i,j,k)*sgsSTy(i,j-1:j,k))
            !             if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+sum(fs%itpi_z(:,i,j,k)*sgsSTz(i,j,k-1:k))
            !          end do
            !       end do
            !    end do
            ! end block STmodel_add
            
            ! Assemble explicit residual
            this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
            this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
            this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW
            
            ! Form implicit residuals
            call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
            
            ! Apply these residuals
            this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
            this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
            this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
            
            ! Apply other boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
            
            ! Solve Poisson equation
            call this%fs%update_laplacian()
            call this%fs%correct_mfr()
            call this%fs%get_div()
            call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
            this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
            this%fs%psolv%sol=0.0_WP
            call this%fs%psolv%solve()
            call this%fs%shift_p(this%fs%psolv%sol)
            
            ! Correct velocity
            call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
            this%fs%P=this%fs%P+this%fs%psolv%sol
            this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
            this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
            this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
            
            ! Increment sub-iteration counter
            this%time%it=this%time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         call this%fs%get_div()
         
         ! Perform volume-fraction-to-droplet transfer
         call transfer_vf_to_drops(this)
         
         ! Output to ensight
         if (this%ens_evt%occurs()) then
            ! Update surfmesh object
            update_smesh: block
               use irl_fortran_interface
               integer :: i,j,k,nplane,np
               ! Transfer polygons to smesh
               call this%vf%update_surfmesh(this%smesh)
               ! Also populate nplane variable
               this%smesh%var(1,:)=1.0_WP
               np=0
               do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
                  do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                     do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh
            ! Update partmesh object
            update_pmesh: block
               integer :: i
               ! Transfer particles to pmesh
               call this%lp%update_partmesh(this%pmesh)
               ! Also populate diameter variable
               do i=1,this%lp%np_
                  this%pmesh%var(1,i)=this%lp%p(i)%d
               end do
            end block update_pmesh
            ! Perform ensight output
            call this%ens_out%write_data(this%time%t)
         end if
         
         ! Perform and output monitoring
         call this%lp%get_max()
         call this%fs%get_max()
         call this%vf%get_max()
         call this%mfile%write()
         call this%cflfile%write()
         call this%sprayfile%write()
         
         ! After we're done clip all VOF at the exit area
         vf_side_clipping: block
            integer :: i,j,k
            do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_
               do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_
                  do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
                     if (i.ge.this%vf%cfg%imax-5) this%vf%VF(i,j,k)=0.0_WP
                  end do
               end do
            end do
         end block vf_side_clipping
         
         ! Finally, see if it's time to save restart files
         if (this%save_evt%occurs()) then
            save_restart: block
               character(len=str_medium) :: dirname,timestamp
               ! Prefix for files
               dirname='restart_'; write(timestamp,'(es12.5)') this%time%t
               ! Prepare a new directory
               if (this%fs%cfg%amRoot) call execute_command_line('mkdir -p '//trim(adjustl(dirname))//trim(adjustl(timestamp)))
               ! Populate df and write it
               call this%df%pushval(name=  't',val=this%time%t )
               call this%df%pushval(name= 'dt',val=this%time%dt)
               call this%df%pushvar(name=  'U',var=this%fs%U   )
               call this%df%pushvar(name=  'V',var=this%fs%V   )
               call this%df%pushvar(name=  'W',var=this%fs%W   )
               call this%df%pushvar(name=  'P',var=this%fs%P   )
               call this%df%pushvar(name='Pjx',var=this%fs%Pjx )
               call this%df%pushvar(name='Pjy',var=this%fs%Pjy )
               call this%df%pushvar(name='Pjz',var=this%fs%Pjz )
               call this%df%pushvar(name= 'LM',var=this%sgs%LM )
               call this%df%pushvar(name= 'MM',var=this%sgs%MM )
               call this%df%pushvar(name= 'VF',var=this%vf%VF  )
               call this%df%write(fdata=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data')
               ! Also output IRL interface
               call this%vf%write_interface(filename=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.irl')
               ! Also output particles
               call this%lp%write(filename=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.lpt')
            end block save_restart
         end if
         
      end do
      
   end subroutine step
   
   
   !> Transfer vf to drops
   subroutine transfer_vf_to_drops(this)
      implicit none
      class(flow) :: this
      
      ! Perform a first pass with simplest CCL
      call this%cc%build_lists(VF=this%vf%VF,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      
      ! Loop through identified detached structs and remove those that are spherical enough
      remove_struct: block
         use mathtools, only: pi
         integer :: m,n,l,i,j,k,np
         real(WP) :: lmin,lmax,eccentricity,diam
         
         ! Loops over film segments contained locally
         do m=1,this%cc%n_meta_struct
            
            ! Test if sphericity is compatible with transfer
            lmin=this%cc%meta_structures_list(m)%lengths(3)
            if (lmin.eq.0.0_WP) lmin=this%cc%meta_structures_list(m)%lengths(2) ! Handle 2D case
            lmax=this%cc%meta_structures_list(m)%lengths(1)
            eccentricity=sqrt(1.0_WP-lmin**2/lmax**2)
            if (eccentricity.gt.this%max_eccentricity) cycle
            
            ! Test if diameter is compatible with transfer
            diam=(6.0_WP*this%cc%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
            if (diam.eq.0.0_WP.or.diam.gt.this%d_threshold) cycle
            
            ! Create drop from available liquid volume - only one root does that
            if (this%cc%cfg%amRoot) then
               ! Make room for new drop
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(0,8)                                                                                 !< Give id (maybe based on break-up model?)
               this%lp%p(np)%dt  =0.0_WP                                                                                   !< Let the drop find it own integration time
               this%lp%p(np)%Acol=0.0_WP                                                                                   !< Give zero collision force (axial)
               this%lp%p(np)%Tcol=0.0_WP                                                                                   !< Give zero collision force (tangential)
               this%lp%p(np)%d   =diam                                                                                     !< Assign diameter to account for full volume
               this%lp%p(np)%pos =[this%cc%meta_structures_list(m)%x,this%cc%meta_structures_list(m)%y,this%cc%meta_structures_list(m)%z] !< Place the drop at the liquid barycenter
               this%lp%p(np)%vel =[this%cc%meta_structures_list(m)%u,this%cc%meta_structures_list(m)%v,this%cc%meta_structures_list(m)%w] !< Assign mean structure velocity as drop velocity
               this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])                !< Place the drop in the proper cell for the lp%cfg
               this%lp%p(np)%flag=0                                                                                        !< Activate it
               ! Increment particle counter
               this%lp%np_=np
            end if
            
            ! Find local structs with matching id
            do n=this%cc%sync_offset+1,this%cc%sync_offset+this%cc%n_struct
               if (this%cc%struct_list(this%cc%struct_map_(n))%parent.ne.this%cc%meta_structures_list(m)%id) cycle
               ! Remove liquid in meta-structure cells
               do l=1,this%cc%struct_list(this%cc%struct_map_(n))%nnode ! Loops over cells within local
                  i=this%cc%struct_list(this%cc%struct_map_(n))%node(1,l)
                  j=this%cc%struct_list(this%cc%struct_map_(n))%node(2,l)
                  k=this%cc%struct_list(this%cc%struct_map_(n))%node(3,l)
                  ! Remove liquid in that cell
                  this%vf%VF(i,j,k)=0.0_WP
               end do
            end do
            
         end do
         
      end block remove_struct
      
      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      
      ! Clean up CCL
      call this%cc%deallocate_lists()
      
      ! Perform more detailed CCL in a second pass
      this%cc%max_interface_planes=2
      call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%cc%get_min_thickness()
      call this%cc%sort_by_thickness()
      
      ! Loop through identified films and remove those that are thin enough
      remove_film: block
         use mathtools, only: pi
         integer :: m,n,i,j,k,np,ip,np_old
         real(WP) :: Vt,Vl,Hl,Vd
         
         ! Loops over film segments contained locally
         do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
            
            ! Skip non-liquid films
            if (this%cc%film_list(this%cc%film_map_(m))%phase.ne.1) cycle
            
            ! Skip films that are still thick enough
            if (this%cc%film_list(this%cc%film_map_(m))%min_thickness.gt.this%min_filmthickness) cycle
            
            ! We are still here: transfer the film to drops
            Vt=0.0_WP      ! Transferred volume
            Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
            np_old=this%lp%np_  ! Remember old number of particles
            do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
               i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
               j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
               k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
               ! Increment liquid volume to remove
               Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
               ! Estimate drop size based on local film thickness in current cell
               Hl=max(this%cc%film_thickness(i,j,k),this%min_filmthickness)
               Vd=pi/6.0_WP*(this%diam_over_filmthickness*Hl)**3
               ! Create drops from available liquid volume
               do while (Vl-Vd.gt.0.0_WP)
                  ! Make room for new drop
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
                  this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                  this%lp%p(np)%acol =0.0_WP                                     !< Give zero collision force
                  this%lp%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
                  this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                          !< Place the drop at the liquid barycenter
                  this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
                  this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
                  this%lp%p(np)%flag=0                                          !< Activate it
                  ! Increment particle counter
                  this%lp%np_=np
                  ! Update tracked volumes
                  Vl=Vl-Vd
                  Vt=Vt+Vd
               end do
               ! Remove liquid in that cell
               this%vf%VF(i,j,k)=0.0_WP
            end do
            
            ! Based on how many particles were created, decide what to do with left-over volume
            if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
               ! Add one last drop for remaining liquid volume
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
               this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
               this%lp%p(np)%acol =0.0_WP                                     !< Give zero collision force
               this%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
               this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                          !< Place the drop at the liquid barycenter
               this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
               this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the lp%cfg
               this%lp%p(np)%flag=0                                          !< Activate it
               ! Increment particle counter
               this%lp%np_=np
            else ! Some particles were created, make them all larger
               do ip=np_old+1,this%lp%np_
                  this%lp%p(ip)%d=this%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
               end do
            end if
         end do
         
      end block remove_film
      
      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      
      ! Clean up CCL
      call this%cc%deallocate_lists()
      
      ! Resync the spray
      call this%lp%sync()
      
   end subroutine transfer_vf_to_drops
   
   
   !> Finalize the NGA2 simulation
   subroutine final(this)
      implicit none
      class(flow), intent(inout) :: this
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%SR,this%sgsSTx,this%sgsSTy,this%sgsSTz)
      
   end subroutine final
   
   
end module flow_class