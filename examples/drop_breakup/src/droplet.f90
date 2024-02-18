!> Various definitions and tools for initializing NGA2 config
!> Modified from NGA2/examples/coupler_tester/src/geometry.f90 for using overset mesh in droplet simulations
module geometry
   use mpi_f08,      only: MPI_Group
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   public :: geometry_init
   public :: cfg_g,grp_g,isInGrp_g
   public :: cfg_l,grp_l,isInGrp_l
   
   !> Two groups and partitions, along with logicals
   integer, dimension(3) :: partition_g,partition_l
   logical :: isInGrp_g,isInGrp_l
   type(MPI_Group) :: grp_g,grp_l
   
   !> These are the two configs
   type(config) :: cfg_g,cfg_l

   ! !> Single config files
   ! type(config), public :: cfg !< This is the atomizing flow config
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      use parallel, only: comm,group,nproc,rank
      use mpi_f08,  only: MPI_Group,MPI_Group_range_incl
      implicit none
      integer, dimension(3,1) :: grange
      integer :: ierr
      ! type(sgrid) :: grid
      
      ! We start by reading in the two partitions
      call param_read('Gas Phase Partition',partition_g)
      call param_read('Liquid Phase Partition',partition_l)

      ! Create an MPI group along with logical for the first grid on the lowest ranks
      grange(:,1)=[0,product(partition_g)-1,1]
      call MPI_Group_range_incl(group,1,grange,grp_g,ierr)
      isInGrp_g=.false.; if (rank.le.product(partition_g)-1) isInGrp_g=.true.
      
      ! Create an MPI group along with logical for the second grid on the highest ranks
      grange(:,1)=[nproc-product(partition_l),nproc-1,1]
      call MPI_Group_range_incl(group,1,grange,grp_l,ierr)
      isInGrp_l=.false.; if (rank.ge.nproc-product(partition_l)) isInGrp_l=.true.
      
      ! Create gas phase grid from input params
      if (isInGrp_g) then
         create_grid_g: block
            use sgrid_class, only: cartesian,sgrid
            type(sgrid) :: grid
            integer :: i,j,k,nx,ny,nz
            real(WP) :: Lx,Ly,Lz
            real(WP), dimension(:), allocatable :: x,y,z
            ! Read in grid definition
            call param_read('Gas Phase Lx',Lx); call param_read('Gas Phase nx',nx); allocate(x(nx+1))
            call param_read('Gas Phase Ly',Ly); call param_read('Gas Phase ny',ny); allocate(y(ny+1))
            call param_read('Gas Phase Lz',Lz); call param_read('Gas Phase nz',nz); allocate(z(nz+1))
            ! Create simple rectilinear grid
            do i=1,nx+1
               x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.25_WP*Lx
            end do
            do j=1,ny+1
               y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
            end do
            do k=1,nz+1
               z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
            end do
            ! General serial grid object with overlap=3 for tpns
            grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='gas_phase')
            ! Use it to create a config
            cfg_g=config(grp=grp_g,decomp=partition_g,grid=grid)
            ! Do not place walls
            cfg_g%VF=1.0_WP
         end block create_grid_g
      end if
      
      ! Create liquid phase grid from input params
      if (isInGrp_l) then
         create_grid_l: block
            use sgrid_class, only: cartesian,sgrid
            type(sgrid) :: grid
            integer :: i,j,k,nx,ny,nz
            real(WP) :: Lx,Ly,Lz
            real(WP), dimension(:), allocatable :: x,y,z
            real(WP), dimension(3) :: center,radii
            ! Read in grid definition
            ! call param_read('Liquid Phase Lx',Lx); 
            call param_read('Liquid Phase nx',nx); allocate(x(nx+1))
            ! call param_read('Liquid Phase Ly',Ly); 
            call param_read('Liquid Phase ny',ny); allocate(y(ny+1))
            ! call param_read('Liquid Phase Lz',Lz); 
            call param_read('Liquid Phase nz',nz); allocate(z(nz+1))
            ! Create simple rectilinear grid
            call param_read('Droplet center',center)
            call param_read('Droplet radii',radii)
            Lx=2.4*radii(1)
            Ly=2.4*radii(2)
            Lz=2.4*radii(3)
            do i=1,nx+1
               x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx+center(1)
            end do
            do j=1,ny+1
               y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly+center(2)
            end do
            do k=1,nz+1
               z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz+center(3)
            end do
            ! General serial grid object with overlap=3 for tpns
            grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='liquid_phase')
            ! Use it to create a config
            cfg_l=config(grp=grp_l,decomp=partition_l,grid=grid)
            ! Do not place walls
            cfg_l%VF=1.0_WP   
         end block create_grid_l
      end if
      
      
      ! ! Create config from grid
      ! create_cfg: block
      !    use parallel, only: group
      !    integer, dimension(3) :: partition
      !    ! Read in partition
      !    call param_read('Partition',partition,short='p')
      !    ! Create partitioned grid
      !    cfg=config(grp=group,decomp=partition,grid=grid)
      !    ! No walls
      !    cfg%VF=1.0_WP
      ! end block create_cfg
      
      
   end subroutine geometry_init
   
   
end module geometry
