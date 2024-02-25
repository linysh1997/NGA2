!> Various definitions and tools for initializing NGA2 config
!> Modified from NGA2/examples/coupler_tester/src/geometry.f90 for using overset mesh in droplet simulations
module geometry
   use mpi_f08,      only: MPI_Group
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   public :: geometry_init
   public :: cfg,grp,isInGrp
   
   !> Two groups and partitions, along with logicals
   integer, dimension(3) :: partition
   logical :: isInGrp
   type(MPI_Group) :: grp
   
   !> These are the two configs
   type(config) :: cfg
   
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
      call param_read('Gas Phase Partition',partition)

      ! Create an MPI group along with logical for the first grid on the lowest ranks
      grange(:,1)=[0,product(partition)-1,1]
      call MPI_Group_range_incl(group,1,grange,grp,ierr)
      isInGrp=.false.; if (rank.le.product(partition)-1) isInGrp=.true.
      
      ! Create an MPI group along with logical for the second grid on the highest ranks
      grange(:,1)=[nproc-product(partition),nproc-1,1]
      call MPI_Group_range_incl(group,1,grange,grp,ierr)
      isInGrp=.false.; if (rank.ge.nproc-product(partition)) isInGrp=.true.
      
      ! Create gas phase grid from input params
      if (isInGrp) then
         createrid: block
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
            cfg=config(grp=grp,decomp=partition,grid=grid)
            ! Do not place walls
            cfg%VF=1.0_WP
         end block createrid
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
