!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,      only: WP
   use flow_class,     only: flow
   use droplet_class,  only: droplet
   use coupler_class,  only: coupler
   implicit none
   private
   
   !> Flow simulation
   type(flow) :: gas
   logical :: isInHITGrp
   
   !> Droplet atomization simulation
   type(droplet) :: drop
   
   !> Coupler from gas to drop
   type(coupler) :: xcpl,ycpl,zcpl
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group
      implicit none
      ! type(MPI_Group) :: flow_group
      
      ! Initialize atomization simulation
      call drop%init()
      call gas%init()
      
      ! ! Create an MPI group using leftmost processors only
      ! create_flow_group: block
      !    use parallel, only: group,comm
      !    use mpi_f08,  only: MPI_Group_incl
      !    integer, dimension(:), allocatable :: ranks
      !    integer, dimension(3) :: coord
      !    integer :: n,ngrp,ierr,ny,nz
      !    ngrp=drop%cfg%npy*drop%cfg%npz
      !    allocate(ranks(ngrp))
      !    ngrp=0
      !    do nz=1,drop%cfg%npz
      !       do ny=1,drop%cfg%npy
      !          ngrp=ngrp+1
      !          coord=[0,ny-1,nz-1]
      !          call MPI_CART_RANK(drop%cfg%comm,coord,ranks(ngrp),ierr)
      !       end do
      !    end do
      !    call MPI_Group_incl(group,ngrp,ranks,flow_group,ierr)
      !    if (drop%cfg%iproc.eq.1) then
      !       isInHITGrp=.true.
      !    else
      !       isInHITGrp=.false.
      !    end if
      ! end block create_flow_group
      
      ! ! Prepare HIT simulation
      ! if (isInHITGrp) then
      !    prepare_flow: block
      !       real(WP) :: dt
      !       ! Initialize HIT
      !       call gas%init(group=flow_group,xend=drop%cfg%x(drop%cfg%imin))
      !       ! Run HIT until t/tau_eddy=20
      !       dt=0.15_WP*gas%cfg%min_meshsize/gas%Urms_tgt !< Estimate maximum stable dt
      !       do while (gas%time%t.lt.20.0_WP*gas%tau_tgt); call gas%step(dt); end do
      !    end block prepare_flow
      ! end if
      
      ! Initialize couplers from gas to drop
      create_coupler: block
         use parallel, only: group
         xcpl=coupler(src_grp=group,dst_grp=group,name='gas2drop')
         ycpl=coupler(src_grp=group,dst_grp=group,name='gas2drop')
         zcpl=coupler(src_grp=group,dst_grp=group,name='gas2drop')
         call xcpl%set_src(gas%cfg,'x')
         call ycpl%set_src(gas%cfg,'y')
         call zcpl%set_src(gas%cfg,'z')
         call xcpl%set_dst(drop%cfg,'x'); call xcpl%initialize()
         call ycpl%set_dst(drop%cfg,'y'); call ycpl%initialize()
         call zcpl%set_dst(drop%cfg,'z'); call zcpl%initialize()
      end block create_coupler
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Atomization drives overall time integration
      do while (.not.drop%time%done())
         ! Advance atomization simulation
         call drop%step()
         call gas%step()
         ! Finish transfer
         pull_velocity: block
            call xcpl%transfer(); call xcpl%pull(drop%resU)
            call ycpl%transfer(); call ycpl%pull(drop%resV)
            call zcpl%transfer(); call zcpl%pull(drop%resW)
         end block pull_velocity
         ! Apply time-dependent Dirichlet condition
         apply_boundary_condition: block
            use tpns_class, only: bcond
            type(bcond), pointer :: mybc
            integer :: n,i,j,k
            call drop%fs%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               drop%fs%U(i  ,j,k)=1.0_WP+drop%resU(i  ,j,k)
               drop%fs%V(i-1,j,k)=       drop%resV(i-1,j,k)
               drop%fs%W(i-1,j,k)=       drop%resW(i-1,j,k)
            end do
         end block apply_boundary_condition
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize atomization simulation
      call drop%final()
      
      ! Finalize HIT simulation
      if (isInHITGrp) call gas%final()
      
   end subroutine simulation_final
   
   
end module simulation