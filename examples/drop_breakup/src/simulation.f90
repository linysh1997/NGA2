!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,      only: WP
   use flow_class,     only: flow
   use droplet_class,  only: droplet
   use coupler_class,  only: coupler
   implicit none
   private
   
   !> Flow simulation
   type(flow) :: airflow
   logical :: isInFlowGrp
   
   !> Droplet simulation
   type(droplet) :: drop
   logical :: isInDropGrp
   
   !> Coupler from airflow to drop
   type(coupler) :: xcpl,ycpl,zcpl
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group,MPI_Group_range_incl
      use param, only: param_read
      implicit none
      type(MPI_Group) :: flow_group
      type(MPI_Group) :: drop_group

      create_flow_group: block
         use parallel, only: group,comm,nproc,rank
         integer :: ierr
         integer, dimension(3,1) :: grange
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Whole Domain Partition',partition)
         grange(:,1)=[0,product(partition)-1,1]
         call MPI_Group_range_incl(group,1,grange,flow_group,ierr)
         isInFlowGrp=.false.; if (rank.le.product(partition)-1) isInFlowGrp=.true.
      end block create_flow_group

      create_drop_group: block
         use parallel, only: group,comm,nproc,rank
         integer :: ierr
         integer, dimension(3,1) :: grange
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Droplet Domain Partition',partition)
         grange(:,1)=[nproc-product(partition),nproc-1,1]
         call MPI_Group_range_incl(group,1,grange,drop_group,ierr)
         isInDropGrp=.false.; if (rank.ge.nproc-product(partition)) isInDropGrp=.true.
      end block create_drop_group

      ! Initialize airflow simulation
      if (isInFlowGrp) call airflow%init(flow_group,isInFlowGrp)
      ! Initialize droplet simulation
      if (isInDropGrp) call drop%init(drop_group,isInDropGrp)

      ! If restarting, the domains could be out of sync, so resync
      ! time by forcing airflow to be at same time as droplet
      airflow%time%t=drop%time%t
      
      ! Initialize couplers from airflow to drop
      create_coupler: block
         use parallel, only: group
         xcpl=coupler(src_grp=flow_group,dst_grp=drop_group,name='airflow2drop')
         ycpl=coupler(src_grp=flow_group,dst_grp=drop_group,name='airflow2drop')
         zcpl=coupler(src_grp=flow_group,dst_grp=drop_group,name='airflow2drop')
         if(isInFlowGrp) then 
            call xcpl%set_src(airflow%cfg,'x')
            call ycpl%set_src(airflow%cfg,'y')
            call zcpl%set_src(airflow%cfg,'z')
         end if
         if(isInDropGrp) then
            call xcpl%set_dst(drop%cfg,'x')
            call ycpl%set_dst(drop%cfg,'y')
            call zcpl%set_dst(drop%cfg,'z')
         end if
         call xcpl%initialize()
         call ycpl%initialize()
         call zcpl%initialize()
      end block create_coupler
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      ! Airflow drives overall time integration
      if (isInFlowGrp) then
         do while (.not.airflow%time%done())
            ! Advance airflow simulation
            call airflow%step()
            ! Advance droplet simulation until it's caught up
            if (isInDropGrp) then
               do while (drop%time%t.le.airflow%time%t)
                  call drop%step()
               end do
            end if  
            
            ! Finish transfer
            ! Handle coupling between air flow and droplet
            coupling_a2d: block
               use tpns_class, only: bcond
               integer :: n,i,j,k
               type(bcond), pointer :: mybc_flow
               ! Exchange data using cplx/y/z couplers
               if(isInFlowGrp) call xcpl%push(airflow%fs%U); call xcpl%transfer(); if(isInDropGrp) call xcpl%pull(drop%resU)
               if(isInFlowGrp) call ycpl%push(airflow%fs%V); call ycpl%transfer(); if(isInDropGrp) call ycpl%pull(drop%resV)
               if(isInFlowGrp) call zcpl%push(airflow%fs%W); call zcpl%transfer(); if(isInDropGrp) call zcpl%pull(drop%resW)
               ! Apply time-varying Dirichlet conditions
               call airflow%fs%get_bcond('inflow',mybc_flow)
               do n=1,mybc_flow%itr%no_
                  i=mybc_flow%itr%map(1,n); j=mybc_flow%itr%map(2,n); k=mybc_flow%itr%map(3,n)
                  airflow%fs%U(i  ,j,k)=airflow%resU(i  ,j,k)*sum(airflow%fs%itpr_x(:,i  ,j,k)*airflow%cfg%VF(i-1:i,    j,    k))
                  airflow%fs%V(i-1,j,k)=airflow%resV(i-1,j,k)*sum(airflow%fs%itpr_y(:,i-1,j,k)*airflow%cfg%VF(i-1  ,j-1:j,    k))
                  airflow%fs%W(i-1,j,k)=airflow%resW(i-1,j,k)*sum(airflow%fs%itpr_z(:,i-1,j,k)*airflow%cfg%VF(i-1  ,j    ,k-1:k))
               end do
            end block coupling_a2d
         end do
   end if
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize droplet simulation
      if (isInDropGrp) call drop%final()
      
      ! Finalize airflow simulation
      if (isInFlowGrp) call airflow%final()
      
   end subroutine simulation_final
   
   
end module simulation