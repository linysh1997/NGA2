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
   logical :: isInHITGrp
   
   !> Droplet atomization simulation
   type(droplet) :: drop
   
   !> Coupler from airflow to drop
   type(coupler) :: xcpl,ycpl,zcpl
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      ! type(MPI_Group) :: flow_group
      
      ! Initialize air flow simulation
      call airflow%init()

      ! Initialize droplet atomization simulation
      call drop%init()

      ! If restarting, the domains could be out of sync, so resync
      ! time by forcing airflow to be at same time as droplet atomization
      airflow%time%t=drop%time%t  
      
      ! Initialize couplers from airflow to drop
      create_coupler: block
         use parallel, only: group
         xcpl=coupler(src_grp=group,dst_grp=group,name='airflow2drop')
         ycpl=coupler(src_grp=group,dst_grp=group,name='airflow2drop')
         zcpl=coupler(src_grp=group,dst_grp=group,name='airflow2drop')
         call xcpl%set_src(airflow%cfg,'x')
         call ycpl%set_src(airflow%cfg,'y')
         call zcpl%set_src(airflow%cfg,'z')
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
         ! Advance airflow simulation until it's caught up
         do while (airflow%time%t.le.drop%time%t)
            call airflow%step()
         end do
         
         ! Finish transfer
         ! Handle coupling between injector and atomization
         coupling_a2d: block
            use tpns_class, only: bcond
            integer :: n,i,j,k
            type(bcond), pointer :: mybc
            ! Exchange data using cplx/y/z couplers
            call xcpl%push(airflow%fs%U); call xcpl%transfer(); call xcpl%pull(drop%resU)
            call ycpl%push(airflow%fs%V); call ycpl%transfer(); call ycpl%pull(drop%resV)
            call zcpl%push(airflow%fs%W); call zcpl%transfer(); call zcpl%pull(drop%resW)
            ! Apply time-varying Dirichlet conditions
            call drop%fs%get_bcond('inflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               drop%fs%U(i  ,j,k)=drop%resU(i  ,j,k)*sum(drop%fs%itpr_x(:,i  ,j,k)*drop%cfg%VF(i-1:i,    j,    k))
               drop%fs%V(i-1,j,k)=drop%resV(i-1,j,k)*sum(drop%fs%itpr_y(:,i-1,j,k)*drop%cfg%VF(i-1  ,j-1:j,    k))
               drop%fs%W(i-1,j,k)=drop%resW(i-1,j,k)*sum(drop%fs%itpr_z(:,i-1,j,k)*drop%cfg%VF(i-1  ,j    ,k-1:k))
            end do
         end block coupling_a2d
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize atomization simulation
      call drop%final()
      
      ! Finalize HIT simulation
      call airflow%final()
      
   end subroutine simulation_final
   
   
end module simulation