program Rough_temp
!code implementing coupled LBM on square lattice with D2Q9 for thermal flow over a rough surface
implicit none

   !subroutines used in the program
   interface
      subroutine streaming(f, c, u_x, u_y, state)
         real(8), dimension(:,:,0:), intent (in out) :: f, c
         real(8), dimension(:,:), intent(in) :: u_x, u_y
         integer, dimension(:,:), intent(in) :: state
      end subroutine streaming
      subroutine collision(f, c, f_eq, c_eq, tau, omega_d, alpha_g, rho, temperature, state)
         real(8), dimension(:,:,0:), intent(in out) :: f, c
         real(8), dimension(:,:,0:), intent(in) :: f_eq, c_eq
         real(8), intent(in) :: tau, omega_d, alpha_g
         real(8), dimension(:,:), intent(in) :: rho, temperature
         integer, dimension(:,:), intent(in) :: state
      end subroutine 
      subroutine equilibrium(rho, temperature, u_x, u_y, f_eq, c_eq, tau, omega_d, force, state)
         real(8), dimension(:,:), intent(in) :: rho, u_x, u_y, temperature
         real(8), dimension(:,:,:), allocatable, intent(out) :: f_eq, c_eq
         real(8), intent(in) :: force, tau, omega_d
         integer, dimension(:,:), intent(in) :: state
      end subroutine
      subroutine moment_update(f, c, rho, u_x, u_y, temperature, force, state)
         real(8), dimension(:,:,0:), intent(in) :: f, c
         integer, dimension(:,:), intent(in) :: state
         real(8), dimension(:,:), allocatable, intent(out) :: rho, u_x, u_y, temperature
         real(8), intent(in) :: force
      end subroutine
      subroutine LB_steps(n_x, n_y, steps, force, tau, omega_d, alpha_g, f, c, t, state)
         integer, intent(in) :: n_x, n_y, steps
         real(8), intent(in) :: force, tau, omega_d, alpha_g
         real(8), dimension(:,:,0:), intent(in out) :: f, c
         integer, intent(in out) :: t
         integer, dimension(:,:), intent(in) :: state
      end subroutine
   end interface

   integer, parameter :: n_x = 802, n_y = 161  !lattice size
   real(8), parameter :: force = 5d-8, tau = 2, viscosity = (tau - 0.5) / 3, omega_d = 1.98, alpha_g= 1d-4 
                         !force is the constant horizontal forcing, tau is the relaxation time for the f populations, omega_d is the
                         !time step size divided by the relaxation time for the c populations, alpha_g is the parametere controlling
                         !gravity
   real(8), dimension(n_x, n_y, 0:8) :: f, c  !populations
   real(8), dimension(:,:), allocatable :: rho, u_x, u_y, temperature  !macroscopic quantities: density, velocity and temperature
   real(8), dimension(:,:,:), allocatable :: f_eq, c_eq !equilibrium populations
   integer, dimension(n_x, n_y) :: state  !array labeling the fluid nodes
   integer :: t, i, j
   real(8) :: tot_temp, tot_density !temperature integrated over the system and total fluid mass
   
   !reading the array labeling the nodes: 0 for the solid nodes, 1 for the fluid nodes in the bulk and the other numbers to
   !distinguish the shape of the wall around boundary nodes
   open(12, file = "state.dat", status = "old")
      read(12, *) state
   close(12)

!   state = 1; state(:,1) = 13   !when this line is uncommented it sets the bottom surface to be flat and at the bottom of the
                                 !lattice
   
   !shaping the macroscopic quantity arrays with lattice size
   allocate(rho(n_x,n_y), u_x(n_x,n_y), u_y(n_x, n_y), temperature(n_x, n_y))

   !initializing the macroscopic quantities
   rho = 1; u_x = 0; u_y = 0; temperature = 0; t = 0
   
   !opening the files on which the data is going to be written
   open(18, file= "concentration_field.dat", status = "replace")
   write(18, *) "# force ", force, " tau ", tau, " omega_d ", omega_d, "alpha_g", alpha_g
   open(19, file= "velocity_field.dat", status = "replace")
   write(19, *) "# force ", force, " tau ", tau, " omega_d ", omega_d, "alpha_g", alpha_g
   open(20, file= "average.dat", status = "replace")
   write(20, *) "# force ", force, " tau ", tau, " omega_d ", omega_d, "alpha_g", alpha_g

   !initializing the populations to be the local equilibrium for the initial values of the macroscopic quantities
   call equilibrium(rho, temperature, u_x, u_y, f_eq, c_eq, tau, omega_d, force, state)
   f = f_eq
   c = c_eq

!   open(17, file = "populations.dat", status = "old")      !used when continuing a previous simulation
!   open(16, file = "concentrations.dat", status = "old")   !reading the populations from the last step of a previous simulations to 
!      read(17, *) f                                        !continue it
!      read(16, *) c                                        ! -//-
!   close(17); close(16)                                    ! -//-

   write(18, *) "  "                                                                          !used only when starting a new simulation
   write(18, *) "  "                                                                          !writing on file the initial values of 
   write(18, *) "# t = ", t                                                                   !the macroscopi quantities
   do i = 1, n_x                                                                              ! -//-
      do j = 1, n_y                                                                           ! -//-
         if (state(i,j) .ne. 0) then                                                          ! -//-
            write(18, *) i, j, 1, temperature(i,j)                                            ! -//-
         end if                                                                               ! -//-
      end do                                                                                  ! -//-
   end do                                                                                     ! -//-

   write(19, *) " "; write(19,*) " "; write(19, *) "# t = ", t                                ! -//-
   do i = 1, n_x                                                                              ! -//-
      do j = 1, n_y                                                                           ! -//-
         if (state(i,j) .ne. 0) then                                                          ! -//-
            write(19, *) i, j, u_x(i, j) + force/(2*rho(i, j)), u_y(i,j)                      ! -//-
         end if                                                                               ! -//-
      end do                                                                                  ! -//-
   end do                                                                                     ! -//-

   !calling the subroutine that implements a number of time steps (1000 in this case) and writes data on file at the end
   !i controls how many times the subroutine is called
   do i = 1, 200
      call LB_steps(n_x, n_y, 1000, force, tau, omega_d, alpha_g, f, c, t, state)
   end do
   close(18)
     
   !write the last values of the populations in order to be able to continue this simulation at a later time
   open(17, file = "populations.dat", status = "replace")
   open(16, file = "concentrations.dat", status = "replace")
      write(17, *) f
      write(16, *) c
   close(17); close(16)

end program Rough_temp

subroutine LB_steps(n_x, n_y, steps, force, tau, omega_d, alpha_g, f, c, t, state)
implicit none
!implementing a number of time steps specified by the variable steps and printing on file the temperature and velocity fields,
!temperature integrated over the system, speed averaged over the system

   integer, intent(in) :: n_x, n_y, steps  !lattice size and number of steps implemented by the subroutine
   real(8), intent(in) :: force, tau, omega_d, alpha_g !force is the constant horizontal forcing, tau is the relaxation time for the f
                                                       !populations, omega_d is the time step size divided by the relaxation time
                                                       !for the c populations, alpha_g is the parameter controlling gravity
   real(8), dimension(:,:,0:), intent(in out) :: f, c !populations
   integer, intent(in out) :: t   !variable tracking time
   integer, dimension(:,:), intent(in) :: state !array labeling the nodes

   real(8), dimension(:,:), allocatable :: rho, u_x, u_y, temperature !macroscopi quantities
   real(8), dimension(:,:,:), allocatable :: f_eq, c_eq  !equilibrium populations
   real(8), dimension(n_x, n_y) :: abso_vel !speed
   real(8) :: average_vel, tot_temp, tot_density !average_vel is the speed averaged over the system, tot_temp is the temperature 
                                                 !integrated over the system, tot_density is the total fluid mass in the 
                                                 !system
   integer :: i, j

   !all the subroutines used in this subroutine
   interface
      subroutine streaming(f, c, u_x, u_y, state)
         real(8), dimension(:,:,0:), intent (in out) :: f, c
         real(8), dimension(:,:), intent(in) :: u_x, u_y
         integer, dimension(:,:), intent(in) :: state
      end subroutine streaming
      subroutine collision(f, c, f_eq, c_eq, tau, omega_d, alpha_g, rho, temperature, state)
         real(8), dimension(:,:,0:), intent(in out) :: f, c
         real(8), dimension(:,:,0:), intent(in) :: f_eq, c_eq
         real(8), intent(in) :: tau, omega_d, alpha_g
         real(8), dimension(:,:), intent(in) :: rho, temperature
         integer, dimension(:,:), intent(in) :: state
      end subroutine 
      subroutine equilibrium(rho, temperature, u_x, u_y, f_eq, c_eq, tau, omega_d, force, state)
         real(8), dimension(:,:), intent(in) :: rho, u_x, u_y, temperature
         real(8), dimension(:,:,:), allocatable, intent(out) :: f_eq, c_eq
         real(8), intent(in) :: force, tau, omega_d
         integer, dimension(:,:), intent(in) :: state
      end subroutine
      subroutine moment_update(f, c, rho, u_x, u_y, temperature, force, state)
         real(8), dimension(:,:,0:), intent(in) :: f, c
         integer, dimension(:,:), intent(in) :: state
         real(8), dimension(:,:), allocatable, intent(out) :: rho, u_x, u_y, temperature
         real(8), intent(in) :: force
      end subroutine
   end interface

   !implement a number of time steps specified by the variable steps
   do i = 1, steps
      call moment_update(f, c, rho, u_x, u_y, temperature, force, state)
      call equilibrium(rho, temperature, u_x, u_y, f_eq, c_eq, tau, omega_d, force, state)
      call collision(f, c, f_eq, c_eq, tau, omega_d, alpha_g, rho, temperature, state)
      call streaming(f, c, u_x, u_y, state)
   end do

   abso_vel = sqrt((u_x + force/(2*rho)) ** 2 + u_y** 2) !calculate the speed for each node
   t = t + steps  !updating time
   
   !calculating integrated temperature and total fluid mass
   tot_temp = 0; tot_density = 0
   do i = 1, n_x
      do j = 1, n_y
         if (state(i,j) .ne. 0) then
            tot_temp = tot_temp + temperature(i,j)
            tot_density = tot_density + rho(i,j)
         end if
      end do
   end do

   !writing the speed and temperature field on file
   write(18, *) "  "; write(18, *) "  "; write(18, *) "# t = ", t
   do i = 1, n_x
      do j = 1, n_y
         if (state(i,j) .ne. 0) then
            write(18, *) i, j, abso_vel(i,j), temperature(i,j)
         end if
      end do
   end do

   !writing the velocity field on file
   write(19, *) "  "; write(19, *) "  "; write(19, *) "# t = ", t
   do i = 1, n_x
      do j = 1, n_y
         if (state(i,j) .ne. 0) then
            write(19, *) i, j, u_x(i, j) + force/(2*rho(i, j)), u_y(i,j)
         end if
      end do
   end do
 
   !calculating the average speed
   average_vel = 0
   do i = 1, n_x
      do j = 1, n_y
         if (state(i,j) .ne. 0) then
            average_vel = average_vel + abso_vel(i,j)
         end if
      end do
   end do
   average_vel = average_vel / (n_x * n_y)

   !writing time, average velocity, integrated temperature and total fluid mass of file
   write(20, *) t, average_vel, tot_temp, tot_density
end subroutine

subroutine streaming(f, c, u_x, u_y, state)
implicit none
!streaming of the populations on a D2Q9 lattice with periodic boundary conditions

   real(8), dimension(0:,0:,0:), intent(in out) :: f, c !populations
   real(8), dimension(0:,0:), intent(in) :: u_x, u_y  !velocity
   integer, dimension(0:,0:), intent(in) :: state  !array labeling the nodes

   real(8), dimension(:,:,:), allocatable :: temp ! temporary array
   integer :: n_x, n_y, i, j !number of lattice sites

   !extracting the lattice size from the population array
   n_x = size(f, 1)
   n_y = size(f, 2)

   !shaping the temporary array with the same size as the population array
   allocate(temp(0:n_x-1,0:n_y-1,0:8))

   temp = f !copying the previous populations on the temporary array to prevent writing over populations that have not yet been 
            !processed

   !streaming the f populations, depending on whether the node is in the bulk or on the bottom surface and the shape of the boundary
   !around the node
   !half-way bounce-back is used to impose no-slip on the bottom surface
   do i = 0, n_x - 1
      do j = 0, n_y - 2
         select case(state(i,j))
            case(1) 
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 2) = temp(i, modulo(j-1, n_y), 2)
               f(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
               f(i, j, 6) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(13)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 2) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 6) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(11)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 2) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 3) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
               f(i, j, 6) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(10)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 2) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 3) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
               f(i, j, 6) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
               f(i, j, 7) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(12)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 2) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 3) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 6) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(3)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 2) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 3) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(i, j, 5)
               f(i, j, 6) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
               f(i, j, 7) = temp(i, j, 7)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(5)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 2) = temp(i, modulo(j-1, n_y), 2)
               f(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 6) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(8)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 2) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 6) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(9)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 2) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 6) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(7)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 2) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 6) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
            case(2)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 2) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 6) = temp(i,j,6)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(i,j,8)
           case(14)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 2) = temp(i, modulo(j-1, n_y), 2)
               f(i, j, 3) = temp(modulo(i, n_x), j, 3)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 6) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
            case(16)
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 2) = temp(i, modulo(j-1, n_y), 2)
               f(i, j, 3) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
               f(i, j, 6) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
               f(i, j, 7) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
            case(6)           
               f(i, j, 0) = temp(i, j, 0)
               f(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
               f(i, j, 2) = temp(i, modulo(j-1, n_y), 2)
               f(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
               f(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
               f(i, j, 5) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
               f(i, j, 6) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
               f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
               f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
         end select
      end do 
   end do
     
   !streaming the f populations on the top wall and imposing free-slip boundary conditions
   j = n_y - 1
   do i = 0, n_x - 1
      f(i, j, 0) = temp(i, j, 0)
      f(i, j, 1) = temp(modulo(i-1,n_x), j, 1)
      f(i, j, 2) = temp(i, modulo(j-1,n_y), 2)
      f(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
      f(i, j, 4) = temp(i, modulo(j-1,n_y), 2)
      f(i, j, 5) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
      f(i, j, 6) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
      f(i, j, 7) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
      f(i, j, 8) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
   end do

   !streaming the f populations at the outlet and imposing zero-gradient conditions
   i = n_x - 1
   do j = 0, n_y - 1
      if (state(i,j) .ne. 0) then
         f(i, j, 0) = f(i-1, j, 0)
         f(i, j, 1) = f(i-1, j, 1)
         f(i, j, 2) = f(i-1, j, 2)
         f(i, j, 3) = f(i-1, j, 3)
         f(i, j, 4) = f(i-1, j, 4)
         f(i, j, 5) = f(i-1, j, 5)
         f(i, j, 6) = f(i-1, j, 6)
         f(i, j, 7) = f(i-1, j, 7)
         f(i, j, 8) = f(i-1, j, 8)
      end if
   end do

   temp = c  !copying the previous c populations on the temperory array

   !streaming the c populations, depending on wheter the node is in the bulk or on the bottom surface and the shape of the 
   !boundary around the node   
   !the constant value 0 is imposed on the bottom surface
   do i = 1, n_x - 1
      do j = 0, n_y - 2
         if (state(i,j) .eq. 1) then
            c(i, j, 0) = temp(i, j, 0)
            c(i, j, 1) = temp(modulo(i-1, n_x), j, 1)
            c(i, j, 2) = temp(i, modulo(j-1, n_y), 2)
            c(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
            c(i, j, 4) = temp(i, modulo(j+1, n_y), 4)
            c(i, j, 5) = temp(modulo(i-1,n_x), modulo(j-1,n_y), 5)
            c(i, j, 6) = temp(modulo(i+1,n_x), modulo(j-1,n_y), 6)
            c(i, j, 7) = temp(modulo(i+1,n_x), modulo(j+1,n_y), 7)
            c(i, j, 8) = temp(modulo(i-1,n_x), modulo(j+1,n_y), 8)
         else if (state(i,j) .ne. 0) then
            c(i, j, 0) = 0 * 4d0 / 9
            c(i, j, 1) = (0 + 3 * u_x(i,j)) / 9
            c(i, j, 2) = (0 + 3 * u_y(i,j)) / 9
            c(i, j, 3) = (0 - 3 * u_x(i,j)) / 9
            c(i, j, 4) = (0 - 3 * u_y(i,j)) / 9
            c(i, j, 5) = (0 + 3 * (u_x(i,j) + u_y(i,j))) / 36
            c(i, j, 6) = (0 - 3 * (u_x(i,j) - u_y(i,j))) / 36
            c(i, j, 7) = (0 - 3 * (u_x(i,j) + u_y(i,j))) / 36
            c(i, j, 8) = (0 + 3 * (u_x(i,j) - u_y(i,j))) / 36
         end if
      end do 
   end do

   !streaming populations on the top wall and imposing free-slip conditions
   j = n_y -1
   do i = 1, n_x - 1
      c(i, j, 0) = temp(i, j, 0)
      c(i, j, 1) = temp(modulo(i-1,n_x), j, 1)
      c(i, j, 2) = temp(i, modulo(j-1,n_y), 2)
      c(i, j, 3) = temp(modulo(i+1, n_x), j, 3)
      c(i, j, 4) = temp(i, modulo(j-1,n_y), 2)
      c(i, j, 5) = temp(modulo(i-1,n_x), modulo(j-2,n_y), 5)
      c(i, j, 6) = temp(modulo(i,n_x), modulo(j-2,n_y), 6)
      c(i, j, 7) = temp(modulo(i,n_x), modulo(j-2,n_y), 6)
      c(i, j, 8) = temp(modulo(i-1,n_x), modulo(j-2,n_y), 5)
   end do

   !streaming populations at the inlet and imposing temperature C=1
   i = 0
   do j = 0, n_y - 1
      if (state(i,j) .ne. 0) then
         c(i, j, 0) = 1 * 4d0 / 9
         c(i, j, 1) = (1 + 3 * u_x(i,j)) / 9
         c(i, j, 2) = (1 + 3 * u_y(i,j)) / 9
         c(i, j, 3) = (1 - 3 * u_x(i,j)) / 9
         c(i, j, 4) = (1 - 3 * u_y(i,j)) / 9
         c(i, j, 5) = (1 + 3 * (u_x(i,j) + u_y(i,j))) / 36
         c(i, j, 6) = (1 - 3 * (u_x(i,j) - u_y(i,j))) / 36
         c(i, j, 7) = (1 - 3 * (u_x(i,j) + u_y(i,j))) / 36
         c(i, j, 8) = (1 + 3 * (u_x(i,j) - u_y(i,j))) / 36
      end if
   end do

   !streaming populations at the outlet and imposing zero-gradient conditions
   i = n_x - 1
   do j = 0, n_y - 1
      if (state(i,j) .ne. 0) then
         c(i, j, 0) = c(i-1, j, 0)
         c(i, j, 1) = c(i-1, j, 1)
         c(i, j, 2) = c(i-1, j, 2)
         c(i, j, 3) = c(i-1, j, 3)
         c(i, j, 4) = c(i-1, j, 4)
         c(i, j, 5) = c(i-1, j, 5)
         c(i, j, 6) = c(i-1, j, 6)
         c(i, j, 7) = c(i-1, j, 7)
         c(i, j, 8) = c(i-1, j, 8)
      end if
   end do

end subroutine streaming

subroutine collision(f, c, f_eq, c_eq, tau, omega_d, alpha_g, rho, temperature, state)
implicit none
!calculates the new populations on a D2Q9 lattice after collision

   real(8), dimension(:,:,0:), intent(in out) :: f, c !populations
   real(8), dimension(:,:,0:), intent(in) :: f_eq, c_eq  !equilibrium distributions
   real(8), intent(in) :: tau, omega_d, alpha_g !tau is the relaxation time for the f populations, omega_d is the time step
                                                        !divided by the relaxation time for the c populations, alpha_g is the parameter 
                                                        !controlling gravity and force is the constant horizontal forcing
   real(8), dimension(:,:), intent(in) :: rho, temperature !density and temperature
   integer, dimension(:,:), intent(in) :: state !array labeling the nodes

   real(8) :: omega, omega_p, omega_d_p  !omega is the time step divided by the relaxation time and omega_p is 1-omega
   integer :: i, j, k

   omega = 1/tau
   omega_p = 1 - omega
   omega_d_p = 1 - omega_d

   !write in f the populations after collision
   do i = 1, size(f_eq, 1)
      do j = 1, size(f_eq, 2)
         if (state(i,j) .ne. 0) then
            f(i, j, 0) = omega_p * f(i,j,0) + omega * f_eq(i,j,0)
            f(i, j, 1) = omega_p * f(i,j,1) + omega * f_eq(i,j,1)
            f(i, j, 2) = omega_p * f(i,j,2) + omega * f_eq(i,j,2) + alpha_g * rho(i,j) * temperature(i,j) / 27
            f(i, j, 3) = omega_p * f(i,j,3) + omega * f_eq(i,j,3)
            f(i, j, 4) = omega_p * f(i,j,4) + omega * f_eq(i,j,4) - alpha_g * rho(i,j) * temperature(i,j) / 27
            f(i, j, 5) = omega_p * f(i,j,5) + omega * f_eq(i,j,5) + alpha_g * rho(i,j) * temperature(i,j) / 108
            f(i, j, 6) = omega_p * f(i,j,6) + omega * f_eq(i,j,6) + alpha_g * rho(i,j) * temperature(i,j) / 108
            f(i, j, 7) = omega_p * f(i,j,7) + omega * f_eq(i,j,7) - alpha_g * rho(i,j) * temperature(i,j) / 108
            f(i, j, 8) = omega_p * f(i,j,8) + omega * f_eq(i,j,8) - alpha_g * rho(i,j) * temperature(i,j) / 108
            do k = 0, 8
               c(i, j, k) = omega_d_p * c(i, j, k) + omega_d * c_eq(i, j, k)
            end do
         end if
     end do
   end do

end subroutine collision

subroutine moment_update(f, c, rho, u_x, u_y, temperature, force, state)
implicit none
!calculate the macroscopic density, velocity and temperature on a D2Q9 lattice

   real(8), dimension(:,:,0:), intent(in) :: f, c   !populations
   integer, dimension(:,:), intent(in) :: state !array labeling the nodes
   real(8), dimension(:,:), allocatable, intent(out) :: rho, u_x, u_y, temperature !macroscopic density, velocity and temperature
   real(8), intent(in) :: force

   integer :: n_x, n_y, i, j ! number of lattice sites

   !extrating the size of the lattice from the populations array
   n_x = size(f,1)
   n_y = size(f,2)

   !shaping the density, velocity and concentration arrays with size of the lattice
   allocate(rho(n_x, n_y), u_x(n_x, n_y), u_y(n_x, n_y), temperature(n_x, n_y))
   rho = 0; u_x = 0; u_y = 0

   !calculate density and velocity at each lattice site
   do i = 1, n_x
      do j = 1, n_y
         if (state(i, j) .ne. 0) then 
            rho(i, j) = f(i, j, 0) + f(i, j, 1) + f(i, j, 2) + f(i, j, 3) + f(i, j, 4) + f(i, j, 5) + f(i, j, 6)&
                + f(i, j, 7) + f(i, j, 8)
            u_x(i, j) = (f(i, j, 1) + f(i, j, 5) + f(i, j, 8) - f(i, j, 3) - f(i, j, 6) - f(i, j, 7)) / rho(i, j)
            u_y(i, j) = (f(i, j, 2) + f(i, j, 5) + f(i, j, 6) - f(i, j, 4) - f(i, j, 7) - f(i, j, 8)) / rho(i, j)
            temperature(i,j) = c(i,j,0) + c(i,j,1) + c(i,j,2) + c(i,j,3) + c(i,j,4) + c(i,j,5) + c(i,j,6) + c(i,j,7) + c(i,j,8)
         end if
      end do
   end do

   !inlet
   u_x(1,:) = force / 2
   u_y(1,:) = 0

end subroutine moment_update

subroutine equilibrium(rho, temperature, u_x, u_y, f_eq, c_eq, tau, omega_d, force, state)
implicit none
!calculate the equilibrium distributions on a D2Q9 latticeusing the given density and distribution

   real(8), dimension(:,:), intent(in) :: rho, u_x, u_y, temperature  !macroscopic density, velocity and temperature
   real(8), dimension(:,:,:), allocatable, intent(out) :: f_eq, c_eq   !equilibrium distributions
   real(8), intent(in) :: force, tau, omega_d
   integer, dimension(:,:), intent(in) :: state

   integer, parameter :: c_s_2 = 3, c_s_4 = 9   !the constants c_s^(-2) and c_s^(-4), where c_s is the sound speed
   integer :: n_x, n_y, i, j  !number of lattice sites
   real(8) :: u_sq, u_x_eq   !modulus squared of the velocity

   !extracting the size of the lattice from the density array
   n_x = size(rho, 1)
   n_y = size(rho, 2)

   !shaping the equilibrium distribution arrays with the size of the lattice and the number of discrete velocities
   allocate(f_eq(n_x, n_y, 0:8), c_eq(n_x, n_y, 0:8))
   f_eq = 0

   
   !calculate and write in f_eq the equilibrium distributions
   do i = 1, n_x
      do j = 1, n_y
         if (state(i, j) .ne. 0) then
            u_x_eq = u_x(i, j) + (tau * force / rho(i,j)) !the forcing is introduced in the equilibrium distribution
            u_sq = (u_x_eq ** 2) + (u_y(i,j) ** 2) !u^2
            f_eq(i, j, 0) = 2 * rho(i,j) * (2 - c_s_2 * u_sq) / 9
            f_eq(i, j, 1) = rho(i, j) * (2 + 2*c_s_2*u_x_eq + c_s_4*(u_x_eq**2) - c_s_2*u_sq) / 18
            f_eq(i, j, 2) = rho(i, j) * (2 + 2*c_s_2*u_y(i,j) + c_s_4*(u_y(i,j)**2) - c_s_2*u_sq) / 18
            f_eq(i, j, 3) = rho(i, j) * (2 - 2*c_s_2*u_x_eq + c_s_4*(u_x_eq**2) - c_s_2*u_sq) / 18
            f_eq(i, j, 4) = rho(i, j) * (2 - 2*c_s_2*u_y(i,j) + c_s_4*(u_y(i,j)**2) - c_s_2*u_sq) / 18
            f_eq(i, j, 5) = rho(i, j) * (1 + c_s_2*(u_x_eq+u_y(i,j)) + c_s_4*u_x_eq*u_y(i,j) + (c_s_4-c_s_2)*u_sq/2) / 36
            f_eq(i, j, 6) = rho(i, j) * (1 - c_s_2*(u_x_eq-u_y(i,j)) - c_s_4*u_x_eq*u_y(i,j) + (c_s_4-c_s_2)*u_sq/2) / 36
            f_eq(i, j, 7) = rho(i, j) * (1 - c_s_2*(u_x_eq+u_y(i,j)) + c_s_4*u_x_eq*u_y(i,j) + (c_s_4-c_s_2)*u_sq/2) / 36
            f_eq(i, j, 8) = rho(i, j) * (1 + c_s_2*(u_x_eq-u_y(i,j)) - c_s_4*u_x_eq*u_y(i,j) + (c_s_4-c_s_2)*u_sq/2) / 36
            u_sq = (u_x(i,j)**2) + (u_y(i,j) ** 2)
            c_eq(i, j, 0) = 2 * temperature(i,j) * (2 - c_s_2 * u_sq) / 9
            c_eq(i, j, 1) = temperature(i, j) * (2 + 2*c_s_2*u_x(i,j) + c_s_4*(u_x(i,j)**2) - c_s_2*u_sq) / 18
            c_eq(i, j, 2) = temperature(i, j) * (2 + 2*c_s_2*u_y(i,j) + c_s_4*(u_y(i,j)**2) - c_s_2*u_sq) / 18
            c_eq(i, j, 3) = temperature(i, j) * (2 - 2*c_s_2*u_x(i,j) + c_s_4*(u_x(i,j)**2) - c_s_2*u_sq) / 18
            c_eq(i, j, 4) = temperature(i, j) * (2 - 2*c_s_2*u_y(i,j) + c_s_4*(u_y(i,j)**2) - c_s_2*u_sq) / 18
            c_eq(i, j, 5) = temperature(i, j) * (1 + c_s_2*(u_x(i,j)+u_y(i,j)) + c_s_4*u_x(i,j)*u_y(i,j) &
               + (c_s_4-c_s_2)*u_sq/2) / 36
            c_eq(i, j, 6) = temperature(i, j) * (1 - c_s_2*(u_x(i,j)-u_y(i,j)) - c_s_4*u_x(i,j)*u_y(i,j) &
               + (c_s_4-c_s_2)*u_sq/2) / 36
            c_eq(i, j, 7) = temperature(i, j) * (1 - c_s_2*(u_x(i,j)+u_y(i,j)) + c_s_4*u_x(i,j)*u_y(i,j) &
               + (c_s_4-c_s_2)*u_sq/2) / 36
            c_eq(i, j, 8) = temperature(i, j) * (1 + c_s_2*(u_x(i,j)-u_y(i,j)) - c_s_4*u_x(i,j)*u_y(i,j) &
               + (c_s_4-c_s_2)*u_sq/2) / 36            
         end if
      end do
   end do

end subroutine equilibrium

