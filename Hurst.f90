program Hurst
!generates rough profile, discretizes it and labels the boundary depending on its local orientation
implicit none
   integer, parameter :: N = 400, scal = 2 !N is the number of points for which the height of the profile is generated
   real(8), parameter :: pi = 3.14159265358979 * 2d0 / N , Hu = 0.8, v_scaling = 1d-1 !Hu is the Hurst exponent of the rough profile
   real(8), dimension(0: N * scal) :: h_sc   !heights generated plus the interpolated ones
   integer, dimension(0: N * scal) :: h_discr !heights approximated to the closest integer
   real(8), dimension(0:N) :: delta, h !delta are the random phases in the Fourier series and h are the generated heights
   real(8), dimension(N) :: h_k  !Fourier coefficients
   logical, dimension(0: N * scal, 0: N * scal) :: border !logical array: true for border sites, false otherwise
   integer, dimension(0: N * scal, 0: N * scal) :: state, copy !temporary arrays to store lattice site labels
   integer, dimension(0: N * scal + 1, 0: N * scal) :: fin_state !final labels for lattice sites
   real(8) :: m, q
   complex(8) :: temp
   integer :: i, k, a

   h_k = (/ ((pi * k) ** (- 5d-1 - Hu), k = 1, N) /) !Fourier coefficients
   call random_number(delta) !random phase

   !generating the profile height for N+1 points with the Fourier transform method in R.F.Voss
   h = 0
   do i = 0, N
      do k = 1, N
         h(i) = h(i) + h_k(k) * sin(k * i * pi + N * pi &
         * delta(k))
      end do
      h(i) = h(i) * v_scaling
   end do

   !interpolating and taking scal-1 points in between previously generated points 
   h_sc(0) = h(0) * scal
   do k = 0, N-1
      m = h(k+1) - h(k)
      q = h(k) - m * k
      do i = 1, scal-1
         h_sc((k * scal) + i) = ((m * dble((k * scal) + i) / scal) + q) * scal
      end do
      h_sc((k * scal) + scal) = h(k + 1) * scal
   end do

   !approximating heights to closest integer
   h_discr = nint(h_sc)

   !border is a logical array of the size of the lattice; true for nodes on the boundary, false otherwise
   border = .false.
   border(0, h_discr(0) - minval(h_discr) + 1) = .true.
   do i = 1, N * scal
      if (abs(h_discr(i) - h_discr(i-1)) > 0) then
         do k = min(h_discr(i), h_discr(i-1)), max(h_discr(i), h_discr(i-1))
            border(i, k - minval(h_discr) +1) = .true.
         end do
      else
         border(i, h_discr(i) - minval(h_discr) + 1) = .true.
      end if
   end do

   !labeling the nodes: 0 for solid nodes, 1 for fluid nodes in the bulk, 2 for boundary nodes
   do k = 0, N * scal
      a = 0
      do i = 0, N * scal
         state(k, i) = a
         if (border(k, i)) then
            state(k, i) = 2
            a = 1
         end if
      end do
   end do

   !relabeling the boundary nodes base on the local orientation of the boundary
   do k = 0, (N * scal)
      do i = 0, N * scal
         if (border(k, i)) then
            if (border(k, i+1)) then
               if (border(k, i-1)) then
                  if ((state(modulo(k+1, (N*scal)+1), i) .eq. 0) .or. (state(modulo(k-1, (N*scal)+1), i) .eq. 1)) then
                     state(k, i) = 16
                  else if ((state(modulo(k-1, (N*scal)+1), i) .eq. 0) .or. (state(modulo(k+1, (N*scal)+1), i) .eq. 1)) then
                     state(k, i) = 14
                  end if
               else if ((border(modulo(k+1, (N*scal)+1), i)) .and. (.not.(border(modulo(k+1, (N*scal)+1), i+1)))) then
                  state(k, i) = 8
               else if (border(modulo(k-1, (N*scal)+1), i)) then 
                  state(k, i) = 11
               else if (border(modulo(k+1, (N*scal)+1), i)) then
                  state(k, i) = 8
               end if
            else if (border(k, i-1)) then
               if ((border(modulo(k-1, (N*scal)+1), i)) .and. ((state(modulo(k-1, (N*scal)+1), i) .eq. 8) .or. (state(modulo(k-1,&
               (N*scal)+1), i) .eq. 13) .or. (state(modulo(k-1, (N*scal)+1), i) .eq. 6))) then 
                  state(k, i) = 5
               else if (border(modulo(k+1, (N*scal)+1), i)) then
                  state(k, i) = 6
               end if
            else if (border(modulo(k+1, (N*scal)+1), i)) then
               if (border(modulo(k-1, (N*scal)+1), i)) then
                  state(k, i) = 13
               end if
            end if
         end if
      end do
   end do
   copy = state
   do k = 0, (N*scal)
      do i = 0, (N * scal)
         if ((copy(k,i) .eq. 11) .and. ((copy(modulo(k-1, N*scal),i) .eq. 13) .or. (copy(modulo(k-1, N*scal),i) .eq. 8))) then
            if (copy(k,i+1) .eq. 16) then
               state(k,i) = 3
            else
               state(k,i) = 12
            end if
         else if ((copy(k,i) .eq. 11) .and. (copy(k,i+1) .eq. 16)) then
            state(k,i) = 10
         else if ((copy(k,i) .eq. 8) .and. ((copy(modulo(k+1, N*scal),i) .eq. 13) .or. (copy(modulo(k+1, N*scal),i) .eq. 11))) then
            if (copy(k,i+1) .eq. 14) then
               state(k,i) = 2
            else
               state(k,i) = 9
            end if
         else if ((copy(k,i) .eq. 8) .and. ((copy(k,i+1) .eq. 14))) then
            state(k,i) = 7
         end if
      end do
   end do

   !transferring the labels on another array copying the last column twice to facilitate zero-gradient boundary conditions
   fin_state(0: N * scal, :) = state
   fin_state(N * scal + 1, :) = state(N * scal, :)

   !printing the labeled lattice on file to be used by the LBM simulation program
   open(15, file = "state.dat", status = "replace")
   write(15, *) fin_state
   close(15)

   !printing the labeled lattice on file to be plotted for visualization
   open(12, file = "labeled_lattice.dat", status = "replace")
      do k = 0, N * scal + 1 
         do i = 0, N * scal
            write(12, *) dble(k) + 1, dble(i) + 1, fin_state(k, i)
         end do
      end do
   close(12)
end program Hurst
