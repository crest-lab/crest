!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Stefan Grimme
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

subroutine fragmentize(nspin,maxsystem, maxmagnat, jab, neigh, ispinsyst, nspinsyst, nsystem)
    implicit none
    !Dummy Arguments:
    integer, intent(in)  :: nspin
    integer, intent(in)  :: maxsystem
    integer, intent(in)  :: maxmagnat
    integer, intent(in)  :: neigh(200,nspin)
    real*8,  intent(in)  :: jab(nspin*(nspin+1)/2)
    integer, intent(out) :: ispinsyst  (10*nspin,maxsystem)
    integer, intent(out) :: nspinsyst  (maxsystem)
    integer, intent(out) :: nsystem

    !Stack
    integer:: i
    integer:: current
    integer:: ass
    integer:: k
    integer:: j
    integer:: lin
    integer:: maxdistatoms(2)
    integer:: max_linkatoms(2)
    integer:: assigned_to_frag(nspin)
    integer:: precessor(nspin)
    real*8:: magdist(nspin,nspin)
    real*8:: rmaxab(nspin,nspin)
    real*8:: shortest_distance
    real*8:: maxdist
    real*8:: cur_dist
    real*8:: max_link
    logical:: visited(nspin)
    logical:: assigned(nspin,nspin)

    real*8, parameter :: thr_j=0.001d0

!    write(*,*) nspin, maxsystem, nsystem
!      write( *, * ) "Coupling-Matrix"
       do i = 1, nspin
!           write(*,'(16(f6.3,1x))') (abs(jab(lin(j,i))),j=1,nspin)
       end do ! End Loop over  i from 1 to nspin



    magdist = huge(1.0d0)
!    write( *, * ) "inverse couplings"
    do i=1,nspin
        do j=1,nspin
            if(i.eq.j.or.abs(jab(lin(j,i))).lt. thr_j) cycle
            !magdist(j,i) = maxab - abs(jab(lin(j,i)))
            magdist(j,i) = 1 / abs(jab(lin(j,i)))
        enddo
!        write(*,'(16(f6.3,1x))') (magdist(j,i),j=1,nspin)
    enddo

    rmaxab = 0
!    write( *, * ) "shortest distances"
    do i = 1, nspin
        do j = 1, nspin
            rmaxab(j,i) = shortest_distance(nspin, i, j, neigh, magdist, visited, precessor)
        end do ! End Loop over  j from 1 to nspin
!        write(*,'(16(f6.3,1x))') (rmaxab(j,i),j=1,nspin)
    end do ! End Loop over  i from 1 to nspin


    ispinsyst = 0
    nspinsyst = 0
    
    assigned_to_frag = 0
    nsystem = 0
    !If spin systems are already separated, keep them:
    if (maxval(rmaxab, mask=assigned) == huge(1.0d0)) then
	assigned = .true.
	do while (maxval(rmaxab, mask=assigned) == huge(1.0d0))
	nsystem = nsystem + 1
	maxdistatoms= maxloc(rmaxab, mask=assigned)
	cur_dist = shortest_distance(nspin, maxdistatoms(1), maxdistatoms(2), neigh, magdist, visited, precessor)
	nspinsyst(nsystem) = count(visited)
	    ispinsyst(:,nsystem) = 0
	    k = 1
	    do i = 1, nspin
		if (visited(i) .eqv. .true.) then
		    ispinsyst(k,nsystem) = i
		    assigned(i,:) = .false.
		    k = k + 1
		end if
	    end do ! End Loop over  i from 1 to count(visited)     
	end do
    else
    !else put everything in one big spinsystem
    nsystem = 1
    nspinsyst(1) = nspin
    do i = 1, nspin
        ispinsyst(i,1) = i
    end do ! End Loop over  i from 1 to nspin
    end if

    

    !Find Sub-Systems with nspinsyst > maxmagnat
    ass = maxloc(nspinsyst,1)
    do while (nspinsyst(ass) > maxmagnat)


        !Set Spins in this list as available for algorithm
        assigned = .false.        
        do i = 1, nspinsyst(ass)
            do j = 1, nspinsyst(ass)
                assigned(ispinsyst(i,ass),ispinsyst(j,ass)) = .true.
            end do ! End Loop over  j from 1 to nspinsyst(ass)
        end do ! End Loop over  i from 1 to nspin

!               do i = 1, nspin
!                   write(*,'(16(L3))') (assigned(j,i),j=1,nspin)
!               end do ! End Loop over  i from 1 to nspin
!        
!                write( *, * ) ""
!                do i = 1, nspin
!                    write(*,'(16(f6.3,1x))') (rmaxab(j,i),j=1,nspin)
!                end do ! End Loop over  i from 1 to nspin


        !Find largest distance: 
        maxdist=maxval(rmaxab, mask=assigned)
        maxdistatoms= maxloc(rmaxab, mask=assigned)
        !        write( *, * ) "maxdistatoms", maxdistatoms

        !If a Path is found between A and B
        if (maxdist < huge(1.0d0)) then

            !get shortest Path from A to B
            cur_dist=shortest_distance(nspin, maxdistatoms(1), maxdistatoms(2), neigh, magdist, visited, precessor)
            current = maxdistatoms(2)

            !loop while A and B are still connected
            do while (cur_dist < huge(1.0d0))

                !find weakest link
                max_link = 0
                max_linkatoms= 0
                do while (precessor(current) /= 0)
                    if (magdist(current, precessor(current)) > max_link) then
                        max_link = magdist(current, precessor(current))
                        max_linkatoms(1) = current
                        max_linkatoms(2) = precessor(current)
                    end if
                    current = precessor(current)
                end do ! End loop: while precessor(current) /= 0


                !Split weakest link, set distance to infinity
                magdist(max_linkatoms(1), max_linkatoms(2)) = huge(1.0d0)
                magdist(max_linkatoms(2), max_linkatoms(1)) = huge(1.0d0)

                !Get next-shortest Path:
                cur_dist=shortest_distance(nspin, maxdistatoms(1), maxdistatoms(2), neigh, magdist, visited, precessor)
                current = maxdistatoms(2)
            end do ! cur_dist < huge(1.0d0

            !A and B are now Seperated:
            rmaxab(maxdistatoms(1),maxdistatoms(2))=huge(1.0d0)
            rmaxab(maxdistatoms(2),maxdistatoms(1))=huge(1.0d0)
        end if ! End if: while maxdist < huge(1.0d0)

        !Split into subsystems:
        !Overwrite old spinsystem:
        
        cur_dist = shortest_distance(nspin, maxdistatoms(1), maxdistatoms(2), neigh, magdist, visited, precessor)
        nspinsyst(ass) = count(visited)
        ispinsyst(:,ass) = 0
        k = 1
        do i = 1, nspin
            if (visited(i)) then
                ispinsyst(k,ass) = i
                k = k + 1
            end if
        end do ! End Loop over  i from 1 to count(visited)

        !add new spinsystem
        cur_dist = shortest_distance(nspin, maxdistatoms(2), maxdistatoms(1), neigh, magdist, visited, precessor)
        nsystem = nsystem + 1
        nspinsyst(nsystem) = count(visited)
        k = 1
        do i = 1, nspin
            if (visited(i)) then
                ispinsyst(k,nsystem) = i
                k = k + 1
            end if
        end do ! End Loop over  i from 1 to count(visited)
        ass = maxloc(nspinsyst,1)

    !        write(*,*)'spinsystem    nuclei (after iteration)'
    !        do i=1,nsystem
    !            write(*,'(i5,5x,20i4)')i,(ispinsyst(j,i),j=1,nspinsyst(i))
    !        end do
    end do ! End loop: while nspinsyst(ass) > maxmagnat

!   rmaxab = 0
!   write( *, * ) "Final map"
!   do i = 1, nspin
!       do j = 1, nspin
!           rmaxab(j,i) = shortest_distance(nspin, i, j, neigh, magdist, visited, precessor)
!       end do ! End Loop over  j from 1 to nspin
!       write(*,'(16(f6.3,1x))') (rmaxab(j,i),j=1,nspin)
!   end do ! End Loop over  i from 1 to nspin




end subroutine 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculates the shortest path between two atoms
!start and goal are integers, determining the index in xyz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function shortest_distance(nspin, start, goal, neighbours, input_distances, visited, precessor)
    implicit none
    !Dummy Arguments:
    integer, intent(in) :: nspin
    integer, intent(in) :: start
    integer, intent(in) :: goal
    integer, intent(in) :: neighbours(200,nspin)
    real*8,  intent(in) :: input_distances(nspin,nspin)
    logical, intent(out):: visited(nspin)
    integer, intent(out):: precessor(nspin)

    !Stack:
    integer:: current
    integer:: neighbour
    integer:: i_neighbours
    real*8:: alt_dist
    integer:: bonds
    real*8:: distance(nspin)

    !logical:: visited(nspin)

    bonds = 0
    precessor = 0
    distance = huge(distance)
    distance(start) = 0
    visited = .false.
    do while (all(visited) .eqv. .false.) !as long there are unvisited nodes
        current = minloc(distance,1, .not. visited)

        !Abort if Fragments are not connected:
        if (distance(current) == huge(distance)) then
            shortest_distance = huge(distance)
            return
        end if ! distance(current) == huge(distance)

        visited(current) = .true.
        if (current == goal) then
            shortest_distance = distance(goal)
            !route to target found
            do while (precessor(current) /= 0)
                bonds = bonds + 1
                current = precessor(current)
            end do ! End loop: while precessor(current) /= 0
            exit
        else
            !loop over all neighbours of current atom
            do i_neighbours = 1, neighbours(200,current)
                neighbour = neighbours(i_neighbours,current)
                if (visited(neighbour) .eqv. .false.) then
                    !distanzupdate
                    alt_dist = distance(current) + input_distances(current,neighbour)
                    !write( *, * ) alt_dist, distance(current),
                    if ((alt_dist < distance(neighbour))) then
                        distance(neighbour) = alt_dist
                        precessor(neighbour) = current
                    end if ! (alt_dist < distance(neighbour))
                endif !(visited(neighbour))
            end do ! End Loop over  i_neighbours from 1 to all_Atoms(current)%n_neighbours
        endif
    end do ! End loop: while sum(visited) /= 0
        !initialize
    shortest_distance = distance(goal)
end function shortest_distance
