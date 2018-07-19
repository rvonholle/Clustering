module MyPoints
   type Point
      real :: x = 0
      real :: y = 0
      integer :: clID = 0
      logical :: isCore = .false.
      logical :: unChecked = .true.
   end type Point
end module MyPoints

module MyClusters
   use MyPoints
   type Cluster
      type(Point), allocatable :: points(:)
   end type Cluster
end module MyClusters

program dbscan_2D
   use MyPoints
   use MyClusters
   implicit none
   integer, parameter :: MAXSIZE = 4000
   integer, parameter :: MINPTS = 4

   type(Point), allocatable :: points(:)
   type(Cluster), allocatable :: clusters(:)
   type(Cluster) :: tempCluster
   character(len=100) :: dataFile, clusteredFile
   integer :: numData, numPoints, numClusters, temp, pointcount, i, j, k
   integer :: n
   type(Point) :: origin
   real :: EPS
   logical :: allPointsChecked = .false.
   logical :: contCluster = .true.

! This program clusters points in a three-dimensional region
   write (*,*)
   write (*,'(A38)') "What file is the data being read from?"
   read *, dataFile
   write (*,'(A42)') "What file will the data be clustered into?"
   read *, clusteredFile

! Read and sort the data points
   numData = 0
   write (*,*)
   write (*,'(A10)') "Reading..."
   open (1, file = dataFile)
   do
      read (1, *, end=10)
      numData = numData + 1
   end do
   10 close(1)
   allocate(points(numData))
   points(:) = origin
   open (1, file = dataFile)
   do i = 1, numData
      read (1,*) points(i)%x, points(i)%y
   end do
   close (1)
   write (*,'(A16)') "Reading complete"
   write (*,*)
   write (*,'(A23)') "Generating parameter..."
   EPS = epsPick(points, MINPTS)
   write (*,'(A19)') "Parameter generated"

! Find the cores
   do i = 1, numData
      temp = 0
      do j = 1, numData
         if (dist(points(i), points(j)) <= EPS) then
            temp = temp + 1
         end if
      end do
      if (temp >= MINPTS) then
         points(i)%isCore = .true.
      end if
   end do

! Cluster the points
   write (*,*)
   write (*,'(A13)') "Clustering..."
   numClusters = 1
   i = 1
   do i = 1, numData
      if (points(i)%isCore .and. points(i)%unChecked .and. points(i)%clID == 0) then
         points(i)%unChecked = .false.
         points(i)%clID = numClusters
         numClusters = numClusters + 1
         call appendPoint(tempCluster%points, points(i))
         do j = 1, numData
            if (points(j)%isCore .and. i /= j .and. points(j)%unChecked .and. &
                  points(j)%clID == 0 .and. dist(points(i), points(j)) <= EPS) then
               points(j)%unchecked = .false.
               points(j)%clID = points(i)%clID
               call appendPoint(tempCluster%points, points(j))
            end if
         end do
         contCluster = .true.
         do while (contCluster)
            temp = size(tempCluster%points)
            do j = 1, temp
               do k = 1, numData
                  if (points(k)%unChecked .and. points(k)%clID==0 .and. &
                        dist(tempCluster%points(j),points(k))<=EPS) then
                     points(k)%unChecked = .false.
                     points(k)%clID = tempCluster%points(j)%clID
                     call appendPoint(tempCluster%points, points(k))
                  end if
               end do
            end do
            if (size(tempCluster%points) == temp) then
               contCluster = .false.
               call appendCluster(clusters, tempCluster)
               deallocate(tempCluster%points)
            end if
         end do
      end if
   end do
   do i = 1, numData
      do j = 1, size(clusters)
         do k = 1, size(clusters(j)%points)
            if (points(i)%unChecked .and. dist(points(i), clusters(j)%points(k)) <= EPS &
                  .and. .not. points(i)%isCore) then
               write(*,*) 'I got here???'
               points(i)%unChecked = .false.
               points(i)%clID = clusters(j)%points(k)%clID
               call appendPoint(clusters(j)%points, points(i))
               exit
            end if
         end do
      end do
   end do
   write (*,'(A19)') "Clustering complete"

! Write the clusters
   write (*,*)
   write (*,'(A10)') "Writing..."
   open (2, file = clusteredFile)
   if (points(1)%x == points(numData)%x .and. points(1)%y == points(numData)%y) then
      write (*,'(A25)') "The data is all the same."
      do i = 1, numData
         if (abs(points(i)%x) < 0.1 .or. abs(points(i)%x) > 1e6) then
            write (2, '(A1,ES13.7,A1)', Advance = 'No') "(", points(i)%x, ","
         else
            write (2, '(A1,F12.7,A1)', Advance = 'No') "(", points(i)%x, ","
         end if
         if (abs(points(i)%y) < 0.1 .or. abs(points(i)%y) > 1e6) then
            write (2, '(ES13.7,A1)') points(i)%y, ")"
         else
            write (2, '(F12.7,A1)') points(i)%y, ")"
         end if
      end do
   else
      do i = 1, size(clusters)
         pointcount = 1
!         write (2,'(A7,I4,A1)') "Cluster", i, ":"
         do j = 1, size(clusters(i)%points)
            if (abs(clusters(i)%points(j)%x) < 0.1 .or. abs(clusters(i)%points(j)%x) > 1e6) then
               write (2, 101, Advance = 'No') pointcount, i, clusters(i)%points(j)%x, ","
            else
               write (2, 102, Advance = 'No') pointcount, i, clusters(i)%points(j)%x, ","
            end if
            if (abs(points(j)%y) < 0.1 .or. abs(points(j)%y) > 1e6) then
               write (2, '(ES13.7)') clusters(i)%points(j)%y
            else
               write (2, '(F12.7)') clusters(i)%points(j)%y
            end if
            pointcount = pointcount + 1
         end do
      end do
!      write (2,'(A9)') "Outliers:"
      pointcount = 1
      do i = 1, numData
         if (points(i)%clID == 0 .and. .not. points(i)%isCore) then
            if (abs(points(i)%x) < 0.1 .or. abs(points(i)%x) > 1e6) then
               write (2, 101, Advance = 'No') pointcount, -1, points(i)%x, ","
            else
               write (2, 102, Advance = 'No') pointcount, -1, points(i)%x, ","
            end if
            if (abs(points(i)%y) < 0.1 .or. abs(points(i)%y) > 1e6) then
               write (2, '(ES13.7)') points(i)%y
            else
               write (2, '(F12.7)') points(i)%y
            end if
            pointcount = pointcount + 1
         end if
      end do
   end if
   write (*,'(A16)') "Writing complete"
   write (*,*)

   101 Format(I5,1X,I5,ES13.7,A1)
   102 Format(I5,1X,I5,F12.7,A1)
contains

subroutine sortReals(realList)
   use MyPoints
   use MyClusters
   implicit none
   real, allocatable :: realList(:)
   real :: temp
   integer :: i, j
   do i = 1, size(realList)
      do j = i, size(realList)
         if (realList(i) > realList(j)) then
            temp = realList(i)
            realList(i) = realList(j)
            realList(j) = temp
         end if
      end do
   end do
end subroutine sortReals

! This subroutine appends a real variable to a list of real variables
subroutine appendReal(list, newReal)
   use MyPoints
   use MyClusters
   implicit none
   integer :: i
   real :: newReal
   real, allocatable :: list(:)
   real, allocatable :: tempList(:)
   if (allocated(list)) then
      allocate(tempList(size(list) + 1))
      do i = 1, size(list)
         tempList(i) = list(i)
      end do
      tempList(size(list) + 1) = newReal
      deallocate(list)
      call move_alloc(tempList, list)
      if (allocated(tempList)) then
         deallocate(tempList)
      end if
   else
      allocate(list(1))
      list(1) = newReal
   end if
end subroutine appendReal

! This subroutine appends a point to a list of points
subroutine appendPoint(list, newPoint)
   use MyPoints
   use MyClusters
   implicit none
   integer :: i
   type(Point) :: newPoint
   type(Point), allocatable :: list(:)
   type(Point), allocatable :: tempList(:)
   if (allocated(list)) then
      allocate(tempList(size(list) + 1))
      do i = 1, size(list)
         tempList(i) = list(i)
      end do
      tempList(size(list) + 1) = newPoint
      deallocate(list)
      call move_alloc(tempList, list)
      if (allocated(tempList)) then
         deallocate(tempList)
      end if
   else
      allocate(list(1))
      list(1) = newPoint
   end if
end subroutine appendPoint

! This subroutine appends a cluster to a list of clusters
subroutine appendCluster(list, newCluster)
   use MyPoints
   use MyClusters
   implicit none
   integer :: i
   type(Cluster) :: newCluster
   type(Cluster), allocatable :: list(:)
   type(Cluster), allocatable :: tempList(:)
   if (allocated(list)) then
      allocate(tempList(size(list) + 1))
      do i = 1, size(list)
         tempList(i) = list(i)
      end do
      tempList(size(list) + 1) = newCluster
      deallocate(list)
      call move_alloc(tempList, list)
      if (allocated(tempList)) then
         deallocate(tempList)
      end if
   else
      allocate(list(1))
      list(1) = newCluster
   end if
end subroutine appendCluster

! This function finds the square of the distance between two points
function dist(pointA, pointB)
   use MyPoints
   use MyClusters
   implicit none
   type(Point), intent(in) :: pointA, pointB
   real :: dist
   dist = ((pointB%x - pointA%x) ** 2)
   dist = dist + ((pointB%y - pointA%y) ** 2)
end function dist

! This function determines the optimal value of epsilon
function epsPick(points, MINPTS)
   use MyPoints
   use MyClusters
   type(Point), allocatable :: points(:)
   integer :: MINPTS
   real :: epsPick
   real, allocatable :: tempDists(:)
   real, allocatable :: dists(:)
   real, allocatable :: distsPrime(:)
   real, allocatable :: distsPrime2(:)
   do i = 1, size(points)
      do j = 1, size(points)
         if (i /= j) then
            call appendReal(tempDists, sqrt(dist(points(i), points(j))))
         end if
      end do
      call sortReals(tempDists)
      call appendReal(dists, tempDists(MINPTS))
      deallocate(tempDists)
   end do
   do i = 2, size(dists) - 1
      call appendReal(distSPrime, (dists(i + 1) - dists(i - 1)) / 2)
   end do
   do i = 2, size(distsPrime) - 1
      call appendReal(distsPrime2, (distsPrime(i + 1) - distsPrime(i - 1)) / 2)
   end do
   do i = 1, size(distsPrime2)
      distsPrime2(i) = abs(distsPrime2(i))
   end do
   epsPick = dists(maxloc(distsPrime2, 1) + 2)
end function epsPick

end program dbscan_2D
