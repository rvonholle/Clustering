module MyPoints
   type Point
      real :: x = 0
      real :: y = 0
      real :: z = 0
      integer :: clID = 0
      integer :: orig_pos = 0
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
   integer, parameter :: MINPTS = 5

   type(Point), allocatable :: points(:)
   type(Cluster), allocatable :: clusters(:)
   type(Cluster) :: tempCluster, tempCluster_2
   real, dimension(:,:), allocatable :: dists(:,:)
   real, allocatable :: tempDists(:)
   character(len=100) :: dataFile, clusteredFile
   integer :: numData, numCores, numClusters, numChecked, temp, pointcount, i, j, k
   integer :: n
   type(Point) :: origin
   real :: EPS
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
      read (1,*) points(i)%x, points(i)%y, points(i)%z
      points(i)%orig_pos = i
   end do
   close (1)
   write (*,'(A16)') "Reading complete"
   write (*,*)
   write (*,'(A31)') "Generating EPSILON parameter..."
   EPS = epsPick(points, MINPTS, dists)
   write (*,'(A27)') "EPSILON parameter generated"
!   write (*,*) EPS

! Find the cores
   numCores = 0
   write (*,*)
   write (*,'(A16)') "Finding cores..."
   do i = 1, numData
      temp = 0
      do j = 1, numData
         if (dists(j,i) <= EPS) then
            temp = temp + 1
         end if
      end do
      if (temp >= MINPTS) then
         points(i)%isCore = .true.
         numCores = numCores + 1
      end if
   end do
   write (*,'(A11)') "Cores found"

! Cluster the points
   write (*,*)
   write (*,'(A13)') "Clustering..."
   numClusters = 1
   numChecked = 0
   do i = 1, numData
      if (points(i)%isCore .and. points(i)%unChecked) then
         points(i)%unChecked = .false.
         numChecked = numChecked + 1
         points(i)%clID = numClusters
         numClusters = numClusters + 1
         call appendPoint(tempCluster%points, points(i))
         do j = 1, i - 1
            if (points(j)%isCore .and. points(j)%unChecked .and. &
                  dists(j,i) <= EPS .and. size(tempCluster%points) < MAXSIZE) then
               points(j)%unchecked = .false.
               numChecked = numChecked + 1
               points(j)%clID = points(i)%clID
               call appendPoint(tempCluster%points, points(j))
            else if (size(tempCluster%points) >= MAXSIZE) then
               exit
            end if
         end do
         do j = i + 1, numData
            if (points(j)%isCore .and. points(j)%unChecked .and. &
                  dists(j,i) <= EPS .and. size(tempCluster%points) < MAXSIZE) then
               points(j)%unchecked = .false.
               numChecked = numChecked + 1
               points(j)%clID = points(i)%clID
               call appendPoint(tempCluster%points, points(j))
            else if (size(tempCluster%points) >= MAXSIZE) then
               exit
            end if
         end do
         contCluster = .true.
         do while (contCluster)
            temp = size(tempCluster%points)
            do j = 1, temp
               do k = 1, numData
                  if (points(k)%unChecked .and. dists(k,tempCluster%points(j)%orig_pos) <= EPS  &
                         .and. size(tempCluster%points) < MAXSIZE) then
                     points(k)%unChecked = .false.
                     numChecked = numChecked + 1
                     points(k)%clID = tempCluster%points(j)%clID
                     call appendPoint(tempCluster%points, points(k))
                  else if (size(tempCluster%points) >= MAXSIZE) then
                     exit
                  end if
               end do
            end do
            if (size(tempCluster%points) == temp .and. temp >= MINPTS) then
               contCluster = .false.
               call appendCluster(clusters, tempCluster)
               deallocate(tempCluster%points)
               write (*,'(A10,I4,A10)') "Clustering", nint(100.0 * numChecked / size(points)), "% complete"
            else if (size(tempCluster%points) == temp) then
               contCluster = .false.
               do j = 1, temp
                  points(j)%isCore = .false.
                  points(j)%clID = 0
               end do
               deallocate(tempCluster%points)
            end if
         end do
      end if
      if (numChecked == numCores) then
         exit
      end if
   end do
   do i = 1, numData
      do j = 1, size(clusters)
         do k = 1, size(clusters(j)%points)
            if (points(i)%unChecked .and. dists(clusters(j)%points(k)%orig_pos,i) <= EPS &
                  .and. .not. points(i)%isCore) then
               points(i)%unChecked = .false.
               points(i)%clID = clusters(j)%points(k)%clID
               call appendPoint(clusters(j)%points, points(i))
               exit
            end if
         end do
      end do
   end do
   write (*,'(A19)') "Clustering complete"

! Add the outliers
   do i = 1, numData
      if (points(i)%clID == 0) then
         allocate(tempDists(numData))
         do j = 1, numData
            tempDists(j) = dists(j,i)
         end do
         do while (points(minloc(tempDists,1))%clID == 0)
            tempDists(minloc(tempDists,1)) = maxval(tempDists,1)
         end do
         points(i)%clID = points(minloc(tempDists,1))%clID
         call appendPoint(clusters(points(i)%clID)%points, points(i))
         deallocate(tempDists)
      end if
   end do

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
            write (2, 101, Advance = 'No') pointcount, i, clusters(i)%points(j)%x
            write (2, '(ES17.7,3X)', Advance = 'No') clusters(i)%points(j)%y
            write (2, '(ES17.7)') clusters(i)%points(j)%z
            pointcount = pointcount + 1
         end do
      end do
!      write (2,'(A9)') "Outliers:"
!      pointcount = 1
!      do i = 1, numData
!         if (points(i)%clID == 0 .and. .not. points(i)%isCore) then
!            if (abs(points(i)%x) < 0.1 .or. abs(points(i)%x) > 1e6) then
!               write (2, 101, Advance = 'No') pointcount, 0, points(i)%x, ","
!            else
!               write (2, 102, Advance = 'No') pointcount, 0, points(i)%x, ","
!            end if
!            if (abs(points(i)%y) < 0.1 .or. abs(points(i)%y) > 1e6) then
!               write (2, '(ES13.7)') points(i)%y
!            else
!               write (2, '(F12.7)') points(i)%y
!            end if
!            pointcount = pointcount + 1
!         end if
!      end do
   end if
   write (*,'(A16)') "Writing complete"
   write (*,*)

   101 Format(I5,1X,I5,3X,ES17.7,3X)
   102 Format(I5,1X,I5,3X,F12.7,3X)

contains

! This subroutine sorts a list of real numbers
recursive subroutine sortReals(list, low, high)
   use MyPoints
   use MyClusters
   implicit none
   real, allocatable :: list(:)
   integer :: low, high
   integer :: partID = 0
   if (low < high) then
      call partition(list, low, high, partID)
      call sortReals(list, low, partID - 1)
      call sortReals(list, partID + 1, high)
   end if
end subroutine sortReals

! This subroutine partitions a list of real numbers
subroutine partition(list, low, high, partID)
   use MyPoints
   use MyClusters
   implicit none
   real, allocatable :: list(:)
   integer :: low, high, partID
   integer :: i, j
   real :: pivot, temp
   pivot = list(high)
   i = (low - 1)
   do j = low, high - 1
      if (list(j) <= pivot) then
         i = i + 1
         temp = list(i)
         list(i) = list(j)
         list(j) = temp
      end if
   end do
   temp = list(i + 1)
   list(i + 1) = list(high)
   list(high) = temp
   partID = i + 1
end subroutine partition

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

! This function finds the distance between two points
function dist(pointA, pointB)
   use MyPoints
   use MyClusters
   implicit none
   type(Point), intent(in) :: pointA, pointB
   real :: dist
   dist = ((pointB%x - pointA%x) ** 2)
   dist = dist + ((pointB%y - pointA%y) ** 2)
   dist = dist + ((pointB%z - pointA%z) ** 2)
   dist = sqrt(dist)
end function dist

! This function determines the optimal value of epsilon
function epsPick(points, MINPTS, allDists)
   use MyPoints
   use MyClusters
   type(Point), allocatable :: points(:)
   integer :: MINPTS
   real :: epsPick
   real, dimension(:,:), allocatable :: allDists(:,:)
   real, allocatable :: tempDists(:)
   real, allocatable :: dists(:)
   real, allocatable :: comps(:)
   integer :: i, j
   real :: lineMaker
   allocate(dists(size(points)))
   allocate(allDists(size(points),size(points)))
   write (*,'(A25)') "Creating distance list..."
   do i = 1, size(points)
      do j = 1, size(points)
         allDists(j,i) = dist(points(j), points(i))
      end do
      if (floor(100.0 * i / size(points)) == 25 .and. &
         floor(100.0 * (i - 1) / size(points)) /= 25) then
         write (*,'(A26)') "Distance list 25% complete"
      else if (floor(100.0 * i / size(points)) == 50 .and. &
         floor(100.0 * (i - 1) / size(points)) /= 50) then
         write (*,'(A26)') "Distance list 50% complete"
      else if (floor(100.0 * i / size(points)) == 75 .and. &
         floor(100.0 * (i - 1) / size(points)) /= 75) then
         write (*,'(A26)') "Distance list 75% complete"
      end if
   end do
   do i = 1, size(points)
      allocate(tempDists(size(points)))
      do j = 1, size(points)
         tempDists(j) = allDists(j,i)
      end do
      do j = 1, MINPTS
         tempDists(minloc(tempDists)) = maxval(tempDists)
      end do
      dists(i) = minval(tempDists)
      deallocate(tempDists)
   end do
   write (*,'(A22)') "Distance list complete"
   write (*,'(A24)') "Sorting distance list..."
   call sortReals(dists, 1, size(dists))
   write (*,'(A20)') "Distance list sorted"
   lineMaker = (dists(size(dists)) - dists(1)) / size(dists)
   allocate(comps(size(dists)))
   do i = 1, size(dists)
      comps(i) = (dists(1) + (i - 1) * lineMaker) - dists(i)
   end do
   epsPick = dists(maxloc(comps, 1))
end function epsPick

end program dbscan_2D
