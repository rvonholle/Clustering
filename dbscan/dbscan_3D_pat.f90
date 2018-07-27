module MyPoints
   type Point
      real :: x = 0
      real :: y = 0
      real :: z = 0
      integer :: clID = 0
      integer :: orig_pos = 0
      integer :: nodeID = 0
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

module FileStorage
   use MyPoints
   use MyClusters
   type titleCard
      character(len=133) :: line1
      character(len=100) :: title
   end type titleCard
   type summaryData
      integer :: item(4)
      integer :: numNodes
      integer :: numElems
      integer :: item2
      integer :: numClusters
      integer :: item3
      character(len=133) line2
   end type summaryData
   type packet2
      integer :: numNodes
      integer :: item
      integer :: clID
      integer :: item2
      real :: item3(3)
      integer, allocatable :: nodes(:)
   end type packet2
   type packet4
      character(len=2) :: packID = "04"
      integer :: clustNum
      integer :: item(7)
      character(len=133) line2
   end type packet4
   type InputFile
      type(titleCard) :: title
      type(summaryData) :: summary
      type(packet2), allocatable :: elems(:)
      type(packet4) :: clusts
      integer :: packet3Lines
   end type InputFile
end module FileStorage      

program dbscan_2D
   use MyPoints
   use MyClusters
   use FileStorage
   implicit none
   integer, parameter :: MAXSIZE = 4000
   integer, parameter :: MINPTS = 8
   integer, parameter :: FILENUMBER = 72

   type(Point), allocatable :: points(:)
   type(Cluster), allocatable :: clusters(:)
   type(Cluster) :: tempCluster, tempCluster_2
   type(InputFile) :: origFile
   real, dimension(:,:), allocatable :: dists(:,:)
   real, allocatable :: tempDists(:)
   character(len=100) :: dataFile, clusteredFile
   character(len=300) :: copyCommand
   integer :: numData, numCores, numClusters, numChecked, temp, pointcount, i, j, k
   integer :: n
   type(Point) :: origin
   real :: EPS
   logical :: contCluster = .true., inCluster = .true.

! This program clusters points in a three-dimensional region
   write (*,*)
   write (*,'(A38)') "What file is the data being read from?"
   read *, dataFile
   write (*,'(A42)') "What file will the data be clustered into?"
   read *, clusteredFile

! Read and sort the data points
   write (*,*)
   write (*,'(A10)') "Reading..."
   call readPatFile(dataFile, origFile, points)
   numData = size(points)
   write (*,'(A16)') "Reading complete"
   write (*,*)
   write (*,'(A31)') "Generating EPSILON parameter..."
   EPS = epsPick(points, MINPTS, dists)
   write (*,'(A27)') "EPSILON parameter generated"

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
         if (temp >= MINPTS) then
            exit
         end if
      end do
      if (temp >= MINPTS) then
         points(i)%isCore = .true.
         numCores = numCores + 1
      end if
   end do
   write (*,'(A11)') "Cores found"
!   write (*,*) numCores

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

! Cluster the elements
   do i = 1, origFile%summary%numElems
      do j = 1, size(clusters)
         inCluster = .true.
         do k = 1, size(origFile%elems(i)%nodes)
            if (.not. any(clusters(j)%points(:)%nodeID == origFile%elems(i)%nodes(k))) then
               inCluster = .false.
            end if
         end do
         if (inCluster) then
            origFile%elems(i)%clID = j
         else
            do k = 1, size(origFile%elems(i)%nodes)
               if (any(clusters(j)%points(:)%nodeID == origFile%elems(i)%nodes(k))) then
                  inCluster = .true.
               end if
            end do
            if (inCluster) then
               origFile%elems(i)%clID = j
            end if
         end if
      end do
   end do

! Write the clusters
   write (*,*)
   write (*,'(A10)') "Writing..."
   call writePatFile(dataFile, clusteredFile, FILENUMBER, clusters, origFile)
   write (*,'(A16)') "Writing complete"
   write (*,*)

contains

! This subroutine sorts a list of real numbers
recursive subroutine sortReals(list, low, high)
   use MyPoints
   use MyClusters
   use FileStorage
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
   use FileStorage
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
   use FileStorage
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
   use FileStorage
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
   use FileStorage
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
   use FileStorage
   type(Point), allocatable :: points(:)
   integer :: MINPTS
   real :: epsPick
   real, dimension(:,:), allocatable :: allDists(:,:)
   real, allocatable :: tempDists(:)
   real, allocatable :: dists(:)
   real, allocatable :: comps(:)
   integer :: temp, i, j
   real :: lineMaker
   allocate(dists(size(points)))
   allocate(allDists(size(points),size(points)))
   write (*,'(A25)') "Creating distance list..."
   temp = floor(25 * size(points) / 100.0)
   do i = 1, temp
      do j = 1, size(points)
         allDists(j,i) = dist(points(j), points(i))
      end do
   end do
   write (*,'(A26)') "Distance list 25% complete"
   temp = floor(50 * size(points) / 100.0)
   do i = floor(25 * size(points) / 100.0) + 1, temp
      do j = 1, size(points)
         allDists(j,i) = dist(points(j), points(i))
      end do
   end do
   write (*,'(A26)') "Distance list 50% complete"
   temp = floor(75 * size(points) / 100.0)
   do i = floor(50 * size(points) / 100.0), temp
      do j = 1, size(points)
         allDists(j,i) = dist(points(j), points(i))
      end do
   end do
   write (*,'(A26)') "Distance list 75% complete"
   do i = temp + 1, size(points)
      do j = 1, size(points)
         allDists(j,i) = dist(points(j), points(i))
      end do
   end do
   write (*,'(A22)') "Distance list complete"
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

! This subroutine reads from the input file
subroutine readPatFile(fileName, patFile, points)
   use MyPoints
   use MyClusters
   use FileStorage
   implicit none
   character(len=100) :: fileName
   type(InputFile) :: patFile
   type(Point), allocatable :: points(:)
   character(len=5) :: skipItem
   integer :: i, j, stat
   open(1,file = fileName)
! Read header
   read (1,'(A)') patFile%title%line1
   read (1,'(A)') patFile%title%title
   read (1,*) patFile%summary%item(1), patFile%summary%item(2), patFile%summary%item(3), &
      patFile%summary%item(4), patFile%summary%numNodes, patFile%summary%numElems, &
      patFile%summary%item2, patFile%summary%numClusters, patFile%summary%item3
   read (1,'(A)') patFile%summary%line2
   allocate(patFile%elems(patFile%summary%numElems))
   allocate(points(patFile%summary%numNodes))
! Read nodes
   do i = 1, patFile%summary%numNodes
      read (1,*) skipItem, points(i)%nodeID
      read (1,*) points(i)%x, points(i)%y, points(i)%z
      read (1,*)
      points(i)%orig_pos = i
   end do
! Read elements
   do i = 1, patFile%summary%numElems
      read (1,*)
      read (1,*) patFile%elems(i)%numNodes, patFile%elems(i)%item, patFile%elems(i)%clID, &
         patFile%elems(i)%item2, patFile%elems(i)%item3(1), patFile%elems(i)%item3(2), patFile%elems(i)%item3(3)
      allocate(patFile%elems(i)%nodes(patFile%elems(i)%numNodes))
      do j = 1, patFile%elems(i)%numNodes - 1
         read (1,'(I8)',Advance='No',iostat=stat) patFile%elems(i)%nodes(j)
      end do
      read (1,*) patFile%elems(i)%nodes(patFile%elems(i)%numNodes)
   end do
! Read packet 3
   read (1,iostat=stat) skipItem, skipItem, skipItem, patFile%packet3Lines
   if (stat == 0) then
      do i = 1, patFile%packet3Lines
         read (1,*)
      end do
! Read packet 4s
      read (1,*) patFile%clusts%packID, skipItem, patFile%clusts%item(1), &
         patFile%clusts%item(2), patFile%clusts%item(3), patFile%clusts%item(4), &
         patFile%clusts%item(5), patFile%clusts%item(6), patFile%clusts%item(7)
      read (1,'(A)') patFile%clusts%line2
   else
      patFile%clusts%item(1) = 1
      patFile%clusts%item(2) = 1
      patFile%clusts%item(3) = 8
      patFile%clusts%item(4) = 8
      patFile%clusts%item(5) = 0
      patFile%clusts%item(6) = 1
      patFile%clusts%item(7) = 0
      patFile%clusts%line2 = "  1.00000000e+00"
   end if
   close(1)
end subroutine readPatFile

! This subroutine writes to the output file
subroutine writePatFile(origFile, outputFile, FILENUMBER, clusters, patFile)
   use MyPoints
   use MyClusters
   use FileStorage
   implicit none
   character(len=100) :: origFile, outputFile
   integer :: FILENUMBER
   type(Cluster), allocatable :: clusters(:)
   type(InputFile) :: patFile
   character(len=133) :: tempLine
   integer :: numLines = 0, linesWritten = 0
   open(1, file = origFile)
   do
      read(1,*,end=10)
      numLines = numLines + 1
   end do
   10 close(1)
   open(1, file = origFile)
   open(FILENUMBER, file=outputFile)
! Write header
   write(FILENUMBER,'(A133)') patFile%title%line1
   write(FILENUMBER,'(A100)') outputFile
   write(FILENUMBER,'(I2,I8,I8,I8,I8,I8,I8,I8,I8)') patFile%summary%item(1), patFile%summary%item(2), &
      patFile%summary%item(3), patFile%summary%item(4), patFile%summary%numNodes, &
      patFile%summary%numElems, patFile%summary%item2, size(clusters), patFile%summary%item3
   write(FILENUMBER,'(A133)') patFile%summary%line2
   linesWritten = 4
   do i = 1, linesWritten
      read (1,*)
   end do
! Write nodes
   do i = 1, patFile%summary%numNodes * 3
      read (1,'(A)') tempLine
      write (FILENUMBER,'(A133)') tempLine
      linesWritten = linesWritten + 1
   end do
! Write elements
   do i = 1, size(patFile%elems)
      read (1,'(A)') tempLine
      write (FILENUMBER,'(A133)') tempLine
      linesWritten = linesWritten + 1
      write (FILENUMBER,'(I8,I8,I8,I8,ES16.8,ES16.8,ES16.8)') patFile%elems(i)%numNodes, &
         patFile%elems(i)%item, patFile%elems(i)%clID, patFile%elems(i)%item2, 0.0, 0.0, 0.0
      linesWritten = linesWritten + 1
      do j = 1, patFile%elems(i)%numNodes - 1
         write (FILENUMBER,'(I8)',Advance='No') patFile%elems(i)%nodes(j)
      end do
      write (FILENUMBER,'(I8)') patFile%elems(i)%nodes(j)
      linesWritten = linesWritten + 1
      read (1,*)
      read (1,*)
   end do
! Write packet 3
   if (patFile%summary%item2 > 0) then
      read (1,'(A)') tempLine
      write (FILENUMBER,'(A133)') tempLine
      linesWritten = linesWritten + 1
      do i = 1, patFile%packet3Lines
         read (1,'(A)') tempLine
         write (FILENUMBER,'(A133)') tempLine
         linesWritten = linesWritten + 1
      end do
   end if
! Write packet 4s
   if (linesWritten < numLines - 1) then
      do i = 1, patFile%summary%numClusters
         read (1,*)
         read (1,*)
         linesWritten = linesWritten + 2
      end do
   end if
   do i = 1, size(clusters)
      write (FILENUMBER,'(A2,I8,I8,I8,I8,I8,I8,I8,I8)') patFile%clusts%packID, i, &
         patFile%clusts%item(1), patFile%clusts%item(2), patFile%clusts%item(3), &
         patFile%clusts%item(4), patFile%clusts%item(5), patFile%clusts%item(6), patFile%clusts%item(7)
      write (FILENUMBER,'(A)') patFile%clusts%line2
   end do
! Write remainder of file
   if (linesWritten < numLines) then
      do i = linesWritten, numLines
         read (1,'(A)', end=20) tempLine
         write (FILENUMBER,'(A133)') tempLine
      end do
20    close(1)
   else if (linesWritten >= numLines) then
      write (FILENUMBER,'(I2,I8,I8,I8,I8,I8,I8,I8,I8)') 99, 0, 0, 1, 0, 0, 0, 0, 0
   end if
   close(FILENUMBER)
end subroutine writePatFile

end program dbscan_2D
