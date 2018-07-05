program k_means_clustering
   implicit none
   type Point
      real :: x_coord
      real :: y_coord
   end type Point
   type Cluster
      type(Point) :: center
      type(Point), allocatable :: clusterPoints(:)
   end type Cluster
   type(Point), allocatable :: dataPoints(:)
   type(Point), dimension(:,:), allocatable :: clusters(:,:)
   type(Point), allocatable :: tempList(:)
   type(Point) :: origin
   character(len=100) :: dataFile, clusteredFile
   character(len=1) :: chooseClusters
   real :: numData_copy, numClusters_copy
   integer :: numData, numClusters, temp, i, j, k, m, n
   origin%x_coord = 0
   origin%y_coord = 0

! This program clusters points on a two-dimensional plane
   write (*,*) "What file is the data being read from?"
   read *, dataFile
   write (*,*) "What file will the data be clustered into?"
   read *, clusteredFile
   write (*,*) "Would you like to choose the number of clusters (y/n)?"
   read *, chooseClusters
   do while (chooseClusters /= 'y' .and. chooseClusters /= 'n')
      write (*,*) "Please enter y (yes) or n (no)."
      write (*,*) "Would you like to choose the number of clusters (y/n)?"
      read *, chooseClusters
   end do
   if (chooseClusters == 'y') then
      write (*,*) "How many clusters will there be?"
      read *, numClusters
   end if

! Read and sort the data points
   numData = 0
   open (1, file = dataFile)
   do
      read (1, *, end=10)
      numData = numData + 1
   end do
   10 close (1)
   allocate(dataPoints(numData))
   dataPoints(:)%x_coord = 0
   dataPoints(:)%y_coord = 0
   open (1, file = dataFile)
   do i = 1, numData
      read (1,*) dataPoints(i)%x_coord, dataPoints(i)%y_coord
   end do
   close (1)
   call sortPoints(dataPoints, origin, numData)

! Determine the number of clusters
   numData_copy = numData
   numClusters_copy = numClusters
   if (chooseClusters == 'n') then
      numClusters_copy = sqrt(numData_copy / 2)
      numClusters = nint(numClusters_copy)
   else
      do while (numClusters >= numData .and. chooseClusters /= 'n')
         do while (numClusters > numData)
            write (*,*) "There are more clusters than there is data."
            write (*,*) "Please choose a smaller number of clusters."
            read *, numClusters
         end do
         do while (numClusters == numData .and. chooseClusters /= 'n')
            write (*,*) "Each data point will be in its own cluster."
            write (*,*) "This can be avoided if the number of clusters is smaller."
            write (*,*) "Would you like to choose a smaller number of clusters (y/n)?"
            read *, chooseClusters
            do while (chooseClusters /= 'y' .and. chooseClusters /= 'n')
               write (*,*) "Please enter y (yes) or n (no)."
               read *, chooseClusters
            end do
            if (chooseClusters == 'y') then
               write (*,*) "How many clusters will there be?"
               read *, numClusters
            end if
         end do
      end do
   end if
   allocate(clusters(numClusters+1,numData+1))
   allocate(tempList(numData+1))
   clusters(:,:)%x_coord = 0
   clusters(:,:)%y_coord = 0
   tempList(:)%x_coord = 0
   tempList(:)%y_coord = 0

! Cluster the points
   if(dataPoints(1)%x_coord/=dataPoints(numData)%x_coord.or.dataPoints(1)%y_coord/=dataPoints(numData)%y_coord)then
      temp = nint(numData_copy / numClusters_copy)
      do i = 1, numClusters
         clusters(i,1) = dataPoints(1 + (temp * (i - 1)))
      end do
      do j = 1, 50
         call pickCenters(dataPoints, clusters, numData, numClusters)
         do i = 1, numClusters
            call setTempList(clusters, tempList, numData, i)
            call findCentroid(tempList, numData)
            call resetCluster(clusters, tempList, numData, i)
         end do
         do i = 1, numData
         end do
      end do
   end if

! Write the clusters
   if(dataPoints(1)%x_coord==dataPoints(numData)%x_coord.and.dataPoints(1)%y_coord==dataPoints(numData)%y_coord)then
      write (*,*) "The data is all the same."
      open (2, file = clusteredFile)
      write (2, *) "Cluster 1:"
      do i = 1, numData
         if (abs(dataPoints(1)%x_coord) < 0.1 .or. abs(dataPoints(1)%x_coord) > 1e7) then
            write (2, '(A1,ES13.2,A1)', Advance = 'No') "(", dataPoints(i)%x_coord, ","
         else
            write (2, '(A1,F12.2,A1)', Advance = 'No') "(", dataPoints(i)%x_coord, ","
         end if
         if (abs(dataPoints(1)%y_coord) < 0.1 .or. abs(dataPoints(1)%y_coord) > 1e7) then
            write (2, '(ES13.2,A1)') dataPoints(i)%y_coord, ")"
         else
            write (2, '(F12.2,A1)') dataPoints(i)%y_coord, ")"
         end if
      end do
   else
      m = 0
      n = 0
      open (2, file = clusteredFile)
      do i = 1, numClusters
         if (m == n) then
            write (2, '(A8,I5,A1)') "Cluster", (i - n), ":"
         else
            m = m + 1
         end if
         call setTempList(clusters, tempList, numData, i)
         k = 0
         do j = 2, numData
            if (tempList(j)%x_coord /= 0 .and. tempList(j)%y_coord /= 0) then
               k = k + 1
            end if
         end do
         if (k > 0) then
            do j = 2, numData
               if (tempList(j)%x_coord /= 0 .and. tempList(j)%y_coord /= 0) then
                  if (abs(tempList(j)%x_coord) < 0.1 .or. abs(tempList(j)%x_coord) >= 1e7) then
                     write (2, '(A1,ES16.2,A1)', Advance = 'No') "(", tempList(j)%x_coord, ","
                  else
                     write (2, '(A1,F12.2,A1)', Advance = 'No') "(", tempList(j)%x_coord, ","
                  end if
                  if (abs(tempList(j)%y_coord) < 0.1 .or. abs(tempList(j)%y_coord) >= 1e7) then
                     write (2, '(ES16.2,A1)') tempList(j)%y_coord, ")"
                  else
                     write (2, '(F12.2,A1)') tempList(j)%y_coord, ")"
                  end if
               end if
            end do
         else
            n = n + 1
         end if
      end do
   end if
   write (*,*)

contains

! This subroutine sorts the data based on distance from a center point
subroutine sortPoints(dataPoints, center, numData)
   implicit none
   type(Point) :: dataPoints(:)
   type(Point), intent(in) :: center
   integer, intent(in) :: numData
   type(Point) :: temp
   integer i, j
   do i = 1, numData
      do j = i, numData
         if (dist(dataPoints(i), center) > dist(dataPoints(j), center)) then
            temp = dataPoints(i)
            dataPoints(i) = dataPoints(j)
            dataPoints(j) = temp
         end if
      end do
   end do
end subroutine sortPoints

! This subroutine determines which center each point is closest to
subroutine pickCenters(dataPoints, clusters, numData, numClusters)
   implicit none
   type(Point), intent(in) :: dataPoints(:)
   type(Point) :: clusters(:,:)
   integer, intent(in) :: numData
   integer, intent(in) :: numClusters
   real :: centerDist(numClusters)
   integer i, j, m, n
   do i = 1, numClusters
      do j = 2, numData
         clusters(i,j)%x_coord = 0
         clusters(i,j)%y_coord = 0
      end do
   end do
   do i = 1, numData
      centerDist(:) = 0
      do j = 1, numClusters
         centerDist(j) = dist(dataPoints(i), clusters(j,1))
      end do
      m = minloc(centerDist, dim=1)
      n = 1
      do while (clusters(m,n)%x_coord /= 0 .and. clusters(m,n)%y_coord /= 0)
         n = n + 1
         if (n > numData + 1) then
            exit
         end if
      end do
      if (n <= numData + 1) then
         clusters(m,n) = dataPoints(i)
      end if
   end do
end subroutine pickCenters

! This function finds the distance between two points
function dist(pointA, pointB)
   implicit none
   type(Point), intent(in) :: pointA, pointB
   real :: dist
   dist = ((pointB%x_coord - pointA%x_coord) ** 2)
   dist = dist + ((pointB%y_coord - pointA%y_coord) ** 2)
   dist = sqrt(dist)
end function dist

! This subroutine finds the centroid of a cluster
subroutine findCentroid(currCluster, numData)
   implicit none
   type(Point) :: currCluster(:)
   integer :: numData
   integer :: i, nonOrigin
   currCluster(1)%x_coord = 0
   currCluster(1)%y_coord = 0
   nonOrigin = 0
   do i = 1, numData
      currCluster(1)%x_coord = currCluster(1)%x_coord + currCluster(i+1)%x_coord
      currCluster(1)%y_coord = currCluster(1)%y_coord + currCluster(i+1)%y_coord
      if (currCluster(i+1)%x_coord /= 0 .and. currCluster(i+1)%y_coord /= 0) then
         nonOrigin = nonOrigin + 1
      end if
   end do
   currCluster(1)%x_coord = currCluster(1)%x_coord / nonOrigin
   currCluster(1)%y_coord = currCluster(1)%y_coord / nonOrigin
end subroutine findCentroid

! This subroutine sets tempList equal to the current cluster
subroutine setTempList(clusters, tempList, numData, clusterID)
   implicit none
   type(Point) :: clusters(:,:)
   type(Point) :: tempList(:)
   integer :: numData, clusterID
   integer :: i
   do i = 1, numData + 1
      tempList(i)%x_coord = 0
      tempList(i)%y_coord = 0
   end do
   do i = 1, numData + 1
      if (clusters(clusterID,i)%x_coord /= 0 .and. clusters(clusterID,i)%y_coord /= 0) then
         tempList(i) = clusters(clusterID,i)
      end if
   end do
end subroutine setTempList

! This subroutine sets the current cluster equal to tempList
subroutine resetCluster(clusters, tempList, numData, clusterID)
   implicit none
   type(Point) :: clusters(:,:)
   type(Point) :: tempList(:)
   integer :: numData, clusterID
   integer :: i
   do i = 2, numData + 1
      clusters(clusterID,i)%x_coord = 0
      clusters(clusterID,i)%y_coord = 0
   end do
   do i = 1, numData + 1
      if (tempList(i)%x_coord /= 0 .and. tempList(i)%y_coord /= 0) then
         clusters(clusterID,i) = tempList(i)
      end if
   end do
end subroutine resetCluster

end program k_means_clustering
