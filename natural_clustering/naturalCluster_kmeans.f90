program naturalCluster_kmeans
   implicit none
   real, dimension(:,:), allocatable :: clusters(:,:)
   real, allocatable :: dataPoints(:)
   real, allocatable :: tempList_1(:)
   real, allocatable :: tempList_2(:)
   real, allocatable :: compList(:)
   real :: numData_copy, numClusters_copy, temp, histRange, rangeLimit
   character(len=100) :: dataFile, clusteredFile
   character(len=1) :: chooseClusters
   integer :: numData, numClusters, numBins, i, j, k, m, n, fivesCounter

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
   if (chooseClusters == "y") then
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
   numData_copy = numData
   allocate(dataPoints(numData))
   open (1, file = dataFile)
   do i = 1, numData
      read (1,*) dataPoints(i)
   end do
   close (1)
   do i = 1, numData
      do j = i, numData
         if (dataPoints(i) > dataPoints(j)) then
            temp = dataPoints(i)
            dataPoints(i) = dataPoints(j)
            dataPoints(j) = temp
         end if
      end do
   end do

! Determine number of clusters and histogram range
   histRange = dataPoints(numData) - dataPoints(1)
   histRange = histRange / numBins
   if (chooseClusters == "n") then
      numClusters_copy = sqrt(numData_copy / 2)
      numClusters = nint(numClusters_copy)
   else
      do while (numClusters >= numData .and. chooseClusters /= 'n')
         do while (numClusters > numData)
            write (*,*) "There are more clusters than there is data."
            write (*,*) "Please choose a smaller number of clusters."
            read *, numClusters
         end do
         do while (numClusters == numData)
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
   allocate(tempList_1(numData+1))
   allocate(tempList_2(numData))
   allocate(compList(numData))

! Cluster the points
   if (dataPoints(1) /= dataPoints(numData)) then
      temp = (dataPoints(numData) - dataPoints(1)) / (numClusters + 1)
      do i = 1, numClusters
         clusters(i,1) = dataPoints(1) + (temp * i)
      end do
      do j = 1, 100
         do i = 1, numData
            tempList_2(i) = 0
         end do
         do i = 1, numClusters
            compList(i) = clusters(1,i+1)
         end do
         call pickCenters(dataPoints, clusters, numData, numClusters)
         do i = 1, numClusters
            call setTempList(clusters, tempList_1, numData, i)
            call findCentroid(tempList_1, numData)
            call resetCluster(clusters, tempList_1, numData, i)
         end do
         do i = 1, numData
            tempList_2(i) = clusters(1,i+1)
         end do
         if (all(compList == tempList_2)) then
            exit
         end if
      end do
   end if

! Write the clusters and the histogram
   if (dataPoints(1) == dataPoints(numData)) then
      write (*,*) "The data is all the same. No histogram will be written."
      open (2, file = clusteredFile)
      write (2, *) "Cluster 1:"
      do i = 1, numData
         if (abs(dataPoints(1)) < 0.1) then
            write (2, '(ES13.7)') dataPoints(i)
         else
            write (2, '(F12.7)') dataPoints(i)
         end if
      end do
   else
      write (*,*) "A histogram of the data will be displayed."
      write (*,*) "How many bins will be in the histogram?"
      read *, numBins
      write (*,*)
      write (*,*) "Key:"
      write (*,*) "# = 5 data points, | = gap between clusters"
      m = 0
      n = 0
      fivesCounter = 0
      rangeLimit = dataPoints(1) + histRange
      if (abs(dataPoints(1)) < 0.1 .or. abs(dataPoints(1)) >= 1e7) then
         write (*, '(A1, ES13.7, A1)', Advance = 'No') "[", dataPoints(1), ","
      else
         write (*, '(A1, F9.7, A1)', Advance = 'No') "[", dataPoints(1), ","
      end if
      if (abs(rangeLimit) < 0.1 .or. abs(rangeLimit) >= 1e7) then
         write (*, '(ES13.7, A2)', Advance = 'No') rangeLimit, "):"
      else
         write (*, '(F9.7, A2)', Advance = 'No') rangeLimit, "):"
      end if
      open (2, file = clusteredFile)
      do i = 1, numClusters
         if (m == n) then
            write (2, '(A8, I5, A1)') "Cluster", (i - n), ":"
            if (i > 1) then
               write (*, '(A1)', Advance = 'No') "|"
            end if
         else
            m = m + 1
         end if
         call setTempList(clusters, tempList_1, numData, i)
         k = 0
         do j = 2, numData
            if (tempList_1(j) /= 0) then
               k = k + 1
            end if
         end do
         if (k > 0) then
            do j = 2, numData
               if (tempList_1(j) /= 0) then
                  if (tempList_1(j) > rangeLimit) then
                     if (fivesCounter >= 1) then
                        write (*, '(A1)') "#"
                     else
                        write (*,*)
                     end if
                     if (abs(rangeLimit) < 0.1 .or. abs(rangeLimit) >= 1e7) then
                        write (*, '(A1, ES13.7, A1)', Advance = 'No') "[", rangeLimit, ","
                     else
                        write (*, '(A1, F9.7, A1)', Advance = 'No') "[", rangeLimit, ","
                     end if
                     rangeLimit = rangeLimit + histRange
                     if (abs(rangeLimit) < 0.1 .or. abs(rangeLimit) >= 1e7) then
                        write (*, '(ES13.7, A2)', Advance = 'No') rangeLimit, "):"
                     else
                        write (*, '(F9.7, A2)', Advance = 'No') rangeLimit, "):"
                     end if
                     fivesCounter = 0
                  end if
                  if (abs(tempList_1(j)) < 0.1 .or. abs(tempList_1(j)) >= 1e7) then
                     write (2, '(ES16.7)') tempList_1(j)
                  else
                     write (2, '(F12.7)') tempList_1(j)
                  end if
                  fivesCounter = fivesCounter + 1
                  if (fivesCounter == 5) then
                     write (*, '(A1)', Advance = 'No') "#"
                     fivesCounter = 0
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

! This subroutine will determine which center each point is closest to
subroutine pickCenters(dataPoints, clusters, numData, numClusters)
   implicit none
   real :: dataPoints(:)
   real :: clusters(:,:)
   integer :: i, j, m, n, numData, numClusters
   real, allocatable :: centerDist(:)
   allocate(centerDist(numClusters))
   do i = 1, numClusters
      do j = 2, numData
         clusters(i,j) = 0
      end do
   end do
   do i = 1, numData
      do j = 1, numClusters
         centerDist(j) = abs(dataPoints(i) - clusters(j,1))
      end do
      m = minloc(centerDist, dim=1)
      n = 0
      do while (clusters(m,n) /= 0)
         n = n + 1
      end do
      clusters(m,n) = dataPoints(i)
   end do
end subroutine pickCenters

! This subroutine finds the centroid of a cluster
subroutine findCentroid(currCluster, numData)
   implicit none
   real :: currCluster(:)
   integer :: numData
   integer :: i, nonZero
   currCluster(1) = 0
   nonZero = 0
   do i = 1, numData
      currCluster(1) = currCluster(1) + currCluster(i+1)
      if (currCluster(i+1) /= 0) then
         nonZero = nonZero + 1
      end if
   end do
   currCluster(1) = currCluster(1) / nonZero
end subroutine findCentroid

! This subroutine sets tempList equal to the current cluster
subroutine setTempList(clusters, tempList, numData, clusterID)
   implicit none
   real :: clusters(:,:)
   real :: tempList(:)
   integer :: numData, clusterID
   integer :: i
   do i = 1, numData + 1
      tempList(i) = 0
   end do
   do i = 1, numData + 1
      tempList(i) = clusters(clusterID, i)
   end do
end subroutine setTempList

! This subroutine sets the current cluster equal to templist
subroutine resetCluster(clusters, templist, numData, clusterID)
   implicit none
   real :: clusters(:,:)
   real :: tempList(:)
   integer :: numData, clusterID
   integer :: i
   do i = 2, numData + 1
      clusters(clusterID, i) = 0
   end do
   do i = 1, numData + 1
      clusters(clusterID, i) = tempList(i)
   end do
end subroutine resetCluster

end program naturalCluster_kmeans
