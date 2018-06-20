      program naturalCluster
      real, allocatable :: dataPoints(:)
      real, allocatable :: pointGap(:)
      real :: metric, temp
      character(len=100) :: dataFile, clusteredFile
      integer :: numData, i, j

c     This program reads a data file and clusters the data into natural clusters
      print *, "What file is the data being read from?"
      read *, dataFile

c     Read and sort the data points
      numData = 0
      open (1, file = dataFile)
      do
         read (1,*, end=10)
         numData = numData + 1
      end do
 10    close (1)
      allocate(dataPoints(numData))
      allocate(pointGap(numData - 1))
      open (1, file = dataFile)
      do i = 1, numData
         read (1,*) dataPoints(i)
      end do
      close (1)
      do i = 1, numData
         do j = i, numData
            if (dataPoints(i) .gt. dataPoints(j)) then
               temp = dataPoints(i)
               dataPoints(i) = dataPoints(j)
               dataPoints(j) = temp
            end if
         end do
      end do

c     Determine the metric
      metric = 0
      do i = 1, (numData - 1)
         pointGap(i) = dataPoints(i + 1) - dataPoints(i)
      end do
      do i = 1, (numData - 1)
         metric = metric + pointGap(i)
      end do
      metric = metric / (numData - 1)

c     Cluster the data
      print *, "What file will the data be clustered into?"
      read *, clusteredFile
      j = 1
      open (2, file = clusteredFile)
      write (2,*) "Cluster", j, ":"
      do i = 1, numData
         if (pointGap(i - 1) .gt. metric) then
            j = j + 1
            write (2,*) "Cluster", j, ":"
            write (2,*) dataPoints(i)
         else
            write (2,*) dataPoints(i)
         end if
      end do
      close (2)
      end
