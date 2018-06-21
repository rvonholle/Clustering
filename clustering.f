      program clustering
      real, allocatable :: dataPoints(:)
      real :: temp
      character(len=100) :: dataFile, clusteredFile
      integer :: numGroups, numData, groupSize, i, j, k

c     This program reads a data file and clusters the data into a certain number of groups
      print *, "What file is the data being read from?"
      read *, dataFile
      print *, "How many groups will the data be clustered into?"
      read *, numGroups

c     Find number of data points in file
      numData = 0
      open (1, file = dataFile)
      do 
         read (1,*, end=10)
         numData = numData + 1
      end do
 10    close (1)
      allocate(dataPoints(numData))

c     Sort the data points
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

c     Cluster the data
      print *, "What file will the data be clustered into?"
      read *, clusteredFile
      groupSize = numData / numGroups
      i = 1
      j = 1
      k = 1
      open (2, file = clusteredFile)
      do i = 1, numGroups
         write (2,*) "Cluster", i, ":"
         do j = 1, groupSize
            if (dataPoints(k) > -0.1 .and. dataPoints(k) < 0.1) then
               write (2, '(ES15.7)') dataPoints(k)
            else
               write (2, '(F11.7)') dataPoints(k)
            end if
            k = k + 1
         end do
      end do
      if (k .lt. numData) then
         write (2,*) "Remainder cluster:"
         do i = k, numData
            write (2,*) dataPoints(i)
         end do
      end if
      close (2)
      end
