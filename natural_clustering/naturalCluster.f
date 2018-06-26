      program naturalCluster
      real, allocatable :: dataPoints(:)
      real, allocatable :: pointGap(:)
      real :: metric, temp, histRange, rangeLimit
      character(len=100) :: dataFile, clusteredFile
      integer :: numData, numBins, i, j, k, fivesCounter

c     This program reads a data file and clusters the data into natural clusters
      print *, "What file is the data being read from?"
      read *, dataFile
      print *, "What file will the data be clustered into?"
      read *, clusteredFile
      print *, "A histogram of the data will be displayed."
      print *, "How many bins will be in the histogram?"
      read *, numBins
      print *,
      print *, "Key:"
      print *, "#=5 data points, |=gap between clusters"

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
            if (dataPoints(i) > dataPoints(j)) then
               temp = dataPoints(i)
               dataPoints(i) = dataPoints(j)
               dataPoints(j) = temp
            end if
         end do
      end do

c     Determine the metric and the histogram range
      metric = 0
      do i = 1, (numData - 1)
         pointGap(i) = dataPoints(i + 1) - dataPoints(i)
      end do
      do i = 1, (numData - 1)
         metric = metric + pointGap(i)
      end do
      metric = metric / (numData - 1)
      metric = metric * 2
      histRange = dataPoints(numData) - dataPoints(1)
      histRange = histRange / numBins

c     Cluster the data and create the histogram
      j = 1
      fivesCounter = 0
      rangeLimit = dataPoints(1) + histRange
      open (2, file = clusteredFile)
      write (2, '(A8, I5, A1)') "Cluster", j, ":"
      if (dataPoints(1) > -0.1 .and. dataPoints(1) < 0.1) then
         write (*, '(A1, ES13.7)', Advance = 'No') "[", dataPoints(1)
         if (rangeLimit > -0.1 .and. rangelimit < 0.1) then
            write(*,'(A1,ES13.7,A2)',Advance='No') ",", rangeLimit, "):"
         else
            write(*,'(A1,F9.7,A2)',Advance='No') ",", rangeLimit, "):"
         end if
      else
         write (*, '(A1, F12.7)', Advance = 'No') "[", dataPoints(1)
         if (rangeLimit > -0.1 .and. rangeLimit < 0.1) then
            write(*,'(A1,ES13.7,A2)',Advance='No') ",", rangeLimit, "):"
         else
            write(*,'(A1,F9.7,A2)',Advance='No') ",", rangeLimit, "):"
         end if
      end if
      do i = 1, numData
         if (dataPoints(i) > rangeLimit) then
            if (fivesCounter > 0) then
               write (*,'(A1)', Advance = 'No') "#"
            end if
            fivesCounter = 0
            write (*,*)
            if (rangeLimit + histRange < dataPoints(numData)) then
              if (rangeLimit > -0.1 .and. rangeLimit < 0.1) then
                write(*,'(A1,ES13.7,A1)',Advance='No')"[",rangeLimit,","
              else
                 write(*,'(A1,F9.7,A2)',Advance='No')"[",rangeLimit,","
              end if
              rangeLimit = rangeLimit + histRange
              if (rangeLimit > -0.1 .and. rangeLimit < 0.1) then
                 write(*,'(ES13.7,A2)',Advance='No') rangeLimit, "):"
              else 
                 write (*, '(F9.7, A2)', Advance='No') rangeLimit, "):"
              end if
            end if
         end if
         fivesCounter = fivesCounter + 1
         if (pointGap(i - 1) > metric) then
            j = j + 1
            write (2, '(A8, I5, A1)') "Cluster", j, ":"
            if (fivesCounter > 0) then
               write (*,'(A1)', Advance = 'No') "#"
            end if
            fivesCounter = 0
            write (*, '(A1)', Advance = 'No') "|"
            if (dataPoints(i) > -0.1 .and. dataPoints(i) < 0.1) then
               write (2, '(ES16.7)') dataPoints(i)
            else
               write (2, '(F12.7)') dataPoints(i)
            end if
            if (fivesCounter == 5) then
               write (*, '(A1)', Advance = 'No') "#"
               fivesCounter = 0
            end if
         else
            if (dataPoints(i) > -0.1 .and. dataPoints(i) < 0.1) then
               write (2, '(ES16.7)') dataPoints(i)
            else
               write (2, '(F12.7)') dataPoints(i)
            end if
            if (fivesCounter == 5) then
               write (*, '(A1)', Advance = 'No') "#"
               fivesCounter = 0
            end if
         end if
      end do
      close (2)
      write (*,*)
      end
