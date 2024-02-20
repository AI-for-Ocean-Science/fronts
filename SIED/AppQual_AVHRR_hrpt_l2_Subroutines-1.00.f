c***********************************************************************
c     Subroutines and functions in this package:
c     
c     function CommonAppQualSubsVersionNo()
c     The above function was replaced by code in the subroutine below.
c
c     subroutine AVHRR_hrpt_l2_GetDateTime(ImageName, Year, Month, 
c     1     YearT, MonthT, Day, Hour, Minute, Second)
c     
c     subroutine AVHRR_hrpt_l2p_GetGeoFilename(ImageName, YearC, MonthC, 
c     1     GeoNameFull)
c     
c     subroutine AVHRR_hrpt_l2_GetSSTOutData( FileNameIn, ConditionedSST, 
c     1     SecondsSince1970)
c     
c     subroutine AVHRR_hrpt_l2_CreateOutputSSTFile( FileNameOut, ncOutID, 
c     1     sstIDout, ClearIDout, FrontsIDout, DayOrNightIDout,
c     2     SecondsSinceID, Start, YearC, MonthC, DayC)
c     
c      subroutine AVHRR_hrpt_l2_GetGeoData( FileNameIn, LatIn, LonIn, 
c     1     SecondsSince1970, FileNameOut, NumClear, YearC, MonthC, Day 
c     2     Hour, Minute, Second, GeoNameFull)
c     
c***********************************************************************
c     Modifications 
c     
c     1.00
c     
c     8/27/13 - PCC - Created file from modified subroutines for AppQual
c     1/7/14 - PCC - Changed the name of this file to 
c       AppQual_AMSRE_v7_l3 to differentiate it from subroutines for 
c       other versions and levels.
c      Changed name: AMSRE_GetGeoFilename to AMSRE_v7_l2_GetGeoFilename
c      Removed CommonAppQualSubsVersionNo because this function has to
c       to be called by all versions of the subroutine package. Instead,
c       the variable CommonAppQualSubsVersionNo has been added to the
c       subroutine that reads in the date and time: xxx_GetDatTime and
c       is returned from that subroutine. The version number for this 
c       group of subroutines is set there REMEMBER TO CHANGE IT WHEN
c       UPDATING THE THIS PACKAGE.
c     12/16/14 - JPS - Modified from AMSRE_v7_l2 to process the L2P gac 
c       AVHRR data
c     
c***********************************************************************
      subroutine AVHRR_hrpt_l2_GetDateTime( ImageName, Year, Month, 
     1     YearT, MonthT, Day, Hour, Minute, Second)
c***********************************************************************
c     This subroutine gets the year, month, day, hour, minute and
c     seconds for the image specified by the image name. These times
c     are approximate and are used to determine if this image is in
c     in the time range specified. They should be considered accurate
c     the nearest hour. The precise time will be determined from 
c     metadata in the input file if available.
c     
      implicit none

c     Functions

      real secnds
      real*8 SecondsSince

c******Parameter statements

      include 'ParameterStatements.f'

c     General variables

c     Day - Day of the satellite pass.
      integer Day 

c     Hour - Hour of the satellite pass.
      integer Hour 

c     ImageName - the filename for the input image.
      character*319 ImageName
      
c     Minute - Minute of the satellite pass.
      integer Minute 

c     Month - This is the month in the loop over months.
      integer Month
c     MonthT - Month of the satellite pass.
      integer MonthT

c     Second - Second of the satellite pass.
      integer Second

c     Year - This is the year in the loop over years.
      integer Year
c     YearT - Year of the satellite pass.
      integer YearT

c---------------------------------start here --------------------------- 
      print *, 'ERROR THIS TEMPLATE SUBROUTINE HAS NOT BEEN EDITED', 
     1     ' YET GO THROUGH LINE BY LINE'
      return
      
      if(Debug .eq. 1) print *, 'GetTimes #100'

c       The next line sets the version of the subroutines in this file.
c       It needs to be incremented whenever this file is modified.

      CommonAppQualSubsVersionNo = '1.01'

c     Get year, month, day. 

      read( ImageName(1:4), '(i4)') YearT
      read( ImageName(5:6), '(i2)') MonthT
      read( ImageName(7:8), '(i2)') Day

      if(Debug .eq. 1) print *, 'GetTimes #110: YearT, ',
     1     'MonthT, Day: ', YearT, MonthT, Day

c     Make sure that the time in the satellite name corresponds to the
c     year and month being read.

      if((Year .ne. YearT) .or. (Month .ne. MonthT)) then
         write(UnitLog, *) 'GetTimes #120. Problem with times. ',
     1        'Year and month in filename: ', YearT, MonthT,
     2        ' does not match that in loops:', Year, Month
         print *, 'GetTimes #120. Problem with times. ',
     1        'Year and month in filename: ', YearT, MonthT,
     2        ' does not match that in loops:', Year, Month
         stop
      endif

c     Set the flag that says that the times written out are GMT.

      LSTorGMT = 1

c     Get seconds since the start of 1970. Will get a more accurate 
c     time later when it reads from the file, but have to decompress
c     and do all that, so skip for now.

      Hour = 0
      Minute = 0
      Second = 0

      if(Debug .eq. 1) print *, 'GetTimes #999'

      end subroutine AVHRR_hrpt_l2_GetDateTime

c***********************************************************************
      subroutine AVHRR_hrpt_l2_GetGeoFilename(ImageName, YearC, 
     1     MonthC, GeoNameFull)
c***********************************************************************
c     
      implicit none

c     Functions

      real secnds
      real*8 SecondsSince

c******Parameter statements

      include 'ParameterStatements.f'

c     General variables

c     ImageName - the filename for the input image.
      character*319 ImageName
      
c---------------------------------start here --------------------------- 
      print *, 'ERROR THIS TEMPLATE SUBROUTINE HAS NOT BEEN EDITED', 
     1     ' YET GO THROUGH LINE BY LINE'
      return
      
      if(Debug .eq. 1) print *, 'GetGeoFileName #100'

c     Build the output geolocation filename. GeoName is written out to
c     the SST file. GeoNameFull is used to create the geospatial file.

      GeoName = BlankFill(:len(GeoName))
      GeoName = 'GeoLoc/' // YearC // '/' // Monthc // '/' // 
     1     ImageName(1:Loc1-1) // '_geolocation.nc'

      if(Debug .eq. 1) print *, 'GetGeoFileName #110: ',
     1     'GeoName::', trim(GeoName), '::'
      
      GeoNameFull = BlankFill(:len(GeoNameFull))
      GeoNameFull = trim(BaseDir) // trim(GeoName)

      if(Debug .eq. 1) print *, 'GetGeoFileName #999'

      end subroutine AVHRR_hrpt_l2_GetGeoFilename

c***********************************************************************
      subroutine AVHRR_hrpt_l2_GetSSTOutData( FileNameIn, sstOut, 
     1     SecondsSince1970)
c***********************************************************************
c     
c     This subroutine will get the input SST data for AMSR-E. 
c     
c     Will also divide the input by 10 as agreed with Pierre for CMS
c     SSTs in centidegrees.
c     
c     For example, if the lowest acceptable value is -100 and the 
c     largest acceptable value is 1400, the Offset would be -100 and
c     the Range 1500.

c******Subroutine arugments
c     
c     sstIn - the input array.
c     sstOut - the output - conditioned - array.

c******General Variables
c     
c     ix, iy - do loop parameters
c     
c     NumberTooSmall - number of pixels that fall below Offset, 
c     . excluding FillValues.
c     NumberTooLarge - the number of pixels that are above Range plus
c     . Offset, excluding FillValues.
c     
c     MinSST - the min value in the input array excluding FillValues
c     MaxSST - the max value in the input array excluding FillValues
c     
      
      use netcdf

      implicit none

c     Function statements

      real*8 SecondsSince

c     Parameter statements

      include 'ParameterStatements.f'
      
c     Commons

c     RsstFillValue - missng value in input arrays.
      integer*2 RsstFillValue

c     General variables

      character*319 FileNameIn

c     ConditionedSST - SST after subtracting Offset and truncating.
      integer*2, allocatable :: ConditionedSST(:,:)

c     ncInID - netCDF identifier for the input file.
      integer ncInID
c     nxDimID - netCDF identifier for the x dimension variable
      integer nxDimID
c     nxN - length of the x dimension
      integer nxN
c     nyDimID - netCDF identifier for the y dimension variable
      integer nyDimID
c     nyN - length of the y dimension
      integer nyN

c     proxIDin - netCDF ID for ProximityConfidence
      integer proxIDin
c     ProximityConfidence - proximity flag in input: 
c     1=Bad, data rejected
c     2=Suspected Bad, data that has any confidence flags bit 0-5 thrown
c     3=Unprocessed proximity confidence flag, should be Good data
c     4=Good data
      byte, allocatable :: ProximityConfidence(:,:,:)

c     RejectionFlag
c     bit 0: near land or ice
c     bit 1: RFI
c     bit 2: sunglint
c     bit 3: rain detected
c     big 4: ic
c     bit 5: wind > 20 m/s or SST out of bounds
c     bit 6: land
c     bit 7: edge of swath " ;
      byte, allocatable :: RejectionFlag(:,:,:)
c     rejIDin - netCDF ID for RejectionFlag
      integer rejIDin
      
c     SecondsIn - reference time, seconds since 1981/1/1 00:00 0.0
      real*4 SecondsIn
c     SecondsSince1970 - seconds since 1970 passed in, which is a rough
c     number. It is updated in this routine and the correct value is
c     returned.
      real*8 SecondsSince1970
c     SecondsSinceTemp - the temporary value of seconds since 1970 used
c     to store the good value for comparison with the values passed in.
      real*8 SecondsSinceTemp
c     sstIn - SST read in.
      real*4, allocatable :: sstIn(:,:)
c     sstIDin - netCDF identifier for the input SST field
      integer sstIDin
c     sstOut - SST written out.
      integer*2 sstOut(LenX, LenY) 
c     sstQualIn - quality for SST on input
      byte, allocatable :: sstQualIn(:,:)
c     sstQualIDin - netCDF identifier for the input SST qualifier field
      integer sstQualIDin

c     TimeIDin - netCDF variable ID for reference time.
      integer TimeIDin
c     TempSST - variable used to read SST, allocated for each file.
      integer*2, allocatable :: TempSST(:,:,:)

      integer*2 MinSST, MaxSST, NumberTooSmall, NumberTooLarge

      integer SecondsSinceID

      integer*4 ix, jy

c----------------------------------start here --------------------------
c     

      print *, 'ERROR THIS TEMPLATE SUBROUTINE HAS NOT BEEN EDITED', 
     1     ' YET GO THROUGH LINE BY LINE'
      return
      
      if(Debug .eq. 1) print *, 'AVHRR_GetSSTData #000: ',
     1     'trim(FileNameIn)::', trim(FileNameIn), '::'

c     Allocate input and output variables.

      allocate( sstIn(LenX,LenY),
     1     sstQualIn(LenX,LenY),
     2     ConditionedSST(LenX,LenY),
     *     stat=ierr)
      
      if (ierr .ne. 0) then
         write(UnitLog,*) 'AVHRR_GetSSTData #100: Allocation error, ',
     1        'exiting.'
         stop 'AVHRR_GetSSTData #100: Allocation error, exiting.'
      endif

c     Open the input file.

      status = nf90_open( FileNameIn, nf90_nowrite, ncInID)
      if(status .ne. nf90_noerr) call handle_err(status)

c     Get the time from inside the file. Will use the reference time, 
c     not, the time at each pixel location. 

      status = nf90_inq_varid( ncInID, 'time', TimeIDin)
      if(status .ne. nf90_noerr) call handle_err(status)

      status = nf90_get_var( ncInID, TimeIDin, SecondsIn)
      if(status .ne. nf90_noerr) call handle_err(status)

      if(Debug .eq. 1) print *, 'AVHRR_GetSSTData #001: ',
     1     'Got Variable: Time'

c     Is the time read within 24 hours of the time calculated earier?

      SecondsSinceTemp = SecondsSince( 1981, 1, 1, 0, 0, 0,
     1     1970, 1, 1, 0, 0, 0) + SecondsIn

      if(abs(SecondsSinceTemp-SecondsSince1970) .gt. 86400*3/2) then
         write(UnitLog, *) 'AVHRR_GetSSTData #110: There is a problem',
     1        'with the time read in: ', SecondsSinceTemp, 
     2        'It is not close to the value determined from the file ',
     3        'name: ', SecondsSince1970
         print *, 'AVHRR_GetSSTData #110: There is a problem',
     1        'with the time read in: ', SecondsSinceTemp, 
     2        'It is not close to the value determined from the file ',
     3        'name: ', SecondsSince1970
         stop
      endif

      SecondsSince1970 = SecondsSinceTemp

c     Get the input variable IDs

      status = nf90_inq_varid( ncInID, 'sea_surface_temperature', 
     1     sstIDin)
      if(status .ne. nf90_noerr) call handle_err(status)
      
      if(Debug .eq. 1) print *, 'AVHRR_GetSSTData #002: ',
     1     'Got Variable: SST'

      status = nf90_inq_varid( ncInID, 'quality_level', proxIDin)
      if(status .ne. nf90_noerr) call handle_err(status)
      
      status = nf90_inq_varid( ncInID, 'quality_level', rejIDin)
      if(status .ne. nf90_noerr) call handle_err(status)

c$$$  status = nf90_inq_varid( ncInID, 'diurnal_amplitude', dwIDin)
c$$$  if(status .ne. nf90_noerr) call handle_err(status)

c     Get the length of the input dimensions to make sure that they
c     are consistent with xLenIn and yLenIn

      status = nf90_inq_dimid( ncInID, 'ni', nxDimID)
      if(status .ne. nf90_noerr) call handle_err(status)

      status = nf90_inq_dimid( ncInID, 'nj', nyDimID)
      if(status .ne. nf90_noerr) call handle_err(status)

      status = nf90_inquire_dimension( ncInID, nxDimID, len=nxN)
      if(status .ne. nf90_noerr) call handle_err(status)

      status = nf90_inquire_dimension( ncInID, nyDimID, len=nyN)
      if(status .ne. nf90_noerr) call handle_err(status)

      if(debug .eq. 1) then
         write(UnitLog,*) 'AVHRR_GetSSTData #120: nxN, nyN: ', nxN, nyN
         print *, 'AVHRR_GetSSTData #120: nxN, nyN: ', nxN, nyN
      endif

      if( nxN .gt. LenXIn) then
         write(UnitLog,*) 'WriteSSTData #130 ERROR ** Long orbit: ', nxN
         print *, 'AVHRR_GetSSTData #130: ERROR ** Long orbit: ', nxN
         stop '******* Dimension error *******'
      endif

      if( nyN .gt. LenYIn) then
         write(UnitLog,*) 'AVHRR_GetSSTData #140: ERROR *** Bad # of ',
     1        'elements: ', nyN, ' should be less than', LenYIn
         print *, 'AVHRR_GetSSTData #140: Error *** Bad # of ',
     1        'elements: ', nyN, ' should be less than', LenYIn
         stop '******* Dimension error *******'
      endif

c     Get the missing values for the input SST, Lat and Lon

      status = nf90_get_att( ncInID, sstIDin, '_FillValue', 
     1     RsstFillValue)
      if(status .ne. nf90_noerr) call handle_err(status)
      
      if(debug .eq. 1) print *, 'AVHRR_GetSSTData #150: SST input  ',
     1     'fill value: ', RsstFillValue

c     Check to make sure that the offset and slope read in by ReadInputs
c      correspond to that in the input file.

      status = nf90_get_att( ncInID, sstIDin, 'add_offset', 
     1     sstOffsetIn)
      if(status .ne. nf90_noerr) call handle_err(status)
      
      status = nf90_get_att( ncInID, sstIDin, 'scale_factor', 
     1     SSTScaleFactorIn)
      if(status .ne. nf90_noerr) call handle_err(status)
      
      if((InputOffset .ne. sstOffsetIn) .or. 
     1     (InputScaleFactor .ne.SSTScaleFactorIn)) then
         print *, 'AVHRR_GetSSTData #155: SST offset from ',
     1        'run_control_file:  ', InputOffset, ' or scale factor ',
     2        InputScaleFactor, ' does not match that in the input ',
     3        'file: ', sstOffsetIn, ' or ', SSTScaleFactorIn
         print *, 'AVHRR_GetSSTData #156: STOPPING.'
      endif

c     Initialize input values to the missing value for that variable 
c     type. Data will be read into a temporary array if if the data
c     read in passes the quality tests it the values in sstIN, LatIn
c     and LonIn will be assigned the value read in.

      do 1033 jy=1,LenY
         do 1034 ix=1,LenX
c$$$  sstIn(ix,jy) = RsstFillValue
            sstIn(ix,jy) = FillValueReal
 1034    continue
 1033 continue

c     Now allocate variable for input SST.

      allocate(TempSST(nxN,nyN,1), 
     7     ProximityConfidence(nxN,nyN,1),
     8     RejectionFlag(nxN,nyN,1),
     9     stat=ierr)
      if (ierr .ne. 0) then
         write(UnitLog,*) 'AVHRR_GetSSTData #160: Allocation error.'
         stop 'AVHRR_GetSSTData #160: Allocation error. STOPPING.'
      endif

c     Now read the data

      status = nf90_get_var( ncInID, sstIDin, TempSST)
      If(status .ne. nf90_noerr) call handle_err(status)
      
      status = nf90_get_var( ncInID, rejIDin, RejectionFlag)
      If(status .ne. nf90_noerr) call handle_err(status)

c      status = nf90_get_var( ncInID, proxIDin, ProximityConfidence)
c      If(status .ne. nf90_noerr) call handle_err(status)

c     Load TempSST into sstIn and apply Proximity Confidence and 
c     Rejection Flag. At this point any RejectionFlag bit that is
c     set will result in the SST value being set to the fill value
c     and if the RejectionFlag is either 1 or 2 (bad or probably
c     bad), SST will also be set to the fill value

      do 1043 jy=1,nyN
         do 1044 ix=1,nxN
            if(RejectionFlag(ix,jy,1) .gt. 3) then 
               sstIn(ix,jy) = TempSST(ix,jy,1)
            endif
 1044    continue
 1043 continue

c     Now condition SST data. ******************************************

      call ConditionInput( sstIn, ConditionedSST)
      
c     Next Median filter the conditioned data.

      call median( ConditionedSST, SSTOut, 3)

c     Now clean up.

      status = nf90_close(ncInID)
      if(status .ne. nf90_noerr) call handle_err(status)

      end subroutine AVHRR_hrpt_l2_GetSSTOutData

c**********************************************************************
      subroutine AVHRR_hrpt_l2_CreateOutputSSTFile( FileNameOut, 
     1     ncOutID,ProgName, sstIDout, ClearIDout, FrontsIDout, 
     2     DayOrNightIDout, SecondsSinceID, Start, YearC, MonthC, DayC)
c**********************************************************************
c     
c     This subroutine will generate the output file.
c     
      use netcdf

      implicit none

c******Parameter statements

      include 'ParameterStatements.f'

c     Functions

c     AttDescription - attribute description when needed.
      character*1000 AttDescription

c     SatName_uc - the upper case version of the satellite name.
      character*20 SatName_uc

      character*200 Title
      character*67 Source
      character*218 Comments

c     Filenames

      character*319 FileNameOut

c     netCDF variables

      integer Dummy

      integer ncOutID, nxDimIDin, nyDimIDin, nxDimIDout, nyDimIDout
      integer sstIDout, SecondsSinceID
      integer ClearIDout, FrontsIDout, DayOrNightIDout

      integer sstDimsOut(2) 

      integer Values(8)
      character*5 Zone
      character*10 Time
      character*8 Date
      character*19 ProcessingDateTime
c      character*3000 SummaryPart1
      
      character*20 OffsetChar, ScaleInputChar, RangeChar

      character*16 Conventions
      character*64 Meta_Conventions
      character*64 standard_name_vocab
      character*64 creator_name
      character*64 creator_email
      character*256 creator_url
      character*1024 license
      character*512 contributor_name
      character*512 contributor_role
      character*256 publisher_name
      character*256 publisher_url
      character*256 publisher_institution
      character*32 cdm_data_type

      character*32 time_coverage_start
      real*4 time_coverage_duration

      integer i

      logical exis

c---------------------------------start here --------------------------


      print *, 'ERROR THIS TEMPLATE SUBROUTINE HAS NOT BEEN EDITED', 
     1     ' YET GO THROUGH LINE BY LINE'
      return
c***********************************************************************
c**** Define lobal attributes particular to this archive, but don't 
c     write them out yet. Will do this in the section below that does 
c     not change from archive-to-archive.
c     
c**** This section will change from archive-to-archive.

      if(Debug .eq. 1) print *, 'CreateMedian: #000. FileNameOut::', 
     1     trim(FileNameOut)

c     AMSR-E or TMI

      SatName_uc = trim(SatName)
      call upper_case(SatName_uc)

      Loc1 = index( SatName, 'AVHRR')
      
      Title = BlankFill(:len(Title))
      Title  = 'Conditioned ' // trim(SatName_uc) // ' '
     1     // trim(ParameterName) // ' fields from the Centre for '
     2     // 'Environmental Data Archival.'

      if(Debug .eq. 1) print *, 'CreateMedian: #100: SatName::', 
     1     trim(SatName), '::, SatName_uc::', trim(SatName_uc),
     2     '::, Title::', trim(Title)

c     Build the summary for this data set in 3 parts. Three parts are
c     needed because of formatting required for the variables in the 
c     middle part.

      SummaryPart1 = BlankFill(1:len(SummaryPart1))
      SummaryPart1 = 'The field in this file was generated from the '
     1     // trim(SatName_uc) // ' ' // trim(ParameterName) 
     2     // ' field obtained from the Centre for Environmental '
     3     // 'Data Archival. '
     4     // trim(ParameterName) // ' values for which the Rejection '
     5     // 'Flag was equal to 0 and the Proximity Confidence was '
     6     // 'greater than 2 (the Rejection Flag and Proximity '
     7     // 'Confidence are defined in the original ' 
     8     // trim(ParameterName) // ' field obtained from RSS) were '
     9     // 'scaled as follows: Counts_out = nint( Counts_in * '

      Source = BlankFill(:len(Source))
      Source = 'Satellite observation' 

      creator_name = BlankFill(:len(creator_name))
      creator_name = 'John Salter'

      creator_url = BlankFill(:len(creator_url))
      creator_url = 'http://www.sstfronts.org'

      creator_email = BlankFill(:len(creator_email))
      creator_email = 'john_salter@my.uri.edu'

      institution = BlankFill(:len(institution))
      institution = 'University of Rhode Island - ' //
     1     'Graduate School of Oceanography'

      project = BlankFill(:len(project))
      project = 'SST Fronts'

      acknowledgement = BlankFill(:len(acknowledgement))
      acknowledgement = 'The Graduate School of Oceanography ' //
     1     'generated this data file with support from X. ' //
     2     'Please acknowledge the use of these data with ' //
     3     'the following: {These data, were ' //
     4     'produced by the University of Rhode Island with ' //
     5     'sponsorship from X from SST fields ' //
     6     'developed by and obtained from the Centre for ' //
     7     ' Environmental Data Archival}  See the license ' //
     8     'attribute for information regarding use and ' //
     9     're-distribution of these data.'

      license = BlankFill(:len(license))
      license = 'The fields associated with the AVHRR data are '
     1     // 'publicly available from the Centre for Environmental '
     2     // 'Data Archival at: '
     3     // 'ftp.ceda.ac.uk/neodc/esacci_sst/data/lt/AVHRR/'
     4     // '. Please acknowledge use of these data with the '
     5     // 'text given in the acknowledgement attribute.'

      contributor_name = BlankFill(:len(contributor_name))
      contributor_name = 'Centre for Environmental Data Archival:' //
     1     ' http://www.ceda.ac.uk/'

      contributor_role = BlankFill(:len(contributor_role))
      if(Loc1 .eq. 1) then
         contributor_role = 'Generated SST fields from the AVHRR ' //
     1        'data stream.'
      else
         contributor_role = 'Generated SST fields from the TMI ' //
     1        'data stream.'
      endif

      publisher_name = BlankFill(:len(publisher_name))
      publisher_name = 'SST Fronts'

      publisher_url = BlankFill(:len(publisher_url))
      publisher_url = 'http://www.sstfronts.org'

      publisher_institution = BlankFill(:len(publisher_institution))
      publisher_institution = 'University of Rhode Island - ' //
     1     'Graduate School of Oceanography'

      cdm_data_type = BlankFill(:len(cdm_data_type))
      cdm_data_type = 'Grid'

      ProgVersionNo = AppQual_l2p_AVHRR_gac_v1_VersionNo

      if(Debug .eq. 1) print *, 'CreateMedian: #101: SatName::',
     1     'Ready to Call CreateOutputSSTFile'
      
      call CreateOutputSSTFile( FileNameOut, ncOutID, ProgName, 
     1     nxDimIDout, nyDimIDout, sstDimsOut, sstIDout, ClearIDout, 
     2     FrontsIDout, DayOrNightIDout, SecondsSinceID, Title, 
     3     SummaryPart1, creator_name, creator_url, creator_email,
     4     license, contributor_name, contributor_role, publisher_name, 
     5     publisher_url, publisher_institution, cdm_data_type)

      return

      end subroutine AVHRR_hrpt_l2_CreateOutputSSTFile

c***********************************************************************
      subroutine AVHRR_hrpt_l2_GetGeoData( FileNameIn, LatIn, LonIn, 
     1     SecondsSince1970, FileNameOut, NumClear, YearC, MonthC, Day, 
     2     Hour, Minute, Second, GeoNameFull)
c***********************************************************************
c     
c     This subroutine will get the input SST data for AVHRR. 
c     
c     Will also divide the input by 10 as agreed with Pierre for CMS
c     SSTs in centidegrees.
c     
c     For example, if the lowest acceptable value is -100 and the 
c     largest acceptable value is 1400, the Offset would be -100 and
c     the Range 1500.

c******Subroutine arugments
c     
c     sstIn - the input array.
c     sstOut - the output - conditioned - array.

c******General Variables
c     
c     ix, iy - do loop parameters
c     
c     NumberTooSmall - number of pixels that fall below Offset, 
c     . excluding FillValues.
c     NumberTooLarge - the number of pixels that are above Range plus
c     . Offset, excluding FillValues.
c     
c     MinSST - the min value in the input array excluding FillValues
c     MaxSST - the max value in the input array excluding FillValues
c     
      
      use netcdf

      implicit none

c     Function statements

      real*8 SecondsSince

c     Parameter statements

      include 'ParameterStatements.f'
      
c     Commons

c     RsstFillValue - missng value in input arrays.
      real*4 RlatFillValue
      real*4 RlonFillValue
      common /MissingValuesLatLon/ RlatFillValue, RlonFillValue

c     General variables

      character*319 FileNameIn, FileNameOut
      character*319 TempFile

c     CountAVHRR - the number of input elements to read.
      integer CountAVHRR(3)

c     geo_exis - .true. if geo file exists, .false. otherwise.
      logical geo_exis

c     i - do loop parameter
      integer i
c     iSave - length of the data read.
      integer iSave
c     iUpper - the upper limit of the 1st dimension do loop for sstIn.
      integer iUpper

c     j - do loop parameter
      integer j

c     LatIDin - netCDF identifier for the latitude variable.
      integer LatIDin
c     LatIn - latitude at each pixel locaion, read in.
      real*4 LatIn(LenX,LenY)
c     LatLonFileIn - input file if lat,lon are separate from the SST 
c      file.
      character*319 LatLonFileIn
c     LonIDin - netCDF identifier for the longitude variable.
      integer LonIDin
c     LonIn - longitude at each pixel locaion, read in.
      real*4 LonIn(LenX,LenY)

c     MinLat - the minumum Lat value read from the ascii file.
      real*4 MinLat
c     Minute - minute of reference time.
      integer Minute

c     ncInID - netCDF identifier for the input file.
      integer ncInID
c     ncOutID - netCDF identifier for the output file.
      integer ncOutID
c     NumClear - number of good pixels in this image.
      integer NumClear
c     nxDimID - netCDF identifier for the x dimension variable
      integer nxDimID
c     nxN - length of the x dimension
      integer nxN
c     nyDimID - netCDF identifier for the y dimension variable
      integer nyDimID
c     nyN - length of the y dimension
      integer nyN
      
c     Second - seconds of reference time.
      integer Second
c     SecondsIn - reference time, seconds since 1981/1/1 00:00 0.0
      real*4 SecondsIn
c     SecondsSince1970 - seconds since 1970 passed in, which is a rough
c     number. It is updated in this routine and the correct value is
c     returned.
      real*8 SecondsSince1970
c     SecondsSinceTemp - the temporary value of seconds since 1970 used
c     to store the good value for comparison with the values passed in.
      real*8 SecondsSinceTemp
c     SeparateLatLon - 1 if separate lat,lon file, 0 if lat, lon in SST 
c     file
      integer SeparateLatLon
c     StartAMSR - the staring row, line and 3rd dimension.
      integer StartAMSR(3)
c     StrideAMSR - the stride to use when reading in data.
      integer StrideAMSR(3)

c     TimeIDin - netCDF variable ID for reference time.
      integer TimeIDin
c     TempFile2 - temporary file for lat,lon in decompress.
      character*319 TempFile2
c     TempLatLon - variable used to read latitude and longitude, 
c     allocated for each file if read from txt file that Chelle
c     sent us.
      real*4, allocatable :: TempLatLon(:,:)

      integer Day, Hour

      integer*2 MinSST, MaxSST, NumberTooSmall, NumberTooLarge

      integer ClearIDout, FrontsIDout, DayOrNightIDout

      integer SecondsSinceID

      integer*4 ix, jy

      logical Compressed

      integer Storage

c----------------------------------start here --------------------------
c     

      print *, 'ERROR THIS TEMPLATE SUBROUTINE HAS NOT BEEN EDITED', 
     1     ' YET GO THROUGH LINE BY LINE'
      return
      
      if(Debug .eq. 1) print *, 'AVHRR_GetGeoData #000: ',
     1     'trim(FileNameIn)::', trim(FileNameIn), '::'

c     Switch here to read the lat,lon data from the input SST file or 
c     from a separate file.

      SeparateLatLon = 0

      if(SeparateLatLon .eq. 1) then

c     Allocate temp variable for latitude and longitude, 2d not 3d.

         allocate(TempLatLon(243,8900), stat=ierr)

         if (ierr .ne. 0) then
            write(UnitLog,*) 'AVHRR_GetGeoData #110: Allocation error.'
            print *,'AVHRR_GetGeoData #110: Allocation error. STOPPING.'
            stop
         endif

c     Need to read a separate lat,lon file. Generate the filename.

         if(Debug .eq. 1) print *, 'AVHRR_GetGeoData #120: ', 
     1        'FileNameIn::', trim(FileNameIn), '::'

         Loc1 = index( FileNameIn, '.dat')
         LatLonFileIn = blankfill(:len(LatLonFileIn))
         LatLonFileIn = '/Volumes/Original/archive/amsre-archive/' //
     1        'Navigation/0' // FileNameIn(Loc1-5:Loc1-1) // '.txt.gz'

         if(Debug .eq. 1) print *, 'AMSRE_GetGeoData #130, Loc1: ',
     2        Loc1, ' LatLonFileIn::', trim(LatLonFileIn), '::'

c     And the temporary filename to use for the decompressed lat,lon 
c     file.

         TempFile2 = blankfill(:len(TempFile2))
         TempFile2 = trim(BaseDir) // 'TmpDir/' // 
     1        FileNameIn(Loc1-5:Loc1-1) // '.txt'

         if(Debug .eq. 1) print *, 'AMSRE_GetGeoData #140, TempFile2::', 
     1        trim(TempFile2), '::'

         Command = BlankFill(:len(Command))
         Command = 'gzip -dc ' // trim(LatLonFileIn) // 
     1        ' > ' // trim(TempFile2)

         if(Debug .eq. 1) print *, 'AMSRE_GetGeoData #150::', 
     1        trim(Command), '::'

         call system(Command)

c     Now open and read the file.

         open(Unit=101, file=TempFile2, form='formatted', 
     1        status='unknown')

c     Get latitude and longitude. First zero out column 1 of the
c     temporary array. Will use this to check on the size of the
c     arrays read in.

         if(Debug .eq. 1) print *, 'AMSRE_GetGeoData #160. Read ',
     1        'Lat/Lon.'

         do 4000 i=1,8900
            TempLatLon(1,i) = -99999.0
 4000    continue

c     Now read the data

         read(101,4100,end=4101) TempLatLon
 4100    format(2138400F12.4)

 4101    continue

c     Get the length of the data read.

         iSave = 0
         do 4300 i=1,8900
            if(TempLatLon(1,i) .eq. -99999.0) go to 4301
            iSave = i
 4300    continue
 4301    continue

         if(( iSave .eq. 0) .or. ( iSave .ne. 8800) .or. 
     1        (Debug .eq. 1) ) print *, 'AMSRE_GetGeoData #180: ',
     2        'iSave=', iSave

         iSave = iSave / 2

c     Next separate the lat and lon fields.

         MinLat = 90.0
         iUpper = min( LenX, iSave)
         do 4400 i=1,iUpper
            do 4410 j=1,243
               Latin(i,j) = TempLatLon(j,i)
               Lonin(i,j) = TempLatLon(j,i+iSave)
               if(Latin(i,j) .lt. MinLat) MinLat = Latin(i,j)
 4410       continue
 4400    continue

         if(Debug .eq. 1) print *, 'AMSRE_GetGeoData #185: MinLat: ', 
     1        MinLat

c     Assuming that the minimum lat value is the fill value set
c      to lat and lon values equal to the minimum value to 
c      FillValueReal. This is because, for some reason, using the 
c      minimum value to flag missing values doesn't work. 

         do 4420 i=1,iUpper
            do 4430 j=1,243
               if(Latin(i,j) .eq. MinLat) then
                  Latin(i,j) = FillValueReal
                  Lonin(i,j) = FillValueReal
               endif
 4430       continue
 4420    continue

c     Close the input file and delete the temporary file.

         close(101)

         Command = BlankFill(:len(Command))
         Command = 'rm ' // trim(TempFile2)

         if(Debug .eq. 1) print *, 'AMSRE_GetGeoData #190::', 
     1        trim(Command), '::'

         call system(Command)

      else

c*************** Here if lat,lon in SST file ***************************

         status = nf90_open( FileNameIn, nf90_nowrite, ncInID)
         if(status .ne. nf90_noerr) call handle_err(status)

c     Get the input variable IDs

         status = nf90_inq_varid( ncInID, 'lat', latIDin)
         if(status .ne. nf90_noerr) call handle_err(status)
         
         status = nf90_inq_varid( ncInID, 'lon', lonIDin)
         if(status .ne. nf90_noerr) call handle_err(status)

c     Get the length of the input dimensions to make sure that they
c     are consistent with xLenIn and yLenIn

         status = nf90_inq_dimid( ncInID, 'ni', nxDimID)
         if(status .ne. nf90_noerr) call handle_err(status)

         status = nf90_inq_dimid( ncInID, 'nj', nyDimID)
         if(status .ne. nf90_noerr) call handle_err(status)

         status = nf90_inquire_dimension( ncInID, nxDimID, len=nxN)
         if(status .ne. nf90_noerr) call handle_err(status)

         status = nf90_inquire_dimension( ncInID, nyDimID, len=nyN)
         if(status .ne. nf90_noerr) call handle_err(status)

         if(debug .eq. 1) then
            write(UnitLog,*) 'AVHRR_GetGeoData #200, nxN, nyN: ', 
     1           nxN, nyN
            print *, 'AVHRR_GetGeoData #200, nxN, nyN: ', nxN, nyN
         endif

         if( nxN .gt. LenXIn) then
            write(UnitLog,*) 'AVHRR_GetSSTData #210: ERROR ** Long ',
     1           'orbit: ', nxN
            print *, 'AVHRR_GetSSTData #220:: ERROR ** Long orbit: ', 
     1           nxN
            stop '******* Dimension error *******'
         endif

         if( nyN .gt. LenYIn) then
            write(UnitLog,*) 'AVHRR_GetSSTData #220: ERROR *** Bad # ',
     1           'of elements: ', nyN
            print *, 'AVHRR_GetSSTData #220: Error *** Bad # of ',
     1           'elements: ', nyN, '. STOPPING'
            stop '******* Dimension error *******'
            return
         endif

c     Get the missing values for the input SST, Lat and Lon

c         status = nf90_get_att( ncInID, LatIDin, '_FillValue', 
c     1        RlatFillValue)
c         if(status .ne. nf90_noerr) call handle_err(status)
c         
c         status = nf90_get_att( ncInID, LonIDin, '_FillValue', 
c     1        RlonFillValue)
c         if(status .ne. nf90_noerr) call handle_err(status)

         RlatFillValue = FillValueReal
         RlonFillValue = FillValueReal
         
         if(debug .eq. 1) print *, 'AVHRR_GetGeoData #230, Lat, Lon ',
     1        'fill values: ', RlatFillValue, RlonFillValue

c     Initialize input values to the missing value for that variable 
c     type. Data will be read into a temporary array if if the data
c     read in passes the quality tests it the values in sstIN, LatIn
c     and LonIn will be assigned the value read in.

         do 1033 jy=1,LenY
            do 1034 ix=1,LenX
               LatIn(ix,jy) = FillValueReal
               LonIn(ix,jy) = FillValueReal
 1034       continue
 1033    continue

c     Allocate temp variable for latitude and longitude, 2d not 3d.

         allocate(TempLatLon(nxN,nyN), stat=ierr)

         if (ierr .ne. 0) then
            write(UnitLog,*) 'AppQual #100B: Allocation error, exiting.'
            stop
         endif

c     And read in latitude and longitude

         if(Debug .eq. 1) print *, 'AVHRR_GetGeoData #240, Got past ', 
     1        'read for SSTin.'

         status = nf90_get_var( ncInID, latIDin, TempLatLon)
         If(status .ne. nf90_noerr) call handle_err(status)

c     Load TempLatLon into LatIn

         do 1053 jy=1,nyN
            do 1054 ix=1,nxN
               LatIn(ix,jy) = TempLatLon(ix,jy)
 1054       continue
 1053    continue

c     Same drill for longitude

         status = nf90_get_var( ncInID, lonIDin, TempLatLon)
         If(status .ne. nf90_noerr) call handle_err(status)
         
c     Load TempLatLon into LonIn

         do 1063 jy=1,nyN
            do 1064 ix=1,nxN
               LonIn(ix,jy) = TempLatLon(ix,jy)
 1064       continue
 1063    continue

c     close the input file.

         status = nf90_close(ncInID)
         if(status .ne. nf90_noerr) call handle_err(status)

      endif

      if(Debug .eq. 1) print *, 'AVHRR_GetGeoData #999'

      end subroutine AVHRR_hrpt_l2_GetGeoData


