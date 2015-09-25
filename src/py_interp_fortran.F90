MODULE routines
IMPLICIT NONE
CONTAINS
!  Subroutines extracted from p_interp.F90 Authors:
!
!  November 2007 - Cindy Bruyere
!  December 2009 - Lluis Fita (Santander Meteorolgy Group, Univ. Cantabria)
!  January 2012 - J. Fernandez (Santander Meteorolgy Group, Univ. Cantabria) -> jf
!
SUBROUTINE interp (data_out, data_in, pres_field, interp_levels, psfc, ter, tk, qv, ix, iy, iz, it, &
                     num_metgrid_levels, LINLOG, extrapolate, GEOPT, MISSING)
    IMPLICIT NONE
    ! Input fields:
    !
    ! ix, iy, iz, it     NetCDF dimensions
    ! num_metgrid_levels Number pf pressure levels.
    ! LINLOG             Interpolation option (linear or log)
    ! data_in            Data of the input variable
    ! pres_field, tk, qv Pressure, temperature and specific humidity fields.
    ! psfc               Surface pressure
    ! ter                Terrain height
    ! interp_levels      Pressure levels to interpolate
    ! GEOPT              Boolean: True if we are interpolating the geopotential
    ! extrapolate        Boolean: Extrapolate below ground
    !
    INTEGER, INTENT(IN)                              :: ix, iy, iz, it    
    INTEGER, INTENT(IN)                              :: num_metgrid_levels, LINLOG
    REAL(8), DIMENSION(ix, iy, iz, it), INTENT(IN)      :: data_in, pres_field, tk, qv 

    REAL(8), DIMENSION(ix, iy, it), INTENT(IN)          :: psfc
    REAL(8), DIMENSION(ix, iy), INTENT(IN)              :: ter
    REAL(8), DIMENSION(num_metgrid_levels), INTENT(IN)  :: interp_levels
    REAL(8)                                             :: ptarget, dp, dpmin, expon, pbot, zbot, tbotextrap, tvbotextrap
    LOGICAL, INTENT(IN)                                 :: GEOPT
    INTEGER, INTENT(IN)                                 :: extrapolate
    REAL(8), INTENT(IN)                                 :: MISSING
    
    REAL(8), DIMENSION(ix, iy, num_metgrid_levels, it), INTENT(OUT)  :: data_out
    INTEGER                                          :: i, j, itt, k, kk, kin, kupper
    REAL(8), DIMENSION(num_metgrid_levels)              :: data_out1D
    REAL(8), DIMENSION(iz)                              :: data_in1D, pres_field1D

!     REAL, DIMENSION(ix, iy, num_metgrid_levels, it)  :: N
!     REAL                                             :: sumA, sumN, AVE_geopt


!    N = 1.0

    expon=287.04*.0065/9.81


    do itt = 1, it
        do j = 1, iy
            do i = 1, ix
               data_in1D(:)    = data_in(i,j,:,itt)
               pres_field1D(:) = pres_field(i,j,:,itt)
               CALL int1D (data_out1D, data_in1D, pres_field1D, interp_levels, iz, num_metgrid_levels, LINLOG, MISSING)
               data_out(i,j,:,itt) = data_out1D(:)
            end do
        end do
    end do


    ! Fill in missing values
    IF ( extrapolate == 0 ) RETURN       !! no extrapolation - we are out of here

    ! First find where about 400 hPa is located
    kk = 0
    find_kk : do k = 1, num_metgrid_levels
        kk = k
        if ( interp_levels(k) <= 40000. ) exit find_kk
    end do find_kk


    IF ( GEOPT ) THEN     !! geopt is treated different below ground

    do itt = 1, it
       do k = 1, kk
          do j = 1, iy
              do i = 1, ix
                 IF ( data_out(i,j,k,itt) == MISSING .AND. interp_levels(k) < psfc(i,j,itt) ) THEN

!                We are below the first model level, but above the ground 

                    data_out(i,j,k,itt) = ((interp_levels(k) - pres_field(i,j,1,itt))*ter(i,j)*9.81 +  &
                                           (psfc(i,j,itt) - interp_levels(k))*data_in(i,j,1,itt) ) /   &
                                          (psfc(i,j,itt) - pres_field(i,j,1,itt))

                 ELSEIF ( data_out(i,j,k,itt) == MISSING ) THEN

!                We are below both the ground and the lowest data level.

!                First, find the model level that is closest to a "target" pressure
!                level, where the "target" pressure is delta-p less that the local
!                value of a horizontally smoothed surface pressure field.  We use
!                delta-p = 150 hPa here. A standard lapse rate temperature profile
!                passing through the temperature at this model level will be used
!                to define the temperature profile below ground.  This is similar
!                to the Benjamin and Miller (1990) method, except that for
!                simplicity, they used 700 hPa everywhere for the "target" pressure.
!                Code similar to what is implemented in RIP4

                    ptarget = (psfc(i,j,itt)*.01) - 150.
                    dpmin=1.e4
                    kupper = 0
                    loop_kIN : do kin=iz,1,-1
                       kupper = kin
                       dp=abs( (pres_field(i,j,kin,itt)*.01) - ptarget )
                       if (dp.gt.dpmin) exit loop_kIN
                       dpmin=min(dpmin,dp)
                    enddo loop_kIN

                    pbot=max(pres_field(i,j,1,itt),psfc(i,j,itt))
                    zbot=min(data_in(i,j,1,itt)/9.81,ter(i,j))

                    tbotextrap=tk(i,j,kupper,itt)*(pbot/pres_field(i,j,kupper,itt))**expon
                    tvbotextrap=virtual(tbotextrap,qv(i,j,1,itt))

                    data_out(i,j,k,itt) = (zbot+tvbotextrap/.0065*(1.-(interp_levels(k)/pbot)**expon))*9.81
               
                 ENDIF
              enddo
          enddo
       enddo
    enddo


    !!! Code for filling missing data with an average - we don't want to do this
    !!do itt = 1, it
       !!loop_levels : do k = 1, num_metgrid_levels
          !!sumA = SUM(data_out(:,:,k,itt), MASK = data_out(:,:,k,itt) /= MISSING)
          !!sumN = SUM(N(:,:,k,itt), MASK = data_out(:,:,k,itt) /= MISSING)
          !!IF ( sumN == 0. ) CYCLE loop_levels
          !!AVE_geopt = sumA/sumN
          !!WHERE ( data_out(:,:,k,itt) == MISSING )
             !!data_out(:,:,k,itt) = AVE_geopt
          !!END WHERE
       !!end do loop_levels
    !!end do

    END IF

    !!! All other fields and geopt at higher levels come here
    do itt = 1, it
    do j = 1, iy
    do i = 1, ix
      do k = 1, kk
         if ( data_out(i,j,k,itt) == MISSING ) data_out(i,j,k,itt) = data_in(i,j,1,itt)
      end do
      do k = kk+1, num_metgrid_levels
         if ( data_out(i,j,k,itt) == MISSING ) data_out(i,j,k,itt) = data_in(i,j,iz,itt)
      end do
    end do
    end do
    end do

END SUBROUTINE interp
 
SUBROUTINE int1D(xxout, xxin, ppin, ppout, npin, npout, LINLOG, MISSING)
    IMPLICIT NONE
!
! Modified from int2p - NCL code
! routine to interpolate from one set of pressure levels
! .   to another set  using linear or ln(p) interpolation
!
! NCL: xout = int2p (pin,xin,pout,linlog)
! This code was originally written for a specific purpose.
! .   Several features were added for incorporation into NCL's
! .   function suite including linear extrapolation.
!
! nomenclature:
!
! .   ppin   - input pressure levels. The pin can be
! .            be in ascending or descending order
! .   xxin   - data at corresponding input pressure levels
! .   npin   - number of input pressure levels >= 2
! .   ppout  - output pressure levels (input by user)
! .            same (ascending or descending) order as pin
! .   xxout  - data at corresponding output pressure levels
! .   npout  - number of output pressure levels
! .   linlog - if abs(linlog)=1 use linear interp in pressure
! .            if abs(linlog)=2 linear interp in ln(pressure)
! .   missing- missing data code. 
!                                                ! input types
    INTEGER, INTENT(IN)     :: npin,npout,linlog
    real(8), INTENT(IN)     :: ppin(npin),xxin(npin),ppout(npout)
    real(8), INTENT(IN)     :: MISSING
    !                                                ! output
    real(8),     INTENT(OUT)  :: xxout(npout)
    INTEGER                :: np,nl,nlmax
    INTEGER                :: nlsave,nlstrt
    real(8)                :: slope,pa,pb,pc

    ! automatic arrays
    real(8)      :: p(npin),x(npin)
    real(8)    :: pout(npout),xout(npout)


    xxout = MISSING
    pout  = ppout
    p     = ppin
    x     = xxin
    nlmax = npin

    ! exact p-level matches
    nlstrt = 1
    nlsave = 1
    do np = 1,npout
        xout(np) = MISSING
        do nl = nlstrt,nlmax
            if (pout(np).eq.p(nl)) then
                xout(np) = x(nl)
                nlsave = nl + 1
                go to 10
            end if
        end do
    10     nlstrt = nlsave
    end do

    if (LINLOG.eq.1) then
        do np = 1,npout
            do nl = 1,nlmax - 1
                if (pout(np).lt.p(nl) .and. pout(np).gt.p(nl+1)) then
                  slope = (x(nl)-x(nl+1))/ (p(nl)-p(nl+1))
                  xout(np) = x(nl+1) + slope* (pout(np)-p(nl+1))
                end if
            end do
        end do
    elseif (LINLOG.eq.2) then
        do np = 1,npout
            do nl = 1,nlmax - 1
            if (pout(np).lt.p(nl) .and. pout(np).gt.p(nl+1)) then
              pa = log(p(nl))
              pb = log(pout(np))
            ! special case: in case someone inadvertently enter p=0.
              if (p(nl+1).gt.0.d0) then
                  pc = log(p(nl+1))
              else
                  pc = log(1.d-4)
              end if

              slope = (x(nl)-x(nl+1))/ (pa-pc)
              xout(np) = x(nl+1) + slope* (pb-pc)
            end if
            end do
        end do
    end if

! place results in the return array;
    xxout = xout

END SUBROUTINE int1D
!------------------------------------------------------------------------------
SUBROUTINE compute_mslp(data_out, pres_field, psfc, ter, tk, qv, ix, iy, iz, it)
! New subroutine to compute mean_sealevel pressure values 

    INTEGER, INTENT(IN)                                           :: ix, iy, iz, it
    REAL(8), DIMENSION(ix, iy, it), INTENT(OUT)                   :: data_out
    REAL(8), DIMENSION(ix, iy, iz, it), INTENT(IN)                :: pres_field, tk, qv
    REAL(8), DIMENSION(ix, iy, it), INTENT(IN)                    :: psfc
    REAL(8), DIMENSION(ix, iy), INTENT(IN)                        :: ter
    INTEGER                                                       :: i, j, itt, kin, kupper
    REAL(8)                                                       :: ptarget, dpmin, dp, tkint, tbotextrap
    REAL(8)                                                       :: tvbotextrap, pbot, expon

!    N = 1.0
    expon = 287.04*.0065/9.81

!    data_out = 0.
!    We are below both the ground and the lowest data level.

!    First, find the model level that is closest to a "target" pressure
!    level, where the "target" pressure is delta-p less that the local
!    value of a horizontally smoothed surface pressure field.  We use
!    delta-p = 150 hPa here. A standard lapse rate temperature profile
!    passing through the temperature at this model level will be used
!    to define the temperature profile below ground.  This is similar
!    to the Benjamin and Miller (1990) method, using  
!    700 hPa everywhere for the "target" pressure.
    do itt = 1, it
        do j = 1, iy
            do i = 1, ix
                ptarget = (psfc(i,j,itt)*.01) - 150.
                dpmin = 1.e4
                kupper = 0
                loop_kIN : do kin = iz,1,-1
                kupper = kin
                dp=abs( (pres_field(i, j, kin, itt)*.01) - ptarget )
                if (dp.gt.dpmin) exit loop_kIN
                dpmin=min(dpmin, dp)
            enddo loop_kIN
            ptarget=ptarget*100.
            !          
            if (pres_field(i, j, kupper + 1, itt) - ptarget .ge. 0) then
                kupper = kupper + 1
            endif

!           kupper = 8
!           ptarget = pres_field(i, j, kupper, itt)
!
!           García-Díez 2012-06
!           The reference level temperature is computed by
!           linear interpolation, so there is no jump when the selected eta level changes.
!
            tkint = (tk(i, j, kupper, itt)*abs(ptarget - pres_field(i, j, kupper, itt)) + &
            tk(i, j, kupper + 1, itt)*abs(ptarget - pres_field(i, j, kupper + & 
            1, itt)))/abs(pres_field(i, j, kupper, itt) - pres_field(i, j, kupper + 1, itt))
!           tkint = tk(i, j, kupper, itt)

            pbot = pres_field(i,j,1,itt)
            
            !         tbotextrap = tkint*(psfc(i, j, itt)/ptarget)**expon
            tkint = tk(i, j, kupper, itt)
            tbotextrap = tkint*(pbot/ptarget)**expon
            tvbotextrap = virtual(tbotextrap, qv(i,j,1,itt))
            data_out(i, j, itt) = pbot*((tvbotextrap + 0.0065*ter(i, j))/tvbotextrap)**(1/expon)

!         IF (i==INT(ix/2) .AND. j==INT(iy/2) ) THEN
!         IF (i==54 .AND. j==134) THEN
!           IF (ter(i,j) > 2300.) THEN
!             PRINT *, 'pkupper - 1:', pres_field(i,j,kupper - 1,itt), 'pkupper:', pres_field(i,j,kupper,itt), 'pkupper + 1:', pres_field(i,j,kupper + 1,itt), 'pkupper + 2:', pres_field(i,j,kupper + 2,itt)
!             PRINT *, 'tk kupper:', tk(i,j,kupper,itt), 'tk kupper + 1:', tk(i,j,kupper + 1,itt), 'tkint:', tkint
!             PRINT *, 'qv kupper:', qv(i,j,kupper,itt), 'qv kupper + 1:', qv(i,j,kupper + 1,itt), 'qvint:', qvint
!             PRINT *,itt,' ptarget',ptarget,'kupper:',kupper
!             PRINT *,'tk:',tk(i,j,kupper,itt),'psfc:',psfc(i,j,itt)
!             PRINT *,'tbot:',tbotextrap,'tvbot:',tvbotextrap,'ter:',ter(i,j)
!             PRINT *,'qv:',qv(i,j,kupper,itt),'mslp:',data_out(i,j,itt,1)
!           ENDIF
            enddo ! i
        enddo ! j
    enddo ! itt

END SUBROUTINE compute_mslp

!---------------------------------------------------------------------

SUBROUTINE clt_sundqvist(dx, dy, dz, dt, cldfra, totcloudfr)
!  Subroutine to compute total cloud cover in base 1. BY LLUIS

    IMPLICIT NONE
    INTEGER,INTENT(IN)                                     :: dx, dy, dz, dt
    REAL(8), DIMENSION(dx,dy,dz,dt), INTENT(IN)           :: cldfra
    REAL(8), DIMENSION(dx,dy,dt),    INTENT(OUT)          :: totcloudfr
    ! Local
    INTEGER                                                 :: k
    REAL(8), DIMENSION(dx,dy,dz,dt)                        :: cldfram1, cldmax

    !!!!!!!!!!!!!! Variables
    ! dx, dy, dz, dt: dimensions of fields
    ! cldfra: cloud fraction at each level
    ! totcfr: total cloud fraction
    !
    ! rcode = nf_inq_varid(ncid, cldfraname, idcldfra)
    ! rcode = nf_get_var_real(ncid, idcldfra, cldfra)

    totcloudfr = 1.
    cldfram1(:,:,2:dz,:) = cldfra(:,:,1:dz-1,:)
    cldfram1(:,:,1,:) = 0.
    cldmax = MAX(cldfram1,cldfra)

    WHERE (cldfram1 == 1.) cldfram1 = 0.99

    vertical_levels: DO k=1, dz 
        totcloudfr=totcloudfr*((1. - cldmax(:,:,k,:))/(1. - cldfram1(:,:,k,:)))
!        PRINT *, "totcloudfr(dx/, dy/2, 2):", totcloudfr(50, 50, 2)
!        PRINT *, "cldfram1(dx/, dy/2, 2):", cldfram1(50,50,k,2)
!        PRINT *, "cldmax(dx/, dy/2, 2):", cldmax(50,50,k,2)
    END DO vertical_levels
    totcloudfr = 1. - totcloudfr

END SUBROUTINE clt_sundqvist

!---------------------------------------------------------------------

SUBROUTINE clt_maxrand(tot_cldfra, cldfra, ix, iy, iz, it)
!
!  Subroutine to compute total cloud cover in base 1 using maximum-random overlapping
!
    IMPLICIT NONE
    INTEGER,INTENT(IN)                                     :: ix, iy, iz, it
    REAL(8), DIMENSION(ix,iy,iz,it), INTENT(IN)           :: cldfra
    REAL(8), DIMENSION(ix,iy,it),    INTENT(OUT)          :: tot_cldfra
    ! Local
    REAL(8)                                                 :: tot_col_cldfra
    INTEGER                                                 :: i, j, t
    !
    ! ix, iy, iz, it: dimensions of fields
    ! cldfra: cloud fraction at each model level
    ! tot_cldfra: total cloud fraction (ix, iy, it) array
    ! tot_col_cldfra: column total cloud fraction real
    ! minpres, maxpres, optional bounds for pressure levels in hPa
    ! pres_field, optional pressure field
    !
    do t = 1, it
        do j = 1, iy
            do i = 1, ix
                call clt_maxrand_column(tot_col_cldfra, cldfra(i,j,:,t), iz)
                tot_cldfra(i,j,t) = tot_col_cldfra
            end do ! i
        end do ! j
    end do ! t
    
END SUBROUTINE clt_maxrand

!---------------------------------------------------------------------

SUBROUTINE clt_maxrand_levels(tot_cldfra, cldfra, pres_field, maxpres, minpres, ix, iy, iz, it)
!
! Subroutine to compute total cloud cover in base 1 using maximum-random overlapping.
! This one is able to compute high, low and medium clouds. Its sepparated from the total
! because optional arguments are not working and then it is not efficient.
!
    IMPLICIT NONE
    INTEGER,INTENT(IN)                                     :: ix, iy, iz, it
    REAL(8), DIMENSION(ix,iy,iz,it), INTENT(IN)           :: cldfra
    REAL(8), DIMENSION(ix,iy,iz,it), INTENT(IN)           :: pres_field
    REAL(8), INTENT(IN)                                    :: minpres, maxpres
    REAL(8), DIMENSION(ix,iy,it),    INTENT(OUT)          :: tot_cldfra
    ! Local
    REAL(8), DIMENSION(ix,iy,iz,it)                        :: masked_cldfra
    REAL(8)                                                 :: tot_col_cldfra
    INTEGER                                                 :: i, j, t
    !
    ! ix, iy, iz, it: dimensions of fields
    ! cldfra: cloud fraction at each model level
    ! tot_cldfra: total cloud fraction (ix, iy, it) array
    ! tot_col_cldfra: column total cloud fraction real
    ! minpres, maxpres, bounds for pressure levels in hPa
    ! pres_field, pressure field
    !
    masked_cldfra = cldfra
    where (pres_field <= minpres*100.) masked_cldfra(:,:,:,:) = 0.
    where (pres_field >= maxpres*100.) masked_cldfra(:,:,:,:) = 0.

    do t = 1, it
        do j = 1, iy
            do i = 1, ix
                call clt_maxrand_column(tot_col_cldfra, masked_cldfra(i,j,:,t), iz)
                tot_cldfra(i,j,t) = tot_col_cldfra
            end do ! i
        end do ! j
    end do ! t
    
END SUBROUTINE clt_maxrand_levels

!---------------------------------------------------------------------

SUBROUTINE clt_maxrand_column(tot_col_cldfra, col_cldfra, iz)
!
! Subroutine for performing maximum-random overlapping in one column.
!
    INTEGER,INTENT(IN)                                     :: iz
    REAL(8), DIMENSION(iz), INTENT(IN)                    :: col_cldfra
    REAL(8), INTENT(OUT)                                   :: tot_col_cldfra
    INTEGER                                                 :: k, nseq
    REAL(8),  ALLOCATABLE, DIMENSION(:)                  :: seq_cldfra
    LOGICAL                                                 :: in_cloudy_section
    !
    ! iz (in): Number of vertical levels.
    ! col_cldfra (in): Cloud fraction in model levels.
    ! tot_col_cldfra (out) : Column total cloud fraction real
    ! k: Vertical level counter
    ! nseq: Cloudy section number (levels with cldfra > 1 sepparated by cloud free levels)
    ! seq_cldfra: Allocatable array to store each sections cloud fraction
    ! in_cloudy_section: Flag that tells us if we are inside a cloudy section or not.
    !
    ! First, we need to find the number of cloudy sections.
    !
    nseq = 0
    in_cloudy_section = .FALSE.
    do k = 1, iz
        if (col_cldfra(k) > 0.) then
            !
            ! Check if we are already in a cloudy section
            !
            if (in_cloudy_section) then
                !
                ! We are inside a section, don't count it
                !
                cycle
            else
                !
                ! We found a new cloudy section
                !
                in_cloudy_section = .TRUE.
                nseq = nseq + 1
            endif
        else
            !
            ! We are outside a cloudy section.
            !
            in_cloudy_section = .FALSE.
        endif
    end do ! k
    !debug
    !print *, "Number of cloudy sections found:", nseq
    !
    ! Allocate the vector with cloud fractions of each section
    !
    allocate(seq_cldfra(nseq))
    !
    ! Then loop again over the vertical axis storing the cloud fraction
    ! of each section, computed assuming maximum overlapping.
    !
    nseq = 0
    in_cloudy_section = .FALSE.
    do k = 1, iz
        !
        ! Check if we are already in a cloudy section
        !
        if (col_cldfra(k) > 0.) then
            if (in_cloudy_section) then
                !
                ! We are inside a section, check if this level cltfra is
                ! larger than previous and store it.
                !
                seq_cldfra(nseq) = max(seq_cldfra(nseq), col_cldfra(k))
            else
                !
                ! We found a new cloudy section. Move the counter and 
                ! save its cldfra.
                !
                in_cloudy_section = .TRUE.
                nseq = nseq + 1
                seq_cldfra(nseq) = col_cldfra(k)
            endif
        else
            !
            ! We are outside a cloudy section.
            !
            in_cloudy_section = .FALSE.
        endif
    end do ! k

    !print *, "Cloud section cldfra:", seq_cldfra
    !print *, "All levels cloud fraction:" , col_cldfra
    !
    ! Compute the total cloud cover assuming random overlapping
    ! between cloudy sections.
    !
    tot_col_cldfra = 1. - product(1. - seq_cldfra)
END SUBROUTINE clt_maxrand_column
!------------------------------------------------------------------------------
FUNCTION virtual (tmp,rmix)
    IMPLICIT NONE
!      This function returns virtual temperature in K, given temperature
!      in K and mixing ratio in kg/kg.

    real(8), intent(IN)                  :: tmp, rmix
    real(8)                              :: virtual

    virtual=tmp*(0.622+rmix)/(0.622*(1.+rmix))

END FUNCTION virtual

end MODULE
