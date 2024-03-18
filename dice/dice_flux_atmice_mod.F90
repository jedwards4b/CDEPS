module dice_flux_atmice_mod

    ! shared constants
  use dshr_methods_mod, only : cdeps_real_kind, cxx, cl, cs
  use shr_const_mod, only : loc_zvir => shr_const_zvir, loc_cpdair => shr_const_cpdair
  use shr_const_mod, only : loc_cpvir => shr_const_cpvir, loc_karman => shr_const_karman
  use shr_const_mod, only : loc_g => shr_const_g, loc_latvap => shr_const_latvap
  use shr_const_mod, only : loc_latice => shr_const_latice, loc_stebol => shr_const_stebol
  use shr_const_mod, only : spval => shr_const_spval
  implicit none

  integer,parameter :: dbug = 0 ! internal debug level

!===============================================================================
contains
!===============================================================================

  subroutine dice_flux_atmice(             &
       mask  ,zbot  ,ubot  ,vbot  ,thbot,  &
       qbot  ,rbot  ,tbot  ,ts    ,sen,    &
       lat   ,lwup  ,evap  ,taux  ,tauy,   &
       tref  ,qref  ,logunit               )

    !-------------------------------------------------------------------------------
    ! PURPOSE:
    !   using atm & ice state variables, compute atm/ice fluxes
    !   and diagnostic 10m air temperature and humidity
    !
    ! NOTE:
    !   o all fluxes are positive downward
    !   o net heat flux = net sw + lw up + lw down + sen + lat
    !   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
    !   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
    !
    ! ASSUME:
    !   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
    !-------------------------------------------------------------------------------

    !--- input arguments --------------------------------
    integer    ,intent(in)  :: mask (:)    ! 0 <=> cell NOT in model domain
    real(cdeps_real_kind)   ,intent(in)  :: zbot (:)    ! atm level height  (m)
    real(cdeps_real_kind)   ,intent(in)  :: ubot (:)    ! atm u wind     (m/s)
    real(cdeps_real_kind)   ,intent(in)  :: vbot (:)    ! atm v wind     (m/s)
    real(cdeps_real_kind)   ,intent(in)  :: thbot(:)    ! atm potential T   (K)
    real(cdeps_real_kind)   ,intent(in)  :: qbot (:)    ! atm specific humidity (kg/kg)
    real(cdeps_real_kind)   ,intent(in)  :: rbot (:)    ! atm air density   (kg/m^3)
    real(cdeps_real_kind)   ,intent(in)  :: tbot (:)    ! atm T       (K)
    real(cdeps_real_kind)   ,intent(in)  :: ts   (:)    ! surface temperature
    integer    ,intent(in)  :: logunit     ! logging unit number

    !--- output arguments -------------------------------
    real(cdeps_real_kind)   ,intent(out) :: sen  (:)    ! sensible      heat flux  (W/m^2)
    real(cdeps_real_kind)   ,intent(out) :: lat  (:)    ! latent        heat flux  (W/m^2)
    real(cdeps_real_kind)   ,intent(out) :: lwup (:)    ! long-wave upward heat flux  (W/m^2)
    real(cdeps_real_kind)   ,intent(out) :: evap (:)    ! evaporative water flux ((kg/s)/m^2)
    real(cdeps_real_kind)   ,intent(out) :: taux (:)    ! x surface stress (N)
    real(cdeps_real_kind)   ,intent(out) :: tauy (:)    ! y surface stress (N)
    real(cdeps_real_kind)   ,intent(out) :: tref (:)    ! 2m reference height temperature
    real(cdeps_real_kind)   ,intent(out) :: qref (:)    ! 2m reference height humidity

    !--- local constants --------------------------------
    real(cdeps_real_kind),parameter :: umin   =  1.0_cdeps_real_kind            ! minimum wind speed (m/s)
    real(cdeps_real_kind),parameter :: zref   = 10.0_cdeps_real_kind            ! ref height           ~ m
    real(cdeps_real_kind),parameter :: ztref  =  2.0_cdeps_real_kind            ! ref height for air T ~ m
    real(cdeps_real_kind),parameter :: zzsice = 0.0005_cdeps_real_kind          ! ice surface roughness

    !--- local variables --------------------------------
    integer     :: lsize  ! array dimensions
    integer     :: n      ! array indicies
    real(cdeps_real_kind)    :: vmag   ! surface wind magnitude   (m/s)
    real(cdeps_real_kind)    :: thvbot ! virtual temperature      (K)
    real(cdeps_real_kind)    :: ssq    ! sea surface humidity     (kg/kg)
    real(cdeps_real_kind)    :: dssqdt ! derivative of ssq wrt Ts (kg/kg/K)
    real(cdeps_real_kind)    :: delt   ! potential T difference   (K)
    real(cdeps_real_kind)    :: delq   ! humidity difference      (kg/kg)
    real(cdeps_real_kind)    :: stable ! stability factor
    real(cdeps_real_kind)    :: rdn    ! sqrt of neutral exchange coefficient (momentum)
    real(cdeps_real_kind)    :: rhn    ! sqrt of neutral exchange coefficient (heat)
    real(cdeps_real_kind)    :: ren    ! sqrt of neutral exchange coefficient (water)
    real(cdeps_real_kind)    :: rd     ! sqrt of exchange coefficient (momentum)
    real(cdeps_real_kind)    :: rh     ! sqrt of exchange coefficient (heat)
    real(cdeps_real_kind)    :: re     ! sqrt of exchange coefficient (water)
    real(cdeps_real_kind)    :: ustar  ! ustar
    real(cdeps_real_kind)    :: qstar  ! qstar
    real(cdeps_real_kind)    :: tstar  ! tstar
    real(cdeps_real_kind)    :: hol    ! H (at zbot) over L
    real(cdeps_real_kind)    :: xsq    ! temporary variable
    real(cdeps_real_kind)    :: xqq    ! temporary variable
    real(cdeps_real_kind)    :: psimh  ! stability function at zbot (momentum)
    real(cdeps_real_kind)    :: psixh  ! stability function at zbot (heat and water)
    real(cdeps_real_kind)    :: alz    ! ln(zbot/z10)
    real(cdeps_real_kind)    :: ltheat ! latent heat for surface
    real(cdeps_real_kind)    :: tau    ! stress at zbot
    real(cdeps_real_kind)    :: cp     ! specific heat of moist air

    real(cdeps_real_kind)    :: bn     ! exchange coef funct for interpolation
    real(cdeps_real_kind)    :: bh     ! exchange coef funct for interpolation
    real(cdeps_real_kind)    :: fac    ! interpolation factor
    real(cdeps_real_kind)    :: ln0    ! log factor for interpolation
    real(cdeps_real_kind)    :: ln3    ! log factor for interpolation

    !--- local functions --------------------------------
    real(cdeps_real_kind)   :: Tk      ! temperature (K)
    real(cdeps_real_kind)   :: qsat    ! the saturation humidity of air (kg/m^3)
    real(cdeps_real_kind)   :: dqsatdt ! derivative of qsat wrt surface temperature
    real(cdeps_real_kind)   :: xd      ! dummy argument
    real(cdeps_real_kind)   :: psimhu  ! unstable part of psimh
    real(cdeps_real_kind)   :: psixhu  ! unstable part of psimx

    qsat(Tk)    = 627572.4_cdeps_real_kind / exp(5107.4_cdeps_real_kind/Tk)
    dqsatdt(Tk) = (5107.4_cdeps_real_kind / Tk**2) * 627572.4_cdeps_real_kind / exp(5107.4_cdeps_real_kind/Tk)
    psimhu(xd)  = log((1.0_cdeps_real_kind+xd*(2.0_cdeps_real_kind+xd))*(1.0_cdeps_real_kind+xd*xd)/8.0_cdeps_real_kind) - 2.0_cdeps_real_kind*atan(xd) + 1.571_cdeps_real_kind
    psixhu(xd)  =  2.0_cdeps_real_kind * log((1.0_cdeps_real_kind + xd*xd)/2.0_cdeps_real_kind)

    !--- formats ----------------------------------------
    character(*),parameter ::    F01 = "('(dice_flux_atmIce) ',a, i7,2x,d21.14)"
    character(*),parameter :: subName =  "(dice_flux_atmIce) "
    !-------------------------------------------------------------------------------

    lsize = size(tbot)

    do n = 1,lsize

       if (mask(n) == 0) then
          sen  (n) = spval
          lat  (n) = spval
          lwup (n) = spval
          evap (n) = spval
          taux (n) = spval
          tauy (n) = spval
          tref (n) = spval
          qref (n) = spval
       else
          !--- define some needed variables ---
          vmag   = max(umin, sqrt(ubot(n)**2+vbot(n)**2))
          thvbot = thbot(n)*(1.0_cdeps_real_kind + loc_zvir * qbot(n)) ! virtual pot temp (K)
          ssq   =  qsat  (ts(n)) / rbot(n)           ! sea surf hum (kg/kg)
          dssqdt = dqsatdt(ts(n)) / rbot(n)           ! deriv of ssq wrt Ts
          delt   = thbot(n) - ts(n)                   ! pot temp diff (K)
          delq   = qbot(n) - ssq                        ! spec hum dif (kg/kg)
          alz    = log(zbot(n)/zref)
          cp     = loc_cpdair*(1.0_cdeps_real_kind + loc_cpvir*ssq)
          ltheat = loc_latvap + loc_latice

          !----------------------------------------------------------
          ! first estimate of Z/L and ustar, tstar and qstar
          !----------------------------------------------------------

          !--- neutral coefficients, z/L = 0.0 ---
          rdn = loc_karman/log(zref/zzsice)
          rhn = rdn
          ren = rdn

          !--- ustar,tstar,qstar ----
          ustar = rdn * vmag
          tstar = rhn * delt
          qstar = ren * delq

          !--- compute stability & evaluate all stability functions ---
          hol    = loc_karman * loc_g * zbot(n) &
               &     * (tstar/thvbot+qstar/(1.0_cdeps_real_kind/loc_zvir+qbot(n))) / ustar**2
          hol    = sign( min(abs(hol),10.0_cdeps_real_kind), hol )
          stable = 0.5_cdeps_real_kind + sign(0.5_cdeps_real_kind , hol)
          xsq    = max(sqrt(abs(1.0_cdeps_real_kind - 16.0_cdeps_real_kind*hol)) , 1.0_cdeps_real_kind)
          xqq    = sqrt(xsq)
          psimh  = -5.0_cdeps_real_kind*hol*stable + (1.0_cdeps_real_kind-stable)*psimhu(xqq)
          psixh  = -5.0_cdeps_real_kind*hol*stable + (1.0_cdeps_real_kind-stable)*psixhu(xqq)

          !--- shift all coeffs to measurement height and stability ---
          rd = rdn / (1.0_cdeps_real_kind+rdn/loc_karman*(alz-psimh))
          rh = rhn / (1.0_cdeps_real_kind+rhn/loc_karman*(alz-psixh))
          re = ren / (1.0_cdeps_real_kind+ren/loc_karman*(alz-psixh))

          !--- update ustar, tstar, qstar w/ updated, shifted coeffs --
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq

          !----------------------------------------------------------
          ! iterate to converge on Z/L, ustar, tstar and qstar
          !----------------------------------------------------------

          !--- compute stability & evaluate all stability functions ---
          hol    = loc_karman * loc_g * zbot(n) &
               &      * (tstar/thvbot+qstar/(1.0_cdeps_real_kind/loc_zvir+qbot(n))) / ustar**2
          hol    = sign( min(abs(hol),10.0_cdeps_real_kind), hol )
          stable = 0.5_cdeps_real_kind + sign(0.5_cdeps_real_kind , hol)
          xsq    = max(sqrt(abs(1.0_cdeps_real_kind - 16.0_cdeps_real_kind*hol)) , 1.0_cdeps_real_kind)
          xqq    = sqrt(xsq)
          psimh  = -5.0_cdeps_real_kind*hol*stable + (1.0_cdeps_real_kind-stable)*psimhu(xqq)
          psixh  = -5.0_cdeps_real_kind*hol*stable + (1.0_cdeps_real_kind-stable)*psixhu(xqq)

          !--- shift all coeffs to measurement height and stability ---
          rd = rdn / (1.0_cdeps_real_kind+rdn/loc_karman*(alz-psimh))
          rh = rhn / (1.0_cdeps_real_kind+rhn/loc_karman*(alz-psixh))
          re = ren / (1.0_cdeps_real_kind+ren/loc_karman*(alz-psixh))

          !--- update ustar, tstar, qstar w/ updated, shifted coeffs --
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq

          !----------------------------------------------------------
          ! compute the fluxes
          !----------------------------------------------------------

          tau = rbot(n) * ustar * ustar

          !--- momentum flux ---
          taux(n) = tau * ubot(n) / vmag
          tauy(n) = tau * vbot(n) / vmag

          !--- heat flux ---
          sen (n) =   cp * tau * tstar / ustar
          lat (n) =  ltheat * tau * qstar / ustar
          lwup(n) = -loc_stebol * ts(n)**4

          !--- water flux ---
          evap(n) = lat(n)/ltheat

          !----------------------------------------------------------
          ! compute diagnostic: 2m reference height temperature
          !----------------------------------------------------------

          !--- Compute function of exchange coefficients. Assume that
          !--- cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore
          !--- 1/sqrt(cn(n))=1/rdn and sqrt(cm(n))/ch(n)=1/rh
          bn = loc_karman/rdn
          bh = loc_karman/rh

          !--- Interpolation factor for stable and unstable cases
          ln0 = log(1.0_cdeps_real_kind + (ztref/zbot(n))*(exp(bn) - 1.0_cdeps_real_kind))
          ln3 = log(1.0_cdeps_real_kind + (ztref/zbot(n))*(exp(bn - bh) - 1.0_cdeps_real_kind))
          fac = (ln0 - ztref/zbot(n)*(bn - bh))/bh * stable &
               &   + (ln0 - ln3)/bh * (1.0_cdeps_real_kind-stable)
          fac = min(max(fac,0.0_cdeps_real_kind),1.0_cdeps_real_kind)

          !--- actual interpolation
          tref(n) = ts(n) + (tbot(n) - ts(n))*fac
          qref(n) = qbot(n) - delq*fac

       endif
    enddo

    if (dbug > 0) then
       do n = 1,lsize
          if (mask(n) /= 0) then
             write(logunit, F01)'n,mask  = ',n,mask(n)
             write(logunit, F01)'n,zbot  = ',n,zbot(n)
             write(logunit, F01)'n,ubot  = ',n,ubot(n)
             write(logunit, F01)'n,vbot  = ',n,vbot(n)
             write(logunit, F01)'n,thbot = ',n,thbot(n)
             write(logunit, F01)'n,qbot  = ',n,qbot(n)
             write(logunit, F01)'n,tbot  = ',n,tbot(n)
             write(logunit, F01)'n,ts    = ',n,ts(n)
             write(logunit, F01)'n,lat   = ',n,lat(n)
             write(logunit, F01)'n,sen   = ',n,sen(n)
             write(logunit, F01)'n,taux  = ',n,taux(n)
             write(logunit, F01)'n,taux  = ',n,tauy(n)
             write(logunit, F01)'n,lwup  = ',n,lwup(n)
             write(logunit, F01)'n,evap  = ',n,evap(n)
             write(logunit, F01)'n,tref  = ',n,tref(n)
             write(logunit, F01)'n,qref  = ',n,qref(n)
          end if
       end do
    end if

  end subroutine dice_flux_atmIce

end module dice_flux_atmice_mod
