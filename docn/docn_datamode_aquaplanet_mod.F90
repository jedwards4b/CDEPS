module docn_datamode_aquaplanet_mod

  use ESMF             , only : ESMF_SUCCESS, ESMF_State, ESMF_Mesh, ESMF_MeshGet, ESMF_LogWrite, ESMF_LOGMSG_INFO
  use NUOPC            , only : NUOPC_Advertise
  use shr_const_mod    , only : shr_const_TkFrz, shr_const_pi
  use shr_sys_mod      , only : shr_sys_abort
  use dshr_methods_mod , only : dshr_state_getfldptr, dshr_fldbun_getfldptr, chkerr
  use dshr_methods_mod , only : cdeps_real_kind, i8, cl, cs
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_add

  implicit none
  private ! except

  public :: docn_datamode_aquaplanet_advertise
  public :: docn_datamode_aquaplanet_init_pointers
  public :: docn_datamode_aquaplanet_advance

  ! export fields
  real(cdeps_real_kind), pointer    :: So_omask(:)  => null()    ! real ocean fraction sent to mediator
  real(cdeps_real_kind), pointer    :: So_t(:)      => null()
  real(cdeps_real_kind), pointer    :: So_u(:)      => null()
  real(cdeps_real_kind), pointer    :: So_v(:)      => null()

  ! model mesh lats and lons in radians
  real(cdeps_real_kind), pointer    :: rlon(:), rlat(:)

  ! parameters
  real(cdeps_real_kind) , parameter :: tkfrz = shr_const_tkfrz ! freezing point, fresh water (kelvin)
  real(cdeps_real_kind) , parameter :: pi    = shr_const_pi    ! 3.14....
  real(cdeps_real_kind) , parameter :: pio180 = shr_const_pi/180._cdeps_real_kind

  ! parameters for zonally symmetric experiments
  real(cdeps_real_kind) , parameter :: t0_max     = 27._cdeps_real_kind
  real(cdeps_real_kind) , parameter :: t0_min     = 0._cdeps_real_kind
  real(cdeps_real_kind) , parameter :: maxlat     = 60._cdeps_real_kind*pio180
  real(cdeps_real_kind) , parameter :: shift      = 5._cdeps_real_kind*pio180
  real(cdeps_real_kind) , parameter :: shift9     = 10._cdeps_real_kind*pio180
  real(cdeps_real_kind) , parameter :: shift10    = 15._cdeps_real_kind*pio180

  ! parameters for zonally asymmetric experiments
  real(cdeps_real_kind) , parameter :: t0_max6    = 1._cdeps_real_kind
  real(cdeps_real_kind) , parameter :: t0_max7    = 3._cdeps_real_kind
  real(cdeps_real_kind) , parameter :: latcen     = 0._cdeps_real_kind*pio180
  real(cdeps_real_kind) , parameter :: loncen     = 0._cdeps_real_kind*pio180
  real(cdeps_real_kind) , parameter :: latrad6    = 15._cdeps_real_kind*pio180
  real(cdeps_real_kind) , parameter :: latrad8    = 30._cdeps_real_kind*pio180
  real(cdeps_real_kind) , parameter :: lonrad     = 30._cdeps_real_kind*pio180

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine docn_datamode_aquaplanet_advertise(exportState, fldsexport, flds_scalar_name, rc)

    ! input/output variables
    type(esmf_State)   , intent(inout) :: exportState
    type(fldlist_type) , pointer       :: fldsexport
    character(len=*)   , intent(in)    :: flds_scalar_name
    integer            , intent(out)   :: rc

    ! local variables
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport, 'So_omask'            )
    call dshr_fldList_add(fldsExport, 'So_t'                )
    call dshr_fldList_add(fldsExport, 'So_u'                )
    call dshr_fldList_add(fldsExport, 'So_v'                )

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(docn_comp_advertise): Fr_ocn'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

  end subroutine docn_datamode_aquaplanet_advertise

  !===============================================================================
  subroutine docn_datamode_aquaplanet_init_pointers(exportState, ocn_fraction, rc)

    ! input/output variables
    type(ESMF_State) , intent(inout) :: exportState
    real(cdeps_real_kind)         , intent(in)    :: ocn_fraction(:)
    integer          , intent(out)   :: rc

    ! local variables
    character(len=*), parameter :: subname='(docn_init_pointers): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call dshr_state_getfldptr(exportState, 'So_omask' , fldptr1=So_omask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'So_t'     , fldptr1=So_t, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'So_u'     , fldptr1=So_u, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'So_v'     , fldptr1=So_v, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    So_u(:)   = 0.0_cdeps_real_kind
    So_v(:)   = 0.0_cdeps_real_kind

    ! Set export state ocean fraction (So_omask)
    So_omask(:) = ocn_fraction(:)

  end subroutine docn_datamode_aquaplanet_init_pointers

  !===============================================================================
  subroutine docn_datamode_aquaplanet_advance(exportState, model_mesh, sst_option, sst_constant_value, rc)

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: exportState
    type(ESMF_Mesh)        , intent(in)    :: model_mesh
    integer , optional     , intent(in)    :: sst_option
    real(cdeps_real_kind), optional     , intent(in)    :: sst_constant_value
    integer                , intent(out)   :: rc

    ! local variables
    logical           :: first_time = .true.
    integer           :: lsize
    integer           :: n,i
    integer           :: spatialDim         ! number of dimension in mesh
    integer           :: numOwnedElements   ! size of mesh
    real(cdeps_real_kind)          :: tmp, tmp1
    real(cdeps_real_kind), pointer :: ownedElemCoords(:) ! mesh lat and lons
    character(len=*), parameter :: subname='(docn_datamode_aquaplanet): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (first_time) then
       ! set model lats and lons
       call ESMF_MeshGet(model_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(ownedElemCoords(spatialDim*numOwnedElements))
       allocate(rlon(numOwnedElements), rlat(numOwnedElements))
       call ESMF_MeshGet(model_mesh, ownedElemCoords=ownedElemCoords)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,numOwnedElements
          rlon(n) = ownedElemCoords(2*n-1) * pio180
          rlat(n) = ownedElemCoords(2*n) * pio180
       end do
       deallocate(ownedElemCoords)

       first_time = .false.
    end if

    ! Determine So_t

    if (present(sst_constant_value)) then

       So_t(:) = sst_constant_value

    else if (present(sst_option)) then

       ! Determine local size
       lsize = size(rlon)

       ! Error checks
       if (sst_option < 1 .or. sst_option > 10) then
          call shr_sys_abort ('docn_prescribed_sst: ERROR: sst_option must be between 1 and 10')
       end if

       ! Determine So_t in degrees C
       if (sst_option == 1 .or. sst_option == 6 .or. sst_option == 7 .or. sst_option == 8) then
          do i = 1,lsize
             if (abs(rlat(i)) > maxlat) then
                So_t(i) = t0_min
             else
                tmp = sin(rlat(i)*pi*0.5_cdeps_real_kind/maxlat)
                tmp = 1._cdeps_real_kind - tmp*tmp
                So_t(i) = tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end if
       if (sst_option == 2) then ! Flat
          do i = 1,lsize
             if (abs(rlat(i)) > maxlat) then
                So_t(i) = t0_min
             else
                tmp = sin(rlat(i)*pi*0.5_cdeps_real_kind/maxlat)
                tmp = 1._cdeps_real_kind - tmp*tmp*tmp*tmp
                So_t(i) = tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end if
       if (sst_option == 3) then ! Qobs
          do i = 1,lsize
             if (abs(rlat(i)) > maxlat) then
                So_t(i) = t0_min
             else
                tmp = sin(rlat(i)*pi*0.5_cdeps_real_kind/maxlat)
                tmp = (2._cdeps_real_kind - tmp*tmp*tmp*tmp - tmp*tmp)*0.5_cdeps_real_kind
                So_t(i) = tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end if
       if (sst_option == 4) then ! Peaked
          do i = 1,lsize
             if (abs(rlat(i)) > maxlat) then
                So_t(i) = t0_min
             else
                tmp = (maxlat - abs(rlat(i)))/maxlat
                tmp1 = 1._cdeps_real_kind - tmp
                So_t(i) = t0_max*tmp + t0_min*tmp1
             end if
          end do
       end if
       if (sst_option == 5) then ! Control-5N
          do i = 1,lsize
             if (abs(rlat(i)) > maxlat) then
                So_t(i) = t0_min
             else if (rlat(i) > shift) then
                tmp = sin((rlat(i)-shift)*pi*0.5_cdeps_real_kind/(maxlat-shift))
                tmp = 1._cdeps_real_kind - tmp*tmp
                So_t(i) = tmp*(t0_max - t0_min) + t0_min
             else
                tmp = sin((rlat(i)-shift)*pi*0.5_cdeps_real_kind/(maxlat+shift))
                tmp = 1._cdeps_real_kind - tmp*tmp
                So_t(i) = tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end if
       if (sst_option == 6) then ! 1KEQ
          do i = 1,lsize
             if (abs(rlat(i)-latcen) <= latrad6) then
                tmp1 = cos((rlat(i)-latcen)*pi*0.5_cdeps_real_kind/latrad6)
                tmp1 = tmp1*tmp1
                tmp = abs(rlon(i)-loncen)
                tmp = min(tmp , 2._cdeps_real_kind*pi-tmp)
                if(tmp <= lonrad) then
                   tmp = cos(tmp*pi*0.5_cdeps_real_kind/lonrad)
                   tmp = tmp*tmp
                   So_t(i) = So_t(i) + t0_max6*tmp*tmp1
                end if
             end if
          end do
       end if
       if (sst_option == 7) then ! 3KEQ
          do i = 1, lsize
             if (abs(rlat(i)-latcen) <= latrad6) then
                tmp1 = cos((rlat(i)-latcen)*pi*0.5_cdeps_real_kind/latrad6)
                tmp1 = tmp1*tmp1
                tmp = abs(rlon(i)-loncen)
                tmp = min(tmp , 2._cdeps_real_kind*pi-tmp)
                if (tmp <= lonrad) then
                   tmp = cos(tmp*pi*0.5_cdeps_real_kind/lonrad)
                   tmp = tmp*tmp
                   So_t(i) = So_t(i) + t0_max7*tmp*tmp1
                end if
             end if
          end do
       end if
       if (sst_option == 8) then ! 3KW1
          do i = 1, lsize
             if (abs(rlat(i)-latcen) <= latrad8) then
                tmp1 = cos((rlat(i)-latcen)*pi*0.5_cdeps_real_kind/latrad8)
                tmp1 = tmp1*tmp1
                tmp = cos(rlon(i)-loncen)
                So_t(i) = So_t(i) + t0_max7*tmp*tmp1
             end if
          end do
       end if
       if (sst_option == 9) then ! Control-10N
          do i = 1, lsize
             if (abs(rlat(i)) > maxlat) then
                So_t(i) = t0_min
             else if (rlat(i) > shift9) then
                tmp = sin((rlat(i)-shift9)*pi*0.5_cdeps_real_kind/(maxlat-shift9))
                tmp = 1._cdeps_real_kind - tmp*tmp
                So_t(i) = tmp*(t0_max - t0_min) + t0_min
             else
                tmp = sin((rlat(i)-shift9)*pi*0.5_cdeps_real_kind/(maxlat+shift9))
                tmp = 1._cdeps_real_kind - tmp*tmp
                So_t(i) = tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end if
       if (sst_option == 10) then ! Control-15N
          do i = 1, lsize
             if (abs(rlat(i)) > maxlat) then
                So_t(i) = t0_min
             else if(rlat(i) > shift10) then
                tmp = sin((rlat(i)-shift10)*pi*0.5_cdeps_real_kind/(maxlat-shift10))
                tmp = 1._cdeps_real_kind - tmp*tmp
                So_t(i) = tmp*(t0_max - t0_min) + t0_min
             else
                tmp = sin((rlat(i)-shift10)*pi*0.5_cdeps_real_kind/(maxlat+shift10))
                tmp = 1._cdeps_real_kind - tmp*tmp
                So_t(i) = tmp*(t0_max - t0_min) + t0_min
             end if
          end do
       end if

       ! Convert to kelvin
       So_t(:) = So_t(:) + TkFrz

    else
       call shr_sys_abort("ERROR: either sst_constant value or sst_option must be input")
    end if

  end subroutine docn_datamode_aquaplanet_advance

end module docn_datamode_aquaplanet_mod
