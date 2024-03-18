module dice_datamode_ssmi_mod

  use ESMF                 , only : ESMF_State, ESMF_LogWrite, ESMF_Array, ESMF_MeshGet
  use ESMF                 , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_DistGrid
  use ESMF                 , only : ESMF_ArrayCreate, ESMF_ArrayDestroy
  use NUOPC                , only : NUOPC_Advertise
  use shr_sys_mod          , only : shr_sys_abort
  use shr_const_mod        , only : shr_const_pi, shr_const_spval, shr_const_tkfrz, shr_const_latice
  use shr_frz_mod          , only : shr_frz_freezetemp
  use dshr_strdata_mod     , only : shr_strdata_get_stream_pointer, shr_strdata_type
  use dshr_methods_mod     , only : dshr_state_getfldptr, dshr_fldbun_getfldptr, dshr_fldbun_regrid, chkerr
  use dshr_methods_mod     , only : cdeps_real_kind, i8, cs, cl
  use dshr_mod             , only : dshr_restart_read, dshr_restart_write
  use dice_flux_atmice_mod , only : dice_flux_atmice
  use dshr_fldlist_mod     , only : fldlist_type, dshr_fldlist_add

  implicit none
  private ! except

  public  :: dice_datamode_ssmi_advertise
  public  :: dice_datamode_ssmi_init_pointers
  public  :: dice_datamode_ssmi_advance
  public  :: dice_datamode_ssmi_restart_write
  public  :: dice_datamode_ssmi_restart_read

  ! restart fields
  real(cdeps_real_kind), pointer, public :: water(:) => null()

  ! internal fields
  real(cdeps_real_kind), pointer :: yc(:)      => null() ! mesh lats (degrees)
  integer , pointer :: imask(:)   => null()
  !real(cdeps_real_kind), pointer:: ifrac0(:)  => null()

  ! export fields
  real(cdeps_real_kind), pointer ::  Si_imask(:)      => null()
  real(cdeps_real_kind), pointer ::  Si_ifrac(:)      => null()
  real(cdeps_real_kind), pointer ::  Si_t(:)          => null()
  real(cdeps_real_kind), pointer ::  Si_tref(:)       => null()
  real(cdeps_real_kind), pointer ::  Si_qref(:)       => null()
  real(cdeps_real_kind), pointer ::  Si_avsdr(:)      => null()
  real(cdeps_real_kind), pointer ::  Si_anidr(:)      => null()
  real(cdeps_real_kind), pointer ::  Si_avsdf(:)      => null()
  real(cdeps_real_kind), pointer ::  Si_anidf(:)      => null()
  real(cdeps_real_kind), pointer ::  Faii_swnet(:)    => null()
  real(cdeps_real_kind), pointer ::  Faii_sen(:)      => null()
  real(cdeps_real_kind), pointer ::  Faii_lat(:)      => null()
  real(cdeps_real_kind), pointer ::  Faii_lwup(:)     => null()
  real(cdeps_real_kind), pointer ::  Faii_evap(:)     => null()
  real(cdeps_real_kind), pointer ::  Faii_taux(:)     => null()
  real(cdeps_real_kind), pointer ::  Faii_tauy(:)     => null()
  real(cdeps_real_kind), pointer ::  Fioi_melth(:)    => null()
  real(cdeps_real_kind), pointer ::  Fioi_meltw(:)    => null()
  real(cdeps_real_kind), pointer ::  Fioi_swpen(:)    => null()
  real(cdeps_real_kind), pointer ::  Fioi_taux(:)     => null()
  real(cdeps_real_kind), pointer ::  Fioi_tauy(:)     => null()
  real(cdeps_real_kind), pointer ::  Fioi_salt(:)     => null()
  real(cdeps_real_kind), pointer ::  Fioi_bcpho(:)    => null()
  real(cdeps_real_kind), pointer ::  Fioi_bcphi(:)    => null()
  real(cdeps_real_kind), pointer ::  Fioi_flxdst(:)   => null()
  real(cdeps_real_kind), pointer ::  Si_ifrac_n(:,:)  => null()
  real(cdeps_real_kind), pointer ::  Fioi_swpen_ifrac_n(:,:) => null()

  ! import fields
  real(cdeps_real_kind), pointer :: Faxa_swvdr(:)    => null()
  real(cdeps_real_kind), pointer :: Faxa_swvdf(:)    => null()
  real(cdeps_real_kind), pointer :: Faxa_swndr(:)    => null()
  real(cdeps_real_kind), pointer :: Faxa_swndf(:)    => null()
  real(cdeps_real_kind), pointer :: Fioo_q(:)        => null()
  real(cdeps_real_kind), pointer :: Sa_z(:)          => null()
  real(cdeps_real_kind), pointer :: Sa_u(:)          => null()
  real(cdeps_real_kind), pointer :: Sa_v(:)          => null()
  real(cdeps_real_kind), pointer :: Sa_ptem(:)       => null()
  real(cdeps_real_kind), pointer :: Sa_shum(:)       => null()
  real(cdeps_real_kind), pointer :: Sa_dens(:)       => null()
  real(cdeps_real_kind), pointer :: Sa_tbot(:)       => null()
  real(cdeps_real_kind), pointer :: So_s(:)          => null()
  real(cdeps_real_kind), pointer :: Faxa_bcph(:,:)   => null()
  real(cdeps_real_kind), pointer :: Faxa_ocph(:,:)   => null()
  real(cdeps_real_kind), pointer :: Faxa_dstdry(:,:) => null()
  real(cdeps_real_kind), pointer :: Faxa_dstwet(:,:) => null()

  ! surface albedo constants
  real(cdeps_real_kind) , parameter :: snwfrac = 0.286_cdeps_real_kind ! snow cover fraction ~ [0,1]
  real(cdeps_real_kind) , parameter :: as_nidf = 0.950_cdeps_real_kind ! albedo: snow,near-infr,diffuse
  real(cdeps_real_kind) , parameter :: as_vsdf = 0.700_cdeps_real_kind ! albedo: snow,visible  ,diffuse
  real(cdeps_real_kind) , parameter :: as_nidr = 0.960_cdeps_real_kind ! albedo: snow,near-infr,direct
  real(cdeps_real_kind) , parameter :: as_vsdr = 0.800_cdeps_real_kind ! albedo: snow,visible  ,direct
  real(cdeps_real_kind) , parameter :: ai_nidf = 0.700_cdeps_real_kind ! albedo: ice, near-infr,diffuse
  real(cdeps_real_kind) , parameter :: ai_vsdf = 0.500_cdeps_real_kind ! albedo: ice, visible  ,diffuse
  real(cdeps_real_kind) , parameter :: ai_nidr = 0.700_cdeps_real_kind ! albedo: ice, near-infr,direct
  real(cdeps_real_kind) , parameter :: ai_vsdr = 0.500_cdeps_real_kind ! albedo: ice, visible  ,direct
  real(cdeps_real_kind) , parameter :: ax_nidf = ai_nidf*(1.0_cdeps_real_kind-snwfrac) + as_nidf*snwfrac
  real(cdeps_real_kind) , parameter :: ax_vsdf = ai_vsdf*(1.0_cdeps_real_kind-snwfrac) + as_vsdf*snwfrac
  real(cdeps_real_kind) , parameter :: ax_nidr = ai_nidr*(1.0_cdeps_real_kind-snwfrac) + as_nidr*snwfrac
  real(cdeps_real_kind) , parameter :: ax_vsdr = ai_vsdr*(1.0_cdeps_real_kind-snwfrac) + as_vsdr*snwfrac

  ! other parameters
  real(cdeps_real_kind) , parameter :: pi       = shr_const_pi     ! pi
  real(cdeps_real_kind) , parameter :: spval    = shr_const_spval  ! flags invalid data
  real(cdeps_real_kind) , parameter :: tFrz     = shr_const_tkfrz  ! temp of freezing
  real(cdeps_real_kind) , parameter :: latice   = shr_const_latice ! latent heat of fusion
  real(cdeps_real_kind) , parameter :: waterMax = 1000.0_cdeps_real_kind        ! wrt iFrac comp & frazil ice (kg/m^2)

  character(*) , parameter :: nullstr = 'null'
  character(*) , parameter :: rpfile  = 'rpointer.ice'
  character(*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dice_datamode_ssmi_advertise(importState, exportState, fldsimport, fldsexport, &
       flds_scalar_name, flds_i2o_per_cat, rc)

    ! --------------------------------------------------------------
    ! determine export and import fields to advertise to mediator
    ! --------------------------------------------------------------

    ! input/output arguments
    type(ESMF_State)   , intent(inout) :: importState
    type(ESMF_State)   , intent(inout) :: exportState
    type(fldlist_type) , pointer       :: fldsimport
    type(fldlist_type) , pointer       :: fldsexport
    character(len=*)   , intent(in)    :: flds_scalar_name
    logical            , intent(in)    :: flds_i2o_per_cat ! .true. if select per ice thickness
    integer            , intent(out)   :: rc

    ! local variables
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    ! Advertise export fields
    call dshr_fldList_add(fldsExport , trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport ,'Si_ifrac'    )
    call dshr_fldList_add(fldsExport ,'Si_imask'    )
    call dshr_fldList_add(fldsExport ,'Si_t'        )
    call dshr_fldList_add(fldsExport ,'Si_tref'     )
    call dshr_fldList_add(fldsExport ,'Si_qref'     )
    call dshr_fldList_add(fldsExport ,'Si_avsdr'    )
    call dshr_fldList_add(fldsExport ,'Si_anidr'    )
    call dshr_fldList_add(fldsExport ,'Si_avsdf'    )
    call dshr_fldList_add(fldsExport ,'Si_anidf'    )
    call dshr_fldList_add(fldsExport ,'Faii_swnet'  )
    call dshr_fldList_add(fldsExport ,'Faii_sen'    )
    call dshr_fldList_add(fldsExport ,'Faii_lat'    )
    call dshr_fldList_add(fldsExport ,'Faii_lwup'   )
    call dshr_fldList_add(fldsExport ,'Faii_evap'   )
    call dshr_fldList_add(fldsExport ,'Faii_taux'   )
    call dshr_fldList_add(fldsExport ,'Faii_tauy'   )
    call dshr_fldList_add(fldsExport ,'Fioi_melth'  )
    call dshr_fldList_add(fldsExport ,'Fioi_meltw'  )
    call dshr_fldList_add(fldsExport ,'Fioi_swpen'  )
    call dshr_fldList_add(fldsExport ,'Fioi_taux'   )
    call dshr_fldList_add(fldsExport ,'Fioi_tauy'   )
    call dshr_fldList_add(fldsExport ,'Fioi_salt'   )
    call dshr_fldList_add(fldsExport ,'Fioi_bcpho'  )
    call dshr_fldList_add(fldsExport ,'Fioi_bcphi'  )
    call dshr_fldList_add(fldsExport ,'Fioi_flxdst' )
    if (flds_i2o_per_cat) then
       call dshr_fldList_add(fldsExport, 'Si_ifrac_n'        , ungridded_lbound=1, ungridded_ubound=1)
       call dshr_fldList_add(fldsExport, 'Fioi_swpen_ifrac_n', ungridded_lbound=1, ungridded_ubound=1)
    end if

    ! Advertise import fields
    call dshr_fldList_add(fldsImport , trim(flds_scalar_name))
    call dshr_fldList_add(fldsImport, 'Faxa_swvdr' )
    call dshr_fldList_add(fldsImport, 'Faxa_swvdf' )
    call dshr_fldList_add(fldsImport, 'Faxa_swndr' )
    call dshr_fldList_add(fldsImport, 'Faxa_swndf' )
    call dshr_fldList_add(fldsImport, 'Fioo_q'     )
    call dshr_fldList_add(fldsImport, 'Sa_z'       )
    call dshr_fldList_add(fldsImport, 'Sa_u'       )
    call dshr_fldList_add(fldsImport, 'Sa_v'       )
    call dshr_fldList_add(fldsImport, 'Sa_ptem'    )
    call dshr_fldList_add(fldsImport, 'Sa_shum'    )
    call dshr_fldList_add(fldsImport, 'Sa_dens'    )
    call dshr_fldList_add(fldsImport, 'Sa_tbot'    )
    call dshr_fldList_add(fldsImport, 'So_s'       )
    call dshr_fldList_add(fldsImport, 'Faxa_bcph'  , ungridded_lbound=1, ungridded_ubound=3)
    call dshr_fldList_add(fldsImport, 'Faxa_ocph'  , ungridded_lbound=1, ungridded_ubound=3)
    call dshr_fldList_add(fldsImport, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)
    call dshr_fldList_add(fldsImport, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(dice_comp_advertise): Fr_ice'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    fldlist => fldsImport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(importState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(dice_comp_advertise): Fr_ice'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

  end subroutine dice_datamode_ssmi_advertise

  !===============================================================================
  subroutine dice_datamode_ssmi_init_pointers(importState, exportState, sdat, flds_i2o_per_cat, rc)

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    type(shr_strdata_type) , intent(in)    :: sdat
    logical                , intent(in)    :: flds_i2o_per_cat ! .true. if select per ice thickness
    integer                , intent(out)   :: rc

    ! local variables
    integer             :: n
    integer             :: lsize
    type(ESMF_DistGrid) :: distGrid
    type(ESMF_Array)    :: elemMaskArray
    integer             :: spatialDim         ! number of dimension in mesh
    integer             :: numOwnedElements   ! size of mesh
    real(cdeps_real_kind), pointer   :: ownedElemCoords(:) ! mesh lat and lons
    character(len=*), parameter :: subname='(dice_init_pointers): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lsize = sdat%model_lsize

    ! Set Si_imask (this corresponds to the ocean mask)
    call dshr_state_getfldptr(exportState, fldname='Si_imask'    , fldptr1=Si_imask    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(imask(sdat%model_lsize))
    call ESMF_MeshGet(sdat%model_mesh, numOwnedElements=numOwnedElements, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    elemMaskArray = ESMF_ArrayCreate(distGrid, imask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(sdat%model_mesh, elemMaskArray=elemMaskArray, rc=rc) ! set the varues of imask
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    Si_imask(:) = real(imask(:), kind=cdeps_real_kind) ! set the mask as real
    call ESMF_ArrayDestroy(elemMaskArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set export state module pointers
    call dshr_state_getfldptr(exportState, fldname='Si_ifrac'    , fldptr1=Si_ifrac    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_t'        , fldptr1=Si_t        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_tref'     , fldptr1=Si_tref     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_qref'     , fldptr1=Si_qref     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_avsdr'    , fldptr1=Si_avsdr    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_anidr'    , fldptr1=Si_anidr    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_avsdf'    , fldptr1=Si_avsdf    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_anidf'    , fldptr1=Si_anidf    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_swnet'  , fldptr1=Faii_swnet  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_sen'    , fldptr1=Faii_sen    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_lat'    , fldptr1=Faii_lat    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_lwup'   , fldptr1=Faii_lwup   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_evap'   , fldptr1=Faii_evap   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_taux'   , fldptr1=Faii_taux   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_tauy'   , fldptr1=Faii_tauy   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_melth'  , fldptr1=Fioi_melth  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_meltw'  , fldptr1=Fioi_meltw  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_swpen'  , fldptr1=Fioi_swpen  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_taux'   , fldptr1=Fioi_taux   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_tauy'   , fldptr1=Fioi_tauy   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_salt'   , fldptr1=Fioi_salt   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_bcpho'  , fldptr1=Fioi_bcpho  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_bcphi'  , fldptr1=Fioi_bcphi  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_flxdst' , fldptr1=Fioi_flxdst , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_i2o_per_cat) then
       call dshr_state_getfldptr(exportState, fldname='Si_ifrac_n', fldptr2=Si_ifrac_n, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(exportState, fldname='Fioi_swpen_ifrac_n', fldptr2=Fioi_swpen_ifrac_n, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Set pointers to importState fields
    call dshr_state_getfldptr(importState, fldname='Faxa_swvdr'  , fldptr1=Faxa_swvdr  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_swvdf'  , fldptr1=Faxa_swvdf  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_swndr'  , fldptr1=Faxa_swndr  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_swndf'  , fldptr1=Faxa_swndf  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_bcph'   , fldptr2=Faxa_bcph   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_ocph'   , fldptr2=Faxa_ocph   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_dstdry' , fldptr2=Faxa_dstdry , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_dstwet' , fldptr2=Faxa_dstwet , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Fioo_q'      , fldptr1=Fioo_q      , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_z'        , fldptr1=Sa_z        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_u'        , fldptr1=Sa_u        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_v'        , fldptr1=Sa_v        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_ptem'     , fldptr1=Sa_ptem     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_tbot'     , fldptr1=Sa_tbot     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_shum'     , fldptr1=Sa_shum     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_dens'     , fldptr1=Sa_dens     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='So_s'        , fldptr1=So_s        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! initialize import arrays
    ! used for the first use in generating the export state and should have no impact on the solution
    Faxa_swvdr(:)    = 0._cdeps_real_kind
    Faxa_swvdf(:)    = 0._cdeps_real_kind
    Faxa_swndr(:)    = 0._cdeps_real_kind
    Faxa_swndf(:)    = 0._cdeps_real_kind
    Faxa_bcph(:,:)   = 0._cdeps_real_kind
    Faxa_ocph(:,:)   = 0._cdeps_real_kind
    Faxa_dstdry(:,:) = 0._cdeps_real_kind
    Faxa_dstwet(:,:) = 0._cdeps_real_kind
    Fioo_q(:)        = 0._cdeps_real_kind
    Sa_z(:)          = 10.0_cdeps_real_kind
    Sa_u(:)          = 5.0_cdeps_real_kind
    Sa_v(:)          = 5.0_cdeps_real_kind
    Sa_ptem(:)       = 260.0_cdeps_real_kind
    Sa_tbot(:)       = 260.0_cdeps_real_kind
    Sa_shum(:)       = 0.0014_cdeps_real_kind
    Sa_dens(:)       = 1.3_cdeps_real_kind
    So_s(:)          = 0._cdeps_real_kind

    ! Determine water from restart if needed
    ! allocate module arrays that are not part of import and export state

    ! Determine model latitudes
    allocate(yc(lsize))
    call ESMF_MeshGet(sdat%model_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    call ESMF_MeshGet(sdat%model_mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       yc(n) = ownedElemCoords(2*n)
    end do
    deallocate(ownedElemCoords)

  end subroutine dice_datamode_ssmi_init_pointers

  !===============================================================================
  subroutine dice_datamode_ssmi_advance(exportState, importState, cosarg, flds_i2o_per_cat, &
       flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0, dt, logunit, restart_read, rc)

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: exportState
    type(ESMF_State)       , intent(inout) :: importState
    real(cdeps_real_kind)               , intent(in)    :: cosarg     ! for setting ice temp pattern
    logical                , intent(in)    :: flds_i2o_per_cat
    real(cdeps_real_kind)               , intent(in)    :: flux_swpf  ! short-wave penatration factor
    real(cdeps_real_kind)               , intent(in)    :: flux_Qmin  ! bound on melt rate
    logical                , intent(in)    :: flux_Qacc  ! activates water accumulation/melt wrt Q
    real(cdeps_real_kind)               , intent(in)    :: flux_Qacc0 ! initial water accumulation value
    real(cdeps_real_kind)               , intent(in)    :: dt
    integer                , intent(in)    :: logunit
    logical                , intent(in)    :: restart_read
    integer                , intent(out)   :: rc

    ! local variables
    logical               :: first_time = .true.
    integer               :: n
    integer               :: lsize
    real(cdeps_real_kind)              :: qmeltall ! q that would melt all accumulated water
    real(cdeps_real_kind), allocatable :: tfreeze(:)
    character(len=*), parameter :: subname='(dice_datamode_ssmi_advance): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lsize = size(Si_ifrac)

    if (first_time) then

       ! initialize water(n)
       if (.not. restart_read) then
          ! Note that if restart is read, water is allocated before
          ! the call to this routine in docn_datamode_ssmi_restart_read
          allocate(water(lsize))
          do n = 1,lsize
             if (Si_ifrac(n) > 0.0_cdeps_real_kind) then
                water(n) = flux_Qacc0
             else
                water(n) = 0.0_cdeps_real_kind
             end if
          end do
          ! iFrac0 = iFrac  ! previous step's ice fraction
       endif

       ! reset first time
       first_time = .false.
    end if

    allocate(tfreeze(lsize))
    do n = 1,lsize
       !--- convert to Kelvin
       tfreeze(n) = shr_frz_freezetemp(So_s(n)) + tFrz

       !--- fix erroneous iFrac ---
       Si_ifrac(n) = min(1.0_cdeps_real_kind,max(0.0_cdeps_real_kind,Si_ifrac(n)))

       !--- fabricate ice surface T, fix erroneous iFrac ---
       if ( yc(n) > 0.0_cdeps_real_kind) then
          Si_t(n) = 260.0_cdeps_real_kind + 10.0_cdeps_real_kind*cos(cosArg)
       else
          Si_t(n) = 260.0_cdeps_real_kind - 10.0_cdeps_real_kind*cos(cosArg)
       end if

       !--- set albedos (constant) ---
       Si_avsdr(n) = ax_vsdr
       Si_anidr(n) = ax_nidr
       Si_avsdf(n) = ax_vsdf
       Si_anidf(n) = ax_nidf

       !--- swnet is sent to cpl as a diagnostic quantity only ---
       !--- newly recv'd swdn goes with previously sent albedo ---
       !--- but albedos are (currently) time invariant         ---
       Faii_swnet(n)   = (1.0_cdeps_real_kind - Si_avsdr(n))*Faxa_swvdr(n) &
                       + (1.0_cdeps_real_kind - Si_anidr(n))*Faxa_swndr(n) &
                       + (1.0_cdeps_real_kind - Si_avsdf(n))*Faxa_swvdf(n) &
                       + (1.0_cdeps_real_kind - Si_anidf(n))*Faxa_swndf(n)

       !--- compute melt/freeze water balance, adjust iFrac  -------------
       if ( .not. flux_Qacc ) then ! Q accumulation option is OFF
          Fioi_melth(n) = min(Fioo_q(n),0.0_cdeps_real_kind )          ! q<0 => melt potential
          Fioi_melth(n) = max(Fioi_melth(n),Flux_Qmin   ) ! limit the melt rate
          Fioi_meltw(n) =    -Fioi_melth(n)/latice        ! corresponding water flux

       else                                 ! Q accumulation option is ON
          !--------------------------------------------------------------
          ! 1a) Q<0 & iFrac > 0  =>  infinite supply of water to melt
          ! 1b) Q<0 & iFrac = 0  =>  melt accumulated water only
          ! 2a) Q>0 & iFrac > 0  =>  zero-out accumulated water
          ! 2b) Q>0 & iFrac = 0  =>  accumulated water
          !--------------------------------------------------------------

          if ( Fioo_q(n) <  0.0_cdeps_real_kind ) then ! Q<0 => melt
             if (Si_ifrac(n) > 0.0_cdeps_real_kind ) then
                Fioi_melth(n) = Si_ifrac(n)*max(Fioo_q(n),Flux_Qmin)
                Fioi_meltw(n) =    -Fioi_melth(n)/latice
                !  water(n) = < don't change this value >
             else
                Qmeltall   = -water(n)*latice/dt
                Fioi_melth(n) = max(Fioo_q(n), Qmeltall, Flux_Qmin )
                Fioi_meltw(n) = -Fioi_melth(n)/latice
                water(n) =  water(n) - Fioi_meltw(n)*dt
             end if
          else                       ! Q>0 => freeze
             if (Si_ifrac(n) > 0.0_cdeps_real_kind ) then
                Fioi_melth(n) = 0.0_cdeps_real_kind
                Fioi_meltw(n) = 0.0_cdeps_real_kind
                water(n) = 0.0_cdeps_real_kind
             else
                Fioi_melth(n) = 0.0_cdeps_real_kind
                Fioi_meltw(n) = 0.0_cdeps_real_kind
                water(n) = water(n) + dt*Fioo_q(n)/latice
             end if
          end if

          if (water(n) < 1.0e-16_cdeps_real_kind ) water(n) = 0.0_cdeps_real_kind

          !--- non-zero water => non-zero iFrac ---
          if (Si_ifrac(n) <= 0.0_cdeps_real_kind  .and.  water(n) > 0.0_cdeps_real_kind) then
             Si_ifrac(n) = min(1.0_cdeps_real_kind,water(n)/waterMax)
             ! Si_t(n) = tfreeze(n)     ! T can be above freezing?!?
          end if

          !--- cpl multiplies Fioi_melth & Fioi_meltw by iFrac ---
          !--- divide by iFrac here => fixed quantity flux (not per area) ---
          if (Si_ifrac(n) > 0.0_cdeps_real_kind) then
             Si_ifrac(n) = max( 0.01_cdeps_real_kind, Si_ifrac(n)) ! min iFrac
             Fioi_melth(n) = Fioi_melth(n)/Si_ifrac(n)
             Fioi_meltw(n) = Fioi_meltw(n)/Si_ifrac(n)
          else
             Fioi_melth(n) = 0.0_cdeps_real_kind
             Fioi_meltw(n) = 0.0_cdeps_real_kind
          end if
       end if

       !--- modify T wrt iFrac: (iFrac -> 0) => (T -> tfreeze) ---
       Si_t(n) = tfreeze(n) + Si_ifrac(n)*(Si_t(n)-tfreeze(n))
    end do

    ! compute ice/ice surface fluxes
    call dice_flux_atmice( &
         imask     ,Sa_z      ,Sa_u      ,Sa_v     , &
         Sa_ptem   ,Sa_shum   ,Sa_dens   ,Sa_tbot  , &
         Si_t      ,Faii_sen  ,Faii_lat  ,Faii_lwup, &
         Faii_evap ,Faii_taux ,Faii_tauy ,Si_tref  , &
         Si_qref   ,logunit )

    ! compute ice/oce surface fluxes (except melth & meltw, see above)
    do n=1,lsize
       if (imask(n) == 0) then
          Fioi_swpen(n) = spval
          Fioi_melth(n) = spval
          Fioi_meltw(n) = spval
          Fioi_salt (n) = spval
          Fioi_taux(n)  = spval
          Fioi_tauy(n)  = spval
          Si_ifrac(n)   = 0.0_cdeps_real_kind
       else
          !--- penetrating short wave ---
          Fioi_swpen(n) = max(0.0_cdeps_real_kind, flux_swpf*Faii_swnet(n) ) ! must be non-negative

          !--- i/o surface stress ( = atm/ice stress) ---
          Fioi_taux(n) = Faii_taux(n)
          Fioi_tauy(n) = Faii_tauy(n)

          !--- salt flux ---
          Fioi_salt(n) = 0.0_cdeps_real_kind
       end if
       ! !--- save ifrac for next timestep
       ! iFrac0(n) = Si_ifrac(n)
    end do

    ! Compute outgoing aerosol fluxes
    do n = 1,lsize
       Fioi_bcpho(n) = Faxa_bcph(2,n)
       Fioi_bcphi(n) = Faxa_bcph(1,n) + Faxa_bcph(3,n)
       Fioi_flxdst(n) =  Faxa_dstdry(1,n) + Faxa_dstwet(1,n) + Faxa_dstdry(2,n) + Faxa_dstwet(2,n) &
                       + Faxa_dstdry(3,n) + Faxa_dstwet(3,n) + Faxa_dstdry(4,n) + Faxa_dstwet(4,n)
    end do

    ! optional per thickness category fields
    if (flds_i2o_per_cat) then
       do n = 1,lsize
          Si_iFrac_n(1,n)         = Si_ifrac(n)
          Fioi_swpen_iFrac_n(1,n) = Fioi_swpen(n) * Si_ifrac(n)
       end do
    end if

    ! Reset first_time
    first_time = .false.

    deallocate(tfreeze)

  end subroutine dice_datamode_ssmi_advance

  !===============================================================================
  subroutine dice_datamode_ssmi_restart_write(case_name, inst_suffix, ymd, tod, &
       logunit, my_task, sdat)

    ! input/output variables
    character(len=*)            , intent(in)    :: case_name
    character(len=*)            , intent(in)    :: inst_suffix
    integer                     , intent(in)    :: ymd       ! model date
    integer                     , intent(in)    :: tod       ! model sec into model date
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    type(shr_strdata_type)      , intent(inout) :: sdat
    !-------------------------------------------------------------------------------

    call dshr_restart_write(rpfile, case_name, 'dice', inst_suffix, ymd, tod, &
         logunit, my_task, sdat, fld=water, fldname='water')

  end subroutine dice_datamode_ssmi_restart_write

  !===============================================================================
  subroutine dice_datamode_ssmi_restart_read(rest_filem, inst_suffix, logunit, my_task, mpicom, sdat)

    ! input/output arguments
    character(len=*)            , intent(inout) :: rest_filem
    character(len=*)            , intent(in)    :: inst_suffix
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    integer                     , intent(in)    :: mpicom
    type(shr_strdata_type)      , intent(inout) :: sdat
    !-------------------------------------------------------------------------------

    ! allocate module memory for restart fields that are read in
    allocate(water(sdat%model_lsize))

    ! read restart
    call dshr_restart_read(rest_filem, rpfile, inst_suffix, nullstr, logunit, my_task, mpicom, sdat, &
         fld=water, fldname='water')

  end subroutine dice_datamode_ssmi_restart_read

end module dice_datamode_ssmi_mod
