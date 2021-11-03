! Wrapper for initialization routine
!
! For documentation, see rrtmg_lw_init.f90
!
! Modifications w/r/to rrtmg_lw_init:
! - change real(kind=rb) to real(8)
!   (f2py has trouble with parameterized types)
subroutine rrtmg_lw_ini_wrapper(cpdair)

  use rrtmg_lw_init, only: rrtmg_lw_ini
  implicit none
  real(8), intent(in) :: cpdair
  
  call rrtmg_lw_ini(cpdair)

end subroutine rrtmg_lw_ini_wrapper

! Wrapper for radiative transfer calculation.
!
! For documentation, see rrtmg_lw_rad.f90
!
! Modifications w/r/to rrtmg_lw_rad:
! - change icld from intent(inout) to intent(in)
!   to clarify that modifications will not be reflected
!   in value of input variable after subroutine returns
!   (avoids different behavior for Python scalars and rank-0 arrays)
! - change real(kind=rb) to real(8)
!   and integer(kind=im) to integer(4)
!   (f2py has trouble with parameterized types)
! - convert intent(out) arrays to intent(inout)
!   (f2py has trouble with intent(out) assumed shape arrays)
! - convert optional arguments to required arguments
!   (not sure how f2py handles optional arguments)
subroutine rrtmg_lw_wrapper(                        &
  ncol, nlay, icld, idrv,                           &
  play, plev, tlay, tlev, tsfc,                     &
  h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,     &
  cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis,      &
  inflglw, iceflglw, liqflglw, cldfr,               &
  taucld, cicewp, cliqwp, reice, reliq,             &
  tauaer,                                           &
  start_band, end_band,                             &
  uflx, dflx, hr, uflxc, dflxc, hrc,                &
  wavenumber_range,                                 &
  duflx_dt, duflxc_dt )

  use rrtmg_lw_rad, only: rrtmg_lw
  implicit none

  integer(4), intent(in) :: ncol            
  integer(4), intent(in) :: nlay            
  integer(4), intent(in) :: icld         
  integer(4), intent(in) :: idrv            
  real(8), intent(in) :: play(:,:)          
  real(8), intent(in) :: plev(:,:)          
  real(8), intent(in) :: tlay(:,:)          
  real(8), intent(in) :: tlev(:,:)          
  real(8), intent(in) :: tsfc(:)            
  real(8), intent(in) :: h2ovmr(:,:)        
  real(8), intent(in) :: o3vmr(:,:)         
  real(8), intent(in) :: co2vmr(:,:)        
  real(8), intent(in) :: ch4vmr(:,:)        
  real(8), intent(in) :: n2ovmr(:,:)        
  real(8), intent(in) :: o2vmr(:,:)         
  real(8), intent(in) :: cfc11vmr(:,:)      
  real(8), intent(in) :: cfc12vmr(:,:)      
  real(8), intent(in) :: cfc22vmr(:,:)      
  real(8), intent(in) :: ccl4vmr(:,:)       
  real(8), intent(in) :: emis(:,:)          
  integer(4), intent(in) :: inflglw         
  integer(4), intent(in) :: iceflglw        
  integer(4), intent(in) :: liqflglw        
  real(8), intent(in) :: cldfr(:,:)         
  real(8), intent(in) :: cicewp(:,:)        
  real(8), intent(in) :: cliqwp(:,:)        
  real(8), intent(in) :: reice(:,:)         
  real(8), intent(in) :: reliq(:,:)         
  real(8), intent(in) :: taucld(:,:,:)      
  real(8), intent(in) :: tauaer(:,:,:)      
  integer(4), intent(in) :: start_band
  integer(4), intent(in) :: end_band
  real(8), intent(inout) :: uflx(:,:)         
  real(8), intent(inout) :: dflx(:,:)         
  real(8), intent(inout) :: hr(:,:)           
  real(8), intent(inout) :: uflxc(:,:)        
  real(8), intent(inout) :: dflxc(:,:)        
  real(8), intent(inout) :: hrc(:,:)
  real(8), intent(inout) :: wavenumber_range(:)  
  real(8), intent(inout) :: duflx_dt(:,:)     
  real(8), intent(inout) :: duflxc_dt(:,:)  

  integer icld_inout
 
  icld_inout = icld
  call rrtmg_lw(                                      &  
    ncol, nlay, icld_inout, idrv,                     &
    play, plev, tlay, tlev, tsfc,                     &
    h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,     &
    cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis,      &
    inflglw, iceflglw, liqflglw, cldfr,               &
    taucld, cicewp, cliqwp, reice, reliq,             &
    tauaer,                                           &
    start_band, end_band,                             &
    uflx, dflx, hr, uflxc, dflxc, hrc,                &
    wavenumber_range,                                 &
    duflx_dt, duflxc_dt )

end subroutine rrtmg_lw_wrapper
