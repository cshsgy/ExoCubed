module ck_opacity_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

  !! Common constants
  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp*pi, fourpi = 4.0_dp*pi
  real(dp), parameter :: half = 1.0_dp/2.0_dp
  real(dp), parameter :: third = 1.0_dp/3.0_dp, twothird = 2.0_dp/3.0_dp

  real(dp), parameter :: amu = 1.66053906660e-24_dp
  real(dp), parameter :: amu_cgs = amu
  real(dp), parameter :: sb  = 5.670374419e-8_dp
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: kb_cgs = kb
  real(dp), parameter :: R_gas = 8.31446261815324_dp
  real(dp), parameter :: hp = 6.62607015e-34_dp
  real(dp), parameter :: c_s = 2.99792458e8_dp
  real(dp), parameter :: kb_si = 1.380649e-23_dp

  type ck_table

    character(len=10) :: sp  ! Species name
    character(len=250) :: path ! Path to table data
    integer :: nwl
    real(dp), allocatable, dimension(:) :: wl ! Central wavelengths of table (um)
    real(dp), allocatable, dimension(:) :: wn ! Central wavenumbers of table (cm-1)

    integer :: iVMR ! VMR index for consituent species

    integer :: nT ! Number of temperature points in table
    real(dp), allocatable, dimension(:) :: T  ! Temperature points in table (K)
    real(dp), allocatable, dimension(:) :: lT  ! Temperature points in table (K)

    integer :: nP ! Number of pressure points in table
    real(dp), allocatable, dimension(:) :: P ! Pressure points in table (bar) [CARE UNITS sometimes!!]
    real(dp), allocatable, dimension(:) :: lP ! Pressure points in table (bar) [CARE UNITS sometimes!!]

    integer :: nG ! Number of g-ordinances (for corr-k)
    real(dp), allocatable, dimension(:) :: Gx, Gw ! Gx ordinances in table, Gw weight of ordinances

    real(dp), allocatable, dimension(:,:,:,:) :: kap ! Kappa values [cm2 molecule-1] (CARE: convert units from table source)

  end type ck_table

  type CIA_table

    character(len=10) :: sp  ! Species name
    logical :: i3 = .False. ! Flag for 3 part special species
    character(len=10), dimension(2) :: sp_con ! Constituent species (for CIA)
    character(len=10), dimension(3) :: sp_con_3 ! Constituent species for 3 part special
    integer, dimension(2) :: iVMR ! VMR index for consituent species
    integer, dimension(3) :: iVMR_3 ! VMR index for 3 part special
    character(len=250) :: path ! Path to table data

    integer :: form

    integer :: nwl
    real(dp), allocatable, dimension(:) :: wl ! Central wavelengths of table (um)
    real(dp), allocatable, dimension(:) :: wn ! Central wavenumbers of table (cm-1)

    integer :: nT ! Number of temperature points in table
    real(dp), allocatable, dimension(:) :: T  ! Temperature points in table (K)

    real(dp), allocatable, dimension(:,:) :: kap ! table values [Usually: cm5 molecule-2] (CARE: Check format units)

  end type CIA_table

  type Ray_table

    character(len=10) :: sp  ! Species name
    integer :: iVMR ! VMR index for consituent species
    character(len=250) :: path ! Path to table data
    integer :: nwl
    real(dp), allocatable, dimension(:) :: kap

  end type Ray_table

  integer :: u_nml

  logical :: first_call = .True.
  integer :: n_wl, n_b, n_g
  real(dp), allocatable, dimension(:) :: wl_e, wl, wn_e, wn, freq
  character(len=250) :: wl_path

  ! CK table data - read in pre-calculated cm2 molecule-1 k-coefficents for each band
  logical :: PM, RORR, AEE
  integer :: ck_interp
  integer :: n_ck
  integer, allocatable, dimension(:) :: ck_form
  type(ck_table), allocatable, dimension(:) :: ck
  character(len=10), allocatable, dimension(:) :: ck_sp
  character(len=250), allocatable, dimension(:) :: ck_paths

    ! CIA table data - read in HITRAN based CIA table then interpolate
  integer :: n_CIA
  integer, allocatable, dimension(:) :: CIA_form
  type(CIA_table), allocatable, dimension(:) :: CIA
  character(len=10), allocatable, dimension(:) :: CIA_sp
  character(len=250), allocatable, dimension(:) :: CIA_paths

  ! Rayleigh table data - read in pre-calcuilated cm2 molecule-1 xsections for each band
  integer :: n_Ray
  integer, allocatable, dimension(:) :: Ray_form
  type(Ray_table), allocatable, dimension(:) :: Ray
  character(len=10), allocatable, dimension(:) :: Ray_sp
  character(len=250), allocatable, dimension(:) :: Ray_paths

  namelist /ck_nml/ PM, RORR, AEE, ck_sp, ck_paths, ck_form, ck_interp
  namelist /CIA_nml/ CIA_sp, CIA_paths, CIA_form
  namelist /Ray_nml/ Ray_sp, Ray_paths, Ray_form


contains

  subroutine ck_opacity(nlay, n_ck, n_CIA, n_Ray, n_b, n_g, wl_in, grav_in, T_in, p_in, pe_in, mu_in, n_sp, sp_list, VMR, &
    & k_tot, a_tot, g_tot)
    implicit none

    ! What's required as input is the layer T, P, VMR
    integer, intent(in) :: nlay, n_ck, n_CIA, n_Ray, n_sp, n_b, n_g
    character(len=10), dimension(n_sp), intent(in) :: sp_list
    real(dp), intent(in) :: grav_in
    real(dp), dimension(nlay), intent(in) :: T_in, p_in, mu_in
    real(dp), dimension(nlay+1), intent(in) :: pe_in
    real(dp), dimension(n_b+1), intent(in) :: wl_in
    real(dp), dimension(n_sp, nlay), intent(in) :: VMR
    real(dp), dimension(n_g,n_b,nlay), intent(out) :: k_tot, a_tot, g_tot

    integer :: b, g, u_nml, k
    real(dp), allocatable, dimension(:,:,:), save :: k_ck_sp
    real(dp), allocatable, dimension(:,:), save :: k_ck
    real(dp), dimension(nlay) :: k_cont, k_Ray

    real(dp) :: grav
    real(dp), dimension(nlay) :: Tl, pl, Nl, RH, mu
    real(dp), dimension(nlay+1) :: pe

    logical, save :: first_call = .True.

    if (first_call .eqv. .True.) then
      call read_gas_data(n_ck, n_CIA, n_Ray, n_sp, n_b, sp_list)
      !! Read input variables from namelist
      allocate(k_ck_sp(n_ck,n_g,nlay))
      allocate(k_ck(n_g,nlay))

      ! Find wavelength centers and wavenumbers, freq etc
      allocate(wl_e(n_b+1), wl(n_b), wn_e(n_b+1), wn(n_b), freq(n_b))
      do b = 1, n_b
        wl_e(b) = wl_in(b)
        wn_e(b) = 1.0_dp/(wl_e(b) * 1.0e-4_dp)
        wl(b) = (wl_in(b+1)+wl_in(b))/2.0_dp
        wn(b) = 1.0_dp/(wl(b) * 1.0e-4_dp)
        freq(b) = c_s / (wl(b) * 1.0e-6_dp)
      end do
      wl_e(n_b+1) = wl_in(n_b+1)
      wn_e(n_b+1) = 1.0_dp/(wl_in(n_b+1) * 1.0e-4_dp)

      first_call = .False.
    end if

    ! Keep everything in CGS units - p in bar
    grav = grav_in*100.0_dp ! gravity in cm s-2
    Tl(:) = T_in(:)
    pe(:) = pe_in(:)*10.0_dp ! pe in dyne
    pl(:) = p_in(:)/1e5_dp
    mu(:) = mu_in(:)
    Nl(:) = (pl(:)*1e6_dp) / (kb_cgs * Tl(:))
    RH(:) = (pl(:)*1e6_dp * mu(:) * amu_cgs) / (kb_cgs * Tl(:))

    ! Loop through each band - we follow various recipes from a mix of CMCRT, CHIMERA and NASA AMES
    do b = 1, n_b

      ! First, bi-linearly interp each species k-tables to mid point T and P
      ! Then perform k-coefficent mixing (Random overlap)
      if (n_ck > 0) then

        !print*, b, 'ck_interp'
        do k = 1, nlay
          if (ck_interp == 1) then
            call interp_ck_sp_TP(n_ck, n_g, Tl(k), pl(k), b, k_ck_sp(:,:,k))
          else if (ck_interp == 2) then
            call interp_ck_sp_TP_Bezier(n_ck, n_g, Tl(k), pl(k), b, k_ck_sp(:,:,k))
          end if
        end do

        if (RORR .eqv. .True.) then
          !print*, b, 'RO'
          do k = 1, nlay
            call random_overlap(n_ck, n_g, n_sp, Nl(k), RH(k), VMR(:,k), k_ck_sp(:,:,k), k_ck(:,k))
          end do
        else if (AEE .eqv. .True.) then
          !print*, b, 'AEE'
          call adap_equiv_extinction(nlay, n_ck, n_g, n_sp, wl(b), grav, &
          &  pe(:), Tl(:), Nl(:), RH(:), VMR(:,:), k_ck_sp(:,:,:), k_ck(:,:))
        else if (PM .eqv. .True.) then
          !print*, b, 'PM'
          k_ck(:,:) = k_ck_sp(1,:,:)
        end if
      else
        k_ck(:,:) = 0.0_dp
      end if

      ! Second, interp to find  continuum opacity
      if (n_CIA > 0) then
        !print*, b, 'conti'
        do k = 1, nlay
          !call interp_conti(n_CIA, n_sp, n_b, b, Tl(k), Nl(k), RH(k), VMR(:,k), k_cont(k))
          call interp_conti_Bezier(n_CIA, n_sp, n_b, b, Tl(k), Nl(k), RH(k), VMR(:,k), k_cont(k))
        end do
      else
        k_cont(:) = 0.0_dp
      end if

      ! Third, calculate Rayleigh scattering
      if (n_Ray > 0) then
        !print*, b, 'Ray'
        do k = 1, nlay
          call calc_Rayleigh(n_Ray, n_sp, b, Nl(k), RH(k), VMR(:,k), k_Ray(k))
        end do
      else
        k_Ray(:) = 0.0_dp
      end if

      ! Fourth, aerosol opacity
      ! (passed to subroutine directly)

      ! Find total values for each g ordinate
      do g = 1, n_g
         k_tot(g,b,:) = k_ck(g,:) + k_cont(:) + k_Ray(:) !+ k_cl_ext(b)  ! Total extinction
         a_tot(g,b,:) = min(k_Ray(:)/k_tot(g,b,:),0.99_dp) !(k_Ray + k_cl_sca(b)) / k_tot(b,g)           ! Effective single scattering albedo
         g_tot(g,b,:) = 0.0_dp !(g_cl(b)*k_cl_sca(b)) / (k_cl_sca(b) + k_Ray) ! Effective asymmetry factor
      end do
      !print*, b, k_ck(1), k_ck(8), k_Ray, k_cont, a_tot(b,1) , a_tot(b,8)
    end do

    ! MKS conversion
    k_tot(:,:,:) = k_tot(:,:,:) * 0.1_dp

    !stop

  end subroutine ck_opacity

  subroutine adap_equiv_extinction(nlay, n_ck, n_g, n_sp, wl_in, grav, pe, Tl, Nl, RH, VMR, k_ck_sp, k_ck)
    implicit none

    integer, intent(in) :: nlay, n_ck, n_g, n_sp
    real(dp), intent(in) :: grav, wl_in
    real(dp), dimension(nlay), intent(in) :: Tl, Nl, RH
    real(dp), dimension(nlay+1), intent(in) :: pe
    real(dp), dimension(n_sp,nlay), intent(in) :: VMR
    real(dp), dimension(n_ck,n_g,nlay), intent(in) :: k_ck_sp

    real(dp), dimension(n_g,nlay), intent(out) :: k_ck

    integer :: k, s, g
    integer, dimension(n_ck) :: max_idx
    real(dp) :: top, bot, dpe, tau_tot, tau_s, tau_all
    real(dp), dimension(nlay) :: bl
    real(dp), dimension(n_ck,nlay) :: k_av
    real(dp), dimension(n_ck) :: tau_av
    logical, dimension(n_ck) :: major

    !! Calculate average kappa in for each species in each layer for thermal component
    !! First just try a general non-weighted average
    !! - seems a good compromise for both stellar and thermal components
    do s = 1, n_ck
      do k = 1, nlay
        !bl(k) = BB(wl_in, Tl(k))
        top = 0.0_dp
        bot = 0.0_dp
        do g = 1, n_g
          top = top + ck(s)%Gw(g)*k_ck_sp(s,g,k) !*bl(k)
          bot = bot + ck(s)%Gw(g) !*bl(k)
        end do
        k_av(s,k) = top/bot * (VMR(ck(s)%iVMR,k) * Nl(k) / RH(k))
      end do
    end do

    !! Find the major absorber species
    !! If we use mean opacity for all species, then using tau is fine
    !! no need to do the transmission functions
    tau_av(:) = 0.0_dp
    tau_tot = 0.0_dp
    do k = 1, nlay
      dpe = pe(k+1) - pe(k)
      do s = 1, n_ck
        tau_s = (k_av(s,k) * dpe)/grav
        tau_av(s) = tau_av(s) + tau_s
      end do
      tau_all = (sum(k_av(:,k)) * dpe)/grav
      tau_tot = tau_tot + tau_all
      !print*, k, pe(k), tau_av(:), tau_tot
      if (tau_tot >= 1.0_dp) then
        exit
      end if
    end do

    !! Find species with highest tau_av
    !! This also checks for the largest tau should the surface be at tau_av < 1
    max_idx = maxloc(tau_av,dim=1)
    !print*, max_idx
    major(:) = .False.
    major(max_idx(1)) = .True.

    !! Combine the grey average opacities of minor species with the k-coefficents of the major species
    k_ck(:,:) = 0.0_dp
    do k = 1, nlay
      do s = 1, n_ck
        if (major(s) .eqv. .True.) then
          ! Combine full k-table
          k_ck(:,k) = k_ck(:,k) + k_ck_sp(s,:,k) * (VMR(ck(s)%iVMR,k) * Nl(k) / RH(k))
        else
          ! Combine average grey value
          k_ck(:,k) = k_ck(:,k) + k_av(s,k)
        end if
      end do
    end do

  end subroutine adap_equiv_extinction

  subroutine random_overlap(n_ck, n_g,  n_sp, N, RH, VMR, k_ck_sp, k_ck)
    implicit none

    integer, intent(in) :: n_ck, n_g, n_sp
    real(dp), intent(in) :: N, RH
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), dimension(n_ck,n_g), intent(in) :: k_ck_sp
    real(dp), dimension(n_g), intent(out) :: k_ck

    integer :: i, j, g, s, n_g2
    real(dp) :: VMR_tot, VMR_cum, intg_sum
    real(dp), dimension(n_g*n_g) :: k_mix, k_mix_sort, logkmix, k_mix_cp
    real(dp), dimension(n_g*n_g) :: wt_mix, wt_mix_sort
    real(dp), dimension(0:n_g*n_g) :: intg, x
    integer, dimension(n_g*n_g) :: sort_indicies
    integer :: loc

    integer :: ix, ix1
    real(dp) :: xval, x0, x1, y0, y1, yval

    !! If 1 ck table, then it's just VMR * ck coefficents
    if (n_ck == 1)then
      k_ck(:) =  VMR(ck(1)%iVMR) * k_ck_sp(1,:) * N / RH
      return
    end if

    !! Convolved size
    n_g2 = n_g * n_g

    !! Start RO procedure
    ! Initial k-table is 1st table * VMR
    VMR_tot = VMR(ck(1)%iVMR)
    k_ck(:) = k_ck_sp(1,:)
    ! Loop over all other tables
    do s = 2, n_ck
      ! Track current cumulative VMR (_cum) and next VMR (_tot)
      VMR_cum = VMR_tot
      VMR_tot = VMR_tot + VMR(ck(s)%iVMR)
      ! Skip this species if very low abundance
      if (VMR(ck(s)%iVMR)*sum(k_ck_sp(s,:)*ck(s)%Gw(:)) < 1.0e-30_dp) then
        cycle
      end if
      ! Loop n_g by n_g, perform random k and weight mixing following Amundsen et al. (2017)
      do i = 1, n_g
        do j = 1, n_g
          k_mix((i-1)*n_g+j) = (VMR_cum*k_ck(i) + VMR(ck(s)%iVMR)*k_ck_sp(s,j)) &
            & / VMR_tot
          wt_mix((i-1)*n_g+j) = ck(s-1)%Gw(i) * ck(s)%Gw(j)
        end do
      end do

      ! Sort the mixed k tables and assosiated weights
      call sort2(n_g2, k_mix, wt_mix)

      ! Now reconstruct the x (g) coordinate
      ! Find cumulative sum of the mixed weights
      intg(0) = 0.0_dp
      intg(1) = wt_mix(1)
      do g = 2, n_g2
        intg(g) = intg(g-1) + wt_mix(g)
      end do
      ! Normalised cumulative sum of weights
      x(:) = intg(:)/maxval(intg) !*2.0_dp - 1.0_dp !(*2 - 1 not needed here, as here weights go 0-1)

      ! Note:, I belive this works in this case, as due to the larger the weight the more
      ! likely the opacity of that x coordinate is to be sampled in a probabilistic sense
      ! (i.e. takes up more range in the x coordinate), so the cumulative weights normalised gives the fraction of the
      ! importance of that g-ordinate to the total opacity distribution

      ! Now interpolate to the origional x grid, this is the mixed k-table.
      do g = 1, n_g
        xval = ck(s)%Gx(g)
        call locate(x(:),xval,ix)
        ix1 = ix + 1
        call linear_log_interp(xval, x(ix), x(ix1), k_mix(ix), k_mix(ix1), k_ck(g))
      end do

    end do

    k_ck(:) = VMR_tot * N * k_ck(:) / RH

  end subroutine random_overlap

  subroutine interp_ck_sp_TP(n_ck,n_g, T, P, b, k_ck_sp)
    implicit none

    integer, intent(in) :: n_ck, n_g, b
    real(dp), intent(in) :: T, P
    real(dp), dimension(n_ck,n_g), intent(out) :: k_ck_sp

    integer :: s, g
    integer :: iT, iT1, iP, iP1
    real(dp) :: lT, lP
    real(dp) :: xval, x0, x1, y0, y1, yval, aval
    real(dp) :: a00, a01, a10, a11

    !! A large if block statement is used to determine if the T-p point
    !! is within the ck table, then interpolates.
    !! If T > Tmax, then T = Tmax, if P > Pmax, P = Pmax and if P < Pmin, P = Pmin
    !! This is done to avoid crashes with k-tables that may not cover the intendend T-p range
    !! IT IS AVDISED TO CREATE k-tables WITHIN APPRORIATE T-p RANGES!!
    !! NOTE: no such check is done for T < Tmin of the table, which will give an error

    lT = log10(T)
    lP = log10(P)

    do s = 1, n_ck

      ! Find temperature grid index
      call locate(ck(s)%T(:),T,iT)
      iT1 = iT + 1
      ! Find pressure grid index
      call locate(ck(s)%P(:),P,iP)
      iP1 = iP + 1

      if (iT == ck(s)%nT) then
        ! Temperature is too high, outside table range - use highest availible T data
        if (iP == ck(s)%nP) then
           ! Pressure is too high, outside table range - use highest availible P data
            k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,iT,iP,:)
        else if (iP == 0) then
           ! Pressure is too low, outside table range - use highest lowest P data
           k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,iT,1,:)
        else
          ! Pressure is within table range, perform linear interpolation at highest T
          do g = 1, ck(s)%nG
            xval = lP
            x0 = ck(s)%lP(iP) ; x1 = ck(s)%lP(iP1)
            y0 = ck(s)%kap(b,iT,iP,g) ; y1 = ck(s)%kap(b,iT,iP1,g)
            call linear_interp(xval, x0, x1, y0, y1, yval)
            !call linear_log_interp(xval, x0, x1, y0, y1, yval)
            k_ck_sp(s,g) = 10.0_dp**yval
          end do
        end if
      else if (iT == 0) then
        ! Temperature is too low, outside table range - use lowest availible T data
        if (iP == ck(s)%nP) then
           ! Pressure is too high, outside table range - use highest availible P data
            k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,1,iP,:)
        else if (iP == 0) then
           ! Pressure is too low, outside table range - use highest lowest P data
           k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,1,1,:)
        else
          ! Pressure is within table range, perform linear interpolation at highest T
          do g = 1, ck(s)%nG
            xval = lP
            x0 = ck(s)%lP(iP) ; x1 = ck(s)%lP(iP1)
            y0 = ck(s)%kap(b,1,iP,g) ; y1 = ck(s)%kap(b,1,iP1,g)
            call linear_interp(xval, x0, x1, y0, y1, yval)
            !call linear_log_interp(xval, x0, x1, y0, y1, yval)
            k_ck_sp(s,g) = 10.0_dp**yval
          end do
        end if
      else
        ! Temperature is within the normal range
        if (iP == ck(s)%nP) then
          ! Pressure is too high, outside table range - use highest availible P
          do g = 1, ck(s)%nG
            xval = lT
            x0 = ck(s)%lT(iT) ; x1 = ck(s)%lT(iT1)
            y0 = ck(s)%kap(b,iT,iP,g) ; y1 = ck(s)%kap(b,iT1,iP,g)
            call linear_interp(xval, x0, x1, y0, y1, yval)
            !call linear_log_interp(xval, x0, x1, y0, y1, yval)
            k_ck_sp(s,g) = 10.0_dp**yval
          end do
        else if (iP == 0) then
          ! Pressure is too low, outside table range - use lowest availible P
          do g = 1, ck(s)%nG
            xval = lT
            x0 = ck(s)%lT(iT) ; x1 = ck(s)%lT(iT1)
            y0 = ck(s)%kap(b,iT,1,g) ; y1 = ck(s)%kap(b,iT1,1,g)
            call linear_interp(xval, x0, x1, y0, y1, yval)
            !call linear_log_interp(xval, x0, x1, y0, y1, yval)
            k_ck_sp(s,g) = 10.0_dp**yval
          end do
        else
          ! Both pressure and temperature are within table bounds, perform bi-linear interpolation
          xval = lT
          yval = lP
          x0 = ck(s)%lT(iT) ; x1 = ck(s)%lT(iT1)
          y0 = ck(s)%lP(iP) ; y1 = ck(s)%lP(iP1)
          do g = 1, ck(s)%nG
            a00 = ck(s)%kap(b,iT,iP,g) ; a01 = ck(s)%kap(b,iT,iP1,g)
            a10 = ck(s)%kap(b,iT1,iP,g) ; a11 = ck(s)%kap(b,iT1,iP1,g)
            call bilinear_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
            !call bilinear_log_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)
            k_ck_sp(s,g) = 10.0_dp**aval
          end do
        end if
      end if

    end do

  end subroutine interp_ck_sp_TP

  subroutine interp_ck_sp_TP_Bezier(n_ck,n_g, T, P, b, k_ck_sp)
    implicit none

    integer, intent(in) :: n_ck, n_g, b
    real(dp), intent(in) :: T, P
    real(dp), dimension(n_ck,n_g), intent(out) :: k_ck_sp

    integer :: s, g
    integer :: iT1, iT2, iT3, iP1, iP2, iP3
    real(dp) :: lT, lP

    real(dp), dimension(3) :: lPa, lTa, lka, lka_ck

    !! A large if block statement is used to determine if the T-p point
    !! is within the ck table, then interpolates.
    !! If T > Tmax, then T = Tmax, if P > Pmax, P = Pmax and if P < Pmin, P = Pmin
    !! This is done to avoid crashes with k-tables that may not cover the intendend T-p range
    !! IT IS AVDISED TO CREATE k-tables WITHIN APPRORIATE T-p RANGES!!

    lT = log10(T)
    lP = log10(P)

    do s = 1, n_ck

      ! Find temperature grid index triplet
      call locate(ck(s)%T(:),T,iT2)
      iT1 = iT2 - 1
      iT3 = iT2 + 1

      if (iT1 <= 0) then
        iT1 = 1
        iT2 = 2
        iT3 = 3
      else if (iT3 > ck(s)%nT) then
        iT1 = ck(s)%nT - 2
        iT2 = ck(s)%nT - 1
        iT3 = ck(s)%nT
      end if

      lTa(1) = ck(s)%lT(iT1)
      lTa(2) = ck(s)%lT(iT2)
      lTa(3) = ck(s)%lT(iT3)

      ! Find pressure grid index triplet
      call locate(ck(s)%P(:),P,iP2)
      iP1 = iP2 - 1
      iP3 = iP2 + 1

      if (iP1 <= 0) then
        iP1 = 1
        iP2 = 2
        iP3 = 3
      else if (iP3 > ck(s)%nP) then
        iP1 = ck(s)%nP - 2
        iP2 = ck(s)%nP - 1
        iP3 = ck(s)%nP
      end if

      lPa(1) = ck(s)%lP(iP1)
      lPa(2) = ck(s)%lP(iP2)
      lPa(3) = ck(s)%lP(iP3)

      if (T >= ck(s)%T(ck(s)%nT)) then
        ! Temperature is too high, outside table range - use highest availible T data
        if (P >= ck(s)%P(ck(s)%nP)) then
           ! Pressure is too high, outside table range - use highest availible P data
            k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,ck(s)%nT,ck(s)%nP,:)
        else if (P <= ck(s)%P(1)) then
           ! Pressure is too low, outside table range - use highest lowest P data
           k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,ck(s)%nT,1,:)
        else
          ! Pressure is within table range, perform Bezier interpolation at highest T
          do g = 1, ck(s)%nG
            lka(1) = ck(s)%kap(b,ck(s)%nT,iP1,g)
            lka(2) = ck(s)%kap(b,ck(s)%nT,iP2,g)
            lka(3) = ck(s)%kap(b,ck(s)%nT,iP3,g)
            call Bezier_interp(lPa(:), lka(:), 3, lP, k_ck_sp(s,g))
            k_ck_sp(s,g) = 10.0_dp**k_ck_sp(s,g)
          end do
        end if
      else if (T <= ck(s)%T(1)) then
        ! Temperature is too low, outside table range - use lowest availible T data
        if (P >= ck(s)%P(ck(s)%nP)) then
           ! Pressure is too high, outside table range - use highest availible P data
            k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,1,ck(s)%nP,:)
        else if (P <= ck(s)%P(1)) then
           ! Pressure is too low, outside table range - use highest lowest P data
           k_ck_sp(s,:) = 10.0_dp**ck(s)%kap(b,1,1,:)
        else
          ! Pressure is within table range, perform linear interpolation at highest T
          do g = 1, ck(s)%nG
            lka(1) = ck(s)%kap(b,1,iP1,g)
            lka(2) = ck(s)%kap(b,1,iP2,g)
            lka(3) = ck(s)%kap(b,1,iP3,g)
            call Bezier_interp(lPa(:), lka(:), 3, lP, k_ck_sp(s,g))
            k_ck_sp(s,g) = 10.0_dp**k_ck_sp(s,g)
          end do
        end if
      else
        ! Temperature is within the normal range
        if (P >= ck(s)%P(ck(s)%nP)) then
          ! Pressure is too high, outside table range - use highest availible P
          do g = 1, ck(s)%nG
            lka(1) = ck(s)%kap(b,iT1,ck(s)%nP,g)
            lka(2) = ck(s)%kap(b,iT2,ck(s)%nP,g)
            lka(3) = ck(s)%kap(b,IT3,ck(s)%nP,g)
            call Bezier_interp(lTa(:), lka(:), 3, lT, k_ck_sp(s,g))
            k_ck_sp(s,g) = 10.0_dp**k_ck_sp(s,g)
          end do
        else if (P <= ck(s)%P(1)) then
          ! Pressure is too low, outside table range - use lowest availible P
          do g = 1, ck(s)%nG
            lka(1) = ck(s)%kap(b,iT1,1,g)
            lka(2) = ck(s)%kap(b,iT2,1,g)
            lka(3) = ck(s)%kap(b,IT3,1,g)
            call Bezier_interp(lTa(:), lka(:), 3, lT, k_ck_sp(s,g))
            k_ck_sp(s,g) = 10.0_dp**k_ck_sp(s,g)
          end do
        else
          ! Both pressure and temperature are within table bounds, perform Bezier interpolation 4 times
          do g = 1, ck(s)%nG
            lka(1) = ck(s)%kap(b,iT1,iP1,g)
            lka(2) = ck(s)%kap(b,iT1,iP2,g)
            lka(3) = ck(s)%kap(b,IT1,iP3,g)
            call Bezier_interp(lPa(:), lka(:), 3, lP, lka_ck(1)) ! Result at T1, P_in
            lka(1) = ck(s)%kap(b,iT2,iP1,g)
            lka(2) = ck(s)%kap(b,iT2,iP2,g)
            lka(3) = ck(s)%kap(b,IT2,iP3,g)
            call Bezier_interp(lPa(:), lka(:), 3, lP, lka_ck(2)) ! Result at T2, P_in
            lka(1) = ck(s)%kap(b,iT3,iP1,g)
            lka(2) = ck(s)%kap(b,iT3,iP2,g)
            lka(3) = ck(s)%kap(b,IT3,iP3,g)
            call Bezier_interp(lPa(:), lka(:), 3, lP, lka_ck(3)) ! Result at T3, P_in
            call Bezier_interp(lTa(:), lka_ck(:), 3, lT, k_ck_sp(s,g)) ! Result at T_in, P_in
            k_ck_sp(s,g) = 10.0_dp**k_ck_sp(s,g)
          end do
        end if
      end if

    end do

  end subroutine interp_ck_sp_TP_Bezier

  subroutine interp_conti(n_CIA, n_sp, n_b, b,  T, N, RH, VMR, k_cont)
    implicit none

    integer, intent(in) :: b, n_CIA, n_sp, n_b
    real(dp), intent(in) :: T, N, RH
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(out) :: k_cont

    integer, allocatable, dimension(:,:), save :: iwn, iwn1
    logical, save :: f_call = .True.

    integer :: s
    integer :: iT, iT1
    real(dp) :: k_Hm, k_ff
    real(dp) :: x0, x1, y0, y1, xval, yval
    real(dp) :: a00, a01, a10, a11, aval

    ! If this is the first call, find the wavenumber index numbers
    ! so we don't have to calculate them each time
    if (f_call .eqv. .True.) then
      if (allocated(iwn) .eqv. .False.) then
        allocate(iwn(n_CIA,n_b), iwn1(n_CIA,n_b))
      end if
      do s = 1, n_CIA
        if ((trim(CIA(s)%sp) == 'H-') .or. (CIA(s)%form == 2)) then
          cycle
        end if
        call locate(CIA(s)%wn(:),wn(b),iwn(s,b))
        iwn1(s,b) = iwn(s,b) + 1
        !print*, wn(b), iwn(s,b), iwn1(s,b), CIA(s)%wn(iwn(s,b)), CIA(s)%wn(iwn1(s,b))
      end do
      if (b == n_b) then
        f_call = .False.
      end if
    end if

    k_cont = 0.0_dp
    do s = 1, n_CIA

      if (trim(CIA(s)%sp) == 'H-') then
        call CIA_Hminus(n_sp, s, b, T, N, VMR, k_Hm)
        k_cont = k_cont + k_Hm
        cycle
      end if
      if (CIA(s)%form == 2) then
        call CIA_ff(n_sp, s, b, T, N, VMR, k_ff)
        k_cont = k_cont + k_ff
        !print*, s, b, k_ff
        cycle
      end if

      ! If band is outside table wavelength range
      if ((iwn1(s,b) > CIA(s)%nwl) .or. (iwn(s,b) < 1)) then
        cycle
      end if

      ! Locate required T indexes in CIA wn array for layer temperature
      call locate(CIA(s)%T(:),T,iT)
      iT1 = iT + 1

      !! Perform temperature edge case check
      if (iT < 1) then
        ! Temperature of layer is outside lower bounds of table
        ! Perform wn linear interp to minval(T)
        xval = wn(b) ; x0 = CIA(s)%wn(iwn(s,b)) ; x1 = CIA(s)%wn(iwn1(s,b))
        y0 = CIA(s)%kap(iwn(s,b),1) ; y1 = CIA(s)%kap(iwn1(s,b),1)

        ! Perform log linear interpolation
        call linear_interp(xval, x0, x1, y0, y1, yval)

        ! Add to result to work variable in units of [cm-1]
        k_cont = k_cont + yval &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N

      else if (iT1 > CIA(s)%nT) then

        ! Temperature of layer is outside upper bounds of table
        ! Perform wn linear interp to maxval(T)
        xval = wn(b) ; x0 = CIA(s)%wn(iwn(s,b)) ; x1 = CIA(s)%wn(iwn1(s,b))
        y0 = CIA(s)%kap(iwn(s,b),CIA(s)%nT) ; y1 = CIA(s)%kap(iwn1(s,b),CIA(s)%nT)

        ! Perform log linear interpolation
        call linear_interp(xval, x0, x1, y0, y1, yval)

        ! Add to result to work variable in units of [cm-1]
        k_cont = k_cont + yval &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N

      else

        !! wn and T are within the table bounds
        xval = wn(b) ; x0 = CIA(s)%wn(iwn(s,b)) ; x1 = CIA(s)%wn(iwn1(s,b))
        yval = T ; y0 = CIA(s)%T(iT) ; y1 = CIA(s)%T(iT1)
        a00 = CIA(s)%kap(iwn(s,b),iT) ; a10 = CIA(s)%kap(iwn1(s,b),iT)
        a01 = CIA(s)%kap(iwn(s,b),iT1) ; a11 = CIA(s)%kap(iwn1(s,b),iT1)

        ! Perform bi-linear interpolation
        call bilinear_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)

        ! Add to result to work variable in units of [cm-1]
        k_cont = k_cont + aval &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N

      end if

    end do

    k_cont = k_cont / RH

  end subroutine interp_conti


  subroutine interp_conti_Bezier(n_CIA, n_sp, n_b, b,  T, N, RH, VMR, k_cont)
    implicit none

    integer, intent(in) :: b, n_CIA, n_sp, n_b
    real(dp), intent(in) :: T, N, RH
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(out) :: k_cont

    integer, allocatable, dimension(:,:), save :: iwn, iwn1
    logical, save :: f_call = .True.

    integer :: s
    integer :: iT, iT1
    real(dp) :: k_Hm, k_ff
    real(dp) :: x0, x1, y0, y1, xval, yval
    real(dp) :: a00, a01, a10, a11, aval

    ! If this is the first call, find the wavenumber index numbers
    ! so we don't have to calculate them each time
    if (f_call .eqv. .True.) then
      if (allocated(iwn) .eqv. .False.) then
        allocate(iwn(n_CIA,n_b), iwn1(n_CIA,n_b))
      end if
      do s = 1, n_CIA
        if ((trim(CIA(s)%sp) == 'H-') .or. (CIA(s)%form == 2)) then
          cycle
        end if
        call locate(CIA(s)%wn(:),wn(b),iwn(s,b))
        iwn1(s,b) = iwn(s,b) + 1
        !print*, wn(b), iwn(s,b), iwn1(s,b), CIA(s)%wn(iwn(s,b)), CIA(s)%wn(iwn1(s,b))
      end do
      if (b == n_b) then
        f_call = .False.
      end if
    end if

    k_cont = 0.0_dp
    do s = 1, n_CIA

      if (trim(CIA(s)%sp) == 'H-') then
        call CIA_Hminus(n_sp, s, b, T, N, VMR, k_Hm)
        k_cont = k_cont + k_Hm
        cycle
      end if
      if (CIA(s)%form == 2) then
        call CIA_ff(n_sp, s, b, T, N, VMR, k_ff)
        k_cont = k_cont + k_ff
        !print*, s, b, k_ff
        cycle
      end if

      ! If band is outside table wavelength range
      if ((iwn1(s,b) > CIA(s)%nwl) .or. (iwn(s,b) < 1)) then
        cycle
      end if

      ! Locate required T indexes in CIA wn array for layer temperature
      call locate(CIA(s)%T(:),T,iT)
      iT1 = iT + 1

      !! Perform temperature edge case check
      if (iT < 1) then
        ! Temperature of layer is outside lower bounds of table
        ! Perform wn linear interp to minval(T)
        xval = wn(b) ; x0 = CIA(s)%wn(iwn(s,b)) ; x1 = CIA(s)%wn(iwn1(s,b))
        y0 = CIA(s)%kap(iwn(s,b),1) ; y1 = CIA(s)%kap(iwn1(s,b),1)

        ! Perform log linear interpolation
        call linear_interp(xval, x0, x1, y0, y1, yval)

        ! Add to result to work variable in units of [cm-1]
        k_cont = k_cont + yval &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N

      else if (iT1 > CIA(s)%nT) then

        ! Temperature of layer is outside upper bounds of table
        ! Perform wn linear interp to maxval(T)
        xval = wn(b) ; x0 = CIA(s)%wn(iwn(s,b)) ; x1 = CIA(s)%wn(iwn1(s,b))
        y0 = CIA(s)%kap(iwn(s,b),CIA(s)%nT) ; y1 = CIA(s)%kap(iwn1(s,b),CIA(s)%nT)

        ! Perform log linear interpolation
        call linear_interp(xval, x0, x1, y0, y1, yval)

        ! Add to result to work variable in units of [cm-1]
        k_cont = k_cont + yval &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N

      else

        !! wn and T are within the table bounds
        xval = wn(b) ; x0 = CIA(s)%wn(iwn(s,b)) ; x1 = CIA(s)%wn(iwn1(s,b))
        yval = T ; y0 = CIA(s)%T(iT) ; y1 = CIA(s)%T(iT1)
        a00 = CIA(s)%kap(iwn(s,b),iT) ; a10 = CIA(s)%kap(iwn1(s,b),iT)
        a01 = CIA(s)%kap(iwn(s,b),iT1) ; a11 = CIA(s)%kap(iwn1(s,b),iT1)

        ! Perform bi-linear interpolation
        call bilinear_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)

        ! Add to result to work variable in units of [cm-1]
        k_cont = k_cont + aval &
          & * VMR(CIA(s)%iVMR(1)) * N &
          & * VMR(CIA(s)%iVMR(2)) * N

      end if

    end do

    k_cont = k_cont / RH

  end subroutine interp_conti_Bezier

  subroutine calc_Rayleigh(n_Ray, n_sp, b, Nl, RH, VMR, k_Ray)
    implicit none

    integer, intent(in) :: b, n_Ray, n_sp
    real(dp), intent(in) :: Nl, RH
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(out) :: k_Ray

    integer :: s

    ! For Rayleigh scattering we just need the pre-tabulated cross section
    ! in cm2 molecule-1 for each Rayleigh scattering species in this band.
    ! The only process that affects the strength of Rayleigh scattering is
    ! changes in the VMR and local total number density

    k_Ray = 0.0_dp
    do s = 1, n_Ray
      k_Ray = k_Ray + VMR(Ray(s)%iVMR) * Ray(s)%kap(b)
      !print*, s, b, k_Ray, VMR(Ray(s)%iVMR), Ray(s)%kap(b)
    end do
    k_Ray = k_Ray * Nl / RH

  end subroutine calc_Rayleigh

  ! Start up module to read in tables and prepare opacities
  subroutine read_gas_data(n_ck, n_CIA, n_Ray, n_sp, n_b,  sp_list)
    implicit none

    integer, intent(in) :: n_ck, n_CIA, n_Ray, n_sp, n_b
    character(len=10), dimension(n_sp), intent(in) :: sp_list

    integer :: u, l, i, j, g, b, n, s, uR
    integer :: n_wlr, ni

    integer :: n_T, n_P, n_G, n_Rays
    real(dp) :: dum1r
    real(dp), allocatable, dimension(:) :: wl_dum

    integer :: stat

    character(len=10) :: name
    integer :: nrec, iwn, iwn1
    real(dp) :: wn_s, wn_e, temp_r, kmax, dum

    logical :: exists
    integer :: iostat_end

    ! Open the namelist file
    open(newunit=u_nml, file='FMS_RC.nml', action='read')


    if (n_ck > 0) then

      print*, 'Reading in ck data'
      !! read in ck data and prep
      allocate(ck_sp(n_ck), ck_paths(n_ck), ck(n_ck), ck_form(n_ck))
      ! Read ck namelist variables
      read(u_nml, nml=ck_nml)
      ck(:)%sp = ck_sp(:)
      ck(:)%path = ck_paths(:)
      ! Find index for VMR
      do s = 1, n_ck
        if (PM .eqv. .True.) then
          cycle
        end if
        exists = .False.
        do j = 1, n_sp
          if (trim(ck(s)%sp) == trim(sp_list(j))) then
            ck(s)%iVMR = j
            exists = .True.
            exit
          end if
        end do
        if (exists .eqv. .False.) then
          print*, 'ERROR - Specifed CK species in namelist not found in species list - STOPPING'
          print*, 'Species: ', ck(s)%sp
          print*, 'List: ', sp_list(:)
          stop
        end if
      end do

      ! Now read in prepared ck_opacity tables
      do s = 1, n_ck
        if (ck_form(s) == 1) then
          print*, 'Reading M. Line ck table: ', trim(ck(s)%sp), ' @ ',  trim(ck(s)%path)
          open(newunit=u, file=trim(ck(s)%path), action='read')
          read(u,*) ck(s)%nwl, ck(s)%nT, ck(s)%nP, ck(s)%nG
          allocate(ck(s)%wl(ck(s)%nwl), ck(s)%wn(ck(s)%nwl))
          allocate(ck(s)%T(ck(s)%nT))
          allocate(ck(s)%P(ck(s)%nP))
          allocate(ck(s)%Gx(ck(s)%nG), ck(s)%Gw(ck(s)%nG))
          allocate(wl_dum(ck(s)%nwl+1))
          do l = 1, ck(s)%nwl+1
            read(u,*) wl_dum(l)
          end do
          do l = 1, ck(s)%nwl
             ck(s)%wl(l) = (wl_dum(l) + wl_dum(l+1)) / 2.0_dp
             ck(s)%wn(l) = 1.0_dp/ (ck(s)%wl(l) * 1e-4_dp)
             !print*, l, ck(s)%wl(l), ck(s)%wn(l)
          end do
          deallocate(wl_dum)
          do g = 1, ck(s)%nG
            read(u,*) ck(s)%Gx(g), ck(s)%Gw(g)
            !print*, g, ck(s)%Gx(g), ck(s)%Gw(g)
          end do
          do i = 1, ck(s)%nT
            read(u,*) ck(s)%T(i)
            !print*, i, ck(s)%T(i)
          end do
          do j = 1, ck(s)%nP
            read(u,*) ck(s)%P(j)
            !print*, j, ck(s)%P(j)
          end do
          read(u,*)
          allocate(ck(s)%kap(ck(s)%nwl,ck(s)%nT,ck(s)%nP,ck(s)%nG))
          do b = 1, ck(s)%nwl
            do i = 1, ck(s)%nT
              do j = 1, ck(s)%nP
                read(u,*) (ck(s)%kap(b,i,j,g), g = 1, ck(s)%nG)
                !print*, b,i,j, (ck(s)%kap(b,i,j,g), g = 1, ck(s)%nG)
              end do
            end do
          end do
        else if (ck_form(s) == 2) then
          print*, 'Reading HELIOS ck table: ', trim(ck(s)%sp), ' @ ',  trim(ck(s)%path)
          open(newunit=u, file=trim(ck(s)%path), action='read')
          read(u,*)
          read(u,*) ck(s)%nT, ck(s)%nP, ck(s)%nwl, ck(s)%nG
          allocate(ck(s)%wl(ck(s)%nwl+1), ck(s)%wn(ck(s)%nwl+1))
          allocate(ck(s)%T(ck(s)%nT),ck(s)%lT(ck(s)%nT))
          allocate(ck(s)%P(ck(s)%nP),ck(s)%lP(ck(s)%nP))
          allocate(ck(s)%Gx(ck(s)%nG), ck(s)%Gw(ck(s)%nG))
          do i = 1, ck(s)%nT
            read(u,*) ck(s)%T(i)
            !print*, i, ck(s)%T(i)
          end do
          ck(s)%lT(:) = log10(ck(s)%T(:))
          do j = 1, ck(s)%nP
            read(u,*) ck(s)%P(j)
            !print*, j, ck(s)%P(j)
          end do
          ck(s)%lP(:) = log10(ck(s)%P(:))
          read(u,*) (ck(s)%wl(g), g = 1, ck(s)%nwl+1)
          !ck(s)%wl(:) =
          ck(s)%wn(:) = 1.0_dp/(ck(s)%wl(:)*1.0e-4_dp)
          read(u,*) (ck(s)%Gx(g), g = 1, ck(s)%nG)
          read(u,*) (ck(s)%Gw(g), g = 1, ck(s)%nG)
          read(u,*)
          allocate(ck(s)%kap(ck(s)%nwl,ck(s)%nT,ck(s)%nP,ck(s)%nG))
          ! Some slight differences between premixed and not premixed opacitiy tables here
          ! Due to GGChem grid ordering in the premixed opacitiy calculation
          if (PM .eqv. .False.) then
            do i = 1, ck(s)%nT
              do j = 1, ck(s)%nP
                do b = ck(s)%nwl, 1, -1
                  read(u,*) (ck(s)%kap(b,i,j,g), g = 1, ck(s)%nG)
                  ck(s)%kap(b,i,j,:) = max(ck(s)%kap(b,i,j,:),1.0e-99_dp)
                  ck(s)%kap(b,i,j,:) = log10(ck(s)%kap(b,i,j,:))
                  !print*, b,i,j, (ck(s)%kap(b,i,j,g), g = 1, ck(s)%nG)
                end do
              end do
            end do
          else if (PM .eqv. .True.) then
            do i = 1, ck(s)%nP
              do j = 1, ck(s)%nT
                do b = ck(s)%nwl, 1, -1
                  read(u,*) (ck(s)%kap(b,j,i,g), g = 1, ck(s)%nG)
                  ck(s)%kap(b,j,i,:) = max(ck(s)%kap(b,j,i,:),1.0e-99_dp)
                  ck(s)%kap(b,j,i,:) = log10(ck(s)%kap(b,j,i,:))
                  !print*, b,i,j, (ck(s)%kap(b,j,i,g), g = 1, ck(s)%nG)
                end do
              end do
            end do
          end if
        else
          print*, 'No ck species, n_ck = ', n_ck
          stop
        end if
      end do
    end if

    if (n_CIA > 0) then

      !! read in ck data and prep
      allocate(CIA_sp(n_CIA), CIA_paths(n_CIA), CIA(n_CIA), CIA_form(n_CIA))

      ! Read the cia namelist
      read(u_nml, nml=cia_nml)

      CIA(:)%sp = CIA_sp(:)
      CIA(:)%path = CIA_paths(:)
      CIA(:)%form = CIA_form(:)

      ! Find the CIA constituents from lookup table
      call find_CIA_consituents(n_CIA)

      do s = 1, n_CIA

        ! Check for 3 species special
        if (CIA(s)%i3 .eqv. .True.) then
          ni = 3
        else
          ni = 2
        end if

        do i = 1, ni
          exists = .False.
          do j = 1, n_sp
            if (ni == 2) then
              if (CIA(s)%sp_con(i) == sp_list(j)) then
                CIA(s)%iVMR(i) = j
                exists = .True.
                exit
              end if
            else if (ni == 3) then
              if (CIA(s)%sp_con_3(i) == sp_list(j)) then
                CIA(s)%iVMR_3(i) = j
                exists = .True.
                exit
              end if
            end if
          end do

          if (exists .eqv. .False.) then
            print*, 'ERROR - Specified CIA species component not found in prf VMR list - STOPPING'
            if (ni == 2) then
              print*, 'Species 2 part: ', CIA(s)%sp, CIA(s)%sp_con(i)
            else if (ni == 3) then
              print*, 'Species 3 part: ', CIA(s)%sp, CIA(s)%sp_con_3(i)
            end if
            stop
          end if

        end do
      end do
    end if

    ! Read in prepared CIA tables
    if (n_CIA > 0) then

      do s = 1, n_CIA

        print*, 'Reading CIA species: ', trim(CIA(s)%sp), ' @ ', trim(CIA(s)%path)

        ! Skip H-
        if ((trim(CIA(s)%sp) == 'H-')) then
          cycle
        end if

        if (CIA(s)%form == 1) then

          open(newunit=u,file=trim(CIA(s)%path),action='read',status='old')

          ! Allocate CIA table temperature arrays
          allocate(CIA(s)%T(CIA(s)%nT))

          ! Read and allocate data until error (end of file)
          do n = 1, CIA(s)%nT
            read(u,*,iostat=stat) name, wn_s, wn_e, nrec, temp_r, kmax, dum
            !print*, n, name, wn_s, wn_e, nrec, temp_r, kmax, dum
            if (n == 1) then
              ! Allocate CIA table wn and table value array
              allocate(CIA(s)%wn(nrec))
              allocate(CIA(s)%kap(nrec,CIA(s)%nT))
              CIA(s)%nwl = nrec
            end if
            ! Check if end of file reached
            if (is_iostat_end(stat)) then
              print*,'Reached end of HITRAN CIA file: ', CIA(s)%sp, CIA(s)%path
              exit
            else
              ! Temperature point of table
              CIA(s)%T(n) = temp_r
              ! Read the record data
              do i = 1, nrec
                read(u,*) CIA(s)%wn(i), CIA(s)%kap(i,n)
                !print*, i, CIA(s)%wn(i), CIA(s)%kap(i,n)
              end do
            end if

          end do
          close(u)

        else if (CIA(s)%form == 2) then

          CIA(s)%nT = 8

          if (CIA(s)%sp == 'He-') then
            CIA(s)%nwl = 16
          else if (CIA(s)%sp == 'H2-') then
            CIA(s)%nwl = 18
          end if

          open(newunit=u,file=trim(CIA(s)%path),action='read',status='old')
          read(u,*)

          allocate(CIA(s)%wl(CIA(s)%nwl))
          allocate(CIA(s)%wn(CIA(s)%nwl))
          allocate(CIA(s)%T(CIA(s)%nT))
          allocate(CIA(s)%kap(CIA(s)%nwl,CIA(s)%nT))

          read(u,*) CIA(s)%T(:)

          do n = 1, CIA(s)%nwl
            read(u,*) CIA(s)%wl(n), CIA(s)%kap(n,:)
            !print*, s, n, CIA(s)%wl(n), CIA(s)%kap(n,:)
          end do
          !CIA(s)%wl(:) = CIA(s)%wl(:) * 1.0e-4_dp
          CIA(s)%wn(:) = 1.0_dp/(CIA(s)%wl(:) * 1.0e-8_dp)


          close(u)
        end if
      end do
    end if

    if (n_Ray > 0) then
      ! Read Rayleigh namelist
      allocate(Ray(n_Ray))
      allocate(Ray_form(n_Ray),Ray_sp(n_Ray),Ray_paths(n_Ray))
      read(u_nml, nml=Ray_nml)

      Ray(:)%sp = Ray_sp(:)
      Ray(:)%path = Ray_paths(:)

      do s = 1, n_Ray
        exists = .False.
        do j = 1, n_sp
          if (Ray(s)%sp == sp_list(j)) then
            Ray(s)%iVMR = j
            exists = .True.
            exit
          end if
        end do

        if (exists .eqv. .False.) then
          print*, 'ERROR - Specified Rayleigh species component not found in prf VMR list - STOPPING'
          print*, 'Species: ', s, Ray(s)%sp
          stop
        end if
      end do

      ! Read in Rayleigh scattering data
      if (n_Ray > 0) then
        do s = 1, n_Ray
          allocate(Ray(s)%kap(n_b))
          print*, 'Reading Rayleigh data: ', trim(Ray(s)%sp), ' @ ',trim(Ray(s)%path)
          open(newunit=uR,file=trim(Ray(s)%path),action='read',status='old')
          read(uR,*)
          read(uR,*) Ray(s)%kap(:)
          !print*, s, Ray(s)%kap(:)
          close(uR)
        end do
      end if
    end if


    print*, 'Finished reading'

  end subroutine read_gas_data

  ! Adapted CMCRT H- calculation
  subroutine CIA_Hminus(n_sp, s, b, T, Nl, VMR, k_cia)
    implicit none
    ! John (1988) paramaters
    real(dp), dimension(6), parameter :: &
      & An_ff1 = (/518.1021_dp, 472.2636_dp, -482.2089_dp, 115.5291_dp, 0.0_dp, 0.0_dp/), &
      & Bn_ff1 = (/-734.8666_dp, 1443.4137_dp, -737.1616_dp, 169.6374_dp, 0.0_dp, 0.0_dp/), &
      & Cn_ff1 = (/1021.1775_dp, -1977.3395_dp, 1096.8827_dp, -245.6490_dp, 0.0_dp, 0.0_dp/), &
      & Dn_ff1 = (/-479.0721_dp, 922.3575_dp, -521.1341_dp, 114.2430_dp, 0.0_dp, 0.0_dp/), &
      & En_ff1 = (/93.1373_dp, -178.9275_dp, 101.7963_dp, -21.9972_dp, 0.0_dp, 0.0_dp/), &
      & Fn_ff1 = (/-6.4285_dp, 12.3600_dp, -7.0571_dp, 1.5097_dp, 0.0_dp, 0.0_dp/)

    real(dp), dimension(6), parameter :: &
      & An_ff2 = (/0.0_dp, 2483.3460_dp, -3449.8890_dp, 2200.0400_dp, -696.2710_dp, 88.2830_dp/), &
      & Bn_ff2 = (/0.0_dp, 285.8270_dp, -1158.3820_dp, 2427.7190_dp, -1841.4000_dp, 444.5170_dp/), &
      & Cn_ff2 = (/0.0_dp, -2054.2910_dp, 8746.5230_dp, -13651.1050_dp, 8642.9700_dp, -1863.8640_dp/), &
      & Dn_ff2 = (/0.0_dp, 2827.7760_dp, -11485.6320_dp, 16755.5240_dp, -10051.5300_dp, 2095.2880_dp/), &
      & En_ff2 = (/0.0_dp, -1341.5370_dp, 5303.6090_dp, -7510.4940_dp, 4400.0670_dp, -901.7880_dp/), &
      & Fn_ff2 = (/0.0_dp, 208.9520_dp, -812.9390_dp, 1132.7380_dp, -655.0200_dp, 132.9850_dp/)

    real(dp), dimension(6), parameter :: &
      & Cn_bf = (/152.519, 49.534, -118.858, 92.536, -34.194, 4.982 /)

    real(dp), parameter :: alf = 1.439e8_dp, lam_0 = 1.6419_dp

    integer, intent(in) :: s, b, n_sp
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(in) :: T, Nl
    real(dp), intent(out) :: k_cia

    integer :: n
    real(dp) :: kff, kbf, fbf, xbf, sff
    real(dp) :: T5040

    T5040 = 5040.0_dp / T

    ! Do bound free calculation
    if (wl(b) > lam_0) then
      xbf = 0.0_dp
      !kbf = 0.0_dp
    else
      fbf = 0.0_dp
      do n = 1, 6
        fbf = fbf + Cn_bf(n) * (1.0_dp/wl(b) - 1.0_dp/lam_0)**((real(n,dp)-1.0_dp)/2.0_dp)
      end do
      xbf = 1.0e-18_dp * wl(b)**3 * (1.0_dp/wl(b) - 1.0_dp/lam_0)**(3.0_dp/2.0_dp) * fbf
      !kbf = 0.750_dp * T**(-5.0_dp/2.0_dp) * exp(alf/(lam_0*T)) * (1.0_dp - exp(-alf/(wl(b)*T))) * xbf
    end if

    ! Do free free calculation
    sff = 0.0_dp
    if (wl(b) >= 0.3645_dp) then
      do n = 1, 6
        sff = sff + T5040**((real(n,dp)+1.0_dp)/2.0_dp) &
                & * (wl(b)**2*An_ff2(n) + Bn_ff2(n)  + Cn_ff2(n)/wl(b) + Dn_ff2(n)/wl(b)**2 + En_ff2(n)/wl(b)**3 &
                & + Fn_ff2(n)/wl(b)**4)
      end do
      kff = 1.0e-29_dp * sff
    else if ((wl(b) < 0.3645_dp) .and. (wl(b) > 0.1823_dp)) then
      do n = 1, 6
        sff = sff + T5040**((real(n,dp)+1.0_dp)/2.0_dp) &
                & * (wl(b)**2*An_ff1(n) + Bn_ff1(n)  + Cn_ff1(n)/wl(b) + Dn_ff1(n)/wl(b)**2 + En_ff1(n)/wl(b)**3 &
                & + Fn_ff1(n)/wl(b)**4)
      end do
      kff = 1.0e-29_dp * sff
    else
      kff = 0.0_dp
    end if

    ! xbf is in [cm2 molecule-1] and kff is in [cm4 dyne-1] - convert to [cm-1] before CMCRT output
    kbf = xbf * VMR(CIA(s)%iVMR_3(1)) * Nl !! * H- [molecule cm-3]
    !kbf = kbf * (VMR_lay(CIA_tab(s)%iVMR_3(2),z) * N_lay(z) * VMR_lay(CIA_tab(s)%iVMR_3(3),z) * N_lay(z)) * kb * T ! * P(e-) [dyne cm-2] * H [cm-3]
    kff = kff * (VMR(CIA(s)%iVMR_3(2)) * Nl * VMR(CIA(s)%iVMR_3(3)) * Nl) * kb_cgs * T !! * P(e-) [dyne cm-2] * H [cm-3]

    k_cia = kbf + kff

  end subroutine CIA_Hminus

  subroutine CIA_Heminus(n_sp, s, b, T, Nl, VMR, k_Hem)
    implicit none

    integer, intent(in) :: s, b, n_sp
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(in) :: T, Nl
    real(dp), intent(out) :: k_Hem

    real(dp) :: aa, bb, cc, kff

    ! Polynomial coefficents, frequency f in [Hz]
    aa = 3.397e-46_dp + (-5.216e-31_dp + 7.039e-15_dp/freq(b))/freq(b)
    bb = -4.116e-42_dp + (1.067e-26_dp + 8.135e-11_dp/freq(b))/freq(b)
    cc = 5.081e-37_dp + (-8.724e-23_dp - 5.659e-8_dp/freq(b))/freq(b)
    kff = aa*T + bb + cc/T

    kff = kff * VMR(CIA(s)%iVMR(1)) * Nl * VMR(CIA(s)%iVMR(2)) * Nl !  He [cm-3] * e- [cm-3]

    ! End units [cm-1]

    k_Hem = kff

  end subroutine CIA_Heminus

  subroutine CIA_ff(n_sp, s, b, T, Nl, VMR, kff)
    implicit none

    integer, intent(in) :: s, b, n_sp
    real(dp), dimension(n_sp), intent(in) :: VMR
    real(dp), intent(in) :: T, Nl
    real(dp), intent(out) :: kff

    integer :: iwl, iwl1, iT, iT1
    real(dp) :: T5040, wlA
    real(dp) :: xval, x0, x1, yval, y0, y1
    real(dp) :: a00, a10, a01, a11, aval

    T5040 = 5040.0_dp / T
    wlA = wl(b) * 1.0e4_dp

    call locate(CIA(s)%wl(:),wlA,iwl)
    iwl1 = iwl + 1

    !print*, s, b, wlA, CIA(s)%wl(iwl), CIA(s)%wl(iwl1)

    if (iwl == 0 .or. iwl1 > CIA(s)%nwl) then
       kff = 0.0_dp
       return
    end if

    call locate(CIA(s)%T(:),T5040,iT)
    iT1 = iT + 1

    !print*, s, b, T5040, CIA(s)%T(iT), CIA(s)%T(iT1)

    if (iT == 0 .or. iT1 > CIA(s)%nT) then
       kff = 0.0_dp
       return
    end if

    !! wn and T are within the table bounds
    xval = wlA ; x0 = CIA(s)%wl(iwl) ; x1 = CIA(s)%wl(iwl1)
    yval = T5040 ; y0 = CIA(s)%T(iT) ; y1 = CIA(s)%T(iT1)
    a00 = CIA(s)%kap(iwl,iT) ; a10 = CIA(s)%kap(iwl1,iT)
    a01 = CIA(s)%kap(iwl,iT1) ; a11 = CIA(s)%kap(iwl1,iT1)

    ! Perform bi-linear interpolation
    call bilinear_interp(xval, yval, x0, x1, y0, y1, a00, a10, a01, a11, aval)

   ! print*, s,b,wl(b), T, aval

    kff = 1.0e-26_dp * aval * (VMR(CIA(s)%iVMR(1)) * Nl * VMR(CIA(s)%iVMR(2)) * Nl) * kb_cgs * T !! * P(e-) [dyne cm-2] * (H2 or He) [cm-3]

  end subroutine CIA_ff


  subroutine find_CIA_consituents(n_CIA)
    implicit none

    integer, intent(in) :: n_CIA
    integer :: s

    do s = 1, n_CIA

      select case(CIA(s)%sp)

      case('H2-H2')
        CIA(s)%sp_con(1) = 'H2'
        CIA(s)%sp_con(2) = 'H2'

        CIA(s)%nT = 113

      case('H2-He','He-H2')
        CIA(s)%sp_con(1) = 'H2'
        CIA(s)%sp_con(2) = 'He'

        CIA(s)%nT = 334

      case('H2-H','H-H2')
        CIA(s)%sp_con(1) = 'H2'
        CIA(s)%sp_con(2) = 'H'

        CIA(s)%nT = 4

      case('H-He','He-H')
        CIA(s)%sp_con(1) = 'He'
        CIA(s)%sp_con(2) = 'H'

        CIA(s)%nT = 10

      case('H-')
        CIA(s)%sp_con_3(1) = 'H-'
        CIA(s)%sp_con_3(2) = 'H'
        CIA(s)%sp_con_3(3) = 'e-'

        CIA(s)%i3 = .True.

      case('He-')
        CIA(s)%sp_con(1) = 'He'
        CIA(s)%sp_con(2) = 'e-'

      case('H2-')
        CIA(s)%sp_con(1) = 'H2'
        CIA(s)%sp_con(2) = 'e-'

      case default
        print*, 'ERROR - CIA species constituents could not be found - STOPPING'
        print*, 'Species: ', CIA(s)%sp
        stop
      end select

    end do

  end subroutine find_CIA_consituents

  subroutine locate(arr, var, idx)
    implicit none

    integer, intent(out) :: idx
    real(dp), dimension(:), intent(in) :: arr
    real(dp),intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section/binary search (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = size(arr)+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if (var > arr(jm)) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

  ! Perform linear interpolation  space
  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = xval
    lx1 = x1; lx2 = x2
    ly1 = y1; ly2 = y2

    norm = 1.0_dp / (lx2 - lx1)

    yval = ((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_interp

  subroutine bilinear_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(dp) :: lxval, lyval, lx1, lx2, ly1, ly2, la11, la21, la12, la22
    real(dp), intent(out) :: aval
    real(dp) :: norm

    lxval = xval ; lyval = yval
    lx1 = x1 ; lx2 = x2 ; ly1 = y1 ; ly2 = y2
    la11 = a11 ; la21 = a21 ; la12 = a12 ; la22 = a22

    norm = 1.0_dp / (lx2 - lx1) / (ly2 - ly1)

    aval = la11 * (lx2 - lxval) * (ly2 - lyval) * norm &
      & + la21 * (lxval - lx1) * (ly2 - lyval) * norm &
      & + la12 * (lx2 - lxval) * (lyval - ly1) * norm &
      & + la22 * (lxval - lx1) * (lyval - ly1) * norm

    aval = aval

  end subroutine bilinear_interp

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2
    real(dp) :: lxval, ly1, ly2, lx1, lx2
    real(dp), intent(out) :: yval
    real(dp) :: norm

    lxval = log10(xval)
    lx1 = log10(x1); lx2 = log10(x2)
    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / (lx2 - lx1)

    yval = 10.0_dp**((ly1 * (lx2 - lxval) + ly2 * (lxval - lx1)) * norm)

  end subroutine linear_log_interp

  subroutine bilinear_log_interp(xval, yval, x1, x2, y1, y2, a11, a21, a12, a22, aval)
    implicit none

    real(dp), intent(in) :: xval, yval, x1, x2, y1, y2, a11, a21, a12, a22
    real(dp) :: lxval, lyval, lx1, lx2, ly1, ly2, la11, la21, la12, la22
    real(dp), intent(out) :: aval
    real(dp) :: norm

    lxval = log10(xval) ; lyval = log10(yval)
    lx1 = log10(x1) ; lx2 = log10(x2) ; ly1 = log10(y1) ; ly2 = log10(y2)
    la11 = log10(a11) ; la21 = log10(a21) ; la12 = log10(a12) ; la22 = log10(a22)

    norm = 1.0_dp / (lx2 - lx1) / (ly2 - ly1)

    aval = la11 * (lx2 - lxval) * (ly2 - lyval) * norm &
      & + la21 * (lxval - lx1) * (ly2 - lyval) * norm &
      & + la12 * (lx2 - lxval) * (lyval - ly1) * norm &
      & + la22 * (lxval - lx1) * (lyval - ly1) * norm

    aval = 10.0_dp**(aval)

  end subroutine bilinear_log_interp

  ! Perform Bezier interpolation
  subroutine Bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: xc, dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x > xi(1) .and. x < xi(2)) then
      ! left hand side interpolation
      !print*,'left'
      w = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine Bezier_interp

  subroutine sort2(N,RA,RB)
    integer, intent(in) :: N
    integer :: L, IR, I, J
    real(dp), dimension(N), intent(inout) :: RA, RB
    real(dp) :: RRA, RRB
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
      L=L-1
      RRA=RA(L)
      RRB=RB(L)
    ELSE
      RRA=RA(IR)
      RRB=RB(IR)
      RA(IR)=RA(1)
      RB(IR)=RB(1)
      IR=IR-1
      IF(IR.EQ.1)THEN
        RA(1)=RRA
        RB(1)=RRB
        RETURN
      ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
      IF(J.LT.IR)THEN
        IF(RA(J).LT.RA(J+1)) THEN
          J=J+1
        ENDIF
      ENDIF
      IF(RRA.LT.RA(J))THEN
        RA(I)=RA(J)
        RB(I)=RB(J)
        I=J
        J=J+J
      ELSE
        J=IR+1
      ENDIF
      GO TO 20
    ENDIF
    RA(I)=RRA
    RB(I)=RRB
    GO TO 10
  end subroutine sort2

  real(dp) function BB(wl_in, T_in)
    implicit none

    ! Planck function in wavelength units

    real(dp), intent(in) :: wl_in, T_in
    real(dp) :: left, right
    real(dp) :: wl_cm

    wl_cm = wl_in * 1.0e-4_dp

    left = (2.0_dp * hp * c_s**2)/(wl_cm)**5
    right = 1.0_dp / (exp((hp * c_s) / (wl_cm * kb * T_in)) - 1.0_dp)
    BB = left * right

  end function BB

end module ck_opacity_mod
