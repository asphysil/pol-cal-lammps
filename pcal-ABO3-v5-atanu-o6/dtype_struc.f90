
MODULE data_structure
    use constant, only : dp 
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndim=3, ntype = 3
    INTEGER, PARAMETER :: n1st_max=12, ndim3 = ndim
  
    !INTEGER, PARAMETER :: natms_max=1000000
    !INTEGER, PARAMETER :: na=natms_max/5
  
    INTEGER, PARAMETER :: no_ab_neigh=8, no_ao_neigh=12, no_bo_neigh=6
    REAL(dp) :: dmin_ab = 4.5
    REAL(dp) :: dmin_ao = 3.7
    REAL(dp) :: dmin_bo = 3.5

    REAL(dp) :: unitvol
    !INTEGER, DIMENSION(:) :: atm_info(27*natms1_max)
    !REAL(dp), DIMENSION(:) ::  d(27*natms1_max)
    !REAL(dp), DIMENSION(:,:) :: ab_neigh_coord(3, n1st_max, na)
    !
    REAL(dp),DIMENSION(ndim,ndim) :: latt_vec, ainv
    REAL(dp), DIMENSION(ndim) :: work, ipiv
    INTEGER :: info
    INTEGER :: tndump
  
    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::   a_atms,b_atms,a_disp, b_disp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  a_pol, b_pol, pol_unit
    REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  o_atms
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: sp
    REAL(dp), DIMENSION(:,:,:),ALLOCATABLE ::   pol_ab 
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     ab_neigh
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     ao12_neigh
    INTEGER, DIMENSION(:, :),ALLOCATABLE ::     bo6_neigh
  END MODULE data_structure
