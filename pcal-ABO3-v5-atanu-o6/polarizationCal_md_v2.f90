! Purpose to calculate polarization from a structure

PROGRAM POLARIZATION_CAL
use constant
use data_structure
use nearest_neighbour_list, only : nearest_neighbour
use relative_to_O_disp, only :  cal_a_site_pol_disp,cal_b_site_pol_disp
use writefiles, only : write_xsf, write_pol, write_dump_ovito, write_pol_ovito
use readfiles, only : read_dumpxyz,  readdump_type1, readdump_type2
IMPLICIT NONE   
!     REAL(dp), DIMENSION(:,:) :: frac_pb(3,na),frac_ti(3,na), frac_O(3,3*na)

     REAL(dp),DIMENSION(ndim3) ::  ptemp

     REAL(dp) :: x,y, start, finish, &
                vol,  nti_dp, born_charge
    REAL(dp) :: dab 

     INTEGER :: ntot,  npb, nti, no, ndump, tot_ndump,  no_neigh
     INTEGER :: i, j, k, m, junk1_int

     INTEGER :: fileID, fileID_dump
     CHARACTER(LEN=40):: filename
     LOGICAl :: fexist

npb=0
nti=0
no=0


DO j=1,3
   DO i=1,3
     latt_vec(i,j)=0.0
ENDDO
ENDDO

!PRINT*, " Enter total number of times dumped the atomic coordinates"
!READ*, tot_ndump

tot_ndump = 1

PRINT*,"***************Warning*****************************"
PRINT*,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++"
PRINT*,"%%%%%%%the infomation in INPUT file should be in the following sequence%%%%%%"
PRINT*, " file should be in xyz FORMAT and structure should be ABO3 pervoskite"
PRINT*," At first, the files should have all A atoms coordinates"
PRINT*," After that, the files should have all B atoms coordinates"
PRINT*,"! Last element should be oxige atom Coordinates"
PRINT*,"*********************************************"

PRINT*, " Enter total number of times dumped the atomic coordinates"
READ*, tot_ndump

!
!PRINT*, "Please change total number of times dumped . default is one"
!tot_ndump=1
!!!!!!!!!!!!!!!!!
OPEN(UNIT=12, FILE='dump.xyz', ACTION='READ')
READ(12,*)
DO i=1, 2
  READ(12,*)
ENDDO
READ(12,*) ntot
close(12)
!!!!!!!!!!!!!!!

npb=ntot/5
nti=npb 
no=3*nti
nti_dp =  DBLE(nti)

! ALLOCATE 
ALLOCATE(a_atms(ndim3, npb),b_atms(ndim3, npb))
ALLOCATE( a_disp(ndim3, npb), b_disp(ndim3, npb))
ALLOCATE( a_pol(ndim3, npb), b_pol(ndim3, npb), pol_unit(ndim3, npb)) 

ALLOCATE(sp(ndim3, npb))

ALLOCATE(o_atms(ndim3, no))
ALLOCATE(pol_ab(ndim3,no_ab_neigh,npb))

ALLOCATE(ab_neigh(no_ab_neigh, npb))
ALLOCATE(ao12_neigh(no_ao_neigh, npb))
ALLOCATE(bo6_neigh(no_bo_neigh, npb))

OPEN(UNIT=fileID_in, FILE='dump.xyz', ACTION='READ')

! 
OPEN(UNIT=fileID_totpol,FILE="Total_P.dat", ACTION="WRITE")   
OPEN(UNIT=fileID_locpol,FILE="Local_P.dat", ACTION="WRITE")

OPEN(UNIT=fileID_ovito,FILE="dump-ovito.xyz", ACTION="WRITE")
OPEN(UNIT=fileID_pol_ovito, FILE="dump-pol-ovito.xyz", ACTION="WRITE")

! reading dump file
CALL CPU_TIME(start)
CALL read_dumpxyz(npb, nti, no)

!CALL readdump_type1(npb, nti, no)
!CALL readdump_type2(npb, nti, no)
!
CALL CPU_TIME(finish)
PRINT*, "  //////////////// "
PRINT*, " Total CPU time to read a single dump =", (finish-start)/60.0, "Minutes"
         !  cell volume
vol = ABS( latt_vec(1,1) * (latt_vec(2,2)*latt_vec(3,3) - latt_vec(2,3)*latt_vec(3,2)) &
          -latt_vec(2,1) * (latt_vec(1,2)*latt_vec(3,3) - latt_vec(1,3)*latt_vec(3,2)) &
          +latt_vec(3,1) * (latt_vec(1,2)*latt_vec(2,3) - latt_vec(1,3)*latt_vec(2,2)) )
        
unitvol=vol/nti
!!for calculating inverse Matrix
DO j=1, 3
    DO i=1, 3
     ainv(i,j)=latt_vec(i,j)
   ENDDO
 ENDDO


!! call DGETRF(ndim, ndim, ainv, ndim, ipiv, info)
!! IF (info /= 0) THEN
!!    PRINT*, 'Matrix is numerically singular!'
!!    STOP
!! ENDIF
!!
!!   call DGETRI(ndim, ainv, ndim, ipiv, work, ndim, info)
!!
!!     IF (info /= 0) THEN
!!        PRINT*, 'Matrix inversion failed!'
!!        STOP
!!     ENDIF
!!
!!! Cartesian to fractional coordinate transformation
!     frac_ti=MATMUL(tib(1:nti,:),ainv)
!     frac_pb=MATMUL(pba(1:npb,:),ainv)
!     frac_O=MATMUL(O(1:no,:),ainv)

!!!! finding nearest neighbour distance
CALL CPU_TIME(start)



INQUIRE(FILE="ab-neighbours.dat", exist=fexist)
IF (fexist) THEN
    PRINT*, "ab-neighbours.dat file exist"
    
    OPEN(UNIT=fileID_in2, FILE="ab-neighbours.dat", ACTION="READ")
    do i=1, nti
     READ(fileID_in2, *) junk1_int, (ab_neigh(j,i), j=1,8)
    enddo
ELSE

  PRINT*, "ab-neighbours.dat file does not exist"
  OPEN(UNIT=fileID_in2, FILE="ab-neighbours.dat", ACTION="WRITE")
  dab = dmin_ab
  no_neigh=no_ab_neigh
  CALL  nearest_neighbour(no_neigh, nti, npb, dab, b_atms(:,1:nti), a_atms(:,1:npb), ab_neigh(:,1:nti))
  
   do i=1, nti
   WRITE(fileID_in2, *) i, (ab_neigh(j,i), j=1,8)
   enddo
  
 ENDIF



INQUIRE(FILE="ao12-neighbours.dat", exist=fexist)
IF (fexist) THEN
   PRINT*, "ao12-neighbours.dat file exist"
   
   OPEN(UNIT=fileID_in3, FILE="ao12-neighbours.dat", ACTION="READ")
   
   do i=1, npb
    READ(fileID_in3, *) junk1_int, (ao12_neigh(j,i), j=1,12)
   enddo

ELSE
PRINT*, "ao12-neighbours.dat file does not exist"
OPEN(UNIT=fileID_in3, FILE="ao12-neighbours.dat", ACTION="WRITE")

dab = dmin_ao
no_neigh=no_ao_neigh
CALL  nearest_neighbour(no_neigh, npb, no, dab, a_atms(:,1:npb), o_atms(:,1:no),ao12_neigh(:,1:npb))
 
 do i=1, npb
 WRITE(fileID_in3, *) i, (ao12_neigh(j,i), j=1,12)
 enddo

 ENDIF

INQUIRE(FILE="bo6-neighbours.dat", exist=fexist)
IF (fexist) THEN
   PRINT*, "bo6-neighbours.dat file exist"

   OPEN(UNIT=fileID_in4, FILE="bo6-neighbours.dat", ACTION="READ")

   do i=1, nti
    READ(fileID_in4, *) junk1_int, (bo6_neigh(j,i), j=1,6)
   enddo

ELSE
PRINT*, "bo6-neighbours.dat file does not exist"
OPEN(UNIT=fileID_in4, FILE="bo6-neighbours.dat", ACTION="WRITE")

dab = dmin_bo
no_neigh=no_bo_neigh
CALL  nearest_neighbour(no_neigh, nti, no, dab, b_atms(:,1:nti), o_atms(:,1:no), bo6_neigh(:,1:nti)) 

 do i=1, nti
 WRITE(fileID_in4, *) i, (bo6_neigh(j,i), j=1,6)
 enddo

 ENDIF

!print*, 'ok'
CALL CPU_TIME(finish)
PRINT*, "  //////////////// "
PRINT*, " Total CPU time for calculating nearest neighbour and distance =", (finish-start)/60.0, "Minutes"
!
! Calculate a-site polarization
born_charge=zqa
CALL cal_a_site_pol_disp(born_charge, npb, no_ao_neigh)


! Calculate b-site polarization  
born_charge=zqb
CALL cal_b_site_pol_disp(born_charge, nti,  no_bo_neigh)
 
! arranging data
DO i=1, nti
  DO k=1, no_ab_neigh
    m = ab_neigh(k, i)
   pol_ab(1:3,k,i)= a_pol(1:3,m)
  ENDDO 
ENDDO

  !!Polarization per unit cell
  DO i=1, nti
    ptemp=(/0.0,0.0,0.0/)
    DO k=1, no_ab_neigh
      DO j= 1, 3 
        ptemp(j) = ptemp(j) +  pol_ab(j,k,i)
      ENDDO
    ENDDO ! x y z
    pol_unit(1:3,i) = (b_pol(1:3,i) + (ptemp(1:3)/8.0))
    !print*, (pol_unit(j,i), j=1,3)
ENDDO

! Write file
fileID=30
filename='xcrysden-init.xsf'
CALL write_xsf(fileID, filename,  npb, ntot)
CALL write_pol(npb)
!CALL write_dump_ovito( npb)
!CALL write_pol_ovito(npb) 

PRINT*, " reading dump >1 "

PRINT*, " $$$$$$$$$ "
PRINT*, " Unit of the Local and Total polarization is in C/m^2"
!PRINT*,  (ptot(j), j = 1, 3), SQRT(ptot(1)**2 + ptot(2)**2 + ptot(3)**2)
CALL CPU_TIME(start)
!!!!!!!! Second itrations!!!!!!!!!!!

DO ndump=2, tot_ndump
! reading dump file
CALL read_dumpxyz(npb, nti, no)

!  CALL readdump_type1(npb, nti, no)
!   CALL readdump_type2(npb, nti, no)
   !!*******************************!!
            !  cell volume
    vol = ABS( latt_vec(1,1) * (latt_vec(2,2)*latt_vec(3,3) - latt_vec(2,3)*latt_vec(3,2)) &
               - latt_vec(2,1) * (latt_vec(1,2)*latt_vec(3,3) - latt_vec(1,3)*latt_vec(3,2)) &
               + latt_vec(3,1) * (latt_vec(1,2)*latt_vec(2,3) - latt_vec(1,3)*latt_vec(2,2)) )
    unitvol=vol/nti
   ! Calculate a-site polarization
   born_charge=zqa
   CALL cal_a_site_pol_disp(born_charge, npb, no_ao_neigh)
            
  ! Calculate b-site polarization
   born_charge=zqb
   CALL cal_b_site_pol_disp(born_charge, nti, no_bo_neigh)
    
   ! arranging data 
   DO i=1, nti
       DO k=1, no_ab_neigh
         m = ab_neigh(k, i)
       pol_ab(1:3,k,i)= a_pol(1:3,m)
       ENDDO 
   ENDDO
   
   !!Polarization per unit cell
    DO i=1, nti
        ptemp=(/0.0,0.0,0.0/)
        DO k=1, no_ab_neigh
           DO j= 1, 3 
           ptemp(j) = ptemp(j) +  pol_ab(j,k,i)
          ENDDO
        ENDDO ! x y z
    pol_unit(1:3,i) = (b_pol(1:3,i) + (ptemp(1:3)/8.0))
   ENDDO
  
 !CALL write_dump_ovito( npb) 
 !CALL write_pol_ovito( npb) 
 CALL write_pol(npb)


    IF ( MOD(ndump,100)==0) THEN
           PRINT*, " number of dumps", ndump, "// in step of 100"
    ENDIF
ENDDO ! ndump

fileID=31
filename='xcrysden-final.xsf'
CALL write_xsf(fileID, filename,  npb, ntot)

CALL CPU_TIME(finish)
PRINT*, "  //////////////// "
PRINT*, " Total CPU time =", (finish-start)/3600.0, "hrs"

close(fileID_in)
close(fileID_ovito)
close(fileID_locpol)
close(fileID_totpol)

END PROGRAM POLARIZATION_CAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
