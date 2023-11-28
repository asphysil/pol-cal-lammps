MODULE writefiles
    use constant, only : dp,& 
             fileID_locpol, fileID_totpol, fileID_ovito,&
             fileID_pol_ovito
    use data_structure, only : latt_vec, a_atms, b_atms, o_atms,&
                              a_pol, b_pol, pol_unit, sp, tndump 
    IMPLICIT NONE 
 
 private    
 public ::  write_xsf, write_pol ,write_dump_ovito, write_pol_ovito 
contains 
SUBROUTINE write_xsf(fileID, filename,  nb, ntot)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: fileID,  nb, ntot 
    CHARACTER(LEN=40), INTENT(IN) :: filename
    INTEGER :: no 
    REAL :: f(3)  
    INTEGER :: i, j 
    INTEGER :: atmnum(3)

    no = 3*nb 
    f = (/0.0, 0.0, 0.0/)
    atmnum(1) = 83 ! Bi
    atmnum(2) = 26 ! Fe
    atmnum(3) = 8  ! O
OPEN(UNIT=fileID, FILE=filename, ACTION="WRITE")
!!!!
!!! Polarization calculation
!!!!!!WRITEING for xcruden .xsf file FORMAT
WRITE(fileID,*) '#A B O'
WRITE(fileID,*)'CRYSTAL'
WRITE(fileID,*)'PRIMVEC'
DO i=1,3
  WRITE(fileID,300) latt_vec(i,1), latt_vec(i,2), latt_vec(i,3)
ENDDO
WRITE(fileID,*)'CONVVEC'
DO i=1,3
  WRITE(fileID,300)latt_vec(i,1), latt_vec(i,2), latt_vec(i,3)
ENDDO
!300 FORMAT(1X, 3F13.8)
WRITE(fileID,*)'PRIMCOORD'
WRITE(fileID,*) ntot, '1'

DO i =1, nb 
WRITE(fileID,301) atmnum(1), (a_atms(j,i),j=1,3), (a_pol(j,i), j=1,3)
ENDDO

DO i =1, nb 
WRITE(fileID,301) atmnum(2), (b_atms(j,i),j=1,3), (b_pol(j,i), j=1,3)
ENDDO

DO i=1, no
  WRITE(fileID,301) atmnum(3), (o_atms(j,i), j=1,3), f(1), f(2), f(3)
ENDDO
300 FORMAT(2X, 3F12.7)
301 FORMAT (2X, I5, 6F12.8)

close(fileID)
END SUBROUTINE write_xsf 

SUBROUTINE write_pol(nb)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) ::  nb 
    REAL(dp) :: ptot(3)  
    REAL(dp) :: nb_dp
    REAL(dp) :: px, py, pz  
    INTEGER :: i, j 



nb_dp = DBLE(nb)

ptot(1)=SUM(pol_unit(1,1:nb))/nb_dp
ptot(2)=SUM(pol_unit(2,1:nb))/nb_dp
ptot(3)=SUM(pol_unit(3,1:nb))/nb_dp

! PRINT*, " ----------------------Warning-------------"
! PRINT*, " Please change total polarization calculation formula --" 
! PRINT*, 'if the following strutural infomations are  not correct'

! PRINT*, " x, y, and z axis are  along a, b and c crystalographic direction, respectively. Axis are orthogonal"
! !PRINT*, " x-axis along [1-10], y-axis along [110] and z-axis along [001]"
! PRINT*, "-----------------------------------------------"
! !px = ptot(1)/sqrt(2.0) + ptot(2)/sqrt(2.0)
! !py =  - ptot(1)/sqrt(2.0) + ptot(2)/sqrt(2.0)
! !pz = ptot(3)
! !
! PRINT*, " x-axis along [100], y-axis along [100] and z-axis along [001]"
!
px = ptot(1)
py = ptot(2)
pz = ptot(3)
WRITE(fileID_totpol, 303) px, py, pz, SQRT(px**2 + py**2 + pz**2)

DO i=1, nb 
    WRITE(fileID_locpol,302) (b_atms(j,i),j=1,3), (pol_unit(j,i),j=1,3)
ENDDO 
302 FORMAT(2X, 3F12.7, 3X, 3F12.7)
303 FORMAT(2X,  3F12.7, 3X, F12.7)

END SUBROUTINE write_pol

SUBROUTINE write_dump_ovito(nb)
INTEGER, INTENT(IN) ::  nb
INTEGER :: i, j
INTEGER :: id=1 
REAL :: x 
REAL(dp) :: ploc 
x=0.0000

write(fileID_ovito,'(A)') 'ITEM: TIMESTEP'
write(fileID_ovito,'(I7)') tndump
write(fileID_ovito,'(A)') 'ITEM: NUMBER OF ATOMS'
write(fileID_ovito, '(I7)') nb 
write(fileID_ovito, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
DO i=1,3
    write(fileID_ovito,'(2F15.8)') x, latt_vec(i,i)
ENDDO
write(fileID_ovito, '(A)') 'ITEM: ATOMS x y z fx fy fz mux muy muz mu' 
DO i = 1, nb
ploc = SQRT(pol_unit(1,i)*pol_unit(1,i) + pol_unit(2,i)*pol_unit(2,i) + pol_unit(3,i)*pol_unit(3,i) )
write(fileID_ovito, 400) (b_atms(j,i), j=1,3), (sp(j,i), j=1,3), (pol_unit(j,i),j=1,3), ploc
ENDDO
400 FORMAT(1X, 3F15.10, 3F15.10, 3F15.10, F15.10)
END SUBROUTINE write_dump_ovito

SUBROUTINE write_pol_ovito(nb)
INTEGER, INTENT(IN) ::  nb
INTEGER :: i, j
INTEGER :: id=1
REAL :: x
REAL(dp) :: ploc
x=0.0000

write(fileID_pol_ovito,'(A)') 'ITEM: TIMESTEP'
write(fileID_pol_ovito,'(I7)') tndump
write(fileID_pol_ovito,'(A)') 'ITEM: NUMBER OF ATOMS'
write(fileID_pol_ovito, '(I7)') nb
write(fileID_pol_ovito, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
DO i=1,3
    write(fileID_pol_ovito,'(2F15.8)') x, latt_vec(i,i)
ENDDO
write(fileID_pol_ovito, '(A)') 'ITEM: ATOMS x y z mux muy muz mu'
DO i = 1, nb
ploc = SQRT(pol_unit(1,i)*pol_unit(1,i) + pol_unit(2,i)*pol_unit(2,i) + pol_unit(3,i)*pol_unit(3,i) )
write(fileID_pol_ovito, 400) (b_atms(j,i), j=1,3), (pol_unit(j,i),j=1,3), ploc
ENDDO
400 FORMAT(1X, 3F15.10, 3F15.10, F15.10)
END SUBROUTINE write_pol_ovito

END MODULE writefiles
