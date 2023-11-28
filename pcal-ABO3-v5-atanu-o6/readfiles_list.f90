MODULE readfiles
 use constant, only : dp, fileID_in 
 use data_structure, only : latt_vec, a_atms, b_atms, o_atms, sp, tndump 
    IMPLICIT NONE 
private
public :: read_dumpxyz, readdump_type1,readdump_type2
contains

SUBROUTINE read_dumpxyz(na, nb, no)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: na, nb, no
 REAL(dp) :: x,y, start, finish
 INTEGER :: i, j, ntot
     !reading .xyz file

!WRITE(15, *) '# Px   Py,   Pz,  and  P = SQRT(px^2 +py^2 + pz^2)'
!!*********Reading files***********!!

   READ(fileID_in,*)
   READ(fileID_in,*) tndump
   READ(fileID_in,*)

   READ(fileID_in,*) ntot
   READ(fileID_in,*)
   DO i=1,3
     READ(fileID_in,*) x, y
     latt_vec(i,i)=ABS(y-x)
   ENDDO

   READ(fileID_in,*)
DO i=1, na
    READ (fileID_in,*) a_atms(1,i), a_atms(2,i), a_atms(3,i)
 ENDDO

 DO i=1, nb
  READ(fileID_in,*) b_atms(1,i), b_atms(2,i), b_atms(3,i)
  !, sp(1,i), sp(2,i), sp(3,i)
 ENDDO

DO i=1, no
   READ(fileID_in,*) o_atms(1,i), o_atms(2,i), o_atms(3,i)
 ENDDO

END SUBROUTINE read_dumpxyz 

SUBROUTINE readdump_type1(na, nb, no)
 IMPLICIT NONE 
 INTEGER, INTENT(IN) :: na, nb, no  
 REAL(dp) :: x,y, start, finish
 INTEGER :: i, j, ntot   
     !reading .xyz file

!WRITE(15, *) '# Px   Py,   Pz,  and  P = SQRT(px^2 +py^2 + pz^2)'
!!*********Reading files***********!!

   READ(fileID_in,*)
   READ(fileID_in,*) tndump
   READ(fileID_in,*)

   READ(fileID_in,*) ntot
   READ(fileID_in,*)
   DO i=1,3
     READ(fileID_in,*) x, y
     latt_vec(i,i)=ABS(y-x)
   ENDDO

   READ(fileID_in,*)
DO i=1, na
    READ (fileID_in,*) a_atms(1,i), a_atms(2,i), a_atms(3,i)
 ENDDO

 DO i=1, nb
  READ(fileID_in,*) b_atms(1,i), b_atms(2,i), b_atms(3,i), sp(1,i), sp(2,i), sp(3,i)
 ENDDO
  
DO i=1, no
   READ(fileID_in,*) o_atms(1,i), o_atms(2,i), o_atms(3,i)
 ENDDO

END SUBROUTINE readdump_type1

SUBROUTINE readdump_type2(na, nb, no)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: na, nb, no  
    REAL(dp) :: x,y, start, finish
    INTEGER :: i, j, ntot 
    INTEGER :: id !, type   
   !reading .xyz file
   !WRITE(15, *) '# Px   Py,   Pz,  and  P = SQRT(px^2 +py^2 + pz^2)'
   !!*********Reading files***********!!
      READ(fileID_in,*)
      READ(fileID_in,*) tndump
      READ(fileID_in,*)
      READ(fileID_in,*) ntot
      READ(fileID_in,*)
      DO i=1,3
        READ(fileID_in,*) x, y
        latt_vec(i,i)=ABS(y-x)
      ENDDO
   
      READ(fileID_in,*)
   DO i=1, na
       READ (fileID_in,*) id,  a_atms(1,i), a_atms(2,i), a_atms(3,i)
    ENDDO
   
    DO i=1, nb
     READ(fileID_in,*) id, b_atms(1,i), b_atms(2,i), b_atms(3,i), sp(1,i), sp(2,i), sp(3,i)
    ENDDO
     
   DO i=1, no
      READ(fileID_in,*) id,  o_atms(1,i), o_atms(2,i), o_atms(3,i)
    ENDDO
   
  
   !!*******************************!  
   END SUBROUTINE readdump_type2
   
END MODULE readfiles
