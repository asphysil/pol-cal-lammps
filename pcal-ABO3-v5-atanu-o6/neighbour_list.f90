MODULE nearest_neighbour_list
use constant, only : dp 
use data_structure, only : ndim3, dp, latt_vec
private 
public :: nearest_neighbour
contains 

SUBROUTINE nearest_neighbour(n1st, n1, n2, dab, atm1, atm2, ab_loc)
 
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1st, n1, n2
    REAL(dp), INTENT(IN) :: dab, atm1(ndim3,n1), atm2(ndim3,n2)
    INTEGER, INTENT(OUT) :: ab_loc(n1st,n1)
    !
    INTEGER :: atm_info(27)
    REAL(dp) ::  d(27), dist, d_sq, dmax
  
    REAL(dp) :: v1(ndim3), vi(ndim3), vf(ndim3)
  !
    INTEGER :: ncount 
    INTEGER :: i, j,  n, l1, l2, l3, loc
  
     dmax=100 ! angs
     
  d_sq = dab*dab 

    DO i=1, n1
        v1(1)=atm1(1,i)
        v1(2)=atm1(2,i)
        v1(3)=atm1(3,i)
        ncount =0 
        DO j=1, n2
          if ( ncount == n1st ) exit

           vi(1)=atm2(1,j)
           vi(2)=atm2(2,j)
           vi(3)=atm2(3,j)
  
           dist = (v1(1)-vi(1))**2 + &
                  (v1(2)-vi(2))**2 + &
                  (v1(3)-vi(3))**2

           if (dist < d_sq) then
            ncount = ncount + 1
            ab_loc(ncount,i) = j 
            !print*, dist 
            cycle 
           else 
           n=0
           DO l1=-1, 1
              DO l2 =-1, 1
                 DO l3=-1, 1
                   n=n+1
                    vf(1)=vi(1)+ l1*latt_vec(1,1) !
                    vf(2)=vi(2)+ l2*latt_vec(2,2) !
                    vf(3)=vi(3)+ l3*latt_vec(3,3)
  
                    atm_info(n)=j
                    d(n)=(v1(1)-vf(1))*(v1(1)-vf(1)) + &
                              (v1(2)-vf(2))*(v1(2)-vf(2)) + &
                              (v1(3)-vf(3))*(v1(3)-vf(3))
                   ENDDO
                 ENDDO
               ENDDO
               loc=MINLOC(d, DIM=1)
               dist = d(loc)
               if (dist < d_sq) then
                 ncount = ncount + 1
                 ab_loc(ncount,i)=atm_info(loc)
               endif 
          endif

          ENDDO 
      ENDDO
  
  END SUBROUTINE nearest_neighbour

  END MODULE nearest_neighbour_list
  