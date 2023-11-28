MODULE relative_to_O_disp
    use constant, only : dp 
    use data_structure, only : ndim3, dp, latt_vec, unitvol, &
    a_atms, o_atms, a_disp, a_pol, dmin_ao,ao12_neigh,&
    b_atms, b_disp, b_pol, dmin_bo,bo6_neigh  
    
    IMPLICIT NONE 
    private 
    public :: cal_a_site_pol_disp,cal_b_site_pol_disp
    contains
    
    

SUBROUTINE cal_a_site_pol_disp(qcharge, n1, n1st)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1st, n1
    REAL(dp), INTENT(IN) ::  qcharge
   
  
    REAL(dp) ::  atm_info(3,27), d(27)
    REAL(dp) :: ab_loc_coord(3, n1st, n1)
    REAL(dp), DIMENSION(:) :: v1(ndim3), vi(ndim3), vf(ndim3)
  
    REAL(dp) :: c_mass(3), x_sum, y_sum, z_sum, n1st_d
    REAL(dp) :: d_sq, dist_temp
  !
  
    INTEGER :: i, j,  n, l1, l2, l3, loc
  
   d_sq = dmin_ao*dmin_ao

     n1st_d = DBLE(n1st)
  
     DO i=1, n1
      DO j =1, n1st 
        ab_loc_coord(1, j, i) = o_atms(1, ao12_neigh(j,i))
        ab_loc_coord(2, j, i) = o_atms(2, ao12_neigh(j,i))
        ab_loc_coord(3, j, i) = o_atms(3, ao12_neigh(j,i))
      ENDDO 
    ENDDO
  
    DO i=1, n1
        v1(1)=a_atms(1,i)
        v1(2)=a_atms(2,i)
        v1(3)=a_atms(3,i)
   
        x_sum = 0.0;
        y_sum = 0.0;
        z_sum = 0.0;
  
        DO j=1, n1st
           vi(1) = ab_loc_coord(1, j, i)
           vi(2) = ab_loc_coord(2, j, i)
           vi(3) = ab_loc_coord(3, j, i)
  
          dist_temp = (v1(1)-vi(1))*(v1(1)-vi(1)) + &
                      (v1(2)-vi(2))*(v1(2)-vi(2)) + &
                      (v1(3)-vi(3))*(v1(3)-vi(3))
  
          IF (dist_temp > d_sq) THEN 
            n=0
             DO l1=-1, 1
                DO l2 =-1, 1
                  DO l3=-1, 1
                    n=n+1
                    vf(1)=vi(1)+ l1*latt_vec(1,1) !
                    vf(2)=vi(2)+ l2*latt_vec(2,2) !
                    vf(3)=vi(3)+ l3*latt_vec(3,3)
                    atm_info(1,n)=vf(1)
                    atm_info(2,n)=vf(2)
                    atm_info(3,n)=vf(3)
                    d(n)= (v1(1)-vf(1))*(v1(1)-vf(1)) + &
                          (v1(2)-vf(2))*(v1(2)-vf(2)) + &
                          (v1(3)-vf(3))*(v1(3)-vf(3))

                   ENDDO
                 ENDDO
               ENDDO
  
                loc=MINLOC(d, DIM=1)
                !ab_loc_coord(1,k,i) = INT(atm_info(1,loc))
               ! ab_loc_coord(1:3,j,i) = atm_info(1:3,loc)
                x_sum = x_sum + atm_info(1,loc)
                y_sum = y_sum + atm_info(2,loc)
                z_sum = z_sum +  atm_info(3,loc)
              ELSE
                x_sum = x_sum + vi(1)
                y_sum = y_sum + vi(2)
                z_sum = z_sum + vi(3)
          ENDIF
        ENDDO
  
                  ! Center of mass
               c_mass(1)=x_sum/n1st_d
               c_mass(2)=y_sum/n1st_d
               c_mass(3)=z_sum/n1st_d
         ! displacement
                a_disp(1,i)=v1(1)-c_mass(1)
                a_disp(2,i)=v1(2)-c_mass(2)
                a_disp(3,i)=v1(3)-c_mass(3)
  
                a_pol(1,i) = qcharge*a_disp(1,i)*16.02/unitvol
                a_pol(2,i) = qcharge*a_disp(2,i)*16.02/unitvol
                a_pol(3,i) = qcharge*a_disp(3,i)*16.02/unitvol
  !print*, a_pol(1,i), a_pol(2,i), a_pol(3,i) 
      ENDDO
  
  END SUBROUTINE cal_a_site_pol_disp
  
  
  SUBROUTINE cal_b_site_pol_disp(qcharge, n1, n1st)
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: n1st, n1
    REAL(dp), INTENT(IN) :: qcharge
   
  
    REAL(dp) ::  atm_info(3,27), d(27)
    REAL(dp) :: ab_loc_coord(3, n1st, n1)
    !
    REAL(dp), DIMENSION(:) :: v1(ndim3), vi(ndim3), vf(ndim3)
  
    REAL(dp) :: c_mass(3), x_sum, y_sum, z_sum, n1st_d
    REAL(dp) :: d_sq, dist_temp
  
    INTEGER :: i, j,  n, l1, l2, l3, loc
  
  
    d_sq = dmin_bo*dmin_bo 

     n1st_d = DBLE(n1st)
  
     DO i=1, n1
      DO j =1, n1st 
        ab_loc_coord(1, j, i) = o_atms(1, bo6_neigh(j,i))
        ab_loc_coord(2, j, i) = o_atms(2, bo6_neigh(j,i))
        ab_loc_coord(3, j, i) = o_atms(3, bo6_neigh(j,i))
      ENDDO 
    ENDDO
  
    DO i=1, n1
        v1(1)=b_atms(1,i)
        v1(2)=b_atms(2,i)
        v1(3)=b_atms(3,i)
   
        x_sum = 0.0;
        y_sum = 0.0;
        z_sum = 0.0;
  
        DO j=1, n1st
           vi(1) = ab_loc_coord(1, j, i)
           vi(2) = ab_loc_coord(2, j, i)
           vi(3) = ab_loc_coord(3, j, i)
  
  
          dist_temp = (v1(1)-vi(1))*(v1(1)-vi(1)) + &
                      (v1(2)-vi(2))*(v1(2)-vi(2)) + &
                      (v1(3)-vi(3))*(v1(3)-vi(3))
         
          IF (dist_temp > d_sq ) THEN 
            n=0
             DO l1=-1, 1
                DO l2 =-1, 1
                  DO l3=-1, 1
                    n=n+1
                    vf(1)=vi(1)+ l1*latt_vec(1,1) !
                    vf(2)=vi(2)+ l2*latt_vec(2,2) !
                    vf(3)=vi(3)+ l3*latt_vec(3,3)
                    atm_info(1,n)=vf(1)
                    atm_info(2,n)=vf(2)
                    atm_info(3,n)=vf(3)
                    d(n)= (v1(1)-vf(1))*(v1(1)-vf(1)) + &
                              (v1(2)-vf(2))*(v1(2)-vf(2)) + &
                              (v1(3)-vf(3))*(v1(3)-vf(3))
                   ENDDO
                 ENDDO
               ENDDO
  
                loc=MINLOC(d, DIM=1)
                !ab_loc_coord(1,k,i) = INT(atm_info(1,loc))
               ! ab_loc_coord(1:3,j,i) = atm_info(1:3,loc)
                x_sum = x_sum + atm_info(1,loc)
                y_sum = y_sum + atm_info(2,loc)
                z_sum = z_sum +  atm_info(3,loc)
              ELSE
                x_sum = x_sum + vi(1)
                y_sum = y_sum + vi(2)
                z_sum = z_sum + vi(3)
          ENDIF
        ENDDO
  
                  ! Center of mass
               c_mass(1)=x_sum/n1st_d
               c_mass(2)=y_sum/n1st_d
               c_mass(3)=z_sum/n1st_d
         ! displacement
                b_disp(1,i)=v1(1)-c_mass(1)
                b_disp(2,i)=v1(2)-c_mass(2)
                b_disp(3,i)=v1(3)-c_mass(3)
  
                b_pol(1,i) = qcharge*b_disp(1,i)*16.02/unitvol
                b_pol(2,i) = qcharge*b_disp(2,i)*16.02/unitvol
                b_pol(3,i) = qcharge*b_disp(3,i)*16.02/unitvol
  
      ENDDO
  
  END SUBROUTINE cal_b_site_pol_disp
  
END MODULE relative_to_O_disp