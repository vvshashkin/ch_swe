program test_standard_atmosphere
    use standard_atmosphere_mod, only : get_standard_atm_temp_z,   &
                                        get_standard_atm_temp_p,   &
                                        get_standard_atm_pres_z,   &
                                        get_standard_atm_height_p, &
                                        get_standard_atm_theta_z,  &
                                        get_standard_atm_theta_p

    implicit none

    real(kind=8) :: z(9) = [0.0_8, 10e3_8, 12e3_8, 17e3_8, 22e3_8, 32e3_8, 50e3_8, 60e3_8, 80e3_8]

    print *, "TEMP", get_standard_atm_temp_z(z(1:9))
    print *, "TEMP(p)-Temp", get_standard_atm_temp_p(get_standard_atm_pres_z(z(1:9))) - &
                            get_standard_atm_temp_z(z(1:9))
    print *, "Pres", get_standard_atm_pres_z(z(1:9))
    print *, "H(p)-H", get_standard_atm_height_p(get_standard_atm_pres_z(z(1:9)))-z(1:9)
    print *, "Theta", get_standard_atm_theta_z(z)


    print *, "T250", get_standard_atm_temp_p(250e2_8)
end program
