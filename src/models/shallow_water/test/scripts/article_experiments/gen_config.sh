gen_config(){
    TEMPLATE=$1
    N=$2
    DT=$3
    STAGGERING=$4
    # staggering=echo $4 | tr '[:upper:]' '[:lower:]'
    SBP_ORDER=$5
    UPSTREAM_ORDER=$6

    Ah_TEMPLATE="
    swm_operator_section = 'operator_conf'\n
    operator_conf = {\n
        swm_operator_type = 'vector_invariant_swm_operator',\n
        v_components_type = 'covariant',\n
        div_op_name = 'div_ah_sbp_proj_%%%SBP_ORDER',\n
        grad_op_name = 'grad_ah_sbp_proj_%%%SBP_ORDER',\n
        coriolis_op_name = 'coriolis_colocated',\n
        curl_op_name = 'curl_div_ah_sbp_proj_%%%SBP_ORDER',\n
        KE_op_name = 'KE_colocated',\n
        co2contra_op_name = 'co2contra_colocated',\n
        massflux_op_name = 'massflux_colocated',\n
        quadrature_name = 'SBP_Ah%%%SBP_ORDER_quadrature',\n
    },\n
    diffusion_section = 'diffusion_Ah'\n
    diffusion_Ah = {\n
        diff_time_scheme = 'explicit_Eul1'\n
        uv_diff_coeff = 0.03,\n
        hordiff_uv_name = 'hordiff_vec_xyz_Ah_sbp_%%%SBP_ORDER_narrow',\n
        h_diff_coeff = 0.01,\n
        hordiff_h_name = 'hordiff_scalar_Ah_sbp_%%%SBP_ORDER_narrow',\n
    },\n
    output_diag_section = 'output_diag_conf_Ah'\n
    output_diag_conf_Ah = {\n
        output_curl        = .true.,\n
        output_div         = .true.,\n
        curl_operator_name = 'curl_div_ah_sbp_proj_%%%SBP_ORDER',\n
        div_operator_name = 'div_ah_sbp_proj_%%%SBP_ORDER'\n
        co2contra_operator_name = 'co2contra_colocated',\n
        h2uv_operator_name = 'interp2d_p2uv_identity',\n
        quadrature_name          = 'SBP_Ah%%%SBP_ORDER_quadrature'\n
    }
    "
    Ch_TEMPLATE="swm_operator_section = 'operator_conf'\n
        operator_conf = { \n
        swm_operator_type        = 'advective_swm_operator'\n
        v_components_type        = 'covariant',\n
        div_op_name              = 'div_ch_sbp_sat_%%%SBP_ORDER_proj',\n
        grad_op_name             = 'grad_ch_sbp_sat_%%%SBP_ORDER',\n
        coriolis_op_name         = 'coriolis_Ch_sbp%%%SBP_ORDER',\n
        co2contra_op_name        = 'co2contra_ch_sbp%%%SBP_ORDER',\n
        vector_advection_op_name = 'vector_advection_Ch_%%%UPSTREAM_ORDER_cov'\n
        massflux_op_name         = 'massflux_ch_%%%UPSTREAM_ORDER',\n
        quadrature_name          = 'SBP_C%%%SBP_ORDER_quadrature',\n
        },\n
        output_diag_section = 'output_diag_conf_Ch'\n
        output_diag_conf_Ch = {\n
            output_curl        = .true.,\n
            output_div         = .true.\n
            curl_operator_name = 'curl_ch_sbp_sat_%%%SBP_ORDER',\n
            div_operator_name = 'div_ch_sbp_sat_%%%SBP_ORDER_proj'\n
            co2contra_operator_name = 'co2contra_ch_sbp%%%SBP_ORDER',\n
            h2uv_operator_name = 'interp2d_p2uv_Ch_sbp%%%SBP_ORDER',\n
            uv2q_operator_name = 'interp2d_uv2qvec_Ch_sbp%%%SBP_ORDER'\n
            quadrature_name          = 'SBP_C%%%SBP_ORDER_quadrature'\n
        }\n
    "
    C_TEMPLATE="swm_operator_section = 'operator_conf'\n
        operator_conf = { \n
        swm_operator_type        = 'advective_swm_operator'\n
        v_components_type        = 'contravariant',\n
        div_op_name              = 'div_c_sbp_sat_%%%SBP_ORDER',\n
        grad_op_name             = 'grad_c_sbp_sat_%%%SBP_ORDER',\n
        coriolis_op_name         = 'coriolis_Cgrid_noncons_sbp%%%SBP_ORDER',\n
        co2contra_op_name        = ''co2contra_c_sbp%%%SBP_ORDER_new',\n
        vector_advection_op_name = 'vector_advection_C_up4_cov'\n
        massflux_op_name         = 'massflux_c_up4',\n
        quadrature_name          = 'SBP_C%%%SBP_ORDER_quadrature',\n
        },\n
        output_diag_section = 'output_diag_conf_C'\n
        output_diag_conf_C = {\n
            output_curl        = .true.,\n
            output_div         = .true.\n
            curl_operator_name = 'curl_c_sbp_sat_%%%SBP_ORDER',\n
            div_operator_name = 'div_c_sbp_sat_%%%SBP_ORDER'\n
            co2contra_operator_name = 'co2contra_c_sbp%%%SBP_ORDER_new',\n
            h2uv_operator_name = 'interp2d_p2uv_C_sbp%%%SBP_ORDER',\n
            uv2q_operator_name = 'interp2d_uv2qvec_C_sbp%%%SBP_ORDER'\n
            quadrature_name          = 'SBP_C%%%SBP_ORDER_quadrature'\n
        }
    "

    if [[ $STAGGERING == "Ah" ]]; then
        TEMPLATE=$TEMPLATE"\n"$Ah_TEMPLATE
    elif  [[ $STAGGERING == "Ch" ]]; then
        TEMPLATE=$TEMPLATE"\n"$Ch_TEMPLATE
    elif  [[ $STAGGERING == "C" ]]; then
        TEMPLATE=$TEMPLATE"\n"$C_TEMPLATE
	fi

	echo -e $TEMPLATE |
             sed "s/%%%N/$N/" |
             sed "s/%%%DT/$DT/" |
             sed "s/%%%STAGGERING/$STAGGERING/" |
             sed "s/%%%SBP_ORDER/$SBP_ORDER/"   |
             sed "s/%%%UPSTREAM_ORDER/$UPSTREAM_ORDER/"
}

