gen_config(){
    TEMPLATE=$1
    N=$2
    DT=$3
    STAGGERING=$4
    SBP_ORDER=$5
    HMEAN=$6

    Ah_TEMPLATE="
    swm_operator_section = 'operator_conf'\n
    operator_conf = {\n
        swm_operator_type = 'linear_swm_operator',\n
        v_components_type = 'covariant',\n
        h_mean = %%%HMEAN,\n
        div_op_name = 'div_ah_sbp_proj_%%%SBP_ORDER',\n
        grad_op_name = 'grad_ah_sbp_proj_%%%SBP_ORDER',\n
        coriolis_op_name = 'coriolis_colocated',\n
        co2contra_op_name = 'co2contra_colocated',\n
        quadrature_name = 'SBP_Ah%%%SBP_ORDER_quadrature',\n
        helm_solver_cfg = {
            solver_name = 'cg',
            max_iter = 100, 
            rel_tol  = 0.00001 
            },\n
    },\n
    output_diag_section = 'output_diag_conf_Ah'\n
    output_diag_conf_Ah = {\n
        output_curl        = .false.,\n
        output_div         = .false.,\n
        curl_operator_name = 'curl_div_ah_sbp_proj_%%%SBP_ORDER',\n
        div_operator_name = 'div_ah_sbp_proj_%%%SBP_ORDER'\n
        co2contra_operator_name = 'co2contra_colocated',\n
        h2uv_operator_name = 'interp2d_p2uv_identity',\n
        quadrature_name          = 'SBP_Ah%%%SBP_ORDER_quadrature'
    }
    "
    Ch_TEMPLATE="swm_operator_section = 'operator_conf'\n
        operator_conf = { \n
        swm_operator_type        = 'linear_swm_operator'\n
        v_components_type        = 'covariant',\n
        h_mean                   = %%%HMEAN,\n
        div_op_name              = 'div_ch_sbp_sat_%%%SBP_ORDER_proj',\n
        grad_op_name             = 'grad_ch_sbp_sat_%%%SBP_ORDER',\n
        coriolis_op_name         = 'coriolis_Ch_sbp%%%SBP_ORDER',\n
        co2contra_op_name        = 'co2contra_ch_sbp%%%SBP_ORDER',\n
        quadrature_name          = 'SBP_C%%%SBP_ORDER_quadrature',\n
        helm_solver_cfg = {
            solver_name = 'cg',
            max_iter = 100, 
            rel_tol  = 0.00001 
            },\n
        },\n
        output_diag_section = 'output_diag_conf_Ch'\n
        output_diag_conf_Ch = {\n
            output_curl        = .false.,\n
            output_div         = .false.\n
            curl_operator_name = 'curl_ch_sbp_sat_%%%SBP_ORDER',\n
            div_operator_name = 'div_ch_sbp_sat_%%%SBP_ORDER_proj'\n
            co2contra_operator_name = 'co2contra_ch_sbp%%%SBP_ORDER',\n
            h2uv_operator_name = 'interp2d_p2uv_Ch_sbp%%%SBP_ORDER',\n
            uv2q_operator_name = 'interp2d_uv2qvec_Ch_sbp%%%SBP_ORDER'\n
            quadrature_name          = 'SBP_C%%%SBP_ORDER_quadrature'\n
        }\n
    "

    if [[ $STAGGERING == "Ah" ]]; then
        TEMPLATE=$TEMPLATE"\n"$Ah_TEMPLATE
    elif  [[ $STAGGERING == "Ch" ]]; then
        TEMPLATE=$TEMPLATE"\n"$Ch_TEMPLATE
    # elif  [[ $STAGGERING == "C" ]]; then
    #     TEMPLATE=$TEMPLATE"\n"$C_TEMPLATE
	fi

	echo -e $TEMPLATE |
             sed "s/%%%N/$N/" |
             sed "s/%%%DT/$DT/" |
             sed "s/%%%HMEAN/$HMEAN/" |
             sed "s/%%%STAGGERING/$STAGGERING/" |
             sed "s/%%%SBP_ORDER/$SBP_ORDER/"
}

