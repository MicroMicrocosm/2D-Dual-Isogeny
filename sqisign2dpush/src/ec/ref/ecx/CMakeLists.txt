set(ECX_SOURCE_FILES
    ${ECX_DIR}/poly-mul.c
    ${ECX_DIR}/poly-redc.c
    ${ECX_DIR}/ec.c
    ${ECX_DIR}/ec_dual.c
    ${ECX_DIR}/curve.c
    ${ECX_DIR}/tedwards.c
    ${ECX_DIR}/kps.c
    ${ECX_DIR}/xisog.c
    ${ECX_DIR}/xeval.c
    ${ECX_DIR}/isog_chains.c
    ${ECX_DIR}/basis.c
)

function(sqisign_add_ec_library variant_upper variant_lower)
    add_library(${LIB_EC_${variant_upper}} ${ECX_SOURCE_FILES})
    target_include_directories(${LIB_EC_${variant_upper}} PRIVATE
        ${INC_INTBIG}
        ${INC_PRECOMP_${variant_upper}}
        ${INC_PUBLIC}
        ${INC_GF_${variant_upper}}
        ${INC_COMMON}
        ${INC_EC}
    )
    target_compile_definitions(${LIB_EC_${variant_upper}} PUBLIC SQISIGN_VARIANT=${variant_lower})
    target_compile_options(${LIB_EC_${variant_upper}} PRIVATE ${C_OPT_FLAGS})

    if(ARGC GREATER 3 AND ARGV3)
        target_link_libraries(${LIB_EC_${variant_upper}} PUBLIC ${LIB_GF_${variant_upper}})
    endif()
endfunction()