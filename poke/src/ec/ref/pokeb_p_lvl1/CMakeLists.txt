set(SOURCE_FILES_EC_${SVARIANT_UPPER}_REF
    ${B_ECX_DIR}/ec.c
    ${B_ECX_DIR}/kps.c
    ${B_ECX_DIR}/xisog.c
    ${B_ECX_DIR}/xeval.c
    ${B_ECX_DIR}/poly-mul.c
    ${B_ECX_DIR}/poly-redc.c
    ${ECX_DIR}/isog_chains.c
    ${B_ECX_DIR}/basis.c
    ${B_ECX_DIR}/biextension.c
    ${B_ECX_DIR}/dac.c
    ${B_ECX_DIR}/b_isog.c
)

add_library(${LIB_EC_${SVARIANT_UPPER}} ${SOURCE_FILES_EC_${SVARIANT_UPPER}_REF})
target_include_directories(${LIB_EC_${SVARIANT_UPPER}} PRIVATE ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_EC})
target_link_libraries(${LIB_EC_${SVARIANT_UPPER}} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_GF_${SVARIANT_UPPER}})

add_subdirectory(test)