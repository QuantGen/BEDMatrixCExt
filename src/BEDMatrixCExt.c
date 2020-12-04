#include "BEDMatrix.h" // requires LinkingTo: BEDMatrix in DESCRIPTION

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP numSamples(SEXP X) {
    struct BEDMatrix *state = R_ExternalPtrAddr(X);
    if (state == NULL) {
        Rf_error("BEDMatrix instance has been unmapped.");
    }
    return Rf_ScalarInteger(state->num_samples);
}

SEXP numVariants(SEXP X) {
    struct BEDMatrix *state = R_ExternalPtrAddr(X);
    if (state == NULL) {
        Rf_error("BEDMatrix instance has been unmapped.");
    }
    return Rf_ScalarInteger(state->num_variants);
}

SEXP extractGenotypeCartesian(SEXP X, SEXP i, SEXP j) {
    struct BEDMatrix *state = R_ExternalPtrAddr(X);
    if (state == NULL) {
        Rf_error("BEDMatrix instance has been unmapped.");
    }
    R_xlen_t ni = Rf_xlength(i);
    int *pi = INTEGER(i);
    R_xlen_t nj = Rf_xlength(j);
    int *pj = INTEGER(j);
    SEXP out = PROTECT(Rf_allocMatrix(INTSXP, ni, nj));
    int *pout = INTEGER(out);
    int num_bytes_per_variant = compute_num_bytes_per_variant(state->num_samples);
    for (R_xlen_t col_idx = 0; col_idx < nj; col_idx++) {
        for (R_xlen_t row_idx = 0; row_idx < ni; row_idx++) {
            int genotype = extract_genotype_cartesian(
                state->data,
                pi[row_idx] - 1,
                pj[col_idx] - 1,
                num_bytes_per_variant
            );
            pout[col_idx * ni + row_idx] = recode_genotype(genotype, NA_INTEGER);
        }
    }
    return out;
}

SEXP extractGenotypeLinear(SEXP X, SEXP k) {
    struct BEDMatrix *state = R_ExternalPtrAddr(X);
    if (state == NULL) {
        Rf_error("BEDMatrix instance has been unmapped.");
    }
    R_xlen_t nk = Rf_xlength(k);
    int *pk = INTEGER(k);
    SEXP out = PROTECT(Rf_allocVector(INTSXP, nk));
    int *pout = INTEGER(out);
    int num_bytes_per_variant = compute_num_bytes_per_variant(state->num_samples);
    for (R_xlen_t idx = 0; idx < nk; idx++) {
        int genotype = extract_genotype_linear(
            state->data,
            pk[idx] - 1,
            state->num_samples,
            num_bytes_per_variant
        );
        pout[idx] = recode_genotype(genotype, NA_INTEGER);
    }
    return out;
}

static const R_CallMethodDef callMethods[] = {
    {"numSamples", (DL_FUNC) &numSamples, 1},
    {"numVariants", (DL_FUNC) &numVariants, 1},
    {"extractGenotypeCartesian", (DL_FUNC) &extractGenotypeCartesian, 1},
    {"extractGenotypeLinear", (DL_FUNC) &extractGenotypeLinear, 1},
    {NULL, NULL, 0}
};

void R_init_BEDMatrixCExt(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
