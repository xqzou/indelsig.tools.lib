// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// segment
Rcpp::DataFrame segment(std::vector<std::string> string, std::vector<std::string> context);
RcppExport SEXP _indelsig_tools_lib_segment(SEXP stringSEXP, SEXP contextSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type string(stringSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type context(contextSEXP);
    rcpp_result_gen = Rcpp::wrap(segment(string, context));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_indelsig_tools_lib_segment", (DL_FUNC) &_indelsig_tools_lib_segment, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_indelsig_tools_lib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
