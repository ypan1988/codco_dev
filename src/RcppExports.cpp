// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// obj_fun
double obj_fun(const arma::mat& A, const arma::vec& U2);
RcppExport SEXP _codco_obj_fun(SEXP ASEXP, SEXP U2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type U2(U2SEXP);
    rcpp_result_gen = Rcpp::wrap(obj_fun(A, U2));
    return rcpp_result_gen;
END_RCPP
}
// c_obj_fun
double c_obj_fun(const arma::mat& M, const arma::vec& C);
RcppExport SEXP _codco_c_obj_fun(SEXP MSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(c_obj_fun(M, C));
    return rcpp_result_gen;
END_RCPP
}
// gen_m
arma::mat gen_m(const arma::mat& X, const arma::mat& A);
RcppExport SEXP _codco_gen_m(SEXP XSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(gen_m(X, A));
    return rcpp_result_gen;
END_RCPP
}
// gen_u
arma::vec gen_u(const arma::mat& M, const arma::mat& X, const arma::vec& C);
RcppExport SEXP _codco_gen_u(SEXP MSEXP, SEXP XSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_u(M, X, C));
    return rcpp_result_gen;
END_RCPP
}
// remove_one
double remove_one(const arma::mat& A, arma::uword i, const arma::vec& u);
RcppExport SEXP _codco_remove_one(SEXP ASEXP, SEXP iSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::uword >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(remove_one(A, i, u));
    return rcpp_result_gen;
END_RCPP
}
// remove_one_mat
arma::mat remove_one_mat(const arma::mat& A, arma::uword i);
RcppExport SEXP _codco_remove_one_mat(SEXP ASEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::uword >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(remove_one_mat(A, i));
    return rcpp_result_gen;
END_RCPP
}
// add_one
double add_one(const arma::mat& A, double sigma_jj, const arma::vec& f, const arma::vec& u);
RcppExport SEXP _codco_add_one(SEXP ASEXP, SEXP sigma_jjSEXP, SEXP fSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type sigma_jj(sigma_jjSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type f(fSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(add_one(A, sigma_jj, f, u));
    return rcpp_result_gen;
END_RCPP
}
// add_one_mat
arma::mat add_one_mat(const arma::mat& A, double sigma_jj, const arma::vec& f);
RcppExport SEXP _codco_add_one_mat(SEXP ASEXP, SEXP sigma_jjSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type sigma_jj(sigma_jjSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(add_one_mat(A, sigma_jj, f));
    return rcpp_result_gen;
END_RCPP
}
// ChooseSwap
double ChooseSwap(arma::uvec idx_in, arma::mat A, arma::mat sig, arma::vec u, arma::uvec& out2, arma::mat& out3);
RcppExport SEXP _codco_ChooseSwap(SEXP idx_inSEXP, SEXP ASEXP, SEXP sigSEXP, SEXP uSEXP, SEXP out2SEXP, SEXP out3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type idx_in(idx_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type out2(out2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type out3(out3SEXP);
    rcpp_result_gen = Rcpp::wrap(ChooseSwap(idx_in, A, sig, u, out2, out3));
    return rcpp_result_gen;
END_RCPP
}
// Grad
arma::uvec Grad(arma::uvec idx_in, arma::mat A, arma::mat sig, arma::vec u, double tol, bool trace);
RcppExport SEXP _codco_Grad(SEXP idx_inSEXP, SEXP ASEXP, SEXP sigSEXP, SEXP uSEXP, SEXP tolSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type idx_in(idx_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(Grad(idx_in, A, sig, u, tol, trace));
    return rcpp_result_gen;
END_RCPP
}
// ChooseSwapRobust
double ChooseSwapRobust(arma::uword nlist, arma::uvec idx_in, const arma::mat& A_list, const arma::mat& sig_list, const arma::vec& u_list, const arma::vec& weights, arma::uvec& out2, arma::mat& out3);
RcppExport SEXP _codco_ChooseSwapRobust(SEXP nlistSEXP, SEXP idx_inSEXP, SEXP A_listSEXP, SEXP sig_listSEXP, SEXP u_listSEXP, SEXP weightsSEXP, SEXP out2SEXP, SEXP out3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type nlist(nlistSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type idx_in(idx_inSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A_list(A_listSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sig_list(sig_listSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u_list(u_listSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type out2(out2SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type out3(out3SEXP);
    rcpp_result_gen = Rcpp::wrap(ChooseSwapRobust(nlist, idx_in, A_list, sig_list, u_list, weights, out2, out3));
    return rcpp_result_gen;
END_RCPP
}
// GradRobust
arma::uvec GradRobust(arma::uword nlist, arma::uvec idx_in, arma::mat A_list, arma::mat M_list, arma::vec C_list, arma::mat X_list, arma::mat sig_list, arma::vec u_list, arma::vec weights, double tol, bool trace);
RcppExport SEXP _codco_GradRobust(SEXP nlistSEXP, SEXP idx_inSEXP, SEXP A_listSEXP, SEXP M_listSEXP, SEXP C_listSEXP, SEXP X_listSEXP, SEXP sig_listSEXP, SEXP u_listSEXP, SEXP weightsSEXP, SEXP tolSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type nlist(nlistSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type idx_in(idx_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A_list(A_listSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type M_list(M_listSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type C_list(C_listSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig_list(sig_listSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u_list(u_listSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(GradRobust(nlist, idx_in, A_list, M_list, C_list, X_list, sig_list, u_list, weights, tol, trace));
    return rcpp_result_gen;
END_RCPP
}
// uvec_minus
arma::uvec uvec_minus(const arma::uvec& v, arma::uword rm_idx);
RcppExport SEXP _codco_uvec_minus(SEXP vSEXP, SEXP rm_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type rm_idx(rm_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(uvec_minus(v, rm_idx));
    return rcpp_result_gen;
END_RCPP
}
// GradRobustStep
Rcpp::List GradRobustStep(arma::uvec idx_in, arma::vec C_list, arma::mat X_list, arma::mat sig_list, arma::vec weights);
RcppExport SEXP _codco_GradRobustStep(SEXP idx_inSEXP, SEXP C_listSEXP, SEXP X_listSEXP, SEXP sig_listSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type idx_in(idx_inSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type C_list(C_listSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig_list(sig_listSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(GradRobustStep(idx_in, C_list, X_list, sig_list, weights));
    return rcpp_result_gen;
END_RCPP
}
// GradRobustAlg1
Rcpp::List GradRobustAlg1(arma::uvec idx_in, arma::vec C_list, arma::mat X_list, arma::mat sig_list, arma::vec weights);
RcppExport SEXP _codco_GradRobustAlg1(SEXP idx_inSEXP, SEXP C_listSEXP, SEXP X_listSEXP, SEXP sig_listSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type idx_in(idx_inSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type C_list(C_listSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_list(X_listSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sig_list(sig_listSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(GradRobustAlg1(idx_in, C_list, X_list, sig_list, weights));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _codco_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _codco_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _codco_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _codco_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_codco_obj_fun", (DL_FUNC) &_codco_obj_fun, 2},
    {"_codco_c_obj_fun", (DL_FUNC) &_codco_c_obj_fun, 2},
    {"_codco_gen_m", (DL_FUNC) &_codco_gen_m, 2},
    {"_codco_gen_u", (DL_FUNC) &_codco_gen_u, 3},
    {"_codco_remove_one", (DL_FUNC) &_codco_remove_one, 3},
    {"_codco_remove_one_mat", (DL_FUNC) &_codco_remove_one_mat, 2},
    {"_codco_add_one", (DL_FUNC) &_codco_add_one, 4},
    {"_codco_add_one_mat", (DL_FUNC) &_codco_add_one_mat, 3},
    {"_codco_ChooseSwap", (DL_FUNC) &_codco_ChooseSwap, 6},
    {"_codco_Grad", (DL_FUNC) &_codco_Grad, 6},
    {"_codco_ChooseSwapRobust", (DL_FUNC) &_codco_ChooseSwapRobust, 8},
    {"_codco_GradRobust", (DL_FUNC) &_codco_GradRobust, 11},
    {"_codco_uvec_minus", (DL_FUNC) &_codco_uvec_minus, 2},
    {"_codco_GradRobustStep", (DL_FUNC) &_codco_GradRobustStep, 5},
    {"_codco_GradRobustAlg1", (DL_FUNC) &_codco_GradRobustAlg1, 5},
    {"_codco_rcpparma_hello_world", (DL_FUNC) &_codco_rcpparma_hello_world, 0},
    {"_codco_rcpparma_outerproduct", (DL_FUNC) &_codco_rcpparma_outerproduct, 1},
    {"_codco_rcpparma_innerproduct", (DL_FUNC) &_codco_rcpparma_innerproduct, 1},
    {"_codco_rcpparma_bothproducts", (DL_FUNC) &_codco_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_codco(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
