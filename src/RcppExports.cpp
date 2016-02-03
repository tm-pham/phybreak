// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// CCphyloconstruct
std::vector<int> CCphyloconstruct(std::vector<int> pars, std::vector<int> dims);
RcppExport SEXP phybreakCO_CCphyloconstruct(SEXP parsSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<int> >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type dims(dimsSEXP);
    __result = Rcpp::wrap(CCphyloconstruct(pars, dims));
    return __result;
END_RCPP
}
// CCphylotree
std::vector<int> CCphylotree(std::vector<int> pars, std::vector<int> dims);
RcppExport SEXP phybreakCO_CCphylotree(SEXP parsSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<int> >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type dims(dimsSEXP);
    __result = Rcpp::wrap(CCphylotree(pars, dims));
    return __result;
END_RCPP
}
// CCtranstreeconstruct
std::vector<int> CCtranstreeconstruct(std::vector<int> pars, std::vector<int> dims);
RcppExport SEXP phybreakCO_CCtranstreeconstruct(SEXP parsSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<int> >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type dims(dimsSEXP);
    __result = Rcpp::wrap(CCtranstreeconstruct(pars, dims));
    return __result;
END_RCPP
}
// CCtranstree
std::vector<int> CCtranstree(std::vector<int> pars, std::vector<int> dims);
RcppExport SEXP phybreakCO_CCtranstree(SEXP parsSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<int> >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type dims(dimsSEXP);
    __result = Rcpp::wrap(CCtranstree(pars, dims));
    return __result;
END_RCPP
}
// likseq
double likseq(CharacterVector SNPs, IntegerVector SNPfreqs, IntegerVector nodeparents, NumericVector nodetimes, double mutrate, int obs);
RcppExport SEXP phybreakCO_likseq(SEXP SNPsSEXP, SEXP SNPfreqsSEXP, SEXP nodeparentsSEXP, SEXP nodetimesSEXP, SEXP mutrateSEXP, SEXP obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< CharacterVector >::type SNPs(SNPsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type SNPfreqs(SNPfreqsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nodeparents(nodeparentsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nodetimes(nodetimesSEXP);
    Rcpp::traits::input_parameter< double >::type mutrate(mutrateSEXP);
    Rcpp::traits::input_parameter< int >::type obs(obsSEXP);
    __result = Rcpp::wrap(likseq(SNPs, SNPfreqs, nodeparents, nodetimes, mutrate, obs));
    return __result;
END_RCPP
}
// likseqenv
double likseqenv(Environment pbenv, IntegerVector nodestochange, IntegerVector tips);
RcppExport SEXP phybreakCO_likseqenv(SEXP pbenvSEXP, SEXP nodestochangeSEXP, SEXP tipsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Environment >::type pbenv(pbenvSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nodestochange(nodestochangeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type tips(tipsSEXP);
    __result = Rcpp::wrap(likseqenv(pbenv, nodestochange, tips));
    return __result;
END_RCPP
}
// MLphylotree_MCC
std::vector<double> MLphylotree_MCC(std::vector<int> pars, std::vector<int> dims);
RcppExport SEXP phybreakCO_MLphylotree_MCC(SEXP parsSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<int> >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type dims(dimsSEXP);
    __result = Rcpp::wrap(MLphylotree_MCC(pars, dims));
    return __result;
END_RCPP
}
// ptr
std::vector<int> ptr(IntegerVector pars, int ID);
RcppExport SEXP phybreakCO_ptr(SEXP parsSEXP, SEXP IDSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< int >::type ID(IDSEXP);
    __result = Rcpp::wrap(ptr(pars, ID));
    return __result;
END_RCPP
}
// sct
NumericVector sct(NumericVector tle);
RcppExport SEXP phybreakCO_sct(SEXP tleSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type tle(tleSEXP);
    __result = Rcpp::wrap(sct(tle));
    return __result;
END_RCPP
}