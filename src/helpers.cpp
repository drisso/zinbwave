# include <Rmath.h>
# include <Rcpp.h>

//'@useDynLib zinbwave
// [[Rcpp::export]]
Rcpp::NumericVector clog1pexp(Rcpp::NumericVector x){
                              // double c0 = -37.0,
                              // double c1 = 18.0,
                              // double c2 = 33.3) {

    Rcpp::NumericVector y(x.length());
    for(int i=0; i<x.length(); i++) {
        y[i] = log1pexp(x[i]);
    }
    return y;
//
//     Rcpp::NumericVector y(x.length());
//
//     for(int i=0; i<y.length(); i++) {
//         if(c0 < x[i] & x[i] <= c1) {
//             y[i] = std::log1p(std::exp(x[i]));
//         }
//         if(c1 < x[i] & x[i] <= c2) {
//             y[i] = x[i] + 1/std::exp(x[i]);
//         }
//         if(c2 < x[i]) {
//             y[i] = x[i];
//         }
//     }
//
//     return y;
}

