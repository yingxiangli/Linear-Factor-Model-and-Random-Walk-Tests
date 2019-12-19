# Linear Factor Model and Random Walk Tests
 Linear Factor Model and Random Walk Tests

This repository includes the julia code to run the following tests:

TEST OF LINEAR FACTOR MODEL (Chapter 12 of Cochrane (2000))
1) time-series test of linear factor models (ap_test_time)
   - Users are able to use iid, time iid, Newey-West and Hansen-Hodrick covariance matrix
   - Users are able to use asymptotic and finite-distribution OLS and GMM tests
2) cross-sectional test of linear factor models (ap_test_cross)
   - Users are able to use iid, time iid, Newey-West and Hansen-Hodrick covariance matrix
   - Users are able to use OLS, GLS, OLS and GLS with Shanken correction, GMM and GLS GMM
3) cross-sectional test of linear factor models with the Fama-MacBeth procedure (ap_test_fm)
* I have extended the tests to multi-factor case

TEST OF RANDOM WALK (Chapter 2 of Campbell, Ho and MacKinley (1997))
1) Box and Pierce (1970) test (rw_test)
2) Ljung and Box (1978) test (rw_test)
3) variance ratio test (rw_test)


The ap_rw_tests.ipynb runs those tests.
