"""
ap_test_time(re, factors; cov_method = "time iid", lag = 0, test_method = "asymptotic")

TIME-SERIES ASSET PRICING TESTS

Asset pricing models predict a restriction on the intercepts of the time-series regressions.
The asset pricing model is valid if the intercepts of time-series regressions are jointly zero.

Regress asset excess returns on factors to obtain betas

INPUTS
`re': T x N matrix of excess returns, where T is the number of periods for each test assets and N is the number of test assets
`factors': T x K matrix of factors. Factors must be excess returns for time-series tests
`cov_method': method for calculating the variance-covariance matrix of the error terms:
            - time iid: error terms are iid across time but cross-sectionally correlated
            - iid: error terms are iiid both across time and assets
            - NW: Newey-West standard errors with m lags
            - HH: Hansen-Hodrick standard errors with m lag
`lag': number of lags used in the Newey-West or Hansen-Hodrick standard errors
`test_method': method of the time-series asset pricing tests:
            - asymptotic: OLS asymptotic distribution
            - finite: OLS finite-sample distribution - Gibbons Ross and Shanken (1989) test statistics
            - GMM asymptotic: asymptotic GMM exploiting the moment conditions in time-series regressions
            - GMM finite: finite-distribution GMM exploiting the moment conditions in time-series regressions

OUTPUTS
`alpha': N x 1 vector of alphas of the test assets
`beta': N x K matrix of the betas of the test assets. There are N test assets and K factors
`resid': T x N matrix of the pricing errors
`cov_resid': N x N matrix of covariance matrix of the pricing errors

`test_statistics': Wald-statistics of the asset pricing tests
`pvalue': pvalue of the test statistics. For asymptotics, the test statistics has a Chi-square distribution. For finite sample, it has a F distribution

"""



function ap_test_time(re, factors, cov_method = "time iid", lag = 0, test_method = "asymptotic")
    (T, K) = size(factors);
    N = size(re, 2);
    Ef = mean(factors, dims = 1)'; # mean values of factors K x 1
    cov_f = cov(factors); # coviance matrix for the factors K x K

    X = [ones(T, 1) factors]; # T x (K+1)
    (alpha, beta, resid, cov_resid) = OLSFn(re, X, cov_method, lag); # time-series regressions


    # Methods of time series test
    if test_method == "asymptotic" # OLS Asymptotic
        test_statistics = T * ((1.0 .+ Ef' * (cov_f \ Ef)) \ (alpha')) * (cov_resid \ alpha);
        pvalue = 1.0 .- cdf.(Chisq(N), test_statistics);

    elseif test_method == "finite" # OLS Finite (Gibbons Ross and Shanken (1989))
        test_statistics = ((T - N - K) / N) * ((1.0 .+ Ef' * (cov_f \ Ef)) \ (alpha')) * (cov_resid \ alpha);
        pvalue = 1.0 .- cdf.(FDist(N, T - N - K), test_statistics);

    elseif test_method == "GMM asymptotic" || test_method == "GMM finite" # GMM Asymptotic
        h = [ones(T, 1) factors];
        f = kron(ones(1, K+1), resid) .* kron(h, ones(1, N)); # Moment Conditions - p.218 in Cochrane (2000)
        d = - kron(h' * h / T, I(N)); # Jacobian of the moment conditions
        dinv = d \ I; # Inverse of the Jacobian

        Ef2 = CovMoment(f, cov_method, lag); # spectral density matrix
        vp = (1.0 / T) * dinv * Ef2 * dinv'; # asymptotic covariance of the regression coefficiets - p.218 equation (172) in Cochrane  
        wt = vp[1:N, 1:N] \ I; # covariance matrix of the intercepts (alphas)
        test_statistics = alpha' * wt * alpha;
        pvalue = 1.0 .- cdf.(Chisq(N), test_statistics);

        if test_method == "GMM finite" # GMM Finite
            test_statistics = test_statistics * ((T - N - K) / N) * ((1.0 .+ Ef' * (cov_f \ Ef)) \ I) / (T * ((1.0 .+ Ef' * (cov_f \ Ef)) \ I));
            pvalue = 1.0 .- cdf.(FDist(N, T - N - K), test_statistics);
        end
    end

    return (alpha, beta, resid, cov_resid, test_statistics, pvalue)
end
