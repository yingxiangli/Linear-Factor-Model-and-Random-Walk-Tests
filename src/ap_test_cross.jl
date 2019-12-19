"""
ap_test_cross(re, factors; cov_method = "time iid", lag = 0, test_method = "Shanken")

CROSS-SECTIONAL ASSET PRICING TESTS

Asset pricing models predict a restriction on the error terms of the second-stage cross-sectional regressions.
The asset pricing model is valid if the error terms of cross-sectional regressions are jointly zero-mean.

Two-Pass Regression Estimation Procedure
- First-stage: time-series regressions of excess returns on factors to obtain estimates of betas for each assets
- Second-stage: cross-sectional regressions (without intercept) of time-series average returns on estiamted betas to obtain risk premium for each assets

INPUTS
`re': T x N matrix of excess returns, where T is the number of periods for each test assets and N is the number of test assets
`factors': T x K matrix of factors. Factors must be excess returns for time-series tests
`cov_method': method for calculating the variance-covariance matrix of the error terms:
            - time iid: error terms are iid across time but cross-sectionally correlated
            - iid: error terms are iiid both across time and assets
            - NW: Newey-West standard errors with m lag
            - HH: Hansen-Hodrick standard errors with m lags
`lag': number of lags used in the Newey-West or Hansen-Hodrick standard errors
`test_method': method of the time-series asset pricing tests:
            - OLS: OLS
            - GLS: GLS
            - Shanken: OLS with Shanken Correction (1992)
            - GLS Shanken: GLS with Shanken Correction (1992)
            - GMM: GMM exploiting the moment conditions from the first-stage time-series regressions and second-stage cross-sectional regressions. The model is overidentified
            - GLS GMM: GLS GMM exploiting the moment conditions from the first-stage time-series regressions and second-stage cross-sectional regressions. The model is overidentified

OUTPUTS
`Ere' : N x 1 vector of expected excess returns for N test assets
`lambda': K x 1 vector of risk premiums for K factors
`alpha': N x 1 vector of alphas for N test assets

`test_statistics': Wald-statistics of the asset pricing tests
`pvalue': pvalue of the test statistics. The test statistics has a Chi-square distribution.


"""


function ap_test_cross(re, factors, cov_method = "time iid", lag = 0, test_method = "Shanken")
    (T, K) = size(factors);
    N = size(re, 2);
    cov_f = cov(factors); # coviance matrix for the factors K x K

    X = [ones(T,1) factors];
    (~, beta, resid, cov_resid) = OLSFn(re, X, cov_method, lag); # First-stage time-series regressions

    Ere = mean(re, dims = 1)'; # N x 1 vector of sample average excess returns for each test assets

    # Methods of cross-sectional tests
    if test_method == "OLS" # OLS
        lambda = (beta' * beta) \ (beta' * Ere); # OLS second stage regressions - Risk premiums K x 1
        alpha = Ere- beta * lambda; # Alaphas N x 1
        Q = I - beta * ((beta' * beta) \ beta');
        cov_alpha = 1.0 ./ T * Q * cov_resid * Q;
        test_statistics = alpha' * (cov_alpha \ alpha);
        pvalue = 1.0 .- cdf.(Chisq(N - K), test_statistics);

    elseif test_method == "GLS" # GLS
        lambda = (beta' * (cov_resid \ beta)) \ beta' * (cov_resid \ Ere); # GLS second stage regressions - Risk premiums K x 1
        alpha = Ere- beta * lambda; # Alaphas N x 1
        test_statistics = T * alpha' * (cov_resid \ alpha);
        pvalue = 1.0 .- cdf.(Chisq(N - K), test_statistics);

    elseif test_method == "Shanken" # OLS Shanken Correction (1992)
        lambda = beta' * beta \ (beta' * Ere); # OLS second stage regressions - Risk premiums K x 1
        alpha = Ere- beta * lambda; # Alaphas N x 1
        Q = I - beta * ((beta' * beta) \ beta');
        cov_alpha = 1.0 ./ T * Q * cov_resid * Q .* (1.0 .+ lambda' * (cov_f \ lambda));
        test_statistics = alpha' * (cov_alpha \ alpha);
        pvalue = 1.0 .- cdf.(Chisq(N - K), test_statistics);

    elseif test_method == "GLS Shanken" # GLS Shanken Correction (1992)
        lambda = (beta' * (cov_resid \ beta)) \ beta' * (cov_resid \ Ere); # GLS second stage regressions - Risk premiums K x 1
        alpha = Ere - beta * lambda; # Alaphas N x 1
        cov_alpha = (1.0 ./ T) * (cov_resid - beta * (beta' * (cov_resid \ beta) \ beta')) .* (1.0 .+ lambda' * (cov_f \ lambda));
        test_statistics = T * (1.0 .+ lambda' * (cov_f \ lambda)) * alpha' * (cov_resid \ alpha);
        pvalue = 1.0 .- cdf.(Chisq(N - K), test_statistics);

    elseif test_method == "GMM" || test_method == "GLS GMM" # OLS GMM and GLS GMM
        lambda = beta' * beta \ (beta' * Ere); # OLS second stage regressions - Risk premiums K x 1
        alpha = Ere- beta * lambda; # Alaphas N x 1
        h = [ones(T, 1) factors];
        f = kron(ones(1, K+1), resid) .* kron(h, ones(1, N)); # Time-series moment conditions generalized to k factors
        f = [f re .- kron(ones(T, 1), (beta * lambda)')]; # Adding cross-sectional moment conditions gernerlized to k factors - p.225 in Cochrane (2000)
        d = - kron(h' * h / T, I(N)); # Jacobian of the time-seires GMM moment conditions
        d = [d zeros((1 + K) * N, K); kron([0 lambda'],I(N)) -beta]; # Adding the Jacobian of the cross-sectional GMM moment conditions

        # Weigthing matrix
        if test_method == "GMM"
            a = [I zeros((K + 1) * N, N); zeros(K, (K + 1) * N) beta']; # weighting matrix for GMM
        else
            a = [I zeros((K + 1) * N, N); zeros(K, (K + 1) * N) beta' * (cov_resid \ I)]; # weighting matrix for GLS GMM
        end

        invad = (a * d) \ I;
        premom = I(size(d * invad * a)[1]) - d * invad * a; #  for the sampling distributionof the moments - p.190 eqn (147) in Cochrane (2000)
        Ef2 = CovMoment(f, cov_method, lag); # spectral density matrix
        varb2 = (1.0 ./ T) * invad * a * Ef2 * a' * invad; # variance of the coefficients - p.190 eqn (146) in Cochrane (2000)
        varmom = (1.0 ./ T) * premom * Ef2 * premom'; # asymptotic covariance matrix of the moments - p.190 eqn (147) in Cochrane (2000)
       
        test_statistics = alpha' * pinv(varmom[N * (K + 1) + 1:end, N * (K + 1) + 1:end]) * alpha; # test the moment conditions that the the error terms(alphas) in the second-stage cross-sectional regressions are jointly zero
        pvalue = 1.0 .- cdf.(Chisq(N - K), test_statistics);
    end

    return(Ere, lambda, alpha, test_statistics, pvalue)
end
