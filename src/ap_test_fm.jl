"""
ap_test_fm(re, factors)

CROSS-SECTIONAL TESTS USING THE FAMA-MACBETH REGRESSION

Fama-MacBeth standard errors do not include corrections for the fact that the betas are also estimated.

INPUTS
`re': T x N matrix of excess returns, where T is the number of periods for each test assets and N is the number of test assets
`factors': T x K matrix of factors. Factors must be excess returns for time-series tests

OUTPUTS
`lambda': T x K vector of risk premiums for the factors
`alpha': T x N vector of alphas for each test assets

`test_statistics': Wald-statistics of the asset pricing test using Fama-MacBeth procedure
`pvalue': pvalue of the test statistics. The test statistics has a Chi-square distribution.

"""

function ap_test_fm(re, factors)
    (T, K) = size(factors);
    N = size(re, 2);
    cov_f = cov(factors); # coviance matrix for the factors K x K

    X = [ones(T,1) factors];
    (a, beta, resid, cov_resid) = OLSFn(re, X) # First-stage time-series regressions to obatin betas N x K

    lambda = fill(NaN, (T, K))
    alpha = fill(NaN, (T, N))
    
    for t = 1:T # run cross-sectional regressions for each period
        lambda[t, :] = ((beta' * beta) \  (beta' * vec(re[t, :])))'; # 1 x K
        alpha[t, :] = re[t, :]' - lambda[t, :]' * beta';
    end

    alpha_mean = mean(alpha, dims = 1)'; # time-series average of alphas N x 1
    alpha_mean_cov = zeros(N, N);
        
    for t = 1:T
        alpha_mean_cov = alpha_mean_cov + (alpha[t, :] - alpha_mean) * (alpha[t, :] - alpha_mean)';
    end
    alpha_mean_cov = alpha_mean_cov ./ T^2;
    
    test_statistics = alpha_mean' * (alpha_mean_cov \ alpha_mean);
    pvalue = 1.0 .- cdf.(Chisq(N - K), test_statistics);

    return (alpha_mean, test_statistics, pvalue)
end
