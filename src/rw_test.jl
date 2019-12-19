"""
rw_test(p; method = "VR", m = 1)

TEST OF RANDOM WALK

Random walk with uncorrelated increments (RW3) is the weakest form version of the random walk.
The other versions are random walk with iid increments (RW1) and random walk with independent increments (RW2).


Box and Pierce (1970) test of RW(3): the autocorrelation of increments at all leads and lags should be zero under random walk.
By summing the squared autocorrelation, the Box-Pierce Q-statistics is designed to detect departures from zero autocorrelations in either directions and at all lags

Ljung and Box (1978) test of RW(3): the Ljung-Box Q-statistics provides the finite-sample correction which yields a better fit to the chi-squared distribution for small sample size.

Variance ratio test of RW(1): the variance of random walk increments must be a linear function of the time interval under random walk. The variance ratios statistics tests such property.

INPUT
`p': a T x 1 vector of asset prices
`method': method used to test random walk hypthesis with uncorrelated increments:
        - BP: the Box and Pierce (1970) statistics
        - LB: the Ljung and Box (1978) statistics
        - VR: variance ratio test
`m': number of autocorrelation lags used for test

OUTPUT
`test_statistics': test statistics of the random walk test
`pvalue': p-value of the test statistics

"""

function rw_test(p, method = "VR", m = 1)
    r = diff(log.(p), dims = 1)
    T = length(r); # length of the time series

    if method == "BP" # Box and Pierce (1970) test of RW(3)
        Qm = 0.0;
        for k = 1:m
            rho_k = autocor(r, [k]; demean = false);
            Qm = Qm .+ rho_k.^2;
        end
        test_statistics = T * Qm;
        pvalue = 1.0 .- cdf.(Chisq(m), test_statistics);
    elseif method == "LB" # Ljung and Box (1978) test of RW(3)
        Qm = 0.0;
        for k = 1:m
            rho_k = autocor(r, [k]; demean = false);
            Qm = Qm .+ rho_k.^2;
        end
        test_statistics = T * (T + 2) / (T - m) * Qm;
        pvalue = 1.0 .- cdf.(Chisq(m), test_statistics);
    elseif method == "VR" # Variance ratio test of RW(1)
        VR = 1.0;
        for k = 1:m-1
            rho_k = autocor(r, [k]; demean = false);
            VR = VR .+ 2 .* (1 .- k / m) .* rho_k;
        end
        test_statistics = sqrt(T * m) * (VR .- 1.0) ./ sqrt(2.0 .* (m - 1));
        pvalue = 1.0 .- cdf.(Normal(0.0, 1.0), test_statistics);
    end

    return (test_statistics, pvalue)
end
