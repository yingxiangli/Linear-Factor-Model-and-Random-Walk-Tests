"""
OLSFn

Run the (first-stage) time-series regressions for the asset pricing model to obtain the alpha and betas for each test assets

INPUT
`Y': T x N matrix of dependent varaible
`X': T x (K + 1) matrix with K explanatory variables and one intercept term
`cov_method': method for calculating the variance-covariance matrix of the error terms:
            - time iid: error terms are iid across time but cross-sectionally correlated
            - iid: error terms are iiid both across time and assets
            - NW: Newey-West covariance matrix with m lags
            - HH: Hansen-Hodrick covariane matrix with m lags
`m': number of lags for the Newey-West or Hansen-Hodrick covariance matrix


OUTPUT
`alpha': N x 1 vector of alphas of the test assets
`beta': N x K matrix of the betas of the test assets. There are N test assets and K factors
`resid': T x N matrix of the pricing errors
`cov_resid': N x N matrix of covariance matrix of the pricing errors
"""


function OLSFn(Y, X, cov_method = "time iid", m = 0)
    beta = (X' * X) \ (X' * Y); # (K+1) x N
    resid = Y - X * beta;
    alpha = beta[1, :]; # alphas N x 1
    beta = beta[2:end, :]'; # betas N x K

    # Method of estimating the covariance matrix of residuals for OLS
    if cov_method == "time iid" # iid across time but cross-sectionally correlated
        cov_resid = cov(resid);
    elseif cov_method == "iid" # iid
        cov_resid = cov(resid);
        cov_resid = diagm(diag(cov_resid));
    elseif cov_method == "NW" # Newey-West
        cov_resid = cov(resid);
        for i = 1:m
            w = 1 - i / (m + 1);
            cov_resid = cov_resid + 2 * w * (1 / (T - 1 - i)) * resid[1:end - i, :]' * resid[1 + i:end, :];
        end
    elseif cov_method == "HH" # Hansen-Hodrick
        cov_resid = cov(resid)
        for i = 1:m
            cov_resid = cov_resid + (2 / (T - 1 - i)) * resid[1:end - i, :]' * resid[1 + i:end, :];
        end
    end

    return (alpha, beta, resid, cov_resid)
end
