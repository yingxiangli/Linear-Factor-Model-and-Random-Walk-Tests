"""
CovMoment(f; m = 0)

Calculate the spectural densitiy matrix - p.190 equation (145) in Cochrane (2000)

INPUT
`f`: T x q array of q moment conditions
`cov_method': method for calculating the variance-covariance matrix of the error terms:
            - time iid: error terms are iid across time but cross-sectionally correlated
            - iid: error terms are iiid both across time and assets
            - NW: Newey-West coveriance matrix with m lags
            - HH: Hansen-Hodrick covariance matrix with m lag
`m': number of lags for the Newey-West or Hansen-Hodrick covariance matrix
OUTPUT
- `Ef2`: q x q covariance matrix
"""

function CovMoment(f, cov_method = "time iid", m = 0)
    Ef2 = (1 / (T-1)) * (f' * f);
    if cov_method == "time iid"
        Ef2 = Ef2; # Spectural densitiy matrix under time iid errors
    elseif cov_method == "iid"
        Ef2 = diagm(diag(Ef2)); # Spectural densitiy matrix under iid errors
    elseif cov_method == "NW"
        for i = 1:m
            w = 1 - i / (m + 1);
            Ef2 = Ef2 + 2 * w * (1 / (T - 1 - i)) * f[1:end - i, :]' * f[1 + i:end, :]; # Spectural densitiy matrix under Newey-West errors
        end
    elseif cov_method == "HH"
        for i = 1:m
            Ef2 = Ef2 + (2 / (T - 1 - i)) * f[1:end - i, :]' * f[1 + i:end, :]; # Spectural densitiy matrix under Hansen-Hodrick errors
        end
    end

    return Ef2

end
