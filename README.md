# Adaptive-LASSO
Subset ARMA order selection via the adaptive Lasso

Model selection is critical to subset autoregressive moving-average (ARMA) modeling. This is commonly
done by subset selection methods, which may be computationally intensive and even impractical when the true ARMA orders of the underlying model are high. Conversely, automatic variable selection methods based on regularization do not directly apply to this problem because the innovation process is latent. To solve this problem, we
propose to identify the optimal subset ARMA model by fitting an adaptive Lasso regression of the time series on its lags and the lags of the residuals from a long autoregression
fitted to the time series data, where the residuals serve as proxies for the innovations
