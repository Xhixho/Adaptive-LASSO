# Subset ARMA order selection via the adaptive Lasso


Model selection is critical to subset autoregressive moving-average (ARMA) modeling. This is commonly
done by subset selection methods, which may be computationally intensive and even impractical when the true ARMA orders of the underlying model are high. Conversely, automatic variable selection methods based on regularization do not directly apply to this problem because the innovation process is latent. To solve this problem, we
propose to identify the optimal subset ARMA model by fitting an adaptive Lasso regression of the time series on its lags and the lags of the residuals from a long autoregression
fitted to the time series data, where the residuals serve as proxies for the innovations.


For more details, you can read the full paper By Chen \& Chan 
https://www.intlpress.com/site/pub/files/_fulltext/journals/sii/2011/0004/0002/SII-2011-0004-0002-a014.pdf


Here you will find the Code about the method proposed from the paper and the Monte Carlo MC  simulations developed by Eleni Dretaki \& George Xhixho
There are also some plug ins, we add some extra features like Adaptive Ridge and Adaptive Elastic-Net.
