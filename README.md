# Adaptive-LASSO
Subset ARMA order selection via the adaptive Lasso

Model selection is critical to subset autoregressive moving-average (ARMA) modeling. This is commonly
done by subset selection methods, which may be computationally intensive and even impractical when the true ARMA orders of the underlying model are high. Conversely, automatic variable selection methods based on regularization do not directly apply to this problem because the innovation process is latent. To solve this problem, we
propose to identify the optimal subset ARMA model by fitting an adaptive Lasso regression of the time series on its lags and the lags of the residuals from a long autoregression
fitted to the time series data, where the residuals serve as proxies for the innovations.


In an ARMA(p,q) model with an $AR(\infty)$ representation, it implies that all coefficients ($\Phi$) are known. Consequently, the entire set of residuals ${\epsilon_t}$ can be derived from the equation $(1-\sum_{j=1}^{\infty}a_jB^j)Y_t=\epsilon_t$.
	
This suggests that for a long-term fitting of the AR(p) process to $Y_t$, one can approximate the unobserved $\epsilon_t$. However, the challenge lies in selecting appropriate model orders. The conventional method, still in use today, involves heuristically determining the model orders by inspecting sample auto-correlation and partial auto-correlation functions, such as ACF and PACF \cite{book3}. Despite its common usage, this traditional approach may result in misspecified models and potentially lead to poor forecasting performance.
	
When dealing with a time series showcasing a seasonal or cyclical phenomenon with a period of $s$, the initial step in choosing an appropriate model order involves defining $m = \max{(p,q)}$. In most cases, $m$ remains small relative to the total observations $T$ ($m < T$) and satisfies $m < s$. A situation where $m \geq s$ can lead to extended gaps in pertinent information for accurate predictions.\\
	
For time series exhibiting periodicity, the Seasonal Autoregressive Moving Average (SARMA) model serves as a natural extension of the ARMA framework. For cases where $s > 1$, the SARMA model takes the form $SARMA(p,q) \times (P,Q)$, where $P$ denotes the order of the seasonal autoregressive component, and $Q$ is the order of the seasonal moving average component. This model can be expressed through Equation (8) as:
	
$\Phi_P(B^s)\phi_p(B) Y_t=\Theta_Q(B^s)\theta_p(B)\epsilon_t$
	
	To simplify the SARMA model, it can be reduced to an $ARMA(p^*,q^*)$ process, where $\Phi_{p}^{}(B)Y_t=\Theta_{q}^{}(B)\epsilon_t$. Here, $m = \max\{p^*,q^*\} = \max\{Ps+p,Qs+q \}$, with $[p,P,q,Q,s] \in \mathbb{N}^5$. Consequently, our focus shifts to fitting the $ARMA(p^*,q^*)$ model for further analysis. \\
	Let p and q to be upper bounds such that $p\geq p^*$ and $q\geq q^*$, for the $(p+1)(q+1)$ different models, the final order selection can be based off minimization of some measure of prediction error. Commonly used metrics for this purpose include information criteria such as AIC or BIC, which effectively penalize model complexity. To facilitate this process, step-wise selection algorithms are typically employed. These algorithms systematically assess and compare various model configurations, iteratively refining the selection based on the chosen criterion. The problem here is that the model selection becomes less effective and very computational costly for large p \& q when searching for $2^{p+q}$ unique subset of ARMA models. \\
	Chen \& Chan \cite{chan} introduce order selection via adaptive Lasso. Let $y = [y_m,...,y_\tau]^{'} ,\epsilon=[\epsilon_m,...,\epsilon_\tau]^{'}$, $\beta = [\phi^{'},\theta^{'}]=[\phi_1,...,\phi_p,\theta_1,...,\theta_q]^{'}$ and
	\[
	X = \begin{bmatrix}
		x_{m}^{'} \\
		x_{m+1}^{'} \\
		\vdots \\
		x_{\tau}^{'}
	\end{bmatrix}
	= \begin{bmatrix}
		y_{m-1} & \cdots & y_{m-p} & \hat{\epsilon}_{m-1} & \cdots & \hat{\epsilon}_{m-q} \\
		y_{m} & \cdots & y_{m-p+1} & \hat{\epsilon}_{m} & \cdots & \hat{\epsilon}_{m-q+1} \\
		\vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\
		y_{\tau -1} & \cdots & y_{\tau-p} & \hat{\epsilon}_{\tau-1} & \cdots & \hat{\epsilon}_{\tau-q}
	\end{bmatrix}
	\]
	\begin{equation}
		y = X\beta + \epsilon
	\end{equation}
	where $\epsilon$, are the residuals from fitted $AR(p^{'})$ models used to estimate unknown innovations and $beta_j$ $j=1,...,p+q$ are nonzero . The goal is to identify the correct subset of nonzero components in the ARMA model. It has been proven that the adaptive Lasso method, in linear regression model can produce asymptotically unbiased estimators fro nonzero coefficients. However, this method does not apply directly to the ARMA model, due to $\epsilon_t$ terms.
	lasso estimator of $\beta$ is given by:
	\begin{equation}
		\hat{\beta}(\lambda) = \underset{\tau}{argmin}\left\{||y-X\beta|| + \lambda_\tau \sum_{j=1}^{p+q} \hat{w}_j |\tau_j| \right \}
	\end{equation}
	Where $\lambda_\tau$ is the tuning parameter controlling the degree of penalization and $\hat{w}$ are the weights. From \cite{chan} the weights are:
	$$\hat{w}=|\widetilde{\tau}|^{-\eta}$$
	$$\widetilde{\tau}=(\hat{X}^T \hat{X})^{-1}\hat{X}^{T}y$$
	which is the least square estimator of $\beta$ based on $\hat{X}$ and $\eta=2$ \cite{zou}.
