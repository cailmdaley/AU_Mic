##emcee
# def lnlike(theta, x, y, yerr):
#     a, b = theta
#     model = a*x**2 + b
#     inv_sigma2 = 1.0/(yerr**2)
#     return -0.5*(np.sum((y-model)**2*inv_sigma2))
# 
# nll = lambda *args: -lnlike(*args)
# result = op.minimize(nll, [2.5, 0.75], args=(x, data, gaussian)) #result["x"]
# 
# def lnprior(theta):
#     a, b = theta
#     if -10. < a < 10. and -10. < b < 10.0:
#         return 0.0
#     return -np.inf
# 
# def lnprob(theta, x, y, yerr):
#     lp = lnprior(theta)
#     if not np.isfinite(lp):
#         return -np.inf
#     return lnlike(theta, x, y, yerr) + lp
# 
# 
# ndim, nwalkers = 2, 100
# #pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
# pos = [[1.5, 5.] + 5*np.random.randn(ndim) for i in range(nwalkers)]
# 
# sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, data, 1.6))
# sampler.run_mcmc(pos, 500)
# 
# 
# #burn in plot
# plt.close()
# fig, axes = plt.subplots(2, 1, sharex=True, figsize=(8, 9))
# axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
# axes[0].yaxis.set_major_locator(MaxNLocator(5))
# axes[0].axhline(3, color="#888888", lw=2)
# axes[0].set_ylabel("$a$")
# 
# axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
# axes[1].yaxis.set_major_locator(MaxNLocator(5))
# axes[1].axhline(1, color="#888888", lw=2)
# axes[1].set_ylabel("$b$")
# axes[1].set_xlabel("step number")
# 
# fig.tight_layout(h_pad=0.0)
# plt.show()
# 
# 
# #corner plot
# samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
# fig = corner.corner(samples, labels=["$a$", "$b$"],truths=[3, 1, np.log(gaussian)])
# plt.show()
