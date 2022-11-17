# epsilon=0.000001
# guess = [5.5e3, 5.5e3, 0.0]
#parameters, covariance = curve_fit(Hcsa1, x, y, p0 = guess, absolute_sigma = False, bounds=[(guess[0]-epsilon, guess[1] -epsilon,0-epsilon), (guess[0]+epsilon, guess[1] +epsilon, 0+epsilon)])
#parameters, covariance = curve_fit(Hcsa1, x, y, p0 = guess, sigma=1/y, absolute_sigma = True, bounds=[(-20000, 10,-1), (-1000, 100000,1)])
# parameters, covariance = curve_fit(Hcsa1, x, y)
# fit_A1 = parameters[0]
# fit_B1 = parameters[1]
# fit_C1 = parameters[2]

#parameters, covariance = curve_fit(HQ2, x, y )
#fit_A2 = parameters[0]
#fit_B2 = parameters[1]
#fit_C2 = parameters[2]
#fit_D2 = parameters[3]
#fit_E2 = parameters[4]

# print('A1 fit value', fit_A1); print('B1 fit value',fit_B1); print('C1 fit value',fit_C1)

# fit_Hcsa1 = Hcsa1(fit_A1, fit_B1, fit_C1, x)
# # Chi square Test
# print(stats.chisquare(f_exp= fit_Hcsa1, f_obs = y))

# #fit_HQ2 = HQ2(fit_A2, fit_B2, fit_C2, fit_D2, fit_E2, x)

# plt.plot(x, y, 'o', label='data')
# plt.plot(x, fit_Hcsa1, '-', label='fit')
# plt.xlabel ("Angle")
# plt.ylabel("Frequency")
# plt.show()
# plt.savefig("HQ_fit_plot")
#plt.plot(x,y, 'o')
#plt.plot(x,fit_HQ2,'-')
#plt.show()