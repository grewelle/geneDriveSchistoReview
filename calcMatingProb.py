import scipy #import scipy packages and math
import math
import numpy as np
from scipy.integrate import quad
import pylab as pl


def integrand(theta, alpha, k):
    return (1-np.cos(theta))/(1+alpha*np.cos(theta))**(1+k)


def main():
    steps = 1
    years =10
    eggshedIndex = []
    numInfectedSnails = []
    comp = 1 #number of body compartments
    m = 40/comp #mean worm burden
    k = 0.24 #clumping parameter (also known as r in negative binomial)
    phi = [] #probability of mating
    w = [] #mean shedding female worms
    sigma = [] #prevalence of at least one mated pair (egg-shedding)
    #b = 9.137*10**-7/1.74622164 #200 eggs per worm per day * 7 days per week * prob of egg reaching water and contacting snail
    b = 8.302489585*10**-6
    #b = 4.41176*10**-5/18.342036000428315

    #b = 3.807106599*10**-6 / 1.94378511811
    #b = 9.517766497 * 10 ** -8 / 3.28311420544
    a1 = 1.066666666666666667*10**-3/comp #50 cercariae shed per day per snail * 7 days per week * prob of cerc contacting human and maturing
    N = 10**4 #number of total snails
    mu_1 = .004 #+ .02666667 #weekly death rate of worms per capita
    mu_2 = 0.25 #weekly death rate of infected snails per capita
    H = 1000*comp #number of humans
    rho1 = 0 #fraction of resistant snails [0-1]
    g = 0./12 #genedrive efficiency divided by weeks in snail generation

    """for i in range(steps):
        alpha = m / (m + k)
        phi.append(0)
        w.append(0)
        part1 = ((1-alpha)**(1+k))/(2*np.pi)
        part2 = quad(integrand, 0, 2*np.pi, args=(alpha,k))
        phi[i] = 1 - part1*part2[0]
        w[i] = 0.5*phi[i]*m
        sigma.append(0)
        sigma[i] = 1-2*(1+m/(2*k))**-k +(1+m/k)**-k
        beta = b * w[i]  # transmission parameter * number of mated worm pairs
        m += 1
    print(phi)
    print(w)
    print(sigma)"""

    def dX_dt(X, t=0):  # specify the initial time point t0=0

        alpha = X[0] / (X[0] + k)
        part1 = ((1 - alpha) ** (1 + k)) / (2 * np.pi)
        part2 = quad(integrand, 0, 2 * np.pi, args=(alpha, k))
        phi = 1 - part1 * part2[0]
        w = 0.5*phi
        #print(w*X[0])
        a = a1
        #a = a1*40/X[0] #include immunity proportional to intensity of worm burden
        #X[2] = rho1*X[0]/40 #include evolving resistance in snails proportional to intensity of worm burden

        #print(w)
        """a = a*np.sin(np.pi/26*t+np.pi/2)+a
        b = (3.807106599*10**-6 / 18.342036000428315) * np.sin(
            np.pi / 26 * t + np.pi / 2) + 3.807106599*10**-6 / 18.342036000428315
        N = (7000)*np.sin(np.pi/26*t+np.pi/2)+10000"""

        beta= b*w
        #beta=b*X[0]
        #beta= b*1.94378511811
        #beta = 9.517766497*10**-8
        #print(beta)
        y = scipy.array([a*N*X[1] - mu_1 * X[0], beta * H * X[0]*(1-X[1]-X[2]) - mu_2 * X[1], g*X[2]*(1-X[2])])
        return y

    if __name__ == "__main__":

        X0 = scipy.array([m, 0.015, rho1])  # initials conditions: x0=10  and y0=5
        S = scipy.array([X0])
        t = scipy.linspace(0, 52*years, 52*years)
        #X, infodict = scipy.integrate.odeint(dX_dt, X0, t, full_output=True)
        """for z in range(years):
            t_i = scipy.linspace(0, 52 +1, 52 +1)  # create 1000 time points stating at t0=0

            X, infodict = scipy.integrate.odeint(dX_dt, X0, t_i, full_output=True)

            S = np.append(S, X[1:], axis=0)
            #print(S)
            X0 = scipy.array([.25*X[-1,0], X[-1,1], X[-1,2]])
        S=S[0:-1]"""
        for u in range(101):
            X, infodict = scipy.integrate.odeint(dX_dt, X0, t, full_output=True)
            #print(100*(40-X[-1, 0])/40) # percent reduction in mean worm burden
            print(X[0,0])
            alpha_ = X[-1, 0] / (X[-1, 0] + k)
            part1_ = ((1 - alpha_) ** (1 + k)) / (2 * np.pi)
            part2_ = quad(integrand, 0, 2 * np.pi, args=(alpha_, k))
            phi_ = 1 - part1_ * part2_[0]
            w = 0.5 * phi_
            #print((18.3420359942-w)/18.3420359942*100) #percent reduction in mated worm pairs
            sigma = 1 - 2 * (1 + X[-1, 0] / (2 * k)) ** -k + (1 + X[-1, 0] / k) ** -k
            #print(100*(0.602598259783-sigma)/0.602598259783) #percent reduction in prevalence
            r0 = H*N*a1*b*w*(1-rho1)/mu_1/mu_2
            #print(r0)
            #print(100*(1.015228426 - r0)/1.015228426) #percent reduction in R0
            rho1 += .01
            X0 = scipy.array([m, 0.015, rho1])
        #	infodict['message']                                   #'Integration successful.'

        #pl.plot(t, X[:, 0], 'r-', t, X[:, 1], 'b-')  # x = X[:,0] and y = X[:,1]
        #pl.plot(t, 1 - 2 * (1 + S[:,0] / (2 * k1)) ** -k1 + (1 + S[:,0] / k1) ** -k1)
        #pl.plot(t, 1 - 2 * (1 + S[:, 0] / (2 * k1*S[t,0]/40)) ** -k1*S[t,0]/40 + (1 + S[:, 0] / k1*S[t,0]/40) ** -k1*S[t,0]/40)

        """for u in range(101):
            alpha_ = X[-1, 0] / (X[-1, 0] + k)
            part1_ = ((1 - alpha_) ** (1 + k)) / (2 * np.pi)
            part2_ = quad(integrand, 0, 2 * np.pi, args=(alpha_, k))
            phi_ = 1 - part1_ * part2_[0]
            w = 0.5 * phi_ * X[-1, 0]
            for z in range(10):
                t_i = scipy.linspace(0, 520 / 10, 1040 / 10)  # create 1000 time points stating at t0=0

                X, infodict = scipy.integrate.odeint(dX_dt, X0, t_i, full_output=True)

                S = np.append(S, X, axis=0)

                X0 = scipy.array([.25 * X[-1, 0], X[-1, 1], X[-1, 2]])
            S = S[1:]
            rho += .01

            flag = True
            for h in t:
                if (1 - 2 * (1 + S[h,0] / (2 * k)) ** -k + (1 + S[h,0] / k) ** -k) < 0.01 and flag == True:
                    print(t[h])
                    flag = False"""
        """for r in range(len(S)):
            k = k1 * X[0] / 40
            alpha = S[r,0] / (S[r,0] + k)
            part1 = ((1 - alpha) ** (1 + k)) / (2 * np.pi)
            part2 = quad(integrand, 0, 2 * np.pi, args=(alpha, k))
            phi = 1 - part1 * part2[0]
            w = 0.5 * phi * S[r,0]
            eggshedIndex.append(0)
            eggshedIndex[r] = w
            numInfectedSnails.append(0)
            numInfectedSnails[r] = S[r,1]*((7000)*np.sin(np.pi/26*r+np.pi/2)+10000)"""


        #pl.plot(t, eggshedIndex)
        #pl.plot(t, S[:,0])
        #pl.plot(t,S[:,1])
        #pl.plot(t,numInfectedSnails)
        #pl.plot(t, X[:,0])
        #pl.plot(t, X[:, 1])
        #pl.plot(t, X[:, 2])
        #pl.legend(['m', 'y'])
        #pl.show()

main()