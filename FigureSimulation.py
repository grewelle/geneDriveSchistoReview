import scipy
import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt
import seaborn as sns; sns.set(style="white", color_codes=True)

SMALL_SIZE = 36
MEDIUM_SIZE = 48
BIGGER_SIZE = 64

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#integration function for neg binomial mating
def integrand(theta, alpha, k):
    return (1-np.cos(theta))/(1+alpha*np.cos(theta))**(1+k)


def main():

    #parameter block
    years = 10
    m = 40  # equilibrium mean worm burden
    k = 0.24  # clumping parameter (also known as r in negative binomial)
    b = 2.075622396 * 10 ** -7
    a = 1.066666666666666667 * 10 ** -3   # 50 cercariae shed per day per snail * 7 days per week * prob of cerc contacting human and maturing
    N = 10 ** 4  # number of total snails
    mu_1 = .004  # + .02666667 #weekly death rate of worms per capita
    mu_2 = 0.25  # weekly death rate of infected snails per capita
    H = 1000  # number of humans
    rho = 0  # fraction of resistant snails [0-1]

    #begin simulation
    X0 = scipy.array([m, 0.015])  # initials conditions: x0=10  and y0=5
    t = scipy.linspace(0, 52 * years, 52 * years)

    def dX_dt(X, t=0):  # specify the initial time point t0=0

        alpha = X[0] / (X[0] + k)
        part1 = ((1 - alpha) ** (1 + k)) / (2 * np.pi)
        part2 = quad(integrand, 0, 2 * np.pi, args=(alpha, k))
        phi = 1 - part1 * part2[0]
        w = 0.5*phi
        beta= b*w
        #print(beta)
        y = scipy.array([a*N*X[1] - mu_1 * X[0], beta * H * X[0]*(1-X[1]-rho) - mu_2 * X[1]])
        return y


    wormBurdenRed = [] # per capita worm burden reduction %
    matedPairsRed = [] # per capita number mated pair reduction %
    prevalenceRed = [] # prevalence reduction %
    reproNumberRed = [] # effective reproductive number reduction %

    for u in range(101):
        X, infodict = scipy.integrate.odeint(dX_dt, X0, t, full_output=True)
        wormBurdenRed.append(100*(40-X[-1, 0])/40)
        alpha_ = X[-1, 0] / (X[-1, 0] + k)
        part1_ = ((1 - alpha_) ** (1 + k)) / (2 * np.pi)
        part2_ = quad(integrand, 0, 2 * np.pi, args=(alpha_, k))
        phi_ = 1 - part1_ * part2_[0]
        w = 0.5 * phi_
        matedPairsRed.append((18.3421134484-w*X[-1,0])/18.3421134484*100)
        sigma = 1 - 2 * (1 + X[-1, 0] / (2 * k)) ** -k + (1 + X[-1, 0] / k) ** -k
        prevalenceRed.append(100*(0.602598639541-sigma)/0.602598639541)
        rt = H * N * a * b * w * (1 - rho) / mu_1 / mu_2
        reproNumberRed.append(100*(1.0152306002 - rt)/1.0152306002)
        rho += .01
        X0 = scipy.array([m, 0.015, rho])

    rhoRange = scipy.linspace(0, 1, num=101, endpoint=True)
    print(len(rhoRange))


    fig = plt.figure(figsize=(18, 18))
    plt.axis((0, 1, 100, 0))
    plt.xlabel('Rho')
    plt.ylabel('Percent Reduction from Baseline')
    #wormBurden = plt.plot(rhoRange, wormBurdenRed,  'r:', label= 'Worm Burden', linewidth= 4)
    matedPairs = plt.plot(rhoRange, matedPairsRed, 'm--', label= 'Infection Intensity', linewidth= 4)
    prevalence = plt.plot(rhoRange, prevalenceRed, 'b-.', label= 'Prevalence', linewidth= 4)
    reproNumber = plt.plot(rhoRange, reproNumberRed, 'k', label= 'Reproductive No.', linewidth= 4)
    currentAxis = plt.gca()
    #plt.legend()

    fig.savefig("GeneDriveReviewFigure1b.png")
    fig.savefig("GeneDriveReviewFigure1b.tiff")

main()

