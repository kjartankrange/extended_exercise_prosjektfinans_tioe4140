
import math 
from scipy.stats import norm
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def save_plot(plotter,name_of_figure):
    plotter.savefig(f"figures/{name}.png")

#A

def black_scholes_model(S,K,T,r,delta,sigma,call_or_put):
    d1 = (math.log(S/K)+(r-delta+(sigma**2)/2)*T) / (sigma*T**(1/2))
    d2 = d1-sigma*(T**(1/2))
    c_or_p = call_or_put*(S*math.e**(-delta*T)*norm.cdf(call_or_put*d1)-K*math.e**(-r*T)*norm.cdf(call_or_put*d2))
    return c_or_p
    

#B, mulig forskjell fra metoden over Sufficent number of steps

def nCr(n,k):
    return math.factorial(n)/(math.factorial(k)*(math.factorial(n-k)))

def european_bionomial_option(S,K,T,r,delta,sigma,h,call_or_put):

    #handle probability cases    
    d = math.e**((r-delta)*h-sigma*h**(1/2)) 
    u = math.e**((r-delta)*h+sigma*h**(1/2))
    p = (math.e**((r-delta)*h)-d)/(u-d)
    possibilites = []
    for us in range(int(T/h),0,-1):
        possibilites.append(us)  #["uuuuu","uuuud","uuudd","uuddd","udddd","ddddd"]

    prices = []
    
    for end_node in possibilites: 
        us = end_node
        price = S*u**(us)*d**(len(possibilites)-us)
        value = call_or_put*(price - K)
        if value > 0:
            probabilty = nCr(len(possibilites),us) * p**us * (1-p)**(len(possibilites)-us)
            prices.append(value*probabilty)
    option_price = sum(prices)*math.e**(-(r)*T)
    return option_price

#Function for showing how binomial price converges towards Black-Sholes with sufficient steps
def plot_binomial_prices(max_number_of_steps,delta_steps):
    at_steps = 5
    x = []
    y = []
    while at_steps <= max_number_of_steps:
        y.append(european_bionomial_option(S,K,T,r,delta,sigma,T/at_steps,call_or_put))
        x.append(at_steps)
        at_steps += delta_steps
    return x,y

#C GREEKS


def DELTA(s_primes,K,delta,sigma,r,T):
    deltas = []
    for s in s_primes:
        DELTA = math.e**(-delta*T)*norm.cdf((math.log(s/K)+(r-delta+(sigma**2)/2)*T) / (sigma*T**(1/2)))
        deltas.append(DELTA)
    return deltas




#TODO: Plot for different T's ? 
def GAMMA(deltas):
    gammas = []
    for i in range(len(deltas)-1):
        gammas.append( (deltas[i+1]-deltas[i]) )
    return gammas



def VEGA(s_primes):
    vegas = []
    for s in s_primes:
        vegas.append( black_scholes_model(s,K,T,r,delta,sigma+0.01,call_or_put) - black_scholes_model(s,K,T,r,delta,sigma,call_or_put))
    return vegas


def THETA(s_primes):
    theta = []
    for s in s_primes:
        theta.append( black_scholes_model(s,K,T-1/365,r,delta,sigma,call_or_put) - black_scholes_model(s,K,T,r,delta,sigma,call_or_put))

    return theta



def RHO(s_primes):
    rho = []
    for s in s_primes:
        rho.append( black_scholes_model(s,K,T,r + 0.01,delta,sigma,call_or_put) - black_scholes_model(s,K,T,r,delta,sigma,call_or_put))# Hvorfor er det bedre med en høyere risk free rate?

    return rho



def PSI(s_prime):
    psi = []
    for s in s_primes:
        psi.append( black_scholes_model(s,K,T,r,delta+0.01,sigma,call_or_put) - black_scholes_model(s,K,T,r,delta,sigma,call_or_put))# Hvorfor er det bedre med en høyere risk free rate?

    return psi




if __name__ == "__main__":
    #Parameters
    S = 110 #Current stock price
    K = 100 #Option exercise price
    T = 5 #Time to maturity (in years)
    r = 0.05 #Annual interest rate
    delta = 0.02 #Annual (continuous) dividend yield
    sigma = .3 #Annualized volatility of stock
    call_or_put = 1 # 1 = call, -1 = put
    h = 1
    
    #––––Run toggles––––
    one_a = 0
    one_b = 0
    one_c = 1
    one_c_graphs = 0 
    #plots for 1C)
    DELTA_plot = 1
    GAMMA_plot = 1
    VEGA_plot = 1
    THETA_plot = 1
    RHO_plot = 1
    PSI_plot = 1
    #––––––––––––––––
    if one_a:
        print("Black Scholes")
        print(black_scholes_model(S,K,T,r,delta,sigma,call_or_put))
    
    if one_b:
        print(f"\n European Binomial option")
        print(european_bionomial_option(S,K,T,r,delta,sigma,h,call_or_put))
        x,y = plot_binomial_prices(1000,10)
        plt.axhline(y=black_scholes_model(S,K,T,r,delta,sigma,call_or_put), color='r', linestyle='-')
        plt.plot(x,y)
        plt.xlabel("Number of steps")
        plt.ylabel("Price USD")
        plt.show()
    if one_c: 
        s_primes = [S]
        if one_c_graphs:
            s_primes = [x for x in range(1,401)]
        #print("\n GREEKS")
        deltas = DELTA(s_primes,K,delta,sigma,r,T)
        print(f"Delta {DELTA(s_primes,K,delta,sigma,r,T)[0]}")
        print(f"Gamma {GAMMA(DELTA(s_primes,K,delta,sigma,r,T))}")
        print(f"Vega {VEGA(s_primes)[0]}")
        print(f"Theta {THETA(s_primes)[0]}")
        print(f"Rho {RHO(s_primes)[0]}")
        print(f"Psi {PSI(s_primes)[0]}")
        if one_c_graphs:
            plt.plot(s_primes,deltas)
            plt.title("Delta")
            if DELTA_plot:
                plt.show()
            
            plt.plot(s_primes[1::],GAMMA(deltas))
            plt.title("Gamma")
            if GAMMA_plot: 
                plt.show()
            
            plt.plot(s_primes,VEGA(s_primes))
            plt.title("Vega")
            if VEGA_plot:
                plt.show()
            
            plt.plot(s_primes,THETA(s_primes)) #Kommentar når prisen er hly er det kjipt at den varer lengre fordi vi vil cashe ut
            plt.title("Theta")
            if THETA_plot:
                plt.show()
            
            plt.plot(s_primes,RHO(s_primes)) #Kommentar når prisen er hly er det kjipt at den varer lengre fordi vi vil cashe ut
            plt.title("Rho")
            if RHO_plot:
                plt.show()
            
            plt.plot(s_primes,PSI(s_primes))
            plt.title("Psi")
            if PSI_plot:
                plt.show()


