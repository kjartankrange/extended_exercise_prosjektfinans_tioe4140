
import math 
from scipy.stats import norm
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


#A

def black_scholes_model(S,K,T,r,delta,sigma,call_or_put):
    

    d1 = (math.log(S/K)+(r-delta+(sigma**2)/2)*T) / (sigma*T**(1/2))
    d2 = d1-sigma*(T**(1/2))
    
    c_or_p = call_or_put*(S*math.e**(-delta*T)*norm.cdf(call_or_put*d1)-K*math.e**(-r*T)*norm.cdf(call_or_put*d2))

    return c_or_p
    
    
S = 110 #Current stock price
K = 100 #Option exercise price
T = 5 #Time to maturity (in years)
r = 0.05 #Annual interest rate
delta = 0.02 #Annual (continuous) dividend yield
sigma = .3 #Annualized volatility of stock
call_or_put = -1 # 1 = call, -1 = put


print("Black Scholes")
print(black_scholes_model(S,K,T,r,delta,sigma,call_or_put))






#B, mulig forskjell fra metoden over Sufficent number of steps

def nCr(n,k):
    return math.factorial(n)/(math.factorial(k)*(math.factorial(n-k)))

def european_bionomial_option(S,K,T,r,delta,sigma,call_or_put):

    #handle probability cases    
    h = 1
    d = math.e**((r-delta)*h-sigma*h**(1/2)) 
    u = math.e**((r-delta)*h+sigma*h**(1/2))
    p = (math.e**((r-delta)*h)-d)/(u-d)
    possibilites = ["uuuuu","uuuud","uuudd","uuddd","udddd","ddddd"]

    prices = []
    
    for end_node in possibilites: 
        us = end_node.count("u")
        price = S*u**(us)*d**(len(possibilites)-1-us)
        value = call_or_put*(price - K)
        if value > 0:
            probabilty = nCr(len(possibilites)-1,5-us) * p**us * (1-p)**(len(possibilites)-1-us)
            prices.append(value*probabilty)
    option_price = sum(prices)*math.e**(-(r)*T)
    return option_price





print(f"\n European Binomial option")
print(european_bionomial_option(S,K,T,r,delta,sigma,call_or_put))






#C GREEKS
print("\n GREEKS")

def DELTA(s_primes,K,delta,sigma,r,T):
    deltas = []
    for s in s_primes:
        DELTA = math.e**(-delta*T)*norm.cdf((math.log(s/K)+(r-delta+(sigma**2)/2)*T) / (sigma*T**(1/2)))
        deltas.append(DELTA)
    return deltas
s_primes = [x for x in range(1,401)]
deltas = DELTA(s_primes,K,delta,sigma,r,T)

plt.plot(s_primes,deltas)
plt.title("Delta")
#plt.show()



#TODO: Plot for different T's ? 
def GAMMA(deltas):
    gammas = []
    for i in range(len(deltas)-1):
        gammas.append( (deltas[i+1]-deltas[i]) )
    return gammas

plt.plot(s_primes[1::],GAMMA(deltas))
plt.title("Gamma")
#plt.show()

def VEGA(s_primes):
    vegas = []
    for s in s_primes:
        vegas.append( black_scholes_model(s,K,T,r,delta,sigma+0.01,call_or_put) - black_scholes_model(s,K,T,r,delta,sigma,call_or_put))
    return vegas
plt.plot(s_primes,VEGA(s_primes))
plt.title("Vega")
#plt.show()

def THETA(s_primes):
    theta = []
    for s in s_primes:
        theta.append( black_scholes_model(s,K,T-1/365,r,delta,sigma,call_or_put) - black_scholes_model(s,K,T,r,delta,sigma,call_or_put))

    return theta

plt.plot(s_primes,THETA(s_primes)) #Kommentar når prisen er hly er det kjipt at den varer lengre fordi vi vil cashe ut
plt.title("Theta")
#plt.show()

def RHO(s_primes):
    rho = []
    for s in s_primes:
        rho.append( black_scholes_model(s,K,T,r + 0.01,delta,sigma,call_or_put) - black_scholes_model(s,K,T,r,delta,sigma,call_or_put))# Hvorfor er det bedre med en høyere risk free rate?

    return rho

plt.plot(s_primes,RHO(s_primes)) #Kommentar når prisen er hly er det kjipt at den varer lengre fordi vi vil cashe ut
plt.title("Rho")
#plt.show()

def PSI(s_prime):
    psi = []
    for s in s_primes:
        psi.append( black_scholes_model(s,K,T,r,delta+0.01,sigma,call_or_put) - black_scholes_model(s,K,T,r,delta,sigma,call_or_put))# Hvorfor er det bedre med en høyere risk free rate?

    return psi

plt.plot(s_primes,PSI(s_primes))
plt.title("Psi")
#plt.show()