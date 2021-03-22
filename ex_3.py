import math 
import numpy as np
from numpy import random
import matplotlib.pyplot as plt





#3A

def probability_of_having_passed(t,S,sigma, S_t, H):
    if t == 0:
        return 0 
    return math.e**((-2/(t*sigma**2))*math.log(S/H)*math.log(S_t/H) )

def nCr(n,k):
    return math.factorial(n)/(math.factorial(k)*(math.factorial(n-k)))


def european_knock_in_option(S,K,T,r,delta,sigma,h,call_or_put,H):
    
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
        
        if price < H: 
            #Multiply by the chance of having passed H on the path, since the new valuation is E(pay|S): P*max(S-K,0) (1-P)*0 
            value = value*probability_of_having_passed(int(T/h),S,sigma,price,H)

        if value > 0:
            probabilty = nCr(len(possibilites),us) * p**us * (1-p)**(len(possibilites)-us)
            prices.append(value*probabilty)

    option_price = sum(prices)*math.e**(-(r)*T)
    return option_price

#Keytakaway: The code works, but we see no difference for option as there is no u and d combination after 5 years which gives a K < price < H and   

#print(f"European_knock_in_option: {european_knock_in_option(S,K,T,r,delta,sigma,call_or_put,H)}")


#3b

def european_monte_carlo(S,K,T,r,delta,sigma,call_or_put,H,h,number_of_simulations):
    #Monte Carlo with random distribution growt
    sims = []
    for i in range(number_of_simulations): 
        sim = []
        t = 0
        s_t = S
        while t <= T:
            alpha = 0.05
            sim.append(s_t)
            s_t = s_t*math.e**((alpha-delta - (sigma**(2))/2)*h+ sigma*h**(1/2)*random.normal(0,1) )
            t += h

        sims.append(sim)
    
    #Collect all exercercise values which have passed the barrier
    exercise_values = []

    for sim in sims:
        for price in sim: 
            if price >= H:
                #sim[-1] is S_T
                exercise_values.append( max(call_or_put*(sim[-1] - K), 0) )
                break                 
    
    
    #discount all cashflows
    cash_flow = 0
    for exercise_value in exercise_values:
        
        cash_flow += exercise_value*math.e**(-r*T)
    
    return cash_flow/number_of_simulations



def meta_simulator(simulation_round_limit,simulation_step_size,simulation_amount,S,K,T,r,delta,sigma,call_or_put,H,h):
    at_simulation = simulation_step_size
    variance = []
    x_axis = [at_simulation]
    
    while at_simulation <= simulation_round_limit:
        variance_i = []
        for i in range(simulation_amount):
            variance_i.append(european_monte_carlo(S,K,T,r,delta,sigma,call_or_put,H,h,at_simulation))
               
        variance.append(np.var(variance_i))

        at_simulation += simulation_step_size
        x_axis.append(at_simulation)
    return x_axis,variance


#3c
def save_and_return(data,value,t,pos_y):
    data[ (t,pos_y) ] = value
    return value
 


def american_bionomial_barrier_option(S,K,t,r,delta,sigma,h,call_or_put,pos_y,dictionary):
    if t == 5: 
        temp = max(call_or_put*(S-K),0)
        if S < 160:
            #print(probability_of_having_passed(t,110,sigma,S,160))
            temp = temp * probability_of_having_passed(t,110,sigma,S,160)
        
        dictionary[ (t,pos_y) ]  = temp
        return temp

    #handle probability cases    
    u = math.e**((r-delta)*h+sigma*h**(1/2))
    d = math.e**((r-delta)*h-sigma*h**(1/2)) 
    p = (math.e**((r-delta)*h)-d)/(u-d)

    return save_and_return(dictionary,max( call_or_put*(S-K)*probability_of_having_passed(t,110,sigma,S,H), math.e**((-r)*h)*american_bionomial_barrier_option(u*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y+h,dictionary)*p+math.e**((-r)*h)*american_bionomial_barrier_option(d*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y-h,dictionary)*(1-p) ),t,pos_y)



if __name__ == "__main__":
    #What simulations to run 
    #––––Run toggles––––
    three_a = 0
    three_b = 1 
    three_b_var = 0 
    three_c = 0
    #––––––––––––––––

    #Params for simulations
    S = 110 #Current stock price
    K = 100 #Option exercise price
    T = 5 #Time to maturity (in years)
    h = 0.05 #stepsize in years
    r = 0.05 #Annual interest rate
    delta = 0.02 #Annual (continuous) dividend yield
    sigma = .3 #Annualized volatility of stock
    call_or_put = 1 # 1 = call, -1 = put
    H = S+50 #barrier which the option has to pass in value
    number_of_simulations = 1000               



    #3a test
    if three_a: 
        print(f"European knock in:  {european_knock_in_option(S,K,T,r,delta,sigma,h,call_or_put,H)}")
       

    #3b test
    if three_b:
        h = 1 #NB running this for low h takes to long
        meta_simulations = 1000  
        
        
        data = []
        for x in range(number_of_simulations):
            data.append(european_monte_carlo(S,K,T,r,delta,sigma,call_or_put,H,h,number_of_simulations))
        print(f"European barrier Monte Carlo:  {np.mean(data)}")

    if three_b_var:
        simulation_round_limit = 1000
        simulation_step_size = 10
        simulation_amount = 10
        numb_simulations, variance = meta_simulator(simulation_round_limit,simulation_step_size,simulation_amount,S,K,T,r,delta,sigma,call_or_put,H,h)
        plt.plot(numb_simulations[1:],variance)
        plt.show()
    
    #3c test
    if three_c:
        tree = {}
        print(american_bionomial_barrier_option(S,K,0,r,delta,sigma,h,call_or_put,0,tree))
        print(tree)
    

