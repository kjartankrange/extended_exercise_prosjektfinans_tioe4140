import math 
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from math import e
from scipy.stats import norm
import numpy as np
from sklearn.linear_model import LinearRegression





#3A 

def probability_of_having_passed(t,S,sigma, S_t, H):
    if t == 0:
        return 0 
    return math.e**( (-2/(t*sigma**2))*math.log(S/H)*math.log(S_t/H) )

def nCr(n,k):
    return math.factorial(n)/(math.factorial(k)*(math.factorial(n-k)))


def european_knock_in_option(S,K,T,r,delta,sigma,h,call_or_put,H,knock_in):
    
    #handle probability cases    
    d = math.e**((r-delta)*h-sigma*h**(1/2)) 
    u = math.e**((r-delta)*h+sigma*h**(1/2))
    p = (math.e**((r-delta)*h)-d)/(u-d)
    
    possibilites = []
    for us in range(int(T/h),0,-1):
        possibilites.append(us)  #["uuuuu","uuuud","uuudd","uuddd","udddd","ddddd"]
    values = []
    #loop to find values in end nodes
    for end_node in possibilites: 
        us = end_node
        price = S*u**(us)*d**(len(possibilites)-us)
        value = max(call_or_put*(price - K),0)
        
        if knock_in: 
            if price < H: 
                value = value*probability_of_having_passed(int(T/h),S,sigma,price,H)
                

        else: #knock_out
            if price > H: 
                value = 0             
            else: 
                value = value*(1-probability_of_having_passed(int(T/h),S,sigma,price,H))
    
        if value > 0:
            probabilty = nCr(len(possibilites),us) * p**us * (1-p)**(len(possibilites)-us)
            values.append(value*probabilty)

    option_price = sum(values)*math.e**(-(r)*T)
    return option_price

#Keytakaway: The code works, but we see no difference for option as there is no u and d combination after 5 years which gives a K < price < H and   

#print(f"European_knock_in_option: {european_knock_in_option(S,K,T,r,delta,sigma,call_or_put,H)}")


#3b

def european_monte_carlo(S,K,T,r,delta,sigma,call_or_put,H,h,number_of_simulations,knock_in):
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
        knocked_out = 0
        for price in sim: 
            if price >= H and knock_in:
                #sim[-1] is S_T
                exercise_values.append( max(call_or_put*(sim[-1] - K), 0) )
                break 
            if price >= H and not knock_in:
                knocked_out = True
                break
        #This is a knock-out option and it hasn't been knocked out   
        if not knocked_out and not knock_in:
            exercise_values.append( max(call_or_put*(sim[-1] - K), 0) )
    
    
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
            variance_i.append(european_monte_carlo(S,K,T,r,delta,sigma,call_or_put,H,h,at_simulation,knock_in))
               
        variance.append(np.var(variance_i))

        at_simulation += simulation_step_size
        x_axis.append(at_simulation)
    return x_axis,variance


#3c
def save_and_return(data,value,t,pos_y):
    data[ (t,pos_y) ] = value
    return value
 
def american_bionomial_barrier_option(S,K,t,r,delta,sigma,h,call_or_put,pos_y,knock_in,dictionary):
    if t == T: 
        temp = max(call_or_put*(S-K),0)
        dictionary[ (t,pos_y) ]  = temp
        return temp
    
    #handle probability cases    
    u = math.e**((r-delta)*h+sigma*h**(1/2))
    d = math.e**((r-delta)*h-sigma*h**(1/2)) 
    p = (math.e**((r-delta)*h)-d)/(u-d)
    if t == T-h and K < S and S < H: 
           p_u = min(probability_of_having_passed(t+h,110,sigma,S*u,H),1) #print(f"Probability of having passed @upwards: {probability_of_having_passed(t+h,110,sigma,S*u,H)}")
           p_d = min(probability_of_having_passed(t+h,110,sigma,S*d,H),1) #print(f"Probability of having passed @downwards {probability_of_having_passed(t+h,110,sigma,S*d,H)}") #This wrongly prints larger than 1
           if not knock_in:
               p_u = 1 - p_u
               p_d = 1 - p_d
           return save_and_return(dictionary,max( call_or_put*(S-K), p_u*math.e**((-r)*h)*american_bionomial_barrier_option(u*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y+h,knock_in,dictionary)*p+p_d*math.e**((-r)*h)*american_bionomial_barrier_option(d*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y-h,knock_in,dictionary)*(1-p) ),t,pos_y)
    
    return save_and_return(dictionary,max( call_or_put*(S-K), math.e**((-r)*h)*american_bionomial_barrier_option(u*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y+h,knock_in,dictionary)*p+math.e**((-r)*h)*american_bionomial_barrier_option(d*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y-h,knock_in,dictionary)*(1-p) ),t,pos_y)

def monte_carlo_simulation(S,T,delta,sigma,r,h,rounds):
    sims = []
    for i in range(rounds): 
        sim = []
        t = 0
        s_t = S
        while t <= T:
            alpha = 0.05
            sim.append(s_t)
            s_t = s_t*e**((alpha-delta - (sigma**(2))/2)*h + sigma*h**(1/2)*random.normal(0,1) )
            
            t += h
        sims.append(sim)
    return sims

def price_options(T,K,h,call_or_put):
    # stock sims = [[13,13.8,14,15,16],[12,13,15,14,10],...]
    # cash_flows = [[13,13.8,14,15,16],[12,13,15,14,10],...]
    rounds = 100000
    stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,rounds)

    cash_flows = []
    stop_flag = []
    extra=int(1/h)
    
    for i in range(len(stock_sims)):
        stop_flag.append([0]*T*extra)
        cash_flows.append([0]*T*extra)
 
    #unique case at t==T
    
    for i in range(len(stock_sims)):
        stock_sims[i] = stock_sims[i][1:]
    

    #print(cash_flows[4])

    #knock-in
    knock_in = []
    H = 160
    
    #price american knock-in and knock out

    if call_or_put==1:
        for i in range(len(stock_sims)):
            knock_in.append([0]*T*extra)
            flag = 0
            for t in range(len(stock_sims[i])):
                if stock_sims[i][t]>=H:
                    flag = 1
                knock_in[i][t]=flag
    else:
        for i in range(len(stock_sims)):
            knock_in.append([1]*T*extra)
            flag = 1
            for t in range(len(stock_sims[i])):
                if stock_sims[i][t]>=H:
                    flag = 0
                knock_in[i][t]=flag

    for i in range(len(stock_sims)):
        if sum(knock_in[i])>=1:
            cash_flows[i][-1]=max(call_or_put*(stock_sims[i][-1] - K), 0)
            if call_or_put*(stock_sims[i][-1] - K)> 0:
                stop_flag[i][-1] = 1

    for t in range(T*extra-1,0,-1):

        regression_table_t = []
        

        for i in range(len(stock_sims)):
            #get the last values
            #print(len(cash_flows[t+1]))
            
            #if the option is in the money

            if call_or_put*(stock_sims[i][t-1] - K) > 0:

                if cash_flows[i][t]!=0:

                    regression_table_t.append( (cash_flows[i][t]*e**(-r*h) , stock_sims[i][t-1] ) )
        
        #print(regression_table_t)
        if regression_table_t!=[]:
            c, c_x = get_regression(regression_table_t)
        #print(c)
        #print(c_x)
        count_ex = 0
        count_co = 0
        for i in range(len(stock_sims)):
            continuation_value = c+c_x*stock_sims[i][t-1]
            exercise_value = max(call_or_put*(stock_sims[i][t-1] - K),0)
            #print(f"{t}")
            #print(continuation_value)
            #print(exercise_value)
            #print("---")
            if exercise_value > continuation_value:
                if exercise_value!=0 and knock_in[i][t]==1:
                    for x in range(len(stop_flag[i])):
                        stop_flag[i][x]=0
                        cash_flows[i][x]=0
                    stop_flag[i][t-1]=1
                    cash_flows[i][t-1]=exercise_value
                    count_ex += 1
            else:
                count_co += 1
                cash_flows[i][t-1]=0
        #print(f"count ex: {count_ex}")
        #print(f"count co: {count_co}")

    #return average
    summ = 0
    for i in range(len(stock_sims)):
        #print(cash_flows[i])
        for t in range(0,T*extra):
            if stop_flag[i][t] == 1:
                summ += cash_flows[i][t]*e**((-r)*(t*h+h))
    #print(stop_flag)
    
    #for E
    plot_hist_MC(stop_flag)
    
    return summ / len(stock_sims)
    
def get_regression(regression_table):
    x = []
    y = []
    #print(len(regression_table))
    for tup in regression_table: 
        x.append(tup[1])
        #print(tup[0])
        y.append(tup[0])
    
    plt.figure(1)
    plt.scatter(x, y)
    
    x = np.array(x).reshape((-1,1)) 
    y = np.array(y)
    #print(x)
    #print(y)

    model = LinearRegression().fit(x,y)
    #print(model.coef_)
    #print(model.intercept_)
    #print(model.intercept_, model.coef_[0])
    plt.plot(x, model.intercept_+ x*model.coef_[0], "r")
    #plt.show()
    return model.intercept_, model.coef_[0]

def plot_hist_MC(table):
    x = []
    for lst in table:
        #print(lst)
        for t in range(len(lst)):
            if (lst[t]==1):
                x.append(t*h+h)
    plt.figure(2)
    plt.hist(x, bins = [0.5,1.5,2.5,3.5,4.5,5.5]) #[0.5,1.5,2.5,3.5,4.5,5.5]
    plt.show()


if __name__ == "__main__":
    #What simulations to run 
    #––––Run toggles––––
    three_a = 1
    three_b = 0 
    three_b_var = 0
    three_c = 0
    three_c_mc = 0
    #––––––––––––––––

    #Params for simulations
    S = 110 #Current stock price
    K = 100 #Option exercise price
    T = 5 #Time to maturity (in years)
    h = 1/4 #stepsize in years
    r = 0.05 #Annual interest rate
    delta = 0.02 #Annual (continuous) dividend yield
    sigma = .3 #Annualized volatility of stock
    call_or_put = 1 # 1 = call, -1 = put
    H = S+50 #barrier which the option has to pass in value
    number_of_simulations = 5_000               
    knock_in = 1 #1 for knock-in 0 for knock-out


    #3a test
    if three_a: 
        h = 0.01
        knock_in = 1
        print(f"European call knock-in:  {european_knock_in_option(S,K,T,r,delta,sigma,h,call_or_put,H,knock_in)}, H = {H}, step size = {h}" )
        knock_in = 0
        print(f"European call knock-out:  {european_knock_in_option(S,K,T,r,delta,sigma,h,call_or_put,H,knock_in)}, H = {H}, step size = {h}" )


    #3b test
    if three_b:
        h = 0.01 #NB running this for low h takes to long
        meta_simulations = 20  
        
        
        data_in = []
        data_out = []
        for x in range(meta_simulations):
            data_in.append(european_monte_carlo(S,K,T,r,delta,sigma,call_or_put,H,h,number_of_simulations,1))
            data_out.append(european_monte_carlo(S,K,T,r,delta,sigma,call_or_put,H,h,number_of_simulations,0))
        print(f"European knock-in Monte Carlo mean:  {np.mean(data_in)} variance: {np.var(data_in)}")
        print(f"European knock-out Monte Carlo mean:  {np.mean(data_out)} variance: {np.var(data_out)}")

    if three_b_var:
        simulation_round_limit = 10_000
        simulation_step_size = 10
        simulation_amount = 10
        plt.xlabel("Number of simulations")
        plt.ylabel("Variance")
        numb_simulations, variance = meta_simulator(simulation_round_limit,simulation_step_size,simulation_amount,S,K,T,r,delta,sigma,call_or_put,H,h)
        plt.plot(numb_simulations[1:],variance)
        
        plt.show()
    
    #3c test
    if three_c:
        delta = 0.0
        tree = {}
        knock_in = 1
        print(f"American call knock-in option: {american_bionomial_barrier_option(S,K,0,r,delta,sigma,h,1,0,knock_in,tree)}, h = {h}")
        knock_in = 0
        print(f"American call knock-out option: {american_bionomial_barrier_option(S,K,0,r,delta,sigma,h,1,0,knock_in,tree)}, h = {h}")
        knock_in = 1
        print(f"American put knock-in option: {american_bionomial_barrier_option(S,K,0,r,delta,sigma,h,-1,0,knock_in,tree)}, h = {h}")
        knock_in = 0
        print(f"American put knock-out option: {american_bionomial_barrier_option(S,K,0,r,delta,sigma,h,-1,0,knock_in,tree)}, h = {h}")
        if three_c_mc:
            if(call_or_put==1):
                print("call")
            else:
                print("put")
            print("MC",price_options(T,K,h,call_or_put))
            call_or_put = -1
            if(call_or_put==1):
                print("call")
            else:
                print("put")
            print("MC",price_options(T,K,h,call_or_put))
    

