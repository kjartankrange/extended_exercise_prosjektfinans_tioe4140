from math import e
from scipy.stats import norm
from numpy import random
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
import sys
sys.setrecursionlimit(10**6)



#2a
def nCr(n,k):
    return math.factorial(n)/(math.factorial(k)*(math.factorial(n-k)))

stock_tree = {}
def save_and_return(S, data,value,t,pos_y):
    data[ (t,pos_y) ] = value
    stock_tree[(t, pos_y)] = S
    return value

full_tree = {}

#fill full tree
def american_bionomial_option(S,K,t,r,delta,sigma,h,call_or_put,pos_y, T):
   
    if t == T: 
        temp = max(call_or_put*(S-K),0)
        full_tree[ (t,pos_y) ] = temp
        stock_tree[ (t,pos_y) ] = S
        return temp

    #handle probability cases    
    u = e**((r-delta)*h+sigma*h**(1/2))
    d = e**((r-delta)*h-sigma*h**(1/2)) 
    p = (e**((r-delta)*h)-d)/(u-d)
    return save_and_return(S, full_tree,max( call_or_put*(S-K), e**((-r)*h)*american_bionomial_option(u*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y+h,T)*p+e**((-r)*h)*american_bionomial_option(d*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y-h, T)*(1-p) ),t,pos_y)





#2b
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


def sims_plotter(sims,T):
    ts = [i for i in range(T+1)]
    for sim in sims:

        plt.plot(ts,sim)
    plt.show()


def price_options(T,K,h,call_or_put, stock_sims_in):
    # stock sims = [[13,13.8,14,15,16],[12,13,15,14,10],...]
    # cash_flows = [[13,13.8,14,15,16],[12,13,15,14,10],...]

    cash_flows = []
    stop_flag = []
    extra=int(1/h)
    
    for i in range(len(stock_sims_in)):
        stop_flag.append([0]*T*extra)
        cash_flows.append([0]*T*extra)
 
    #unique case at t==T
    
    for i in range(len(stock_sims_in)):
        stock_sims[i] = stock_sims_in[i][1:]
    
   
    for i in range(len(stock_sims)):
        cash_flows[i][-1]=max(call_or_put*(stock_sims[i][-1] - K), 0)
        if call_or_put*(stock_sims[i][-1] - K)> 0:
            stop_flag[i][-1] = 1
    

    for t in range(T*extra-1,0,-1):
        regression_table_t = []
        

        for i in range(len(stock_sims)):
            #get the last values
            
            #if the option is in the money

            if call_or_put*(stock_sims[i][t-1] - K) > 0:

                if cash_flows[i][t]!=0:

                    regression_table_t.append( (cash_flows[i][t]*e**(-r*h) , stock_sims[i][t-1] ) )
        
        if regression_table_t!=[]:
            c, c_x, c_x_2 = get_regression(regression_table_t)

        for i in range(len(stock_sims)):
            continuation_value = c+c_x*stock_sims[i][t-1]+c_x_2*stock_sims[i][t-1]**2
            exercise_value = max(call_or_put*(stock_sims[i][t-1] - K),0)

            if exercise_value > continuation_value:
                if exercise_value!=0:
                    for x in range(T*extra):
                        stop_flag[i][x]=0
                        cash_flows[i][x]=0
                    stop_flag[i][t-1]=1
                    cash_flows[i][t-1]=exercise_value
            else:
                cash_flows[i][t-1]=0


    #return average
    summ = 0
    for i in range(len(stock_sims)):
        for t in range(0,T*extra):
            if stop_flag[i][t] == 1:
                summ += cash_flows[i][t]*e**((-r)*(t*h+h))

    #for F
    if task_2f:
        plot_hist_MC(stop_flag)
    
    return summ / len(stock_sims)
    
def get_regression(regression_table):
    x = []
    y = []
    z = []
    for tup in regression_table: 
        x.append([tup[1], tup[1]**2])
        z.append(tup[1])
        y.append(tup[0])
    
    #plt.figure(1)
    #plt.scatter(z, y)
    
    x = np.array(x) #.reshape((-1,1)) 
    y = np.array(y)
    z = np.array(z).reshape((-1,1))

    model = LinearRegression().fit(x,y)
    """
    fx = []
    for i in range(len(z)):
        fx.append(model.intercept_+ z[i]*model.coef_[0]+model.coef_[1]*z[i]**2)

    plt.plot(z,fx, "r")
    plt.show()
    """
    return model.intercept_, model.coef_[0], model.coef_[1]
    
#2c

def starting_point(r, delta, sigma, K, call_or_put):
    h_1 = 1/2 - (r-delta)/sigma**2 + (((r-delta)/sigma**2-1/2)**2 + 2*r / sigma**2)**(1/2)
    h_2 = 1/2 - (r-delta)/sigma**2 - (((r-delta)/sigma**2-1/2)**2 + 2*r / sigma**2)**(1/2)
    if call_or_put == 1:
        H_c = K *(h_1/(h_1-1))
        return H_c
    else:
        H_p = K *(h_2/(h_2-1))
        return H_p


def exercise_boundary(full_tree, call_or_put,h):
    stock_price = []
    time = []
    y_pos = 0
    max_y = T*(int(1/h))
    start = starting_point(r, delta, sigma, K, call_or_put)
    z = american_bionomial_option(start,K, -T,r,delta,sigma,h,call_or_put,0, T)
    boundary = call_or_put * 1000
    for x_pos in range(T*(int(1/h))+1):
        y_pos =  - max_y
        max_y+=1
        while y_pos<=max_y:
            S = stock_tree[(x_pos*h,y_pos*h)]
            if (call_or_put==1 and S<boundary) or (call_or_put==-1 and S>boundary):
                if (call_or_put*(S-K)>e**((-r)*h)*full_tree[x_pos*h,y_pos*h]):
                    if boundary==call_or_put*1000:
                        stock_price.append(S)
                        time.append(0)
                    stock_price.append(S)
                    time.append(x_pos*h)
                    boundary = S
            y_pos+=2
    if task_2c:
        plt.figure(1)
        plt.scatter(time, stock_price, marker="|")
        plt.plot(time, stock_price)
        plt.xlim(0,T)
        if call_or_put==1:
            plt.ylim(100,500)
        else:
            plt.ylim(20,100)
        plt.show()
    return time, stock_price


#2e

def make_granular_function(excersice_boundaries,h,times):
    y_data = []*T*int(1/h)
    #[0,1,2,2.66,3.33,4.7]
    
    growth_lst = []
    for i in range(len(excersice_boundaries)-1):
        growth = (excersice_boundaries[i+1]-excersice_boundaries[i]) / ( (times[i+1]-times[i])/h)
        growth_lst.append(growth)
    
    growth_extenscive = []    

    t = 0
    while t <  T:
        #add growth in time period
        for i in range(len(times)-1): 
            if t < times[i+1]:
                growth_extenscive.append(growth_lst[i]) 
                break
        t += h
               
    
    value_t = excersice_boundaries[0]
    granular_exercise_lst = [value_t]
    
    for growth in growth_extenscive:
        value_t += growth
        granular_exercise_lst.append(value_t)

    return granular_exercise_lst


def ex_times(stock_sims, time, ex_bounderies):
    gran_ex_list = make_granular_function(ex_bounderies,h, time)
    ex_times = []
    for sim in stock_sims:
        for i in range(len(sim)): 
            if call_or_put==1:
                if sim[i] >= gran_ex_list[i]:
                    ex_times.append(i*h)
                    break
            else:
                if sim[i] <= gran_ex_list[i]:
                    ex_times.append(i*h)
                    break
    #percentage before T
    sum_total = 0
    sum_before_T = 0
    for time in ex_times:
        if time<T:
            sum_before_T+=1
        sum_total+=1
    print(f"Percentage of paths excercised before T = {T}: ",sum_before_T/sum_total*100)

    plt.figure(3)
    plt.hist(ex_times, bins = T*int(1/h)) #[0.5,1.5,2.5,3.5,4.5,5.5]
    plt.show()


#F
def plot_hist_MC(table):
    x = []
    for lst in table:
        for t in range(len(lst)):
            if (lst[t]==1):
                x.append(t*h+h)
    plt.figure(2)
    plt.hist(x, bins = T*int(1/h)) #[0.5,1.5,2.5,3.5,4.5,5.5]
    plt.show()



if __name__ == "__main__":
    #What simulations to run 
    #––––Run toggles––––
    task_2a = 0
    task_2b = 0
    task_2c = 1
    task_2d = 0
    task_2e = 0
    task_2f = 0
    #––––––––––––––––

    #Params for simulations
    S = 110 #Current stock price
    K = 100 #Option exercise price
    T = 5 #Time to maturity (in years)
    h = 1/2 #stepsize in years
    r = 0.05 #Annual interest rate
    delta = 0.02 #Annual (continuous) dividend yield
    sigma = .3 #Annualized volatility of stock
    call_or_put = 1 # 1 = call, -1 = put
    number_of_simulations = 20000               

    #2a

    if task_2a:
        call_or_put = 1
        print("American call option (binomial): ",american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0, T))
        call_or_put = -1
        print("American put option (binomial): ",american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0, T))

    #2b
    if task_2b:
        call_or_put = 1
        stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,number_of_simulations)
        print("American call option (MC): ",price_options(T,K,h,call_or_put, stock_sims))
        call_or_put = -1
        stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,number_of_simulations)
        print("American put option (MC): ",price_options(T,K,h,call_or_put, stock_sims))

    #2c
    #TODO: do we need to use the hint here? Add more st devs?
    if task_2c:

        stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,number_of_simulations)
        call_or_put=1
        x,y= exercise_boundary(full_tree, call_or_put,h)
        stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,number_of_simulations)
        call_or_put=-1
        x,y=exercise_boundary(full_tree, call_or_put,h)

    #2d - discussion

    #2e and 2f need same price paths
    if task_2e or task_2f:
        stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,number_of_simulations)

    #2e 
    if task_2e:
        call_or_put = -1
        price = american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0)
        time, ex_bounderies = exercise_boundary(full_tree, call_or_put,h)
        ex_times(stock_sims, time, ex_bounderies)
    
    #2f
    if task_2f:
        call_or_put = -1
        price  = price_options(T,K,h,call_or_put, stock_sims)
