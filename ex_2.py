from math import e

from scipy.stats import norm
from numpy import random
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression


def nCr(n,k):
    return math.factorial(n)/(math.factorial(k)*(math.factorial(n-k)))

def save_and_return(data,value,t,pos_y):
    data[ (t,pos_y) ] = value
    return value

#A
full_tree = {}
#fill ful
def american_bionomial_option(S,K,t,r,delta,sigma,h,call_or_put,pos_y):
   
    if t == 5: 
        temp = max(call_or_put*(S-K),0)
        full_tree[ (t,pos_y) ] = temp
        return temp

    #handle probability cases    
    u = e**((r-delta)*h+sigma*h**(1/2))
    d = e**((r-delta)*h-sigma*h**(1/2)) 
    p = (e**((r-delta)*h)-d)/(u-d)

    return save_and_return(full_tree,max( call_or_put*(S-K), e**((-r)*h)*american_bionomial_option(u*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y+h)*p+e**((-r)*h)*american_bionomial_option(d*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y-h)*(1-p) ),t,pos_y)

S = 110 #Current stock price
K = 100 #Option exercise price
T = 5 #Time to maturity (in years)
r = 0.05 #Annual interest rate
delta = 0.02 #Annual (continuous) dividend yield
sigma = .3 #Annualized volatility of stock
call_or_put = 1 # 1 = call, -1 = put
h = 1 #step size in years

#print(american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0))

#print(full_tree)

#B
def monte_carlo_simulation(S,T,delta,sigma,r,h,rounds):
    sims = []
    for i in range(rounds): 
        sim = []
        t = 0
        s_t = S
        while t <= T:
            alpha = 0.09
            sim.append(s_t)
            s_t = s_t*e**((alpha-delta - (sigma**(2))/2)*h+ sigma*h**(1/2)*random.normal(0,1) )
            
            t += h
        sims.append(sim)
    return sims





def sims_plotter(sims,T):
    ts = [i for i in range(T+1)]
    for sim in sims:

        plt.plot(ts,sim)
    plt.show()

def price_options(T,K,h,call_or_put):
    # stock sims = [[13,13.8,14,15,16],[12,13,15,14,10],...]
    # cash_flows = [[13,13.8,14,15,16],[12,13,15,14,10],...]
    rounds = 10_000
    stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,rounds)
    cash_flows = []
    stop_flag = []

    
    for i in range(len(stock_sims)):
        stop_flag.append([0,0,0,0,0])
        cash_flows.append([0,0,0,0,0])
 
    #unique case at t==T
    
    for i in range(len(stock_sims)):
        stock_sims[i] = stock_sims[i][1:]
    
   
    for i in range(len(stock_sims)):
        for x in range(T):
            cash_flows[i][4]=max(call_or_put*(stock_sims[i][4] - K), 0)
            if call_or_put*(stock_sims[i][4] - K)> 0:
                stop_flag[i][4] = 1
    
    #print(cash_flows[4])

    for t in range(T-1,0,-h):
        #print(t)
        regression_table_t = []
        

        for i in range(len(stock_sims)):
            #get the last values
            #print(len(cash_flows[t+1]))
            
            #if the option is in the money

            #if call_or_put*(stock_sims[i][t-1] - K) > 0:
                
            regression_table_t.append( (cash_flows[i][t]*e**(-r+delta)*h , stock_sims[i][t-1] ) ) 
        
        
        #print(regression_table_t)
        if True:
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
                if exercise_value!=0:
                    for x in range(T):
                        stop_flag[i][x]=0
                        cash_flows[i][x]=0
                    stop_flag[i][t-1]=1
                    count_ex += 1
                    cash_flows[i][t-1]=exercise_value
            else:
                count_co += 1
                cash_flows[i][t-1]=0
        #print(f"count ex: {count_ex}")
        #print(f"count co: {count_co}")

    #return average
    summ = 0
    for i in range(len(stock_sims)):
        #print(cash_flows[i])
        for t in range(T):
            if stop_flag[i][t] == 1:
                summ += cash_flows[i][t]*e**((r+delta)*(t+1))
    #print(stop_flag)
    
    #for E
    #plot_hist_MC(stop_flag)
    
    return summ / len(stock_sims)
    
def get_regression(regression_table):
    x = []
    y = []
    #print(len(regression_table))
    for tup in regression_table: 
        x.append(tup[1])
        #print(tup[0])
        y.append(tup[0])
    
    #plt.figure(1)
    #plt.scatter(x, y)
    
    x = np.array(x).reshape((-1,1)) 
    y = np.array(y)
    #print(x)
    #print(y)

    model = LinearRegression().fit(x,y)
    #print(model.coef_)
    #print(model.intercept_)
    #print(model.intercept_, model.coef_[0])
    #plt.plot(x, model.intercept_+ x*model.coef_[0], "r")
    #plt.show()
    return model.intercept_, model.coef_[0]
    

#C
def exercise_bounderies(full_tree,call_or_put):
    h_1 = 1/2
    
    pass 

def make_granular_function(excersice_boundaries,h,times):
    y_data = []*T*int(1/h)
    #[0,1,2,2.66,3.33,4.7]
    
    growth_lst = []
    for i in range(len(excersice_boundaries)-1):
        growth = (excersice_boundaries[i+1]-excersice_boundaries[i]) / ( (times[i+1]-times[i])/h)
        growth_lst.append(growth)
    #print(growth_lst)
    
     
    growth_extenscive = []
    
    #for i in range(len(growth_lst)-1):
    #   growth_lst[i] = growth_lst[i]/( (times[i+1]-times[i])/h )
    

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



    

"""
def stock_sim():
    new_exercise_boundry_delta = T/len(exercise_bounderies)
    new_boundaries = []
    for i in range(T):
        new_boundaries.append(i*new_exercise_boundry_delta)
    
    for sim in stock_sims:
        for i in range(len(stock_sim)):
            last_exercise = i/T
            newest_exercise = 0
            for bound in new_boundaries:
                if last_exercise < bound:
                    newest_exercise = bound
                    break
            boundary = exercise_bounderies[new_exercise]



    exercise_bounderies[i]
"""
#E

def plot_hist_MC(table):
    x = []

    for lst in table:
        for t in range(len(lst)):
            if (lst[t]==1):
                x.append(t+1)
    plt.figure(2)
    plt.hist(x, bins = [0.5,1.5,2.5,3.5,4.5,5.5])
    plt.show()


if __name__ == "__main__":
    last = []
    first = []
    #for stock_sim in stock_sims:
    #    last.append(stock_sim[-1])
    #    first.append(stock_sim[1])

    #plt.hist(last,bins=100)
    #plt.show()
    #sims_plotter(stock_sims,T)

    run_all = False
    h = 0.1
    plt.plot( [h*i for i in range(int(T/h)+2)], make_granular_function([283.00947164471603, 283.00947164471603, 249.13137450756523, 222.62313534145832, 195.9736802534565, 169.9459364846473, 147.3750010220282, 127.80176670011116],h,[0, 1.5, 2.25, 3.5, 4.25, 4.5, 4.75, 5.0]))
    
    plt.show()

    if run_all: 
        if(call_or_put==1):
            print("call")
        else:
            print("put")
        print("MCTS",price_options(T,K,h,call_or_put))
        print("Tree",american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0))

        call_or_put = -1
        if(call_or_put==1):
            print("call")
        else:
            print("put")
        print("MCTS",price_options(T,K,h,call_or_put))
        print("Tree",american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0))

    


