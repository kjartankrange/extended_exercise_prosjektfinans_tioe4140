from math import e

from scipy.stats import norm
from numpy import random
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

S = 110 #Current stock price
K = 100 #Option exercise price
T = 5 #Time to maturity (in years)
r = 0.05 #Annual interest rate
delta = 0.02 #Annual (continuous) dividend yield
sigma = .3 #Annualized volatility of stock
call_or_put = -1 # 1 = call, -1 = put
h = 1/4 #step size in years



def nCr(n,k):
    return math.factorial(n)/(math.factorial(k)*(math.factorial(n-k)))

stock_tree = {}
def save_and_return(S, data,value,t,pos_y):
    data[ (t,pos_y) ] = value
    stock_tree[(t, pos_y)] = S
    return value

#A
full_tree = {}

#fill ful
def american_bionomial_option(S,K,t,r,delta,sigma,h,call_or_put,pos_y):
   
    if t == 5: 
        temp = max(call_or_put*(S-K),0)
        full_tree[ (t,pos_y) ] = temp
        stock_tree[ (t,pos_y) ] = S
        return temp

    #handle probability cases    
    u = e**((r-delta)*h+sigma*h**(1/2))
    d = e**((r-delta)*h-sigma*h**(1/2)) 
    p = (e**((r-delta)*h)-d)/(u-d)
    return save_and_return(S, full_tree,max( call_or_put*(S-K), e**((-r)*h)*american_bionomial_option(u*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y+h)*p+e**((-r)*h)*american_bionomial_option(d*S,K,t+h,r,delta,sigma,h,call_or_put,pos_y-h)*(1-p) ),t,pos_y)





def exercise_boundary(full_tree, call_or_put,h):
    stock_price = []
    time = []
    y_pos = 0
    max_y = 0
    boundary = call_or_put* 1000
    for x_pos in range(T*(int(1/h))+1):
        y_pos = -max_y
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

#print(american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0))

#print(full_tree)

#print(stock_tree)

#exercise_boundary(full_tree, call_or_put,h)






#B
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

rounds = 10000
stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,rounds)

def sims_plotter(sims,T):
    ts = [i for i in range(T+1)]
    for sim in sims:

        plt.plot(ts,sim)
    plt.show()


def price_options(T,K,h,call_or_put):
    # stock sims = [[13,13.8,14,15,16],[12,13,15,14,10],...]
    # cash_flows = [[13,13.8,14,15,16],[12,13,15,14,10],...]

    cash_flows = []
    stop_flag = []
    extra=int(1/h)
    
    for i in range(len(stock_sims)):
        stop_flag.append([0]*T*extra)
        cash_flows.append([0]*T*extra)
 
    #unique case at t==T
    
    for i in range(len(stock_sims)):
        stock_sims[i] = stock_sims[i][1:]
    
   
    for i in range(len(stock_sims)):
        cash_flows[i][-1]=max(call_or_put*(stock_sims[i][-1] - K), 0)
        if call_or_put*(stock_sims[i][-1] - K)> 0:
            stop_flag[i][-1] = 1
    
    #print(cash_flows[4])

    for t in range(T*extra-1,0,-1):
        #print(t)
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
                if exercise_value!=0:
                    for x in range(T*extra):
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
    

#C
def exercise_bounderies(full_tree,call_or_put):
    h_1 = 1/2
    
    pass 

#E

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



print(american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0))
time, ex_bounderies = exercise_boundary(full_tree, call_or_put,h)


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

    plt.figure(3)
    plt.hist(ex_times, bins = 5) #[0.5,1.5,2.5,3.5,4.5,5.5]
    plt.show()

stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,10000)
ex_times(stock_sims, time, ex_bounderies)


#F
def plot_hist_MC(table):
    x = []
    for lst in table:
        for t in range(len(lst)):
            if (lst[t]==1):
                x.append(t*h+h)
    plt.figure(2)
    plt.hist(x, bins = [0.5,1.5,2.5,3.5,4.5,5.5]) #[0.5,1.5,2.5,3.5,4.5,5.5]
    plt.show()

last = []
first = []

#for stock_sim in stock_sims:
#    last.append(stock_sim[-1])
#    first.append(stock_sim[1])

#plt.hist(last,bins=100)
#plt.show()
#sims_plotter(stock_sims,T)

if __name__ == "__main__":

    #2c
    task_2c = False
    if task_2c:
        call_or_put=1
        print(american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0))
        x,y= exercise_boundary(full_tree, call_or_put,h)
        call_or_put=-1
        print(american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0))
        x,y=exercise_boundary(full_tree, call_or_put,h)

    #print(full_tree)

    #print(stock_tree)


    if(call_or_put==1):
        print("call")
    else:
        print("put")
    #print("MCTS",price_options(T,K,h,call_or_put))
    #print("Tree",american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0))
    #print(full_tree)

    call_or_put = -1
    if(call_or_put==1):
        print("call")
    else:
        print("put")
    print("MCTS",price_options(T,K,h,call_or_put))
    #print("Tree",american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0))


