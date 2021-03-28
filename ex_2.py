from math import e
from scipy.stats import norm
from numpy import random
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
import sys
from numpy.polynomial import Polynomial

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
    ts = [i for i in range(T*int(1/h_sim)+1)]
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

            if  call_or_put*(stock_sims[i][t-1] - K) > 0:

                if cash_flows[i][t]!=0:

                    regression_table_t.append( (cash_flows[i][t]*e**(-r*h*(t-1)) , stock_sims[i][t-1] ) )
        
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
                summ += cash_flows[i][t]*e**((-r)*((t)*h))

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
    #c, stats = np.polynomial.polynomial.polyfit(x,y,2,full=True)

    """
    fx = []
    for i in range(len(z)):
        fx.append(model.intercept_+ z[i]*model.coef_[0]+model.coef_[1]*z[i]**2)

    plt.plot(z,fx, "r")
    plt.show()
    """
    return model.intercept_, model.coef_[0], model.coef_[1]
    #return c[0], c[1], c[2]
    
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
    max_y =  T*(int(1/h))
    start = starting_point(r, delta, sigma, K, call_or_put)
    z = american_bionomial_option(start,K, -T,r,delta,sigma,h,call_or_put,0, T)
    boundary = call_or_put * 1000
    for x_pos in range(T*(int(1/h))+1):
        y_pos =  - max_y*call_or_put
        max_y+=1
        while (y_pos<=max_y and call_or_put==1) or (y_pos>=-max_y and call_or_put==-1):
            S = stock_tree[(x_pos*h,y_pos*h)]
            if (call_or_put==1 and S<boundary) or (call_or_put==-1 and S>boundary):
                if (call_or_put*(S-K)>e**((-r)*h)*full_tree[x_pos*h,y_pos*h]):
                    if boundary==call_or_put*1000:
                        stock_price.append(S)
                        time.append(0)
                    stock_price.append(S)
                    time.append(x_pos*h)
                    boundary = S
                    break
            y_pos=y_pos+call_or_put*2
    stock_price.append(100)
    time.append(5)
    if task_2c:
        plt.figure(1)
        plt.scatter(time, stock_price, marker="|", s = 12*6)
        plt.plot(time, stock_price, linewidth=3.0, label = f"Sigma = {sigma}")
        plt.xlabel("Time in Years")
        plt.ylabel("Exercise Boundary")
        plt.xlim(0,T)
        plt.legend()
        if call_or_put==1:
            plt.ylim(100,500)
        else:
            plt.ylim(20,100)
        #plt.show()
    return time, stock_price



#2e

def make_granular_function(excersice_boundaries,h,times):
    y_data = []*T*int(1/h)
    #[0,1,2,2.66,3.33,4.7]
    
    growth_lst = []
    for i in range(len(excersice_boundaries)-1):
        if times[i+1]!=times[i]:
            growth = (excersice_boundaries[i+1]-excersice_boundaries[i]) / ( (times[i+1]-times[i])/h_sim)
            growth_lst.append(growth)
        else:
            growth_lst.append(0)
    
    growth_extenscive = []    

    t = 0
    while t <  T:
        #add growth in time period
        for i in range(len(times)-1): 
            if t < times[i+1]:
                growth_extenscive.append(growth_lst[i]) 
                break
        t += h_sim
               
    
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
                    ex_times.append(i*h_sim)
                    break
            else:
                if sim[i] <= gran_ex_list[i]:
                    ex_times.append(i*h_sim)
                    break
    #percentage before T
    sum_total = 0
    sum_before_T = 0
    sum_before_T_app = 0
    for time in ex_times:
        if time<T:
            sum_before_T+=1
        if time<4.5:
            sum_before_T_app+=1
        sum_total+=1
    print(f"Percentage of paths excercised before T  = {T} (Ex. boundary): ",sum_before_T/len(ex_times)*100)
    print(f"Percentage of paths excercised before T  = 4.5 (Ex. boundary): ",sum_before_T_app/sum_total*100)

    plt.figure(3)
    plt.title("Using Exercise Boundary")
    plt.xlabel("Exercise Time in Years")
    plt.ylabel("Number of Options exercised")
    plt.hist(ex_times, bins = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]) #[0,0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25]
    plt.show()



#F
def plot_hist_MC(table):
    x = []
    #percentage before T
    sum_total = 0
    sum_before_T = 0
    for lst in table:
        for t in range(len(lst)):
            if (lst[t]==1):
                x.append(t*h_sim+h_sim)
                if t!=(T-1):
                    sum_before_T+=1

        sum_total+=1
    print(f"Percentage of paths excercised before T  = {T} (MC)         : ",sum_before_T/sum_total*100)


    plt.figure(2)
    plt.title("Using the Least Squares Monte Carlo-method")
    plt.xlabel("Exercise Time in Years")
    plt.ylabel("Number of Options exercised")
    plt.hist(x, bins = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]) #[0.5,1.5,2.5,3.5,4.5,5.5], T*int(1/h)
    plt.show()

#TASK 2

if __name__ == "__main__":
    #What simulations to run 
    #––––Run toggles––––
    task_2a = 0
    task_2b = 0
    task_2c = 0
    task_2d = 0
    task_2e = 1
    task_2f = 1
    #––––––––––––––––

    #Params for simulations
    S = 110 #Current stock price
    K = 100 #Option exercise price
    T = 5 #Time to maturity (in years)
    h = 0.5 #stepsize in years
    r = 0.05 #Annual interest rate
    delta = 0.02 #Annual (continuous) dividend yield
    sigma = .3 #Annualized volatility of stock
    call_or_put = 1 # 1 = call, -1 = put
    number_of_simulations = 10000               

    #2a

    if task_2a:
        call_or_put = 1
        print("American call option (binomial): ",american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0, T))
        call_or_put = -1
        print("American put option (binomial): ",american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0, T))

    #2b
    if task_2b:
        puts = []
        calls = []
        for x in range(1):
            print(x)
            call_or_put = 1
            stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,number_of_simulations)
            price = price_options(T,K,h,call_or_put, stock_sims)
                
            calls.append(price)
            print("American call option (MC): ",price)
            
            call_or_put = -1
            stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h,number_of_simulations)
            price = price_options(T,K,h,call_or_put, stock_sims)
            puts.append(price)
            print("American put option (MC): ",price)
            
        print(calls)
        print(puts)
        variance_put = np.var(puts)
        variance_call = np.var(calls)
        mean_put = np.mean(puts)
        mean_call = np.mean(calls)
        print("varput", variance_put)
        print("varcall", variance_call)
        print("meancall", mean_call)
        print("meanput", mean_put)




    #2c
    if task_2c:

        call_or_put=1
        sigma = 0.1
        x,y= exercise_boundary(full_tree, call_or_put,h)
        sigma = 0.3
        x,y= exercise_boundary(full_tree, call_or_put,h)
        sigma = 0.5
        x,y= exercise_boundary(full_tree, call_or_put,h)
        plt.show()
        call_or_put=-1
        sigma = 0.1
        x,y= exercise_boundary(full_tree, call_or_put,h)
        sigma = 0.3
        x,y= exercise_boundary(full_tree, call_or_put,h)
        sigma = 0.5
        x,y= exercise_boundary(full_tree, call_or_put,h)
        plt.show()

    #2d - discussion

    #2e and 2f need same price paths
    h_sim=h
    if task_2e or task_2f:

        stock_sims  = monte_carlo_simulation(S,T,delta,sigma,r,h_sim,number_of_simulations)

    #2e 
    if task_2e:
        call_or_put = -1
        price = american_bionomial_option(S,K,0,r,delta,sigma,h,call_or_put,0, T)
        time, ex_bounderies = exercise_boundary(full_tree, call_or_put,h)
        ex_times(stock_sims, time, ex_bounderies)
    
    #2f
    if task_2f:
        call_or_put = -1
        price  = price_options(T,K,h_sim,call_or_put, stock_sims)
