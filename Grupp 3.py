
#Modellering över våren 2020, startdatum 31 januari (första fallet) slutdatum 30 juni

#%%
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# describe the model
def deriv(y, t, N, beta, gamma, delta):
    S, E, I, R = y
    
    if t < 29+14:               #Mellan 31 januari och 14 mars
        beta = beta_start
        
    elif 29+14 < t < 29+17:     #Mellan 14 mars och 17 mars
        beta = beta_res
        
    
    elif 29+17 < t < 29+27:     #Mellan 17 mars och 27 mars
        beta = beta_skol
        
        
    elif 29+27 < t:             #Efter 27 mars
        beta = beta_50
        
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - delta * E 
    dIdt = delta * E - gamma * I
    dRdt = gamma * I      
        
        
        
    return dSdt, dEdt, dIdt, dRdt




# describe the parameters
#N =  1113               # population
N = 10000000             #Sveriges befolkning
gamma=1/7                #sju sjukdagar enligt https://www.folkhalsomyndigheten.se/smittskydd-beredskap/utbrott/aktuella-utbrott/covid-19/fragor-och-svar/  
beta_start = 3*gamma
beta_res = 2.8*gamma     #reserestriktioner
beta_skol = 1.5*gamma    #nedstängning av skolor
beta_50 = 0.6*gamma      #gräns på 50 personer  
delta = 1/5              #fem inkubationsdagar enligt https://www.folkhalsomyndigheten.se/smittskydd-beredskap/utbrott/aktuella-utbrott/covid-19/fragor-och-svar/      

# initial conditions: one infected, rest susceptible
S0 = N-70
E0 = 70
I0 = 0
R0 = 0



t = np.linspace(0, 152, 153) # Grid of time points (in days)
y0 = S0, E0, I0, R0 # Initial conditions vector


# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta))
S, E, I, R = ret.T #delar upp resultatet i fyra separata delar




def plotsir(t, S, E, I, R):
  f, ax = plt.subplots(1,1,figsize=(10,4))
  # ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')     #Visas inte för att man ska se skillnaderna i de andra graferna
  ax.plot(t, E, 'y', alpha=0.7, linewidth=2, label='Exposed')
  ax.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
  ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')

  ax.set_xlabel('Time (days)')

  ax.yaxis.set_tick_params(length=0)
  ax.xaxis.set_tick_params(length=0)
  ax.grid(b=True, which='major', c='w', lw=2, ls='-')
  legend = ax.legend()
  legend.get_frame().set_alpha(0.5)
  for spine in ('top', 'right', 'bottom', 'left'):
      ax.spines[spine].set_visible(False)
    
  # plt.savefig("Plot.png")  
  plt.show();
  
plotsir(t, S, E, I, R)
#%%

