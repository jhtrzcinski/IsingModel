"""
Jacob Trzcinski
Computational Physics
Dr. Matteo Luisi
Semester Project
"""
#%%
### Import some handy packages ###
import numpy as np
import matplotlib.pyplot as plt

### Make a definition of a hot start so we dont have to repeat it all the time ###
#   A hot start is one where it is an array with complete randomness of +1 and -1 states
def HotStart(N):
    system = np.random.choice([1,-1], size=(N,N))
    
    return system

k = 1.381e-23   # I'm going to define the Boltzman constant here so it's global

### make a definition to find the neighbors for a given cell in the array ###
def getNeighbors(arr, i, j):
    ### Boundary Conditions ###
    # boundary conditions to make topography of a torus for contiunity
    # clockwise starting from 12 (12,3,6,9)
    # n is a list of the 4 neighbors states
    # i >
    # j V

    N = np.shape(arr)[0]

    if i == 0:
        if j == 0:   # top left corner
            n = [arr[0,-1], arr[1,0], arr[0,1], arr[-1,0]]
        elif j == N-1: # bottom left corner
            n = [arr[0,-2], arr[1,-1], arr[0,0], arr[-1,-1]]
        else:   # everything else on the left edge
            n = [arr[0,j-1], arr[1,j], arr[0,j+1], arr[-1,j]]
    elif i == N-1:  
        if j == 0:  # top right corner
            n = [arr[-1,-1], arr[0,0], arr[-1,1], arr[-2,0]]
        elif j == N-1:   # bottom right corner
            n = [arr[-1,-2], arr[0,-1], arr[-1,0], arr[-2,-1]]
        else:   # everything else on the right edge
            n = [arr[-1,j-1], arr[0,j], arr[-1,j+1], arr[-2,j]]
    else: 
        if j == 0:  # top edge
            n = [arr[i,-1], arr[i+1,0], arr[i,1], arr[i-1,j]]
        elif j == N-1:   # bottom edge
                n = [arr[i,-2], arr[i+1,-1], arr[i,0], arr[i-1,-1]]
        else:   # everything inside
                n = [arr[i,j-1], arr[i+1,j], arr[i,j+1], arr[i-1,j]]

    return n

### definition of energy: E = -s_i,j * (sum of neighbor's states) ###
def getEnergy(arr, i, j):
    n = getNeighbors(arr, i, j)
    Energy = -1*arr[i,j] * sum(n)

    return Energy

def Ising(arr, T, n, N):
    for s in range(n):
        # Select random coordinates
        i = np.random.randint(0,N)
        j = np.random.randint(0,N)

        # find energy of those coordinates
        E = getEnergy(arr, i, j)

        # determine if the state will flip
        if E > 0:
            arr[i,j] = -1*arr[i,j]
        elif E <= 0:
            r = np.random.rand()
            if r < np.exp(2*E/(k*T)):
                arr[i,j] = -1*arr[i,j]
        ### else do nothing, no switch ###

    return arr

"""
Below is the routine to show the model go from hot start
all the way down to low T, showing the model at each temp increment
and the magnetization of the model as a function of temp
"""

T = 500   # initial temperature, no significant changes if it's any higher
M = []  # make a list to store the magnetization of the model
t = []  # make a list to store the temperatures used

IsingArray = HotStart(50)  # initialize the array with a hot start
# show that new array
plt.figure()
plt.imshow(IsingArray)
plt.title('Initialized Model, Hot Start')
plt.colorbar()
plt.show()

# loop to increment the temperature from T = initial to T = 0.5
while T > 0:
    t.append(T) # add T to list
    # run the ising method to flip states
    IsingArray = Ising(IsingArray, T, 5000, 50) 

    M.append(np.mean(IsingArray))   # add magnetization of model to list

    # show the pretty model and have it show each temp in the title
    plt.figure()
    plt.imshow(IsingArray)
    plt.title('Model at T = %.0f' %T)
    plt.colorbar()
    plt.show()
    if T > 25 and T <= 100:
        T -= 25
    elif T <= 25 and T > 10:
        T -= 10
    elif T <= 5 and T > 1:
        T -= 1
    elif T == 1:
        T -= 0.5
    else:
        T -= 200    # update temperature

# when the loop is done, go one step further to T = 0.1
T = 0.1
t.append(T) # add T to list
IsingArray = Ising(IsingArray, T, 5000, 50) #update model
M.append(np.mean(IsingArray))   # add magnetization to list
    
# show the updated array
plt.figure()
plt.imshow(IsingArray)
plt.title('Model at T = %.1f' %T)
plt.colorbar()
plt.show()

# show the magnetization of the model as a fucntion of temperature
plt.figure()
plt.plot(t,M)
plt.xlim(4,0)
plt.title('Mean Magnetization vs. Temperature')
plt.xlabel('Temperature')
plt.ylabel('Magnetization')
plt.show()
    
    # %%
