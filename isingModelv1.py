#Importing the necessary libraries and modules required for the program's functionality:

import numpy as np #Support for efficient numerical operations and array manipulation.
import matplotlib.pyplot as plt #To create charts, plots, and visual props.
import random #This provides multiple functions for generating random numbers, shuffling sequences, and making random selections.
import copy #Functions for creating shallow or deep copies of objects.
from matplotlib.colors import ListedColormap #For the black and white colors of the initial and final configurations.
cmap = ListedColormap(['k', 'w'])

#Define the parameters of the system:

rdmsedd = random.seed(24)#Random number seed.
rdd = np.random.seed(24) #Seed for NumPy.
N = 15 #Red size.
J = 1 #Constant J.
kb=1 #Boltzmann's constant.
nsteps = 200000 #Number of Monte Carlo iterations.
#Creating the initial lattice of random spins:
spins = np.random.choice([-1, 1], size=(N, N))
initial_conf = copy.deepcopy(spins) #Copy of the initial configuration to display later.
#Defining the function that calculates the energy:
def energy(spins):
    """Calculates the energy of the spin lattice"""
    #Calculating the sum of the products of adjacent spins:
    neighbors = np.roll(spins, 1, axis=0) + np.roll(spins, -1, axis=0) + np.roll(spins, 1, axis=1) + np.roll(spins, -1, axis=1)
    energy = -J * np.sum(spins * neighbors/2)
    return energy

#Defining the magnetization function:
def magnetization(spins):
    return np.sum(spins)

#Temperature list:
temperatures = np.linspace(0.1, 5.0, 90)

#Creating empty lists to store the average energies, heat capacities, average magnetizations, and susceptibility of magnetization:
avg_energies = [] #List of average energies for each T.
heat_cap = [] #List of average heat capacities for each T.
avg_mag = [] #List of average magnetizations for each T.
sus_mag = [] #List of average magnetic susceptibilities for each T.
for T in temperatures:
    #Monte Carlo energies for this temperature:
    energies = []
    #Monte Carlo magnetizations for this temperature:
    magnetizations = []
    for step in range(nsteps):
        #Select a random spin:
        i, j = random.randint(0, N-1), random.randint(0, N-1)
        #Calculate the energy of the current state:
        E = energy(spins)
        #Calculate the magnetization of the current state:
        M = magnetization(spins)
        #Changing the spin:
        spins[i, j] *= -1
        #Calculate the energy of the new state:
        Enew = energy(spins)
        #Calculate magnetization of the new state:
        Mnew = magnetization(spins)
        #Calculate the energy difference:
        deltaE = Enew - E
        #Accept or reject the change with the appropriate probability based on the Metropolis criterion:
        if deltaE <= 0 or np.exp(-deltaE /(kb*T)) > random.uniform(0, 1):
            energies.append(Enew)
            magnetizations.append(Mnew)
        else:
            #Revert the change:
            spins[i, j] *= -1
            energies.append(E)
            magnetizations.append(M)
    #Average energy for this temperature:
    avg_energies.append(2*np.mean(energies)/N**2)
    #Average magnetization for this temperature:
    avg_mag.append(np.mean(magnetizations)/N**2)
    sqrt_energies = [ener**2 for ener in energies] #The square of the energies.
    Eprom = np.mean(energies) #<E> for this T.
    E2prom = np.mean(sqrt_energies) #<E^2> for this T.
    #Average heat capacity for this temperature:
    heat_cap.append((E2prom-Eprom**2)/(N**2*kb*T**2))
    sqrt_mag = [mags**2 for mags in magnetizations] #Square of the magnetizations.
    AbsMag = [abs(mm) for mm in magnetizations] #Absolute value of the magnetizations.
    Mprom = np.mean(AbsMag) #<|M|>
    M2prom = np.mean(sqrt_mag) #<M^2>
    #Average magnetic susceptibility for this temperature:
    sus_mag.append((M2prom-Mprom**2)/(N**2*kb*T))

    print("para T= ", T, " Eprom = ", 2*np.mean(energies)/N**2, " Mprom = ", np.mean(magnetizations)/N**2 , " capcal= ", (E2prom-Eprom**2)/(N**2*kb*T**2), "sus = ", (M2prom-Mprom**2)/(N**2*kb*T))
heat_cap[0]=0 #In order to obtain a correct graph, we replace the first value of heatcap since it is very large.
sus_mag[0]=0 #Likewise the susceptibility.
final_conf=copy.deepcopy(spins) #Making a copy of the final configuration.
#To display the plots together:
fig = plt.figure()
fig.clf()
ax = fig.subplots(2,3)

ax[0,0].scatter(temperatures, avg_energies, color="black", s=10)
ax[0,0].set_xlabel('Temperature')
ax[0,0].set_ylabel('<E>')
ax[0,0].set_title('<E> vs T')

ax[0,1].scatter(temperatures, heat_cap, color="black", s=10)
ax[0,1].set_xlabel('Temperature')
ax[0,1].set_ylabel('<C>')
ax[0,1].set_title('<C> vs T')

ax[1,0].scatter(temperatures, avg_mag, color="black", s=10)
ax[1,0].set_xlabel('Temperature')
ax[1,0].set_ylabel('<M>')
ax[1,0].set_title('<M> vs T')

ax[1,1].scatter(temperatures, sus_mag, color="black", s=10)
ax[1,1].set_xlabel('Temperature')
ax[1,1].set_ylabel('Magnetic susceptibility')
ax[1,1].set_title('Magnetic susceptibility vs Temperature')

#Plot the configuration, using "plt.matshow": it visualizes a 2D array or matrix as a color-coded grid, where each cell's
#color represents the corresponding value in the matrix.

ax[0,2].matshow(initial_conf, cmap=cmap)
ax[0,2].set_title('Initial configuration')

ax[1,2].matshow(final_conf, cmap=cmap)
ax[1,2].set_title('Final configuration')

fig.tight_layout()
fig.show()
plt.show()

#To display each plot individually:
plt.scatter(temperatures, avg_energies, color="black", s=10)
plt.xlabel('Temperature')
plt.ylabel('<E>')
plt.title('<E> vs T')
plt.show()

plt.scatter(temperatures, heat_cap, color="black", s=10)
plt.xlabel('Temperature')
plt.ylabel('<C>')
plt.title('<C> vs T')
plt.show()

plt.scatter(temperatures, avg_mag, color="black", s=10)
plt.xlabel('Temperature')
plt.ylabel('<M>')
plt.title('<M> vs T')
plt.show()

plt.scatter(temperatures, sus_mag, color="black", s=10)
plt.xlabel('Temperature')
plt.ylabel('Magnetic susceptibility')
plt.title('Magnetic susceptibility vs Temperature')
plt.show()

plt.matshow(initial_conf, cmap=cmap)
plt.title('Initial configuration')
plt.show()

plt.matshow(final_conf, cmap=cmap)
plt.title('Final configuration')
plt.show()
