import numpy as np
import matplotlib.pyplot as plt
import random
import copy
from matplotlib.colors import ListedColormap #para colores blanco y negro de conf inicial y final
cmap = ListedColormap(['k', 'w']) #para colores blanco y negro de conf inicial y final
rdmsedd = random.seed(24)#semilla de numero aleatorios
rdd = np.random.seed(24) #semilla para numpy
N = 15 # Tamaño de la red
J = 1 # Constante J
kb=1 #constante de Boltzmann
nsteps = 200000 # Número de iteraciones de Monte Carlo
# Crear la red inicial de espines aleatorios
spins = np.random.choice([-1, 1], size=(N, N))
initial_conf = copy.deepcopy(spins) #copia de configuracion inicial para luego mostrar
#definir la funcion que calcula la energia
def energy(spins):
    """Calcula la energía de la red de espines"""
    # Calcular la suma de los productos de los espines adyacentes
    neighbors = np.roll(spins, 1, axis=0) + np.roll(spins, -1, axis=0) + np.roll(spins, 1, axis=1) + np.roll(spins, -1, axis=1)
    energy = -J * np.sum(spins * neighbors/2)
    return energy

# Definir la función de magnetización
def magnetization(spins):
    return np.sum(spins)

# Lista de temperaturas
temperatures = np.linspace(0.1, 5.0, 90)

avg_energies = [] # Lista de energías promedio para cada T
heat_cap = [] # Lista de capacidades calorificas promedio para cada T
avg_mag = [] # Lista de magnetizaciones promedio para cada T
sus_mag = [] # Lista de suceptibilidades magneticas promedio para cada T
for T in temperatures:
    # Energías de Monte Carlo para esta temperatura
    energies = []
    # Magnetizaciones de monte carlo para esta temperatura
    magnetizations = []
    for step in range(nsteps):
        # Seleccionar un espín aleatorio
        i, j = random.randint(0, N-1), random.randint(0, N-1)
        # Calcular la energía del estado actual
        E = energy(spins)
        #Calcular la magnetizacion del estado actual
        M = magnetization(spins)
        # Cambiar el espín
        spins[i, j] *= -1
        # Calcular la energía del nuevo estado
        Enew = energy(spins)
        #Caclular magnetizacion del nuevo estado
        Mnew = magnetization(spins)
        # Calcular la diferencia de energía
        deltaE = Enew - E
        # Aceptar o rechazar el cambio con la probabilidad adecuada
        if deltaE <= 0 or np.exp(-deltaE /(kb*T)) > random.uniform(0, 1):
            energies.append(Enew)
            magnetizations.append(Mnew)
        else:
            #Revertir el cambio
            spins[i, j] *= -1
            energies.append(E)
            magnetizations.append(M)
    # Energía promedio para esta temperatura
    avg_energies.append(2*np.mean(energies)/N**2)
    # Magnetizacion promedio para esta temperatura
    avg_mag.append(np.mean(magnetizations)/N**2)
    sqrt_energies = [ener**2 for ener in energies] # cuadrado de las energias
    Eprom = np.mean(energies) # <E> para esta T
    E2prom = np.mean(sqrt_energies) # <E^2> para esta T
    # capacidad calorifica promedio para esta temperatura
    heat_cap.append((E2prom-Eprom**2)/(N**2*kb*T**2))
    sqrt_mag = [mags**2 for mags in magnetizations] # cuadrado de las magnetizaciones
    AbsMag = [abs(mm) for mm in magnetizations] # valor absoluto de las magnetizaciones
    Mprom = np.mean(AbsMag) # <|M|>
    M2prom = np.mean(sqrt_mag) # <M^2>
    # susceptibilidad magnetica promedio para esta temperatura
    sus_mag.append((M2prom-Mprom**2)/(N**2*kb*T))

    print("para T= ", T, " Eprom = ", 2*np.mean(energies)/N**2, " Mprom = ", np.mean(magnetizations)/N**2 , " capcal= ", (E2prom-Eprom**2)/(N**2*kb*T**2), "sus = ", (M2prom-Mprom**2)/(N**2*kb*T))
heat_cap[0]=0 # para que salga bien la grafica remplazamos el primer valor de heatcap ya que es muy grande
sus_mag[0]=0 # de igual manera la susceptibilidad
final_conf=copy.deepcopy(spins) # hacemos una copia de la configuracion final
# Para mostrar los graficos juntos
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

ax[0,2].matshow(initial_conf, cmap=cmap)
ax[0,2].set_title('Initial configuration')

ax[1,2].matshow(final_conf, cmap=cmap)
ax[1,2].set_title('Final configuration')

fig.tight_layout()
fig.show()
plt.show()

# Para mostrar cada grafica por individual
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
