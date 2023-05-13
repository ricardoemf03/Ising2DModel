import numpy as np
import matplotlib.pyplot as plt
import random
import copy
#semilla de numero aleatorios
rdmsedd = random.seed(24)
rdd = np.random.seed(24) #semilla para numpy
# Tamaño de la red
N = 15
# Constante J
J = 1
#constante de Boltzmann
kb=1
# Número de iteraciones de Monte Carlo
nsteps = 200000
# Crear la red inicial de espines aleatorios
spins = np.random.choice([-1, 1], size=(N, N))
initial_conf = copy.deepcopy(spins)
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
temperatures = np.linspace(0.1, 5.0, 50)

# Lista de energías promedio
avg_energies = []
heat_cap = []
avg_mag = []
sus_mag = []
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
    avg_mag.append(np.mean(magnetizations)/N**2)
    sqrt_energies = [ener**2 for ener in energies]
    Eprom = np.mean(energies)
    E2prom = np.mean(sqrt_energies)
    heat_cap.append((E2prom-Eprom**2)/(N**2*kb*T**2))
    sqrt_mag = [mags**2 for mags in magnetizations]
    AbsMag = [abs(mm) for mm in magnetizations]
    Mprom = np.mean(AbsMag)
    M2prom = np.mean(sqrt_mag)
    sus_mag.append((M2prom-Mprom**2)/(N**2*kb*T))

    print("para T= ", T, " Eprom = ", 2*np.mean(energies)/N**2, " Mprom = ", np.mean(magnetizations)/N**2 , " capcal= ", (E2prom-Eprom**2)/(N**2*kb*T**2), "sus = ", (M2prom-Mprom**2)/(N**2*kb*T))
heat_cap[0]=0 #para que salga bien la grafica remplazamos el primer valor de heatcap ya que es muy grande
#heat_cap[10]=0
sus_mag[0]=0
final_conf=copy.deepcopy(spins)
#Mostrar los graficos
fig = plt.figure()
fig.clf()
ax = fig.subplots(2,3)

ax[0,0].plot(temperatures, avg_energies)
ax[0,0].set_xlabel('Temperature')
ax[0,0].set_ylabel('<E>')
ax[0,0].set_title('<E> vs T')

ax[0,1].plot(temperatures, heat_cap)
ax[0,1].set_xlabel('Temperature')
ax[0,1].set_ylabel('<C>')
ax[0,1].set_title('<C> vs T')

ax[1,0].plot(temperatures, avg_mag)
ax[1,0].set_xlabel('Temperature')
ax[1,0].set_ylabel('<M>')
ax[1,0].set_title('<M> vs T')

ax[1,1].plot(temperatures, sus_mag)
ax[1,1].set_xlabel('Temperature')
ax[1,1].set_ylabel('Magnetic susceptibility')
ax[1,1].set_title('Magnetic susceptibility vs Temperature')

ax[0,2].matshow(initial_conf)
ax[0,2].set_title('Initial configuration')

ax[1,2].matshow(final_conf)
ax[1,2].set_title('Final configuration')

fig.tight_layout()
fig.show()
plt.show()
#graficas individuales
plt.plot(temperatures, avg_energies)
plt.xlabel('Temperature')
plt.ylabel('<E>')
plt.title('<E> vs T')
plt.show()

plt.plot(temperatures, heat_cap)
plt.xlabel('Temperature')
plt.ylabel('<C>')
plt.title('<C> vs T')
plt.show()

plt.plot(temperatures, avg_mag)
plt.xlabel('Temperature')
plt.ylabel('<M>')
plt.title('<M> vs T')
plt.show()

plt.plot(temperatures, sus_mag)
plt.xlabel('Temperature')
plt.ylabel('Magnetic susceptibility')
plt.title('Magnetic susceptibility vs Temperature')
plt.show()

plt.matshow(initial_conf)
plt.title('Initial configuration')
plt.show()

plt.matshow(final_conf)
plt.title('Final configuration')
plt.show()
