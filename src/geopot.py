import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import lpmn, factorial

class Geopot():
    '''Class to calculate the potential arround the Earth'''

    def __init__(self, body):
        self.mu = body.mu * 1e9  # Gravitational parameter [m^3/s^2]
        self.a = 6378136.3         # Earth equatorial radius [m]

    def gravitational_potential(self, r, elevation, azimuth):
        '''Gravitational potential 4.42'''
        
        self.r = r                 # Radius [m]
        self.elev = elevation      # Elevation angle [rad]
        self.azi = azimuth         # Azimuth angle [rad]
        
        self.grav = np.array([-self.mu / self.r**2, 0, 0])
        self.pot = -self.mu / self.r
        self.pot_ideal = self.pot

        # Read the data
        data = np.loadtxt('data\egm96_to360.ascii.txt')
        n, m = np.shape(data)  # Number of rows and columns

        # Max column 1 and 2
        nn = np.max(data[:, 0])
        mm = np.max(data[:, 1])
        N_max = int(np.max([nn, mm]))
#        N_max = 2

        # Legendre polynomials
        P, Pd = self.legendre(N_max, np.sin(self.elev))

        for i in range(N_max):

            n = int(data[i, 0])
            m = int(data[i, 1])
            C = data[i, 2]
            S = data[i, 3]

            Pn = P[m, n]
            Pnd = Pd[m, n]
            

            self.pot = self.pot + self.mu /self.r * (self.a / self.r)**n * Pn * (C * np.cos(m * self.azi) + S * np.sin(m * self.azi))

            self.grav[0] = self.grav[0] - self.a**n * (n+1) * self.mu * Pn * (C*np.cos(m * self.azi) + S * np.sin(m * self.azi)) / self.r**(n+2)
            self.grav[1] = self.grav[1] + self.a**n*self.mu*Pnd*np.cos(self.elev)*(C*np.cos(m * self.azi) + S * np.sin(m * self.azi)) / self.r**(n+2)
            self.grav[2] = self.grav[2] + self.mu / self.r**(n+2) / np.sin(self.elev) * self.a**n * Pn * m *(S * np.cos(m * self.azi) - C * np.sin(m * self.azi))


            rot = np.array([[np.cos(self.azi)*np.sin(self.elev), np.cos(self.azi)*np.cos(self.elev), -np.sin(self.azi)],
                            [np.sin(self.azi)*np.sin(self.elev), np.sin(self.azi)*np.cos(self.elev), np.cos(self.azi)],
                            [np.cos(self.elev), -np.sin(self.elev), 0]])

            self.grav = np.dot(rot, self.grav)

            #print(self.pot)
            #print(self.grav)

        print(f'Gravedad en valor absoluto: {np.linalg.norm(self.grav)}')
        print(f'Potencial: {self.pot}')

    def legendre(self, N_max, x):
        '''Legendre polynomials
        N_max: degree max
        m: order
        x: argument
        Returns: 
        Legendre polynomial normalized evaluated at x
        '''
        P, Pd= lpmn(N_max, N_max, x)

        # Normalization

        for m in range(N_max):       # Order
            for n in range(N_max):   # Degree
                if n >= m:
                    factor = (2*n + 1)*factorial(n - m)/factorial(n + m)
                        
                    if m != 0:
                        factor = factor*2 
                        
                    factor = np.sqrt(factor)

                    P[m, n] = P[m, n] * factor
                    Pd[m, n] = Pd[m, n] * factor
                    

                # If it is nan or inf, replace it with 0
                if np.isnan(P[m, n]) or np.isinf(P[m, n]):
                        P[m, n] = 0
                if np.isnan(Pd[m, n]) or np.isinf(Pd[m, n]):
                        Pd[m, n] = 0
        
        return P, Pd

    def plot_potential(self, resolucion = 180, type='2D'):
        '''Plot the potential arround the Earth'''
        error = 1

        # Create some longitude, latitude, and third data
        longitudes = np.linspace(-180, 180, resolucion)
        latitudes = np.linspace(-90, 90, resolucion)

        third_data = np.zeros((len(latitudes), len(longitudes))) # rows, columns
        # Calculate the third data

        for i, lat in enumerate(latitudes):
            for j, lon in enumerate(longitudes):
                print(f'Calculando {i} de {len(latitudes)} y {j} de {len(longitudes)}')
                self.gravitational_potential(self.r ,lon, lat)
                third_data[i, j] = self.pot - self.pot_ideal
                print(f'Acc: {third_data[i, j]}')
        while error == 1:   # Check if the input is correct
            if type == '2D':
                # Load the map image
                img = plt.imread('data\mapa.png')

                # Create a figure and an axis
                fig, ax = plt.subplots()

                # Show the image in the background of the axis
                ax.imshow(img, extent=[-180, 180, -90, 90], zorder=0)

                # Plot the longitude, latitude, and third data with 50% transparency
                sc = ax.scatter(longitudes, latitudes, c=third_data, cmap='viridis', alpha=0.5, zorder=1)

                # Add a color bar for the third data
                cbar = plt.colorbar(sc)
                cbar.set_label('Third data')

                # Set the limits and label of the axis
                ax.set_xlim(-180, 180)
                ax.set_ylim(-90, 90)
                ax.set_xlabel('Longitude')
                ax.set_ylabel('Latitude')

                # Show the plot
                plt.show()
                error = 0

            elif type == '3D':
                
                x = np.sin(latitudes / 180.0 * np.pi) * np.cos(longitudes / 180.0 * np.pi)
                y = np.sin(latitudes / 180.0 * np.pi) * np.sin(longitudes / 180.0 * np.pi)
                z = np.cos(latitudes / 180.0 * np.pi)

                # Crea una figura y un eje 3D
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')

                # Dibuja una esfera como fondo
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 100)
                x_sphere = np.outer(np.cos(u), np.sin(v))
                y_sphere = np.outer(np.sin(u), np.sin(v))
                z_sphere = np.outer(np.ones(np.size(u)), np.cos(v))
                ax.plot_surface(x_sphere, y_sphere, z_sphere, color='gray', alpha=0.1, rstride=10, cstride=10)

                # Dibuja los datos de longitud, latitud y el tercer dato con una transparencia del 50%
                sc = ax.scatter(x, y, z, c=third_data, cmap='viridis', alpha=0.5)

                # Añade una barra de color para el tercer dato
                cbar = plt.colorbar(sc)
                cbar.set_label('Tercer dato')

                # Configura los límites del eje y la etiqueta
                ax.set_xlim(-1, 1)
                ax.set_ylim(-1, 1)
                ax.set_zlim(-1, 1)
                ax.set_xlabel('Eje X')
                ax.set_ylabel('Eje Y')
                ax.set_zlabel('Eje Z')

                error = 0
            else:
                print('Type not valid')
                error = 1


if __name__ == '__main__':
    from src.CelestialBodies import CelestialBodies

    Tierra = CelestialBodies()
    Tierra.earth()

    Potencial = Geopot(Tierra)

    Potencial.gravitational_potential(7000000 ,1 , 1)

    Potencial.plot_potential()