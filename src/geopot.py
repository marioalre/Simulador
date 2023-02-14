import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import lpn, factorial, lpmn

class Geopot():
    '''Class to calculate the potential arround the Earth'''

    def __init__(self, body):
        self.mu = body.mu * 1e9  # Gravitational parameter [m^3/s^2]
        self.a = 6378136.3         # Earth equatorial radius [m]

    def gravitational_potential(self, r, elevation, azimuth, order=None):
        '''Gravitational potential 4.42
        Args:
            r: radius [m]
            elevation: elevation angle [rad]
            azimuth: azimuth angle [rad]
            order: order of the expansion
        Returns:
            Gravitational potential [m^2/s^2]
            Gravitational acceleration [m/s^2]
        '''
        
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

        if order is None:
            N_max = int(np.max([nn, mm]))
        else:
            N_max = order

        # Legendre polynomials
        P, Pd = self.legendre(N_max+1, np.sin(self.elev))

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

    def legendre(self, N, x):
        '''Legendre polynomials
        N: degree max
        m: order
        x: argument
        Returns: 
        Legendre polynomial normalized evaluated at x
        '''
        L, Ld= lpn(N, x)

        # P, Pd = lpmn(N, N, x)
        
        
        P = np.zeros((N, N))
        Pd = np.zeros((N, N))

        full_derivs = np.zeros((N+1,N))

        P[0, :] = L[:-1]
        Pd[0, :] = Ld[:-1]

        full_derivs[0, :] = Pd[0, :]

 
        for j in range(1, N+1):
            for i in range(2, N):
                n = i - 1
                full_derivs[j, i] = ((2*n + 1)*((j+1)*full_derivs[j-1, i-1] + x*full_derivs[j, i-1]) 
                                    - n*full_derivs[j, i-2])/(n + 1)

        for m in range(1, N):
            for n in range(m, N):
                '''
           #     if norm_type not in ["positive", "egm96"]:
           #         factor = (-1)**(m-1)
           #     else:
                '''
                factor = 1
                P[m, n] = factor*(1 - x**2)**((m)/2)*full_derivs[m-1, n]
                Pd[m, n] = -factor*((m)/2)*(1 - x**2)**(((m)/2)-1)*2*x + \
                        factor*(1-x**2)**((m)/2)*full_derivs[m, n]
        
        # Normalization

        for n in range(N):       # Order
            for m in range(N):   # Degree
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

    def calculate(self, resolucion=90, order=20, option='potential'):
        '''Calculate the potential and gravity field
        Parameters:
        ----------
        resolucion: int
            Resolution of the grid
        order: int
            Order of the expansion
        option: str
            Option to calculate the potential or the gravity field
        Returns:
        -------
        values: list
            List with the values of the potential or gravity field
        '''

        # Create some longitude, latitude, and third data
        elevacion = np.linspace(-180, 180, resolucion)
        elev_rad = elevacion* np.pi / 180
        azimuth = np.linspace(-90, 90, resolucion)
        azi_rad = azimuth* np.pi / 180

        values = []

        cont = 0

        for i, lat in enumerate(elev_rad):
            for j, lon in enumerate(azi_rad):
                print(f'Calculando {i} de {len(elevacion)} y {j} de {len(azimuth)}')
                self.gravitational_potential(self.a ,lon, lat, order=order)

                if option == 'potential':
                    third_data = self.pot
                    print(f'Potential: {third_data} m^2/s^2')
                elif option == 'gravity':
                    third_data = np.linalg.norm(self.grav) - self.mu / self.a**2
                    print(f'Acc: {third_data} m/s^2')
                else:
                    Warning('Option not valid, by default it will calculate the potential')
                    third_data = self.pot
                    print(f'Potential: {third_data} m^2/s^2')


                values.append([self.elev, self.azi, third_data]) 
                cont += 1


        elev , azi, third_data = [i[0] for i in values], [i[1] for i in values], [i[2] for i in values]
        # Numpy arrays
        elev = np.array(elev) # radians
        azi = np.array(azi)   # radians
        third_data = np.array(third_data)  

        # data = [elev, azi, third_data]

        return elev, azi, third_data

    def plot_potential(self, data, dtype='2D'):
        '''Plot the potential arround the Earth
        Parameters:
        ----------
        data: list
            List with the data to plot
        type: str
            Type of plot. 2D or 3D
        Returns:
        -------
        Plot of the potential
        '''

        elev, azi, third_data = data[0], data[1], data[2]

        # Check if it is a numpy array
        if type(elev) != np.ndarray:
            elev = np.array(elev)
        if type(azi) != np.ndarray:
            azi = np.array(azi)
        if type(third_data) != np.ndarray:
            third_data = np.array(third_data)

        error = 1

        while error == 1:   # Check if the input is correct
            if dtype == '2D':
                # Load the map image
                img = plt.imread('data\mapa.png')

                # Create a figure and an axis
                fig, ax = plt.subplots()

                # Show the image in the background of the axis
                ax.imshow(img, extent=[-180, 180, -90, 90], zorder=0)

                # Plot the longitude, latitude, and third data with 50% transparency
                sc = ax.scatter(azi*180/np.pi, elev*180/np.pi, c=third_data, cmap='viridis', alpha=0.5, zorder=1)

                # Add a color bar for the third data
                cbar = plt.colorbar(sc)
                cbar.set_label('Third data')

                # Set the limits and label of the axis
                ax.set_xlim(-180, 180)
                ax.set_ylim(-90, 90)
                ax.set_xlabel('Elevation')
                ax.set_ylabel('Azimuth')

                # Show the plot
                plt.show()
                error = 0

            elif dtype == '3D':

                x = np.cos(elev) * np.cos(azi)
                y = np.cos(elev) * np.sin(azi)
                z = np.sin(elev)

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

                # Muestra el gráfico
                plt.show()

                error = 0
            else:
                print('Type not valid')
                error = 1

    def write_data(self, name, list1, list2, list3):
        '''Write the data in a csv file
        Parameters
        ----------
        name : str
            Name of the file
        list1 : list
            List with the first data
        list2 : list
            List with the second data
        list3 : list
            List with the third data

        Returns
        -------
        None
        '''

        str_name = "data/" + str(name) + ".txt"

        with open('data/data.txt', 'w') as f:
            for i in range(len(list1)):
                line = str(list1[i]) + ',' + str(list2[i]) + ',' + str(list3[i]) + '\n'
                f.write(line)

        # Close the file
        f.close()

    def read_data(self, name):
        '''Read the data from a csv file
        Parameters
        ----------
        name : str
            Name of the file
        Returns
        -------
        data : list
            List with the data
        '''

        path = "data/" + str(name) + ".txt"

        with open(path, 'r') as f:
            data = f.readlines()
            
            f.close()

        # Remove the \n
        data = [i.replace('\n', '') for i in data]

        # Split the data and distribute it in 3 lists
        data = [i.split(',') for i in data]

        # change the data shape
        data = [[float(i[0]), float(i[1]), float(i[2])] for i in data]

        elev, azi, third_data = [i[0] for i in data], [i[1] for i in data], [i[2] for i in data]

        data = [elev, azi, third_data]

        return data


if __name__ == '__main__':
    from src.CelestialBodies import CelestialBodies

    Tierra = CelestialBodies()
    Tierra.earth()

    Potencial = Geopot(Tierra)

    # Potencial.gravitational_potential(7000000 ,1 , 1, 20)

    #data1, data2, data3 = Potencial.calculate(resolucion=90, order=20, option='gravity')

    #Potencial.write_data('gravity', data1, data2, data3)

    data = Potencial.read_data('data')
    #data = [data1, data2, data3]

    Potencial.plot_potential(data, dtype='2D')
    Potencial.plot_potential(data, dtype='3D')