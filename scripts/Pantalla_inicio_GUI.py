from tkinter import *
import webbrowser
from GUI import *

def oppen_web():
    webbrowser.open_new("linkedin.com/in/mario-álvarez-redondo")

# crear la primera ventana
window = Tk()

# Personalizar la ventana
window.title("Mi primera ventana")
window.geometry("1080x720")
# window.maxsize(1080,720)
window.minsize(1080,360)

window.iconbitmap("scripts\images\escudoEIIIA.ico")

window.config(bg="#41B77F")
# Crear frame
frame = Frame(window, bg="#41B77F")

# texto de bienvenida
label_title = Label(frame, text="Simulador para el análisis de misiones espaciales", font=("Arial", 35), fg="white", bg="#41B77F")
label_title.pack()

# Segundo texto
label_subtitle = Label(frame, text="Autor: Mario Álvarez Redondo", font=("Arial", 18), fg="white", bg="#41B77F")
label_subtitle.pack()

frame.pack(expand=YES)

# Crear un botón
button = Button(frame, text="Contacto", command=oppen_web, fg="#41B77F", bg="white", font=("Arial", 18))
button.pack(side=BOTTOM, pady=20, padx=400)

# Crear un botón

def funcion():
    run()
    window.iconify()


button = Button(frame, text="Ejecutar", fg="#41B77F", bg="white", font=("Arial", 14), command=funcion)
button.pack(side=BOTTOM, pady=15, expand=400)

# Insertar imagen
image = PhotoImage(file="scripts\images\escudoEIIIA.png").subsample(2)
label_image = Label(frame, image=image, bg="#41B77F")
label_image.pack()

# Meter en el bucle principal
window.mainloop()