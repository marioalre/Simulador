# This is a GUI for the orbit class. It is not finished yet.

# App class for the GUI
# Path: scripts\GUI.py

import tkinter as tk
import webbrowser
from tkinter import *

import utils as U

def linkedin():
    webbrowser.open_new("linkedin.com/in/mario-álvarez-redondo")

def github():
    webbrowser.open_new("https://github.com/marioalre")


class App(tk.Tk):
    def __init__(self):
        super().__init__()

        ## Setting up Initial Things
        self.title("GUI UniLeón")
        self.geometry("720x550")
        self.resizable(True, True)
        self.iconphoto(True, tk.PhotoImage(file="scripts/images/escudoEIIIA.png"))
    
        ## Creating a container
        container = tk.Frame(self) # bg white
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(10, weight=1)
        container.grid_columnconfigure(10, weight=1)


        ## Initialize Frames
        self.frames = {}
        self.HomePage = HomePage
        self.PreDetOrb = PreDetOrb
        self.Mission = Mission

        ## Defining Frames and Packing it
        for F in {HomePage, PreDetOrb, Mission}:
            frame = F(self, container)
            self.frames[F] = frame
            frame.grid(row=10, column=10, sticky="nsew")    
           
        self.show_frame(HomePage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        menubar = frame.create_menubar(self)
        self.configure(menu=menubar)
        frame.tkraise()                         ## This line will put the frame on front
 



#---------------------------------------- HOME PAGE FRAME / CONTAINER ------------------------------------------------------------------------

class HomePage(tk.Frame):
    def __init__(self, parent, container):
        super().__init__(container)

        label = tk.Label(self, text="Homme page", font=('Times', '20'))
        label.pack()

        ## ADD CODE HERE TO DESIGN THIS PAGE

    def create_menubar(self, parent):
        menubar = Menu(parent, bd=3, relief=RAISED, activebackground="#80B9DC")

        ## Filemenu
        filemenu = Menu(menubar, tearoff=0, relief=RAISED, activebackground="#026AA9")
        menubar.add_cascade(label="Menú", menu=filemenu)
        filemenu.add_command(label="Mission", command=lambda: parent.show_frame(parent.Mission))
        filemenu.add_command(label="Órbitas", command=lambda: parent.show_frame(parent.PreDetOrb))
        filemenu.add_command(label="Close", command=lambda: parent.show_frame(parent.HomePage))
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=parent.quit)  

        ## proccessing menu
        processing_menu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Validation", menu=processing_menu)
        processing_menu.add_command(label="validate")
        processing_menu.add_separator()

        ## help menu
        help_menu = Menu(menubar, tearoff=0, relief=RAISED, activebackground="#026AA9")
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command = github)
        help_menu.add_command(label="Contancto", command= linkedin)
        help_menu.add_separator()

        return menubar


#---------------------------------------- PreDetOrb / CONTAINER ------------------------------------------------------------------------

class PreDetOrb(tk.Frame):
    def __init__(self, parent, container):
        super().__init__(container)

        label = tk.Label(self, text="Determinnación preliminar de órbitas", font=('Times', '20'))
        label.pack(pady=0,padx=0)

        ## ADD CODE HERE TO DESIGN THIS PAGE

    def create_menubar(self, parent):
        menubar = Menu(parent, bd=3, relief=RAISED, activebackground="#80B9DC")

        ## Filemenu
        filemenu = Menu(menubar, tearoff=0, relief=RAISED, activebackground="#026AA9")
        menubar.add_cascade(label="Menú", menu=filemenu)
        filemenu.add_command(label="Mission", command=lambda: parent.show_frame(parent.Mission))
        filemenu.add_command(label="Close", command=lambda: parent.show_frame(parent.HomePage))
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=parent.quit)  

        ## proccessing menu
        processing_menu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Validation", menu=processing_menu)
        processing_menu.add_command(label="validate")
        processing_menu.add_separator()

        ## help menu
        help_menu = Menu(menubar, tearoff=0, relief=RAISED, activebackground="#026AA9")
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command = github)
        help_menu.add_command(label="Contancto", command= linkedin)
        help_menu.add_separator()

        return menubar
    
#---------------------------------------- Mission PAGE FRAME / CONTAINER ------------------------------------------------------------------------
    
class Mission(tk.Frame):
    def __init__(self, parent, container):
        super().__init__(container)

        label = tk.Label(self, text="Configuración de la misión", font=('Times', '20'))
        label.pack(pady=0,padx=0)

        ## ADD CODE HERE TO DESIGN THIS PAGE

    def create_menubar(self, parent):
        menubar = Menu(parent, bd=3, relief=RAISED, activebackground="#80B9DC")

        ## Filemenu
        filemenu = Menu(menubar, tearoff=0, relief=RAISED, activebackground="#026AA9")
        menubar.add_cascade(label="Menú", menu=filemenu)
        filemenu.add_command(label="Órbitas", command=lambda: parent.show_frame(parent.PreDetOrb))
        filemenu.add_command(label="Close", command=lambda: parent.show_frame(parent.HomePage))
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=parent.quit)  

        ## proccessing menu
        processing_menu = Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Validation", menu=processing_menu)
        processing_menu.add_command(label="validate")
        processing_menu.add_separator()

        ## help menu
        help_menu = Menu(menubar, tearoff=0, relief=RAISED, activebackground="#026AA9")
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command = github)
        help_menu.add_command(label="Contancto", command= linkedin)
        help_menu.add_separator()

        return menubar

def run():
    app = App()
    app.mainloop()

if __name__ == "__main__":
    app = App()
    app.mainloop()



'''
Enlaces útiles:
Diseño de tablas:
    https://pythonguides.com/python-tkinter-table-tutorial/

'''
