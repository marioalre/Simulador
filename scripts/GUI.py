# This is a GUI for the orbit class. It is not finished yet.

# Path: scripts\GUI.py

import tkinter as tk
from tkinter import ttk

class GUI:
    def __init__(self, master):
        self.master = master
        master.title("Orbit GUI")

        self.label = tk.Label(master, text="Orbit GUI")
        self.label.pack()

        self.greet_button = tk.Button(master, text="Greet", command=self.greet)
        self.greet_button.pack()

        self.close_button = tk.Button(master, text="Close", command=master.quit)
        self.close_button.pack()


    def greet(self):
        print('Hola')

root = tk.Tk()
my_gui = GUI(root)
root.mainloop()
