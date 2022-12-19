# This is a GUI for the orbit class. It is not finished yet.

# App class for the GUI
# Path: scripts\GUI.py

import tkinter as tk
from tkinter import ttk
import webbrowser

def oppen_web():
    webbrowser.open_new("linkedin.com/in/mario-álvarez-redondo")


class GUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Orbit")
        self.geometry("1080x720")
        self.resizable(True, True)
        self.config(bg="white")
        self.iconbitmap("scripts\images\escudoEIIIA.ico")

        self.create_widgets()

    def create_widgets(self):
        self.create_menu()
        self.create_notebook()

    def create_menu(self):
        self.menu = tk.Menu(self)
        self.config(menu=self.menu)

        self.file_menu = tk.Menu(self.menu, tearoff=0)
        self.menu.add_cascade(label="File", menu=self.file_menu)
        self.file_menu.add_command(label="New")
        self.file_menu.add_command(label="Open")
        self.file_menu.add_command(label="Save")
        self.file_menu.add_command(label="Save as...")
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command=self.destroy)

        self.edit_menu = tk.Menu(self.menu, tearoff=0)
        self.menu.add_cascade(label="Edit", menu=self.edit_menu)
        self.edit_menu.add_command(label="Undo")
        self.edit_menu.add_command(label="Redo")
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Cut")
        self.edit_menu.add_command(label="Copy")
        self.edit_menu.add_command(label="Paste")
        self.edit_menu.add_command(label="Delete")
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Select All")

        self.help_menu = tk.Menu(self.menu, tearoff=0)
        self.menu.add_cascade(label="Help", menu=self.help_menu)
        self.help_menu.add_command(label="Help Index")
        self.help_menu.add_command(label="Contacto", command=oppen_web)

    def create_notebook(self):
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True)

        self.create_tab1()
        self.create_tab2()
        self.create_tab3()

    def create_tab1(self):
        self.tab1 = ttk.Frame(self.notebook)
        self.notebook.add(self.tab1, text="Estudio de óbitas")

    def create_tab2(self):
        self.tab2 = ttk.Frame(self.notebook)
        self.notebook.add(self.tab2, text="Misiones espaciales")

    def create_tab3(self):
        self.tab3 = ttk.Frame(self.notebook)
        self.notebook.add(self.tab3, text="Herraientas")

if __name__ == "__main__":
    gui = GUI()
    gui.mainloop()
