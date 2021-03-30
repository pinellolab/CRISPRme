'''
Various MessageBox functions
'''
import os, shutil
from os.path import isdir
try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True
from tkinter import messagebox

def deleteResultConfirm(jobID='test'):
    current_working_directory = os.getcwd() + '/'
    root = tk.Tk()
    root.withdraw()
    MsgBox = tk.messagebox.askquestion ('Delete Result','Are you sure you want to delete the ' + jobID + ' result',icon = 'warning')
    removed = False
    if MsgBox == 'yes':
        if not isdir(current_working_directory + 'Results/' + jobID):
            tk.messagebox.showerror('ERROR', 'Directory ' + jobID + ' not found.')
        else:
            try:
                shutil.rmtree(current_working_directory + 'Results/' + jobID)
            except PermissionError:
                tk.messagebox.showerror('ERROR', 'Permission denied.')
            except FileNotFoundError:
                tk.messagebox.showerror('ERROR', 'Directory ' + jobID + ' not found.')
            else:
                tk.messagebox.showinfo('Result Deleted','The selected result was removed.')
                removed = True
    root.destroy() 
    return removed  

if __name__ == '__main__':
    deleteResultConfirm()