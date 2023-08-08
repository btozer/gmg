import os
import subprocess 
from sys import platform

def launch():
    ## GET PATH TO gmgpy DIR
    path = os.path.dirname(os.path.abspath(__file__))

    ## RUN APP
    ## IF ON WINDOWS OR LINUX
    if platform == "win32" or platform == "linux" or platform == "linux2":
        subprocess.run(["pythonw", path+"/gmg.py"])
    ## IF ON MAC
    elif platform == "darwin":
        subprocess.run(["pythonw", path+"/gmg.py"])
    else:
        print("OS not supported")
        exit()

if __name__ == '__main__':
    launcher()