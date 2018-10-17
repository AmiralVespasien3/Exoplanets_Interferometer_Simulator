# Exoplanet Imaging Yield Simulator

This code is written in Python 3.5 and computes the list of known exoplanets that would be observable with an imaging or interferometric instrument based on assumptions on the wavelength, inner working angle and contrast. It also produces overview plots as shown below.

The planetery flux is the sum of the Planck thermal flux and the reflected flux assuming a bond albedo of 0.3 (Earth or Jupiter values) 

![equattion_flux](https://user-images.githubusercontent.com/43030278/45496411-ac070d00-b775-11e8-9b6d-5d5546ca4472.png)

Before launching the code, please download the most recent planets.csv file on the website : https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=planets

The instructions to download the code are as follows:
1/ Click on the 'Download Table ' button at the top left of the page
2/ Tick the boxes : CSV format, Download All Rows, Download All columns, Values Only (...)
3/ Click on 'Download Table'

To launch the code, please locate yourself in the directory containing the 'main.ipynb' file and use the command : jupyter-notebook

![fig1](https://user-images.githubusercontent.com/43030278/45622134-3e175a00-ba83-11e8-8a94-1d35d35bd594.png)
![figure3](https://user-images.githubusercontent.com/43030278/45497382-effb1180-b777-11e8-8836-67a8abe94be1.png)
![figure4](https://user-images.githubusercontent.com/43030278/45350412-06f10680-b5b4-11e8-9282-579457a1ea6e.png)
