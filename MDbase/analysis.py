# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('/Users/takato/Simulation/MD/Temperatures.csv', names=['Ta', 'Tb'], header = 1)
MDstep = range(10000)
Ta = df.loc[:, "Ta"]
Tb = df.loc[:, "Tb"]

Tb.head()


plt.plot(MDstep, Tb) 
plt.show()