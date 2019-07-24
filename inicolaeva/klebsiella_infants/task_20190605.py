#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

matrix = np.random.randint(100, size=(100, 10))

named_matrix = dict()
counter = 0
for column in matrix:
    counter += 1
    named_matrix["row_" + str(counter)] = column.copy()

dataframe = pd.DataFrame(named_matrix)

fig, ax = plt.subplots()
ax.boxplot(dataframe)
ax.set_xticklabels(dataframe.columns)
plt.xlabel("Columns")
plt.ylabel("A random value")
plt.show()
