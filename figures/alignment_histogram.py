# -*- coding: utf-8 -*-
"""
Created on Wed May 28 15:41:52 2025

@author: lafields2
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

cv_report = r"D:\Manuscripts\2025_NP_Database_Update\Csap_Genome_Components\cNP_db_matches.csv"
xls = pd.read_csv(cv_report)

plt.rcParams['font.family'] = 'Arial'

# Assume the last column in each sheet is the CV column
cv3 = xls.iloc[:, -1].dropna()

# Plot histogram for 3-plex CV with requested styling
plt.figure()
plt.hist(cv3, bins=40,
         color='#A4D7C6', alpha=0.5,
         edgecolor='#A4D7C6', linewidth=3)
plt.xlim(0, 100)
plt.xlabel('Similarity Score (%)')
plt.ylabel('Frequency')
#plt.title('Histogram of CV for 3-plex')
plt_path = r"D:\Manuscripts\2025_NP_Database_Update\Csap_Genome_Components\similarity_score_histogram.png"
plt.savefig(plt_path, dpi=300, bbox_inches='tight')
