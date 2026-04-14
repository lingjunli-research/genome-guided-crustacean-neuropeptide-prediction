# -*- coding: utf-8 -*-
"""
Created on Wed Nov  5 10:45:55 2025

@author: lafields2
"""

import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib import rcParams

# --- Make SVG text editable in Illustrator ---
rcParams['svg.fonttype'] = 'none'     # keep text as text (not outlines)
rcParams['font.family'] = 'sans-serif'  # use a common, AI-friendly font

# ---- Paths to your two text files ----
file1 = r"D:\Manuscripts\2025_inSilicoPrediction\Manuscript\Figure_Generation\VennDiagram\All_Cbor_Predicted.txt"
file2 = r"D:\Manuscripts\2025_inSilicoPrediction\Manuscript\Figure_Generation\VennDiagram\All_Csap_Predicted.txt"

# ---- Read lists from files ----
with open(file1, 'r', encoding='utf-8') as f:
    list1 = set(line.strip() for line in f if line.strip())

with open(file2, 'r', encoding='utf-8') as f:
    list2 = set(line.strip() for line in f if line.strip())

# ---- Create Venn diagram ----
plt.figure(figsize=(6,6))
v = venn2([list1, list2], set_labels=('C. borealis', 'C. sapidus'))
plt.title('Venn Diagram of Two Lists')

# ---- Export as SVG (Illustrator-friendly) ----
# transparent=False gives a white background; set to True if you want transparency
plt.tight_layout()
plt.savefig(r'D:\Manuscripts\2025_inSilicoPrediction\Manuscript\Figure_Generation\VennDiagram\venn_diagram.svg', format='svg', transparent=False)
plt.close()
