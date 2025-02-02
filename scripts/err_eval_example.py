import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from expsweep.expsweep import merge_columns


d = pd.DataFrame([
    {'model': 'Fully Dense', 'Accuracy': 34, 'Precision': 40},
    {'model': 'Nested Sph. Harm.', 'Accuracy': 38, 'Precision': 45},
    {'model': 'Another Model', 'Accuracy': 24, 'Precision': 49},
])

d = merge_columns(d, ['Accuracy', 'Precision'], 'Metric', 'Error')

# ----- plotting -----

plt.rcParams.update({
    'text.usetex' : True,
    'font.size' : 16,
    'font.family' : 'lmodern',
    'font.weight' : 'bold',
    'figure.facecolor': (0, 0, 0, 0),
    'axes.facecolor': (1, 1, 1, 1),
    'figure.dpi': 150
})
plt.figure(figsize=(8, 4))


sns.barplot(data=d, x='model', y='Error', hue='Metric')
plt.axhline(y=20, color='blue', linestyle='--', label='Accuracy')
# plt.text(x=0, y=21, s='Worst Allowable Accuracy', color='blue')
plt.text(x=0, y=21, s='Worst Allowable Accuracy')
plt.axhline(y=50, color='darkorange', linestyle='--', label='Precision')
# plt.text(x=0, y=51, s='Worst Allowable Precision', color='darkorange')
plt.text(x=0, y=51, s='Worst Allowable Precision')

plt.ylim([0, 55])

plt.xlabel(None)
plt.title('Model Validation Error')
plt.tight_layout()

plt.savefig('err_eval_example.svg')
plt.close()
