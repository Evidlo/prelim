#!/usr/bin/env python3

from glide.science.model_sph import *
from glide.science.forward_sph import *
from glide.common_components.camera import *
from glide.common_components.generate_view_geom import *

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import seaborn as sns

cams = [CameraL1BNFI(), CameraL1BWFI()]
cams[0].mode.n_frames = cams[1].mode.n_frames = 1
sc = gen_mission(num_obs=1, cams=cams)
grid = DefaultGrid()
m = Zoennchen24Model(grid=grid)
g = grid = m.grid
mr = SphHarmModel(grid, max_l=3)
density = m()

f = ForwardSph(sc, grid, use_aniso=True, use_albedo=True)
y = f.noise(density)

# %% plot

import matplotlib
matplotlib.use('Agg')

# plt.style.use('seaborn-v0_8-ticks') # I personally prefer seaborn for the graph style, but you may choose whichever you want.
params = {
    "ytick.color" : "black",
    "xtick.color" : "black",
    "axes.labelcolor" : "black",
    "axes.edgecolor" : "black",
    "text.usetex" : True,
    "font.family" : "serif",
    "font.serif" : ["Computer Modern Serif"],
    "font.size": 12
}
plt.rcParams.update(params)

plt.close()
plt.margins(0)
plt.figure(figsize=(4, 3), dpi=200)
plt.imshow(f.noisies[0])
plt.colorbar()
plt.tight_layout()
plt.savefig(out:=f'../figures/nfi_{cams[0].mode.n_frames}.png')
print(out)

plt.close()
plt.figure(figsize=(4, 3), dpi=200)
plt.imshow(f.noisies[1])
plt.colorbar()
plt.tight_layout()
plt.savefig(out:=f'../figures/wfi_{cams[1].mode.n_frames}.png')
print(out)

# plt.close()
# plt.figure(figsize=(3, 2), layout='constrained', dpi=200)
# plt.imshow(f.noisies[0])
# plt.colorbar()
# plt.savefig(out:=f'../figures/nfi_{cams[0].mode.n_frames}.png')
# print(out)