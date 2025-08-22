#!/usr/bin/env python3

import matplotlib.pyplot as plt
import torch as t
import matplotlib
matplotlib.use('Agg')
plt.close('all')

from glide.common_components.generate_view_geom import gen_mission
from glide.common_components.orbits import carruthers_orbit
from glide.common_components.camera import CameraL1BWFI, CameraL1BNFI
from glide.science.forward_sph import sc2vg



def add_headers(
    fig,
    *,
    row_headers=None,
    col_headers=None,
    row_pad=1,
    col_pad=5,
    rotate_row_headers=True,
    **text_kwargs
):
    # Based on https://stackoverflow.com/a/25814386

    axes = fig.get_axes()

    for ax in axes:
        sbs = ax.get_subplotspec()

        # Putting headers on cols
        if (col_headers is not None) and sbs.is_first_row():
            ax.annotate(
                col_headers[sbs.colspan.start],
                xy=(0.5, 1),
                xytext=(0, col_pad),
                xycoords="axes fraction",
                textcoords="offset points",
                ha="center",
                va="baseline",
                **text_kwargs,
            )

        # Putting headers on rows
        if (row_headers is not None) and sbs.is_first_col():
            ax.annotate(
                row_headers[sbs.rowspan.start],
                xy=(0, 0.5),
                xytext=(-ax.yaxis.labelpad - row_pad, 0),
                xycoords=ax.yaxis.label,
                textcoords="offset points",
                ha="right",
                va="center",
                rotation=rotate_row_headers * 90,
                **text_kwargs,
            )



intervals = 1, 15, 30, 90
seasons = ('spring', 'summer', 'fall', 'winter')
cams = (CameraL1BWFI(), CameraL1BNFI())

figs = []

for interval in intervals:
    fig, axs = plt.subplots(2, 4, figsize=(10, 5), dpi=150, constrained_layout=True)

    fig_row = []
    # ax_cam - plot axes for each camera
    # cam - camera
    # d - LOS downsample factor
    for ax_cam, cam, d in zip(axs, cams, (16, 64)):

        print(f'{interval=}, {cam=}')

        for ax, season in zip(ax_cam, seasons):

            # generate LOS plot with only 2 observations
            sc = carruthers_orbit(start=season, duration=interval, cam=cam, num_obs=2)
            vg = sc2vg(sc)
            # compute start and end of all LOS
            start = vg.ray_starts
            end = vg.ray_starts + 2 * vg.rays * t.linalg.norm(vg.rays, dim=-1, keepdim=True)
            # shape of camera
            N = cam.npix

            for c, s, e in zip(['b', 'r'], start, end):
                # plot every dth pixel of middle row of pixels
                for e_d in e[cam.npix//2, ::d, :]:
                    ax.axline(s[0, 0, :2], e_d[:2], color=c)

                # ax.axis('equal')
                width = 25
                ax.set_xlim([-width, width])
                ax.set_ylim([-width, width])
                ax.grid(True
                        )
                ax.minorticks_on()

    for ax in axs[1]:
        ax.set_xlabel('Re')

    add_headers(
        fig,
        row_headers=list(map(lambda c: c.camID, cams)),
        col_headers=list(map(str.capitalize, seasons))
    )
    # plt.tight_layout()
    plt.suptitle(f'Interval: {interval}d')

    plt.savefig(f'/srv/www/los{interval}.png')
