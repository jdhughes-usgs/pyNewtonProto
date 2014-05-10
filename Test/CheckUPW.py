__author__ = 'JosephHughes'

import os
import numpy as np
import matplotlib.pyplot as plt

def quadScaling(bratio, quadsfactor=1e-5):
    if bratio <= 0.0:
        sf = 0.
    else:
        a = 1. / (1 - quadsfactor)
        if bratio <= quadsfactor:
            sf = 0.5 * a * (bratio**2) / quadsfactor
        elif bratio > quadsfactor and bratio <= (1. - quadsfactor):
            sf = a * bratio + 0.5 * (1. - a)
        elif bratio > (1. - quadsfactor) and bratio < 1.:
            sf = 1. - ((0.5 * a * ((1. - bratio)**2)) / quadsfactor)
        else:
            sf = 1.
    return sf

def anald(bratio, quadsfactor=1e-5):
    tup = 1.
    a = 1. / (1. - quadsfactor)
    if bratio <= 0.0:
        sf = 0.
    elif bratio > 0. and bratio <= quadsfactor:
        #sf = a * bratio / (quadsfactor * tup)
        sf = a * bratio / (quadsfactor)
    elif bratio > quadsfactor and bratio <= (1. - quadsfactor):
        #sf = a / tup
        sf = a
    elif bratio > (1. - quadsfactor) and bratio < 1.:
        #x = 1. - bratio
        #y = -a * x / (quadsfactor * tup)
        #sf = 1. - y
        ##sf = 1. - ((-a * (1. - bratio)) / (quadsfactor * tup))
        sf = a * (1. - bratio) / quadsfactor
    else:
        sf = 0.
    return sf


def main():
    newparams = {'legend.fontsize':8, 'axes.labelsize':10, 'xtick.labelsize':8, 'ytick.labelsize':8}
    plt.rcParams.update(newparams)

    clr = ['black', 'blue', 'green', 'red']
    dr = 0.0005
    bratios = np.arange(0., 1.+dr, dr)
    v = np.zeros((bratios.shape), np.float)
    quadfactors = [0., 1.e-2, 5.e-2, 1.e-1]

    #--quadratic scaling
    fig = plt.figure(figsize=(6, 6), facecolor='#EDDA74', edgecolor='#EDDA74')
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
    ax2 = fig.add_axes([0.6, 0.2, 0.25, 0.25]) # inset axes
    for iclr, quad in enumerate(quadfactors):
        for idx, b in enumerate(bratios):
            v[idx] = quadScaling(b, quadsfactor=quad)
        ax1.plot(bratios, v, color=clr[iclr], linewidth=1.0, label='{:5.3g}'.format(quad))
        ax2.plot(bratios, v, color=clr[iclr], linewidth=1.0)
    ax1.plot([0.0, 0.1, 0.1], [0.1, 0.1, 0.], color='0.5', linewidth=0.9, label='_None')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_xlabel('Aquifer saturation')
    ax1.set_ylabel('Scaling factor')
    ax2.set_xlim(0, 0.1)
    ax2.set_ylim(0, 0.1)
    ax2.set_xlabel('Aquifer saturation')
    ax2.set_ylabel('Scaling factor')
    ax1.legend(loc='upper left')
    fig.savefig(os.path.join('..', 'Figures', 'quadscaling.png'), dpi=300)

    #--analytical derivatives
    fig = plt.figure(figsize=(6, 6), facecolor='#EDDA74', edgecolor='#EDDA74')
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
    ax2 = fig.add_axes([0.375, 0.2, 0.25, 0.25]) # inset axes
    for iclr, quad in enumerate(quadfactors):
        for idx, b in enumerate(bratios):
            v[idx] = anald(b, quadsfactor=quad)
        ax1.plot(bratios, v, color=clr[iclr], linewidth=1.0, label='{:5.3g}'.format(quad))
        ax2.plot(bratios, v, color=clr[iclr], linewidth=1.0)
    ax1.plot([0.0, 0.1, 0.1], [1.2, 1.2, 0.], color='0.5', linewidth=0.9, label='_None')
    ax1.set_xlim(0, 1)
    #ax1.set_ylim(0, 1)
    ax1.set_xlabel('Aquifer saturation')
    ax1.set_ylabel('Analytical derivative')
    ax2.set_xlim(0, 0.1)
    ax2.set_ylim(0, 1.2)
    ax2.set_xlabel('Aquifer saturation')
    ax2.set_ylabel('Analytical derivative')
    ax1.legend(loc='lower right')
    fig.savefig(os.path.join('..', 'Figures', 'anald.png'), dpi=300)

    #--numerical derivatives
    hups = 1.e-7
    fig = plt.figure(figsize=(6, 6), facecolor='#EDDA74', edgecolor='#EDDA74')
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8]) # main axes
    ax2 = fig.add_axes([0.375, 0.2, 0.25, 0.25]) # inset axes
    for iclr, quad in enumerate(quadfactors):
        for idx, b in enumerate(bratios):
            v1 = quadScaling(b, quadsfactor=quad)
            v2 = quadScaling(b+hups, quadsfactor=quad)
            v[idx] = (v2 - v1) / hups
        ax1.plot(bratios, v, color=clr[iclr], linewidth=1.0, label='{:5.3g}'.format(quad))
        ax2.plot(bratios, v, color=clr[iclr], linewidth=1.0)
    ax1.plot([0.0, 0.1, 0.1], [1.2, 1.2, 0.], color='0.5', linewidth=0.9, label='_None')
    ax1.set_xlim(0, 1)
    #ax1.set_ylim(0, 1)
    ax1.set_xlabel('Aquifer saturation')
    ax1.set_ylabel('Numerical derivative')
    ax2.set_xlim(0, 0.1)
    ax2.set_ylim(0, 1.2)
    ax2.set_xlabel('Aquifer saturation')
    ax2.set_ylabel('Numerical derivative')
    ax1.legend(loc='lower right')
    fig.savefig(os.path.join('..', 'Figures', 'numerd.png'), dpi=300)

if __name__ == '__main__':

    main()
