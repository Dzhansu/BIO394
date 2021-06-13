#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Import libraries
from typing import TextIO

import numpy as np
import random
import math
#from bokeh.plotting import figure, output_file, show
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # Set values:
    grid_boundaries = 50
    # Total number of cells
    number_cells = math.floor(grid_boundaries * grid_boundaries * 0.7)
    # Total number of macrophages
    number_macrophage = math.floor(grid_boundaries * grid_boundaries * 0.1)
    # Initial number of plaques
    number_plaques = 11
    # Initial number of drug molecules
    init_number_molecule = 50000
    # Total number of steps a molecule is allowed to move
    steps_molecule = 500
    # Total number of steps a macrophage is allowed to move
    steps_macrophage = 300
    # How many molecules are needed to remove a plaque
    removal_threshold = 40

    cell_x = np.random.randint(0, grid_boundaries, size=number_cells)
    cell_y = np.random.randint(0, grid_boundaries, size=number_cells)
    cell_map = list(map(list, zip(cell_x, cell_y)))

    print('number of cells: {}'.format(len(cell_map)))

    p.circle(cell_x, cell_y, size=20, color="navy", alpha=0.5)

    # Initialize array for plaque
    # [0] = x
    # [1] = y
    # [2] = number of bound molecules
    # [3] = status of the plaque: 0 means not removed, 1 means removable, 2 means removed
    plaque_map = [[0 for i in range(4)] for j in
                  range(number_plaques)]  # also written [[0] * 3] * number_plaques, i.e. [[0] * cols] * rows
    # Create an extra list to later compare with the position of the macrophages
    plaque_list = []
    plot_plaque_x = []
    plot_plaque_y = []
    for i in range(number_plaques):
        while True:
            plaque_map[i][0] = int(random.randint(0, grid_boundaries))
            plaque_map[i][1] = int(random.randint(0, grid_boundaries))
            if [plaque_map[i][0], plaque_map[i][1]] not in cell_map:
                break
        plaque_list.append([plaque_map[i][0], plaque_map[i][1]])
        plot_plaque_x.append(plaque_map[i][0])
        plot_plaque_y.append(plaque_map[i][1])
        plaque_map[i][2] = 0
        plaque_map[i][3] = 0

    print('number of plaques: {}'.format(len(plaque_map)))

    p.circle(plot_plaque_x, plot_plaque_y, size=10, color="gold", alpha=0.5)

    # Initialize array for macrophage
    # [0] = x
    # [1] = y
    # [2] = status of the macrophage: 0 if empty, 1 if containing a plaque
    macrophage_map = [[0 for i in range(3)] for j in
                      range(number_macrophage)]
    plot_macrophage_x = []
    plot_macrophage_y = []
    for i in range(number_macrophage):
        while True:
            macrophage_map[i][0] = int(random.randint(0, grid_boundaries))
            macrophage_map[i][1] = int(random.randint(0, grid_boundaries))
            if [macrophage_map[i][0], macrophage_map[i][1]] not in cell_map \
                    and [macrophage_map[i][0], macrophage_map[i][1]] not in plaque_list:
                break
        macrophage_map[i][2] = 0
        plot_macrophage_x.append(macrophage_map[i][0])
        plot_macrophage_y.append(macrophage_map[i][1])

    print('number of macrophage: {}'.format(len(macrophage_map)))

    p.circle(plot_macrophage_x, plot_macrophage_y, size=10, color="green", alpha=0.5)
    show(p)

    # here we simulate how many molecules we would need to reduce
    # the concentration of plaques to physiological level (23%)
    molecule = init_number_molecule
    print('Initial number of molecules: {}\n'.format(molecule))
    print('Start simulating molecule diffusion...\n')
    bound = 0
    molecule_plt = []

    for i in range(molecule):

        if i % 1000 == 0 and i > 0:
            print('{} molecules iteration'.format(i))

        #  Set where the molecule will enter the brain tissue
     
        x = np.random.randint(0, grid_boundaries)
        y = np.random.randint(0, grid_boundaries)

        # Here wer let the molecule move randomly for a given number of steps. If it doesn't bind
        # it will be considered as lost. Here, the molecule is not moving if it bumps into a cell
        free = 0
        number_of_steps = 0
        for v in range(steps_molecule):
            bound_1=0
            stuck = 0
            # We store [x,y] in a temporary variable in case the move ends in a cell. In that case we recover
            # the previous position
            coord_x, coord_y = x, y
            # Here we choose a direction
            choice = np.random.randint(4)
            if choice == 0:
                x += 1
            elif choice == 1:
                y += 1
            elif choice == 2:
                x -= 1
            elif choice == 3:
                y -= 1
            # Here we check the boundary conditions
            if x > grid_boundaries:
                x = 0
            if x < 0:
                x = grid_boundaries
            if y > grid_boundaries:
                y = 0
            if y < 0:
                y = grid_boundaries

            # If the molecule bumps into a cell it goes back where it was
            if [x, y] in cell_map:
                x, y = coord_x, coord_y
                stuck = 1

            # If co-localization with a plaque add one to that plaque.
            # We only go into this loop if the move was successful
            if stuck == 0:
                for k in range(len(plaque_map)):
                    if plaque_map[k][0] == x and plaque_map[k][1] == y:
                        plaque_map[k][2] += 1
                        bound_1 += 1
                        bound+=1
                        # the molecule is not free anymore
                        free = 1
                        break

            # If it did not bind it will continue to move
            if free == 0:
                number_of_steps += 1
            else:
                break

        molecule = molecule - bound_1

    # Checking which plaques can be removed
    # Here we change the status of the plaque and we count how many of them will be removed
    remove = 0
    for i in range(len(plaque_map)):
        if plaque_map[i][2] >= removal_threshold:
            plaque_map[i][3] = 1
            remove += 1
    print('{} plaque(s) can be removed'.format(remove))

    # Let the macrophage diffuse and find plaques marked to be removed
    print('Start simulating macrophage diffusion...\n')
    for i in range(number_macrophage):

        if i % 100 == 0 and i > 0:
            print('{} macrophages iteration'.format(i))

        # Here wer let the macrophage move randomly for a given number of steps. If it doesn't bind
        # it will be considered as lost. Here, again, the macrophage is not moving if it reaches the
        # boundaries
        free = 0
        number_of_steps = 0
        for v in range(steps_macrophage):

            stuck = 0
            # Again, we store [x,y] in a temporary variable in case the move ends in a cell. In that case we recover
            # the previous position
            coord_x, coord_y = x, y
            # Here we choose a direction
            choice = np.random.randint(4)
            if choice == 0:
                x += 1
            elif choice == 1:
                y += 1
            elif choice == 2:
                x -= 1
            elif choice == 3:
                y -= 1

            # Here we check the boundary conditions
            if x > grid_boundaries:
                x = 0
            if x < 0:
                x = grid_boundaries
            if y > grid_boundaries:
                y = 0
            if y < 0:
                y = grid_boundaries

            # If the macrophage bumps into a cell it goes back where it was
            if [x, y] in cell_map:
                x, y = coord_x, coord_y
                stuck = 1

            # If co-localization with a plaque change the status of the macrophage.
            # We only go into this loop if the move was successful and the macrophage
            # is empty (right now it is always the case).
            if stuck == 0 and macrophage_map[i][2] == 0:
                for k in range(len(plaque_map)):
                    if plaque_map[k][0] == x and plaque_map[k][1] == y and plaque_map[k][3] == 1:
                        plaque_map[k][3] = 2
                        macrophage_map[i][2] = 1
                        print('Bingo! Plaque {} was removed'.format(k))
                        break

            # If the macrophage did not bind it will continue to move
            if macrophage_map[i][2] == 0:
                number_of_steps += 1
            else:
                break

    # Count the number of removed plaques
    remove = 0
    for i in range(len(plaque_map)):
        if plaque_map[i][3] == 2:
            remove += 1
    print('{} plaque(s) were removed'.format(remove))

    print('Simulation ended properly!')


    exit(0)
