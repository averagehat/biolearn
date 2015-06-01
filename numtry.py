from functools import partial
import numpy as np
from func import compose

#next argument is the size
simulate_prizedoor = partial(np.random.randint, 0, 3)
random_col_vals = partial(np.apply_along_axis, np.random.choice, 1)
simulate_guess = np.ones
rowchoice = compose(np.random.choice, np.ma.compressed)
RUNS = 1000


def goat_doors(pzs, gss):
    grid = np.repeat(np.ma.arange(3), RUNS).reshape(3, RUNS)
    unpicked_matrix = (grid == pzs) | (grid == gss)
    grid.mask = unpicked_matrix
    return np.array(map(rowchoice, grid.T))

switch_guess = goat_doors

def win_percentage(pzs, gss):
    #return (pzs == gss).sum()/float(len(gss))
    return 100*(pzs == gss).mean()

def sim_game(switch=False):
    pzs, gss = simulate_prizedoor(RUNS), simulate_guess(RUNS)
    goats = goat_doors(pzs, gss)
    picks = switch_guess(gss, goats) if switch else gss
    return win_percentage(picks, pzs)

print sim_game(True)


