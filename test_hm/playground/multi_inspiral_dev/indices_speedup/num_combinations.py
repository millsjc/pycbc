import argparse
import numpy as np

def number_tail(t2_coinc_window, t3_coinc_window):
    n_tail = sum(
        (t2_coinc_window + 1 + np.arange(0, t2_coinc_window + 1))
        * (t3_coinc_window + 1 + np.arange(0, t2_coinc_window + 1))
    ) + sum(
        (2 * t2_coinc_window + 1)
        * (t3_coinc_window + 1 + np.arange(t2_coinc_window + 1, t3_coinc_window + 1))
    )
    return n_tail

def number_coincident_combinations(tlen, t2_coinc_window, t3_coinc_window):
    """Assumes t3_coinc_window>t2_coinc_window, and tlen is the number
    of time samples.
    """
    n_tail = number_tail(t2_coinc_window, t3_coinc_window)
    T_2 = 2 * t2_coinc_window + 1
    T_3 = 2 * t3_coinc_window + 1
    n_middle = (tlen - 2 * (t3_coinc_window + 1)) * T_2 * T_3
    n_c = 2 * n_tail + n_middle
    return n_c

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute 3-IFO indices.')
    parser.add_argument('tlen', type=int, help='length of time series')
    parser.add_argument('t2_coinc_window', type=int, help='coincidence window for T2')
    parser.add_argument('t3_coinc_window', type=int, help='coincidence window for T3')
    args = parser.parse_args()

    nind = number_coincident_combinations(args.tlen, args.t2_coinc_window, args.t3_coinc_window)
    print(nind) 
