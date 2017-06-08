import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import convolve

def getMACF(lc):

    return convolve(lc,lc)

def main():
    getMACF()


if __name__ == '__main__':
    main()
