from concurrent import futures
from glob import glob
from astropy.io import fits
import os
import numpy as np
from time import time
from concurrent.futures import ThreadPoolExecutor

A, B = 2, 2


def split(data: np.ndarray, i: int) -> np.ndarray:
    h, w = data.shape
    x, y = h//A, w//B
    n, m = i//B, i % B
    data = data[n*x:(n+1)*x,
                m*y:(m+1)*y].copy()
    return data


def concatenate(data_list):
    data_u = np.concatenate(data_list[:B], axis=1)
    data_d = np.concatenate(data_list[B:], axis=1)
    data = np.concatenate([data_u, data_d], axis=0)
    return data


def array_combine(ndarray, mode="mean") -> np.ndarray:
    """ Function to combine 3-D data array

    Parameters
    ----------
    ndarray: array, input data cube (3D)
    model: mean, median, sum, mean_clip, median_clip, default is mean
    """
    if mode == "median":
        array = np.median(ndarray, axis=0)
    elif mode == "median_clip":
        ndarray = np.sort(ndarray, axis=0)[1:-1]
        array = np.median(ndarray, axis=0)
    elif mode == "sum":
        array = np.sum(ndarray, axis=0)
    elif mode == "mean":
        array = np.mean(ndarray, axis=0)
    elif mode == "mean_clip":
        ndarray = np.sort(ndarray, axis=0)[1:-1]
        array = np.mean(ndarray, axis=0)
    return array


def load_bias(path: str, i: int) -> np.ndarray:
    with fits.open(path) as hdul:
        du = hdul[1].data
    du = du.astype(int)
    du = split(du, i)
    return du


def load_dark(path: str, i: int, bias: np.ndarray) -> np.ndarray:
    with fits.open(path) as hdul:
        du = hdul[1].data
        hu = hdul[0].header
    du = du.astype(int)
    du = du - bias
    du = du / hu["EXPTIME"]
    du = split(du, i)
    return du


def load_flat(path: str, i: int, bias: np.ndarray, dark: np.ndarray) -> np.ndarray:
    with fits.open(path) as hdul:
        du = hdul[1].data
        hu = hdul[0].header
    du = du.astype(int)
    du = du - bias - dark * hu["EXPTIME"]
    du = du / hu["EXPTIME"]
    du = du / np.median(du)
    du = split(du, i)
    return du


def save_fits(func, mode: str, path_list, save_path, *args) -> np.ndarray:
    ch_list = []
    with ThreadPoolExecutor() as tpe:
        for i in range(A*B):
            futures = [tpe.submit(func, path, i, *args) for path in path_list]
            du_list = [future.result() for future in futures]
            du = array_combine(du_list, mode)
            ch_list.append(du)
    du = concatenate(ch_list)
    du_output = os.path.basename(path_list[0])
    du_output = du_output.replace("raw", "combine")
    du_output = os.path.join(save_path, du_output)
    du = du.astype(np.float32)
    du_fits = fits.HDUList([fits.PrimaryHDU(data=du)])
    du_fits.writeto(du_output, overwrite=True)
    return du


def main(input_path : str, save_path : str, number_list : list):
    for number in number_list:
        print(number, end=' ')
        bias_path = glob(os.path.join(input_path, "MSC*/*CLB*_" + number + '_*'))
        bias = save_fits(load_bias, "median", bias_path, save_path)
        print('bias finish', end=' ')
        dark_path = glob(os.path.join(input_path, "MSC*/*CLD*_" + number + '_*'))
        dark = save_fits(load_dark, "median", dark_path, save_path, bias)
        print('dark finish', end=' ')
        flat_path = glob(os.path.join(input_path, "MSC*/*CLF*_" + number + '_*'))
        flat = save_fits(load_flat, "median", flat_path, save_path, bias, dark)
        print('flat finish')


if __name__ == "__main__":
    input_path = "/data/test20211012/ref/"
    save_path = "/data/test20211012/ref/"
    number_list = ['06', '07', '08', '09', '11', '12', '13', '14',
                   '15', '16', '17', '18', '19', '20', '22', '23', '24', '25']
    main(input_path, save_path, number_list)

