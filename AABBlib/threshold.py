#!/usr/bin/env python3
from scipy import ndimage
from multiprocessing import Process, Manager
import numpy as np


class Threshold(object):
    """docstring for Threshold"""
    def __init__(self, img, threshold=0):
        self.img = img
        self.threshold = threshold

    def get_img(self):
        return self.img

    def blur(self):
        pix_blur = (6, 6)
        k = np.ones(pix_blur) / float(pix_blur[0] * pix_blur[1])
        self.img = ndimage.convolve(self.img, k, mode='mirror')

        return self

    def _partial_thresholding(self, procnum, return_dict, h_from, lines_no):
        """Do thresholding only on certain rows.intensity

        The function is appropriate for multiprocessing where each process
        can thresholds pixels on different lines.

        :param procnum: number of the process
        :type procnum: int
        :param return_dict: Shared dict among processes
        :param h_from: starting index of a row
        :param lines_no: amount of lines which will be processed by this process
        :type lines_no: int
        :param t: threshold
        :type t: int
        :param self.img: 2D array containing pixels
        :type self.img: numpy.Array
        """
        intensity_array = np.zeros((lines_no, self.img.shape[1]))
        for h in range(0, lines_no):
            for w in range(0, intensity_array.shape[1]):  # size[1] = width
                intensity = self.img[h + h_from][w]
                if (intensity <= self.threshold):
                    x = 0
                else:
                    x = 255
                intensity_array[h][w] = x
        return_dict[procnum] = intensity_array

    def threshold_img(self):
        """Applies given threshold on the given self.img

        :param t: threshold
        :type t: int
        :type self.img: ImageFile
        :return: self.img on which threshold was applied
        :rtype: ImageFile
        """
        # let's devide height to 4 parts
        height_total = self.img.shape[0]
        quarter = int(height_total / 4)

        h_from = 0
        is_total_height = 0
        return_dict = Manager().dict()
        jobs = []

        # let's create 4 processes
        for i in range(0, 4):
            # if creating the last process, give it the rest of the lines
            if i == 3:
                quarter = height_total - is_total_height
            is_total_height += quarter
            p = Process(target=self._partial_thresholding,
                        args=(i, return_dict, h_from, quarter))
            h_from += quarter
            jobs.append(p)
            p.start()

        # join processes
        for p in jobs:
            p.join()

        # let's put data together
        self.img = np.concatenate((return_dict[0], return_dict[1],
                                   return_dict[2], return_dict[3]), axis=0)

        return self

    def otsu(self):
        """Otsu method for finding threshold"""
        total = self.img.size
        current_max = 0
        threshold = 0
        sumT = 0
        sumF = 0
        sumB = 0

        hist = np.histogram(self.img, range(0, 257))
        for i in range(0, 256):
            sumT += i * hist[0][i]

        weightB = 0
        weightF = 0
        varBetween = 0
        meanB = 0
        meanF = 0

        for i in range(0, 256):
            weightB += hist[0][i]
            weightF = total - weightB
            if weightF == 0:
                break
            sumB += i * hist[0][i]
            sumF = sumT - sumB

            if weightB != 0:
                meanB = sumB / weightB
            if weightF != 0:
                meanF = sumF / weightF
            varBetween = weightB * weightF
            varBetween *= (meanB - meanF) ** 2
            if varBetween > current_max:
                current_max = varBetween
                threshold = i

        self.threshold = threshold - 15
        return self
