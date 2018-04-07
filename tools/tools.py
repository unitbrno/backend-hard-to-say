#!/usr/bin/env python 3
import cv2


def mask_thresh(orig_img_path, threshed_img_path):
    im = cv2.imread(orig_img_path)
    thresh = cv2.imread(threshed_img_path)

    # to grayscale
    im_gray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    thresh_gray = cv2.cvtColor(thresh, cv2.COLOR_BGR2GRAY)

    dst = cv2.addWeighted(im_gray, 0.9, thresh_gray, 0.1, 0)

    cv2.imshow('masked', dst)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
