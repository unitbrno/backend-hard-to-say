#!/usr/bin/env python3
"""
File: main.py
Authors: Martin Kopec <xkopec42@gmail.com>
         Maros Kopec <xkopec44@vutbr.cz>
         Patrik Segedy <xseged00@vutbr.cz>
         Tomas Sykora <xsykor25>
"""

from AABBlib import detection
from AABBlib import threshold
import argparse
import csv
import cv2


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--image-path", required=True,
                        help="Path to an image to be processed")
    parser.add_argument("--csv-path", required=True,
                        help="Path where csv file will be stored")
    parser.add_argument("--resize", required=False, default=100,
                        help="Percentage to scale picture down")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    img = cv2.imread(args.image_path, 0)

    if args.resize != 100:
        resize_percentage = int(args.resize) / 100
        img = cv2.resize(img,
                        (int(img.shape[1] * resize_percentage),
                         int(img.shape[0] * resize_percentage)))

    thresh = threshold.Threshold(img)
    thresh.blur()
    thresh.otsu()
    thresh.threshold_img()

    detector = detection.Detector(thresh.get_img())

    bboxes = detector.get_bounded_boxes()

    box_width = []
    box_height = []
    max_thick = []
    max_len = []
    max_points = []
    edge_list = []

    # calculate a coefficient for changig lengths
    # based on resize of input picture
    k = int((1 / int(args.resize)) * 100)

    for bbox in bboxes:
        edge_list.append(detector.convex_hull(bbox))
        box_height.append(bbox.shape[0] * k)
        box_width.append(bbox.shape[1] * k)

    for edges in edge_list:
        max_l, max_p = detector.max_length(edges)
        max_len.append(max_l * k)
        max_points.append(max_p)

    for edges, bbox, point in zip(edge_list, bboxes, max_points):
        max_thick.append(round(detector.max_thickness(point, edges, bbox), 2))

    zipped = zip(range(1, len(max_len) + 1),
                 box_width, box_height, max_len, max_thick)

    with open(args.csv_path, 'w') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(['Part #', 'Width', 'Height',
                         'Max Length', 'Thickness'])
        writer.writerows(zipped)

    cv2.imwrite('thresh.tif', thresh.get_img())    # dump threshed img
