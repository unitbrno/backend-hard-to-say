#!/usr/bin/env python3
from fractions import Fraction
from scipy import ndimage
import math


class Detector(object):
    """Class Detector
    Detect bounded boxes
    """
    def __init__(self, img):
        self.img = img

    def get_bounded_boxes(self):
        """Get bounded boxes from thresholded image"""
        s = ndimage.generate_binary_structure(2, 2)
        labeled_arr, num_objects = ndimage.label(self.img, structure=s)

        dots = ndimage.find_objects(labeled_arr)

        bboxes = []
        for i, j in enumerate(dots):

            if (dots[i][0].start != 0 and
                    dots[i][1].start != 0 and
                    dots[i][0].stop < self.img.shape[0] and
                    dots[i][1].stop < self.img.shape[1] and
                    labeled_arr[j].shape[0] > 10 and
                    labeled_arr[j].shape[1] > 10):

                garbage_arr, num_garbage = ndimage.label(labeled_arr[j],
                                                         structure=s)
                garbage = ndimage.find_objects(garbage_arr)

                if len(garbage) > 1:
                    for k, l in enumerate(garbage):
                        if (not (garbage[k][0].start == 0 and
                                 garbage[k][1].start == 0 and
                                 garbage[k][0].stop == labeled_arr[j].shape[0] and
                                 garbage[k][1].stop == labeled_arr[j].shape[1])):
                            garbage_arr[l] = 0
                    bboxes.append(garbage_arr)
                else:
                    bboxes.append(labeled_arr[j])

        return bboxes

    def _get_border_from_left(self, row):
        for i in range(0, len(row)):
            if row[i] != 0:
                # return index of the column where the
                # color is not background color (black)
                return i

    def _get_border_from_right(self, row):
        length = len(row)
        for i in range(length - 1, -1, -1):
            if row[i] != 0:
                return i

    def _get_border_from_top(self, c, matrix):
        for i in range(0, len(matrix)):
            if matrix[i, c] != 0:
                return i

    def _get_border_from_bottom(self, c, matrix):
        length = len(matrix)
        for i in range(length - 1, -1, -1):
            if matrix[i, c] != 0:
                return i

    def _append_if_not_in(self, what, to):
        if what not in to:
            to.append(what)
        return to

    def convex_hull(self, bbox):
        borders = []
        r_length = len(bbox)
        c_length = len(bbox[0])  # all rows has the same length

        for c in range(0, c_length):
            r = self._get_border_from_top(c, bbox)
            coordinates = [r, c]
            borders = self._append_if_not_in(coordinates, borders)

            r = self._get_border_from_bottom(c, bbox)
            coordinates = [r, c]
            borders = self._append_if_not_in(coordinates, borders)

        for r in range(0, r_length):
            c = self._get_border_from_left(bbox[r])
            coordinates = [r, c]
            borders = self._append_if_not_in(coordinates, borders)

            c = self._get_border_from_right(bbox[r])
            coordinates = [r, c]
            borders = self._append_if_not_in(coordinates, borders)

        return borders

    def _get_len(self, p1, p2):
        """Get distance between 2 points"""
        return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

    def _get_lengths(self, edge_list):
        """Get [{'length' : length, 'points': [[x1, y1], [x2, y2]]}]"""
        return [dict([('length', self._get_len(i, j)), ('points', [i, j])]) for i in edge_list for j in edge_list]

    def max_length(self, edge_list):
        lengths = self._get_lengths(edge_list)

        max_len = 0
        max_points = []
        for d in lengths:
            if d['length'] > max_len:
                max_len = d['length']
                max_points = d['points']

        return round(max_len, 2), max_points

    def _get_line_eq(self, points):
        """ Compute line equation from given two points
        """
        #    A   x +    B   y +      C      = 0
        # (y1-y2)x + (x2-x1)y + (x1y2-x2y1) = 0
        # points = [[x1, y1], [x2, y2]]
        x1 = points[0] [0]
        y1 = points[0] [1]
        x2 = points[1] [0]
        y2 = points[1] [1]
        return {'a':(y1-y2), 'b':(x2-x1), 'c':(x1*y2-x2*y1)}

    def _is_below_line(self, line_eq, point):
        result = line_eq['a'] * point[0] + line_eq['b'] * point[1] + line_eq['c']
        if result > 0:
            return False
        else:
            return True

    def _is_normal(self, line_eq, normal_eq):
        """ True if lines are orthogonal
        """
        if (line_eq['a'] * normal_eq['a'] + line_eq['a'] * normal_eq['b'] + line_eq['b'] * normal_eq['a'] + line_eq['b'] * normal_eq['b']) == 0:
            return True
        else:
            return False

    def _get_normal(self, line_eq, point):
        lamb = (-line_eq['c'] - (line_eq['a']*point[0]) - (line_eq['b']*point[1])) / (line_eq['a']*line_eq['a'] + line_eq['b']*line_eq['b'])
        x = point[0] + lamb * line_eq['a']
        y = point[1] + lamb * line_eq['b']
        return self._get_line_eq([point,[x, y]])

    def _is_uninterrupted(self, box, equation, x1, y1, x2, y2):
        normal = self._draw_line(x1, y1, x2, y2)
        for point in normal:
            if box[point[1]][point[0]] == 0:
                return False
        return True

    def _draw_line(self, x0, y0, x1, y1):
        points = []
        rev = reversed
        if abs(y1 - y0) <= abs(x1 - x0):
            x0, y0, x1, y1 = y0, x0, y1, x1
            rev = lambda x: x
        if x1 < x0:
            x0, y0, x1, y1 = x1, y1, x0, y0
        leny = abs(y1 - y0)
        for i in range(leny + 1):
            points.append([*rev((round(Fraction(i, leny) * (x1 - x0)) + x0, (1 if y1 > y0 else -1) * i + y0))])
        return points

    def _split_along_line(self, line_eq, points):
        points_below = []
        points_above = []
        for point in points:
            if self._is_below_line(line_eq, point):
                points_below.append(point)
            else:
                points_above.append(point)
        return points_below, points_above

    def max_thickness(self, points, edge, box):
        max_value = 0
        value = None

        line_eq = self._get_line_eq(points)

        edge_below, edge_above = self._split_along_line(line_eq, edge)

        normal_eq = self._get_normal(line_eq, edge_above[int(len(edge_above)/2)])

        normal_above_left, normal_above_rigth = self._split_along_line(normal_eq, edge_above)
        normal_below_left, normal_below_rigth = self._split_along_line(normal_eq, edge_below)

        for a in normal_above_left:
            for b in normal_below_left:
                propsed_normal = self._get_line_eq([a, b])
                if self._is_normal(line_eq, propsed_normal):
                    if self._is_uninterrupted(box, propsed_normal, a[0], a[1], b[0], b[1]):
                        value = self._get_len(a, b)
                        if value > max_value:
                            max_value = value

        for a in normal_above_rigth:
            for b in normal_below_rigth:
                propsed_normal = self._get_line_eq([a, b])
                if self._is_normal(line_eq, propsed_normal):
                    if self._is_uninterrupted(box, propsed_normal, a[0], a[1], b[0], b[1]):
                        value = self._get_len(a, b)
                        if value > max_value:
                            max_value = value

        return max_value
