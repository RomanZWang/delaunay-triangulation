from math import sqrt
import numpy as np

from random import randint

def swap(A, i, j):
    temp = A[i]
    A[i] = A[j]
    A[j] = temp

def partition(A, left, right, pIndex):

    pivot = A[pIndex]

    swap(A, pIndex, right)

    pIndex = left

    for i in range(left, right):
        if A[i] <= pivot:
            swap(A, i, pIndex)
            pIndex = pIndex + 1

    swap(A, pIndex, right)

    return pIndex

def quickSelectHelper(A, left, right, k):

    if left == right:
        return A[left]

    pIndex = randint(left, right)

    pIndex = partition(A, left, right, pIndex)

    if k == pIndex:
        return A[k]

    elif k < pIndex:
        return quickSelectHelper(A, left, pIndex - 1, k)

    else:
        return quickSelectHelper(A, pIndex + 1, right, k)

def quickSelect(A, k):
    return quickSelectHelper(A, 0, len(A)-1, k)

def lexigraphic_sort(points):
    points_sorted = sorted(points, key=lambda k: [k[0], k[1]])
    return points_sorted

# Sorts based on y then x coordinates
def lexigraphic_swap_sort(points):
    points_sorted = sorted(points, key=lambda k: [k[1], k[0]])
    return points_sorted

def in_circle(a, b, c, d):

    result = np.array([
        [a[0], a[1], a[0]**2+a[1]**2, 1],
        [b[0], b[1], b[0]**2+b[1]**2, 1],
        [c[0], c[1], c[0]**2+c[1]**2, 1],
        [d[0], d[1], d[0]**2+d[1]**2, 1]
    ])

    return np.linalg.det(result) >= 0

def is_right_turn(p, q, r):
    if len(p) != 2 or len(q) != 2 or len(r) != 2 :
        raise ValueError("Not enough coords for an input")

    px, py = p
    qx, qy = q
    rx, ry = r

    checker_array = np.array(
        [[1, px, py],
        [1, qx, qy],
        [1, rx, ry]]
    )

    return np.linalg.det(checker_array) > 0

def is_left_turn(p, q, r):
    if len(p) != 2 or len(q) != 2 or len(r) != 2 :
        raise ValueError("Not enough coords for an input")

    px, py = p
    qx, qy = q
    rx, ry = r

    checker_array = np.array(
        [[1, px, py],
        [1, qx, qy],
        [1, rx, ry]]
    )

    return np.linalg.det(checker_array) < 0
