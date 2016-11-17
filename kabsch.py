"""
points_one and points_two are each a list of tuples. Each tuple is of the form (x, y, z). They have the same number of points.

We want to minimize the RMSD.

Returns a list of tuples, that is points_one moved so as to minimize the RMSD between itself and points_two.

Credit to this: http://nghiaho.com/?page_id=671 for explaining kabsch's algorithm to me.
"""
import numpy as np
def kabsch(points_one, points_two):
    centroid_one = (0.0, 0.0, 0.0)
    for x, y, z in points_one:
        centroid_one[0] += x
        centroid_one[1] += y
        centroid_one[2] += z

    centroid_one[0] = centroid_one[0]/len(points_one)
    centroid_one[1] = centroid_one[1]/len(points_one)
    centroid_one[2] = centroid_one[2]/len(points_one)

    centroid_two = (0.0, 0.0, 0.0)
    for x, y, z in points_two:
        centroid_two[0] += x
        centroid_two[1] += y
        centroid_two[2] += z

    centroid_two[0] = centroid_two[0]/len(points_two)
    centroid_two[1] = centroid_two[1]/len(points_two)
    centroid_two[2] = centroid_two[2]/len(points_two)

    #now, translate each set of points so the centroids are at the origin
    translate_one = np.array([(x - centroid_one[0], y - centroid_one[1], z - centroid_one[2]) for x,y,z in points_one])
    translate_two = np.array([(x - centroid_two[0], y - centroid_two[1], z - centroid_two[2]) for x,y,z in points_two])
    covariance = np.dot(translate_one.transpose(), translate_two)
    u, s, v = np.linalg.svd(covariance)
    rotation_matrix = np.dot(v, u.transpose())
    if np.linalg.det(rotation_matrix) < 0:
        rotation_matrix = np.dot(rotation_matrix, np.array([(1, 0, 0), (0, 1, 0), (0, 0, -1)]))
    translation_matrix = np.dot(-1*rotation_matrix, centroid_one) + centroid_two
    
