import calfem.geometry as cfg
import calfem.mesh as cfm
from constants import left_element, bottomright_element, top_element, half_circle, circle1, circle2, circle3, upper_length, lower_length, height, thickness, radius, circle_distance

class model():

    # --Defining geometry--
    def geometry():
        g = cfg.Geometry()
        g.point([0.0,0.0])   # point 0
        g.point([lower_length, 0.0])  # point 1
        g.point([upper_length, height]) # point 2
        g.point([0.0, height]) # point 3

        g.spline([0,1], marker = bottomright_element) # line 0
        g.spline([1,2], marker = bottomright_element) # line 1
        g.spline([2,3], marker = top_element) # line 2

        # Half circle and lower boundary
        g.point([0.0, height/2]) # point 4
        g.point([0.0, height/2 - radius]) # point 5
        g.point([radius, height/2]) # point 6
        g.point([0.0, height/2 + radius]) # point 7

        g.spline([3,7], marker=left_element) # line 3
        g.spline([5,0], marker=left_element) # line 4
        g.circle([5,4,6], marker=half_circle) # line 5
        g.circle([6,4,7], marker=half_circle) # line 6

        # circle 1
        g.point([circle_distance, height/2])         # point 8
        g.point([circle_distance - radius, height/2])  # point 9
        g.point([circle_distance, height/2 - radius]) # point 10
        g.point([circle_distance + radius, height/2])  # point 11
        g.point([circle_distance, height/2 + radius]) # point 12

        g.circle([9,8,10], marker=circle1)   # line 7
        g.circle([10,8,11], marker=circle1)  # line 8
        g.circle([11,8,12], marker=circle1)  # line 9
        g.circle([12,8,9], marker=circle1)   # line 10

        # circle 2
        g.point([2*circle_distance, height/2])  # point 13
        g.point([2*circle_distance - radius, height/2]) # point 14
        g.point([2*circle_distance, height/2 - radius]) # point 15
        g.point([2*circle_distance + radius, height/2]) # point 16
        g.point([2*circle_distance, height/2 + radius]) # point 17

        g.circle([14,13,15], marker=circle2) # line 11
        g.circle([15,13,16], marker=circle2) # line 12
        g.circle([16,13,17], marker=circle2) # line 13
        g.circle([17,13,14], marker=circle2) # line 14

        # circle 3
        g.point([3*circle_distance, height/2])  # point 18
        g.point([3*circle_distance - radius, height/2]) # point 19
        g.point([3*circle_distance, height/2 - radius]) # point 20
        g.point([3*circle_distance + radius, height/2]) # point 21
        g.point([3*circle_distance, height/2 + radius]) # point 22

        g.circle([19,18,20], marker=circle3) # line 15
        g.circle([20,18,21], marker=circle3) # line 16
        g.circle([21,18,22], marker=circle3) # line 17
        g.circle([22,18,19], marker=circle3) # line 18

        g.surface([0,1,2,3,6,5,4], [[7,8,9,10],[11,12,13,14],[15,16,17,18]])

        return g

    # --Create the mesh--
    def create_mesh(geometry, elType, dofsPerNode, elSizeFactor):
        mesh = cfm.GmshMesh(geometry, return_boundary_elements=True)

        mesh.elType = elType # Type of meshing element.
        mesh.dofsPerNode = dofsPerNode # Degress of freedom per node.
        mesh.elSizeFactor = elSizeFactor # Element size Factor

        return mesh