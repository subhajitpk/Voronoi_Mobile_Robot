#########SUBHAJIT##############
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from numpy import sqrt
from sympy import Line, Point, Segment
from shapely.geometry import Point, Polygon
#from google.colab import files
from matplotlib.backends.backend_pdf import PdfFile, PdfPages



pdfFile = PdfPages("output.pdf")

# make up data points
n = 10
#points = np.random.randint(2, 18, size=(n,2))
#print("Points whose voronoi diagram to be found")
#print(points)
points = []
# points.append([0,8])
# points.append([2,0])
# points.append([5,20])
# points.append([7,4])
# points.append([9,1])
# points.append([11,10])
# points.append([13,16])
# points.append([6,1])
# points.append([4,0])
# points.append([17,18])

# points.append([5,6])
# points.append([6,5])
# points.append([6,6])
# points.append([5,5])

# points.append([17,1])
# points.append([3,18])
# points.append([15,18])
# points.append([10,11])
# points.append([12,19])
# points.append([6,15])

# points.append([0,2])
# points.append([1,11])
# points.append([15,0])
# points.append([16,8])
# points.append([18,4])
# points.append([19,13])
# points.append([20,15])
# points.append([11,8])
# points.append([17,5])
# points.append([20,4])

# points.append([1,16])
# points.append([8,18])
# points.append([7,16])
# points.append([11,3])
# points.append([15,10])
# points.append([13,15])




points.append([19,8])
points.append([17,5])
points.append([20,4])

points.append([19,16])
points.append([19,18])
points.append([17,16])
points.append([18,3])
points.append([18,10])
points.append([18,15])
a = [-1, -1]
b = [21, -1]
c = [21, 21]
d = [-1, 21]
s1 = Segment(a, b)
s2 = Segment(b, c)
s3 = Segment(c, d)
s4 = Segment(d, a)
w1=[[-1,3],[4,-1]]
#print("vor",detect_corner(w1,a,b,c,d))




############################## Driver Function ############################################################
def centroid_func(points):
  # compute Voronoi diagrams for points (regions can be infinite)
  vor = Voronoi(points)

  # plot
  myfig = voronoi_plot_2d(vor, show_vertices=False, line_colors='black',line_width=2, line_alpha=1, point_size=10)

  # colorize
  for region in vor.regions:
      if not -1 in region:
          polygon = [vor.vertices[i] for i in region]
          plt.fill(*zip(*polygon),alpha=0.4)

  # fix the range of axes to show in the output
  plt.xlim([-1,21]), plt.ylim([-1,21])
  # plt.show()
  #print(vertices)
  #print("Voronoi vertices of finite polygons")
  #print(vor.vertices)


  #Find intersection where the line between two voronoi finite vertex and infinite vertex intersect the the sides our square region
  
  #Function to compute lexicographic order between two points. 1 means p1 is bigger, 0 means p2 is bigger 
  def lexico(p1, p2):
    if ((p1[0]> p2[0]) or (p1[0] == p2[0] and p1[1] > p2[1])):
      return 1
    else:
      return 0
    
  #### Whether point p is inside the square a-b-c-d or not. 1 means outside and 0 means inside
  def verify(p,a,b,c,d): 
    p0 = p[0]
    p1 = p[1]
    if (p0 < a[0] or p0 > b[0] or p1 < a[1] or p1 > c[1]):
      return 1 
    else:
      return 0  
    
  #### Finding the intersection between the line segment between two point
  #### p1, p2 and the sides of the square
  #### return the intersection point on the side
  def find_intersection(p1, p2, s1, s2, s3, s4):
    l1 = Line(p1, p2)
    coord = []
    if l1.intersection(s1) != []:
      coord.append(l1.intersection(s1)[0][0])
      coord.append(l1.intersection(s1)[0][1])
    if l1.intersection(s2) != []:
      coord.append(l1.intersection(s2)[0][0])
      coord.append(l1.intersection(s2)[0][1])
    if l1.intersection(s3) != []:
      coord.append(l1.intersection(s3)[0][0])
      coord.append(l1.intersection(s3)[0][1])
    if l1.intersection(s4) != []:
      coord.append(l1.intersection(s4)[0][0])
      coord.append(l1.intersection(s4)[0][1])
    #print(coord,"shsi")
    r1 = []
    r2 = []
   # if r1 != []:
    r1.append(coord[0])
    r1.append(coord[1])
   # if r2 != []:
    r2.append(coord[2])
    r2.append(coord[3])
    #print(r1, r2)
    if lexico(p1, p2) == 1:
      if lexico(p1, r1) == 1:
        t = r1
      else:
        t = r2
    if lexico(p1, p2) == 0:
      if lexico(p1, r1) == 0:
        t = r1
      else:
        t = r2
    return t

  #Finite regions of the voronoi polygons

  def voronoi_finite_polygons_2d(vor, radius=None):
      """
      Reconstruct infinite voronoi regions in a 2D diagram to finite
      regions.

      Parameters
      ----------
      vor : Voronoi
          Input diagram
      radius : float, optional
          Distance to 'points at infinity'.

      Returns
      -------
      regions : list of tuples
          Indices of vertices in each revised Voronoi regions.
      vertices : list of tuples
          Coordinates for revised Voronoi vertices. Same as coordinates
          of input vertices, with 'points at infinity' appended to the
          end.

      """

      if vor.points.shape[1] != 2:
          raise ValueError("Requires 2D input")

      new_regions = []
      new_vertices = vor.vertices.tolist()

      center = vor.points.mean(axis=0)
      if radius is None:
          radius = vor.points.ptp().max()

      # Construct a map containing all ridges for a given point
      all_ridges = {}
      for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
          all_ridges.setdefault(p1, []).append((p2, v1, v2))
          all_ridges.setdefault(p2, []).append((p1, v1, v2))

      # Reconstruct infinite regions
      for p1, region in enumerate(vor.point_region):
          vertices = vor.regions[region]

          if all(v >= 0 for v in vertices):
              # finite region
              new_regions.append(vertices)
              continue

          # reconstruct a non-finite region
          ridges = all_ridges[p1]
          new_region = [v for v in vertices if v >= 0]

          for p2, v1, v2 in ridges:
              if v2 < 0:
                  v1, v2 = v2, v1
              if v1 >= 0:
                  # finite ridge: already in the region
                  continue

              # Compute the missing endpoint of an infinite ridge

              t = vor.points[p2] - vor.points[p1] # tangent
              t /= np.linalg.norm(t)
              n = np.array([-t[1], t[0]])  # normal

              midpoint = vor.points[[p1, p2]].mean(axis=0)
              direction = np.sign(np.dot(midpoint - center, n)) * n
              far_point = vor.vertices[v2] + direction * radius

              new_region.append(len(new_vertices))
              new_vertices.append(far_point.tolist())

          # sort region counterclockwise
          vs = np.asarray([new_vertices[v] for v in new_region])
          c = vs.mean(axis=0)
          angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
          new_region = np.array(new_region)[np.argsort(angles)]

          # finish
          new_regions.append(new_region.tolist())

      return new_regions, np.asarray(new_vertices)

 
  regions, vertices = voronoi_finite_polygons_2d(vor)
  #print("Regions are:", regions)
  #print("All Vertices are:", vertices, sep = "\n")
  total_region =len(regions)          ### number of vornoi regions
  len_finite_vor = len(vor.vertices)  #### number of finite vertices of the voronoi diagram
  total_vertices = len(vertices)      #### total vertices of the vornoi diagram
  
  #### Eucledian distance between two points a and b
  def distance(a,b):
    p1 = (a[0]-b[0])*(a[0]-b[0])
    p2 = (a[1]-b[1])*(a[1]-b[1])
    p3 = p1 + p2
    p4 = p3**0.5
    return p4
  
  ##### Whether the point c is in the line segment between a and b
  ##### return True if c is in the line segment ab
  def is_between(a,c,b):
      return distance(a,c) + distance(c,b) == distance(a,b)


  ##### Return the corner point of the square, if t11 and t12 lies on the adjacent sides of the square
  def detect_corner(t11, t12, a, b, c, d):
    if ((is_between(a, t11, b) and is_between(b, t12, c)) or (is_between(a, t12, b) and is_between(b, t11, c))):
      return b
    if ((is_between(b, t11, c) and is_between(c, t12, d)) or (is_between(b, t12, c) and is_between(c, t11, d))):
      return c
    if ((is_between(c, t11, d) and is_between(d, t12, a)) or (is_between(c, t12, d) and is_between(d, t11, a))):
      return d
    if ((is_between(d, t11, a) and is_between(a, t12, b)) or (is_between(d, t12, a) and is_between(a, t11, b))):
      return a
    else:
      return []
  

  new_region = []
  for i in range(total_region):
     new_region.append([])           #### Contain all the modified vertices for the region i
  #### Modifying all the vertices.
  for i in range(total_region):
     len_of_reg_i = len(regions[i])
     for j in range(len_of_reg_i):
        p = vertices[regions[i][j]]  
        k = (j-1) % len_of_reg_i
        q = vertices[regions[i][k]]   #### Previous point of p
        l = (j+1) % len_of_reg_i
        r = vertices[regions[i][l]]   #### Next point of p 
        if verify(p,a,b,c,d) == 0:    #### p is inside the square
           new_region[i].append(p)    #### no modification needed
        else:                         #### p is outside the square
           if verify(q,a,b,c,d) == 0: #### q is inside the square
              new_region[i].append(find_intersection(q,p,s1,s2,s3,s4))  #### Find the intersection of p and q with the side 
           if verify(r,a,b,c,d) == 0:   ### r is inside the square
              new_region[i].append(find_intersection(r,p,s1,s2,s3,s4))  #### Find the intersection of p and r with the side
            ### otherwise we ignore the point
  
  
  
  #### inserting the corners in the polygon
  new_finite_region = []
  for i in range(total_region):
     new_finite_region.append([])
  for i in range(total_region):
     len_of_finite_region_i = len(new_region[i])
     for j in range(len_of_finite_region_i):
        p = new_region[i][j]
        k = (j-1) % len_of_finite_region_i
        q = new_region[i][k]
        r = detect_corner(p,q,a,b,c,d)
        if r != []:
           new_finite_region[i].append(r)
        new_finite_region[i].append(p)
         



  # for i in range(total_vertices):
  #    print(i, '=', vertices[i])
  # for i in range(len(new_region)):
  #    print(i,'=', new_region[i])
  # for i in range(len(new_finite_region)):
  #    print(i,'=', new_finite_region[i])
  
  for i in range(len(new_finite_region)):
     p =Polygon(new_finite_region[i])
     centroid = p.centroid
     points[i] = [centroid.x, centroid.y]   
  pdfFile.savefig(myfig)  
    
for i in range(3):
  
  centroid_func(points) 

pdfFile.close()
