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
n = 9
#points = np.random.randint(2, 18, size=(n,2))
#print("Points whose voronoi diagram to be found")
#print(points)
points = []
points.append([0,8])
points.append([2,0])
points.append([5,20])
points.append([7,4])
points.append([9,1])
points.append([11,10])
points.append([13,16])
points.append([6,1])
# points.append([4,0])
points.append([17,18])
# points.append([5,6])
# points.append([6,5])
# points.append([6,6])
# points.append([5,5])
# add 4 distant dummy points
#points = np.append(points, [[10,0], [0,10], [0,0], [10,10]], axis = 0)
#points = np.append(points, [[999,999], [-999,999], [999,-999], [-999,-999]], axis = 0)




############################## Driver Function ############################################################
def centroid_func(points):
  # compute Voronoi tesselation
  vor = Voronoi(points)


  # plot
  myfig = voronoi_plot_2d(vor, show_vertices=False, line_colors='black',line_width=2, line_alpha=1, point_size=10)

  # colorize
  for region in vor.regions:
      if not -1 in region:
          polygon = [vor.vertices[i] for i in region]
          plt.fill(*zip(*polygon),alpha=0.4)

  # fix the range of axes subhajit
  plt.xlim([-1,21]), plt.ylim([-1,21])

  #print(vertices)
  #print("Voronoi vertices of finite polygons")
  #print(vor.vertices)

  
  plt.show()
  #Convex hull of the points 
  hull = ConvexHull(points)
  #plt.plot(points[:,0], points[:,1], '*')


  #for simplex in hull.simplices:

    # plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

  #plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)

  #plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')

  #Find intersection where the line between two voronoi finite vertex and intfinite vertex intersect the the sides our square region

  def lexico(p1, p2):
    if ((p1[0]> p2[0]) or (p1[0] == p2[0] and p1[1] > p2[1])):
      return 1
    else:
      return 0

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
  total_region =len(regions)
  len_finite_vor = len(vor.vertices)

  def distance(a,b):
    p1 = (a[0]-b[0])*(a[0]-b[0])
    p2 = (a[1]-b[1])*(a[1]-b[1])
    p3 = p1 + p2
    p4 = p3**0.5
    return p4
  def is_between(a,c,b):
      return distance(a,c) + distance(c,b) == distance(a,b)



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

  a = [-1, -1]
  b = [21, -1]
  c = [21, 21]
  d = [-1, 21]
  w1=[[-1,3],[4,-1]]
  #print("vor",detect_corner(w1,a,b,c,d))
  s1 = Segment(a, b)
  s2 = Segment(b, c)
  s3 = Segment(c, d)
  s4 = Segment(d, a)
  total_vertices = len(vertices)

  cor_vertices = []

  t1=[]
  count = 0
  for i in range(total_region):
    temp11 = 0  ####region i is a finite region#####
    for l in range(len(regions[i])):
      if(regions[i][l] >= len_finite_vor):
        temp11 = 1 ####region i is infinite region#####
    len_of_i = len(regions[i]) ###no. of vertices in region i ####
  #  t1 = []
    for j in range(len_of_i-1): 
      k=j+1
      if ((regions[i][j] < len_finite_vor and regions[i][k] >= len_finite_vor)):  ### vornoi line is infinite####
        v1 = vertices[regions[i][j]]
        v2 = vertices[regions[i][k]]
        t = find_intersection(v1, v2, s1, s2, s3, s4)
        vertices[regions[i][k]] = t
        t1.append(vertices[regions[i][k]])
        count = count + 1
      if ((regions[i][j] >= len_finite_vor and regions[i][k] < len_finite_vor)):
        v3 = vertices[regions[i][j]]
        v4 = vertices[regions[i][k]]
        t = find_intersection(v4, v3, s1, s2, s3, s4)
        vertices[regions[i][j]] = t
        t1.append(vertices[regions[i][j]])
        count = count + 1
    if ((regions[i][len_of_i-1] < len_finite_vor and regions[i][0] >= len_finite_vor)):
        v1 = vertices[regions[i][len_of_i-1]]
        v2 = vertices[regions[i][0]]
        t = find_intersection(v1, v2, s1, s2, s3, s4)
        vertices[regions[i][0]] = t
        t1.append(vertices[regions[i][0]])
        count = count + 1
    if ((regions[i][len_of_i-1] >= len_finite_vor and regions[i][0] < len_finite_vor)):
        v3 = vertices[regions[i][len_of_i-1]]
        v4 = vertices[regions[i][0]]
        t = find_intersection(v4, v3, s1, s2, s3, s4)
        vertices[regions[i][len_of_i-1]] = t
        t1.append(vertices[regions[i][len_of_i-1]])
        count = count + 1
    if temp11 == 1:
      if detect_corner(t1[count-2], t1[count-1], a, b, c, d) == a:
        # region[i].append(total_vertices)
        cor_vertices.append(a)
        #print("subhahjit", cor_vertices, a)
      if detect_corner(t1[count-2], t1[count-1], a, b, c, d) == b:
        # region[i].append(total_vertices+1)
        cor_vertices.append(b)
      if detect_corner(t1[count-2], t1[count-1], a, b, c, d) == c:
        # region[i].append(total_vertices+2)
        cor_vertices.append(c)
      if detect_corner(t1[count-2], t1[count-1], a, b, c, d) == d:
        # region[i].append(total_vertices+3)
        cor_vertices.append(d)
      if detect_corner(t1[count-2], t1[count-1], a, b, c, d) == []:
        #print("subhahjit p")
        # region[i].append(total_vertices+3)
        cor_vertices.append([])
    else:
      cor_vertices.append([])

  #print("Modified Vertices are:", vertices, sep = "\n")
  def verify(p,a,b,c,d): ########## whether point p is inside the square a-b-c-d or not ######################
    p0 = p[0]
    p1 = p[1]
    if (p0 < a[0] or p0 > b[0] or p1 < a[1] or p1 > c[1]):
      return 1
    else:
      return 0  
  def side_intersection(seg1, s1, s2, s3, s4):
    if seg1.intersection(s1) != []: 
      point = seg1.intersection(s1)
    elif seg1.intersection(s2) != []: 
      point= seg1.intersection(s2)
    elif seg1.intersection(s3) != []: 
      point= seg1.intersection(s3)
    elif seg1.intersection(s4) != []: 
      point = seg1.intersection(s4)
    else:
      point =  []
    point1 =[]
    #print(point[0][0])
    #print(point[0][1])
    point1.append(point[0][0])
    point1.append(point[0][1])
    return point1



  for i in range(total_region):
    temp1_list = []
    len_of_i = len(regions[i])
    if cor_vertices[i] != []:   ####### For infinte regions ############################
      for j in range(len_of_i-1):
        if verify(vertices[regions[i][j]],a,b,c,d) == 1:
          p00 = vertices[regions[i][j]]
          if j == 0:
            q00 = vertices[regions[i][1]]
            q11 = vertices[regions[i][len_of_i-1]]
            seg0 = Segment(q00, p00)
            seg1 = Segment(q11, p00)
            r00 = side_intersection(seg0, s1, s2, s3, s4)
            r11 = side_intersection(seg1, s1, s2, s3, s4)
            if r00 != []:
              t00 = r00
            else:
              t00 = r11
          else:
            q00 = vertices[regions[i][j+1]]
            q11 = vertices[regions[i][j-1]]
            seg0 = Segment(q00, p00)
            seg1 = Segment(q11, p00)
            r00 = side_intersection(seg0, s1, s2, s3, s4)
            r11 = side_intersection(seg1, s1, s2, s3, s4)
            if r00 != []:
              t00 = r00
            else:
              t00 = r11
          temp1_list.append(t00)  
        else:  
          temp1_list.append(vertices[regions[i][j]])
      if (regions[i][j] >= len_finite_vor and regions[i][j+1] >= len_finite_vor):
        temp1_list.append(cor_vertices[i])
      if verify(vertices[regions[i][len_of_i-1]],a,b,c,d) == 1:
        q00 = vertices[regions[i][0]]
        q11 = vertices[regions[i][len_of_i-2]]
        seg0 = Segment(q00, p00)
        seg1 = Segment(q11, p00)
        r00 = side_intersection(seg0, s1, s2, s3, s4)
        r11 = side_intersection(seg1, s1, s2, s3, s4)
        if r00 != []:
          t00 = r00
        else:
          t00 = r11
        temp1_list.append(t00)
      else:
        temp1_list.append(vertices[regions[i][len_of_i-1]])
      if (regions[i][len_of_i-1] >= len_finite_vor and regions[i][0] >= len_finite_vor):
        temp1_list.append(cor_vertices[i])
    else:                       ############## for finite regions ######################
      for j in range(len_of_i):
        if verify(vertices[regions[i][j]],a,b,c,d) == 1: ###### whether the point is outside the region ######
          p00 = vertices[regions[i][j]]
          if j == 0:
            q00 = vertices[regions[i][1]]
            q11 = vertices[regions[i][len_of_i-1]]
            seg0 = Segment(q00, p00)
            seg1 = Segment(q11, p00)
            r00 = side_intersection(seg0, s1, s2, s3, s4)
            r11 = side_intersection(seg1, s1, s2, s3, s4)
            t00 = detect_corner(r00, r11, a,b,c,d)
          elif j == len_of_i-1:
            q00 = vertices[regions[i][0]]
            q11 = vertices[regions[i][len_of_i-2]]
            seg0 = Segment(q00, p00)
            seg1 = Segment(q11, p00)
            r00 = side_intersection(seg0, s1, s2, s3, s4)
            r11 = side_intersection(seg1, s1, s2, s3, s4)
            t00 = detect_corner(r00, r11, a,b,c,d)
          else:
            q00 = vertices[regions[i][j+1]]
            q11 = vertices[regions[i][j-1]]
            seg0 = Segment(q00, p00)
            seg1 = Segment(q11, p00)
            r00 = side_intersection(seg0, s1, s2, s3, s4)
            r11 = side_intersection(seg1, s1, s2, s3, s4)
            t00 = detect_corner(r00, r11, a,b,c,d)
          temp1_list.append(r11)
          if t00 != []:
            temp1_list.append(t00)
          temp1_list.append(r00)
        else:             ########## if the point is inside the region #######
          temp1_list.append(vertices[regions[i][j]])

    
      
    #print(i, temp1_list, sep ="\n")
    p = Polygon(temp1_list)
    centroid = p.centroid
    points[i] = [centroid.x, centroid.y]
  #print("Centroid", points)  
  pdfFile.savefig(myfig)

for i in range(3):
  centroid_func(points)
  
  # files.download('saved.jpeg') 

pdfFile.close()
