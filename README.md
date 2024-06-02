Imagine the points as mobile robots. They are autonomous and deployed on a bounded region (for exapmle: rectangle, square)
<br>
They need to partition the region. We are using voronoi partitioning in such a way that in each round, every robot computes the voronoi diagram of the region and moves to the C.G of the corresponding cell. 