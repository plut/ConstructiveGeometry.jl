using ConstructiveGeometry

t1 = Surface([[0,0,0],[2,0,0],[0,2,0],[0,0,2]],
	[[1,3,2],[1,4,3],[1,2,4],[2,3,4]])
t2 = Surface([[2,0,0],[0,0,0],[2,0,2],[2,2,0]],
	[[1,3,2],[1,4,3],[1,2,4],[2,3,4]])
println("""
\e[31;1mUnion of two tetrahedra:\e[m
$(mesh(t1 ∪ t2))
\e[31;1mIntersection of two tetrahedra:\e[m
$(mesh(t1 ∩ t2))
\e[31;1mDifference of two tetrahedra:\e[m
$(mesh(t1 \ t2))
\e[31;1mConvex hull of two tetrahedra:\e[m
$(mesh(hull(t1, t2)))
""")
# file `tetrahedra.png` contains images of these four constructions
