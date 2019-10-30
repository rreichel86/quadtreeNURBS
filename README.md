## Quadtree Decomposition of Heterogeneous Materials in NURBS-boundary Representation

Composite materials have been for decades widely used in civil engineering and related fields.
Some approaches for treating numerical models, which deal with different materials, involve the
employment of uniform meshes. This might be expensive in terms of computational resources.
High accuracy of the solution is required; therefore, a fine mesh is needed. At the same time, the
number of elements produced by this kind of approaches must be reduced. This dilemma could
be resolved by creating an automatic mesh that suits the geometry of the model, creating more
elements only where it is needed.
A novel meshing technique for composite materials is presented. The description of interfaces
is done with NURBS. This method makes use of the quadtree decomposition procedure, which
generates non-uniform interface-conforming meshes. This allows us to create an automatic nonuniform
mesh, based on the complexity of the geometry, and to maintain an exact representation
of the boundary.