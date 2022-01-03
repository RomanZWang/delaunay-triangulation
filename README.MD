# Installation
```
conda create --name cgeofinal python=3.7
conda activate cgeofinal
pip install numpy pandas matplotlib ipykernel
python -m ipykernel install --user --name=cgeofinal
cd ./delaunay-triangulation
jupyter notebook
# Open notebook
# Select kernel cgeofinal
# Run code blocks
```
# Background
This project idea is based on the [final project of CS274](https://people.eecs.berkeley.edu/~jrs/274/proj.html) from Berkeley. File is implementation of the incremental-insertion Delaunay Triangulation algorithm found in *Primitives for the Manipulation of General Subdivisions and the Computation of Voronoi Diagrams* by Guibas and Stolfi; and a divide-and-conquer algorithm found in *A faster divide-and-conquer algorithm for constructing delaunay triangulations* by Dwyer.

# References
[Quickselect algorithm](https://www.techiedelight.com/quickselect-algorithm/)
[Triangulation data structure](https://github.com/V-Hill/delaunay-triangulation/blob/main/delaunay_triangulation/triangulation_core/edge_topology.py)
[Reference code for clarifications for ambiguities in original papers](https://github.com/V-Hill/delaunay-triangulation/)