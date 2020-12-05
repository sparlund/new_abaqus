### new_abaqus
new_abaqus is a hobby project for creating an FE-solver reading an FE-mesh following the syntax of the commercial software abaqus. 

The name is a joke from the character YinYang in Silicon Valley, who has the idea for "new netflix".

### Example
The program takes 1 argument, an input file containing the mesh and load case definition. On linux:
```
./new_abaqus example_runfiles/ex3.inp
```
This will save an output file of the results called `ex3.vtk` and a logfile called `ex3.log`. The results can be viewed in the open source post processor ParaView.

This example is in 2D, and contains ~800 elements and ~900 nodes, with a mix of quad and trias.  Result comparison:

<img src="example_runfiles/ex3.png" width="75%"/>

| FE-solver      | Load node deflection (red arrow in figure above!) |
| ----------- | ----------- |
| abaqus (c)      | 0.0348       |
| new_abaqus   | 0.0346        |


### Element types available
CPS3, CPS4, C3D10, C3D8
<img src="elements_available.png" width="75%"/>


### To-do & Features implemented
- [x] Implement logic and structure for reading abaqus input files
  - [ ] Disregard unused nodes
- [x] Set up classes and functions for nodes, elements, properties and materials
- [ ] Implement logic for different elements
  - [x] 2D tria (S2)
  - [x] 2D quadrilateral (CPS4)
  - [ ] 3D second order tetra (C3D10)
  - [X] 3D hexahedron (C3D8)
  - [ ] 3D second order hexahedron (C3D20)
- [x] Assemble mass matrix  
- [x] Assemble stiffness matrix
  - [x] Modify stiffness matrix and load vector to account for boundary conditions
- [x] Add support for more keywords
  - [x] *BOUNDARY
  - [x] *CLOAD
  - [X] *MATERIAL (=mid)
  - [X] *SECTION (=pid)
  - [ ] *STEP
  - [ ] *STATIC
  - [ ] *EIGENFREQUENCY
  - [ ] *OUTPUT
- [x] Point load keyword and assemble global load vector
- [x] Solve Ka=f for linear problems
- [ ] Calculate scalar values on elements
  - [ ] Stresses and strains
  - [ ] von Mises stress
- [ ] Solve eigen value problem
  - [ ] Compute mass matrix for elements
  - [ ] Calculate mass matrix in element construction
- [ ] Export results to VTK format to view results in ParaView
  - [x] Nodal displacement
  - [ ] Stresses and strains  
- [x] Re-direct output to a log file for debugging  
  - [x] Print timing for each step in logfile  
- [ ] makefile
- [ ] Automate test cases for comparison solution against abaqus or hand calculations




