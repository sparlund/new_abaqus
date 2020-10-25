## new_abaqus
new_abaqus is a hobby project for creating an FE-solver reading an FE-mesh following the syntax of the commercial software abaqus. 

The name is a joke from the character YinYang in Silicon Valley, who has the idea for "new netflix".


### To-do & Features implemented
- [x] Implement logic and structure for reading abaqus input files
  - [ ] Disregard unused nodes
- [x] Set up classes and functions for nodes, elements, properties and materials
- [ ] Implement logic for different elements
  - [x] 2D tria (S2)
  - [ ] 3D tria (S3)
  - [ ] 3D brick (C3D8)
- [x] Assemble stiffness matrix
  - [x] Modify stiffness matrix and load vector to account for boundary conditions
- [x] Add support for more keywords
  - [x] *BOUNDARY
  - [x] *CLOAD
  - [ ] *MATERIAL
  - [ ] *STEP
  - [ ] *OUTPUT
- [x] Point load keyword and assemble global load vector
- [x] Solve Ka=f for linear problems
- [ ] Calculate scalar values on elements
  - [ ] Stresses and strains
  - [ ] von Mises stress
- [ ] Solve eigen value problem
  - [ ] Calculate mass matrix in element construction
- [ ] Export results to VTK format to view results in ParaView
  - [x] Nodal displacement
  - [ ] Stresses and strains  
- [ ] makefile
- [ ] Automate test cases for comparison solution against abaqus




