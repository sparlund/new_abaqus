## new_abaqus
new_abaqus is a hobby project for creating an FE-solver reading an FE-mesh following the syntax of the commercial software abaqus. 

The name is a joke from the character YinYang in Silicon Valley, who has the idea for "new netflix".

### TODO & Features implemented
- [x] Implement logic and structure for reading abaqus input files
- [x] Set up classes and functions for nodes, elements, properties and materials
- [ ] Implement logic for different elements
  - [x] 2D tria (S2)
  - [ ] 3D tria (S3)
  - [ ] 3D brick (C3D8)

- [x] Assemble stiffness matrix
- [x] Add support for more keywords
  - [x] *BOUNDARY
  - [x] *CLOAD
  - [ ] *MATERIAL
  - [ ] *STEP
  - [ ] *OUTPUT
- [x] Point load keyword and assemble global load vector
- [x] Solve Ka=f system
  - [x] Modify stiffness matrix and load vector to account for boundary conditions
- [ ] Calculate scalar values on elements
  - [ ] von Mises stress
- [ ] Solve eigen value problem
  - [ ] Calculate mass matrix in element construction
- [ ] Figure out smart way to store results  
  - [ ] Implement the format used in open source result viewer ParaView
- [ ] makefile
- [ ] Automate test cases for comparison solution against abaqus




