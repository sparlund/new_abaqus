## new_abaqus
new_abaqus is a hobby project for creating an FE-solver reading an FE-mesh following the syntax of the commercial software abaqus. 

The name is a joke from the character YinYang in Silicon Valley, who has the idea for "new netflix".

### Features implemented


   
### TODO
- [x] Implement logic and structure for reading abaqus input files
- [x] Set up classes and functions for nodes, elements, properties and materials
- [x] Implement stiffness matrix contirbution for 2D tria element
- [ ] Implement logic for 3D shell element
- [x] Assemble stiffness matrix
- [x] Add support for more keywords
  - [x] *BOUNDARY
  - [x] *CLOAD
  - [ ] *STEP
  - [ ] *OUTPUT
- [x] Point load keyword and assemble global load vector
- [ ] Solve for displacement
  - [ ] Partition stiffness matrix
- [ ] Solve eigen value problem
  - [ ] Calculate mass matrix in element construction
- [ ] makefile
- [ ] Automate test cases for comparison solution against CALFEM and or abaqus




