### Input file
- Receptor coordinate file: 3H8K_A.pdb
- Residues to be blocked by helix: 3H8K.pdb_A
Define residues one in a line as: residue_name chain_id sequence_number
### command
`./globalhepo 3H8K_A.pdb 3H8K.pdb_A 1`
For mpi mode (recommended):
`mpirun -np 2 ./globalhepo 3H8K_A.pdb 3H8K.pdb_A 1`
A file like "3H8K_A.print" are expected.
### Showing helix axis coordinates
`./printhelix 3H8K_A.print`
A file like "3H8K_A_hp.pdb" using oxgen atoms presenting the helix are expected.
