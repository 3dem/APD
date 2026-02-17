# APD

The calculation of the APD for any pair of amyloid structures (in PDB or CIF format) consists of a two-step procedure, each of which is implemented as a separate python script. A third script is auxiliary. The scripts require Biopython.  

## Step 1: contacts.py calculates contacts for an amyloid structure
The first script, contacts.py, generates a list of contacts for a given input structure. Because the main chain within a given layer of the amyloid filament may meander up and/or down the helical axis (which is assumed to be along the Z-direction), the closest contacts between pairs of residues may be between different layers of the amyloid. The input structure needs to contain sufficient layers to capture all these contacts. To facilitate the generation of structures with multiple identical layers, an auxiliary python script called helix.py was implemented to apply helical symmetry operators to a specified subset of chains of an input structure.

The contacts.py script identifies protofilaments in the input structure based on the distances between CA atoms of residues with identical residue numbers in the Z-direction. For each protofilament, it then identifies the chain corresponding to the middle layer. For every residue in the middle chain of each protofilament, the script identifies pairwise residue-residue contacts. Two residues are considered to form a contact if at least one pair of their non-hydrogen side-chain atoms (or CA atoms in the case of glycine) are within 6.5 Å. When multiple atom pairs between the same two residues satisfy this criterion, only the shortest distance is stored. Contacts between residues that are adjacent to each other, or separated by one residue within the same chain, are excluded. Moreover, to avoid counting contacts between side chains that are on opposite sides of the amino acid backbone, the script defines a reference plane for each residue. This plane is defined by the positions of the N and C atoms of that residue, and the position of the N atom of the corresponding residue in the next layer. Then, contacts are only kept if, for both of its residues, the two contacting atoms are on the same side of the plane as its CB atom. For glycines, which lack a CB atom, the planes are ignored.

Contacts are classified as intra- or inter-protofilament, depending on whether the two residues involved belong to the same protofilament. To keep track of contacts between different layers of the amyloid, for each contact the script also stores a ca_offset value, which is the difference in Z-coordinates of the CA atoms of the two residues in the middle layer of the amyloid, divided by the average (~4.75 Å) distance between layers within the amyloid.

Finally, the contacts.py script calculates for each residue in the middle layer of the input amyloid structure whether its CB atom is on the left or the right side of its N-N-C plane. The information about all contacts and the left or right-orientation of all middle-layer residues is written in a CSV file.

## Step 2: compare.py calculates the APD from lists of contacts for two structures

The second script, compare.py, takes the CSV files from two different amyloid structures as inputs to perform the comparison. The two input structures may differ in the extent of the ordered filament core. The script will calculate the number of residues that are ordered in both input structures (Ncommon), as well as the number of residues that are only present in the input structure with the largest ordered core (Nextra).

All contacts in the input CSV files that are between residues that are ordered in both input structures are classified to be either unique to one of the input structures, or in common between the two input structures. Contacts are considered as common contacts if they occur between the same two residues in both input structures, with a distance between non-hydrogen side-chain atoms (or CA atoms for glycines) below 4.5 Å in one of the input structures and below 6.5 Å in the other input structure. The more relaxed distance cutoff for the other input structure provides robustness to side-chain conformations in the calculation of common contacts. Contacts are considered as unique contacts if they have a distance below 4.5 Å in one structure and they are not present in the other structure. 

The APD is based on Ndifferent, which is the number of residues that are ordered in both structures and that participate in a unique contact and/or have a different right/left orientation between the two input structures. The APD is then calculated as:

APD = (N_different + N_extra) / (N_common + N_extra) * 100%

The N_extra terms in the calculation of the APD will penalise structures with distinct ordered cores, effectively considering all of them as different. The compare.py script outputs two different distances: the XY-APD considers contacts to be common regardless of whether a contact happens between the same or different amyloid layers. The XYZ-APD considers contacts to be common, only if the residues that make up the contact have the same relative offsets in the Z-direction. The Z-offsets are integer values calculated from the Z coordinates of the CA atoms, rounded to multiples of 4.75 Å. Because the XY-APD ignores the Z-offsets, the XYZ-APD of a given comparison cannot be lower than its XY-APD. APD values can be calculated for the comparison of individual protofilament folds, as well as for inter-protofilament interfaces. In the latter case, only residues involved in contacts between protofilaments in both structures are considered.

The compare.py script calculates XY-APD and XYZ-APD values for all protofilaments and inter-protofilament interfaces in the input structures. In addition, it generates a schematic representation of the comparison. Residues that are in common among the two input structures are represented by a white circle with a black outline and its one-letter amino acid code in black; residues that are only present in one of the structures are represented with a white circle with a grey outline and a grey one-letter code; connections between the circles represent the main chain. The circles of residues with different left/right orientations between the input structures are filled orange; the circles of residues that are mutated between the two structures are filled red. Contacts that are in common between the two input structures are shown in light grey lines; contacts that are unique are shown in dark orange for intra-protofilament contacts and in light orange for inter-protofilament contacts; contacts that are only unique due to differences in Z-offsets are shown in marine; and contacts that involve residues that are only present in one of the input structures are shown in dark yellow. 

## An example

The below commands run the ADP comparison on two structures of mouse-adapted prion strain:

  python helix.py 7QIG.pdb --chains A --twist -0.64 --rise 4.82 --box_and_pixel 384 1.067 --n_copies 4 --output 7QIG-A.pdb 
  
  python helix.py 8EFU.pdb --chains A --twist -0.565 --rise 4.75 --box_and_pixel 384 1.045 --n_copies 4 --output 8EFU-A.pdb 

  python contacts.py 7QIG-A.pdb
  
  python contacts.py 7LNA-A.pdb

  python compare.py 7QIG-A_contacts.csv 8EFU-A_contacts.csv --flip1 --flip2 --rot1 180 --rot2 70

The helix.py script applies helical symmetry to chain A of each PDB file using the supplied twist (in degrees) and rise (in Å). The --n_copies argument generates 4 additional copies of the chain above and below the original chain. If the atomic model was built into a map generated by helical reconstruction in RELION (He and Scheres, 2017), the --box_and_pixel argument calculates the XY coordinates of the helical axis from the box size of the map (in pixels) and the pixel size (in Å). Alternatively, the helical axis coordinates may be specified using the --axis_xy argument. The contacts.py script requires no arguments, other than the coordinate file generated by the helix.py script. The compare.py script includes arguments to control the orientation of the schematic representation (--flip1 and --flip2), while the --rot1 180 and --rot 70 arguments apply clockwise rotations of 180 and 70 degrees to the respective representations. Note that the mirror and rotation operations of the script do not affect the APD calculation (which is rotation and mirror invariant) and are merely used to orient the schematic representations with the structures for convenience of visual comparison.

