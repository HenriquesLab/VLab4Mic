import numpy as np
import yaml
import os
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder, CaPPBuilder
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import matplotlib.pyplot as plt
import random

from ..utils.transform import cif_builder  # Verified
from ..utils.visualisation.matplotlib_plots import (
    draw1nomral_segment,
    add_ax_scatter,
    draw_nomral_segments,
)  # verified

from .labels import Label

from ..utils.data_format.visualisation import (
    format_coordinates,
    set_colorplot,
)  # verified

from ..utils.data_format.structural_format import builder_format  # verified

from ..utils.transform.normals import normals_by_scaling  # verified
from ..utils.sample import arrays

class MolecularStructureParser:
    def __init__(self):
        """
        Initialize a MolecularStructureParser object with default attributes.
        """
        self.id = ""
        self.chain_builder = "CaPPBuilder"
        self.ch_builder_max_dist = 5
        self.ppgen = None
        self.CIFdictionary = None
        self.protein_names = None
        self.assembly_operations = None
        self.chains_dict = None
        self.assembly_refpt = None
        self.assembly_atoms = None
        self.num_assembly_atoms = None
        self.store_atoms_assembly = True
        self.assymetric_defined = None
        self.label_targets = dict()
        self.label_normals = None
        # label_targets is a dictionary contining coordinates and parameters
        # for each labelling used. Keys are the name of the label
        self.plotting_params = dict()
        self.scale = 1e-10  # in meters
        self.axis = dict(pivot=None, direction=None)

    def _clear_labels(self):
        """
        Clear all label targets and normals, resetting plotting parameters except for assembly atoms.
        """
        #print("labels cleared")
        self.label_targets = dict()
        self.label_normals = None
        temp_ = self.plotting_params["assemblyatoms"]
        self.plotting_params = dict()
        self.plotting_params["assemblyatoms"] = temp_

    def initialise_parsers(self, dictionary: dict):
        """
        Create molecular structure of the PDB/CIF based on its format and associated parameters.

        Parameters
        ----------
        dictionary : dict
            Dictionary containing file path, ID, format, and optionally title.
        """
        print("Parsing structure. This might take a few seconds...")
        self.source_file = dictionary["file"]
        self.id = dictionary["ID"]
        self.format = dictionary["format"]
        try:
            self.identifier = dictionary["title"]
        except:
            self.identifier = ""
        if self.format == "PDB":
            parser = PDBParser()
            self.struct = parser.get_structure(self.identifier, self.source_file)
        elif self.format == "CIF":
            parser = MMCIFParser()
            self.struct = parser.get_structure(self.identifier, self.source_file)
        else:
            print(f"{self.format} is not a valid format. Valid options are PDB or CIF.")
        if self.chain_builder == "PPBuilder":
            self.ppgen = PPBuilder(self.ch_builder_max_dist)
        elif self.chain_builder == "CaPPBuilder":
            self.ppgen = CaPPBuilder(self.ch_builder_max_dist)

    def __str__(self):
        return f"Structure object for {self.identifier} of {self.format} format"

    # BUILD STRUCTURE
    def build_structure(self):
        """
        Build the structure by generating chains, assembly operations, and reference point.
        """
        print(f"Building structure for: {self.id}: {self.identifier}...")
        if self.chains_dict is None:
            self.generate_chains_sequences()
        if self.assembly_operations is None:
            self.generate_assemmbly_operations()
        if self.assembly_refpt is None:
            self.generate_assembly_reference_point()

    def generate_MMCIF_dictionary(self):
        """
        Generate a dictionary with all fields parsed from the CIF file.
        """
        # generate the dictionary with all field parsed from CIF file
        self.CIFdictionary = MMCIF2Dict(self.source_file)
        self._gen_protein_names()

    def get_atoms_infile(self):
        """
        Obtain all atom coordinates written in the PDB/CIF file.

        Returns
        -------
        numpy.ndarray
            Array of atom coordinates (N, 3).
        """
        coordinates = []
        for chains in self.struct:
            for chain in chains:
                for residue in chain:
                    for atom in residue:
                        coordinates.append(atom.get_vector().get_array())
        # get array is a method for the class vector of biopython.
        # np.array(atom.get_vector()) does not return the coordinates from the vector
        return np.array(coordinates)

    def generate_chains_sequences(self):
        """
        Construct a dictionary with chain IDs as keys and peptide sequences as values.
        """
        seq_chains_dict = {}
        for chains in self.struct:
            for chain in chains:
                pps = self.ppgen.build_peptides(chain)
                sequence = ""
                for peptide in pps:
                    sequence += peptide.get_sequence()
                seq_chains_dict[chain.get_id()] = str(sequence)
        self.chains_dict = seq_chains_dict
        # nchains = len(list(seq_chains_dict.keys()))
        # print(f"file contains {nchains} chains")
    
    def _gen_protein_names(self):
        """
        Generate a dictionary of protein names and associated chain/strand information.
        """
        protein_name = dict()
        chain_name = self.CIFdictionary["_entity.pdbx_description"]
        chain_number = self.CIFdictionary["_entity.id"]
        strand_id = self.CIFdictionary["_entity_poly.pdbx_strand_id"]
        for  name, number, strands in zip(chain_name, chain_number, strand_id):
            protein_name[name] = {"number": number, "strand_id": strands}
        self.protein_names = protein_name
    
    def list_protein_names(self):
        """
        List all protein names in the structure.

        Returns
        -------
        list of str
            List of protein names.
        """
        if self.protein_names is not None:
            if len(list(self.protein_names.keys())) > 0:
                return list(self.protein_names.keys())
        else:
            return []

    def _random_substring(string, size=5):
        """
        Generate a random substring of a given size from a string.

        Parameters
        ----------
        string : str
            Input string.
        size : int, optional
            Length of the substring. Default is 5.

        Returns
        -------
        str
            Random substring.
        """
        if len(string) > size:
            start = random.randint(0, len(string) - (size + 1))
            end = start + size
            return string[start:end]

    def _sequence_substring(self, string, size=5, position="random"):
        """
        Extract a substring from a sequence at a specified position.

        Parameters
        ----------
        string : str
            Input sequence.
        size : int, optional
            Length of the substring. Default is 5.
        position : str, optional
            Position to extract from ("random", "cterminal", "nterminal").

        Returns
        -------
        str
            Extracted substring.
        """
        if position == "random":
            return self._random_substring(string, size)
        elif position == "cterminal":
            start = len(string) - (size + 1)
            return string[start:(len(string)-1)]
        elif position == "nterminal":
            end = size 
            return string[0:end]

    def get_peptide_motif(self, chain_name=None, chain_id=None, size = 5, position="cterminal"):
        """
        Get a peptide motif from a specified chain and position.

        Parameters
        ----------
        chain_name : str, optional
            Name of the chain.
        chain_id : str, optional
            Chain ID.
        size : int, optional
            Length of the motif. Default is 5.
        position : str, optional
            Position in the chain ("cterminal", "nterminal", "random").

        Returns
        -------
        tuple
            (chain_name, chain_id, position, motif)
        """
        if chain_name is None:
            chains_in_structure = list(self.protein_names.keys())
            chain_name = random.choice(chains_in_structure)
        if chain_id is None:
            chain_id = random.choice(self.protein_names[chain_name]["strand_id"].split(","))
        chain_sequence = self.chains_dict[chain_id]
        return chain_name, chain_id, position, self._sequence_substring(chain_sequence, size=size, position=position)

    def generate_assemmbly_operations(self):
        """
        Parse the rotation/translation operations needed to construct a molecular assembly from an asymmetric unit.
        """
        if self.format == "PDB":
            print("Assembly operations are only defined for CIF format files.")
            self.assymetric_defined = False
        else:
            if self.CIFdictionary is None:  # check if already created
                self.generate_MMCIF_dictionary()
            if "AlphaFoldDB" in self.CIFdictionary["_database_2.database_id"]:
                self.assymetric_defined = False
            else:
                # get to know how many transformations there are
                transformation_ids = self.CIFdictionary["_pdbx_struct_oper_list.id"]
                # each element of the list is pair of matrix-vector needed to transform data
                # iterate over each of the transformations
                assembly_transformations = []
                for tr in range(len(transformation_ids)):
                    row1 = [
                        self.CIFdictionary["_pdbx_struct_oper_list.matrix[1][1]"][tr],
                        self.CIFdictionary["_pdbx_struct_oper_list.matrix[1][2]"][tr],
                        self.CIFdictionary["_pdbx_struct_oper_list.matrix[1][3]"][tr],
                    ]

                    row2 = [
                        self.CIFdictionary["_pdbx_struct_oper_list.matrix[2][1]"][tr],
                        self.CIFdictionary["_pdbx_struct_oper_list.matrix[2][2]"][tr],
                        self.CIFdictionary["_pdbx_struct_oper_list.matrix[2][3]"][tr],
                    ]

                    row3 = [
                        self.CIFdictionary["_pdbx_struct_oper_list.matrix[3][1]"][tr],
                        self.CIFdictionary["_pdbx_struct_oper_list.matrix[3][2]"][tr],
                        self.CIFdictionary["_pdbx_struct_oper_list.matrix[3][3]"][tr],
                    ]
                    rotation_matrix = np.array([row1, row2, row3], dtype="float")

                    vector = np.array(
                        [
                            self.CIFdictionary["_pdbx_struct_oper_list.vector[1]"][tr],
                            self.CIFdictionary["_pdbx_struct_oper_list.vector[2]"][tr],
                            self.CIFdictionary["_pdbx_struct_oper_list.vector[3]"][tr],
                        ],
                        dtype="float",
                    )
                    assembly_transformations.append([rotation_matrix, vector])
                struct_oper_dictionary = dict(zip(transformation_ids, assembly_transformations))
                if len(transformation_ids) > 1:
                    self.assymetric_defined = True
                    # print(
                    #    "This model is defined with more than one symmetric transformation.
                    # Will consider the assembly as assymetric defined"
                    # )
                else:
                    # when no assembly unit exist, there is only 1 transform expected
                    self.assymetric_defined = False
                    # print("This model is defined with only one symmetric transformation")
                self.assembly_operations = struct_oper_dictionary

    def generate_assembly_reference_point(self):
        """
        Generate a reference point for the whole structure by constructing a full assembly and averaging all atom coordinates.
        """
        # this function could use a small portion of the atoms defined only
        allatoms_defined = self.get_atoms_infile()  # Parse all atoms in file
        self.atoms_in_file = allatoms_defined
        if self.assembly_operations is None:
            self.generate_assemmbly_operations()
        if self.assymetric_defined:
            # print(
            #    "Strcture is defined by assymetric units.
            # Generating symmetry partners with the following transformations."
            # )
            # hardcoded: use only transforms with integers as id
            assemblyatoms = cif_builder.sym_transforms_numeric_only(
                self.assembly_operations, allatoms_defined, show_transforms_used=False
            )
            assembly_reference_point = np.mean(assemblyatoms, axis=0)
            if self.store_atoms_assembly:
                # print("Storing assembly atoms")
                self.assembly_refpt = assembly_reference_point
                self.assembly_atoms = assemblyatoms
                self.num_assembly_atoms = np.shape(assemblyatoms)[0]
                assembly_plotting_params = dict(
                    plotsize=20, plotalpha=0.01, plotmarker="o", plotcolour="#000000"
                )
                self.plotting_params["assemblyatoms"] = assembly_plotting_params
            else:
                self.assembly_refpt = assembly_reference_point
            # define axis
            # print("Defining axis from rotations")
            vect = self._central_axis_from_rotations()
            self.set_axis_with_vector(vect)
        else:
            # print(
            #    "Strcture is not defined by assymetric units. Using only atoms in file"
            # )
            if self.store_atoms_assembly:
                # print("Storing assembly atoms")
                self.assembly_refpt = np.mean(allatoms_defined, axis=0)
                assembly_plotting_params = dict(
                    plotsize=20, plotalpha=0.01, plotmarker="o", plotcolour="#000000"
                )
                self.assembly_atoms = allatoms_defined
                self.num_assembly_atoms = np.shape(allatoms_defined)[0]
                self.plotting_params["assemblyatoms"] = assembly_plotting_params
            else:
                self.assembly_refpt = np.mean(allatoms_defined, axis=0)
            self.set_axis_with_vector()  # default [0,0,1]

    def _get_assemblyatoms(self):
        return self.assembly_atoms["coordinates"]

    def set_axis_from_point(self, axis_defining_point):
        """
        Set the axis of the structure using a point.

        Parameters
        ----------
        axis_defining_point : numpy.ndarray
            Point to define the axis direction.
        """
        direction = axis_defining_point - self.assembly_refpt
        self.axis = dict(pivot=self.assembly_refpt, direction=direction)

    # central axis definition

    def set_axis_with_vector(self, vector: np.array = np.array([0, 0, 1])):
        """
        Set the axis of the structure using a direction vector.

        Parameters
        ----------
        vector : numpy.ndarray, optional
            Direction vector. Default is [0, 0, 1].
        """
        self.axis = dict(pivot=self.assembly_refpt, direction=vector)

    def _central_axis_from_rotations(self):
        """
        Calculate the central axis orientation from the assembly operations.

        Returns
        -------
        numpy.ndarray
            Central axis vector.
        """
        average = np.identity(3)
        keys = list(self.assembly_operations.keys())
        for i in range(1, len(self.assembly_operations)):
            average = np.matmul(average, self.assembly_operations[keys[i]][0])
        eigenvalues, eigenvectors = np.linalg.eig(average)
        central_axis = eigenvectors[
            :, np.argmax(np.isclose(eigenvalues, 1.0))
        ]  # central axis is the eigenvector that corresponds to the eigenvalue = 1
        if central_axis.dtype == "complex128":
            central_axis = np.array([0, 0, 1])
        return central_axis

    # Parsing atoms

    def _get_atom_res_chain(self, chainnames: list, resnames: list, atomnames: list):
        """
        Parse the structure and return atom coordinates defined by atom name, residue ID, and chain ID.

        Parameters
        ----------
        chainnames : list of str
            List of chain IDs.
        resnames : list of str
            List of residue names.
        atomnames : list of str
            List of atom names.

        Returns
        -------
        numpy.ndarray
            Array of atom coordinates.
        """
        myatoms = []
        for model in self.struct:
            for chain in model:
                if chain.id in chainnames:
                    for residue in chain:
                        # this conditional checks if the residue is
                        # defined in resname list
                        if residue.resname.strip() in resnames:
                            for atom in residue:
                                if atom.name.strip() in atomnames:
                                    myatoms.append(atom.get_coord())
        return np.array(myatoms)

    def _get_site_specific(
        self, chainnames: list, resnames: list, atomnames: list, position: int
    ):
        """
        Parse the structure and return atom coordinates for a specific site.

        Parameters
        ----------
        chainnames : list of str
            List of chain IDs.
        resnames : list of str
            List of residue names.
        atomnames : list of str
            List of atom names.
        position : int
            Position in the chain.

        Returns
        -------
        numpy.ndarray
            Array of atom coordinates for the specific site.
        """
        print(f"looking in: {chainnames}, {resnames}, {atomnames}, {position}")
        myatoms = []
        for model in self.struct:
            for chain in model:
                if chain.id in chainnames:
                    for residue in chain:
                        #pos = pos + 1
                        pos = residue.id[1]
                        # this conditional checks if the residue is
                        # defined in resname list
                        print(
                            f"resname {pos}: {residue.get_resname()}, "
                        )
                        if residue.resname.strip() in resnames and position == pos:
                            for atom in residue:
                                if atom.name.strip() in atomnames:
                                    myatoms.append(atom.get_coord())
                                    print(
                                        f"site speficic found: {pos}, "
                                        f"{residue.resname.strip()}, "
                                        f"{atom.name.strip()} "
                                    )
        print(np.array(myatoms))
        return np.array(myatoms)

    def get_atom_res_chains_assembly(
        self, residues: list, atoms: list, position: int = None, chains: list[str]=None,
    ):
        """
        Parse the structure and return atom coordinates for specified residues and atoms, including assembly partners if defined.

        Parameters
        ----------
        residues : list of str
            List of residue names.
        atoms : list of str
            List of atom names.
        position : int, optional
            Position in the chain.

        Returns
        -------
        numpy.ndarray
            Array of atom coordinates.
        """
        if self.chains_dict is None:
            self.generate_chains_sequences()  # Initialized as none
        if self.assembly_refpt is None:
            self.generate_assembly_reference_point()  # Initialized as none
        if self.assembly_operations is None:
            self.generate_assemmbly_operations()
        if chains is None:
            list_of_chains = list(self.chains_dict.keys())
        else:
            list_of_chains = chains
        # chains in structure on which to look for
        # print(f"Using chains {list_of_chains}")
        if position is None:
            atom_res_chain = self._get_atom_res_chain(list_of_chains, residues, atoms)
        else:
            atom_res_chain = self._get_site_specific(
                list_of_chains, residues, atoms, position
            )
        # atom_res_chain only considers a signle asymmetric unit
        # if the CIF defines an assembly, the asymmetric partners
        # or atom_res_chain are generated as well
        if self.assymetric_defined:
            # print("Generating symmetry partners for atoms parsed")
            atom_res_chain_assembly = cif_builder.sym_transforms_numeric_only(
                self.assembly_operations, atom_res_chain, show_transforms_used=False
            )
            return atom_res_chain_assembly
        else:
            return atom_res_chain

    def gen_targets_by_atoms(self,
                             label_name = "atoms", 
                             residues = None,
                             atoms = None,
                             position = None,
                             chains = None,
                             **kwargs):
        """
        Generate target locations for labelling defined as direct labeling.

        Parameters
        ----------
        label_name : str, optional
            Name for the label. Default is "atoms".
        residues : list of str, optional
            List of residue names.
        atoms : list of str, optional
            List of atom names.
        position : int, optional
            Position in the chain.
        **kwargs
            Additional keyword arguments (e.g., fluorophore, labeling_efficiency).
        """
        try:
            fluorophore = kwargs["fluorophore"]
        except:
            fluorophore = None
        try:
            labeling_efficiency = kwargs["labeling_efficiency"]
        except:
            labeling_efficiency = None
        colour = set_colorplot(self.plotting_params)
        parsed = self.get_atom_res_chains_assembly(
            residues,
            atoms,
            position,
            chains
        )
        self.label_targets[label_name] = self._wrapup_label_dictionary(
            parsed, labeling_efficiency, fluorophore
        )
        label_plotting_params = dict(
            plotsize=20, plotalpha=1, plotmarker="o", plotcolour=colour
        )
        self.plotting_params[label_name] = label_plotting_params

    # Parsing secuences

    def gen_targets_by_sequence(self, target_name, sequence, **kwargs):
        """
        Generate target locations for labelling defined as indirect labeling.

        Parameters
        ----------
        target_name : str
            Name for the target.
        sequence : str
            Epitope sequence to search for.
        **kwargs
            Additional keyword arguments (e.g., fluorophore, labeling_efficiency, method).
        """
        #  only the first appearance of the epitope in every chain is retrieved
        #  if a chain has more than one epitope, only the first one is returned
        target_seq = sequence
        method = "average"
        if "method" in kwargs.keys():
            method = kwargs["method"]
        label_name = target_name
        try:
            fluorophore = kwargs["fluorophore"]
        except:
            fluorophore = None
        try:
            labeling_efficiency = kwargs["labeling_efficiency"]
        except:
            labeling_efficiency = None
        print(f"Searching for sequence: {target_seq}")
        if self.chains_dict is None:
            self.generate_chains_sequences()
        if self.assembly_refpt is None:
            self.generate_assembly_reference_point()
        if self.assembly_operations is None:
            self.generate_assemmbly_operations()
        epitopes_list = []
        for chains in self.struct:
            for chain in chains:
                pps = self.ppgen.build_peptides(chain)
                sequence = ""
                for peptide in pps:
                    sequence += peptide.get_sequence()
                start = sequence.find(target_seq)
                end = start + len(target_seq)
                if start != -1:
                    coords_chain = []
                    for residue, i in zip(chain, range(len(sequence))):
                        if i >= start and i <= end:
                            for atom in residue:
                                coords_chain.append(atom.get_vector().get_array())
                    coords_chain = np.array(coords_chain)
                    epitopes_list.append(coords_chain)
        # print(len(epitopes_list))
        # print(f"Summarizing epitopes location with method {method}")
        summarised_epitopes = cif_builder.summarize_epitope_atoms(
            epitopes_list, method=method
        )
        colour = set_colorplot(self.plotting_params)
        if self.assymetric_defined:
            print("Generating symmetry partners for epitopes parsed")
            epitopes_locations_assembly = cif_builder.sym_transforms_numeric_only(
                self.assembly_operations,
                summarised_epitopes,
                show_transforms_used=False,
            )
            self.label_targets[label_name] = self._wrapup_label_dictionary(
                epitopes_locations_assembly, labeling_efficiency, fluorophore
            )
            label_plotting_params = dict(
                plotsize=20, plotalpha=1, plotmarker="o", plotcolour=colour
            )
            self.plotting_params[label_name] = label_plotting_params
        else:
            self.label_targets[label_name] = self._wrapup_label_dictionary(
                summarised_epitopes, labeling_efficiency, fluorophore
            )
            label_plotting_params = dict(
                plotsize=20, plotalpha=1, plotmarker="o", plotcolour=colour
            )
            self.plotting_params[label_name] = label_plotting_params

    def _wrapup_label_dictionary(
        self, coordinates, labeling_efficiency, fluorophore, normals=None
    ):
        label_dictionary = dict(
            coordinates=coordinates,  # this makes the difference
            labeling_efficiency=labeling_efficiency,
            fluorophore=fluorophore,
            normals=normals,
        )
        return label_dictionary

    # generate targets from label info
    def gen_Targets(self, target_name,  target_type, target_value, **kwargs):
        """
        Generate targets for labelling based on the type of target.

        Parameters
        ----------
        target_name : str
            Name for the target.
        target_type : str
            Type of target ("Atom_residue" or "Sequence").
        target_value : Any
            Value for the target (e.g., atom/residue info or sequence).
        **kwargs
            Additional keyword arguments.
        """
        # call the generator depending on the target type
        if target_type == "Atom_residue":
            #print(target_value)
            self.gen_targets_by_atoms(target_name, **target_value)
        elif target_type == "Sequence":
            self.gen_targets_by_sequence(
                target_name=target_name,
                sequence=target_value)
        else:
            print(f"Label type {target_type} is not a valid label")

    def get_target_coords(self, targetid: str):
        """
        Get coordinates for a specific target.

        Parameters
        ----------
        targetid : str
            Target identifier.

        Returns
        -------
        numpy.ndarray
            Coordinates of the target.
        """
        return self.label_targets[targetid]["coordinates"]

    def get_target_normals(self, targetid: str):
        """
        Get normals for a specific target.

        Parameters
        ----------
        targetid : str
            Target identifier.

        Returns
        -------
        numpy.ndarray
            Normals of the target.
        """
        return self.label_targets[targetid]["normals"]

    def get_target_colour(self, targetid: str):
        """
        Get plotting colour for a specific target.

        Parameters
        ----------
        targetid : str
            Target identifier.

        Returns
        -------
        str
            Colour code.
        """
        return self.plotting_params[targetid]["plotcolour"]

    def set_labefficiency(self, targetname: str, newefficiency: float):
        """
        Set the labelling efficiency for a target.

        Parameters
        ----------
        targetname : str
            Target identifier.
        newefficiency : float
            New labelling efficiency value.
        """
        self.label_targets[targetname]["labeling_efficiency"] = newefficiency

    # visualisation
    def show_target_labels(
        self,
        labelnames: str = None,
        with_assembly_atoms=False,
        assembly_fraction=0.01,
        reference_point=True,
        view_init=[0, 0, 0],
        axesoff=True,
        show_axis=True,
        with_normals=False,
        return_plot=False,
        axis_object=None,
        target_size = None,
        target_plotcolour = None,
        atoms_size = None,
        atoms_alpha = None,
    ):
        """
        Visualize target labels and optionally assembly atoms and reference point.

        Parameters
        ----------
        labelnames : str, optional
            Name of the label to show. If None, show all.
        with_assembly_atoms : bool, optional
            If True, show assembly atoms. Default is False.
        assembly_fraction : float, optional
            Fraction of assembly atoms to show. Default is 0.01.
        reference_point : bool, optional
            If True, show the reference point. Default is True.
        view_init : list of int, optional
            Initial view angles [elev, azim, roll]. Default is [0, 0, 0].
        axesoff : bool, optional
            If True, hide axes. Default is True.
        show_axis : bool, optional
            If True, show the axis. Default is True.
        with_normals : bool, optional
            If True, show normals. Default is False.
        return_plot : bool, optional
            If True, return the matplotlib figure. Default is False.

        Returns
        -------
        matplotlib.figure.Figure or None
            The figure if return_plot is True, otherwise None.
        """
        if axis_object is not None:
            ax = axis_object
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
        if labelnames is None:
            for trgt in list(self.label_targets.keys()):
                if target_size is not None:
                    self.plotting_params[trgt]["plotsize"] = target_size
                if target_plotcolour is not None:
                    self.plotting_params[trgt]["plotcolour"] = target_plotcolour
                if with_normals:
                    draw_nomral_segments(
                        [self.get_target_normals(trgt), self.get_target_coords(trgt)],
                        ax,
                        colors=["grey", self.get_target_colour(trgt)],
                    )
                else:
                    add_ax_scatter(
                        ax,
                        format_coordinates(
                            self.get_target_coords(trgt), **self.plotting_params[trgt]
                        ),
                    )
        if with_assembly_atoms:
            print(f"Showing {assembly_fraction*100}% of the total atoms")
            atoms_subset = arrays.sample_array(
                self.assembly_atoms, assembly_fraction
            )
            if atoms_size is not None:
                self.plotting_params["assemblyatoms"]["plotsize"] = atoms_size
            if atoms_alpha is not None:
                self.plotting_params["assemblyatoms"]["plotalpha"] = atoms_alpha
            add_ax_scatter(
                ax,
                format_coordinates(
                    atoms_subset, **self.plotting_params["assemblyatoms"]
                ),
            )
        if reference_point:
            print(f"Reference point at {self.assembly_refpt}")
            ax.scatter(
                self.assembly_refpt[0],
                self.assembly_refpt[1],
                self.assembly_refpt[2],
                c="k",
                label="Reoriented",
                s=20,
                marker="x",
            )
            self.assembly_refpt
        if show_axis:
            draw1nomral_segment(self.axis, ax, lenght=150, colors=["g", "y"])
        ax.view_init(elev=view_init[0], azim=view_init[1], roll=view_init[2])
        ax.set_box_aspect(
            [ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")]
        )
        if axesoff:
            ax.set_axis_off()
        if axis_object is not None:
            return ax
        else:
            if return_plot:
                plt.close()
                return fig
            else:
                fig.show

    def show_assembly_atoms(
        self, assembly_fraction=0.01, view_init=[0, 0, 0], axesoff=True,return_plot=False
    ):
        """
        Visualize a fraction of the assembly atoms.

        Parameters
        ----------
        assembly_fraction : float, optional
            Fraction of assembly atoms to show. Default is 0.01.
        view_init : list of int, optional
            Initial view angles [elev, azim, roll]. Default is [0, 0, 0].
        axesoff : bool, optional
            If True, hide axes. Default is True.

        Returns
        -------
        matplotlib.figure.Figure
            The figure object.
        """
        print(f"Showing {assembly_fraction*100}% of the total atoms")
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        atoms_subset = arrays.sample_array(
            self.assembly_atoms, assembly_fraction
        )
        add_ax_scatter(
            ax,
            format_coordinates(atoms_subset, **self.plotting_params["assemblyatoms"]),
        )
        ax.view_init(elev=view_init[0], azim=view_init[1], roll=view_init[2])
        ax.set_box_aspect(
            [ub - lb for lb, ub in (getattr(ax, f"get_{a}lim")() for a in "xyz")]
        )
        if axesoff:
            ax.set_axis_off()
        if return_plot:
            plt.close()
            return fig
        else:
            fig.show()

    def create_instance_builder(self, write=False, savingdir=None):
        """
        Create a dictionary with information from target sites, their normals, axis, reference point, etc.
        This dictionary can be used as source for labelled structures in the simulation. 

        Parameters
        ----------
        write : bool, optional
            If True, write the builder to a YAML file. Default is False.
        savingdir : str, optional
            Directory to save the YAML file. If None, uses current working directory.

        Returns
        -------
        dict
            Instance builder dictionary.
        """
        targets = dict()
        for target_name, value in self.label_targets.items():
            targets[target_name] = dict(
                coordinates=value["coordinates"], normals=value["normals"]
            )
        instance_builder = builder_format(
            targets, self.assembly_refpt, self.scale, self.axis, self.identifier
        )
        if write:
            if savingdir is None:
                savingdir = os.getcwd()
            filename = savingdir + "/instance_builder_" + str(self.identifier) + ".yaml"
            with open(filename, "w") as file:
                print(f"Builder saved in: {filename}")
                yaml.dump(instance_builder, file)
        return instance_builder


class MolecularReplicates(MolecularStructureParser):
    def __init__(self, pdbxinfo):
        """
        Initialize a MolecularReplicates object from structure information.

        Parameters
        ----------
        pdbxinfo : dict
            Dictionary with structure file, title, format, and ID.
        """
        MolecularStructureParser.__init__(self)
        self.initialise_parsers(pdbxinfo)
        self.replicates = 1
        self.label_names = []
        self.plotcolours = {}
        self.label_fluorophore = {}

    def add_label(self, labelobj: Label):
        """
        Add a label to the molecular structure and generate its targets.

        Parameters
        ----------
        labelobj : Label
            Label object to add.
        """
        #prepare target value from label object
        self.gen_Targets(
            target_name=labelobj.get_name(),
            target_type=labelobj.get_target_type(),
            target_value=labelobj.params["target"]["value"]
            )
        self.label_names.append(labelobj.get_name())
        self.plotcolours[labelobj.get_name()] = labelobj.get_plotcolour()
        self.label_fluorophore[labelobj.get_name()] = labelobj.get_fluorophore()

    def assign_normals2targets(self, mode="scaling", target=None):
        """
        Assign normals to targets using a specified method.

        Parameters
        ----------
        mode : str, optional
            Method for assigning normals ("scaling"). Default is "scaling".
        target : str, optional
            Specific target to assign normals to. If None, assign to all.
        """
        print(f"Assigning normals to targets with method: {mode}")
        if target is None:
            for target_name, value in self.label_targets.items():
                if mode == "scaling":
                    normals = normals_by_scaling(
                        self.label_targets[target_name]["coordinates"]
                    )
                    self.label_targets[target_name]["normals"] = normals
        else:
            if mode == "scaling":
                normals = normals_by_scaling(self.label_targets[target]["coordinates"])
                self.label_targets[target]["normals"] = normals


def build_structure_cif(
    cif_file: str, struct_title: str = "", cif_id: str = "", format_type="CIF"
):
    """
    Load and parse a PDB/CIF file and build a structure object.

    Parameters
    ----------
    cif_file : str
        Absolute path of the CIF file.
    struct_title : str, optional
        Title of the structure.
    cif_id : str, optional
        Structure ID.
    format_type : str, optional
        File format type ("CIF" or "PDB"). Default is "CIF".

    Returns
    -------
    MolecularReplicates
        Structure object with parsed atoms and information.
    """
    structure_dictionary = {
        "file": cif_file,
        "title": struct_title,
        "format": format_type,
        "ID": cif_id,
    }
    Molecularstructure = MolecularReplicates(structure_dictionary)
    Molecularstructure.build_structure()
    return Molecularstructure
