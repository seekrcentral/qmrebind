"""
run_qmrebind_amber.py

Run qmrebind for Amber inputs. Qmrebind reparametrizes the partial charges
on atoms using a QMMM ONIOM calculation on the ligand in the bound state of
the receptor.
"""
import os
import time
import shutil
import argparse

import qmrebind.qmrebind_base as base
import qmrebind.defaults as defaults
import qmrebind.preparation as preparation
import qmrebind.orca as orca
import qmrebind.postprocessing as postprocessing
import qmrebind.check as check 

def move_output(forcefield_file, output=None):
    if output is None:
        print("The reparametrized forcefield file is located at:", 
              os.path.abspath(forcefield_file))
    else:
        print("The reparametrized forcefield file is being written to:", output)
        shutil.copyfile(forcefield_file, output)
    return

def run_qmrebind_amber(
        input_pdb, forcefield_file, ligand_indices=None, ligand_resname="", output=None,
        cut_off_distance=3.0, nprocs=1, maxiter=2000, qm_method="B3LYP", 
        qm_basis_set="6-311G", qm_charge_scheme="CHELPG", qm_charge=None, 
        qm_mult=1, qm2_method="XTB", qm2_charge_scheme="CHELPG", qm2_charge=None, 
        qm2_mult=1, orca_dir_pwd=None, work_dir=None, skip_checks=False,
        keep_solvent_molecules=False):
    """
    Run a full qmrebind calculation on AMBER inputs.
    """
    starttime = time.time()
    
    if output is not None:
        output = os.path.abspath(output)
    
    if orca_dir_pwd is None:
        orca_path = shutil.which("orca")
        orca_dir_pwd = os.path.dirname(orca_path)
        print("Using ORCA at:", orca_path)
    
    #[input_pdb, forcefield_file] = base.make_work_dir(
    #    [input_pdb, forcefield_file], work_dir, overwrite=False, keep_old=True)
    [input_pdb, forcefield_file] = base.make_work_dir(
        [input_pdb, forcefield_file], work_dir, overwrite=True, keep_old=False)
    
    # Getting started with the ORCA calculation using the modified intial PDB 
    # file
    if keep_solvent_molecules:
        preparation.convert_to_TIP3P(forcefield_file, input_pdb)
    preparation.prepare_pdb(input_pdb=input_pdb, 
                            keep_solvent_molecules=keep_solvent_molecules)
    preparation.strip_topology(forcefield_file=forcefield_file,
                               keep_solvent_molecules=keep_solvent_molecules)
    if ligand_indices is None:
        assert ligand_resname != "", \
            "If ligand_indices are not provided, ligand_resname must be."
        qm_region_atom_indices = base.get_indices_qm_region(
            input_pdb=input_pdb, ligand_resname=ligand_resname)
        
    else:
        assert ligand_resname == "", \
            "If ligand_indices are provided, ligand_resname must not be."
        qm_region_atom_indices = ligand_indices
    
    preparation.get_ligand_pdb(
        input_pdb=input_pdb, ligand_pdb=defaults.ligand_pdb, 
        ligand_indices=qm_region_atom_indices)
    preparation.get_receptor_pdb(
        input_pdb=input_pdb, receptor_pdb=defaults.receptor_pdb, 
        ligand_indices=qm_region_atom_indices,
        keep_solvent_molecules=keep_solvent_molecules)
    print(f"The indices for the atoms in the QM region are: "
          f"{qm_region_atom_indices}, and the number of atoms is: "
          f"{len(qm_region_atom_indices)}.")
    assert len(qm_region_atom_indices) > 0, \
        "No atoms in qm region. Incorrect selection?"
    qm2_region_residue_indices, qm2_region_atom_indices \
        = base.get_indices_qm2_region(
            ligand_pdb=defaults.ligand_pdb, input_pdb=input_pdb, 
            cut_off_distance=cut_off_distance, 
            ligand_indices=qm_region_atom_indices)
    print(f"The indices for atoms in the QM2 region are: "
          f"{qm2_region_atom_indices}, and the number of atoms are: "
          f"{len(qm2_region_atom_indices)}.")
    assert len(qm2_region_atom_indices) > 0, \
        "No atoms in qm2 region. Incorrect selection?"
    print(f"The indices for residues in the QM2 region are: "
          f"{qm2_region_residue_indices}, and the number of residues are: "
          f"{len(qm2_region_residue_indices)}.")
    
    if qm_charge is None:
        qm_charge = int(round(base.get_region_charge(
            forcefield_file, qm_region_atom_indices)))
        print(f"Assigning QM region to have charge {qm_charge}.")
        
    if qm2_charge is None:
        qm2_charge = int(round(base.get_region_charge(
            forcefield_file, qm2_region_atom_indices)))
        print(f"Assigning QM2 region to have charge {qm2_charge}.")
    
    base.run_check(check.check_region_charge(
        forcefield_file, qm_region_atom_indices, qm_charge, region="QM"), 
                   skip_checks)
    
    base.run_check(check.check_region_charge(
        forcefield_file, qm2_region_atom_indices, qm2_charge, region="QM2"), 
                   skip_checks)
    
    orca.prepare_orca_pdb(
        input_pdb=input_pdb,
        ligand_pdb=defaults.ligand_pdb,
        orca_pdb=defaults.orca_pdb,
        ligand_indices=qm_region_atom_indices,
        receptor_pdb=defaults.receptor_pdb,
        cut_off_distance=cut_off_distance,
    )
    orca.get_amber_to_orca_prms(forcefield_file=forcefield_file)
    
    # ORCA calculation
    orca.get_orca_input(
        nprocs=nprocs,
        maxiter=maxiter,
        qm_method=qm_method,
        qm_basis_set=qm_basis_set,
        qm2_method=qm2_method,
        optimization=defaults.optimization,
        frequency_calculation=defaults.frequency_calculation,
        qm_charge_scheme=qm_charge_scheme,
        qm2_charge_scheme=qm2_charge_scheme,
        qm2_charge=qm2_charge,
        qm2_mult=qm2_mult,
        forcefield_file=forcefield_file,
        input_pdb=input_pdb,
        ligand_indices=qm_region_atom_indices,
        orca_pdb=defaults.orca_pdb,
        orca_input_file=defaults.orca_input_file,
        ligand_pdb=defaults.ligand_pdb,
        receptor_pdb=defaults.receptor_pdb,
        cut_off_distance=cut_off_distance,
        qm_charge=qm_charge,
        qm_mult=qm_mult,
    )
    
    """ # TODO: marked for removal
    add_xtb_inputs(
        etemp=etemp,
        solvation=solvation,
        solvent=solvent,
        accuracy=accuracy,
        xtb_memory=xtb_memory,
        xtb_nprocs=xtb_nprocs,
        orca_input_file=orca_input_file,
        XTB_add_inputs=XTB_add_inputs,
    )
    """
    #base.run_check(check.check_ligand_same_molecule(
    #    defaults.orca_pdb, qm_region_atom_indices), skip_checks)
    base.run_check(check.check_ligand_same_molecule(
        input_pdb, forcefield_file, qm_region_atom_indices), skip_checks)

    orca.run_orca_qmmm(
        orca_dir_pwd=orca_dir_pwd,
        orca_input_file=defaults.orca_input_file,
        orca_out_file=defaults.orca_out_file,
    )
    
    # Post calculation
    postprocessing.get_qm_charges(
        orca_out_file=defaults.orca_out_file,
        qm_charge_file=defaults.qm_charge_file,
        input_pdb=input_pdb,
        ligand_indices=qm_region_atom_indices,
        qm_charge_scheme=qm_charge_scheme,
    )
    postprocessing.get_ff_charges(
        forcefield_file=forcefield_file,
        ff_charges_file=defaults.ff_charges_file,
        input_pdb=input_pdb,
    )
    postprocessing.get_ff_qm_charges(
        qm_charge_file=defaults.qm_charge_file,
        ff_charges_file=defaults.ff_charges_file,
        ff_charges_qm_fmt_file=defaults.ff_charges_qm_fmt_file,
        input_pdb=input_pdb,
        ligand_indices=qm_region_atom_indices,
    )
    postprocessing.get_qmrebind_parm(
        forcefield_file=forcefield_file,
        input_pdb=input_pdb,
        ff_charges_qm_fmt_file=defaults.ff_charges_qm_fmt_file,
    )
    
    # TODO: Check the charges from the ORCA output directly with the parm7 file
    
    # Post Analysis
    postprocessing.get_qmrebind_parm_solvent(
        input_pdb=input_pdb,
        forcefield_file=forcefield_file,
        ff_charges_file=defaults.ff_charges_file,
    )
    new_no_solvent_pdb_name = base.rename_receptorligand_pdb(
        input_pdb=input_pdb)
    if not skip_checks:
        ff_base = os.path.splitext(forcefield_file)[0]
        ff_ext = os.path.splitext(forcefield_file)[1]
        forcefield_file_before_qm_no_solvent \
            = f"{ff_base}_before_charge_replacement{ff_ext}"
        forcefield_file_after_qm_no_solvent \
            = f"{ff_base}_no_solvent{ff_ext}"
        print("Energy differences without solvent:")
        check.get_energy_diff(
            forcefield_file=forcefield_file_after_qm_no_solvent, 
            forcefield_file_before_qm=forcefield_file_before_qm_no_solvent, 
            pdb_file_before_qm=new_no_solvent_pdb_name)
        input_pdb_base = os.path.splitext(input_pdb)[0]
        ff_base = os.path.splitext(forcefield_file)[0]
        ff_ext = os.path.splitext(forcefield_file)[1]
        forcefield_file_before_qm_solvent = f"{ff_base}_before_qmmm{ff_ext}"
        print("Energy differences with solvent:")
        check.get_energy_diff(
            forcefield_file=forcefield_file, 
            forcefield_file_before_qm=forcefield_file_before_qm_solvent, 
            pdb_file_before_qm=input_pdb)
        print("Running OpenMM simulation to test stability.")
        check.run_openmm_sim(
            input_pdb=input_pdb, forcefield_file=forcefield_file, 
            sim_steps=defaults.sim_steps, T=defaults.T)
        check.get_charge_diff_file(
            forcefield_file=forcefield_file,
            ligand_pdb=defaults.ligand_pdb,
            ligand_charge_diff_file=defaults.ligand_charge_diff_file,
        )
    
    total_time = time.time()-starttime
    print(f"QMREBIND CALCULATION FINISHED! Time: {total_time:.3f} s")
    move_output(forcefield_file, output)
    return forcefield_file
    

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "input_pdb", metavar="INPUT_PDB", type=str,
        help="The PDB file of the system to be reparametrized. The PDB file "\
        "defines the atomic positions of the system. It also must contain "\
        "box vector information in the form of a CRYST line.")
    argparser.add_argument(
        "forcefield_file", metavar="FORCEFIELD_FILE", type=str, 
        help="The name of the forcefield file whose atomic point charges "\
        "will be refined. For the AMBER forcefield, this would be a .parm7 "\
        "or a .prmtop file.")
    argparser.add_argument(
        "-l", "--ligand_indices", dest="ligand_indices", 
        metavar="LIGAND_INDICES", type=str, default="",
        help="A comma-separated list of integers defining ligand within the "\
        "input_pdb structure. Ex: -l '1,2,0'. Either the '-l' or '-L' "\
        "arguments must be included.")
    argparser.add_argument(
        "-L", "--ligand_resname", dest="ligand_resname", 
        metavar="LIGAND_RESNAME", type=str, default="",
        help="The residue name of the ligand molecule for automatic index "\
        "selection. Either the '-l' or '-L' arguments must be included.")
    argparser.add_argument(
        "-o", "--output", dest="output", default=None,
        help="A path to an output file name for the forcefield file. If left "\
        "at the default of None, a file will be generated in the work "\
        "directory with the same name as the input forcefield file. "\
        "Default: None.", type=str)
    argparser.add_argument(
        "-c", "--cut_off_distance", dest="cut_off_distance", default=3.0,
        help="The cut-off distance (in Angstroms) used to define the QM2 "\
        "region of the ONIOM calculation. Default: 3.0.", type=float)
    argparser.add_argument(
        "-n", "--nprocs", dest="nprocs", default=1,
        help="The number of processors to use for ORCA calculations. "\
        "Default: 1.", type=int)
    argparser.add_argument(
        "-i", "--max_iterations", dest="max_iterations", default=2000,
        help="The maximum number of iterations for ORCA convergence. "\
        "Default: 2000.", type=int)
    argparser.add_argument(
        "-m", "--qm_method", dest="qm_method", default="B3LYP",
        help="The method to use for the QM region of the ONIOM calculation. "\
        "Please see the file orca_methods_basis_sets.pdf for all possible "\
        "options. Default: B3LYP.", type=str)
    argparser.add_argument(
        "-b", "--qm_basis_set", dest="qm_basis_set", default="6-311G",
        help="The basis set to use for the QM region of the ONIOM "\
        "calculation. Please see the file orca_methods_basis_sets.pdf for all "\
        "possible options. Default: 6-311G", type=str)
    argparser.add_argument(
        "-s", "--qm_charge_scheme", dest="qm_charge_scheme", default="CHELPG",
        help="The charge scheme to use for the QM region of the ONIOM "
        "calculation. Please see the file orca_methods_basis_sets.pdf for all "\
        "possible options. Default: CHELPG.", type=str)
    argparser.add_argument(
        "-q", "--qm_charge", dest="qm_charge", default=None,
        help="The total charge of the QM region of the ONIOM calculation. "\
        "If set to None, the quantity will be automatically computed. "\
        "Default: 0.", type=int)
    argparser.add_argument(
        "-u", "--qm_multiplicity", dest="qm_multiplicity", default=1,
        help="The multiplicity of the QM region of the ONIOM calculation. "\
        "Default: 1.", type=int)
    argparser.add_argument(
        "-M", "--qm2_method", dest="qm2_method", default="XTB",
        help="The method to use for the QM2 region of the ONIOM calculation. "\
        "Please see the file orca_methods_basis_sets.pdf for all possible "\
        "options. Default: XTB.", type=str)
    argparser.add_argument(
        "-S", "--qm2_charge_scheme", dest="qm2_charge_scheme", default="CHELPG",
        help="The charge scheme to use for the QM2 region of the ONIOM "
        "calculation. Please see the file orca_methods_basis_sets.pdf for all "\
        "possible options. Default: CHELPG.", type=str)
    argparser.add_argument(
        "-Q", "--qm2_charge", dest="qm2_charge", default=None,
        help="The total charge of the QM2 region of the ONIOM calculation. "\
        "If set to None, the quantity will be automatically computed. "\
        "Default: None.", type=int)
    argparser.add_argument(
        "-U", "--qm2_multiplicity", dest="qm2_multiplicity", default=1,
        help="The multiplicity of the QM2 region of the ONIOM calculation. "\
        "Default: 1.", type=int)
    argparser.add_argument(
        "-O", "--orca_path", dest="orca_path", default=None,
        help="An absolute path to the ORCA program. If not specified, ORCA "\
        "will be found from the shutil.which() command. Default: None.", 
        type=str)
    argparser.add_argument(
        "-w", "--work_dir", dest="work_dir", default=None,
        help="A working directory for all temporary files, log files, etc. "\
        "If left at None, a work directory will be automatically generated "\
        "in the current directory. Default: None.", type=str)
    argparser.add_argument(
        "-x", "--skip_checks", dest="skip_checks", default=False, 
        help="By default, checks will be run at various stages of Qmrebind, "\
        "and if the checks fail, the calculation will not proceed. This "\
        "argument bypasses those checks and allows the calculation to "\
        "proceed anyways. Default: False.", action="store_true")
    #argparser.add_argument(
    #    "-I", "--use_implicit", dest="use_implicit", default=None,
    #    help="Use the Universal Solvation Model (SMD) in ORCA to account for "\
    #    "solvation implicitly. Make sure to specify the solvent itself. "\
    #    "Example: '--use_implicit WATER'.")
    argparser.add_argument(
        "-k", "--keep_solvent_molecules", dest="keep_solvent_molecules", 
        default=False, help="Toggle to keep solvent molecules in the QM2 "\
        "region of the calculation.", action="store_true")
    
    args = argparser.parse_args()
    args = vars(args)
    input_pdb = args["input_pdb"]
    forcefield_file = args["forcefield_file"]
    ligand_indices = args["ligand_indices"]
    if ligand_indices != "":
        ligand_indices = base.initialize_indices(ligand_indices)
    else:
        ligand_indices = None
    ligand_resname = args["ligand_resname"]
    output = args["output"]
    cut_off_distance = args["cut_off_distance"]
    nprocs = args["nprocs"]
    max_iterations = args["max_iterations"]
    qm_method = args["qm_method"]
    qm_basis_set = args["qm_basis_set"]
    qm_charge_scheme = args["qm_charge_scheme"]
    qm_charge = args["qm_charge"]
    qm_multiplicity = args["qm_multiplicity"]
    qm2_method = args["qm2_method"]
    qm2_charge_scheme = args["qm2_charge_scheme"]
    qm2_charge = args["qm2_charge"]
    qm2_mult = args["qm2_multiplicity"]
    orca_path = args["orca_path"]
    work_dir = args["work_dir"]
    skip_checks = args["skip_checks"]
    keep_solvent_molecules = args["keep_solvent_molecules"]
    
    run_qmrebind_amber(
        input_pdb, forcefield_file, ligand_indices, ligand_resname, output=output,
        cut_off_distance=cut_off_distance, nprocs=nprocs, 
        maxiter=max_iterations, qm_method=qm_method, qm_basis_set=qm_basis_set,
        qm_charge_scheme=qm_charge_scheme, qm_charge=qm_charge, 
        qm_mult=qm_multiplicity, qm2_method=qm2_method, 
        qm2_charge_scheme=qm2_charge_scheme, qm2_charge=qm2_charge, 
        qm2_mult=qm2_mult, orca_dir_pwd=orca_path, work_dir=work_dir,
        skip_checks=skip_checks, keep_solvent_molecules=keep_solvent_molecules)
    
    # TODO: extract QM1 and QM2 charges from existing parm7 files.