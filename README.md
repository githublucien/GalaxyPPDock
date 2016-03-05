# GalaxyPPDock: A program for protein-protein docking based on cluster-guieded conformational space annealing.

##0.Remark
The GalaxyPPDock distribution version supports only **Linux 64-bit** OS and binary files compiled with MPI option.
**Linux 32-bit** OS or binary files compiled with serial option are not supported.
GalaxyPPDock takes receptor pdb file, ligand pdb file, and ZDOCK output file.
Optionally, you can use PSSM(Position-Specific Scoring Matrix) for more accurate protein-protein docking.

##1. Installation
GalaxyPPDock is only working **Linux 64-bit** system and binary files compiled with MPI option.
Please install mpich2 first before running GalaxyPPDock.

#####################################################################################
1-1. Install mpich2
####################
1. Download mpich from:
   * http://www.mpich.org/static/downloads/3.2/hydra-3.2.tar.gz
2. Unzip and place the download files
   * tar -xvzf hydra-3.2.tar.gz
3. Move hydra MPI directory
   * cd hydra-3.2
4. Add environment variables in $HOME/.bashrc
   * export LD=ld
5. Configure mpich
   * ./configure
6. Install mpich
   * make

You can find **mpiexec.hydra** file in the current directory.

#####################################################################################
Running GalaxyPPDock
####################
1. Add environment variables in $HOME/.bashrc (or $HOME/.bashrc_profile)
    export GALAXYPPDOCK_HOME=$HOME/galaxyppdock
    export EXEC_MPI=$HOME/hydra-3.2/mpiexec.hydra

    (example)
    export GALAXYPPDOCK_HOME=/home/emiliar/galaxyppdock
    export EXEC_MPI=/home/emiliar/hydra-3.2/mpiexec.hydra

2. cd examples/1cgi (or 1eaw, TA12)
3. ../../run.py receptor.pdb ligand.pdb zdock.out (default option)
   ../../run.py receptor.pdb ligand.pdb zdock.out -pssm receptor.pssm ligand.pssm (using PSSM file)
   ../../run.py receptor.pdb ligand.pdb zdock.out -n 50 (using multiple core, default=2)

model01.pdb ~ model10.pdb files are complex structure models of GalaxyPPDock.

