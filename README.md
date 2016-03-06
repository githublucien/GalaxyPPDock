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
2-1. Running GalaxyPPDock
####################
1. Download the GalaxyPPDock program
   * Download a copy of GalaxyPPDock
      * https://github.com/seoklab/GalaxyPPDock/zipball/master/seoklab-GalaxyPPDock-e5ca21.tar.gz
      * https://github.com/seoklab/GalaxyPPDock/zipball/master/seoklab-GalaxyPPDock-e5ca21.zip
2. Unzip and place the downloaded files
   * unzip seoklab-GalaxyPPDock-e5ca21.zip
   * mv seoklab-GalaxyPPDock-e5ca21 $GALAXYPPDOCK_HOME
      (*example*: GALAXYPPDOCK_HOME=/applic/galaxyppdock)

3. Check the downloaded files
   * There should exist:
      * bin: directory for executables  
      There should be add_miss, energy_ligprotein, ppdock_opt.mpi, and ppdock_csa.mpi
      * data: directory for data files
      * examples: directory for example files

4. Add environment variables in $HOME/.bashrc
   * export GALAXYPPDOCK_HOME=$HOME/galaxyppdock
   * export EXEC_MPI=$HOME/hydra-3.2/mpiexec.hydra
   * (*example*: export GALAXYPPDOCK_HOME=/applic/galaxyppdock)
   * (*example*: export EXEC_MPI=/applic/hydra-3.2/mpiexec.hydra)

5. Move GalaxyPPDock example directory
   * cd $GALAXYPPDOCK_HOME/examples/1cgi (or 1eaw, TA12)

6. Running GalaxyPPDock
   * ../../run.py receptor.pdb ligand.pdb zdock.out (default option)
   * ../../run.py receptor.pdb ligand.pdb zdock.out -pssm receptor.pssm ligand.pssm (using PSSM file)
   * ../../run.py receptor.pdb ligand.pdb zdock.out -n 50 (using multiple core, default=2)

7. Output files of GalaxyPPDock
   * model01.pdb ~ model10.pdb files are complex structure models of GalaxyPPDock.

##3. Release log
   * 06 Mar 2016: The first release of GalaxyPPDock

##4. References
   * H. Lee, H. Park, J. Ko, W.-H. Shin, and C. Seok, GalaxyPPDock: Protein-protein docking by cluster-guided conformational space annealing, in preparation.

##5. Contact
   * chaok@snu.ac.kr
   * compbio.galaxy@gmail.com
