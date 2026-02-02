# Decoherence of Majorana Qubits from 1/f Noise

Code for all numerical calculations in *A. Alase, M. C. Goffage, M. C. Cassidy, S. N. Coppersmith, Decoherence of Majorana Qubits from 1/f Noise, arXiv preprint: arXiv:2506.22394.* For Version 3 on the arXiv. 

# Implementation of QPP_Library for: Decoherence in Majorana Qubits by 1/f Noise

This repository contains the code used to generate all numerical results presented in our paper *A. Alase, M. C. Goffage, M. C. Cassidy, S. N. Coppersmith, Decoherence in Majorana Qubits by 1/f Noise*. This code uses the covariance matrix method to calculate the time evolution of the tetron qubit, composed of two Kitaev chains, in the presence of a two-level-fluctuator (TLF) and returns the probability of exciting a quasiparticle pair in a single Kitaev chain of the qubit. See the Methods section and Sec. 4 of the Supplementary Information of our paper for further details on the calculations.  

For a more in-depth overview of our numerical package, called **QPP_Library**, please see Appendix B of our earlier paper, *M. C. Goffage, A. Alase, M. C. Cassidy, S. N. Coppersmith, Leakage at zero temperature from changes in chemical potential in Majorana qubits, arXiv:2504.17485 (2025)*, which employs the same library as the submitted paper.  

Note that since our numerical results are statistical averages over multiple two-level-fluctuator (TLF) noise realisations, we expect minor differences in the plots between consecutive runs.  

---

## Environment and Dependencies
Requires MATLAB 2024a or more recent.  

---

## How to Run the Code
1. Click the green **"Code"** drop-down menu at the top of this page.  
2. Select **"Download ZIP"**.  
3. Unzip the contents into a local directory.  
4. Open **`Run_Code_Demo.m`** or **`Run_all_Figures.m`** in MATLAB.  
5. In MATLAB, go to the **Editor** tab, select the **Editor** tab and then select **Run**.  
6. If prompted to *Change Folder* to the current directory, choose **Change Folder** (highlighted option).  

---

## Program Details – Code Demo
**`Run_Code_Demo.m`** generates the yellow line in Fig. 2b, which is $P_{QPP}$ versus time for a 3 µm nanowire. This is a key numerical result in the paper.  

**`Run_Code_Demo.m`** run time (tested on a MacBook Pro): ~4 minutes.

---

## Program Details – Run All Figures
**`Run_all_Figures.m`** generates all numerical results and corresponding figures presented in the paper. It produces the figures in the following order: Fig. 2c, Fig. 2b, Fig. S1, Fig. 3c, Fig. 3d, Fig. S2, and Fig. S4.  

**`Run_all_Figures.m`** run time: >72 hours.  

You may also run individual figures by opening the **`code`** directory, then opening the figure run script of interest in MATLAB (for example, `Run_Figure2b.m`). You may then run the code for just that figure by selecting **Run** in MATLAB. Note that `Run_Figure3cd.m` requires you to first execute `Run_FigureS1.m`, and `Run_FigureS2.m` requires you to first execute **`Run_Figure2c.m`**.  

---

## Output
All figures are output into the folder **`results`** and are titled according to their corresponding figure and panel in the main text or Supplementary Information of our paper. Additionally, the corresponding MATLAB workspaces generated for each figure are saved in the **`results`** folder.  

---

## Troubleshooting
If you run into any issues or have any questions, please contact m.goffage@unsw.edu.au.  
