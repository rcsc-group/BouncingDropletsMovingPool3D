# Bouncing droplets onto a moving pool

Direct numerical simulation code infrastructure for **three-dimensional drop impact onto a moving pool**, supporting collaborative work with the Harris Lab at Brown. The code complements the [associated arXiv preprint](https://arxiv.org/abs/2510.02220), and provides tools for parameter sweeps and visualisation.

<img width="100%" alt="3D_SimulationSnapshot" src="https://github.com/user-attachments/assets/86635a45-fe4e-474c-997e-93bfe12b7299" />

---

## 📌 Features

✅ Three-dimensional Navier-Stokes solver for drop-pool impact scenarios  
✅ Parameter sweep support: resolution, drop velocity, pool velocity  
✅ Two-phase, non-coalescing fluid volume implementation  
✅ High-resolution output and animation capabilities  
✅ Supplementary movies for publication and presentation  

---

## 🛠️ Installation

### 1. Requirements

- [Basilisk](http://basilisk.fr/) (compiled with `qcc`)
- C compiler
- Visualization tools (e.g., ffmpeg, imagemagick, gnuplot)
- Recommended:  
  ```bash
  sudo apt install ffmpeg imagemagick gnuplot
  ```

### 2. Clone the repository

```bash
git clone https://github.com/rcsc-group/BouncingDropletsMovingPool3D
cd BouncingDropletsMovingPool3D
```

### 3. Install Basilisk

See Basilisk’s [installation page](http://basilisk.fr/src/INSTALL) for instructions.

### 4. Running the code

After Basilisk is set up, run the driver code using the provided shell script:

```bash
sh run_master_example.sh
```

The script organizes outputs into folders with summaries and movies.

---

## ⚙️ Key Simulation Options

Edit parameters across various layers to control:

- Drop velocity $U_{\textrm{drop}}$ (via the shell script)
- Pool velocity $U_{\textrm{pool}}$ (via the shell script)
- Resolution and domain size (via the shell script)
- Output frequency and visualization options (inside the driver code itself)

Two-phase non-coalescing fluid volume implementation is based on [V. Sanjay’s code](https://github.com/VatsalSy/Lifting-a-sessile-drop/blob/master/CaseI/two-phaseDOD.h).

---

## 📁 Folder Structure

```bash
.
├── DriverCode/                  # Main driver scripts and code for simulations
│   ├── 3D_Impact_Bounce.c       # Main simulation source code
│   ├── run_master_example.sh    # Shell script to run parameter sweeps and manage outputs
│   └── two-phaseDOD.h           # Header file for two-phase, non-coalescing fluid volume implementation
│
├── SupplementaryMovies/         # Supplementary movies and captions referenced in the manuscript
│   ├── SM1.mp4                  # Movie 1
│   ├── SM2.mp4                  # Movie 2
│   ├── SM3.mp4                  # Movie 3
│   ├── SM4.mp4                  # Movie 4
│   └── SM_Captions.txt          # Captions and descriptions for supplementary movies
│
├── LICENSE                      # License information
└── README.md                    # Project documentation
```

---

## 📊 Outputs

Generates:

Summary files and movies in organized output folders:
- `.mp4` animations for simulation data
- Data files (interface data coordinates) for further analysis

Visualisation can be toggled depending on your architecture and needs.

---

## 📚 Citation

If you use this code or data in your work, please cite the associated preprint (for now):

> Harris, D. M., Alventosa, L. F., Sand, O., Silver, E., Mohammadi, A., Sykes, T. C., ... & Cimpeanu, R. (2025). Bouncing to coalescence transition for droplet impact onto moving liquid pools. arXiv:2510.02220.

BibTeX:
```bibtex
@article{Harris2025Bouncing,
  title={Bouncing to coalescence transition for droplet impact onto moving liquid pools},
  author={Harris, D. M. and Alventosa, L. F. and Sand, O. and Silver, E. and Mohammadi, A. and Sykes, T. C. and ... and Cimpeanu, R.},
  journal={arXiv preprint arXiv:2510.02220},
  year={2025}
}
```

---

## 🧑 Contributing

Feel free to:

- Fork this repo
- Open issues
- Submit pull requests

---

