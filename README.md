# AD9361 Catalina Fast-Lock Profile Generator

**Reduced-Calibration Fast Lock Profile Generation for the AD9361 by Trained Min-Max NUFFT and Segmented Chebyshev Interpolation**.

Fourth Year Engineering Project  
**Ben-Gurion University of the Negev** | Dept. of Electrical and Computer Engineering  
**Authors:** Maor Buzaglo & Halel Yosef[cite: 13]

---

## 📌 Overview
Modern RF systems must switch between operating frequencies rapidly while maintaining stable and accurate signal generation[cite: 13]. The AD9361 transceiver achieves this using "Fast Lock" profiles, avoiding slow real-time calibrations[cite: 13]. However, obtaining these profiles normally requires exhaustive laboratory calibration at every required frequency, which does not scale well across the 70 MHz to 6 GHz band[cite: 13].

This project provides a software-based algorithmic engine to drastically reduce the required number of physical calibration measurements[cite: 13]. By treating each Fast Lock register field as an unknown continuous function of frequency, we reconstruct the full register tables from a highly sparse sample set[cite: 13].

## 🧠 Algorithmic Approach
Instead of basic polynomial fitting—which is highly prone to Runge's phenomenon on uniform grids[cite: 13]—this system implements an advanced mathematical pipeline:
1. **Spectral Representation:** Data is transformed using a Non-Uniform Fast Fourier Transform (NUFFT)[cite: 13].
2. **Min-Max Optimization:** The scaling vector of the NUFFT is blindly optimized (`fminsearch`) against the theoretical worst-case error bound to guarantee an error of `< 0.5 LSB`[cite: 4, 13]. This makes the model immune to measurement noise[cite: 4].
3. **Barycentric Chebyshev Interpolation:** The spectral representation is evaluated at Chebyshev nodes and stored as a highly stable barycentric polynomial[cite: 13].
4. **Physical Segmentation & Classification:** The frequency axis is dynamically segmented at RFPLL divider boundaries[cite: 4, 13]. Registers are classified by their physical nature (closed-form synthesizer words, discrete transition lists, or continuous smooth fields) to prevent inappropriate interpolation[cite: 13].

## 📂 Repository Structure

The complete MATLAB source code is divided into modular responsibilities[cite: 13]:

| File | Role |
| :--- | :--- |
| `dataprocess.m` | Cleans, merges, and extracts bit-fields from the raw calibration output[cite: 13]. |
| `pll_registers.m` | Computes closed-form synthesizer registers directly from frequency[cite: 13]. |
| `Catalina_Project.m` | The core engine: Scans, classifies, segments, optimizes, builds, and validates the model[cite: 13]. |
| `Catalina_MakeProfiles.m` | Takes the trained model + target frequencies and generates the 16 hex setup words[cite: 13]. |
| `Catalina_InspectProfiles.m` | Decodes and verifies generated profiles (Reverse check)[cite: 13]. |
| `Catalina_LoadScript.m` | Converts the setup words into a `SPIWrite / WAIT` device script[cite: 13]. |

## 🚀 Usage

The preprocessing script is run once per measurement campaign to produce the register–frequency table[cite: 13]. The remaining modules operate on that table. 

Typical usage of the modeling chain in MATLAB[cite: 13]:

```matlab
% Set cfg.fref in Catalina_Project/default_config first
Catalina_Project                 % 1 = SCAN, 2 = BUILD, 3 = VALIDATE
Catalina_MakeProfiles('selftest')
Catalina_MakeProfiles('model.mat', 'targets.xlsx', 'profiles.xlsx')
Catalina_InspectProfiles('profiles.xlsx', 40)
Catalina_LoadScript('profiles.xlsx', 'fastlock_load.txt')
