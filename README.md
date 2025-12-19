# A-novel-frequency-direction-post-processing-method
The corresponding paper, A Frequency-Direction Post-Processing Method for Improving Energy Concentration and Signal Reconstruction, has been submitted to ICASSP. All source code for this work has been made publicly available on GitHub.

# File Descriptions

## 📁 Main Demonstration Scripts

**`demo_bat.m`** – Bat Echolocation Signal Analysis & Method Comparison  
*Primary reproduction script for the paper's core experiment (Section 3). Loads real bat data (`batdata1.mat`), executes six Time‑Frequency Analysis (TFA) methods (STFT, SST, SET, MSST, SRT, and the proposed method), calculates Renyi entropy for concentration comparison, generates all Time‑Frequency Representations (TFRs) (Figs. 2‑3), performs ridge extraction (Fig. 5), and evaluates signal reconstruction accuracy (Figs. 6‑7). This script fully reproduces the quantitative and visual results presented in the paper.*

**`demo_simulated.m`** – Analysis of Simulated Multi‑Component Signals  
*Supplementary demonstration of TFA performance on a synthetic signal containing a fixed‑frequency component, a Linear Frequency Modulated (LFM) component, and a Quadratic Frequency Modulated (QFM) component. *

**`demo.m`** – Illustration of STFT Limitations  
*Visualizes the two fundamental shortcomings of the Short‑Time Fourier Transform in Section 2.1: energy diffusion in the time‑frequency plane and amplitude distortion for non‑stationary components.*

## ⚙️ Core Algorithm Functions

**`PM_W.m`** – Proposed post‑processing Method 
*Implements the novel post‑processing algorithm central to this work (Section 2.2). It reassigns dispersed time‑frequency coefficients toward the nearest local amplitude ridge along the frequency axis, based on the gradient of the STFT amplitude (Eq. 12). This substantially improves energy concentration and enables accurate ridge-only signal reconstruction (Section 2.3).*

## 📊 Data Files

**`batdata1.mat`** – Bat Echolocation Signal Dataset  
*Contains the real bat echolocation signal used in the main experiment (Section 3). The signal has a sampling frequency of ≈142.857 kHz (1 MHz/7). The script uses a 1060‑point segment corresponding to the interval 1.3‑7.1 ms.*

**`MyColormap2.mat`** – Custom Color Map for Visualization  
*Provides a customized color map used for TFRs in this paper.*

## 📚 Project Overview
This repository contains the complete MATLAB implementation accompanying the paper **“A Frequency‑Direction Post‑Processing Method for Improving Energy Concentration and Signal Reconstruction.”** The code enables full reproduction of the experiments, which comprehensively compare six time‑frequency analysis methods. The core contribution is the novel post‑processing method (`PM_W.m`), which introduces a novel frequency‑direction reassignment strategy to overcome energy diffusion and amplitude loss in TFR, thereby significantly enhancing concentration and reconstruction accuracy for non‑stationary signal analysis.

