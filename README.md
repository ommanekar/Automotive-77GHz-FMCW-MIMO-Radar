# Automotive 77GHz FMCW MIMO Radar

This repository contains a MATLAB simulation of an automotive 77 GHz FMCW MIMO radar system. The project demonstrates how modern automotive radars estimate range, velocity, and angle of arrival (AoA) using signal processing techniques.

## Overview

The radar operates at 77 GHz, the standard frequency band for automotive and ADAS applications. A TDM-MIMO configuration is used to form a 16-element virtual uniform linear array, enabling high-resolution angle estimation.

## Key Features

- 77 GHz FMCW radar waveform simulation  
- TDM-MIMO radar with virtual antenna array  
- Range–Doppler processing using 2D FFT  
- 2D CA-CFAR target detection  
- Angle-of-arrival estimation using:
  - FFT-based beamforming
  - MUSIC
  - Spatial Smoothing MUSIC (SS-MUSIC)  
- Comparison between conventional and high-resolution AoA methods

## Algorithms Used

- Fast Fourier Transform (FFT)
- Range–Doppler Processing
- Cell-Averaging CFAR (2D)
- Beamforming
- MUSIC Algorithm
- Spatial Smoothing MUSIC (SS-MUSIC)

