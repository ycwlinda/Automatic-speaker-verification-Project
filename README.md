# Automatic-speaker-verification-Project
## EE214A Digital Speech Processing Course Project 

The objective of this project is to find a set of acoustic features and algorithms that predict whether two speech segments are uttered by the same speaker or not,  which represented by 0 or 1. 
The database we use includes 50 male speakers saying the same sentences for 5 times. 

Methods used: 
- Gaussian Mixture Model adapted from a Universal Background Model(GMM-UBM)
- Support Vector Machines using GMM SuperVectors (SVM-GSV)
- Dynamic Time Warping (DTW).

Features used:
- MFCCs
- LPCs
- LPCCs
- SSC
- F0.

> Our system yields satisfying performance with both FPR and TPR less than 15% under clean conditions and less than 35% under noisy conditions, where 10dB babble noise is added to the original speech.
