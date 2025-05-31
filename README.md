# TFFT: The "Twiddless" Fast Fourier Transform Algorithm
This repository is the first public implementation of the "Twiddless" Fast Fourier Transform (TFFT) algorithm (draft in [0]).


TFFT introduces a completely novel divide-and-conquer strategy to decompose an N-point Discrete Fourier Transform (DFT) (detailed ahead) into two N/2-point DFTs at the expense of O(N) arithmetic operations. As we briefly discuss ahead, TFFT has unique characteristics that turn it competitive in practice, namely,

* **Computational complexity and numerical stability**: TFFT runs at the same O(NlogN) time complexity of other classicl FFT algorithms. For power-of-two N, algorithms like radix-2 FFT might outperform TFFT's complexity and stability by a constant factor. However, in TFFT N follows a more generic pattern c*2^k, for k≥0 and c being ANY non-power-of-two constant (assumng the base case N=c is handled in the algorithm, in this version, TFFT handles N=c=3 and N=c=5, in addition to N=2 and N=1). This enables TFFT to dispense zero-padding for an infinite number of cases not covered by classic FFTs, providing it with lower complexity and better stability. Preliminary results show that better stability for several non-power-of-two cases (e.g., N=20, 12, 5*2^6, 5*2^7, 5*2^8, etc)
* **Parallelization**: TFFT is highly parallelizable. As we discuss next, the "division" step of TFFT's decomposition approach consists of N-1 complex additions (per stage of recursion) to decimate the signal to N/2 points. All these N/2 points can be computed independently in parallel. Moreover, both  N/2-point smaller DFTs can also be processed in parallel, as they are independent one another.
* **Novel paradigm of DFT Decomposition**: In TFFT, the combination step constitutes of only index re-positioning rather than multiplications by twiddle factors, thereby "butterflies" are dispensed. The "divide" approach builds on recent works that demonstrate how to decimate a signal (in either time of frequency), while preserving the DFT coefficients proportionally between the original and the decimated signal. For example, if the signal is decimated by a factor of 2, then N/2 DFT coefficients are preserved in the decimated signal. With TFFT, we exploit this property recursively, ensuring that the i-th recursive call readily provides N/2^i DFT coefficients of the original N-point signal. In particular, leveraging on [1], the rectangular index coefficient (RIC) [2] compresses an N-point signal x (for even N=CL) to a C-point signal xhat at the expense of N-1 complex additions and no complex multiplications. Denoting X and Xhat as the DFTs of x and xhat, respectively, the correspondence is given by X[c]=Xhat[cL] for c=0,1,...,C-1. In this case, the signal is analyzed in what we refer to as the compressed domain. For example, for C=N/2, Xhat gives the even indices of X. 

TFFT exploits the RIC property recursively to compute all even-indexed DFT coefficients of the input signal x. To obtain the odd-indexed DFT bins, TFFT modulates x in order to shift its frequency bins by r=1 positions, ensuring that X[2c+1]=Xhat[c], c=0,1,...,C-1. This results in a complexity T(N)=2T(N/2)+O(N)=O(Nlog_2(N)), ensuring TFFT in the category of "fast" DFT algorithms. Precisely, T(N)=10N log_2 N - 4N + c, for a constant c>0. Therefore, TFFT's leading constant is higher than classic FFTs like radix-2, whose leading constant is 5. Notwithstanding this, our TFFT presents novel unique characteristics,
* **Potentially simpler circuit**. The most computationally complex module of TFFT is the frequency shift modulation to compute the odd-indexed DFT bins (as explained above). However, the number of positions to be shift is always 1 (i.e., x[n]exp{-j*2*pi*1*n/N}) favouring circuit simplification. 

  
Future work include
* **Optimization**: This TFFT implementation does not benefit from well known techniques adopted by sophisticated libraries (e.g. fftw), like codelets and hardware tailoring. Future implementations can incorporate these techniques.
* Deeper study of stability and parallelization
* **Extension for sparse FFTs**: radix-2 FFT needs structural modifications to achieve faster runtime when computing (very, very) sparse signals. The TFFT paradigm naturally embodies the notion of sparsity. For example, a 50% compression (RIC compression for C=N/2) followed by a DFT (using any DFT algorithm), gives 50% of the spectrum with a true complexity reduction in nearly 50%.
* **Comprehensive performance evaluation study**: Considering stability and non power of two N.
* **Enhance N to broader pattern**

 
  [0] Queiroz, Saulo. Fast Compressed-Domain Discrete Fourier Transform: The ``Twiddless'' FFT Algorithm. draft available online in https://arxiv.org/abs/2505.23718.
  [1]  S. Queiroz, J. P. Vilela, and E. Monteiro, “Fast computation of the discrete
       fourier transform square index coefficients,” IEEE Signal Processing Magazine
       (Tips & Tricks), vol. 42, issue 2, 2025.
  [2] Saulo Queiroz, João P. Vilela, Benjamin Koon Kei Ng, Chan-Tong Lam, Edmundo Monteiro. “Fast computation of the discrete
       fourier transform rectangular index coefficients,” (under review by IEEE Signal Processing Letters). available online in https://arxiv.org/abs/2504.12551
