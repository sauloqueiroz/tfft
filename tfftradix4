function X = tfft_radix4(x)
    % tfft_radix4: Recursive tfft-style algorithm with N/4 compression
    % with flexibility of mixed-radix input length N=c*4^k with the efficiency
    % of radix-4.
    % Input: x - input signal (vector), length N.
    %            N=c*4^k, k>=0 and c>0 is a non-power-of-four constant with base case N=c defined in the code. 
    % Output: X - DFT of x

    N = length(x);
    % Base Cases
    if N <= 5
        X = fft(x);
        return;
    end

    x = x(:).';

    N4 = N / 4; % Length of sub-sequences

    % Or direct indexing for clarity (p implicitly loops from 0 to N4-1 via vector ops)
    idx0 = 1:N4;
    idx1 = (N4 + 1) : (2 * N4);
    idx2 = (2 * N4 + 1) : (3 * N4);
    idx3 = (3 * N4 + 1) : N;

    x0p = x(idx0);
    x1p = x(idx1);
    x2p = x(idx2);
    x3p = x(idx3);

    % Butterfly-like computations (sums/differences for radix-4 DIT). However,
    % These are the inputs to the sub-DFTs *BEFORE* multiplication by W_N^kp factors. THERE IS NO COMBINATION AT ALL
    sum1 = x0p + x2p;
    sum2 = x1p + x3p;
    diff1 = x0p - x2p;
    diff2 = x1p - x3p;

    xhat0_input = sum1 + sum2;
    xhat1_input_unscaled = diff1 - 1j * diff2; % Corresponds to (x0 - jx1 - x2 + jx3)
    xhat2_input_unscaled = sum1 - sum2;         % Corresponds to (x0 - x1 + x2 - x3)
    xhat3_input_unscaled = diff1 + 1j * diff2; % Corresponds to (x0 + jx1 - x2 - jx3)

    % PREPROCESSING Twiddle factors for p = 0, 1, ..., N/4 - 1
    p_vec = 0:(N4-1); % Row vector

    % W_N^p terms
    twiddle_factors_p1 = exp(-1j * 2 * pi * p_vec / N);
    % W_N^2p terms
    twiddle_factors_p2 = exp(-1j * 2 * pi * 2 * p_vec / N);
    % W_N^3p terms
    twiddle_factors_p3 = exp(-1j * 2 * pi * 3 * p_vec / N);

    % Apply twiddle factors to get the inputs for recursive calls
    % xhat0 is not multiplied by a twiddle factor (or W_N^0p = 1)
    xhat0 = xhat0_input;
    xhat1 = xhat1_input_unscaled .* twiddle_factors_p1;
    xhat2 = xhat2_input_unscaled .* twiddle_factors_p2;
    xhat3 = xhat3_input_unscaled .* twiddle_factors_p3;

    % Recursive calls for N/4 length DFTs
    Y0 = tfft_radix4(xhat0);
    Y1 = tfft_radix4(xhat1);
    Y2 = tfft_radix4(xhat2);
    Y3 = tfft_radix4(xhat3);

    % Combine results (interleaving for radix-4)
    X = zeros(1, N);
    X(1:4:N) = Y0;          % X_k for k = 4m
    X(2:4:N) = Y1;          % X_k for k = 4m+1
    X(3:4:N) = Y2;          % X_k for k = 4m+2
    X(4:4:N) = Y3;          % X_k for k = 4m+3

end
N = 5*4^5;
x = randn(1, N) + 1j * randn(1, N);

ytfft = tfft_radix4(x);
yfft = fft(x);

norm(ytfft - yfft)
