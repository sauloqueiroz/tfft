% The Twiddless Fast Fourier transform algorithm.
% Saulo Queiroz, Federal University of Technology, Ponta Grossa, PR, Brazil.
%
% Copyright 2025 Saulo Jorge Beltrao de Queiroz
%
% * If this work is useful to you, please cite:
%********* Saulo Queiroz. Fast Compressed-Domain Discrete Fourier Transform: The ``Twiddless'' FFT Algorithm.
%  Draft available online in https://arxiv.org/abs/2505.23718.
%********** S. Queiroz, J. P. Vilela, and E. Monteiro, “Fast computation of the discrete
%           fourier transform square index coefficients,” IEEE Signal Process. Mag.
%           (Tips & Tricks), vol. 42, issue 2, 2025.
%********** S. Queiroz, J. P. Vilela, and E. Monteiro, “Fast computation of the discrete
% *
% * Licensed under the Apache License, Version 2.0 (the "License");
% * you may not use this file except in compliance with the License.
% * You may obtain a copy of the License at
% *
% *     http://www.apache.org/licenses/LICENSE-2.0
% *
% * Unless required by applicable law or agreed to in writing, software
% * distributed under the License is distributed on an "AS IS" BASIS,
% * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% * See the License for the specific language governing permissions and
% * limitations under the License.
% */
%
% TFFT algorithm. Assume N=c2^k for c=3 or c=5 and k>=0.
function X = tfft(x)
  N = length(x);
if N <= 5
    X = fft(x);
    return;
end
    C = N / 2;

    % Compression for even-indexed frequencies (sum), 
    % references:
    % [1] S. Queiroz, J. P. Vilela, and E. Monteiro, “Fast computation of the discrete Fourier 
    % transform square index coefficients,” IEEE Signal Process. Mag. (Tips & Tricks), vol. 42, issue 2, 2025.
    % [2]  https://arxiv.org/abs/2407.00182  
    x_even = x(1:C) + x(C+1:N);

    % Compression for odd-indexed frequencies (difference)
    x_diff = x(1:C) - x(C+1:N);

    % Apply modulation AFTER compression
    c = 0:C-1;
    phase = exp(-1j * 2 * pi * c / N);
    x_odd = x_diff .* phase;

    % Recursive calls
    X_even = tfft(x_even);
    X_odd = tfft(x_odd);

    % Combine results
    X = zeros(1, N);
    X(1:2:end) = X_even;   
    X(2:2:end) = X_odd;   
end
% example
N=20;
x=rand(1,N)+rand(1,N)*i;
y=fft(x)
ytfft=tfft(x);

