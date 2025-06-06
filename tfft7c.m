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


function xhat = compactric(x, C)
  L=length(x)/C;
  xhat = zeros(1, C); % Initialize output DFT vector
  for c = 0:C-1
    xhat(c+1) = 0;
    for l = 0:L-1
      xhat(c+1) = xhat(c+1) + x(c + 1 + l*C);
    end
  end
end

% frequency shift modulation
function x_mod = shift_factor(x, r)
  % Frequency shift modulation
  % x_mod = x .* exp(-j * 2 * pi * r * n / N)
  N = length(x);
  n = 0:N-1;
  phase = exp(-1j * 2 * pi * r * n / N);
  x_mod = x(:).' .* phase;
end

% TFFT algorithm. Assume N=k7^k, i.e. C=N/7
%eg N=98=2*7^2 gives 7 DFTs of 14 points each.
function X = tfft(x)
    N = length(x);
    if N <=5
      X=fft(x);  
     return;
    end
    L = 7;
    C = N / L;
    %x_row = x(:).'; 
    % Decomposition
    % references:
    % [1] S. Queiroz, J. P. Vilela, and E. Monteiro, “Fast computation of the discrete Fourier 
    % transform square index coefficients,” IEEE Signal Process. Mag. (Tips & Tricks), vol. 42, issue 2, 2025.
    x_1_7 = x(1:C:N); %compactric(x, C); % optimize this.

    % Frequency modulation step
    x_mod = shift_factor(x, 1); % time samples for the 7c+1 frequencies
    x_2_7 = x(2:C:N); % compactric(x_mod, C); % optimize this.

    x_mod = shift_factor(x, 2); % time samples for the 7c+2 frequencies. Note, it could be shift_factor(x_mod, 1);
    x_3_7 = x(3:C:N); %compactric(x_mod, C); % optimize this.

    x_mod = shift_factor(x, 3); % time samples for the 7c+3 frequencies.
    x_4_7 = x(4:C:N); % compactric(x_mod, C); % optimize this.

    x_mod = shift_factor(x, 4); % time samples for the 7c+4 frequencies.
    x_5_7 = x(5:C:N); %compactric(x_mod, C); % optimize this.

    x_mod = shift_factor(x, 5); % time samples for the 7c+4 frequencies.
    x_6_7 = x(6:C:N); %compactric(x_mod, C); % optimize this.
  
    X_1_7 = tfft(x_1_7);
    X_2_7 = tfft(x_2_7);
    X_3_7 = tfft(x_3_7);
    X_4_7 = tfft(x_4_7);
    X_5_7 = tfft(x_5_7);
    X_6_7 = tfft(x_6_7);


    X = zeros(1, N); 
    X(1:L:N) = X_1_7;
    X(2:L:N) = X_2_7;
    X(3:L:N) = X_3_7;
    X(4:L:N) = X_4_7;
    X(5:L:N) = X_5_7;
    X(6:L:N) = X_6_7;
end

% example
N=4802; %pattern N=k*7^2, c<=5. 98=2*7^2 (k=2)
C=14;
L=7;
x=rand(1,N)+rand(1,N)*i;
yfft=fft(x)
ytfft=tfft(x);


ytfft(31)
yfft(31)
tfft7k.m
Exibindo tfft7k.m.
