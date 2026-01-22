clc; clear; close all;

c  = 3e8;
fc = 77e9;
lambda = c/fc;

B = 1e9;              % Hz
Tchirp = 50e-6;        % s
slope = B/Tchirp;      % Hz/s

Ns = 256;              % Number of Samples per Chirp
Nchirp = 128;          % total chirps
fs = Ns/Tchirp;        % ADC sampling rate

Ntx = 2;
Nrx = 8;
Nv  = Ntx*Nrx;         % virtual elements

% Array spacing
d = lambda/2;

t_fast = (0:Ns-1)/fs;
PRT = Tchirp;                  % chirp repetition time
NchirpPerTx = Nchirp/Ntx;      % chirps per Tx stream
if mod(Nchirp, Ntx) ~= 0
    error('Error');
end

% For per-Tx processing
f_r = (0:Ns-1)*(fs/Ns);
R_axis = (c*f_r)/(2*slope);

fd_axis = (-NchirpPerTx/2:NchirpPerTx/2-1)/(NchirpPerTx*PRT);
v_axis = fd_axis*lambda/2;

fprintf('--- Radar settings ---\n');
fprintf('fc=%.1f GHz, B=%.1f MHz, Tchirp=%.1f us, Ns=%d, Nchirp=%d\n', fc/1e9, B/1e6, Tchirp*1e6, Ns, Nchirp);

% RF / Noise Model

kB = 1.38064852e-23;
Tsys = 290;        % K
NF_dB = 4;         % dB (receiver NF)
NF = 10^(NF_dB/10);

Pn = kB*Tsys*B*NF;                 
fprintf('Noise power over B: Pn=%.3e W\n', Pn);

% Radar equation
Pt_dBm = 20;      % dBm
G_dB   = 20;      % dB (Tx/Rx gain lumped)
Pt = 10^((Pt_dBm-30)/10);   % W
G  = 10^(G_dB/10);

% Scenario

% Two close-angle targets at same range to stress AoA
targets(1).R = 18;     targets(1).v = 0.0;   targets(1).ang = 12;  targets(1).sigma = 0.50;
targets(2).R = 18;     targets(2).v = 0.0;   targets(2).ang = 22;  targets(2).sigma = 0.50;

fprintf('\nVirtual array: Ntx=%d, Nrx=%d => Nv=%d virtual elements\n', Ntx, Nrx, Nv);
fprintf('\n--- Targets (Radar Equation amplitude) ---\n');

% Virtual positions
pos_virtual = (0:Nv-1)*d;

% For reporting + amplitude per target
for i=1:numel(targets)
    R = targets(i).R;
    sigma = targets(i).sigma;

    Pr = (Pt * (G^2) * (lambda^2) * sigma) / ((4*pi)^3 * (R^4));   % very simplified
    A  = sqrt(Pr);  % amplitude scale for complex baseband
    targets(i).A = A;

    fprintf('Target %d: R=%.1f m, v=%.1f m/s, ang=%.1f deg, sigma=%.2f m^2 -> Pr=%.3e W\n',...
        i, targets(i).R, targets(i).v, targets(i).ang, targets(i).sigma, Pr);
end

Nframes = 32;

% Data container per Tx stream:
sigTx = cell(1, Ntx);
for tx=1:Ntx
    sigTx{tx} = complex(zeros(Ns, NchirpPerTx, Nrx));
end

Ysnap = [];

for f=1:Nframes
    frameTx = cell(1, Ntx);
    for tx=1:Ntx
        frameTx{tx} = complex(zeros(Ns, NchirpPerTx, Nrx));
    end

    phi0 = 2*pi*rand(numel(targets),1);

    for tx=1:Ntx
        k_global = tx:Ntx:Nchirp;     % length NchirpPerTx
        k_local  = 1:NchirpPerTx;

        for rx=1:Nrx
            % build signal as sum of targets
            x = complex(zeros(Ns, NchirpPerTx));

            for ti=1:numel(targets)
                R   = targets(ti).R;
                v   = targets(ti).v;
                ang = targets(ti).ang*pi/180;
                A   = targets(ti).A;

                % Beat frequency (range)
                fb = 2*slope*R/c;

                % Doppler frequency
                fd = 2*v/lambda;

                % Virtual element index for this (tx,rx)
                vidx = (tx-1)*Nrx + rx;                 % 1..Nv
                p = pos_virtual(vidx);

                % Spatial phase across array
                ph_spatial = 2*pi*(p*sin(ang))/lambda;

                % Slow-time phase progression (per local chirp)
                ph_dopp = 2*pi*fd*(k_local-1)*PRT;

                % Fast-time phase
                ph_fast = 2*pi*fb*t_fast(:);

                % Add random initial phase each frame
                phase = ph_fast + ph_dopp + (phi0(ti) + ph_spatial);

                x = x + A*exp(1j*phase);
            end

            noiseVar = Pn; % keep consistent with your earlier scaling
            n = sqrt(noiseVar/2)*(randn(Ns,NchirpPerTx) + 1j*randn(Ns,NchirpPerTx));

            frameTx{tx}(:,:,rx) = x + n;
        end
    end

    for tx=1:Ntx
        sigTx{tx} = sigTx{tx} + frameTx{tx};
    end

end

for tx=1:Ntx
    sigTx{tx} = sigTx{tx}/Nframes;
end

wR = hann(Ns);
wD = hann(NchirpPerTx).';

RDpow_sum = zeros(Ns, NchirpPerTx);

RDcube = complex(zeros(Ns, NchirpPerTx, Nrx, Ntx));

for tx=1:Ntx
    for rx=1:Nrx
        x = sigTx{tx}(:,:,rx);

        % Range FFT
        xw = x .* (wR);                        
        Xr = fft(xw, Ns, 1);

        % Doppler FFT
        Xr = Xr .* (wD);                       
        Xrd = fft(Xr, NchirpPerTx, 2);
        Xrd = fftshift(Xrd, 2);

        RDcube(:,:,rx,tx) = Xrd;
        RDpow_sum = RDpow_sum + abs(Xrd).^2;
    end
end

RDdB = 10*log10(RDpow_sum + 1e-30);

figure('Name','RD Map');
imagesc(v_axis, R_axis, RDdB);
axis xy;
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('Range-Doppler Map (TDM-MIMO: sum over Rx & Tx)');
colorbar;
grid on;
ylim([0 40]);

% Simple CA-CFAR (2D)
Tr = 10; Gr = 4;    
Td = 8;  Gd = 3;   
Pfa = 1e-5;

[Nr, Nd] = size(RDpow_sum);
detMask = false(Nr, Nd);

for r = (Tr+Gr+1):(Nr-(Tr+Gr))
    for d2 = (Td+Gd+1):(Nd-(Td+Gd))
        rIdx = (r-(Tr+Gr)):(r+(Tr+Gr));
        dIdx = (d2-(Td+Gd)):(d2+(Td+Gd));

        rG = (r-Gr):(r+Gr);
        dG = (d2-Gd):(d2+Gd);

        block = RDpow_sum(rIdx, dIdx);
        guard = RDpow_sum(rG, dG);

        mask = true(size(block));
        rG_local = (Tr+1):(Tr+2*Gr+1);
        dG_local = (Td+1):(Td+2*Gd+1);
        mask(rG_local, dG_local) = false;

        trainingCells = block(mask);
        noiseLevel = mean(trainingCells);

        Ntrain = numel(trainingCells);
        alpha = Ntrain*(Pfa^(-1/Ntrain)-1);

        threshold = alpha*noiseLevel;

        if RDpow_sum(r,d2) > threshold
            detMask(r,d2) = true;
        end
    end
end

[rDet, dDet] = find(detMask);
fprintf('\nDetections found: %d\n', numel(rDet));

% plot detections on RD map
figure('Name','RD + CFAR');
imagesc(v_axis, R_axis, RDdB); axis xy; hold on;
plot(v_axis(dDet), R_axis(rDet), 'r.', 'MarkerSize', 18);
xlabel('Velocity (m/s)'); ylabel('Range (m)');
title('RD Map + CA-CFAR Detections (TDM-MIMO)');
colorbar; grid on; ylim([0 40]);

if isempty(rDet)
    [~, idxMax] = max(RDpow_sum(:));
    [r0, d0] = ind2sub(size(RDpow_sum), idxMax);
else
    lin = sub2ind(size(RDpow_sum), rDet, dDet);
    [~, im] = max(RDpow_sum(lin));
    r0 = rDet(im); d0 = dDet(im);
end

R0 = R_axis(r0);
v0 = v_axis(d0);
fprintf('\n[Angle comparison at strongest RD cell]\n');
fprintf('Using RD cell at R=%.2f m, v=%.2f m/s\n', R0, v0);

y = complex(zeros(Nv,1));
for tx=1:Ntx
    for rx=1:Nrx
        vidx = (tx-1)*Nrx + rx;
        y(vidx) = RDcube(r0, d0, rx, tx);
    end
end

angGrid = -90:0.1:90;
angRad  = angGrid*pi/180;

pos_ula8 = (0:Nrx-1)*d;

AF_ula = zeros(size(angRad));
AF_virt = zeros(size(angRad));

for k=1:numel(angRad)
    a8 = exp(1j*2*pi/lambda*pos_ula8(:)*sin(angRad(k)));
    av = exp(1j*2*pi/lambda*pos_virtual(:)*sin(angRad(k)));
    AF_ula(k)  = abs(sum(a8))/Nrx;
    AF_virt(k) = abs(sum(av))/Nv;
end

AF_ula_dB  = 20*log10(AF_ula/max(AF_ula));
AF_virt_dB = 20*log10(AF_virt/max(AF_virt));

figure('Name','Beam pattern');
plot(angGrid, AF_ula_dB, 'LineWidth', 1.5); hold on;
plot(angGrid, AF_virt_dB, 'LineWidth', 1.5);
grid on; xlabel('Angle (deg)'); ylabel('Array Factor (dB)');
title('Beam Pattern: ULA-8 vs Virtual-16');
legend('ULA (8)','Virtual MIMO (16)','Location','northeast');
ylim([-40 0]);

[~, i0] = min(abs(angGrid));
thrNull = -35;

iL = find(AF_ula_dB(1:i0) < thrNull, 1, 'last');
iR = i0-1 + find(AF_ula_dB(i0:end) < thrNull, 1, 'first');
if isempty(iL), iL = 1; end
if isempty(iR), iR = numel(angGrid); end
bw_ula = angGrid(iR) - angGrid(iL);

iL = find(AF_virt_dB(1:i0) < thrNull, 1, 'last');
iR = i0-1 + find(AF_virt_dB(i0:end) < thrNull, 1, 'first');
if isempty(iL), iL = 1; end
if isempty(iR), iR = numel(angGrid); end
bw_virt = angGrid(iR) - angGrid(iL);

fprintf('\nApprox null-to-null beamwidth:\n');
fprintf('  ULA-8  : %.2f deg\n', bw_ula);
fprintf('  Virt-16: %.2f deg\n', bw_virt);

NfftAng = 2048;
wA = hann(Nv);
y_win = y(:).*wA;

Yfft = fftshift(fft(y_win, NfftAng));
Pfft = abs(Yfft).^2;
Pfft = Pfft/max(Pfft);
Pfft_dB = 10*log10(Pfft + 1e-12);

% Spatial frequency bins u in [-1,1] approx => sin(theta) ~ u
u = linspace(-1,1,NfftAng);
theta_fft = asind(max(-1,min(1,u)));  % degrees

figure('Name','FFT AoA');
plot(theta_fft, Pfft_dB, 'LineWidth', 1.5);
grid on; xlabel('Angle (deg)'); ylabel('Normalized spectrum (dB)');
title('FFT Angle Spectrum at strongest (Range,Doppler)');
xlim([-90 90]); ylim([-50 0]);


dSpan = 2;  % small neighborhood
snapList = [];
for dd = (d0-dSpan):(d0+dSpan)
    if dd>=1 && dd<=Nd
        ytmp = complex(zeros(Nv,1));
        for tx=1:Ntx
            for rx=1:Nrx
                vidx = (tx-1)*Nrx + rx;
                ytmp(vidx) = RDcube(r0, dd, rx, tx);
            end
        end
        snapList = [snapList, ytmp]; %#ok<AGROW>
    end
end

% If still too few snapshots, also use range neighborhood
rSpan = 1;
for rr = (r0-rSpan):(r0+rSpan)
    if rr>=1 && rr<=Nr
        ytmp = complex(zeros(Nv,1));
        for tx=1:Ntx
            for rx=1:Nrx
                vidx = (tx-1)*Nrx + rx;
                ytmp(vidx) = RDcube(rr, d0, rx, tx);
            end
        end
        snapList = [snapList, ytmp]; %#ok<AGROW>
    end
end

% Covariance
Rxx = (snapList*snapList')/size(snapList,2);

% Forward-Backward averaging (often helps)
J = fliplr(eye(Nv));
Rfb = 0.5*(Rxx + J*conj(Rxx)*J);

% MUSIC settings
nSrc = 2;  % we expect 2 close targets
[EV, ED] = eig(Rfb);
[evals, idx] = sort(real(diag(ED)), 'descend');
EV = EV(:, idx);

En = EV(:, nSrc+1:end); % noise subspace

% MUSIC spectrum
Pm = zeros(size(angRad));
for k=1:numel(angRad)
    a = exp(1j*2*pi/lambda*pos_virtual(:)*sin(angRad(k)));
    Pm(k) = 1/real(a'*(En*En')*a);
end
Pm = Pm/max(Pm);
Pm_dB = 10*log10(Pm + 1e-12);

figure('Name','MUSIC AoA');
plot(angGrid, Pm_dB, 'LineWidth', 1.5);
grid on; xlabel('Angle (deg)'); ylabel('Normalized MUSIC spectrum (dB)');
title('MUSIC Angle Spectrum at strongest (Range,Doppler)');
xlim([-90 90]); ylim([-50 0]);

%% ------------------ SS-MUSIC (Spatial Smoothing) ------------------
% Spatial smoothing helps resolve coherent sources (same range/doppler)
L = 10;                  % subarray length (choose <= Nv)
K = Nv - L + 1;          % number of subarrays

Rss = zeros(L,L);
for k=1:K
    Rss = Rss + Rfb(k:k+L-1, k:k+L-1);
end
Rss = Rss/K;

[EVs, EDs] = eig(Rss);
[evals2, idx2] = sort(real(diag(EDs)), 'descend');
EVs = EVs(:, idx2);
Ens = EVs(:, nSrc+1:end);

pos_ss = (0:L-1)*d;

Pss = zeros(size(angRad));
for k=1:numel(angRad)
    a = exp(1j*2*pi/lambda*pos_ss(:)*sin(angRad(k)));
    Pss(k) = 1/real(a'*(Ens*Ens')*a);
end
Pss = Pss/max(Pss);
Pss_dB = 10*log10(Pss + 1e-12);

figure('Name','SS-MUSIC AoA');
plot(angGrid, Pss_dB, 'LineWidth', 1.5);
grid on; xlabel('Angle (deg)'); ylabel('Normalized MUSIC spectrum (dB)');
title('SS-MUSIC Angle Spectrum at strongest (Range,Doppler)');
xlim([-90 90]); ylim([-50 0]);

%% ===================== Print Top Detections with Angle Est (SS-MUSIC peak) =====================
% estimate angle by SS-MUSIC peak
[~, imax] = max(Pss);
ang_est = angGrid(imax);

% print top detections (range/vel) and same angle estimate as simple output
fprintf('\nTop detections (Range, Velocity, Angle_est) [Virtual Array]:\n');
if isempty(rDet)
    fprintf('No CFAR detections; used global max.\n');
    fprintf('R=%.2f m, v=%.2f m/s, angle~%.1f deg\n', R0, v0, ang_est);
else
    % sort detections by power
    lin = sub2ind(size(RDpow_sum), rDet, dDet);
    [~, order] = sort(RDpow_sum(lin), 'descend');
    topN = min(10, numel(order));
    for k=1:topN
        rr = rDet(order(k)); dd = dDet(order(k));
        fprintf('R=%.2f m, v=%.2f m/s, angle~%.1f deg\n', R_axis(rr), v_axis(dd), ang_est);
    end
end

disp('Done.');
