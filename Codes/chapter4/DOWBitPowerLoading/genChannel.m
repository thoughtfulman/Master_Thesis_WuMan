function H = genChannel()
    tDS = 18;
    tapSpace = 0.75; % tap spacing is 0.75ns.
    nTap = 64;       % number of taps.
    vector_time=0:tapSpace:tapSpace*(nTap-1);
    if randi(2,1) == 1
        h_Exp = 1 ./ tDS .* exp(-1*vector_time ./ tDS);
        h_Exp = h_Exp * (1/sum(h_Exp.^2))^(1/2);
        hrand = sqrt(3) * rand(1,nTap) .* h_Exp;   % as the sample size is 3 ns, while coarseness of h is 0.75 ns
    else
        a=12*(11/13)^0.5*tDS;
        h_Ceil = 6.*a.^6./(vector_time+a).^7; 
        h_Ceil = h_Ceil * (1/sum(h_Ceil.^2))^(1/2);
        hrand = sqrt(3) * rand(1,nTap) .* h_Ceil;
    end
    H = fft(hrand,256);
    H = H(1:128);
end