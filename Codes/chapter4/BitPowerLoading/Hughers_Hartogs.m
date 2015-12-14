function [BitLoaded,PowerLoaded]=Hughers_Hartogs(SNR)
    Sigma = 1/(10^(SNR/10)); % sigma^2
    %Hughes-Hartogs Algorithm
    BitLoaded = zeros(nSC,1);
    Ptotal = nSC;
    Rtarget = (nSC-nIdleLF-nIdleHF)*6;
    P = 0;
    R = 0;
    while(R<Rtarget)    
    Pmin = 10000;
    index = 0;
    for i=nIdleLF+1:nSC-nIdleHF
        p = Gamma*Sigma/abs(H(i))^2*2^bitLoaded(i);
        if(p<Pmin)
            Pmin = p;
            index = i;
        end

    end
    bitLoaded(index) = bitLoaded(index)+1;
    P = P+p;
    R = R+1;
    end
end