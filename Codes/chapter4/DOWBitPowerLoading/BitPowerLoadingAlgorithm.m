classdef BitPowerLoadingAlgorithm
    methods(Static)
		function [nSC,nIdleLF,nIdleHF,Gamma,Sigma]=initFunc(SNR)
            nSC=evalin('base','nSC');
            nIdleLF=evalin('base','nIdleLF');
            nIdleHF=evalin('base','nIdleHF');
			targetBER = 10^-4;
            Gamma = -log(5*targetBER)/1.5;
            Sigma = 1/(10^(SNR/10)); % sigma^2
        end
        %% Hughers_Hartogs Algorithm
        function [loadedBit,loadedPower] =...
                Hughers_Hartogs(H,SNR,totalPower,targetRate)
            %Hughes-Hartogs Algorithm
			[nSC,nIdleLF,nIdleHF,Gamma,Sigma]=...
                BitPowerLoadingAlgorithm.initFunc(SNR);
            loadedBit = zeros(nSC,1);
            loadedPower = zeros(nSC,1);
            usedPower = 0;
            archivedRate = 0;
            while(archivedRate<targetRate)    
                minDiffPower = 10000;
                index = 0;
                for i=nIdleLF+1:nSC-nIdleHF
                    p = Gamma*Sigma/abs(H(i))^2*2^loadedBit(i);
                    if(p<minDiffPower && loadedBit(i)<=10)
                        minDiffPower = p;
                        index = i;
                    end

                end
                loadedBit(index) = loadedBit(index)+1;
                usedPower = usedPower+minDiffPower;
                archivedRate = archivedRate+1;
            end
            for i = nIdleLF+1:nSC-nIdleHF
                loadedPower(i) = Gamma*Sigma/abs(H(i))^2*...
                    (2^loadedBit(i)-1) * totalPower/usedPower;
            end
        end
        
        %% Chow Algorithm
        function [loadedBit,loadedPower] =...
           		 Chow(H,SNR,totalPower,targetRate)
             
           	[nSC,nIdleLF,nIdleHF,Gamma,Sigma]=...
                BitPowerLoadingAlgorithm.initFunc(SNR);
			snrOnCarriers = abs(H).^2/Sigma;
			gammaMargin = 0;
			maxCount = 10;
			iterateCount = 0;
			usedCarriers = nSC-nIdleLF-nIdleHF;
			loadedBit = zeros(nSC,1);
			loadedPower = zeros(nSC,1);
			diffBit = zeros(nSC,1);
            archivedRate = 0;
            while(iterateCount<maxCount && archivedRate~=targetRate) 
                archivedRate = 0;
                for i = nIdleLF+1:nSC-nIdleHF
                    b = log2(1+snrOnCarriers(i)/...
                        (Gamma*10^(gammaMargin/10)));                    
                    if(round(b)>10)
                        loadedBit(i) = 10;
                    else
                        loadedBit(i) = round(b);
                    end
                    diffBit(i) = b-loadedBit(i);
                    if(loadedBit(i)==0)
                        usedCarriers = usedCarriers-1;
                    end
                    archivedRate = archivedRate + loadedBit(i);
                end

                if(archivedRate==0)
                    error('The channel is too bad to transmit');
                end
                iterateCount = iterateCount + 1;
                gammaMargin = gammaMargin+10*...
                    log10(2^((archivedRate-targetRate)/usedCarriers));
                
                usedCarriers = nSC-nIdleLF-nIdleHF;
                
            end
            
            while(archivedRate>targetRate)
                index = 0;
                minDiffBit = 10000;
                for i = nIdleLF+1:nSC-nIdleHF
                    if(diffBit(i)<minDiffBit && loadedBit(i)>0)
                        minDiffBit = diffBit(i);
                        index = i;
                    end
                end
                loadedBit(index) = loadedBit(index)-1;
                diffBit(index) = diffBit(index)+1;
                archivedRate = archivedRate-1;
            end
            
            while(archivedRate<targetRate)
                index = 0;
                maxDiffBit = -10000;
                for i = nIdleLF+1:nSC-nIdleHF
                    if(diffBit(i)>maxDiffBit && loadedBit(i)<10)
                        maxDiffBit = diffBit(i);
                        index = i;
                    end
                end
                loadedBit(index) = loadedBit(index)+1;
                diffBit(index) = diffBit(index)-1;
                archivedRate = archivedRate+1;
            end
            
            usedPower = 0;
            for i = nIdleLF+1:nSC-nIdleHF
                loadedPower(i) = Gamma/snrOnCarriers(i)*...
                    (2^loadedBit(i)-1);
                usedPower = usedPower+loadedPower(i);
            end
            loadedPower = loadedPower*totalPower/usedPower;
        end
        %% Fischer Algorithm
        function [loadedBit,loadedPower] =...
                Fischer(H,SNR,totalPower,targetRate)
        
            [nSC,nIdleLF,nIdleHF,Gamma,Sigma]=...
                BitPowerLoadingAlgorithm.initFunc(SNR);
			normNoiseOnCarriers = Sigma ./ abs(H).^2;
            LDN = log2(normNoiseOnCarriers(1:nSC));
            usedCarriers = nSC-nIdleLF-nIdleHF;
            loadedBit = zeros(nSC,1);
            isIdleCarriers = zeros(nSC,1);
            loadedPower = zeros(nSC,1);
            flag = 1;
            while(flag)
                flag = 0;
                for i = nIdleLF+1:nSC-nIdleHF
                    if(isIdleCarriers(i)==0)
                        b = (targetRate + sum(LDN))/usedCarriers-LDN(i);
                        if(b<=0)
                            loadedBit(i) = 0;
                            flag = 1;
                            isIdleCarriers(i) = 1;
                            usedCarriers = usedCarriers-1;
                            LDN(i) = 0;
                        else
                            loadedBit(i) = b;
                        end
                    end
                end
            end
            diffBit = loadedBit - (round(loadedBit).*(round(loadedBit)<=10)...
                +10*(round(loadedBit)>10));
            loadedBit = (round(loadedBit).*(round(loadedBit)<=10)...
                +10*(round(loadedBit)>10));
            archivedRate = sum(loadedBit);
            while(archivedRate>targetRate)
                minDiffBit = 10000;
                index = 0;
                for i = nIdleLF+1:nSC-nIdleHF
                    if(isIdleCarriers(i)==0)
                        if(diffBit(i)<minDiffBit && loadedBit(i)>0)
                            minDiffBit = diffBit(i);
                            index = i;
                        end
                    end
                end
                loadedBit(index) = loadedBit(index)-1;
                diffBit(index) = diffBit(index)+1;
                archivedRate = archivedRate -1;
            end
            while(archivedRate<targetRate)
                maxDiffBit = -10000;
                index = 0;
                for i = nIdleLF+1:nSC-nIdleHF
                    if(isIdleCarriers(i)==0)
                        if(diffBit(i)>maxDiffBit && loadedBit(i)<10)
                            maxDiffBit = diffBit(i);
                            index = i;
                        end
                    end
                end
                loadedBit(index) = loadedBit(index)+1;
                diffBit(index) = diffBit(index)-1;
                archivedRate = archivedRate +1;
            end
            usedPower = 0;
            for i = nIdleLF+1:nSC-nIdleHF
                if(isIdleCarriers(i)==0)
                    usedPower = usedPower + ...
                        normNoiseOnCarriers(i)*2^loadedBit(i);
                end
            end
            for i = nIdleLF+1:nSC-nIdleHF
                if(isIdleCarriers(i)==0)
                    loadedPower(i) = totalPower/usedPower*...
                        normNoiseOnCarriers(i)*2^loadedBit(i);
                end               
            end
        end
    end
end
