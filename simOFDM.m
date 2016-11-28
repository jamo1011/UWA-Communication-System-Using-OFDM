function [T] = simOFDM(mod,noise,propagationOffset,offsetAmplitude)

SNR = noise;   %bit to noise
modMethod = mod;  %1=BPSK, 2=8QPSK, 3=16QAM, 4=64QAM
offset = propagationOffset;   %Amount of delay from multipath propagation
amplitude = offsetAmplitude;   %Magnitude of amplitude change of delayed symbols(Should be less than one)

nSym = 3;  %Number of OFDM symbols to transmit
N = 64;  %FFT Size/Number of subcarriers
Nsd = 64;  %Number of Data subcarriers
Nsp = 0;   %Number of Pilot subcarriers
ofdmBW = 20*10^6;   %Bandwidth


deltaF = ofdmBW/N;
Tfft = 1/deltaF;   %FFT period
Tgi = round(Tfft/4);   %Gaurd interval duration
Tsignal = Tgi + Tfft;   %Symbol duration
Ncp = 20; %N*Tgi/Tfft;  %Number of symbols for cyclic prefix
Nst = Nsd + Nsp;  %Number of total used subcarriers
previousSymbol = zeros(1,Ncp+Nst);    
nBitsPerSym = zeros(1,length(modMethod));
modNames = {'BPSK';'8QPSK';'16QAM';'64QAM'};


%Finds the bits per symbol for each modulation method tested
for i=1:length(modMethod)
    switch modMethod(i)
        case 1
            nBitsPerSym(i)=1;
        case 2 
            nBitsPerSym(i)=3;
        case 3 
            nBitsPerSym(i)=4;
        case 4 
            nBitsPerSym(i)=6;
    end
end
    
    
errors = zeros(length(SNR),length(modMethod));

for m=1:length(modMethod)
    
    for i=1:length(SNR)
        for j=1:nSym
            %% Transmitter
            %Creates random binary data
            s = round(rand(1,Nst*nBitsPerSym(m)));
            
            %Modulator
            switch modMethod(m)
                case 1    %BPSK
                    sMod = zeros(1,length(s)/nBitsPerSym(m));
                    for k=1:length(sMod)
                        if s(k)==1
                            sMod(k) = 1;
                        else
                            sMod(k) = -1;
                        end
                    end
                case 2   %8QPSK
                    sMod = zeros(1,length(s)/nBitsPerSym(m));
                    C = cos(pi/8);
                    S = sin(pi/8);
                    for k=1:length(sMod)
                        switch num2str(s(3*k-2:3*k))
                            case '0  0  0'
                                sMod(k)=C+S*sqrt(-1);
                            case '0  0  1'
                                sMod(k)=S+C*sqrt(-1);
                            case '0  1  1'
                                sMod(k)=-S+C*sqrt(-1);
                            case '0  1  0'
                                sMod(k)=-C+S*sqrt(-1);
                            case '1  1  0'
                                sMod(k)=-C-S*sqrt(-1);
                            case '1  1  1'
                                sMod(k)=-S-C*sqrt(-1);
                            case '1  0  1'
                                sMod(k)=S-C*sqrt(-1);
                            case '1  0  0'
                                sMod(k)=C-S*sqrt(-1);
                            otherwise 
                                sMod(k)=1;
                        end
                    end
                  
                case 3   %16QAM
                    sMod = zeros(1,length(s)/nBitsPerSym(m));
                    C = 1;
                    S = .33333;
                    for k=1:length(sMod)
                        switch num2str(s(4*k-3:4*k))
                            case '0  0  0  0'
                                sMod(k)=-C+C*sqrt(-1);
                            case '0  0  0  1'
                                sMod(k)=-S+C*sqrt(-1);
                            case '0  0  1  0'
                                sMod(k)=S+C*sqrt(-1);
                            case '0  0  1  1'
                                sMod(k)=C+C*sqrt(-1);
                            case '0  1  0  0'
                                sMod(k)=-C+S*sqrt(-1);
                            case '0  1  0  1'
                                sMod(k)=-S+S*sqrt(-1);
                            case '0  1  1  0'
                                sMod(k)=S+S*sqrt(-1);
                            case '0  1  1  1'
                                sMod(k)=C+S*sqrt(-1);
                            case '1  0  0  0'
                                sMod(k)=-C-S*sqrt(-1);
                            case '1  0  0  1'
                                sMod(k)=-S-S*sqrt(-1);
                            case '1  0  1  0'
                                sMod(k)=S-S*sqrt(-1);
                            case '1  0  1  1'
                                sMod(k)=C-S*sqrt(-1);
                            case '1  1  0  0'
                                sMod(k)=-C-C*sqrt(-1);
                            case '1  1  0  1'
                                sMod(k)=-S-C*sqrt(-1);
                            case '1  1  1  0'
                                sMod(k)=S-C*sqrt(-1);
                            case '1  1  1  1'
                                sMod(k)=C-C*sqrt(-1);
                            otherwise 
                                sMod(k)=0;
                        end
                    end
                    
                    
                case 4   %64QAM
                    %TODO
              
            end
            
            %IFFT
            xTime = N/sqrt(Nst)*ifft(sMod);
            
            
            %Add cyclic Prefix
            ofdmSignal = [xTime(N-Ncp+1:N) xTime];
            
            
            %% Channel
            r = awgn(ofdmSignal,SNR(i));
            
            r = isi(r,previousSymbol,offset,amplitude);
                    
            
            
            previousSymbol = ofdmSignal;   
            %% Reciever
            %Removes cyclic prefix
            rParallel = r(Ncp+1:(N+Ncp));
            
            %FFT
            r_Freq = sqrt(Nst)/N*(fft(rParallel));
            
            
            
            %Demodulator
            switch modMethod(m)
                case 1    %BPSK
                    sDemod = zeros(1,length(r_Freq)*nBitsPerSym(m));
                    for k=1:length(r_Freq)
                        if r_Freq(k)>=0
                            sDemod(k) = 1;
                        else
                            sDemod(k) = 0;
                        end
                    end
                case 2   %QPSK
                    sDemod = zeros(1,length(r_Freq)*nBitsPerSym(m));
                    phase = angle(r_Freq);
                    for k=1:length(phase)
                        if phase(k)<0
                            phase(k) = phase(k)+2*pi;
                        end
                        
                        if phase(k)>=0 && phase(k)<pi/4
                            sDemod(3*k-2:3*k) = [0 0 0];
                        elseif phase(k)>=pi/4 && phase(k)<pi/2
                            sDemod(3*k-2:3*k) = [0 0 1];
                        elseif phase(k)>=pi/2 && phase(k)<3*pi/4
                            sDemod(3*k-2:3*k) = [0 1 1];    
                        elseif phase(k)>=3*pi/4 && phase(k)<pi
                            sDemod(3*k-2:3*k) = [0 1 0];
                        elseif phase(k)>=pi && phase(k)<5*pi/4
                            sDemod(3*k-2:3*k) = [1 1 0];
                        elseif phase(k)>=5*pi/4 && phase(k)<3*pi/2
                            sDemod(3*k-2:3*k) = [1 1 1];
                        elseif phase(k)>=3*pi/2 && phase(k)<7*pi/4
                            sDemod(3*k-2:3*k) = [1 0 1];
                        elseif phase(k)>=7*pi/4 && phase(k)<2*pi
                            sDemod(3*k-2:3*k) = [1 0 0];
                        end
                    end                   
                case 3   %16QAM
                    sDemod = zeros(1,length(r_Freq)*nBitsPerSym(m));
                    for k=1:length(r_Freq)
                        im=imag(r_Freq(k));
                        re=real(r_Freq(k));
                        if im>=.6667
                            sDemod(k*4-3:k*4-2) = [0 0];
                        elseif im<.6667 && im>0
                            sDemod(k*4-3:k*4-2) = [0 1];
                        elseif im<0 && im>-.6667
                            sDemod(k*4-3:k*4-2) = [1 0];
                        elseif im<-.6667
                            sDemod(k*4-3:k*4-2) = [1 1];
                        end    
                        
                        if re>=.6667
                            sDemod(k*4-1:k*4) = [1 1];
                        elseif re<.6667 && re>0
                            sDemod(k*4-1:k*4) = [1 0];
                        elseif re<0 && re>-.6667
                            sDemod(k*4-1:k*4) = [0 1];
                        elseif re<-.6667
                            sDemod(k*4-1:k*4) = [0 0];
                        end    
                        
                        
                    end
                    
                case 4   %64QAM
                    %TODO
            end
            
            %Error Calculator
            numErrors = sum(abs(sDemod-s));
            errors(i,m) = errors(i,m)+numErrors;
            
        end
        
        if length(SNR)==1
            BPSK = [-1+0*sqrt(-1) 1+0*sqrt(-1)];
            scatterplot(r_Freq);
            hold on;
            scatter(real(BPSK),imag(BPSK),'+r');
            title(modNames(m));
            %hold off;
        end
    end
    
end

bitsSent = nSym*Nst*nBitsPerSym;
[s1 s2] = size(errors);
simulatedBER = []
for i=1:s2
    simulatedBER = [simulatedBER errors(:,i)./bitsSent(i)*100];
end
    


figure;
plot(SNR,simulatedBER(:,1),'r-o',SNR,simulatedBER(:,2),'b-*',SNR,simulatedBER(:,3),'s-y');
title('Add Title');


if length(modMethod)==3 && length(SNR)==1
    BER = [simulatedBER';0];
    T = table(BER,'RowNames',modNames);
    
end

grid on;

end

