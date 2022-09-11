clear all;
close all;
%TP1 traitement du signal S6

Fe = 24e3;
Te = 1/Fe;
Rb = 3000;

%Génération d'un signal
nb_bits = 100;
bits = randi([0,1],1,nb_bits);

%%  Modulateur 1    %%
%Mappage 
Symboles1 = 2*bits-1;

M1 = 2;
Ts1 = log2(M1)*1/Rb;
Ns1 = Ts1/Te;

%Surechantillonage 
Suite_diracs1 = kron(Symboles1, [1 zeros(1,Ns1-1)]);

%Filtrage
h1  =  ones(1,round(Ns1));
x1 = filter(h1,1,Suite_diracs1);
figure ; plot(x1); title("x1");
figure ;  pwelch(x1,[],[],[],1/Te,'twosided'); title("DSP x1");

%DSP théorique 
f = linspace(-Fe/2,Fe/2,length(x1));
DSP_theorique = Ts1*power(sinc(f*Ts1),2);
figure;
semilogy(f,DSP_theorique);
xlabel("Fréquence (Hz)");
title("DSP théorique x1");


%%  Modulateur 2    %%
x2 = reshape(bits,2,length(bits)/2);
Symboles2 = bit2int(x2,2) - 1.5;

M2 = 4;
Ts2 = log2(M2)*1/Rb;
Ns2 = Ts2/Te;

%Surechantillonage 
Suite_diracs2 = kron(Symboles2, [1 zeros(1,Ns2-1)]);

%Filtrage
h2  =  ones(1,round(Ns2));
x2 = filter(h2,1,Suite_diracs2);
figure ; plot(x2); title("x2");
figure ;  pwelch(x2,[],[],[],1/Te,'twosided'); title("DSP x2");

%DSP théorique 
DSP_theorique = Ts2*power(sinc(f*Ts2),2);
figure;
semilogy(f,DSP_theorique);
xlabel("Fréquence (Hz)");
title("DSP théorique x2");


%%  Modulateur 3    %%
Symboles3 = 2*bits-1;

M3 = 2;
Ts3 = log2(M3)*1/Rb;
Ns3 = Ts3/Te;

%Surechantillonage 
Suite_diracs3 = kron(Symboles3, [1 zeros(1,Ns3-1)]);

%Filtrage
h3  =  rcosdesign(0.5, 5, Ns3);
x3 = filter(h3,1,Suite_diracs3);
figure ; plot(x3) ; title("x3");
figure ;  pwelch(x3,[],[],[],1/Te,'twosided') ; title("DSP x3");

DSP_theorique=zeros(length(f));
v1=(1-0.6)/(2*Ts3);
v2=(1+0.6)/(2*Ts3);
for i=1:length(f)
    if  abs(f(i))<v1
        DSP_theorique(i)=var(x3);
    elseif abs(f(i))<v2
        DSP_theorique(i)=Ts3/2*(1+cos((pi*Ts3/0.6)*(abs(f(i))-v1)))*(var(x3)/Ts3);
    else 
        DSP_theorique(i)=0;
    end
end

figure;
semilogy(f,DSP_theorique);
xlabel("Fréquence (Hz)");
title("DSP théorique x3");

% Ensemble des DSP
figure;

dspx1 = pwelch(x1,[],[],[],1/Te,'twosided') ; 
dspx2 = pwelch(x2,[],[],[],1/Te,'twosided') ; 
dspx3 = pwelch(x3,[],[],[],1/Te,'twosided') ; title("DSP (les trois)");

semilogy(linspace(-Fe/2,Fe/2,length(dspx1)),dspx1, color = [1,0,0]);
hold on;
semilogy(linspace(-Fe/2,Fe/2,length(dspx1)),dspx2, color = [0,1,0]);
semilogy(linspace(-Fe/2,Fe/2,length(dspx1)),dspx3, color = [0,0,1]);
hold off;
