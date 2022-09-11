close all;
clear all;

%tiledlayout(2,5);

%TP1 traitement du signal S6

Fe = 24e3;
Te = 1/Fe;
Rb = 3000;

%Génération d'un signal
nb_bits = 1000;
bits = randi([0,1],1,nb_bits);

%%  Modulation   %%
%Mappage 
Symboles1 = 2*bits-1;

M1 = 2;
Ts1 = log2(M1)*1/Rb;
Ns1 = Ts1/Te;

%Surechantillonage 
Suite_diracs1 = kron(Symboles1, [1 zeros(1,Ns1-1)]);

%Filtrage
h1  =  ones(1,round(Ns1));
%x1 = filter(h1,1,Suite_diracs1);

%%  démodulation   %%

z = filter(conv(h1,h1),1,Suite_diracs1);
nexttile;
plot(linspace(0,Te,length(z)),z); title("sortie filtre de z(t)"); %on retombe sur le signal en sortie de mappage (normal car il n'y a pas de bruit)
xlabel("temps(s)");

%Génération d'une impulsion de dirac

h1convh1 = conv(h1,h1);
FFTh1convh1 = fft(h1convh1,128);

nexttile;
plot(h1convh1); title("réponse impulsionnelle g"); % on choisit n0 = 8 (critère de Nynquist)
xlabel("indice");
nexttile;
plot(reshape([z z],[],length(z)/Ns1)); title("diagramme de l'oeil"); % on retrouve bien que tous les n0 notre signal ne peut prendre que deux valeurs différentes
xlabel("indice");

%échantillonage de z avec n0 = 8
z_echan = z([8:Ns1:length(z)])/Ns1;
bits_reconstruits = (sign(z_echan) + 1)/2;
erreur = length(find(bits ~= bits_reconstruits))/length(bits);
fprintf("\nValeur Taux d'erreur binaire avec n0 = 8 = %0.3e\n", erreur);

%échantillonage de z avec n0 = 3
z_echan = z([3:Ns1:length(z)])/Ns1;
bits_reconstruits = (sign(z_echan) + 1)/2;
erreur = length(find(bits ~= bits_reconstruits))/length(bits);
fprintf("\nValeur Taux d'erreur binaire avec n0 = 3 = %0.3e\n", erreur);

%étude avec canal de propagation sans bruit 1
fc = 8000;
hc = (2*fc/Fe)*sinc(2*(fc/Fe)*[-(100 - 1)/2 : (100 - 1)/2]);

g1 = conv(h1,conv(hc,h1));

nexttile;
plot(g1); title("réponse impulsionnelle g(t) fc = 8e3Hz");
xlabel("Indices");

z2 = filter(g1,1,Suite_diracs1);
nexttile;
plot(reshape([z2 z2],[],length(z2)/Ns1)); title("diagramme de l'oeil fc = 8e3Hz"); 
xlabel("indice");


nexttile;
plot(linspace(-Fe/2,Fe/2,length(FFTh1convh1)), abs(FFTh1convh1/max(FFTh1convh1)));
hold on;
plot(linspace(-Fe/2,Fe/2,length(fft(hc,128))), abs(fft(hc,128)/max(fft(hc,128)))); 
title("Rép fr Hc (o) et Hr.H (b) fc = 8e3Hz");
xlabel("fréquence(Hz)");
hold off;

z_echan = z2([58:Ns1:length(z2)])/Ns1; %50 + 8 decalage du aux filtres 
bits_reconstruits = (sign(z_echan) + 1)/2;
erreur_avec_cannal1 = length(find(bits(1:end-7) ~= bits_reconstruits))/length(bits);
fprintf("\nValeur Taux d'erreur binaire avec canal (fc = 8000Hz) = %0.3e\n", erreur_avec_cannal1);

%étude avec canal de propagation sans bruit 2
fc = 1000; 
hc = (2*fc/Fe)*sinc(2*(fc/Fe)*[-(100 - 1)/2 : (100 - 1)/2]);

g2 = conv(h1,conv(hc,h1));

nexttile;
plot(g2); title("Réponse impulsionnelle g(t) fc = e3Hz");
xlabel("Indices");

z2 = filter(g2,1,Suite_diracs1);
nexttile;
plot(reshape([z2 z2],[],length(z2)/Ns1)); title("diagramme de l'oeil");
xlabel("indice");

nexttile;
plot(linspace(-Fe/2,Fe/2,length(FFTh1convh1)), abs(FFTh1convh1/max(FFTh1convh1)));
hold on;
plot(linspace(-Fe/2,Fe/2,length(fft(hc,128))), abs(fft(hc,128)/max(fft(hc,128))));
title("Rép fr Hc (o) et Hr.H (b) fc = 1e3Hz"); 
xlabel("fréquence (Hz)");
hold off;

z_echan = z2([58:Ns1:length(z2)]); %50 + 8 decalage du aux filtres 
bits_reconstruits = (sign(z_echan) + 1)/2;
erreur_avec_cannal2 = length(find(bits(1:end-7) ~= bits_reconstruits))/length(bits);
fprintf("\nValeur Taux d'erreur binaire avec canal (fc = 1e3Hz) = %0.3e\n", erreur_avec_cannal2);





