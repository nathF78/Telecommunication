close all;
clear all;

tiledlayout(3,4);

%TP1 traitement du signal S6

Fe = 24e3;
Te = 1/Fe;
Rb = 3000;
%Génération d'un signal
nb_bits = 10000;
bits = randi([0,1],1,nb_bits);

%%  Modulation   %%
%Mappage 
Symboles1 = 2*bits-1;

M1 = 2;
Ts1 = log2(M1)*1/Rb;
Ns1 = Ts1/Te;

%Surechantillonage 
Suite_diracs1 = kron(Symboles1, [1 zeros(1,Ns1-1)]);

%% Chaine de Référence

%Filtrage
h0  =  ones(1,Ns1);
x0 = filter(h0,1,Suite_diracs1);

%  démodulation   %
znonbruite0 = filter(conv(h0,h0),1,Suite_diracs1);

nexttile;
plot(reshape([znonbruite0 znonbruite0],[],length(znonbruite0)/Ns1));
title("diagramme de l'oeil chaine de référence"); %On prend n0 = 8
n00 = 8;

% ajout bruit 
Px=mean(abs(x0).^2);
Eb_N0 = [0:1:8];
teb0 = [];
teb_th = [];

for i = 1:length(Eb_N0) 
    Eb_N0_i = Eb_N0(i);
    Eb_N0_i = 10^(Eb_N0_i/10);
    sigma = sqrt((Px*Ns1) / (2*log2(M1)*Eb_N0_i)) ;

    bruit = sigma * randn(1,length(x0));
    x0_bruite = x0 + bruit;

    z0_bruite = filter(fliplr(h0),1,x0_bruite);
    z0_echan_bruite = z0_bruite([n00:Ns1:length(z0_bruite)]);
    
    z0_sign = sign(z0_echan_bruite);
    bits_reconstruits = (z0_sign + 1)/2;
    
    teb_i = length(find(bits ~= bits_reconstruits))/length(bits);
    teb0 = [teb0 teb_i];
    
    teb_th = [teb_th qfunc(sqrt(2*Eb_N0_i))/log2(M1)];

end

nexttile;
semilogy(Eb_N0, teb0);
hold on;
semilogy(Eb_N0, teb_th);
xlabel('Eb/N0');
legend('TEB (b)','TEB théorique (o)');
title('TEB chaine de référence');


%% 1er chaine à étudier 

%Filtrage
h1  =  ones(1,Ns1);
x1 = conv(Suite_diracs1,h1);

%  démodulation   %
hr1  = ones(1,Ns1/2);
znonbruite1 = filter(conv(h1,hr1),1,Suite_diracs1);
%znonbruite1 = [znonbruite1 0 0 0 0 0 0];


%étude sans bruit 5.3.1
nexttile;
plot(reshape([znonbruite1 znonbruite1],[],length(znonbruite1)/Ns1)); % on choisit n0 = 8 (critère de Nynquist)
title("diagramme de l'oeil avec canal 1er chaine"); 

z1_echan = znonbruite1([8:Ns1:length(znonbruite1)]);
z_sign = sign(z1_echan);
bits_reconstruits = (z_sign + 1)/2;

erreur_chaine1_sans_bruit = length(find(bits ~= bits_reconstruits))/length(bits);
fprintf("\nValeur Taux d'erreur binaire avec canal 1 sans bruit = %0.3e\n", erreur_chaine1_sans_bruit);

%étude avec bruit 5.3.2

%calcul des TEB pour différents rapports 

Px = mean(abs(x1).^2);
teb1 = [];

for i = 1:length(Eb_N0) 
    Eb_N0_i = Eb_N0(i);
    Eb_N0_i = 10^(Eb_N0_i/10);
    sigma = sqrt((Px*Ns1) / (2*log2(M1)*Eb_N0_i)) ;

    bruit = sigma * randn(1,length(x1));
    x1_bruite = x1 + bruit;

    z1_bruite = filter(fliplr(hr1),1,x1_bruite);
    z1_echan_bruite = z1_bruite([n00:Ns1:length(z1_bruite)]);

    z1_sign = sign(z1_echan_bruite);
    bits_reconstruits = (z1_sign + 1)/2;
    
    teb_i = length(find(bits ~= bits_reconstruits))/length(bits);
    teb1 = [teb1 teb_i];
end

nexttile;
semilogy(Eb_N0, teb1)
hold on 
semilogy(Eb_N0, teb_th)
xlabel('Eb/N0')
legend('TEB (b)','TEB théorique (o)')
title('TEB chaine canal 1er chaine')

%5.3.2.4

nexttile;
semilogy(Eb_N0, teb1);
hold on ;
semilogy(Eb_N0, teb0);
xlabel('Eb/N0');
legend('TEB (b)','TEB chaine de ref (o)');
title('TEB 1er chaine');

%5.3.2.5
nexttile;
pwelch(conv(h0,h0),[],[],[],1/Te,'twosided') ; 
title("DSP chaine de référence");

nexttile;
pwelch(conv(h1,hr1),[],[],[],1/Te,'twosided') ; 
title("DSP 1er chaine étudiée");

%chaine de transmission 1 plus éfficace car la puissance est répartie sur
%moins de fréquences (donc moins d'informations perturbées par le bruit)

%% 5.4 Deuxième chaine à étudier 

% Modulation
%Mappage

Symboles2 = (2*bi2de(reshape(bits, 2, length(bits)/2).') - 3).';
M2 = 4;
Ts2 = log2(M2)*1/Rb;
Ns2 = Ts2/Te;

%Surechantillonage 
Suite_diracs2 = kron(Symboles2, [1 zeros(1,Ns2-1)]);
Suite_diracs2 = reshape(Suite_diracs2', 1, []);

%Filtrage
h2  =  ones(1,Ns2);
x2 = conv(h2, Suite_diracs2);
%znonbruite2 = conv(h2,x2);
znonbruite2 = filter(conv(h2,h2),1,Suite_diracs2);

%diagramme de l'oeil

nexttile;
plot(reshape([znonbruite2 znonbruite2],[],length(znonbruite2)/Ns2));
title("diagramme de l'oeil avec canal 2eme chaine"); %On prend n02 = 16
n02 = 16;

%échantillonage de z avec n02 = 16
z2_echan = floor(znonbruite2([n02:Ns2:length(znonbruite2)])/Ns2);

BitsDecides2_sansbruit = reshape(de2bi((z2_echan + 3)/2).', 1, length(bits));

erreur_chaine2_sans_bruit = length(find(abs(sign(BitsDecides2_sansbruit) - bits) > 0))/length(bits);
fprintf("\nValeur Taux d'erreur binaire avec canal 2 sans bruit = %0.3e\n", erreur_chaine2_sans_bruit);

%calcul des TEB pour différents rapports 

Px = mean(abs(x2).^2);
teb2 = [];
tes2 = [];
tes_th = [];
teb_th = [];

for i = 1:length(Eb_N0) 
    Eb_N0_i = Eb_N0(i);
    Eb_N0_i = 10^(Eb_N0_i/10);
    sigma = sqrt((Px*Ns2) / (2*log2(M2)*Eb_N0_i)) ;

    bruit = sigma * randn(1,length(x2));
    x2_bruite = x2 + bruit;

    z2_bruite = filter(h2,1,x2_bruite);
    z2_echan_bruite = floor(z2_bruite([n02:Ns2:length(z2_bruite)])/Ns2);
    z2_echan_bruite(find(0 <= z2_echan_bruite & z2_echan_bruite <= 2)) = 1;
    z2_echan_bruite(find(0 > z2_echan_bruite & z2_echan_bruite >= -2)) = -1;
    z2_echan_bruite(find(z2_echan_bruite > 2)) = 3;
    z2_echan_bruite(find(z2_echan_bruite < -2)) = -3;
    bits_reconstruits = reshape(de2bi((z2_echan_bruite + 3)/2).',1,length(bits));

    tes_i = length(find(Symboles2 ~= z2_echan_bruite))/length(Symboles2);
    tes2 = [tes2 tes_i];
    tes_th_i = (3/2)*qfunc(sqrt((4/5)*Eb_N0_i));
    tes_th = [tes_th tes_th_i];
    
    teb_i = length(find(bits ~= bits_reconstruits))/length(bits);
    teb2 = [teb2 teb_i];
    teb_th_i = (3/4)*qfunc(sqrt((4/5)*Eb_N0_i));
    teb_th = [teb_th teb_th_i];
end

nexttile;
semilogy(Eb_N0, tes2);
hold on ;
semilogy(Eb_N0, tes_th);
xlabel('Eb/N0');
legend('TES (b)','TES théorique (o)');
title(['TES 2eme chaine']);

nexttile;
semilogy(Eb_N0, teb2);
hold on ;
semilogy(Eb_N0, teb0);
xlabel('Eb/N0');
legend('TEB (b)','TEB chaine de ref (o)');
title('TEB 2eme chaine');

nexttile;
pwelch(conv(h2,h2),[],[],[],1/Te,'twosided') ; 
title("DSP 2eme chaine étudiée");



