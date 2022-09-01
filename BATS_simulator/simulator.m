%% Boston ASL Template and Simulator (BATS)

clear all

% Data loading
load t1.mat 
load m0.mat 
load mask.mat 
load pw.mat 
load att.mat

% Fixed parameters for CBF calculation
a=0.6;
lambda=0.9;
lblt=2.2;
PLD=1.8;
t1b=1.6;

% One compartment-model (Alsop et al MRM 2015 White Paper)
f1c=(6000*lambda*pw*exp(PLD/t1b))./(2*a*t1b*m0*(1-exp(-lblt/t1b)));
f1c(isinf(f1c)|isnan(f1c))=0;
f1c(f1c>150)=0;
f1c=f1c.*mask;

% Two-compartment model - Fractional extra-ATT (From Alsop and Detre 1996,
% modified for finite labeling bolus)

att2=att+att.*0.6;
t1b=1.6;
alpha=0.6;
lambda=0.9;
pld=1.8;
lbl=2.2;

E1=exp(-att2.*((1/t1b)-(1./t1)));
E2=exp(-(max(pld,att2))./t1);
E3=exp(-(max(pld+lbl,att2))./t1);

A=t1.*E1.*(E2-E3);

E4=exp(-(max(pld,att)/t1b));
E5=exp(-(max(pld+lbl,att)/t1b));
E6=exp(-(max(pld,att2)/t1b));
E7=exp(-(max(pld+lbl,att2)/t1b));

B=t1b.*((E4-E5)-(E6-E7));

f2c_frac=(pw./(m0.*(A+B))).*(lambda/(2.*alpha))*6000;

f2c_frac(isinf(f2c_frac)|isnan(f2c_frac))=0;
f2c_frac=abs(f2c_frac);
f2c_frac(f2c_frac>200)=0;
f2c_frac=f2c_frac.*mask;

% Two-compartment model - Fixed extra-ATT (From Alsop and Detre 1996,
% modified for finite labeling bolus)

att2=att+1;
t1b=1.6;
alpha=0.6;
lambda=0.9;
pld=1.8;
lbl=2.2;

E1=exp(-att2.*((1/t1b)-(1./t1)));
E2=exp(-(max(pld,att2))./t1);
E3=exp(-(max(pld+lbl,att2))./t1);

A=t1.*E1.*(E2-E3);

E4=exp(-(max(pld,att)/t1b));
E5=exp(-(max(pld+lbl,att)/t1b));
E6=exp(-(max(pld,att2)/t1b));
E7=exp(-(max(pld+lbl,att2)/t1b));

B=t1b.*((E4-E5)-(E6-E7));

f2c_fix=(pw./(m0.*(A+B))).*(lambda/(2.*alpha))*6000;

f2c_fix(isinf(f2c_fix)|isnan(f2c_fix))=0;
f2c_fix=abs(f2c_fix);
f2c_fix(f2c_fix>200)=0;
f2c_fix=f2c_fix.*mask;

%% PLD simulation - one compartment - T1 tissue = t1b

clear dm1ctmp dm1c

PLDlist=[0.5:0.2:2.5];        % Range of PLDs to simulate
lblt=1.8;                   % Fixed single-labeling duration

for i=1:length(PLDlist)
    
    PLD=PLDlist(i)
    dm1ctmp=2.*a.*m0.*f1c.*t1b.*exp(-att./t1b).*(exp(-(max(PLD-att,0)./t1b))-exp(-(max(lblt+PLD-att,0)./t1b)))./(lambda*6000);
    dm1c(:,:,:,i)=permute(dm1ctmp,[2 1 3]);
    
end

%% PLD simulation - two-comp - Fractional ATT

att2=att+att.*0.6;
t1b=1.6;
alpha=0.6;
lambda=0.9;
lbl=2.2;

PLDlist=[0.5:0.1:4];        % Range of PLDs to simulate

for i=1:length(PLDlist)
    
    pld=PLDlist(i)
    
    E1=exp(-att2.*((1/t1b)-(1./t1)));
    E2=exp(-(max(pld,att2))./t1);
    E3=exp(-(max(pld+lbl,att2))./t1);
    A=t1.*E1.*(E2-E3);
    
    E4=exp(-(max(pld,att)/t1b));
    E5=exp(-(max(pld+lbl,att)/t1b));
    E6=exp(-(max(pld,att2)/t1b));
    E7=exp(-(max(pld+lbl,att2)/t1b));
    B=t1b.*((E4-E5)-(E6-E7));
    
    dm2c_frac(:,:,:,i)=(f2c_frac.*m0.*(A+B).*2.*a)./(6000*lambda);   
    
end

%% PLD simulation - two-comp - Fixed extra-ATT

att2=att+1;
t1b=1.6;
alpha=0.6;
lambda=0.9;
lbl=2.2;

PLDlist=[0.5:0.1:4];        % Range of PLDs to simulate

for i=1:length(PLDlist)
    
    pld=PLDlist(i)
    
    E1=exp(-att2.*((1/t1b)-(1./t1)));
    E2=exp(-(max(pld,att2))./t1);
    E3=exp(-(max(pld+lbl,att2))./t1);
    A=t1.*E1.*(E2-E3);
    
    E4=exp(-(max(pld,att)/t1b));
    E5=exp(-(max(pld+lbl,att)/t1b));
    E6=exp(-(max(pld,att2)/t1b));
    E7=exp(-(max(pld+lbl,att2)/t1b));
    B=t1b.*((E4-E5)-(E6-E7));
    
    dm2c_fix(:,:,:,i)=(f2c_fix.*m0.*(A+B).*2*a)./(6000*lambda);   
end

%% Paper examples

% PLD simulation - one compartment - globally prolonged ATT

PLDlist=[0.5:0.2:2.5];        % Range of PLDs to simulate
lblt=1.8;                   % Fixed single-labeling duration
attb=att.*1.5;               % Globally prolonged ATT

for i=1:length(PLDlist)
    
    PLD=PLDlist(i)
    dm1ctmp=2.*a.*m0.*f1c.*t1b.*exp(-attb./t1b).*(exp(-(max(PLD-attb,0)./t1b))-exp(-(max(lblt+PLD-attb,0)./t1b)))./(lambda*6000);
    dm1c(:,:,:,i)=permute(dm1ctmp,[2 1 3]);
    
end

figure
colormap(jet)
subplot(311);
imagesc([dm1c(:,:,55,1) dm1c(:,:,55,2) dm1c(:,:,55,3) dm1c(:,:,55,4) dm1c(:,:,55,5) dm1c(:,:,55,6) dm1c(:,:,55,7) dm1c(:,:,55,8) dm1c(:,:,55,9) dm1c(:,:,55,10)])
subplot(312);
imagesc([dm1c(:,:,73,1) dm1c(:,:,73,2) dm1c(:,:,73,3) dm1c(:,:,73,4) dm1c(:,:,73,5) dm1c(:,:,73,6) dm1c(:,:,73,7) dm1c(:,:,73,8) dm1c(:,:,73,9) dm1c(:,:,73,10)])
subplot(313);
imagesc([dm1c(:,:,120,1) dm1c(:,:,120,2) dm1c(:,:,120,3) dm1c(:,:,120,4) dm1c(:,:,120,5) dm1c(:,:,120,6) dm1c(:,:,120,7) dm1c(:,:,120,8) dm1c(:,:,120,8) dm1c(:,:,120,10)])

% PLD simulation - one compartment - local prolonged ATT

PLDlist=[0.5:0.2:2.5];        % Range of PLDs to simulate
lblt=1.8;                   % Fixed single-labeling duration

% Mask simulating a local impairment leading to prolonged ATT (e.g. CVD)

testatt=att;
maskatt=load_nii('attmask.nii.gz');
maskatt=double(maskatt.img);
maskatt=imgaussfilt(maskatt,4);
testatt(maskatt~=0)=testatt(maskatt~=0).*1.5;
attc=imgaussfilt3(testatt,4);

for i=1:length(PLDlist)
    
    PLD=PLDlist(i)
    dm1ctmp=2.*a.*m0.*f1c.*t1b.*exp(-attc./t1b).*(exp(-(max(PLD-attc,0)./t1b))-exp(-(max(lblt+PLD-attc,0)./t1b)))./(lambda*6000);
    dm1c(:,:,:,i)=permute(dm1ctmp,[2 1 3]);
    
end

figure
colormap(jet)
subplot(311);
imagesc([dm1c(:,:,55,1) dm1c(:,:,55,2) dm1c(:,:,55,3) dm1c(:,:,55,4) dm1c(:,:,55,5) dm1c(:,:,55,6) dm1c(:,:,55,7) dm1c(:,:,55,8) dm1c(:,:,55,9) dm1c(:,:,55,10)])
subplot(312);
imagesc([dm1c(:,:,73,1) dm1c(:,:,73,2) dm1c(:,:,73,3) dm1c(:,:,73,4) dm1c(:,:,73,5) dm1c(:,:,73,6) dm1c(:,:,73,7) dm1c(:,:,73,8) dm1c(:,:,73,9) dm1c(:,:,73,10)])
subplot(313);
imagesc([dm1c(:,:,120,1) dm1c(:,:,120,2) dm1c(:,:,120,3) dm1c(:,:,120,4) dm1c(:,:,120,5) dm1c(:,:,120,6) dm1c(:,:,120,7) dm1c(:,:,120,8) dm1c(:,:,120,8) dm1c(:,:,120,10)])


%%
figure ; subplot(131) ; imagesc(flipud(f1c(:,:,72)')) ; subplot(132) ; imagesc(flipud(f2c_frac(:,:,72)')) ; subplot (133) ; imagesc(flipud(f2c_fix(:,:,72)'))
