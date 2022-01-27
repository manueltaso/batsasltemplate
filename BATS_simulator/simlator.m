%% BATS - Boston ASL simulator
clear all

% Data loading
pw=load_nii('BATS_asl_pw_brain_noN4.nii.gz');
pw=pw.img;

att=load_nii('BATS_asl_att_brain.nii.gz');
att=(att.img)./1000;
m0=load_nii('BATS_m0_brain_noN4.nii.gz');
m0=m0.img;
m0=imgaussfilt3(m0,4);

t1=load_nii('BATS_t1q_brain.nii.gz');
t1=(t1.img)./1000;
t1=imgaussfilt3(t1,4);

% Fixed parameters for CBF calculation
a=0.6;
lambda=0.9;
lblt=2.2;
PLD=1.8;
T1b=1.6;

% One compartment-models
CBF_1c=(6000*lambda*pw*exp(PLD/T1b))./(2*a*T1b*m0*(1-exp(-lblt/T1b)));
CBF_1c(isinf(CBF_1c)|isnan(CBF_1c))=0;
CBF_1c(CBF_1c>150)=0;

CBF_1c_att_t1b=(6000*(pw./m0).*((lambda./(2*a.*T1b))).*1./((exp(-(att/T1b))).*(exp(-(max(PLD-att,0))./T1b)-exp(-(max(lblt+PLD-att,0))./T1b))));
CBF_1c_att_t1b(isinf(CBF_1c_att_t1b)|isnan(CBF_1c_att_t1b))=0;
CBF_1c_att_t1b(CBF_1c_att_t1b>150)=0;

CBF_1c_att_t1t=(6000*(pw./m0).*((lambda./(2*a.*t1))).*1./((exp(-(att/T1b))).*(exp(-(max(PLD-att,0))./t1)-exp(-(max(lblt+PLD-att,0))./t1))));
CBF_1c_att_t1t(isinf(CBF_1c_att_t1t)|isnan(CBF_1c_att_t1t))=0;
CBF_1c_att_t1t(CBF_1c_att_t1t>150)=0;

%% PLD simulation - one compartment - T1 tissue = T1b
 
clear dMtmp dM

PLDlist=[0.5:0.1:4];        % Range of PLDs to simulate
lblt=2.2;                   % Fixed single-labeling duration

for i=1:length(PLDlist)
     
    PLD=PLDlist(i)
    dMtmp=2.*a.*m0.*CBF_1c_att_t1b.*T1b.*exp(-att./T1b).*(exp(-(max(PLD-att,0)./T1b))-exp(-(max(lblt+PLD-att,0)./T1b)))./(lambda*6000);
    dM1ct1b(:,:,:,i)=permute(dMtmp,[2 1 3]);
    
end


%% PLD simulation - one compartment - T1 tissue from T1 map

clear dMtmp dM

PLDlist=[0.5:0.1:4];        % Range of PLDs to simulate
lblt=2.2;                   % Fixed single-labeling duration

for i=1:length(PLDlist)
     
    PLD=PLDlist(i)
    dMtmp=2.*a.*m0.*CBF_1c_att_t1t.*t1.*exp(-att./T1b).*(exp(-(max(PLD-att,0)./t1))-exp(-(max(lblt+PLD-att,0)./t1)))./(lambda*6000);
    dM1ct1t(:,:,:,i)=permute(dMtmp,[2 1 3]);
    
end


