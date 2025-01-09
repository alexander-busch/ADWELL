
% Cleaning data

%% Get rid of zeros
Hole(Hole==0)=[];
Pipe (Pipe==0)=[];
Length (Length==0)=[];
Inclination (Inclination==0)=[];
Eccentricity_R (Eccentricity_R==0)=[];
Eccentricity_Theta (Eccentricity_Theta==0)=[];
Flow (Flow==0)=[];
YP (YP==0)=[];
PV (PV==0)=[];
FluidDens (FluidDens==0)=[];
PartDia (PartDia==0)=[];
PartDens (PartDens==0)=[];
Porosity (Porosity==0)=[];
BotHolTemp (BotHolTemp==0)=[];
rpm (rpm==0)=[];
ROP (ROP==0)=[];

% Variables = who;



%% Convert to SI & precalculations
D_H=(Hole-Pipe).*25.4./1000;
CrosSec=pi.*((Hole.*25.4./1000).^2-(Pipe.*25.4./1000).^2)./4;
PV=PV.*1e-3;
Flow = Flow./1000./60;

BulkVel=kron(Flow,1./CrosSec);
BulkVel = [min(BulkVel); mean(BulkVel); max(BulkVel)];

KinVisc = kron(FluidDens,1./PV);
KinVisc = [min(KinVisc); mean(KinVisc); max(KinVisc)];

Re_pre = kron(BulkVel,KinVisc);
Re_pre = [min(Re_pre); mean(Re_pre); max(Re_pre)];
