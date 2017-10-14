clc;
clear;
close all;

%% Problem Definition 

   global nTask;
   global R;
   global MaxSize;
   global nResource;
   global TaskSize;
    
     
%% problem definition  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nTask=32;
nResource=5;
MinSize=20;   %% minimum size of tasks
MaxSize=100;  %% maximum size of tasks
MinTimeOfResponse=1;     %% min Time of Response
MaxTimeOfResponse=5;    %% max Time of Response
MinSpeedOfResource=2;     %% min speed of resources
MaxSpeedOfResource=10;    %% max speed of resources
MinRelyOfSystem=20;   %% min Relaiability of System
MaxRelyOfSystem=100;  %% max Relaiability of System


R1=(MaxSpeedOfResource-MinSpeedOfResource).*rand(1,nResource)+MinSpeedOfResource;   %%speed of cpu
R1=sort(R1,'descend');
R2=(MaxSize-MinSize).*rand(1,nResource)+MinSize; 
R3=(MaxRelyOfSystem-MinRelyOfSystem).*rand(1,nResource)+MinRelyOfSystem;
R3=sort(R3,'descend');
R=[R1
   R2
   R3];
xMin=0;
xMax=4;
  


%% MOPSO Settings

nPop=100;   % Population Size

nRep=100;   % Repository Size

MaxIt=100;  % Maximum Number of Iterations

phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));

w=chi;              % Inertia Weight
wdamp=1;            % Inertia Weight Damping Ratio
c1=chi*phi1;        % Personal Learning Coefficient
c2=chi*phi2;        % Global Learning Coefficient

alpha=0.1;  % Grid Inflation Parameter

nGrid=10;   % Number of Grids per each Dimension

beta=4;     % Leader Selection Pressure Parameter

gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure

VelMax=0.8;

%% Initialization
for ii=1:1
    
TaskSize=round(unifrnd(1,MaxSize,[1 nTask]));
particle=CreateEmptyParticle(nPop);

for i=1:nPop
    particle(i).Velocity=0;
    particle(i).Position=xMin+(xMax-xMin)*(unifrnd(0,1,[1 nTask]));

    particle(i).Cost=MyCost1(particle(i).Position);

    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
end

particle=DetermineDomination(particle);

rep=GetNonDominatedParticles(particle);

rep_costs=GetCosts(rep);
G=CreateHypercubes(rep_costs,nGrid,alpha);

for i=1:numel(rep)
    [rep(i).GridIndex rep(i).GridSubIndex]=GetGridIndex(rep(i),G);
end
    
%% MOPSO Main Loop

for it=1:MaxIt
    for i=1:nPop
        rep_h=SelectLeader(rep,beta);

        particle(i).Velocity=w*particle(i).Velocity...
        +c1*rand*((particle(i).Best.Position)-(particle(i).Position))...
        +c2*rand*((rep_h.Position)-(particle(i).Position));
       

        particle(i).Velocity=min(max(particle(i).Velocity,-VelMax),+VelMax);

        particle(i).Position=particle(i).Position +particle(i).Velocity;

        flag=(particle(i).Position<1 | particle(i).Position>nResource);
        particle(i).Velocity(flag)=-particle(i).Velocity(flag);
        
        %particle(i).Position=min(max(particle(i).Position,-4),4);

        particle(i).Cost=MyCost1(particle(i).Position);

        if Dominates(particle(i),particle(i).Best)
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
        elseif ~Dominates(particle(i).Best,particle(i))
            if rand<0.5
                particle(i).Best.Position=particle(i).Position;
                particle(i).Best.Cost=particle(i).Cost;
            end
        end

    end
    
    particle=DetermineDomination(particle);
    nd_particle=GetNonDominatedParticles(particle);
    
    rep=[rep
         nd_particle];
    
    rep=DetermineDomination(rep);
    rep=GetNonDominatedParticles(rep);
    
    for i=1:numel(rep)
        [rep(i).GridIndex rep(i).GridSubIndex]=GetGridIndex(rep(i),G);
    end
    
    if numel(rep)>nRep
        EXTRA=numel(rep)-nRep;
        rep=DeleteFromRep(rep,EXTRA,gamma);
        
        rep_costs=GetCosts(rep);
        G=CreateHypercubes(rep_costs,nGrid,alpha);
        
    end
   
    disp(['Iteration ' num2str(it) ': Number of Repository Particles = ' num2str(numel(rep))]);
    
    w=w*wdamp;
end

%% Results

costs=GetCosts(particle);
rep_costs=GetCosts(rep);

    [Cost1 Index1]=min(rep_costs(1,:));
    [Cost2 Index2]=min(rep_costs(2,:));
  
    disp(' ');
  
disp(' ');
disp('best MakeSpan');
BestMakespan=rep(Index1,:);
Best_makespan(ii)=BestMakespan.Cost(1);
 disp(['Best MakeSpan : ' num2str(BestMakespan.Cost(1))]);
 disp(['makespan : ' num2str(BestMakespan.Cost(2))]);

disp(' ');
disp('best Reliability');
BestReliability=rep(Index2,:);
Best_Reliability(ii)=BestReliability.Cost(2);
 disp(['MakeSpan : ' num2str(BestReliability.Cost(1))]);
 disp(['Reliability : ' num2str(BestReliability.Cost(2))]);



end
% mean_Makespan=mean(Best_Makespan);
% mean_Reliability=mean(Best_Reliability);

figure(1);

%plot3(costs(1,:),costs(2,:),costs(3,:),'r*');
hold on;
grid on;
plot(rep_costs(1,:),rep_costs(2,:));
xlabel('f_1 : Makespan');
ylabel('f_2 : Reliability');
