%% Contingency Table:
%           Rev        Irrev
%Reg        19         7
%Non Reg    13         2

x = table([19; 13],[7; 2]);

%% Fisher test: use a right-tailed Fischer's exact test 
[h,p,stats] = fishertest(x,'Tail','right','Alpha',0.01)

%% Enrichment Analysis

% get data 

% logGamma = xlsread('CCMthermodynamicsAndRegulationContingency.xlsx',2,'O2:O42');
regulation = xlsread('CCMthermodynamicsAndRegulationContingency.xlsx',2,'P2:P42');

GSEA(regulation, length(regulation), true);

