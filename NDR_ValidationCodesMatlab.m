function [Statistics, Points] = NDR_ValidationCodesMatlab(CombinedExportData)
BootMax = 10000;
BootPoints = 100;


Test = get(CombinedExportData);
%Numbers_to_do = [9,10,11];  %Nitrogen
%ValiNumber = 7; %Nitrogen
%HectaNumber = 5;
 
%Numbers_to_do = [7, 8]; % Water
%ValiNumber = 5; % Water
%HectaNumber =6; % Water

Numbers_to_do = [14]; % Phosporus
ValiNumber = 8; % Phosporus
HectaNumber = 5;

Names = cellstr(Test.VarNames(1,Numbers_to_do));
Offset = Numbers_to_do(1)-1;
Statistics = dataset(zeros(length(Names),1),'ObsNames',Names,'varnames',{'InverseDeviance'});
Points.Winsor = dataset(CombinedExportData.DH_ID,'varnames',{'DH_ID'});
Points.Deviance = Points.Winsor;
Points.NonWinsor = Points.Winsor;
Points.Ratio =  Points.Winsor;
Points.RankingDeviation = Points.Winsor;

% Validation data
Hectares = double(CombinedExportData(:,HectaNumber));
InVarVali = double(CombinedExportData(:,ValiNumber))./Hectares; 
InVarValiLog = log10(InVarVali+1); 
Points.NonWinsor.EIONValidation = InVarVali;
[TestArray(:,1)] = (WinsorFunction(((InVarValiLog)),2.5));
Points.Winsor.EIOValidation = TestArray(:,1);
for i = 1:1:length(Numbers_to_do)
    InVar = double(CombinedExportData(:,i+Offset))./Hectares; %#ok<*NODEF>
    InVarLog = log10(InVar+1); %#ok<*NODEF>
    [TestArray(:,2)] = (WinsorFunction(((InVarLog)),2.5));
    Sheds = Points.Winsor.DH_ID;
    [Outputs] = Accuracy_statistics_Validation(TestArray,Sheds);
    Points.Deviance.(genvarname(char(Names(i)))) = reshape(Outputs.deviation_point,[],1);
    Points.Ratio.(genvarname(char(Names(i)))) = reshape(Outputs.Ratio_point,[],1);
    Points.RankingDeviation.(genvarname(char(Names(i)))) = reshape(Outputs.RankingDeviation,[],1);
    save("all.mat")
    Points.NonWinsor.(genvarname(char(Names(i)))) = InVar;
    Points.Winsor.(genvarname(char(Names(i)))) = TestArray(:,2);

    SelectedPoints = randi(length(InVar),BootPoints,BootMax);
    clear Outputs
    parfor boot = 1:BootMax
        Sheds = Points.Winsor.DH_ID(SelectedPoints(:,boot),1); 
        BootVali = TestArray(SelectedPoints(:,boot),1); %#ok<PFBNS>
        BootVar = TestArray(SelectedPoints(:,boot),2);
        BootArray =[BootVali,BootVar];
        [Outputs] = Accuracy_statistics_Validation(BootArray,Sheds);
        StatisticsbootInverseDeviance(boot,1) = Outputs.mean_double_deviation;
        StatisticsbootRho(boot,1) = Outputs.RHO;
        StatisticsbootRhoPVal(boot,1) = Outputs.PVAL;
        StatisticsbootMedianValueWinsor(boot,1) = nanmedian(BootArray(:,1));
        StatisticsbootMedianValueValiWinsor(boot,1) =  nanmedian(BootArray(:,2));
    end
     Statistics.InverseDeviance(i,1) = nanmedian(StatisticsbootInverseDeviance);
     Statistics.InverseDevianceSTD(i,1) = nanstd(StatisticsbootInverseDeviance);
     Statistics.Rho(i,1) = nanmedian(StatisticsbootRho);
     Statistics.RhoSTD(i,1) = nanstd(StatisticsbootRho);
     Statistics.RhoPVal(i,1) = nanmedian(StatisticsbootRhoPVal);
     Statistics.RhoPValSRD(i,1) = nanstd(StatisticsbootRhoPVal);
     Statistics.MedianValue(i,1) = nanmedian(InVar);
     Statistics.MedianValueVali(i,1) =  nanmedian(InVarVali);
     Statistics.Perc05Value(i,1) = prctile(InVar,5);
     Statistics.Perc95Value(i,1) = prctile(InVar,95);
     Statistics.Perc05ValueVali(i,1) = prctile(InVarVali,5);
     Statistics.Perc95ValueVali(i,1) = prctile(InVarVali,95);
     Statistics.MedianValueWinsor(i,1) = nanmedian(StatisticsbootMedianValueWinsor);
     Statistics.MedianValueValiWinsor(i,1) =  nanmedian(StatisticsbootMedianValueValiWinsor);
     Statistics.MedianValueWinsorSTD(i,1) = nanstd(StatisticsbootMedianValueWinsor);
     Statistics.MedianValueValiWinsorSTD(i,1) =  nanstd(StatisticsbootMedianValueValiWinsor);
end
end

function [Outputs] = Accuracy_statistics_Validation(testArray,Sheds)
% clean data set to determine true N and where top put NaN;
testArray(isinf(testArray)==1) = NaN;
[testArray]  = CleanOutNaN(testArray);

%% Make log if applicable, but not for Ensembles
testVar = testArray;
Outputs = CorrFunc(testVar,Sheds);

end

function  Outputs = CorrFunc(testSet,Sheds)
    Precision = 1./(0.00001);
    c1=find((isnan(testSet(:,1))==1));
    d=find((isnan(testSet(:,2))==1));
    alL = [c1;d];
    A(:,1) = unique(alL);
    testSet(A,:) = [];  %#ok<*FNDSB>
    Datapoint = size(testSet,1);
    [Outs.RHO,Outs.PVAL] = corr(testSet(:,1),testSet(:,2),'type','Spearman');
    Outputs.RankingDeviation = RankingDeviation(testSet,Sheds);
    Outputs.PVAL = single(Outs.PVAL);
    Outputs.RHO = (round(Outs.RHO.*Precision))./Precision;

    % Inverse deviance against a 1:1 line
    clear x_range y_range
    x_range = testSet(:,1);
    y_range = testSet(:,2);
    Deviation_point = abs(y_range-x_range);
    RatioPoint = abs(y_range./x_range);
    Outputs.deviation_point= Deviation_point; %% Accuracy per point
    Outputs.Ratio_point= RatioPoint; 
    Outputs.mean_double_deviation = 1- (nansum(Deviation_point)/Datapoint); %%Accuracy overall
    Outputs.mean_double_deviation = (round( Outputs.mean_double_deviation.*Precision))./Precision;
end

function RankingDevOut = RankingDeviation(Input,Sheds)
    Testset(:,1) = Sheds;
    Testset(:,2) = Input(:,1);
    Testset(:,3) = Input(:,2);
    Testset = sortrows(Testset,2);
    Testset(:,4) = 1:length(Sheds);
    Testset = sortrows(Testset,3);
    Testset(:,5) = 1:length(Sheds);
    Testset = sortrows(Testset,1);
    RankingDevOut = (abs(Testset(:,4) - Testset(:,5)))./length(Sheds);

end

function [ArrayOut] = CleanOutNaN(ArrayIn)
OutVarOrg = 1:1:length(ArrayIn(:,1));
OutVarOrg = OutVarOrg';
a1=find((isnan(ArrayIn(:,1))==1));
b=find((isnan(ArrayIn(:,2))==1));
all = [a1;b];
a(:,1) = unique(all);
ArrayOut = ArrayIn;
ArrayOut(a,:) = [];  
end

function  [OutVar] = WinsorFunction(InVar,percLow)
InVar = reshape(InVar,((size(InVar,1)).*(size(InVar,2))),1);
InVar(InVar<0) = NaN;
prct =  prctile(InVar,percLow);
if prct < 0
    display('zero Values present; correct this first')
    cccc
end
InVar_norm = InVar - prct;
clear InVar
InVar_norm( InVar_norm<0) = 0;
prct(2) = prctile(InVar_norm,(100-percLow));
if exist('InVar_org') == 1; %#ok<EXIST>
    InVar_norm = InVar_org;
end
OutVar = (InVar_norm./prct(2));
OutVar = OutVar + 0.0001;
clear InVar_norm InVar_org testlist testtmp upboud
OutVar(OutVar>1) = 1;
end

