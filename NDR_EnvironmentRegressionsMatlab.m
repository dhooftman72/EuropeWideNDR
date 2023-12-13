function Regressions = NDR_EnvironmentRegressionsMatlab(CombinedData)
btsfact = 2.5; %2.5 degrees maximum interaction, else autocorrelation coefficient = 0
[wij] = Morans(CombinedData.Longitude,CombinedData.Lattitude,btsfact);
save('wij.mat','wij');

Transfer.DFtarget = 225;
Transfer.Bootmax = 10000;
Transfer.Size = size(CombinedData,1);
Transfer.NameList = [{'Slope_Terrain'};{'Annual_Rain'};{'Erosivity'};{'Prop_BaseFlow'};{'WyieldM3'};{'N_CV'};...
            {'P_CV'};{'PerWoodland'};{'PerCropLand'};{'PerImprovedGrassland'};{'PerNotImprovedGrass'};...
            {'PerUrbanAtrificial'};{'PerHighFerti'};{'PerLowFerti'}];
for j = 1:6
    warning off
    if j == 1
        VarIn = double(CombinedData.N_Deviance);
    elseif j == 2
        VarIn = double(CombinedData.N_Ratio_Deviance);
    elseif j == 3
        VarIn = double(CombinedData.N_Ranking_Deviance);
    elseif j == 4
        VarIn = double(CombinedData.P_Deviance);
    elseif j == 5
        VarIn = double(CombinedData.P_Ratio_Deviance);
    elseif j == 6
        VarIn = double(CombinedData.P_Ranking_Deviance);
    end
    Transfer.Type = j;
    AutoCorrelationCoef = AutoCorFunc(VarIn,Transfer.Size);
    [AnovaTable] = InteractionModel(VarIn,CombinedData,AutoCorrelationCoef,Transfer);
    if j == 1
        Regressions.N_Deviance = AnovaTable;
    elseif j == 2 
        Regressions.N_Ratio_Deviance = AnovaTable;
    elseif j == 3 
        Regressions.N_Ranking_Deviance = AnovaTable;
    elseif j == 4
        Regressions.P_Deviance = AnovaTable;
    elseif j == 5 
        Regressions.P_Ratio_Deviance = AnovaTable;
    elseif j == 6 
        Regressions.P_Ranking_Deviance = AnovaTable;
    end
end
end

%%
function [Anova_Table] =  InteractionModel(VarIn,Variables,AutoIn,Transfer)
Transfer.NrVariables = length(Transfer.NameList);
Anova_Table = dataset(zeros(Transfer.NrVariables,1),'ObsNames',Transfer.NameList,'varnames',{'SumSqAuto'});

Variable(:,1)  = double(Variables.Slope_Terrain);
Variable(:,2)  = log10(double(Variables.Annual_Rain)+1);
Variable(:,3) =   log10(double(Variables.Erosivity)+1);   
Variable(:,4) =  asin(sqrt((double(Variables.Prop_BaseFlow))));   
Variable(:,5)  = log10(double(Variables.WyieldM3)+1);
Variable(:,6)  = double(Variables.N_CV);
Variable(:,7)  = double(Variables.P_CV);
Variable(:,8) =  asin(sqrt((double(Variables.PerWoodland)))); 
Variable(:,9) =  asin(sqrt((double(Variables.PerCropLand)))); 
Variable(:,10) =  asin(sqrt((double(Variables.PerImprovedGrassland)))); 
Variable(:,11) =  asin(sqrt((double(Variables.PerNotImprovedGrass)))); 
Variable(:,12) =  asin(sqrt((double(Variables.PerUrbanAtrificial)))); 
Variable(:,13) =  asin(sqrt((double(Variables.PerHighFerti)))); 
Variable(:,14) =  asin(sqrt((double(Variables.PerLowFerti)))); 

for NrVariable = 1: (Transfer.NrVariables)
    display(NrVariable)
    display(Transfer.Type)
    parfor boot = 1:1:Transfer.Bootmax
        % pick random selection
        tmper = randperm(Transfer.Size);
        tmperList = sort(tmper(1:Transfer.DFtarget));
        VarInRun = VarIn(tmperList);
        AutoRun = AutoIn(tmperList);
        VariableRun = Variable(tmperList,NrVariable);

        % [~,outstmp,statstmp] = anovan(VarInRun,AutoRun,'sstype',1,...
        %     'model',[1],'continuous', [1],'display', 'off','varnames', {'AutoCorrelation'});
        % ATableSumSqAuto(boot,1) = outstmp(2,2);
        % ATableDFsAuto(boot,1) = outstmp(2,3);
        % ATableSumSqAutoError(boot,1) = outstmp(3,2);
        % ATableDFsAutoError(boot,1) = outstmp(3,3);
        % SecondTierVar = statstmp.resid';

        SecondTierVar = VarInRun;
        [~,outs,stats] = anovan(SecondTierVar,VariableRun,'sstype',1,...
                'model',[1],'continuous', [1], 'display', 'off',...
                'varnames', {'Variable'});
        ATableSumSq(boot,1) = outs(2,2);
        ATableDFs(boot,1) = outs(2,3);
        ATableSumSqError(boot,1) = outs(3,2);
        ATableDFsError(boot,1) = outs(3,3);
        coeffs(boot,1) = stats.coeffs(2);
        Points_Corrected =  reshape(SecondTierVar,[],1) - reshape(stats.resid,[],1);    
        RsquaredLoop(boot,1) = RsquaredFunc(Points_Corrected,reshape(SecondTierVar,[],1));
    end

    % Anova_Table.SumSqAuto(NrVariable,1) = nanmedian(cell2mat(ATableSumSqAuto));
    % Anova_Table.DFAuto(NrVariable,1) = nanmedian(cell2mat(ATableDFsAuto));
    % Anova_Table.SumSqAutoError(NrVariable,1) = nanmedian(cell2mat(ATableSumSqAutoError));
    % Anova_Table.DFAutoError(NrVariable,1) = nanmedian(cell2mat(ATableDFsAutoError));
    % 
    % Anova_Table.MeanSqAuto(NrVariable,1) = Anova_Table.SumSqAuto(NrVariable,1)./Anova_Table.DFAuto(NrVariable,1); 
    % Anova_Table.MeanSqAutoError(NrVariable,1) = Anova_Table.SumSqAutoError(NrVariable,1)./Anova_Table.DFAutoError(NrVariable,1); 
    % 
    % F_ValueAuto = Anova_Table.MeanSqAuto(NrVariable,1)./Anova_Table.MeanSqAutoError(NrVariable,1);
    % P_ValueAuto =  1- fcdf(F_ValueAuto,Anova_Table.DFAuto(NrVariable,1),Anova_Table.DFAutoError(NrVariable,1));
    % Anova_Table.FAuto(NrVariable,1)  = F_ValueAuto;
    % Anova_Table.PvalueAuto(NrVariable,1)  = P_ValueAuto;

    Anova_Table.SumSq(NrVariable,1) = nanmedian(cell2mat(ATableSumSq));
    Anova_Table.DF(NrVariable,1) = nanmedian(cell2mat(ATableDFs));
    Anova_Table.SumSqError(NrVariable,1) = nanmedian(cell2mat(ATableSumSqError));
    Anova_Table.DFError(NrVariable,1) = nanmedian(cell2mat(ATableDFsError));
    Anova_Table.MeanSq(NrVariable,1) = Anova_Table.SumSq(NrVariable,1)./Anova_Table.DF(NrVariable,1); 
    Anova_Table.MeanSqError(NrVariable,1) = Anova_Table.SumSqError(NrVariable,1)./Anova_Table.DFError(NrVariable,1); 

    F_Value = Anova_Table.MeanSq(NrVariable,1)./Anova_Table.MeanSqError(NrVariable,1);
    P_Value =  1- fcdf(F_Value,Anova_Table.DF(NrVariable,1),Anova_Table.DFError(NrVariable,1));
    Anova_Table.F(NrVariable,1) = F_Value;
    Anova_Table.Pvalue(NrVariable,1) = P_Value;
    Anova_Table.Coef(NrVariable,1) = nanmedian(coeffs);
    Anova_Table.RSquared(NrVariable,1) = nanmedian(cell2mat(RsquaredLoop));
end
end
%%
function [wij] = Morans(X,Y,MaxD) %#ok<INUSD>
if length(X) ~= length(Y)
    display('Fatal Error, unequal grid size')
    return
end
Size = length(X);
parfor i = 1:Size
    for j = 1:Size
        if i == j
            wij(i,j) = 0;  %#ok<*AGROW>
        else
            Dist =  log10((sqrt ((((X(i)-X(j))^2) + (Y(i)-Y(j))^2)))+1);
            wijtmp = (log10(MaxD)-Dist)./log10(MaxD);
            wijtmp(wijtmp<0) = 0;
            wij(i,j) = wijtmp;
        end
    end
end
end
%%
function AutoCorrelationCoef = AutoCorFunc(InVariation,Size)
wij = load('wij.mat');
parfor i = 1:Size
    Autot = 0;
    for j = 1:Size
        if i ~= j
            if isnan(InVariation(j))~= 1
                Autot = Autot + (wij.wij(i,j).*InVariation(j));  %#ok<*NODEF,*PFIIN>
            end
        end
    end
    if Autot > 0
        AutoCorrelationCoef(i,1) = Autot./nansum(wij.wij(i,:)); %#ok<*PFBNS,*AGROW>
    else
        AutoCorrelationCoef(i,1) = 0;
    end
end
end
%%
function Rsquared = RsquaredFunc(ExpecIn,ObservedIn)
% ExpecIn = ExpecIn - (min(ExpecIn));
% ObservedIn =ObservedIn - (min(ObservedIn));
List = find((isnan(ExpecIn)== 1));
List2 = find((isnan(ObservedIn)== 1));
List = [(reshape(List,[],1));(reshape(List2,[],1))];
ExpecIn(List) = [];
ObservedIn(List) = [];
meanY = mean(ObservedIn);
Size = length(ExpecIn);
for t = 1:1:Size
    ssres(t) = ((ExpecIn(t)-ObservedIn(t)).^2);
    sstot(t) =  ((ObservedIn(t)-meanY).^2);
end
Rsquared= {1- ((sum(ssres)) /(sum(sstot)))};
end

