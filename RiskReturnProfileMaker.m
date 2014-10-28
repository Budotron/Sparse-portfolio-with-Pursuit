
avMarVar=sum(MarkowitzVar)/20;
avMarRet3months=sum(realMarkowitzReturns_year)/20;
avOneOnNvar=sum(oneOnNVar)/20;
avOMPVar3months=sum(OMPVar)/20;
avOneOnNRet3months=sum(realOneOnNReturns_year )/20;
avOMPRet3months=sum(realOMPReturns_year)/20;
avMPVar3months=sum(MPVar)/20;
avMPRet3months=sum(realMPReturns_year)/20;
avLSOMPVar3months=sum(LSOMPVar)/20;
avLSOMPRet3months=sum(realLSOMPReturns_year)/20;
avthrVar3months=sum(thrVar)/20;
avthrRet3months=sum(realthrReturns_year)/20;
avForVar3months=sum(ForwardsVar)/20;
avForRet3months=sum(realForwardsReturns_year)/20;
avBackVar3months=sum(BackwardsVar)/20;
avBackRet3months=sum(realBackwardsReturns_year)/20;
for i=1:3
    figure
    hold on
    plot(avMarVar, avMarRet3months, 'rv', 'LineWidth', 2)
    plot(avOneOnNvar, avOneOnNRet3months, 'r*','LineWidth', 2)
    plot(avOMPVar3months(:, i), avOMPRet3months(:, i), 'b*','LineWidth', 2)
    plot(avMPVar3months(:, i), avMPRet3months(:, i), 'c*','LineWidth', 2)
    plot(avLSOMPVar3months(:, i), avLSOMPRet3months(:, i), 'm*','LineWidth', 2)
    plot(avthrVar3months(:, i), avthrRet3months(:, i), 'k*','LineWidth', 2)
    plot(avForVar3months(:, i), avForRet3months(:, i), 'g*','LineWidth', 2)
    plot(avBackVar3months(:, i), avBackRet3months(:, i), 'rx','LineWidth', 2)
    xlabel('Risk')
    ylabel('Return')
    title('Average (100days) Risk-ReturnProfiles')
end

