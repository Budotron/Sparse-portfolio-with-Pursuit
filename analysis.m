clc;
% load oneOnNPortfolio;
% load ExpectedOneOnNReturns
% load realOneOnNReturns_month
% load realOneOnNReturns_3months
% load realOneOnNReturns_year
% load oneOnNVar
% 
% load MarkowitzPortfolio
% load realMarkowitzReturns_month
% load realMarkowitzReturns_3months
% load realMarkowitzReturns_year
% load MarkowitzVar
% 
% load backwardsPort_all
% load BackwardsVar   
% load realBackwardsReturns_month
% load realBackwardsReturns_3months
% load realBackwardsReturns_year
% load backwardsSelectedAssets
% 
% load FowardsPort_all
% load ForwardsVar
% load realForwardsReturns_month
% load realForwardsReturns_3months
% load realForwardsReturns_year
% load ForwardsSelectedAssets
% 
% load OMPX_all
% load OMPVar
% load realOMPReturns_month
% load realOMPReturns_3months
% load realOMPReturns_year
% load OMPSelectedAssets
% 
% load MPX_all
% load MPVar
% load realMPReturns_month
% load realMPReturns_3months
% load realMPReturns_year
% load MPSelectedAssets
% 
% load LSOMPX_all
% load LSOMPVar
% load realLSOMPReturns_month
% load realLSOMPReturns_3months
% load realLSOMPReturns_year
% load LSOMPSelectedAssets
% 
% load thrX_all
% load thrVar
% load realthrReturns_month
% load realthrReturns_3months
% load realthrReturns_year
% load thrSelectedAssets
% 
load names; load datesSets;
periods=0:1:9
periods2=0:10;
% sanity=0
% for i=1:9
%     sanity=sanity+1
% %     tempdates=datesSets{i};
% %     dates=datevec(zeros(length(tempdates),1));
% %     for j=1:length(tempdates)
% %         dates=datestr(datenum('tempdates(j)', 'yyyymmdd'))
% %     end
%     figure(i)
%     hold on
%     plot(periods, [0 realOneOnNReturns_year], '*r');
%     plot(periods, [0, realMarkowitzReturns_year])
%     plot(periods, realMPReturns_year(i, :))
%     legend({'1/N', 'MM', 'MP'});
%     xlabel('time horizon: one year')
%     ylabel('realized return(%)')
%     title(sprintf('RR at cardinality %d', i*10))
%     
% end
    
    for cardinality=1:10
        figure(cardinality)
        hold on
        plot(periods, [0, realOneOnNReturns_year], '*r')
        plot(periods, [0, realMarkowitzReturns_year], 'g')
        plot(periods, [0, realMPReturns_year(:, cardinality)'], 'c')
        plot(periods, [0, realOMPReturns_year(:, cardinality)'], 'b')
        plot(periods, [0, realLSOMPReturns_year(:, cardinality)'], 'm')
        plot(periods, [0, realthrReturns_year(:, cardinality)'], 'k')
        plot(periods, [0, realForwardsReturns_year(:, cardinality)'], '-.r')
        plot(periods, [0, realBackwardsReturns_year(:, cardinality)'], ':r')
        legend({'1/N', 'MM', 'MP', 'OMP', 'LSOMP', 'thr', 'Forwards', 'Backwards'})
        xlabel('time horizon: one year')
        ylabel('realized return(%)')
        title(sprintf('RR at cardinality %d', cardinality*10));
    end