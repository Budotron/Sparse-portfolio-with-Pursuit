horizon=[0, 1, 3, 12]
for i=1:9
        oneOnNvec=[realOneOnNReturns_month(i), realOneOnNReturns_3months(i), realOneOnNReturns_year(i)];
        MarkowitzVec=[realMarkowitzReturns_month(i), realMarkowitzReturns_3months(i), realMarkowitzReturns_year(i)];
        %for card=1:10
            %MPV=[realMPReturns_month(i, card), realMPReturns_3months(i, card), realMPReturns_year(i, card)];
            MPV1=[realMPReturns_month(i, 5), realMPReturns_3months(i, 5), realMPReturns_year(i, 5)];
            OMPV1=[realOMPReturns_month(i, 5), realOMPReturns_3months(i, 5), realOMPReturns_year(i, 5)];
            LSOMPV1=[realLSOMPReturns_month(i, 5), realLSOMPReturns_3months(i, 5), realLSOMPReturns_year(i, 5)];
            thrV1=[realthrReturns_month(i, 5), realthrReturns_3months(i, 5), realthrReturns_year(i, 5)];
            figure (i)  
            cc=hsv(10);
            hold all
            plot(horizon, [0, oneOnNvec], 'r*')
            plot(horizon, [0, MarkowitzVec], 'rv-')
            plot(horizon, [0, MPV1], 'c')
            plot(horizon, [0, OMPV1], 'b')
            plot(horizon, [0, LSOMPV1], 'm')
            plot(horizon, [0, thrV1], 'k')
            legend({'1/N', 'MM', 'MP', 'OMP', 'LSOMP', 'thr'})
            xlabel('Return Horizon(M)')
            ylabel('Realized Returns')
            title(sprintf('Realized returns of 50-sparse porfolios purchased in period %d ', i-1))
        %end
end


