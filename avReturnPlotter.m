avOneOnN_month=sum(realOneOnNReturns_month)/length(realOneOnNReturns_month);
avMarkowitz_month=sum(realMarkowitzReturns_month)/length(realMarkowitzReturns_month);
avForward_month=sum(realForwardsReturns_month)/9;
avBackward_month=sum(realBackwardsReturns_month)/9;
avOMP_month=sum(realOMPReturns_month)/9;
avMP_month=sum(realMPReturns_month)/9;
avLSOMP_month=sum(realLSOMPReturns_month)/9;
avth_month=sum(realthrReturns_month)/9;
card=10:10:100;
figure
hold on
plot(card, avMP_month, 'c');
plot(card, avOneOnN_month, 'r*');
plot(card, avMarkowitz_month, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of MPS portfolios (cardinality K) vs unconstrainedbenchmarks over 1 month')
h=legend({'MPS', '1/N', 'MM'})

figure
hold on
plot(card, avOMP_month, 'b');
plot(card, avOneOnN_month, 'r*');
plot(card, avMarkowitz_month, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of 0MPS portfolios (cardinality K) vs unconstrainedbenchmarks over 1 month')
h=legend({'OMPS', '1/N', 'MM'})

figure
hold on
plot(card, avLSOMP_month, 'm');
plot(card, avOneOnN_month, 'r*');
plot(card, avMarkowitz_month, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of LSS portfolios (cardinality K) vs unconstrainedbenchmarks over 1 month')
h=legend({'LSS', '1/N', 'MM'})

figure
hold on
plot(card, avth_month, 'k');
plot(card, avOneOnN_month, 'r*');
plot(card, avMarkowitz_month, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of TS portfolios (cardinality K) vs unconstrainedbenchmarks over 1 month')
h=legend({'TS', '1/N', 'MM'})

avOneOnN_3months=sum(realOneOnNReturns_3months)/length(realOneOnNReturns_3months);
avMarkowitz_3months=sum(realMarkowitzReturns_3months)/length(realMarkowitzReturns_3months);
avForward_3months=sum(realForwardsReturns_3months)/9;
avBackward_3months=sum(realBackwardsReturns_3months)/9;
avOMP_3months=sum(realOMPReturns_3months)/9;
avMP_3months=sum(realMPReturns_3months)/9;
avLSOMP_3months=sum(realLSOMPReturns_3months)/9;
avth_3months=sum(realthrReturns_3months)/9;
card=10:10:100;
figure
hold on
plot(card, avMP_3months, 'c');
plot(card, avOneOnN_3months, 'r*');
plot(card, avMarkowitz_3months, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of MPS portfolios (cardinality K) vs unconstrainedbenchmarks over 3 months')
h=legend({'MPS', '1/N', 'MM'})

figure
hold on
plot(card, avOMP_3months, 'b');
plot(card, avOneOnN_3months, 'r*');
plot(card, avMarkowitz_3months, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of 0MPS portfolios (cardinality K) vs unconstrainedbenchmarks over 3 months')
h=legend({'OMPS', '1/N', 'MM'})

figure
hold on
plot(card, avLSOMP_3months, 'm');
plot(card, avOneOnN_3months, 'r*');
plot(card, avMarkowitz_3months, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of LSS portfolios (cardinality K) vs unconstrainedbenchmarks over 3 months')
h=legend({'LSS', '1/N', 'MM'})

figure
hold on
plot(card, avth_3months, 'k');
plot(card, avOneOnN_3months, 'r*');
plot(card, avMarkowitz_3months, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of TS portfolios (cardinality K) vs unconstrainedbenchmarks over 3 months')
h=legend({'TS', '1/N', 'MM'})

avOneOnN_year=sum(realOneOnNReturns_year)/length(realOneOnNReturns_year);
avMarkowitz_year=sum(realMarkowitzReturns_year)/length(realMarkowitzReturns_year);
avForward_year=sum(realForwardsReturns_year)/9;
avBackward_year=sum(realBackwardsReturns_year)/9;
avOMP_year=sum(realOMPReturns_year)/9;
avMP_year=sum(realMPReturns_year)/9;
avLSOMP_year=sum(realLSOMPReturns_year)/9;
avth_year=sum(realthrReturns_year)/9;
card=10:10:100;
figure
hold on
plot(card, avMP_year, 'c');
plot(card, avOneOnN_year, 'r*');
plot(card, avMarkowitz_year, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of MPS portfolios (cardinality K) vs unconstrainedbenchmarks over 12 months')
h=legend({'MPS', '1/N', 'MM'})

figure
hold on
plot(card, avOMP_year, 'b');
plot(card, avOneOnN_year, 'r*');
plot(card, avMarkowitz_year, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of 0MPS portfolios (cardinality K) vs unconstrainedbenchmarks over 12 months')
h=legend({'OMPS', '1/N', 'MM'})

figure
hold on
plot(card, avLSOMP_year, 'm');
plot(card, avOneOnN_year, 'r*');
plot(card, avMarkowitz_year, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of LSS portfolios (cardinality K) vs unconstrainedbenchmarks over 12 months')
h=legend({'LSS', '1/N', 'MM'})

figure
hold on
plot(card, avth_year, 'k');
plot(card, avOneOnN_year, 'r*');
plot(card, avMarkowitz_year, 'rv');
xlabel('Cardinality of constructed portfolio')
ylabel('Average realized return (10 years)')
title('Average return of TS portfolios (cardinality K) vs unconstrainedbenchmarks over 12 months')
h=legend({'TS', '1/N', 'MM'})







