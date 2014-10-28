load datasets;
A=datasets{1, 1}
lastTenth=1/10*numDays;
forwardSelectionPurchaseDayPrices=dailyPrices(numDays-(lastTenth+1), :)
retmat=makeReturnsMatrix(A);