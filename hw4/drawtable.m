tables=zeros(5,5);
cd order1
run tvdm1.m
tables(:,1)=n(1:5)';
tables(:,2)=error4(1:5)';
tables(2:5,3)=rate4(1:4)';
cd ../order2
run tvdm1.m
tables(:,4)=error4(1:5)';
tables(2:5,5)=rate4(1:4)';