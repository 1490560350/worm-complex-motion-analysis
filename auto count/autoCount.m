clc
close all
clear all
%the path where the file is located
trainPath='E:\';        
theFiles  = dir([trainPath '*.png']);
disp(length(theFiles));
count = 1;
train_num = length(theFiles);
sort_nat_name=sort_nat({theFiles.name}); 

%read data
HT = xlsread(fullfile(trainPath, 'analysis result', 'HeadTailReg.csv'));
Ph = xlsread(fullfile(trainPath, 'analysis result', 'Pharynx.csv'));
PP = xlsread(fullfile(trainPath, 'analysis result', 'PeakPoints.csv'));
IP = xlsread(fullfile(trainPath, 'analysis result', 'InflectionPoints.csv'));

head = HT(:,1:2);
tail = HT(:,3:4);
pharynx = Ph(:,1:2);

[m,n] = size(PP);
[p,q] = size(IP);

Data = zeros(100,5);

anglechange=[];
maxDist=[];
Datadist = [];

for k = 1:train_num
    fullFileName = sort_nat_name{k};
    I = imread([trainPath fullFileName]);
    count = count +1;

    set(0,'DefaultFigureVisible', 'off')
    
    PPnum = 0;
    IPnum = 0;
    
    for n1=1:n
        if isnan( PP(k,n1) )
            PPnum = PPnum+1;
        end
    end
    for n2=1:q
        if isnan( IP(k,n2) )
            IPnum = IPnum+1;
        end
    end
    
    PPnum1 = n - PPnum;
    IPnum1 = q - IPnum;
 
    dists = [];

    for i = 1:2:n-1
        if isnan( PP(k,i) )
            break;
        end
        
        flag = calculateFzdPosition(pharynx(k,1),pharynx(k,2),tail(k,1),tail(k,2),PP(k,i),PP(k,i+1));         
        dist = Dist( PP(k,i),PP(k,i+1),pharynx(k,1),pharynx(k,2),tail(k,1),tail(k,2) ); 
        dist = dist * flag;
        dists = [dists dist];       
    
    end
    b=length(dists);
    max = abs( dists(1) );
    index1 = 1;
    for t=1:b
        if max < abs( dists(t) )
            max = abs( dists(t) );
            index1 = t;
        end
    end

    maxDist = [maxDist dists(index1)];
    
     for i = 1:2:q-1
        if isnan( IP(k,i) )
            break;
        end
     end   
  
end

bodyBendnum = 0;
max1 = 0;
max2 = 0;
t = 1;
i = 2;
while i <= length(maxDist) - 2
    if sign(maxDist(t)) * sign(maxDist(i)) == 1 || (sign(maxDist(t)) * sign(maxDist(i)) == -1 && sign(maxDist(i+1)) * sign(maxDist(i)) == -1) 
        i = i + 1; 
    elseif sign(maxDist(t)) * sign(maxDist(i)) == -1 && sign(maxDist(i+1)) * sign(maxDist(i)) == 1
        t = i;
        [positiveCount, negativeCount, maxPositive, maxNegative, max1, max2] = deal(0, 0, 0, 0, 0, 0);

        while i + 2 <= length(maxDist) 
            if maxDist(i) > 0
                positiveCount = positiveCount + 1;
                if maxDist(i) > maxPositive
                    maxPositive = maxDist(i);
                    max1 = i;
                end
            elseif maxDist(i) < 0
                negativeCount = negativeCount + 1;  
                if abs(maxDist(i)) > abs(maxNegative)
                    maxNegative = maxDist(i);
                    max2 = i;
                end
            end

            if sign(maxDist(t)) * sign(maxDist(i+1)) == -1
                if (i+1-t > 1) && (sign(maxDist(t)) * sign(maxDist(i+2)) == -1)
                    break;
                end
            end
            i = i + 1;    
        end

        if i + 1 == length(maxDist) 
            if maxDist(i) > 0
                positiveCount = positiveCount + 1;
                if maxDist(i) > maxPositive
                    maxPositive = maxDist(i);
                    max1 = i;
                end
            elseif maxDist(i) < 0
                negativeCount = negativeCount + 1;  
                if abs(maxDist(i)) > abs(maxNegative)
                    maxNegative = maxDist(i);
                    max2 = i;
                end
            end

            if maxDist(i+1) > 0
                positiveCount = positiveCount + 1;
                if maxDist(i+1) > maxPositive
                    maxPositive = maxDist(i+1);
                    max1 = i+1;
                end
            elseif maxDist(i+1) < 0
                negativeCount = negativeCount + 1;  
                if abs(maxDist(i+1)) > abs(maxNegative)
                    maxNegative = maxDist(i+1);
                    max2 = i+1;
                end
            end
        end
        if i + 1 - t >= 2
            if sign(maxPositive) * sign(maxDist(t)) == 1 && (positiveCount >= negativeCount)
                bodyBendnum = bodyBendnum + 1;
                fprintf('Body bend detected in file: %s\n', sort_nat_name{max1});
            elseif sign(maxNegative) * sign(maxDist(t)) == 1 && (negativeCount >= positiveCount)
                bodyBendnum = bodyBendnum + 1;
                fprintf('Body bend detected in file: %s\n', sort_nat_name{max2});
            end
        end
        i = i + 1; 
    end      
end
fprintf('Dist: the number of body bends are %d\n', bodyBendnum);

