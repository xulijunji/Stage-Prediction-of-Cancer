function test_label_predict=kNN(train_data,train_label,test_data,k_value)
% k-NN classifier
%   input data :
%        k_value is the value of "K"
%        train_data: row is sample, column is gene
%        train_label: is column vector
%        test_data: row is sample, column is gene
%   output data:
%        test_label_predict is row vector,equivalent length of "number of test samples"
[nTrainSample,nGene]=size(train_data);

[nTestSample,nGene]=size(test_data);
number_Class=length(unique(train_label));

test_label_predict=zeros(nTestSample,1);

for id_TestSample=1:nTestSample
    x_TestSample=test_data(id_TestSample,:);
    dis_TrainSample_label=zeros(nTrainSample,3);
    dis_TrainSample_label(:,3)=train_label;
    dis_TrainSample_label(:,2)=[1:nTrainSample]';    
    
    for id_TrainSample=1:nTrainSample
        y_TrainSample=train_data(id_TrainSample,:);
        dis_TrainSample_label(id_TrainSample,1)=norm(x_TestSample-y_TrainSample,2);
    end
    dis_TrainSample_label=sortrows(dis_TrainSample_label,1);
    dis_TrainSample_label=dis_TrainSample_label(1:k_value,:); % just considering samples with "K-th" smallest distance 

    nSample_each_label=zeros(number_Class,1);
    for k=1:number_Class
        nSample_each_label(k,1)=length(find(dis_TrainSample_label(:,3)==k));
    end
    max_number_in_each_label=max(nSample_each_label);
    
    k=find(nSample_each_label<max_number_in_each_label);
    [y,k]=setdiff(dis_TrainSample_label(:,3),k);
    
    dis_TrainSample_label=dis_TrainSample_label(k,:);   % delete the labels,whose sample number less than "max_number"
      
    dis_TrainSample_label=sortrows(dis_TrainSample_label,1);
    
    test_label_predict(id_TestSample,1)=dis_TrainSample_label(1,3);
end
%
return
% function is OVER

