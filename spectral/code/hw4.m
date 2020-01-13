clear all
close all
%% Load data
load lynxdata
load sunspotdata

%Add search path for the provided custom MATLAB scripts
addpath ..\code

%% Choose method and data
%method='hoyw';
%method='music';
%method='minnorm';
method='esprit'; 

%y=sunspot;
y=lynx;
%y=loglynx;

%Remove mean and compute number of samples
yMean=mean(y);
y_m=y-yMean;
N=length(y_m);

%% Set method parameters

%Orders
n=4:2:40; %Number of sinusoids is n/2
K=length(n);

%User parameters
Lmax=fix(3*N/4); %Max of search interval for HOYW
mMax=fix(N/2); %Max of seach interval for the subspace methods

%% Explore in 2D (model order and method parameter)

switch method
    case 'hoyw'
        %Allocate matrices and cells
        W=cell(K,Lmax); %Storage for frequency vectors (varying lengths)
        rmseMat=nan(K,Lmax); %Storage for relative MSE values
        yHat=zeros(N,K,Lmax); %Storage for reconstructed signal vectors (fixed length N)
        for k=1:K %Loop for orders, n(k)
            for p=n(k)+1:Lmax %Loop for parameter "p" (here p=L=M)
                %Compute hoyw frequency estimates, compute amplitudes and phases
                %reconstruct the signal, compute relative MSE
                W{k,p} = hoyw(y_m,n(k),p,p);
                [~,rmseMat(k,p),yHat(:,k,p)] = lsa(y_m,W{k,p});
            end
        end
    case 'music'
        W=cell(K,mMax);
        rmseMat=nan(K,mMax);
        yHat=zeros(N,K,mMax); 
        for k=1:K
            for p=(n(k)+1):mMax %Loop for parameter "p" (here p=m)
                %Compute music frequency estimates, compute amplitudes and phases
                %reconstruct the signal, compute relative MSE
                W{k,p} = music(y_m,n(k),p);
                [~,rmseMat(k,p),yHat(:,k,p)] = lsa(y_m,W{k,p}); %Can be obtained from lsa() directly
            end
        end
    case 'minnorm'
        W=cell(K,mMax);
        rmseMat=nan(K,mMax);
        yHat=zeros(N,K,mMax); 
        for k=1:K
            for p=(n(k)+1):mMax %Loop for parameter "p" (here p=m)
                %Compute music frequency estimates, compute amplitudes and phases
                %reconstruct the signal, compute relative MSE
                W{k,p} = minnorm(y_m,n(k),p);
                [~,rmseMat(k,p),yHat(:,k,p)] = lsa(y_m,W{k,p}); %Can be obtained from lsa() directly
            end
        end
    case 'esprit'
        W=cell(K,mMax);
        rmseMat=nan(K,mMax);
        yHat=zeros(N,K,mMax); 
        for k=1:K
            for p=(n(k)+1):mMax %Loop for parameter "p" (here p=m)
                %Compute music frequency estimates, compute amplitudes and phases
                %reconstruct the signal, compute relative MSE
                W{k,p} = esprit(y_m,n(k),p);
                [~,rmseMat(k,p),yHat(:,k,p)] = lsa(y_m,W{k,p}); %Can be obtained from lsa() directly
            end
        end
    otherwise
        error(['No such method:' method])
end

%Extract the best relative MSE for each order (across all parameter choices)
rmseVec = min(rmseMat');

%Compute MSE (not relative), note that the model/methods assume zero mean, so this
%must be used in the order selection too
mseVec=(y_m'*y_m)*rmseVec/N;

%Compute the estimated order using BIC (sinorder() takes the number of
%real-valued sinusoids as the inputed order)
orders = sinorder(n,mseVec,N,4);

%% Plotting

%Plot the relative MSE for all parameter combinations and the resulting
%best fit to the data
[bestMin,bestInd]=min(rmseMat(:));
[i,j]=ind2sub(size(rmseMat),bestInd);
if strcmp(method,'hoyw') %Plot for hoyw (only different titles and labels)
    figure,h=axes;
    imagesc(1:Lmax,n/2,rmseMat);
    set(h,'YTick',n/2)
    xlabel('Parameter L=M (>n)')
    ylabel('Number of sinusoids (n/2)')
    colorbar
    
    figure, plot([y_m yHat(:,i,j)]+yMean) %Add the mean again
    title(['Sinusoids: ' num2str(n(i)/2) ', L=M= ' num2str(j) ', rel. MSE = ' num2str(bestMin)])
    xlabel('Samples')
    ylabel('Signal amplitude')
    legend('Data','Fit')
else %Plot for the subspace methods (only different titles and labels)
    figure,h=axes;
    imagesc(1:mMax,n/2,rmseMat)
    set(h,'YTick',n/2)
    xlabel('Parameter m (>n)')
    ylabel('Number of sinusoids (n/2)')
    colorbar
    
    figure, plot([y_m yHat(:,i,j)]+yMean) %Add the mean again
    title(['Sinusoids: ' num2str(n(i)/2) ', m = ' num2str(j) ', rel. MSE = ' num2str(bestMin)])
    xlabel('Samples')
    ylabel('Signal amplitude')
    legend('Data','Fit')
end

%Plot the best solution (frequencies and amplitudes)
betaHat=lsa(y_m,W{i,j});
figure,stem(W{i,j},abs(betaHat))
title('Estimated frequencies and amplitudes for best fit')
ylabel('Amplitude')
xlabel('Normalized frequency')
xlim([0 pi])


%Plot the BIC solution (frequencies and amplitudes)
k=find(n==2*orders(4));
[~,BICind]=min(rmseMat(k,:));
betaHat=lsa(y_m,W{k,BICind});
figure,stem(W{k,BICind},abs(betaHat))
title('Estimated frequencies and amplitudes for the BIC choice')
ylabel('Amplitude')
xlabel('Normalized frequency')
xlim([0 pi])


%Plot the best relative MSE vs model order (for order selection)
figure,h=axes;
stem(n/2,rmseVec);
hold on
stem(orders(4),rmseVec(k),'g')
title('Best relative MSE vs model order')
xlabel('Number of sinusoids (n/2)')
ylabel('rel. MSE')
set(h,'XTick',n/2)
legend('Best rel. MSE','BIC choice')