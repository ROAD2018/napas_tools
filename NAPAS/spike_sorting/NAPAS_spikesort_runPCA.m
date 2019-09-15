% 2019-09-15 Dylan Royston
% Adapted 2018-04-09 Dylan Royston
%
% Copied from hst2/Analysis Code/unitID/PCA sorting/testsort.m
%
% Function to perform PCA on provided spike snippets
% Copied for cleanup and optimization for Somatomapping data
%
% === INPUTS === 
% InitParams:   initial parameters (mu, sigma)
% Spikes:       spike snippets (time x num_snippets)
%
%
% === OUTPUTS ===  
% muK:          cluster means
% Sigma:        cluster covariances
% ppi:          ?
% K:            ?
% sorted:       ?
%
% === UPDATES ===
%
%%

function [muK,Sigma,ppi,K,sorted] = NAPAS_spikesort_runPCA(InitParams,Spikes, peak_loc, color_flag)

colorspec = jet(48);
% calculate average spike snippet
mu = zeros(48,1);

for n = 1:size(Spikes,2)
    mu = mu + double(Spikes(:,n));
end
mu = mu./size(Spikes,2);

% calculate covariance of spike snippets
sig = zeros(48,48);

for n = 1:size(Spikes,2)
    sig = sig + (double(Spikes(:,n))-mu)*((double(Spikes(:,n))-mu)');
end
sig = sig./size(Spikes,2);

[U S] = eig(sig);

% % optional figure, plots first three eigenvectors
% figure(); %for some reason the end of the matrix had the higest eigenvalues
% plot(1:48,U(:,48),'r',1:48,U(:,47),'g',1:48,U(:,46),'b');xlabel('Sample Number');
% ylabel('Magnitude');title('First Three Eigenvectors');
% legend('1st','2nd','3rd');

z = zeros(2,size(Spikes,2));
U2(:,1) = U(:,48);  %puts eigenvectors in proper order
U2(:,2) = U(:,47);

% figure();hold on;
for n = 1:size(Spikes,2)
    z(:,n) = U2'*(double(Spikes(:,n))-mu);
    %    plot(z(1,n),z(2,n),'k.');
end
% xlabel('PC1 score');ylabel('PC2 score');title('PC Score Scatter Plot');
% hold off;


% lik = zeros(1,8);

figure(1); hold off;

if color_flag == 1
    scatter(z(1,:),z(2,:),2,colorspec(peak_loc, :), 'filled'); hold on;
    set(gca, 'Color', 'k');
else
    
    %     for n = 1:size(Spikes,2)   %plots spikes in 2-d PC
    %         plot(z(1,n),z(2,n),'k.'); hold on;
    %     end
    plot(z(1,:),z(2,:),'k.'); hold on;
end
set(gcf, 'Position', [1929 373 852 623]);



% default format is to enter cluster data into Matlab command window, requires main Matlab window to be on top
commandwindow;
K = input('How many units do you see? ');


while isempty(K) || K > 5 || K < 1 || rem(K,1) ~= 0
    if K == 9
        muK = [];
        Sigma = [];
        % 2018-04-17 Royston: have to assign all output values (as empty) to avoid error
        ppi = [];
        K = 9;
        sorted = [];
        hold off
        return
    else
        K = input('Must have 1-5 units, or 9 to go back to last sortable channel: ');
    end
end

% user input to click on cluster centers in PCA window
if K > 1
    figure(1);
%     [x,y] = ginput(K);
    [x,y] = ginputc(K, 'Color', [1 1 1]);
    KParams.mu = [x';y'];
else
    KParams.mu =  InitParams.mu(:,1:K);
end% IF, K > 1

KParams.pi(1:K) = 1/K; %Initializes values for GMM
for i = 1:K
    KParams.sigma(:,:,i) = InitParams.Sigma;
end% FOR, i=1:K

% attempt to perform EM algorithm around cluster centers
try
    disp('*** RUNNING EM ***');
    [muK,Sigma,ppi] = FUNC_run_GaussianMM(KParams,z); %Runs EM algorithm
catch ME
   disp('Error sorting:')
   disp(ME)
   disp('Try again.')
   K = 9;
   muK = [];
   Sigma = [];
   sorted = [];
   hold off
   return
end

%generate probability distributioons for ellipse plotting
disp('*** GENERATING ELLIPSES ***');
p = zeros(size(Spikes,2),K);
for i = 1:K
    p(:,i) = mvnpdf(z',muK(:,i)',Sigma(:,:,i)).*ppi(K);
end

[~, sorted] = max(p');


disp('*** PLOTTING ELLIPSES ***');
for k = 1:K     %Plots means and STD ellipses over data
    func_plotEllipse(squeeze(muK(:,k)),squeeze(Sigma(:,:,k)));
end

% 2018-04-17 Royston: re-plot PCA space points with cluster assignments for visual inspection

disp('** REPLOTTING CLUSTERS ***');

if color_flag == 1
    cluster_colors = {'w', 'g', 'b', 'c', 'm'};
else
    cluster_colors = {'k', 'g', 'b', 'c', 'm'};
end

% 2018-08-09 Royston: replaced the inefficient "plot every individual point" method with a sensible scatter
% tic;
for k = 1 : K
    current_points = find(sorted == k);
    scatter(z(1,current_points), z(2, current_points), 2, cluster_colors{k}, 'filled');
end


hold off;

end% FUNCTION

