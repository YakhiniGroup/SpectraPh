function [path,clusterid,logLikelihood,confidenceLifts]=ViterbiClustering(V, singleManifoldOpt, VisDistances)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% [path,logLikelihood]=ViterbiClustering(V)
%%%%% Takes as input a k dimensional embedding of two equi-sized manifolds V(2n x k)
%%%%% We assume 2 disjoint curves and decide per loci to which curve it belongs to.
%%%%%
%%%%% We utilize hmmViterbiC from https://github.com/probml/pmtk3's
%%%%% mex implementation of Viterbi. Credit goes to them.
%%%%%
%%%%% Output: 
%%%%% path          = the most likely 
%%%%% logLikelihood = 
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('singleManifoldOpt', 'var')
        singleManifoldOpt = true;
    end
    if ~exist('VisDistances', 'var')
        VisDistances = range(size(V))==0 && isequal(V,V');
    end
    metric='euclidean';
    n=size(V,1)/2;

    %Identify trellis similarity matrix coords
    top2top=num2cell([(1:n-1)' (2:n)'],2);
    top2bot=num2cell([(1:n-1)' ((n+2):size(V,1))'],2);
    bot2bot=num2cell([((n+1):size(V,1)-1)' ((n+2):size(V,1))'],2);
    bot2top=num2cell([((n+1):size(V,1)-1)' (2:n)'],2);
    
    %Compute outgoing probabilities
    if (VisDistances)
        topProbs = [cellfun(@(x) V(x(1),x(2))+eps, top2top)';...
                    cellfun(@(x) V(x(1),x(2))+eps, top2bot)'];
        botProbs = [cellfun(@(x) V(x(1),x(2))+eps, bot2top)';...
                    cellfun(@(x) V(x(1),x(2))+eps, bot2bot)'];
    else
        topProbs = [cellfun(@(x) (pdist2(V(x(1),:),V(x(2),:),metric))+eps, top2top)';...
                    cellfun(@(x) (pdist2(V(x(1),:),V(x(2),:),metric))+eps, top2bot)'];
        botProbs = [cellfun(@(x) (pdist2(V(x(1),:),V(x(2),:),metric))+eps, bot2top)';...
                    cellfun(@(x) (pdist2(V(x(1),:),V(x(2),:),metric))+eps, bot2bot)'];
    end
    topProbs = topProbs-min([topProbs(:);0]); % remove negatives
    botProbs = botProbs-min([botProbs(:);0]); % remove negatives
    normTop =  1-(topProbs./[sum(topProbs);sum(topProbs)]);
    normBot =  1-(botProbs./[sum(botProbs);sum(botProbs)]);
    
    TransitionPr=[normTop; normBot];
    logTr=log(TransitionPr);
    if(singleManifoldOpt)
        [clusterid, logLikelihood,confidenceLifts]=DPOptPathWithTies(logTr);
    else %optimize both paths
        PathProbs=[logTr(1,:)+logTr(4,:); logTr(2,:)+logTr(3,:)];        
        clusterid=(PathProbs(1,:)>PathProbs(2,:))+1;
        logLikelihood=sum(PathProbs(sub2ind(size(PathProbs),clusterid,1:n-1)));
        clusterid(PathProbs(1,:)==PathProbs(2,:))=0;
        confidenceLifts=PathProbs(1,:)'-PathProbs(2,:)';
    end
    path=zeros(1,n-1);
    path(clusterid==1)=find(clusterid==1);
    path(clusterid==2)=find(clusterid==2)+n;
end
