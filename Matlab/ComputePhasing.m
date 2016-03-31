%% Phase embeddings
%% Configuration
addpath(horzcat(pwd,'\Helpers\'));
param.cov=1;
param.diagrem=1;
param.embedd=1;
chromosomes=23;
%%
unphasedCount=zeros(1,chromosomes); phasedCount=zeros(1,chromosomes);
results=cell(1,chromosomes); embeddTypeMDscale=false(1,chromosomes); numComponents=zeros(1,chromosomes);
coverage=cell(1,chromosomes); compsizes=cell(1,chromosomes);
weightedContrib=zeros(1,chromosomes); wresults=zeros(1,chromosomes); wstress=zeros(1,chromosomes);
for chrid = 1:chromosomes
    filename=horzcat('ChrCovEmbedd.',num2str(chrid),'.',num2str(param.resolution),'.',num2str(param.cov),num2str(param.diagrem),...
        num2str(param.embedd),'.mat');
    while(~exist(filename,'file'))
        display(horzcat('Waiting on file ',filename))
        pause(10)
    end
    display(horzcat('Loading ',filename))
    load(filename);
    contigs=load(horzcat('contigMap.',num2str(chrid),'.',num2str(param.resolution),'.csv'));
    unphasedCount(chrid)=sum(cellfun(@length,ensamble.coords(~ensamble.embeded)));
    phasedCount(chrid)=sum(cellfun(@length,ensamble.coords(ensamble.embeded)));
    compsizes{chrid}=cellfun(@length,ensamble.mds(~(cellfun(@isempty,ensamble.mds))));
    coverage{chrid}=zeros(1,sum(~(cellfun(@isempty,ensamble.mds))));
    results{chrid}=zeros(1,sum(~(cellfun(@isempty,ensamble.mds))));

    compid=0;
    for comp=find(~(cellfun(@isempty,ensamble.mds)))
        compid=compid+1;
        n1=size(ensamble.mds{comp},1)/2;
        ord=randi(2,1,n1)-1;
        prm=[(1:n1) + n1.*(ord) (1:n1) + n1.*(~ord)];
    
        [path,cids,loglik,confidence]=ViterbiClustering(ensamble.mds{comp}(prm,:),false);
        prmcids=cids; prmcids(diff(ord)~=0)=~(cids(diff(ord)~=0)-1)+1;
        prmconfidence=confidence'; prmconfidence(diff(ord)~=0)=-prmconfidence(diff(ord)~=0);
        coverage{chrid}(compid)=(length(ensamble.coords{comp})/(2*length(contigs)));
        wstress(chrid)=wstress(chrid)+coverage{chrid}(compid)*ensamble.stress(comp);
        weightedContrib(chrid)=weightedContrib(chrid)+coverage{chrid}(compid);
        numComponents(chrid)=numComponents(chrid)+1;
        
        Results.Chromosome(chrid).HTblockIds(compid) = ensamble.coords(comp);
        Results.Chromosome(chrid).Haplotype(compid) = prmcids;
        
        display(horzcat('Chr ',num2str(chrid), ' Component #',num2str(comp),'H = '));
        display(prmcids);
    end
    weightedContrib(chrid)=weightedContrib(chrid)-numComponents(chrid)/length(tG{chrid,chrid});
end

Results
